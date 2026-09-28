#!/usr/bin/env Rscript

# STEP 01c — Extended robustness analysis for the GEDI NMF–k-means hierarchy
#
# Run only after STEP 01b has completed. This script reads the definitive
# STEP 01b checkpoint/model and produces eight numbered CSV outputs:
#   01 NMF restart stability
#   02 complete-procedure subsample stability
#   03 component-specific k-selection stability
#   04 sampling-design sensitivity
#   05 pairwise class-profile distinctness
#   06 final-class silhouette widths
#   07 blocked spatial holdout
#   08 training-versus-population prevalence (sampling-design audit)
#
# Every test is restartable: an existing output is retained unless overwrite
# is TRUE. The expensive conventional sample is also checkpointed.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(RcppML)
  library(clue)
  library(cluster)
})

cfg <- list(
  step01a_dir = "~/data/step01a_classification_prescreen",
  step01b_dir = "~/data/step01b_final_classification",
  input_files = file.path(
    "~/data",
    sprintf("GEDI_RFLD_sampling_frame_%d.parquet", 2019:2023)
  ),
  output_dir = "~/data/step01c_extended_robustness",
  seed = 42L,
  overwrite = FALSE,
  
  # Tests 1–3
  comparison_n = 100000L,
  repetitions = 20L,
  subsample_fraction = 0.80,
  rank_restarts_per_repetition = 1L,
  
  # Test 4: conventional sampling balanced only by year and forest class
  conventional_training_n = 100000L,
  evaluation_n = 30000L,
  
  # Tests 6–7
  silhouette_n = 10000L,
  spatial_training_n = 100000L,
  spatial_evaluation_n = 30000L,
  spatial_cell_size_m = 10000,
  spatial_holdout_fraction = 0.20,
  
  # Fitting controls (must remain compatible with Steps 01a/01b)
  nmf_maxit = 200L,
  nmf_tol = 1e-5,
  k_candidates = 1:6,
  kmeans_nstart = 50L,
  kmeans_iter_max = 500L,
  kmeans_algorithm = "Lloyd",
  rcppml_threads = max(1L, parallel::detectCores(logical = FALSE) - 1L)
)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(cfg$output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
RcppML::setRcppMLthreads(cfg$rcppml_threads)
set.seed(cfg$seed)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), ..., "\n", sep = "")
  flush.console()
}

out_file <- function(name) file.path(cfg$output_dir, name)
must_run <- function(path) cfg$overwrite || !file.exists(path)
rh_names <- paste0("rh", 0:98)

model_path <- file.path(cfg$step01b_dir, "models", "nmf_shape_height_kmeans_model.rds")
checkpoint_path <- file.path(cfg$step01b_dir, "checkpoints", "03_model_and_training_classes.rds")
assignment_dir <- file.path(cfg$step01b_dir, "assignments_dataset")
required_paths <- c(model_path, checkpoint_path, cfg$input_files)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Missing required inputs:\n", paste(missing_paths, collapse = "\n"))

reference_model <- readRDS(model_path)
o <- readRDS(checkpoint_path)
rank_selected <- as.integer(reference_model$selected_rank)
fixed_k <- vapply(reference_model$kmeans, function(z) as.integer(z$k), integer(1))
fixed_k <- fixed_k[order(as.integer(names(fixed_k)))]

if (rank_selected != length(fixed_k)) {
  stop("Reference model rank and number of component-specific k values disagree")
}

shape_normalize <- function(rh) {
  lo <- apply(rh, 1L, min, na.rm = TRUE)
  hi <- rh[, 99L]
  span <- hi - lo
  if (any(!is.finite(span) | span <= 0)) stop("Invalid RH range in retained profiles")
  z <- (rh - lo) / span
  z[z < 0 & z > -1e-6] <- 0
  z[z > 1 & z < 1 + 1e-6] <- 1
  pmin(pmax(z, 0), 1)
}

as_sparse <- function(x) methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
rowmax <- function(x) max.col(x, ties.method = "first")

ari <- function(x, y) {
  z <- table(x, y)
  c2 <- function(a) a * (a - 1) / 2
  n <- sum(z)
  if (n < 2) return(NA_real_)
  a <- sum(c2(rowSums(z))); b <- sum(c2(colSums(z))); observed <- sum(c2(z))
  expected <- a * b / c2(n)
  denominator <- 0.5 * (a + b) - expected
  if (denominator == 0) return(if (identical(as.integer(x), as.integer(y))) 1 else 0)
  (observed - expected) / denominator
}

nmi <- function(x, y) {
  z <- table(x, y)
  pxy <- z / sum(z); px <- rowSums(pxy); py <- colSums(pxy)
  ii <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(pxy[ii] * log(pxy[ii] / (px[ii[, 1L]] * py[ii[, 2L]])))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 && hy == 0) return(1)
  if (hx == 0 || hy == 0) return(0)
  mi / sqrt(hx * hy)
}

choose_knee <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) <= 2L) return(x[[which.min(y)]])
  xs <- (x - min(x)) / diff(range(x)); ys <- (y - min(y)) / diff(range(y))
  # Maximum perpendicular distance from the line joining the endpoints.
  x1 <- xs[1L]; y1 <- ys[1L]; x2 <- xs[length(xs)]; y2 <- ys[length(ys)]
  distance <- abs((y2-y1)*xs - (x2-x1)*ys + x2*y1 - y2*x1) /
    sqrt((y2-y1)^2 + (x2-x1)^2)
  x[[which.max(distance)]]
}

fit_best_nmf <- function(A, k, restarts, offset) {
  best <- NULL; best_mse <- Inf; runs <- vector("list", restarts)
  for (j in seq_len(restarts)) {
    fit <- RcppML::nmf(
      A, k = k, tol = cfg$nmf_tol, maxit = cfg$nmf_maxit,
      seed = cfg$seed + offset + j, verbose = FALSE
    )
    mse <- RcppML::mse(A, fit$w, fit$d, fit$h)
    runs[[j]] <- data.table(restart = j, mse = mse, iterations = fit$iter)
    if (mse < best_mse) { best <- fit; best_mse <- mse }
  }
  list(model = best, mse = best_mse, runs = rbindlist(runs))
}

match_components <- function(reference_w, candidate_w) {
  correlation <- cor(reference_w, candidate_w)
  correlation[!is.finite(correlation)] <- -1
  # solve_LSAP() requires nonnegative entries. Adding the same constant to
  # every cell preserves the assignment that maximizes the correlation sum.
  assignment_score <- correlation - min(correlation, na.rm = TRUE)
  reference_to_candidate <- as.integer(
    clue::solve_LSAP(assignment_score, maximum = TRUE)
  )
  candidate_to_reference <- integer(length(reference_to_candidate))
  candidate_to_reference[reference_to_candidate] <- seq_along(reference_to_candidate)
  matched <- correlation[cbind(seq_along(reference_to_candidate), reference_to_candidate)]
  list(
    reference_to_candidate = reference_to_candidate,
    candidate_to_reference = candidate_to_reference,
    matched_correlations = matched,
    mean_matched_correlation = mean(matched)
  )
}

height_features <- function(raw_rh, center = NULL, scale = NULL) {
  x <- cbind(rh98 = raw_rh[, 99L], negative_depth = pmax(0, -raw_rh[, 1L]))
  if (is.null(center)) center <- colMeans(x)
  if (is.null(scale)) scale <- apply(x, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  list(x = z, center = center, scale = scale)
}

fit_kmeans_free <- function(component, x, offset = 0L) {
  clusters <- integer(nrow(x)); fits <- list(); rows <- list()
  for (comp in sort(unique(component))) {
    ii <- which(component == comp)
    candidates <- cfg$k_candidates[cfg$k_candidates <= nrow(unique(as.data.table(x[ii, , drop = FALSE])))]
    local <- list(); wcss <- numeric(length(candidates))
    for (j in seq_along(candidates)) {
      k <- candidates[[j]]
      set.seed(cfg$seed + offset + comp * 100L + k)
      local[[as.character(k)]] <- kmeans(
        x[ii, , drop = FALSE], centers = k, nstart = cfg$kmeans_nstart,
        iter.max = cfg$kmeans_iter_max, algorithm = cfg$kmeans_algorithm
      )
      wcss[[j]] <- local[[as.character(k)]]$tot.withinss
    }
    selected_k <- choose_knee(candidates, wcss)
    fit <- local[[as.character(selected_k)]]
    ordering <- order(fit$centers[, 1L], fit$centers[, 2L])
    relabel <- integer(selected_k); relabel[ordering] <- seq_len(selected_k)
    clusters[ii] <- relabel[fit$cluster]
    fits[[as.character(comp)]] <- list(k = selected_k, centers = fit$centers[ordering, , drop = FALSE])
    rows[[as.character(comp)]] <- data.table(
      nmf_component = comp, k = candidates, wcss = wcss,
      selected = candidates == selected_k
    )
  }
  list(cluster = clusters, fits = fits, diagnostics = rbindlist(rows))
}

fit_kmeans_fixed <- function(component, x, k_by_component, offset = 0L) {
  clusters <- integer(nrow(x)); fits <- list()
  for (comp in sort(unique(component))) {
    ii <- which(component == comp); k <- as.integer(k_by_component[[as.character(comp)]])
    set.seed(cfg$seed + offset + comp * 100L + k)
    fit <- kmeans(
      x[ii, , drop = FALSE], centers = k, nstart = cfg$kmeans_nstart,
      iter.max = cfg$kmeans_iter_max, algorithm = cfg$kmeans_algorithm
    )
    ordering <- order(fit$centers[, 1L], fit$centers[, 2L])
    relabel <- integer(k); relabel[ordering] <- seq_len(k)
    clusters[ii] <- relabel[fit$cluster]
    fits[[as.character(comp)]] <- list(k = k, centers = fit$centers[ordering, , drop = FALSE])
  }
  list(cluster = clusters, fits = fits)
}

project_nmf <- function(rhs, nmf_model) {
  w <- sweep(nmf_model$w, 2L, nmf_model$d, "*")
  h <- RcppML::nnls(crossprod(w), crossprod(w, t(rhs)))
  if (is.null(dim(h))) h <- matrix(h, nrow = ncol(w))
  h
}

assign_by_centers <- function(component, x, fits) {
  local <- integer(nrow(x))
  for (comp in sort(unique(component))) {
    ii <- which(component == comp); centers <- fits[[as.character(comp)]]$centers
    distances <- sapply(seq_len(nrow(centers)), function(j) {
      rowSums((x[ii, , drop = FALSE] - matrix(
        centers[j, ], nrow = length(ii), ncol = ncol(x), byrow = TRUE
      ))^2)
    })
    if (is.null(dim(distances))) distances <- matrix(distances, ncol = 1L)
    local[ii] <- max.col(-distances, ties.method = "first")
  }
  local
}

sample_pool <- function(n) {
  set.seed(cfg$seed + 9000L)
  sample.int(nrow(o$rhs), min(as.integer(n), nrow(o$rhs)))
}

comparison_pool <- sample_pool(cfg$comparison_n)
reference_w <- reference_model$nmf$w

# ------------------------------------------------------------------------- #
# 01. NMF stability across random initializations
# ------------------------------------------------------------------------- #
file01 <- out_file("01_NMF_restart_stability.csv")
if (must_run(file01)) {
  logmsg("01: NMF restart stability")
  A <- as_sparse(t(o$rhs[comparison_pool, , drop = FALSE]))
  rows <- vector("list", cfg$repetitions)
  for (r in seq_len(cfg$repetitions)) {
    fit <- fit_best_nmf(A, rank_selected, 1L, 100000L + r * 100L)
    match <- match_components(reference_w, fit$model$w)
    component <- match$candidate_to_reference[rowmax(t(fit$model$h))]
    rows[[r]] <- data.table(
      restart = r, n = length(comparison_pool), mse = fit$mse,
      iterations = fit$model$iter,
      mean_matched_basis_correlation = match$mean_matched_correlation,
      minimum_matched_basis_correlation = min(match$matched_correlations),
      reference_component_ARI = ari(o$component[comparison_pool], component),
      reference_component_NMI = nmi(o$component[comparison_pool], component)
    )
  }
  fwrite(rbindlist(rows), file01)
  rm(A); gc(FALSE)
}

# ------------------------------------------------------------------------- #
# 02–03. Repeated complete classification on 80% subsamples
# ------------------------------------------------------------------------- #
file02 <- out_file("02_restart_and_subsample_stability.csv")
file03 <- out_file("03_k_selection_stability.csv")
if (must_run(file02) || must_run(file03)) {
  logmsg("02–03: complete-procedure subsampling and k-selection stability")
  rows02 <- list(); rows03 <- list()
  for (r in seq_len(cfg$repetitions)) {
    set.seed(cfg$seed + 200000L + r)
    ii <- sample(comparison_pool, floor(length(comparison_pool) * cfg$subsample_fraction))
    fit <- fit_best_nmf(
      as_sparse(t(o$rhs[ii, , drop = FALSE])), rank_selected,
      cfg$rank_restarts_per_repetition, 200000L + r * 100L
    )
    match <- match_components(reference_w, fit$model$w)
    component <- match$candidate_to_reference[rowmax(t(fit$model$h))]
    km <- fit_kmeans_free(component, o$height_features[ii, , drop = FALSE], 300000L + r * 100L)
    final_key <- paste(component, km$cluster, sep = "__")
    reference_key <- paste(o$component[ii], o$local_cluster[ii], sep = "__")
    selected_k <- km$diagnostics[selected == TRUE][order(nmf_component)]
    rows02[[r]] <- data.table(
      replicate = r, n = length(ii), mse = fit$mse,
      selected_final_class_count = uniqueN(final_key),
      mean_matched_basis_correlation = match$mean_matched_correlation,
      reference_component_ARI = ari(o$component[ii], component),
      reference_component_NMI = nmi(o$component[ii], component),
      final_class_ARI = ari(reference_key, final_key),
      final_class_NMI = nmi(reference_key, final_key)
    )
    rows03[[r]] <- selected_k[, .(
      replicate = r, nmf_component, selected_k = k,
      selected_wcss = wcss,
      reference_k = as.integer(fixed_k[as.character(nmf_component)]),
      matches_reference_k = k == as.integer(fixed_k[as.character(nmf_component)])
    )]
  }
  dt02 <- rbindlist(rows02); dt03 <- rbindlist(rows03)
  dt03[, complete_hierarchy := all(matches_reference_k), by = replicate]
  if (must_run(file02)) fwrite(dt02, file02)
  if (must_run(file03)) fwrite(dt03, file03)
}

# ------------------------------------------------------------------------- #
# Create/reuse conventional year × forest-class balanced training sample.
# This is a streaming priority reservoir, so the full population is never
# collected into memory.
# ------------------------------------------------------------------------- #
conventional_checkpoint <- file.path(checkpoint_dir, "conventional_year_forest_sample.rds")

make_conventional_sample <- function() {
  ds <- arrow::open_dataset(cfg$input_files, format = "parquet")
  mb_col <- intersect(c("mapbiomas_class", "rfld_mapbiomas_class"), names(ds))
  if (!length(mb_col)) stop("No MapBiomas forest-class column found in input files")
  mb_col <- mb_col[[1L]]
  cols <- c("sample_id", "year", mb_col, rh_names)
  absent <- setdiff(cols, names(ds)); if (length(absent)) stop("Missing: ", paste(absent, collapse = ", "))
  reader <- arrow::as_record_batch_reader(ds %>% dplyr::select(dplyr::all_of(cols)))
  strata_n <- 5L * 4L
  quota <- ceiling(cfg$conventional_training_n / strata_n)
  retained <- NULL; batch_id <- 0L
  repeat {
    batch <- reader$read_next_batch(); if (is.null(batch)) break
    batch_id <- batch_id + 1L
    d <- as.data.table(as.data.frame(batch)); setnames(d, mb_col, "mapbiomas_class")
    d <- d[year %in% 2019:2023 & mapbiomas_class %in% c(3L, 4L, 5L, 6L)]
    if (!nrow(d)) next
    complete <- complete.cases(d[, ..rh_names])
    d <- d[complete]; if (!nrow(d)) next
    # Deterministic random priority given cfg$seed and streaming order.
    set.seed(cfg$seed + 400000L + batch_id)
    d[, priority__ := runif(.N)]
    retained <- rbindlist(list(retained, d), use.names = TRUE, fill = TRUE)
    retained <- retained[order(priority__), head(.SD, quota), by = .(year, mapbiomas_class)]
    if (batch_id %% 50L == 0L) logmsg("Conventional reservoir: batch ", batch_id,
                                      "; retained ", format(nrow(retained), big.mark = ","))
  }
  setorder(retained, year, mapbiomas_class, priority__)
  retained <- retained[, head(.SD, quota), by = .(year, mapbiomas_class)]
  if (nrow(retained) > cfg$conventional_training_n) retained <- retained[seq_len(cfg$conventional_training_n)]
  retained[, priority__ := NULL]
  retained
}

# ------------------------------------------------------------------------- #
# 04. Sensitivity to conventional versus landscape-enriched sampling
# ------------------------------------------------------------------------- #
file04 <- out_file("04_sampling_design_sensitivity.csv")
if (must_run(file04)) {
  logmsg("04: sampling-design sensitivity")
  if (file.exists(conventional_checkpoint) && !cfg$overwrite) {
    conventional <- readRDS(conventional_checkpoint)
  } else {
    conventional <- make_conventional_sample()
    saveRDS(conventional, conventional_checkpoint, compress = FALSE)
  }
  conventional_raw <- as.matrix(conventional[, ..rh_names]); storage.mode(conventional_raw) <- "double"
  conventional_span <- conventional_raw[, 99L] - apply(conventional_raw, 1L, min)
  valid_conventional <- is.finite(conventional_span) & conventional_span > 0
  if (!all(valid_conventional)) {
    conventional <- conventional[valid_conventional]
    conventional_raw <- conventional_raw[valid_conventional, , drop = FALSE]
  }
  conventional_rhs <- shape_normalize(conventional_raw)
  conventional_hf <- height_features(conventional_raw)
  fit <- fit_best_nmf(as_sparse(t(conventional_rhs)), rank_selected, 3L, 410000L)
  match <- match_components(reference_w, fit$model$w)
  component_train <- match$candidate_to_reference[rowmax(t(fit$model$h))]
  km <- fit_kmeans_free(component_train, conventional_hf$x, 420000L)
  
  set.seed(cfg$seed + 430000L)
  eval_idx <- sample.int(nrow(o$rhs), min(cfg$evaluation_n, nrow(o$rhs)))
  h_eval <- project_nmf(o$rhs[eval_idx, , drop = FALSE], fit$model)
  component_eval <- match$candidate_to_reference[rowmax(t(h_eval))]
  eval_hf <- height_features(o$raw_rh[eval_idx, , drop = FALSE],
                             conventional_hf$center, conventional_hf$scale)$x
  local_eval <- assign_by_centers(component_eval, eval_hf, km$fits)
  alternative_key <- paste(component_eval, local_eval, sep = "__")
  reference_key <- paste(o$component[eval_idx], o$local_cluster[eval_idx], sep = "__")
  selected <- km$diagnostics[selected == TRUE][order(nmf_component)]
  result04 <- data.table(
    comparison = "landscape_enriched_reference_vs_year_forest_balanced",
    conventional_training_n = nrow(conventional), evaluation_n = length(eval_idx),
    alternative_mse = fit$mse,
    mean_matched_basis_correlation = match$mean_matched_correlation,
    component_ARI = ari(o$component[eval_idx], component_eval),
    component_NMI = nmi(o$component[eval_idx], component_eval),
    final_class_ARI = ari(reference_key, alternative_key),
    final_class_NMI = nmi(reference_key, alternative_key),
    alternative_final_class_count = uniqueN(alternative_key),
    k_by_component = paste0(selected$nmf_component, "=", selected$k, collapse = ";")
  )
  fwrite(result04, file04)
  rm(conventional_raw, conventional_rhs, h_eval); gc(FALSE)
}

# ------------------------------------------------------------------------- #
# 05. Pairwise distinctness of final class-average original RH profiles
# ------------------------------------------------------------------------- #
file05 <- out_file("05_pairwise_class_profile_distinctness.csv")
if (must_run(file05)) {
  logmsg("05: pairwise class-profile distinctness")
  classes <- sort(unique(o$architecture_class))
  means <- lapply(classes, function(cl) colMeans(o$raw_rh[o$architecture_class == cl, , drop = FALSE]))
  names(means) <- as.character(classes)
  pairs <- t(combn(classes, 2L))
  rows <- lapply(seq_len(nrow(pairs)), function(j) {
    a <- pairs[j, 1L]; b <- pairs[j, 2L]; difference <- means[[as.character(a)]] - means[[as.character(b)]]
    data.table(
      class_a = a, class_b = b,
      n_a = sum(o$architecture_class == a), n_b = sum(o$architecture_class == b),
      profile_correlation = cor(means[[as.character(a)]], means[[as.character(b)]]),
      profile_RMSE_m = sqrt(mean(difference^2)),
      profile_mean_absolute_difference_m = mean(abs(difference)),
      maximum_absolute_difference_m = max(abs(difference))
    )
  })
  fwrite(rbindlist(rows), file05)
}

# ------------------------------------------------------------------------- #
# 06. Silhouette widths in combined NMF-contribution + height-feature space
# ------------------------------------------------------------------------- #
file06 <- out_file("06_silhouette_widths.csv")
if (must_run(file06)) {
  logmsg("06: final-class silhouette widths")
  set.seed(cfg$seed + 600000L)
  ii <- sample.int(nrow(o$rhs), min(cfg$silhouette_n, nrow(o$rhs)))
  scores <- t(o$model$nmf$h[, ii, drop = FALSE])
  scores <- scores / pmax(rowSums(scores), .Machine$double.eps)
  features <- cbind(scores, o$height_features[ii, , drop = FALSE])
  features <- scale(features)
  sil <- cluster::silhouette(as.integer(factor(o$architecture_class[ii])), dist(features))
  sil_dt <- data.table(
    sample_index = ii,
    architecture_class = o$architecture_class[ii],
    silhouette_width = sil[, "sil_width"]
  )
  summary_class <- sil_dt[, .(
    level = "class", n = .N, mean_silhouette = mean(silhouette_width),
    median_silhouette = median(silhouette_width),
    sd_silhouette = sd(silhouette_width),
    proportion_negative = mean(silhouette_width < 0)
  ), by = .(architecture_class)]
  summary_all <- sil_dt[, .(
    architecture_class = NA_integer_, level = "overall", n = .N,
    mean_silhouette = mean(silhouette_width), median_silhouette = median(silhouette_width),
    sd_silhouette = sd(silhouette_width), proportion_negative = mean(silhouette_width < 0)
  )]
  fwrite(rbindlist(list(summary_all, summary_class), use.names = TRUE), file06)
}

# ------------------------------------------------------------------------- #
# 07. Blocked 10-km spatial holdout
# ------------------------------------------------------------------------- #
file07 <- out_file("07_spatial_holdout.csv")
if (must_run(file07)) {
  logmsg("07: blocked spatial holdout")
  latitude <- o$meta$latitude; longitude <- o$meta$longitude
  mean_lat <- mean(latitude, na.rm = TRUE) * pi / 180
  x_m <- longitude * 111320 * cos(mean_lat); y_m <- latitude * 110540
  cell <- paste(floor(x_m / cfg$spatial_cell_size_m), floor(y_m / cfg$spatial_cell_size_m), sep = "_")
  occupied <- unique(cell)
  set.seed(cfg$seed + 700000L)
  holdout_cells <- sample(occupied, max(1L, floor(length(occupied) * cfg$spatial_holdout_fraction)))
  train_available <- which(!cell %in% holdout_cells); test_available <- which(cell %in% holdout_cells)
  train_idx <- sample(train_available, min(cfg$spatial_training_n, length(train_available)))
  test_idx <- sample(test_available, min(cfg$spatial_evaluation_n, length(test_available)))
  
  fit <- fit_best_nmf(as_sparse(t(o$rhs[train_idx, , drop = FALSE])), rank_selected, 3L, 710000L)
  match <- match_components(reference_w, fit$model$w)
  component_train <- match$candidate_to_reference[rowmax(t(fit$model$h))]
  hf_train <- height_features(o$raw_rh[train_idx, , drop = FALSE])
  km <- fit_kmeans_fixed(component_train, hf_train$x, fixed_k, 720000L)
  h_test <- project_nmf(o$rhs[test_idx, , drop = FALSE], fit$model)
  component_test <- match$candidate_to_reference[rowmax(t(h_test))]
  hf_test <- height_features(o$raw_rh[test_idx, , drop = FALSE], hf_train$center, hf_train$scale)$x
  local_test <- assign_by_centers(component_test, hf_test, km$fits)
  predicted_key <- paste(component_test, local_test, sep = "__")
  reference_key <- paste(o$component[test_idx], o$local_cluster[test_idx], sep = "__")
  result07 <- data.table(
    spatial_cell_size_m = cfg$spatial_cell_size_m,
    occupied_cells = length(occupied), heldout_cells = length(holdout_cells),
    heldout_cell_fraction = length(holdout_cells) / length(occupied),
    training_n = length(train_idx), evaluation_n = length(test_idx),
    mse = fit$mse,
    mean_matched_basis_correlation = match$mean_matched_correlation,
    component_ARI = ari(o$component[test_idx], component_test),
    component_NMI = nmi(o$component[test_idx], component_test),
    final_class_ARI = ari(reference_key, predicted_key),
    final_class_NMI = nmi(reference_key, predicted_key)
  )
  fwrite(result07, file07)
}

# ------------------------------------------------------------------------- #
# 08. Discovery-sample versus classified-population prevalence audit
# ------------------------------------------------------------------------- #
file08 <- out_file("08_training_vs_population_prevalence.csv")
if (must_run(file08)) {
  logmsg("08: training-versus-population prevalence")
  if (!dir.exists(assignment_dir)) stop("Assignment dataset directory missing: ", assignment_dir)
  ads <- arrow::open_dataset(assignment_dir, format = "parquet")
  population <- as.data.table(
    ads %>%
      dplyr::filter(!is.na(architecture_class)) %>%
      dplyr::group_by(architecture_class) %>%
      dplyr::summarise(population_n = dplyr::n()) %>%
      dplyr::collect()
  )
  training <- data.table(architecture_class = o$architecture_class)[, .(training_n = .N), by = architecture_class]
  prevalence <- merge(training, population, by = "architecture_class", all = TRUE)
  prevalence[is.na(training_n), training_n := 0L]
  prevalence[is.na(population_n), population_n := 0L]
  prevalence[, `:=`(
    training_proportion = training_n / sum(training_n),
    population_proportion = population_n / sum(population_n)
  )]
  prevalence[, `:=`(
    percentage_point_difference = 100 * (training_proportion - population_proportion),
    enrichment_ratio = fifelse(population_proportion > 0,
                               training_proportion / population_proportion, NA_real_)
  )]
  setorder(prevalence, architecture_class)
  fwrite(prevalence, file08)
}

# Compact manifest for the Supporting Material.
manifest <- data.table(
  output = sprintf("%02d", 1:8),
  evaluation = c(
    "NMF random-initialization stability",
    "Complete-procedure 80% subsampling stability",
    "Component-specific k-selection stability",
    "Sampling-design sensitivity",
    "Pairwise final-class profile distinctness",
    "Final-class silhouette separation",
    "Blocked 10-km spatial holdout",
    "Training-versus-population prevalence audit"
  ),
  file = c(
    basename(file01), basename(file02), basename(file03), basename(file04),
    basename(file05), basename(file06), basename(file07), basename(file08)
  ),
  status = ifelse(file.exists(c(file01,file02,file03,file04,file05,file06,file07,file08)),
                  "complete", "missing")
)
fwrite(manifest, out_file("00_robustness_manifest.csv"))
writeLines(capture.output(sessionInfo()), out_file("sessionInfo.txt"))
saveRDS(cfg, out_file("robustness_configuration.rds"))
logmsg("STEP 01c complete: ", cfg$output_dir)
