#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
prepared_path <- if (length(args) >= 1L) args[[1L]] else "/scratch/mr3882/temp/Xingu_GAM/Xingu_GAM_prepared.rds"
output_dir <- if (length(args) >= 2L) args[[2L]] else "/scratch/mr3882/temp/Xingu_GAM/fits_reduced_basis"
task_id <- if (length(args) >= 3L) as.integer(args[[3L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
overwrite <- identical(tolower(Sys.getenv("OVERWRITE", "false")), "true")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

specs <- data.table(
  model_id = c("L", "S", "E", "LS", "LE", "SE", "LSE_lagged", "LSE_recent", "LSE_average"),
  representation = c(rep("lagged", 7L), "recent", "average"),
  use_lulc = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  use_landscape = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  use_environment = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
)
architecture_classes <- 1:7
grid <- rbindlist(list(
  CJ(focal_class = architecture_classes, weighting = "weighted", model_row = seq_len(nrow(specs))),
  CJ(focal_class = architecture_classes, weighting = "unweighted",
     model_row = which(specs$model_id == "LSE_lagged"))
))
setorder(grid, focal_class, weighting, model_row)
if (is.na(task_id) || task_id < 0L || task_id >= nrow(grid)) {
  stop("task_id must be between 0 and ", nrow(grid) - 1L)
}
job <- grid[task_id + 1L]
spec <- specs[job$model_row]

tag <- sprintf("class%d_%s_%s", job$focal_class, job$weighting, spec$model_id)
job_dir <- file.path(output_dir, sprintf("class_%d", job$focal_class), job$weighting)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
model_path <- file.path(job_dir, paste0(spec$model_id, ".rds"))
meta_path <- file.path(job_dir, paste0(spec$model_id, "_metadata.rds"))
metric_path <- file.path(job_dir, paste0(spec$model_id, "_metrics.csv"))
prediction_path <- file.path(job_dir, paste0(spec$model_id, "_validation_predictions.csv"))
if (file.exists(model_path) && file.exists(metric_path) && !overwrite) {
  message("Completed task already exists: ", tag)
  quit(save = "no", status = 0L)
}

message("Loading prepared data for ", tag)
p <- readRDS(prepared_path)
b <- as.data.table(p$base)
eligible <- which(b$complete_full & b$analysis_split %chin% c("model", "validation"))
if (!length(eligible)) stop("No complete full-history observations")
train_flag <- b$analysis_split[eligible] == "model"
validation_flag <- b$analysis_split[eligible] == "validation"

scale_vector <- function(x, train) {
  center <- mean(x[train], na.rm = TRUE)
  spread <- sd(x[train], na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- 1
  list(value = (x - center) / spread, center = center, scale = spread)
}
scale_matrix <- function(x, train) {
  center <- mean(x[train, , drop = FALSE], na.rm = TRUE)
  spread <- sd(as.numeric(x[train, , drop = FALSE]), na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- 1
  list(value = (x - center) / spread, center = center, scale = spread)
}

d <- data.frame(
  response = as.integer(b$architecture_class[eligible] == job$focal_class),
  sample_id = b$sample_id[eligible],
  analysis_split = b$analysis_split[eligible],
  mapbiomas_class_current = droplevels(b$mapbiomas_class_current[eligible]),
  tmf_class_current = droplevels(b$tmf_class_current[eligible]),
  stringsAsFactors = FALSE
)
scalers <- list()

add_scaled_vector <- function(name, values) {
  z <- scale_vector(values[eligible], train_flag)
  d[[name]] <<- z$value
  scalers[[name]] <<- z[c("center", "scale")]
}
add_scaled_matrix <- function(name, values) {
  z <- scale_matrix(values[eligible, , drop = FALSE], train_flag)
  d[[name]] <<- I(z$value)
  scalers[[name]] <<- z[c("center", "scale")]
}

if (spec$use_lulc) {
  lulc_names <- setdiff(names(p$lulc), "other_reference")
  for (nm in lulc_names) {
    x <- p$lulc[[nm]]
    if (spec$representation == "recent") {
      x <- x[, c(1L, 21L, 41L, 61L), drop = FALSE]
    } else if (spec$representation == "average") {
      x <- do.call(cbind, lapply(0:3, function(j) {
        rowMeans(x[, (j * 20L + 1L):(j * 20L + 20L), drop = FALSE])
      }))
    }
    add_scaled_matrix(paste0("lulc_", nm), x)
  }
  if (spec$representation == "lagged") {
    d[["lulc_lag"]] <- I(p$axes$lulc_lag[eligible, , drop = FALSE])
    d[["lulc_radius"]] <- I(p$axes$lulc_radius[eligible, , drop = FALSE])
  } else {
    d[["lulc_radius4"]] <- I(p$axes$radius4[eligible, , drop = FALSE])
  }
}

if (spec$use_landscape) {
  for (nm in names(p$landscape)) {
    x <- p$landscape[[nm]]
    if (spec$representation == "recent") x <- x[, 1L]
    if (spec$representation == "average") x <- rowMeans(x)
    if (is.matrix(x)) add_scaled_matrix(paste0("land_", nm), x)
    else add_scaled_vector(paste0("land_", nm), x)
  }
  if (spec$representation == "lagged") {
    d[["landscape_lag"]] <- I(p$axes$landscape_lag[eligible, , drop = FALSE])
  }
}

if (spec$use_environment) {
  for (nm in p$environment_names) add_scaled_vector(nm, b[[nm]])
}

lulc_terms <- character()
if (spec$use_lulc) {
  vars <- paste0("lulc_", setdiff(names(p$lulc), "other_reference"))
  if (spec$representation == "lagged") {
    lulc_terms <- sprintf(
      "te(%s,lulc_radius,lulc_lag,bs=c('ts','ts','ts'),k=c(3,3,4))", vars
    )
  } else {
    lulc_terms <- sprintf(
      "te(%s,lulc_radius4,bs=c('ts','ts'),k=c(3,3))", vars
    )
  }
  lulc_terms <- c(lulc_terms, "mapbiomas_class_current", "tmf_class_current")
}

land_terms <- character()
if (spec$use_landscape) {
  vars <- paste0("land_", names(p$landscape))
  if (spec$representation == "lagged") {
    land_terms <- sprintf(
      "te(%s,landscape_lag,bs=c('ts','ts'),k=c(4,4))", vars
    )
  } else {
    land_terms <- sprintf("s(%s,bs='ts',k=6)", vars)
  }
}

env_terms <- character()
if (spec$use_environment) {
  env_terms <- c(
    "te(longitude,latitude,bs=c('ts','ts'),k=c(10,10))",
    sprintf("s(%s,bs='ts',k=6)", setdiff(p$environment_names, c("longitude", "latitude")))
  )
}
rhs <- c(lulc_terms, land_terms, env_terms)
formula_text <- paste("response ~", paste(rhs, collapse = " + "))
form <- as.formula(formula_text)

train <- d[train_flag, , drop = FALSE]
validation <- d[validation_flag, , drop = FALSE]
weights_train <- if (job$weighting == "weighted") {
  w <- b$sampling_weight_normalized[eligible][train_flag]
  w / mean(w)
} else rep(1, sum(train_flag))

nthreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
message("Fitting ", tag, "; training n=", nrow(train), "; validation n=", nrow(validation))
fit <- bam(
  form, family = binomial(link = "logit"), data = train,
  weights = weights_train, method = "fREML", discrete = TRUE,
  nthreads = max(1L, nthreads), na.action = na.fail, gc.level = 1L
)

balanced_accuracy <- function(y, predicted) {
  sensitivity <- if (sum(y == 1L) > 0L) mean(predicted[y == 1L] == 1L) else NA_real_
  specificity <- if (sum(y == 0L) > 0L) mean(predicted[y == 0L] == 0L) else NA_real_
  mean(c(sensitivity, specificity), na.rm = TRUE)
}
train_prob <- predict(fit, newdata = train, type = "response")
threshold_grid <- seq(0.001, 0.999, by = 0.002)
scores <- vapply(threshold_grid, function(z) {
  balanced_accuracy(train$response, as.integer(train_prob >= z))
}, numeric(1))
threshold <- threshold_grid[which.max(scores)]

metrics <- function(y, prob, threshold, split) {
  pred <- as.integer(prob >= threshold)
  tp <- sum(y == 1L & pred == 1L); tn <- sum(y == 0L & pred == 0L)
  fp <- sum(y == 0L & pred == 1L); fn <- sum(y == 1L & pred == 0L)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(sensitivity) &&
            (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else NA_real_
  data.table(
    split = split, n = length(y), positives = sum(y == 1L), negatives = sum(y == 0L),
    threshold = threshold, balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    f1 = f1, sensitivity = sensitivity, specificity = specificity,
    precision = precision, tp = tp, tn = tn, fp = fp, fn = fn
  )
}

validation_prob <- predict(fit, newdata = validation, type = "response")
result <- rbindlist(list(
  metrics(train$response, train_prob, threshold, "training"),
  metrics(validation$response, validation_prob, threshold, "validation")
))
s <- summary(fit)
result[, `:=`(
  focal_class = job$focal_class,
  weighting = job$weighting,
  model_id = spec$model_id,
  representation = spec$representation,
  use_lulc = spec$use_lulc,
  use_landscape = spec$use_landscape,
  use_environment = spec$use_environment,
  AIC = AIC(fit),
  explained_deviance = s$dev.expl,
  adjusted_r_squared = s$r.sq,
  total_edf = sum(fit$edf),
  converged = isTRUE(fit$converged)
)]
setcolorder(result, c("focal_class", "weighting", "model_id", "representation"))

saveRDS(fit, model_path, compress = FALSE)
saveRDS(list(
  task_id = task_id, formula = formula_text, scalers = scalers,
  threshold = threshold, eligible_row_indices = eligible,
  training_sample_ids = train$sample_id,
  validation_sample_ids = validation$sample_id,
  mgcv_version = as.character(packageVersion("mgcv")),
  prepared_path = prepared_path
), meta_path, compress = "xz")
fwrite(result, metric_path)
fwrite(data.table(
  sample_id = validation$sample_id,
  observed = validation$response,
  predicted_probability = validation_prob,
  predicted_class = as.integer(validation_prob >= threshold),
  threshold = threshold
), prediction_path)
message("Completed ", tag)
