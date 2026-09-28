#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
prepared_path <- if (length(args) >= 1L) args[[1L]] else
  "/scratch/mr3882/temp/Xingu_GAM_v3_random/Xingu_GAM_prepared.rds"
fits_root <- if (length(args) >= 2L) args[[2L]] else
  "/scratch/mr3882/temp/Xingu_GAM_v3_random/fits_reduced_basis"
output_root <- if (length(args) >= 3L) args[[3L]] else
  "/scratch/mr3882/temp/Xingu_GAM_v3_random/figures/LSE_all_effects"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
for (d in c("landscape", "lulc", "environment", "categorical", "spatial", "tables")) {
  dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
}

message("Reading prepared object: ", prepared_path)
p <- readRDS(prepared_path)
b <- as.data.table(p$base)
classes <- sort(unique(as.integer(b$architecture_class)))
classes <- classes[is.finite(classes)]
train <- b$analysis_split == "model" & b$complete_full
if (!any(train)) stop("No complete model-training rows")

class_colours <- c(
  "1" = "#EFC000", "2" = "#E76F00", "3" = "#8E2C62",
  "4" = "#44236D", "5" = "#3FA3A5", "6" = "#8FD694",
  "7" = "#0072B2"
)

find_bundle <- function(cl) {
  weighted <- file.path(fits_root, sprintf("class_%d", cl), "weighted")
  unweighted <- file.path(fits_root, sprintf("class_%d", cl), "unweighted")
  model_path <- file.path(weighted, "LSE_lagged.rds")
  meta_path <- file.path(weighted, "LSE_lagged_metadata.rds")
  source <- "weighted"
  if (!file.exists(model_path) || !file.exists(meta_path)) {
    model_path <- file.path(unweighted, "LSE_lagged.rds")
    meta_path <- file.path(unweighted, "LSE_lagged_metadata.rds")
    source <- "unweighted_fallback"
  }
  if (!file.exists(model_path) || !file.exists(meta_path)) {
    stop("Missing LSE model/metadata for class ", cl)
  }
  list(
    fit = readRDS(model_path), metadata = readRDS(meta_path),
    weighting_source = source, model_path = model_path
  )
}

bundles <- setNames(lapply(classes, find_bundle), classes)

finite_quantile <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

choose_contrast <- function(x) {
  q <- finite_quantile(x, c(.10, .25, .50, .75, .90))
  low <- q[[2L]]; high <- q[[4L]]; rule <- "P75 versus P25"
  if (!is.finite(low) || !is.finite(high) || high <= low) {
    low <- q[[1L]]; high <- q[[5L]]; rule <- "P90 versus P10"
  }
  if (!is.finite(low) || !is.finite(high) || high <= low) {
    xx <- as.numeric(x); xx <- xx[is.finite(xx)]
    low <- if (length(xx)) min(xx) else NA_real_
    high <- if (length(xx)) max(xx) else NA_real_
    rule <- "maximum versus minimum"
  }
  list(low = low, high = high, median = q[[3L]], rule = rule)
}

choose_sparse_lulc_contrast <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  ans <- choose_contrast(x)
  if (is.finite(ans$low) && is.finite(ans$high) && ans$high > ans$low &&
      ans$rule != "maximum versus minimum") return(ans)
  positive <- x[x > 0]
  if (length(positive)) {
    return(list(
      low = 0,
      high = as.numeric(quantile(positive, .75, names = FALSE, type = 8)),
      median = median(x),
      rule = "positive-value P75 versus absence"
    ))
  }
  ans
}

get_smooth <- function(fit, pattern) {
  labels <- vapply(fit$smooth, function(s) s$label, character(1))
  hit <- which(grepl(pattern, labels, fixed = TRUE))
  if (length(hit) != 1L) {
    stop("Expected one smooth matching '", pattern, "'; found: ",
         paste(labels[hit], collapse = " | "))
  }
  fit$smooth[[hit]]
}

smooth_difference <- function(fit, smooth, high_data, low_data) {
  xh <- PredictMat(smooth, high_data)
  xl <- PredictMat(smooth, low_data)
  xd <- xh - xl
  ii <- smooth$first.para:smooth$last.para
  beta <- coef(fit)[ii]
  covariance <- fit$Vp[ii, ii, drop = FALSE]
  estimate <- drop(xd %*% beta)
  se <- sqrt(pmax(0, rowSums((xd %*% covariance) * xd)))
  data.table(
    log_odds_ratio = estimate,
    se = se,
    lower_log_or = estimate - 1.96 * se,
    upper_log_or = estimate + 1.96 * se,
    odds_ratio = exp(estimate),
    lower_or = exp(estimate - 1.96 * se),
    upper_or = exp(estimate + 1.96 * se)
  )
}

smooth_value <- function(fit, smooth, newdata) {
  x <- PredictMat(smooth, newdata)
  ii <- smooth$first.para:smooth$last.para
  beta <- coef(fit)[ii]
  covariance <- fit$Vp[ii, ii, drop = FALSE]
  estimate <- drop(x %*% beta)
  se <- sqrt(pmax(0, rowSums((x %*% covariance) * x)))
  data.table(effect = estimate, se = se,
             lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
}

scale_with <- function(x, scaler) {
  (x - as.numeric(scaler$center)) / as.numeric(scaler$scale)
}

# -------------------------------------------------------------------------
# Smooth relevance table
# -------------------------------------------------------------------------
relevance <- rbindlist(lapply(classes, function(cl) {
  fit <- bundles[[as.character(cl)]]$fit
  tab <- summary(fit)$s.table
  if (is.null(tab) || !nrow(tab)) return(NULL)
  statistic_name <- intersect(c("Chi.sq", "F"), colnames(tab))[[1L]]
  p_name <- grep("p-value", colnames(tab), value = TRUE)[[1L]]
  data.table(
    architecture_class = cl,
    weighting_source = bundles[[as.character(cl)]]$weighting_source,
    term = rownames(tab),
    edf = as.numeric(tab[, "edf"]),
    reference_df = if ("Ref.df" %in% colnames(tab)) as.numeric(tab[, "Ref.df"]) else NA_real_,
    statistic = as.numeric(tab[, statistic_name]),
    p_value = as.numeric(tab[, p_name])
  )
}), fill = TRUE)
relevance[, p_adjusted_within_class := p.adjust(p_value, method = "BH"), by = architecture_class]
relevance[, p_adjusted_global := p.adjust(p_value, method = "BH")]
fwrite(relevance, file.path(output_root, "tables", "LSE_smooth_relevance.csv"))

# -------------------------------------------------------------------------
# Landscape temporal contrasts
# -------------------------------------------------------------------------
land_labels <- c(
  d1 = "Distance from external forest edge",
  d2 = "Distance to another forest patch",
  focal_patch_area = "Focal forest-patch area",
  focal_patch_compactness = "Focal patch compactness",
  nearest_nonforest_patch_area = "Adjacent non-forest patch area"
)

land_contrasts <- rbindlist(lapply(names(land_labels), function(term) {
  z <- p$landscape[[term]][train, , drop = FALSE]
  contrast <- choose_contrast(z)
  rbindlist(lapply(classes, function(cl) {
    bundle <- bundles[[as.character(cl)]]
    fit <- bundle$fit; meta <- bundle$metadata
    key <- paste0("land_", term)
    sm <- get_smooth(fit, paste0("te(", key, ",landscape_lag)"))
    lag <- 0:19
    hi <- data.frame(value = rep(scale_with(contrast$high, meta$scalers[[key]]), 20), lag = lag)
    lo <- data.frame(value = rep(scale_with(contrast$low, meta$scalers[[key]]), 20), lag = lag)
    names(hi) <- names(lo) <- c(key, "landscape_lag")
    ans <- smooth_difference(fit, sm, hi, lo)
    ans[, `:=`(
      predictor = term, predictor_label = land_labels[[term]],
      architecture_class = cl, temporal_lag = lag,
      contrast_low = contrast$low, contrast_high = contrast$high,
      contrast_rule = contrast$rule,
      weighting_source = bundle$weighting_source
    )]
    ans
  }))
}))
fwrite(land_contrasts, file.path(output_root, "tables", "LSE_landscape_temporal_effects.csv"))

p_land <- ggplot(
  land_contrasts,
  aes(temporal_lag, odds_ratio, colour = factor(architecture_class),
      fill = factor(architecture_class))
) +
  geom_hline(yintercept = 1, linewidth = .35, colour = "grey40") +
  geom_ribbon(aes(ymin = lower_or, ymax = upper_or), alpha = .10,
              colour = NA) +
  geom_line(linewidth = .9) +
  facet_wrap(~predictor_label, scales = "free_y", ncol = 2) +
  scale_y_log10() +
  scale_colour_manual(values = class_colours) +
  scale_fill_manual(values = class_colours) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 19)) +
  labs(x = "Temporal lag (years)", y = "Odds ratio: high versus low exposure",
       colour = "Architecture class", fill = "Architecture class") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(output_root, "landscape", "LSE_landscape_all_temporal_effects.png"),
       p_land, width = 11, height = 9, dpi = 350, bg = "white")
ggsave(file.path(output_root, "landscape", "LSE_landscape_all_temporal_effects.pdf"),
       p_land, width = 11, height = 9)

for (term in names(land_labels)) {
  pp <- p_land %+% land_contrasts[predictor == term] + facet_wrap(~predictor_label)
  ggsave(file.path(output_root, "landscape", paste0("LSE_landscape_", term, ".png")),
         pp, width = 8, height = 5.5, dpi = 350, bg = "white")
}

# Full predictor-by-temporal-lag surfaces.  At each lag the effect is expressed
# relative to the training median of that predictor, so exp(effect) is an odds
# ratio rather than an unanchored smooth value.
land_axis_value <- function(term, x) {
  if (term %in% c("d1", "d2")) return(expm1(x) / 1000) # log1p(m) -> km
  if (term %in% c("focal_patch_area", "nearest_nonforest_patch_area")) {
    return(expm1(x) / 100) # log1p(ha) -> km2
  }
  x
}
land_model_value <- function(term, x) {
  if (term %in% c("d1", "d2")) return(log1p(x * 1000)) # km -> log1p(m)
  if (term %in% c("focal_patch_area", "nearest_nonforest_patch_area")) {
    return(log1p(x * 100)) # km2 -> log1p(ha)
  }
  x
}
land_axis_label <- c(
  d1 = "Distance from external forest edge (km)",
  d2 = "Distance to another forest patch (km)",
  focal_patch_area = "Focal forest-patch area (km²)",
  focal_patch_compactness = "Focal patch compactness",
  nearest_nonforest_patch_area = "Adjacent non-forest patch area (km²)"
)

land_surfaces <- rbindlist(lapply(names(land_labels), function(term) {
  observed <- as.numeric(p$landscape[[term]][train, , drop = FALSE])
  bounds <- finite_quantile(observed, c(.02, .50, .98))
  display_bounds <- land_axis_value(term, bounds[c(1L, 3L)])
  display_values <- seq(display_bounds[[1L]], display_bounds[[2L]], length.out = 120L)
  term_grid <- CJ(
    predictor_value = display_values,
    temporal_lag = seq(0, 19, length.out = 100L)
  )
  grid_values <- land_model_value(term, term_grid$predictor_value)
  rbindlist(lapply(classes, function(cl) {
    bundle <- bundles[[as.character(cl)]]
    fit <- bundle$fit; meta <- bundle$metadata
    key <- paste0("land_", term)
    sm <- get_smooth(fit, paste0("te(", key, ",landscape_lag)"))
    hi <- data.frame(
      value = scale_with(grid_values, meta$scalers[[key]]),
      lag = term_grid$temporal_lag
    )
    lo <- data.frame(
      value = rep(scale_with(bounds[[2L]], meta$scalers[[key]]), nrow(term_grid)),
      lag = term_grid$temporal_lag
    )
    names(hi) <- names(lo) <- c(key, "landscape_lag")
    ans <- smooth_difference(fit, sm, hi, lo)
    ans[, `:=`(
      predictor = term, predictor_label = land_labels[[term]],
      predictor_value = term_grid$predictor_value,
      reference_value = land_axis_value(term, bounds[[2L]]),
      temporal_lag = term_grid$temporal_lag,
      architecture_class = cl,
      weighting_source = bundle$weighting_source
    )]
    ans
  }))
}))
fwrite(land_surfaces, file.path(output_root, "tables", "LSE_landscape_lag_surfaces.csv"))

for (term in names(land_labels)) {
  dd <- land_surfaces[predictor == term]
  lim <- max(abs(quantile(dd$log_odds_ratio, c(.01, .99), na.rm = TRUE)), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  subtitle <- sprintf(
    "Fill = log odds ratio relative to median (%.3f); predictor range P02–P98",
    dd$reference_value[[1L]]
  )
  pp <- ggplot(dd, aes(temporal_lag, predictor_value, fill = log_odds_ratio)) +
    geom_raster(interpolate = FALSE) +
    facet_wrap(~paste("Class", architecture_class), nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
      name = "log(OR)"
    ) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 19)) +
    labs(
      title = land_labels[[term]], subtitle = subtitle,
      x = "Temporal lag (years)", y = land_axis_label[[term]]
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(
    file.path(output_root, "landscape", paste0("LSE_landscape_", term, "_lag_heatmaps.png")),
    pp, width = 15, height = 4.5, dpi = 350, bg = "white"
  )
  ggsave(
    file.path(output_root, "landscape", paste0("LSE_landscape_", term, "_lag_heatmaps.pdf")),
    pp, width = 15, height = 4.5
  )
  for (cl in classes) {
    pc <- pp %+% dd[architecture_class == cl] + facet_null()
    ggsave(
      file.path(output_root, "landscape", sprintf("LSE_landscape_%s_class_%d_heatmap.png", term, cl)),
      pc, width = 6.5, height = 5.2, dpi = 350, bg = "white"
    )
  }
}

# -------------------------------------------------------------------------
# LULC proportion x radius x temporal-lag contrasts
# -------------------------------------------------------------------------
radii_grid <- exp(seq(log(100), log(3000), length.out = 60L))
lag_grid <- seq(0, 19, length.out = 80L)
surface_grid <- CJ(radius_m = radii_grid, temporal_lag = lag_grid)
lulc_names <- setdiff(names(p$lulc), "other_reference")

# RFLD MapBiomas/TMF proportion bands are uint8 scaled fractions:
# 0 = 0%, 254 = 100%, 255 = NoData. Aggregated MapBiomas groups can differ
# slightly from 254 because their constituent bands are rounded independently.
# The GAM used standardized values, so division by 254 does not alter fitted
# effects; it is applied here only for human-readable labels and tables.
lulc_scale_factor <- 254
single_band_groups <- intersect(
  c("wetland", "grassland", "pasture", "agriculture", "urban", "mining", "water"),
  lulc_names
)
nodata_audit <- rbindlist(lapply(lulc_names, function(nm) {
  z <- p$lulc[[nm]][train, , drop = FALSE]
  data.table(
    predictor = nm,
    finite_n = sum(is.finite(z)),
    exact_255_n = sum(z == 255, na.rm = TRUE),
    minimum_stored = min(z, na.rm = TRUE),
    maximum_stored = max(z, na.rm = TRUE),
    single_band_group = nm %in% single_band_groups
  )
}))
fwrite(nodata_audit, file.path(output_root, "tables", "LSE_LULC_uint8_audit.csv"))
bad_nodata <- nodata_audit[single_band_group & exact_255_n > 0]
if (nrow(bad_nodata)) {
  warning(
    "Stored value 255 occurs in single-band LULC predictors and may represent leaked NoData: ",
    paste(bad_nodata$predictor, collapse = ", ")
  )
}
message("Using documented RFLD uint8 proportion scale: stored value / 254")

lulc_effects <- rbindlist(lapply(lulc_names, function(term) {
  z <- p$lulc[[term]][train, , drop = FALSE]
  contrast <- choose_sparse_lulc_contrast(z)
  if (!is.finite(contrast$low) || !is.finite(contrast$high) || contrast$high <= contrast$low) {
    warning("Skipping constant LULC term: ", term)
    return(NULL)
  }
  rbindlist(lapply(classes, function(cl) {
    bundle <- bundles[[as.character(cl)]]
    fit <- bundle$fit; meta <- bundle$metadata
    key <- paste0("lulc_", term)
    sm <- get_smooth(fit, paste0("te(", key, ",lulc_radius,lulc_lag)"))
    hi <- data.frame(
      value = rep(scale_with(contrast$high, meta$scalers[[key]]), nrow(surface_grid)),
      radius = log(surface_grid$radius_m), lag = surface_grid$temporal_lag
    )
    lo <- hi
    lo$value <- scale_with(contrast$low, meta$scalers[[key]])
    names(hi) <- names(lo) <- c(key, "lulc_radius", "lulc_lag")
    ans <- smooth_difference(fit, sm, hi, lo)
    ans[, `:=`(
      predictor = term, architecture_class = cl,
      radius_m = surface_grid$radius_m,
      temporal_lag = surface_grid$temporal_lag,
      contrast_low = contrast$low, contrast_high = contrast$high,
      contrast_low_proportion = contrast$low / lulc_scale_factor,
      contrast_high_proportion = contrast$high / lulc_scale_factor,
      lulc_scale_factor = lulc_scale_factor,
      contrast_rule = contrast$rule,
      weighting_source = bundle$weighting_source
    )]
    ans
  }))
}), fill = TRUE)
fwrite(lulc_effects, file.path(output_root, "tables", "LSE_LULC_space_time_effects.csv"))

for (term in unique(lulc_effects$predictor)) {
  dd <- lulc_effects[predictor == term]
  lim <- max(abs(quantile(dd$log_odds_ratio, c(.01, .99), na.rm = TRUE)), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  subtitle <- sprintf(
    "%s: %.4f to %.4f proportion (%.2f%% to %.2f%%); fill = log odds ratio",
    dd$contrast_rule[[1L]],
    dd$contrast_low_proportion[[1L]], dd$contrast_high_proportion[[1L]],
    100 * dd$contrast_low_proportion[[1L]],
    100 * dd$contrast_high_proportion[[1L]]
  )
  pp <- ggplot(dd, aes(radius_m, temporal_lag, fill = log_odds_ratio)) +
    geom_raster(interpolate = FALSE) +
    facet_wrap(~paste("Class", architecture_class), nrow = 1) +
    geom_vline(xintercept = c(100, 500, 1000, 3000),
               linewidth = .18, colour = "white", alpha = .5) +
    scale_x_log10(breaks = c(100, 500, 1000, 3000),
                  labels = c("100", "500", "1,000", "3,000")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                         name = "log(OR)") +
    scale_y_continuous(breaks = c(0, 5, 10, 15, 19)) +
    labs(title = paste("LULC:", gsub("_", " ", term)), subtitle = subtitle,
         x = "Spatial radius (m)", y = "Temporal lag (years)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(file.path(output_root, "lulc", paste0("LSE_LULC_", term, "_heatmaps.png")),
         pp, width = 15, height = 4.2, dpi = 350, bg = "white")
  ggsave(file.path(output_root, "lulc", paste0("LSE_LULC_", term, "_heatmaps.pdf")),
         pp, width = 15, height = 4.2)
  for (cl in classes) {
    pc <- pp %+% dd[architecture_class == cl] + facet_null()
    ggsave(file.path(output_root, "lulc", sprintf("LSE_LULC_%s_class_%d.png", term, cl)),
           pc, width = 6.5, height = 5, dpi = 350, bg = "white")
  }
}

# -------------------------------------------------------------------------
# Non-spatial environmental smooths, relative to the predictor median
# -------------------------------------------------------------------------
env_names <- setdiff(p$environment_names, c("longitude", "latitude"))
env_effects <- rbindlist(lapply(env_names, function(term) {
  observed <- b[[term]][train]
  bounds <- finite_quantile(observed, c(.02, .50, .98))
  grid <- seq(bounds[[1L]], bounds[[3L]], length.out = 150L)
  rbindlist(lapply(classes, function(cl) {
    bundle <- bundles[[as.character(cl)]]
    fit <- bundle$fit; meta <- bundle$metadata
    sm <- get_smooth(fit, paste0("s(", term, ")"))
    hi <- setNames(data.frame(scale_with(grid, meta$scalers[[term]])), term)
    lo <- setNames(data.frame(rep(scale_with(bounds[[2L]], meta$scalers[[term]]), length(grid))), term)
    ans <- smooth_difference(fit, sm, hi, lo)
    ans[, `:=`(
      predictor = term, predictor_value = grid,
      reference_value = bounds[[2L]], architecture_class = cl,
      weighting_source = bundle$weighting_source
    )]
    ans
  }))
}))
fwrite(env_effects, file.path(output_root, "tables", "LSE_environment_effects.csv"))

for (term in env_names) {
  dd <- env_effects[predictor == term]
  pp <- ggplot(dd, aes(predictor_value, odds_ratio,
                       colour = factor(architecture_class),
                       fill = factor(architecture_class))) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = .35) +
    geom_ribbon(aes(ymin = lower_or, ymax = upper_or), alpha = .10,
                colour = NA) +
    geom_line(linewidth = .9) +
    scale_y_log10() +
    scale_colour_manual(values = class_colours) +
    scale_fill_manual(values = class_colours) +
    labs(title = gsub("_", " ", term),
         subtitle = sprintf("Odds ratio relative to median = %.3f", dd$reference_value[[1L]]),
         x = term, y = "Odds ratio", colour = "Class", fill = "Class") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_root, "environment", paste0("LSE_environment_", term, ".png")),
         pp, width = 8, height = 5.5, dpi = 350, bg = "white")
}

# -------------------------------------------------------------------------
# Current categorical-class coefficient tables and plots
# -------------------------------------------------------------------------
categorical <- rbindlist(lapply(classes, function(cl) {
  bundle <- bundles[[as.character(cl)]]
  fit <- bundle$fit
  cf <- coef(fit); vv <- diag(fit$Vp)
  idx <- grep("^(mapbiomas_class_current|tmf_class_current)", names(cf))
  if (!length(idx)) return(NULL)
  data.table(
    architecture_class = cl,
    weighting_source = bundle$weighting_source,
    coefficient = names(cf)[idx],
    log_odds_ratio = as.numeric(cf[idx]),
    se = sqrt(pmax(0, vv[idx])),
    odds_ratio = exp(as.numeric(cf[idx])),
    lower_or = exp(as.numeric(cf[idx]) - 1.96 * sqrt(pmax(0, vv[idx]))),
    upper_or = exp(as.numeric(cf[idx]) + 1.96 * sqrt(pmax(0, vv[idx])))
  )
}), fill = TRUE)
categorical[, dataset := fifelse(grepl("^mapbiomas", coefficient), "MapBiomas", "TMF")]
categorical[, level := sub("^(mapbiomas_class_current|tmf_class_current)", "", coefficient)]
fwrite(categorical, file.path(output_root, "tables", "LSE_current_class_fixed_effects.csv"))

for (dataset_name in unique(categorical$dataset)) {
  dd <- categorical[dataset == dataset_name]
  pp <- ggplot(dd, aes(level, odds_ratio, colour = factor(architecture_class),
                       group = architecture_class)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = .35) +
    geom_errorbar(aes(ymin = lower_or, ymax = upper_or),
                  position = position_dodge(width = .65), width = .15) +
    geom_point(position = position_dodge(width = .65), size = 2) +
    scale_y_log10() + scale_colour_manual(values = class_colours) +
    labs(title = paste(dataset_name, "current-year class effects"),
         x = "Class level (relative to model reference level)", y = "Odds ratio",
         colour = "Architecture class") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(output_root, "categorical", paste0("LSE_current_", dataset_name, ".png")),
         pp, width = 9, height = 5.5, dpi = 350, bg = "white")
}

# -------------------------------------------------------------------------
# Longitude-latitude smooths (diagnostic, not an ecological exposure)
# -------------------------------------------------------------------------
lon_bounds <- finite_quantile(b$longitude[train], c(.02, .98))
lat_bounds <- finite_quantile(b$latitude[train], c(.02, .98))
spatial_grid <- CJ(
  longitude = seq(lon_bounds[[1L]], lon_bounds[[2L]], length.out = 100L),
  latitude = seq(lat_bounds[[1L]], lat_bounds[[2L]], length.out = 100L)
)
spatial_effects <- rbindlist(lapply(classes, function(cl) {
  bundle <- bundles[[as.character(cl)]]
  fit <- bundle$fit; meta <- bundle$metadata
  sm <- get_smooth(fit, "te(longitude,latitude)")
  nd <- data.frame(
    longitude = scale_with(spatial_grid$longitude, meta$scalers$longitude),
    latitude = scale_with(spatial_grid$latitude, meta$scalers$latitude)
  )
  ans <- smooth_value(fit, sm, nd)
  ans[, `:=`(
    longitude = spatial_grid$longitude, latitude = spatial_grid$latitude,
    architecture_class = cl, weighting_source = bundle$weighting_source
  )]
  ans
}))
fwrite(spatial_effects, file.path(output_root, "tables", "LSE_spatial_smooths.csv"))

p_spatial <- ggplot(spatial_effects, aes(longitude, latitude, fill = effect)) +
  geom_raster() + facet_wrap(~paste("Class", architecture_class)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Partial\nlog odds") +
  coord_equal() + labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank())
ggsave(file.path(output_root, "spatial", "LSE_longitude_latitude_smooths.png"),
       p_spatial, width = 12, height = 7, dpi = 350, bg = "white")

# Manifest and contrast audit.
contrast_audit <- rbindlist(list(
  unique(land_contrasts[, .(group = "landscape", predictor,
                            contrast_rule, contrast_low, contrast_high)]),
  unique(lulc_effects[, .(group = "LULC", predictor,
                          contrast_rule, contrast_low, contrast_high,
                          contrast_low_proportion, contrast_high_proportion,
                          lulc_scale_factor)])
), fill = TRUE)
fwrite(contrast_audit, file.path(output_root, "tables", "LSE_effect_contrasts.csv"))
fwrite(data.table(
  architecture_class = classes,
  weighting_source = vapply(bundles, `[[`, character(1), "weighting_source"),
  model_path = vapply(bundles, `[[`, character(1), "model_path")
), file.path(output_root, "tables", "LSE_effect_model_manifest.csv"))

message("All LSE effect products written to: ", output_root)
