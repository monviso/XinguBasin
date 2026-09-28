#!/usr/bin/env Rscript

# STEP 08 — Primary external ALS characterization of the original GEDI
# architecture classes. Conventional ALS metrics are compared among GEDI
# classes; ALS is not forced through the GEDI NMF/k-means model.

options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  input_file = "E:/Xingu_rev/data/step07_ALS_profiles_v3/ALS_footprint_profiles.csv",
  output_dir = "E:/Xingu_rev/data/step08_ALS_characterization",
  minimum_class_n_for_tests = 5L,
  fdr_method = "BH",
  overwrite = FALSE
)

required <- c("data.table", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  descriptive = file.path(cfg$output_dir, "ALS_metrics_by_GEDI_class.csv"),
  omnibus = file.path(cfg$output_dir, "ALS_Kruskal_Wallis_tests.csv"),
  pairwise = file.path(cfg$output_dir, "ALS_pairwise_class_tests.csv"),
  sample = file.path(cfg$output_dir, "ALS_class_sample_sizes.csv"),
  profiles = file.path(cfg$output_dir, "ALS_vegetation_profiles_by_GEDI_class.csv"),
  boxplot = file.path(cfg$output_dir, "SM_ALS_metrics_by_GEDI_class.png"),
  heatmap = file.path(cfg$output_dir, "SM_ALS_metric_signatures.png"),
  profile_plot = file.path(cfg$output_dir, "SM_ALS_vegetation_profiles_by_GEDI_class.png")
)
if (!cfg$overwrite && any(file.exists(unlist(paths)))) {
  stop("Step 08 outputs exist. Set cfg$overwrite=TRUE to rebuild them.")
}

dt <- fread(cfg$input_file, colClasses = list(character = "shot_number"))
class_col <- "architecture_class"
metrics <- c(
  "vegetation_height_p50", "vegetation_height_p75", "vegetation_height_p90",
  "vegetation_height_p95", "vegetation_height_p98", "mean_vegetation_height_m",
  "sd_vegetation_height_m", "canopy_cover_gt2m", "canopy_cover_gt5m",
  "canopy_cover_gt10m", "canopy_cover_gt20m", "canopy_cover_gt30m",
  "gap_fraction_below_2m", "vegetation_return_fraction",
  "foliage_height_diversity"
)
absent <- setdiff(c(class_col, metrics), names(dt))
if (length(absent)) {
  stop("Step 07 must be rerun with the ALS-metric patch. Missing: ",
       paste(absent, collapse = ", "))
}
dt[, architecture_class := as.integer(architecture_class)]
dt <- dt[!is.na(architecture_class)]

labels <- c(
  vegetation_height_p50 = "Vegetation P50 (m)",
  vegetation_height_p75 = "Vegetation P75 (m)",
  vegetation_height_p90 = "Vegetation P90 (m)",
  vegetation_height_p95 = "Vegetation P95 (m)",
  vegetation_height_p98 = "Vegetation P98 (m)",
  mean_vegetation_height_m = "Mean vegetation height (m)",
  sd_vegetation_height_m = "Vertical rugosity (SD, m)",
  canopy_cover_gt2m = "Canopy cover >2 m",
  canopy_cover_gt5m = "Canopy cover >5 m",
  canopy_cover_gt10m = "Canopy cover >10 m",
  canopy_cover_gt20m = "Canopy cover >20 m",
  canopy_cover_gt30m = "Canopy cover >30 m",
  gap_fraction_below_2m = "Gap fraction (<2 m)",
  vegetation_return_fraction = "Vegetation-return fraction",
  foliage_height_diversity = "Foliage-height diversity"
)

sample_sizes <- dt[, .N, by = architecture_class][order(architecture_class)]
sample_sizes[, included_in_inference := N >= cfg$minimum_class_n_for_tests]
fwrite(sample_sizes, paths$sample)
eligible_classes <- sample_sizes[included_in_inference == TRUE, architecture_class]

long <- melt(dt, id.vars = c("shot_number", "architecture_class", "als_dataset",
                             "als_year", "gedi_year", "absolute_year_difference"),
             measure.vars = metrics, variable.name = "metric", value.name = "value")
long <- long[is.finite(value)]
long[, metric_label := unname(labels[metric])]

descriptive <- long[, .(
  n = .N, mean = mean(value), sd = sd(value), median = median(value),
  q25 = as.numeric(quantile(value, .25, names = FALSE)),
  q75 = as.numeric(quantile(value, .75, names = FALSE)),
  minimum = min(value), maximum = max(value)
), by = .(metric, metric_label, architecture_class)]
setorder(descriptive, metric, architecture_class)
fwrite(descriptive, paths$descriptive)

epsilon_squared <- function(H, n, k) {
  if (n <= k) return(NA_real_)
  max(0, as.numeric((H - k + 1) / (n - k)))
}
cliffs_delta <- function(x, y) {
  cmp <- outer(x, y, `-`)
  (sum(cmp > 0) - sum(cmp < 0)) / length(cmp)
}

test_long <- long[architecture_class %in% eligible_classes]
omnibus <- rbindlist(lapply(metrics, function(m) {
  x <- test_long[metric == m]
  k <- uniqueN(x$architecture_class)
  if (k < 2L) return(data.table(metric = m, n = nrow(x), classes_tested = k,
                                statistic = NA_real_, df = NA_real_,
                                p_value = NA_real_, epsilon_squared = NA_real_))
  fit <- kruskal.test(value ~ factor(architecture_class), data = x)
  data.table(metric = m, n = nrow(x), classes_tested = k,
             statistic = as.numeric(fit$statistic), df = as.numeric(fit$parameter),
             p_value = as.numeric(fit$p.value),
             epsilon_squared = epsilon_squared(fit$statistic, nrow(x), k))
}))
omnibus[, metric_label := unname(labels[metric])]
omnibus[, p_fdr := p.adjust(p_value, method = cfg$fdr_method)]
omnibus[, significant_fdr_005 := !is.na(p_fdr) & p_fdr < .05]
setcolorder(omnibus, c("metric", "metric_label", setdiff(names(omnibus), c("metric", "metric_label"))))
fwrite(omnibus, paths$omnibus)

pairs <- combn(eligible_classes, 2L, simplify = FALSE)
pairwise <- rbindlist(lapply(metrics, function(m) {
  rbindlist(lapply(pairs, function(pair) {
    x <- test_long[metric == m & architecture_class == pair[1L], value]
    y <- test_long[metric == m & architecture_class == pair[2L], value]
    test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    data.table(metric = m, class_1 = pair[1L], class_2 = pair[2L],
               n_1 = length(x), n_2 = length(y),
               median_1 = median(x), median_2 = median(y),
               median_difference = median(x) - median(y),
               cliffs_delta = cliffs_delta(x, y), p_value = test$p.value)
  }))
}))
pairwise[, metric_label := unname(labels[metric])]
pairwise[, p_fdr := p.adjust(p_value, method = cfg$fdr_method), by = metric]
pairwise[, significant_fdr_005 := p_fdr < .05]
setorder(pairwise, metric, class_1, class_2)
fwrite(pairwise, paths$pairwise)

selected_metrics <- c("vegetation_height_p95", "sd_vegetation_height_m",
                      "canopy_cover_gt10m", "gap_fraction_below_2m",
                      "foliage_height_diversity")
plot_dt <- long[metric %in% selected_metrics]
p_box <- ggplot(plot_dt, aes(factor(architecture_class), value,
                             fill = factor(architecture_class))) +
  geom_boxplot(width = .72, outlier.alpha = .25) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 3) +
  scale_fill_viridis_d(option = "D", end = .9) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(title = "Independent ALS structure among GEDI architecture classes",
       subtitle = "Footprint-level ALS metrics; boxes show median and interquartile range",
       x = "Original GEDI architecture class", y = NULL)
p_box
ggsave(paths$boxplot, p_box, width = 11, height = 7.3, dpi = 320)

heat <- descriptive[, .(metric, metric_label, architecture_class, median)]
heat[, standardized_median := {
  s <- sd(median)
  if (!is.finite(s) || s == 0) rep(0, .N) else (median - mean(median)) / s
}, by = metric]
p_heat <- ggplot(heat, aes(factor(architecture_class), metric_label,
                           fill = standardized_median)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", standardized_median)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8)) +
  labs(title = "ALS structural signatures of GEDI architecture classes",
       subtitle = "Standardized class medians; standardization is within each metric",
       x = "Original GEDI architecture class", y = NULL,
       fill = "Standardized\nmedian")
p_heat
ggsave(paths$heatmap, p_heat, width = 9.2, height = 7.8, dpi = 320)

rh_names <- sprintf("als_rh%02d", 0:98)
if (all(rh_names %in% names(dt))) {
  rh_long <- melt(dt, id.vars = c("shot_number", "architecture_class"),
                  measure.vars = rh_names, variable.name = "percentile",
                  value.name = "height_m")
  rh_long[, percentile := as.integer(sub("als_rh", "", percentile))]
  profile_summary <- rh_long[, .(
    n = uniqueN(shot_number), median_height_m = median(height_m),
    q25_height_m = as.numeric(quantile(height_m, .25, names = FALSE)),
    q75_height_m = as.numeric(quantile(height_m, .75, names = FALSE))
  ), by = .(architecture_class, percentile)]
  fwrite(profile_summary, paths$profiles)
  p_profile <- ggplot(profile_summary,
                      aes(percentile, median_height_m, colour = factor(architecture_class),
                          fill = factor(architecture_class))) +
    geom_ribbon(aes(ymin = q25_height_m, ymax = q75_height_m),
                alpha = .10, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = .8) +
    scale_colour_viridis_d(option = "D", end = .9) +
    scale_fill_viridis_d(option = "D", end = .9) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = "ALS vegetation-height profiles by original GEDI class",
         subtitle = "Vegetation returns only; median and interquartile range",
         x = "ALS vegetation-height percentile", y = "Height above ground (m)",
         colour = "GEDI class")
  ggsave(paths$profile_plot, p_profile, width = 9, height = 6.5, dpi = 320)
}

message("STEP 08 COMPLETE: ", cfg$output_dir)
message("Classes included in inference: ", paste(eligible_classes, collapse = ", "))
excluded <- sample_sizes[included_in_inference == FALSE, architecture_class]
if (length(excluded)) message("Descriptive only because n<", cfg$minimum_class_n_for_tests,
                              ": class ", paste(excluded, collapse = ", "))
