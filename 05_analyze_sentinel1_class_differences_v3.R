#!/usr/bin/env Rscript

# STEP 05 — Summarize OPERA RTC-S1 time series and test whether the selected
# GEDI architecture classes have different microwave signatures.

options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  observation_csv =
    "~/data/step04_sentinel1_opera_v3/sentinel1_observations.csv",
  output_dir = "~/data/step05_sentinel1_analysis_v3",
  minimum_acquisitions = 6L,
  dry_months = 6:9,
  wet_months = c(1:3, 11:12),
  p_adjust_method = "BH"
)

required <- c("data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(cfg$observation_csv)) stop("Not found: ", cfg$observation_csv)
dt <- fread(cfg$observation_csv)
needed <- c("sample_id", "year", "architecture_class", "spatial_block_5km",
            "month", "VV_dB", "VH_dB", "VH_minus_VV_dB", "RVI")
absent <- setdiff(needed, names(dt))
if (length(absent)) stop("Observation file is missing: ", paste(absent, collapse = ", "))

annual <- dt[, {
  valid <- is.finite(VV_dB) & is.finite(VH_dB)
  vv_dry <- VV_dB[valid & month %in% cfg$dry_months]
  vh_dry <- VH_dB[valid & month %in% cfg$dry_months]
  vv_wet <- VV_dB[valid & month %in% cfg$wet_months]
  vh_wet <- VH_dB[valid & month %in% cfg$wet_months]
  .(
    acquisitions = sum(valid),
    VV_median_dB = median(VV_dB[valid]),
    VH_median_dB = median(VH_dB[valid]),
    ratio_median_dB = median(VH_minus_VV_dB[valid]),
    RVI_median = median(RVI[valid]),
    VV_IQR_dB = IQR(VV_dB[valid]),
    VH_IQR_dB = IQR(VH_dB[valid]),
    VV_amplitude_dB = diff(quantile(VV_dB[valid], c(.1, .9))),
    VH_amplitude_dB = diff(quantile(VH_dB[valid], c(.1, .9))),
    VV_wet_minus_dry_dB = if (length(vv_wet) && length(vv_dry))
      median(vv_wet) - median(vv_dry) else NA_real_,
    VH_wet_minus_dry_dB = if (length(vh_wet) && length(vh_dry))
      median(vh_wet) - median(vh_dry) else NA_real_
  )
}, by = .(sample_id, year, architecture_class, spatial_block_5km, burst_id)]

metrics <- setdiff(names(annual),
                   c("sample_id", "year", "architecture_class", "spatial_block_5km",
                     "burst_id", "acquisitions"))
for (metric in metrics) {
  annual[acquisitions < cfg$minimum_acquisitions | !is.finite(get(metric)),
         (metric) := NA_real_]
}
fwrite(annual, file.path(cfg$output_dir, "annual_sentinel1_metrics.csv"))

qa <- annual[, .(
  sampled_points = .N,
  valid_points = sum(acquisitions >= cfg$minimum_acquisitions),
  valid_proportion = mean(acquisitions >= cfg$minimum_acquisitions),
  median_acquisitions = median(acquisitions)
), by = .(year, architecture_class)]
fwrite(qa, file.path(cfg$output_dir, "sentinel1_extraction_QA.csv"))

long <- melt(annual,
             id.vars = c("sample_id", "year", "architecture_class",
                         "spatial_block_5km", "burst_id"),
             measure.vars = metrics, variable.name = "metric", value.name = "value",
             variable.factor = FALSE)[is.finite(value)]
# Statistical tests compare classes after removing the common median of each
# acquisition year. Raw values remain in summaries and plots.
long[, value_year_centered := value - median(value), by = .(metric, year)]

summary_table <- long[, .(
  n = .N, mean = mean(value), sd = sd(value), median = median(value),
  q25 = quantile(value, .25), q75 = quantile(value, .75)
), by = .(metric, architecture_class)]
fwrite(summary_table, file.path(cfg$output_dir, "sentinel1_class_summary.csv"))

global_tests <- long[, {
  test <- kruskal.test(value_year_centered ~ factor(architecture_class))
  k <- uniqueN(architecture_class)
  epsilon2 <- max(0, (as.numeric(test$statistic) - k + 1) / (.N - k))
  .(n = .N, statistic = as.numeric(test$statistic),
    df = as.numeric(test$parameter), p_value = test$p.value,
    epsilon_squared = epsilon2)
}, by = metric]
global_tests[, p_adjusted_BH := p.adjust(p_value, cfg$p_adjust_method)]
global_tests[, effect_magnitude := fifelse(epsilon_squared < .01, "negligible",
                                           fifelse(epsilon_squared < .08, "small",
                                                   fifelse(epsilon_squared < .26, "moderate", "large")))]
fwrite(global_tests, file.path(cfg$output_dir, "sentinel1_global_tests.csv"))

cliffs_delta <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u / (n1 * n2) - 1
}
pairs <- combn(sort(unique(long$architecture_class)), 2L, simplify = FALSE)
rows <- list(); row_id <- 0L
for (metric_name in metrics) {
  md <- long[metric == metric_name]
  for (pair in pairs) {
    x <- md[architecture_class == pair[1L], value_year_centered]
    y <- md[architecture_class == pair[2L], value_year_centered]
    x_raw <- md[architecture_class == pair[1L], value]
    y_raw <- md[architecture_class == pair[2L], value]
    if (length(x) < 2L || length(y) < 2L) next
    test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    row_id <- row_id + 1L
    rows[[row_id]] <- data.table(
      metric = metric_name, class_1 = pair[1L], class_2 = pair[2L],
      n_1 = length(x), n_2 = length(y),
      median_1_raw = median(x_raw), median_2_raw = median(y_raw),
      median_difference_raw_1_minus_2 = median(x_raw) - median(y_raw),
      cliffs_delta_1_vs_2 = cliffs_delta(x, y),
      statistic = as.numeric(test$statistic), p_value = test$p.value)
  }
}
pairwise <- rbindlist(rows)
pairwise[, p_adjusted_within_metric := p.adjust(p_value, cfg$p_adjust_method),
         by = metric]
pairwise[, p_adjusted_all_tests := p.adjust(p_value, cfg$p_adjust_method)]
pairwise[, effect_magnitude := fifelse(abs(cliffs_delta_1_vs_2) < .147,
                                       "negligible", fifelse(abs(cliffs_delta_1_vs_2) < .33, "small",
                                                             fifelse(abs(cliffs_delta_1_vs_2) < .474, "moderate", "large")))]
fwrite(pairwise, file.path(cfg$output_dir, "sentinel1_pairwise_tests.csv"))

heat <- copy(summary_table)
heat[, standardized_median := {
  s <- sd(median)
  if (!is.finite(s) || s == 0) rep(0, .N) else (median - mean(median)) / s
}, by = metric]
p_heat <- ggplot(heat,
                 aes(factor(architecture_class), factor(metric, levels = rev(metrics)),
                     fill = standardized_median)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", standardized_median)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  theme_bw() +
  labs(x = "Vertical architecture class", y = NULL,
       fill = "Standardized\nmedian",
       title = "Sentinel-1 microwave signatures of GEDI classes")
ggsave(file.path(cfg$output_dir, "SM_sentinel1_class_heatmap.png"),
       p_heat, width = 9, height = 8, dpi = 300)

p_box <- ggplot(long,
                aes(factor(architecture_class), value, fill = factor(architecture_class))) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~metric, scales = "free_y") +
  guides(fill = "none") + theme_bw() +
  labs(x = "Vertical architecture class", y = "Annual SAR metric",
       title = "Sentinel-1 distributions among GEDI classes")
ggsave(file.path(cfg$output_dir, "SM_sentinel1_class_boxplots.png"),
       p_box, width = 13, height = 9, dpi = 300)

message("STEP 05 COMPLETE: ", cfg$output_dir)

