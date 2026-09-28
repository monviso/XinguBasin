#!/usr/bin/env Rscript

# STEP 03 — Build annual HLS indices and test radiometric class differences.
options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  sample_csv = "E:/Xingu_rev/data/step02_radiometry_v3/gedi_radiometric_sample.csv",
  appeears_download_dir =
    "E:/Xingu_rev/data/step02_radiometry_v3/appeears_downloads",
  output_dir = "E:/Xingu_rev/data/step03_radiometric_analysis_v3",
  minimum_clear_observations = 3L,
  p_adjust_method = "BH"
)
required <- c("data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

sample_dt <- fread(cfg$sample_csv)
sample_dt[, sample_id := trimws(as.character(sample_id))]
result_files <- list.files(cfg$appeears_download_dir,
                           pattern = "results\\.csv$", recursive = TRUE, full.names = TRUE,
                           ignore.case = TRUE)
if (!length(result_files)) stop("No AppEEARS results CSV files found")
message("Reading ", length(result_files), " AppEEARS files")
raw <- rbindlist(lapply(result_files, fread), use.names = TRUE, fill = TRUE)

find_column <- function(patterns) {
  for (pattern in patterns) {
    hit <- names(raw)[grepl(pattern, names(raw), ignore.case = TRUE)]
    if (length(hit)) return(hit[1L])
  }
  stop("Could not find column matching: ", paste(patterns, collapse = " | "))
}
# AppEEARS schemas vary across products/package versions. Select the identifier
# column empirically rather than assuming that a column named ID is `subtask`.
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  # Also recover IDs if AppEEARS embeds the subtask in a longer label.
  has_embedded <- grepl("XG_[0-9]{6}", x, perl = TRUE)
  x[has_embedded] <- sub(".*?(XG_[0-9]{6}).*", "\\1",
                         x[has_embedded], perl = TRUE)
  x
}
candidate_id_cols <- unique(c(
  names(raw)[grepl("(^ID$|subtask|category|sample_id)", names(raw),
                   ignore.case = TRUE)],
  names(raw)[vapply(raw, function(x) is.character(x) || is.factor(x), logical(1))]
))
if (!length(candidate_id_cols)) stop("No plausible identifier columns found")
id_scores <- vapply(candidate_id_cols, function(column) {
  sum(unique(normalize_sample_id(raw[[column]])) %chin% sample_dt$sample_id)
}, integer(1))
id_col <- candidate_id_cols[which.max(id_scores)]
if (max(id_scores) == 0L) {
  fwrite(data.table(column = candidate_id_cols, overlap = id_scores),
         file.path(cfg$output_dir, "appeears_identifier_column_audit.csv"))
  stop("No AppEEARS identifier column overlaps sample_id. See ",
       file.path(cfg$output_dir, "appeears_identifier_column_audit.csv"))
}
message("Using AppEEARS identifier column: ", id_col,
        "; matched unique IDs=", max(id_scores))
date_col <- find_column(c("^Date$", "date"))
band_names <- c("B02", "B03", "B04", "B05", "B8A", "B11", "B12", "Fmask")
band_cols <- setNames(vapply(band_names, function(b) {
  find_column(c(paste0("(^|[_.])", b, "$"), paste0("HLSS30.*", b)))
}, character(1)), band_names)

dt <- data.table(sample_id = normalize_sample_id(raw[[id_col]]),
                 date = as.Date(raw[[date_col]]))
for (band in band_names) {
  dt[, (band) := suppressWarnings(as.numeric(raw[[band_cols[[band]]]]))]
}
dt <- merge(dt, sample_dt, by = "sample_id", all.x = TRUE)
if (anyNA(dt$architecture_class)) {
  unmatched <- unique(dt[is.na(architecture_class), .(sample_id)])
  fwrite(unmatched, file.path(cfg$output_dir, "unmatched_appeears_ids.csv"))
  matched_n <- dt[!is.na(architecture_class), .N]
  if (matched_n == 0L) {
    stop("No AppEEARS rows matched the sample. See unmatched_appeears_ids.csv")
  }
  warning(nrow(unmatched), " unmatched AppEEARS identifier(s) were excluded; ",
          "see unmatched_appeears_ids.csv")
  dt <- dt[!is.na(architecture_class)]
}

reflectance <- c("B02", "B03", "B04", "B05", "B8A", "B11", "B12")
for (band in reflectance) {
  dt[!is.finite(get(band)) | get(band) <= -9000, (band) := NA_real_]
  if (median(abs(dt[[band]]), na.rm = TRUE) > 2) {
    dt[, (band) := get(band) * 0.0001]
  }
}
dt[!is.finite(Fmask) | Fmask < 0, Fmask := NA_real_]
# HLS Fmask bits 0-3: cirrus, cloud, adjacent cloud, and cloud shadow.
dt[, clear := !is.na(Fmask) & bitwAnd(as.integer(Fmask), 15L) == 0L]
dt[clear == FALSE, (reflectance) := NA_real_]

safe_ratio <- function(a, b) {
  out <- a / b
  out[!is.finite(out) | abs(b) < 1e-8] <- NA_real_
  out
}
dt[, `:=`(
  NDVI = safe_ratio(B8A - B04, B8A + B04),
  EVI = safe_ratio(2.5 * (B8A - B04), B8A + 6 * B04 - 7.5 * B02 + 1),
  NDMI = safe_ratio(B8A - B11, B8A + B11),
  NBR = safe_ratio(B8A - B12, B8A + B12),
  NDRE = safe_ratio(B8A - B05, B8A + B05)
)]
index_metrics <- c("NDVI", "EVI", "NDMI", "NBR", "NDRE")
band_metrics <- reflectance
metrics <- c(index_metrics, band_metrics)

annual <- dt[, c(
  list(clear_observations = sum(clear & is.finite(NDVI), na.rm = TRUE)),
  lapply(.SD, median, na.rm = TRUE)
), by = .(sample_id, year, architecture_class, spatial_block_5km),
.SDcols = metrics]
for (metric in metrics) {
  annual[clear_observations < cfg$minimum_clear_observations |
           !is.finite(get(metric)), (metric) := NA_real_]
}
fwrite(annual, file.path(cfg$output_dir, "annual_HLS_metrics.csv"))

qa <- annual[, .(
  n = .N,
  valid_n = sum(clear_observations >= cfg$minimum_clear_observations),
  valid_proportion = mean(clear_observations >= cfg$minimum_clear_observations),
  median_clear_observations = median(clear_observations)
), by = .(year, architecture_class)]
fwrite(qa, file.path(cfg$output_dir, "radiometric_extraction_QA.csv"))

long <- melt(annual,
             id.vars = c("sample_id", "year", "architecture_class", "spatial_block_5km"),
             measure.vars = metrics, variable.name = "metric", value.name = "value",
             variable.factor = FALSE)[is.finite(value)]
long[, metric_type := fifelse(metric %chin% index_metrics, "Spectral index",
                              "Reflectance band")]
# Remove common annual differences before inference.
long[, value_year_centered := value - median(value), by = .(metric, year)]

summary_table <- long[, .(
  n = .N, mean = mean(value), sd = sd(value), median = median(value),
  q25 = quantile(value, .25), q75 = quantile(value, .75)
), by = .(metric_type, metric, architecture_class)]
fwrite(summary_table, file.path(cfg$output_dir, "radiometric_class_summary.csv"))

global_tests <- long[, {
  test <- kruskal.test(value_year_centered ~ factor(architecture_class))
  k <- uniqueN(architecture_class)
  epsilon2 <- max(0, (as.numeric(test$statistic) - k + 1) / (.N - k))
  .(n = .N, statistic = as.numeric(test$statistic), df = as.numeric(test$parameter),
    p_value = test$p.value, epsilon_squared = epsilon2)
}, by = .(metric_type, metric)]
global_tests[, p_adjusted_BH := p.adjust(p_value, method = cfg$p_adjust_method)]
global_tests[, effect_magnitude := fifelse(epsilon_squared < .01, "negligible",
                                           fifelse(epsilon_squared < .08, "small",
                                                   fifelse(epsilon_squared < .26, "moderate", "large")))]
fwrite(global_tests, file.path(cfg$output_dir, "radiometric_global_tests.csv"))

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
      metric_type = if (metric_name %chin% index_metrics) "Spectral index"
      else "Reflectance band",
      metric = metric_name, class_1 = pair[1L], class_2 = pair[2L],
      n_1 = length(x), n_2 = length(y),
      median_1_raw = median(x_raw), median_2_raw = median(y_raw),
      median_difference_raw_1_minus_2 = median(x_raw) - median(y_raw),
      cliffs_delta_1_vs_2 = cliffs_delta(x, y),
      statistic = as.numeric(test$statistic), p_value = test$p.value)
  }
}
pairwise_tests <- rbindlist(rows)
pairwise_tests[, p_adjusted_within_metric := p.adjust(p_value, cfg$p_adjust_method),
               by = metric]
pairwise_tests[, p_adjusted_all_tests := p.adjust(p_value, cfg$p_adjust_method)]
pairwise_tests[, effect_magnitude := fifelse(abs(cliffs_delta_1_vs_2) < .147,
                                             "negligible", fifelse(abs(cliffs_delta_1_vs_2) < .33, "small",
                                                                   fifelse(abs(cliffs_delta_1_vs_2) < .474, "moderate", "large")))]
fwrite(pairwise_tests, file.path(cfg$output_dir, "radiometric_pairwise_tests.csv"))

heat <- copy(summary_table[metric_type == "Spectral index"])
heat[, standardized_median := {
  s <- sd(median)
  if (!is.finite(s) || s == 0) rep(0, .N) else (median - mean(median)) / s
}, by = metric]
p_heat <- ggplot(heat,
                 aes(factor(architecture_class), factor(metric, levels = rev(index_metrics)),
                     fill = standardized_median)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", standardized_median)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() + theme_bw() +
  labs(x = "Vertical architecture class", y = NULL,
       fill = "Standardized\nmedian",
       title = "HLS radiometric signatures of GEDI classes")
p_heat
ggsave(file.path(cfg$output_dir, "SM_radiometric_class_heatmap.png"),
       p_heat, width = 8, height = 6, dpi = 300)

p_box <- ggplot(long[metric_type == "Spectral index"],
                aes(factor(architecture_class), value,
                    fill = factor(architecture_class))) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~metric, scales = "free_y") +
  guides(fill = "none") + theme_bw() +
  labs(x = "Vertical architecture class", y = "Annual median index",
       title = "HLS radiometric distributions among GEDI classes")
p_box
ggsave(file.path(cfg$output_dir, "SM_radiometric_class_boxplots.png"),
       p_box, width = 11, height = 7, dpi = 300)

# Raw reflectance bands can retain brightness and shadow differences that are
# intentionally suppressed by normalized indices. Keep them in separate plots
# because their interpretation and scales differ from vegetation indices.
band_heat <- copy(summary_table[metric_type == "Reflectance band"])
band_heat[, standardized_median := {
  s <- sd(median)
  if (!is.finite(s) || s == 0) rep(0, .N) else (median - mean(median)) / s
}, by = metric]

p_band_heat <- ggplot(band_heat,
                      aes(factor(architecture_class), factor(metric, levels = rev(band_metrics)),
                          fill = standardized_median)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", standardized_median)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() + theme_bw() +
  labs(x = "Vertical architecture class", y = NULL,
       fill = "Standardized\nmedian",
       title = "HLS reflectance signatures of GEDI classes")
p_band_heat
ggsave(file.path(cfg$output_dir, "SM_radiometric_band_heatmap.png"),
       p_band_heat, width = 8, height = 7, dpi = 300)

p_band_box <- ggplot(long[metric_type == "Reflectance band"],
                     aes(factor(architecture_class), value, fill = factor(architecture_class))) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~factor(metric, levels = band_metrics), scales = "free_y") +
  guides(fill = "none") + theme_bw() +
  labs(x = "Vertical architecture class", y = "Annual median reflectance",
       title = "HLS reflectance distributions among GEDI classes")
p_band_box
ggsave(file.path(cfg$output_dir, "SM_radiometric_band_boxplots.png"),
       p_band_box, width = 12, height = 7, dpi = 300)

wavelengths_nm <- c(B02 = 490, B03 = 560, B04 = 665, B05 = 705,
                    B8A = 865, B11 = 1610, B12 = 2190)
spectral_profile <- summary_table[
  metric_type == "Reflectance band",
  .(architecture_class, band = metric, median_reflectance = median)]
spectral_profile[, wavelength_nm := wavelengths_nm[band]]
setorder(spectral_profile, architecture_class, wavelength_nm)
fwrite(spectral_profile,
       file.path(cfg$output_dir, "radiometric_class_spectral_profiles.csv"))
p_spectral <- ggplot(spectral_profile,
                     aes(wavelength_nm, median_reflectance,
                         colour = factor(architecture_class), group = architecture_class)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) + theme_bw() +
  scale_x_continuous(breaks = wavelengths_nm,
                     labels = paste0(names(wavelengths_nm), "\n",
                                     wavelengths_nm)) +
  labs(x = "HLS band and approximate centre wavelength (nm)",
       y = "Median annual surface reflectance", colour = "Class",
       title = "Median HLS spectral profiles of GEDI classes")
p_spectral
ggsave(file.path(cfg$output_dir, "SM_radiometric_class_spectral_profiles.png"),
       p_spectral, width = 10, height = 6, dpi = 300)
message("STEP 03 COMPLETE: ", cfg$output_dir)

