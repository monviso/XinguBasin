#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

cfg <- list(
  history = "~/RFLD_full_history_all_with_nonforest_area.parquet",
  climate = "~/predictor_sample_terrain_climate.csv",
  output_dir = "~/Xingu_GAM",
  overwrite = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) cfg$history <- args[[1L]]
if (length(args) >= 2L) cfg$climate <- args[[2L]]
if (length(args) >= 3L) cfg$output_dir <- args[[3L]]

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(cfg$output_dir, "Xingu_GAM_prepared.rds")
if (file.exists(out_rds) && !cfg$overwrite) {
  message("Prepared object already exists: ", out_rds)
  quit(save = "no", status = 0L)
}

message("Reading RFLD history")
h <- as.data.table(read_parquet(cfg$history))
message("Reading terrain and climate")
e <- fread(cfg$climate)

required_h <- c(
  "sample_id", "architecture_class", "analysis_split", "lag_year",
  "observation_year", "history_year", "history_forest_present",
  "mapbiomas_class", "tmf_class", "external_edge_distance_m",
  "nearest_other_forest_distance_m", "area_ha", "compactness",
  "nearest_nonforest_patch_area_ha"
)
required_e <- c(
  "sample_id", "longitude", "latitude", "sampling_weight_normalized",
  "elevation_m", "slope_deg", "aspect_deg"
)
stopifnot(all(required_h %in% names(h)), all(required_e %in% names(e)))
if (h[, anyDuplicated(paste(sample_id, lag_year))] > 0L) {
  stop("Duplicated sample_id x lag_year rows")
}
if (h[, .N] != uniqueN(h$sample_id) * 20L) stop("History is not 20 rows per sample")

base <- unique(h[, .(
  sample_id, architecture_class, analysis_split, observation_year
)], by = "sample_id")
setnames(base, "observation_year", "year")
base <- merge(base, e, by = c("sample_id", "architecture_class", "analysis_split", "year"),
              all.x = TRUE, sort = FALSE)
if (nrow(base) != uniqueN(h$sample_id)) stop("Climate merge changed sample count")
base[, sample_index := .I]

# Align history exactly to base order and lag 0,...,19.
index <- CJ(sample_index = base$sample_index, lag_year = 0:19)
index[, sample_id := base$sample_id[sample_index]]
h[, lag_year := as.integer(lag_year)]
h <- merge(index, h, by = c("sample_id", "lag_year"), all.x = TRUE, sort = FALSE)
setorder(h, sample_index, lag_year)
if (nrow(h) != nrow(base) * 20L) stop("Aligned history has wrong row count")

radii <- c(100L, 500L, 1000L, 3000L)
mb_groups <- list(
  forest = c(3L, 4L, 5L, 6L),
  wetland = 11L,
  grassland = 12L,
  pasture = 15L,
  agriculture = 18L,
  silviculture_oil_palm = c(9L, 35L),
  urban = 24L,
  mining = 30L,
  water = 33L
)
all_selected <- unique(unlist(mb_groups, use.names = FALSE))

matrix20 <- function(x) {
  matrix(as.numeric(x), nrow = nrow(base), ncol = 20L, byrow = TRUE)
}

prop_column <- function(radius, code) {
  sprintf("mapbiomas_prop_%dm_class_%d", radius, code)
}

lulc <- list()
for (group_name in names(mb_groups)) {
  mats <- vector("list", length(radii))
  for (j in seq_along(radii)) {
    cols <- vapply(mb_groups[[group_name]], function(z) prop_column(radii[[j]], z), "")
    missing <- setdiff(cols, names(h))
    if (length(missing)) stop("Missing MapBiomas columns: ", paste(missing, collapse = ", "))
    value <- rowSums(as.matrix(h[, ..cols]), na.rm = FALSE)
    mats[[j]] <- matrix20(value)
  }
  # Radius-major columns: r100 lag0:19, r500 lag0:19, etc.
  lulc[[group_name]] <- do.call(cbind, mats)
}

# Retain the residual composition for QA/descriptions, but it is the omitted
# reference and therefore receives no tensor term.
lulc_other <- vector("list", length(radii))
for (j in seq_along(radii)) {
  pattern <- sprintf("^mapbiomas_prop_%dm_class_[0-9]+$", radii[[j]])
  all_cols <- grep(pattern, names(h), value = TRUE)
  selected_cols <- vapply(all_selected, function(z) prop_column(radii[[j]], z), "")
  other_cols <- setdiff(all_cols, selected_cols)
  lulc_other[[j]] <- matrix20(rowSums(as.matrix(h[, ..other_cols]), na.rm = FALSE))
}
lulc$other_reference <- do.call(cbind, lulc_other)

landscape <- list(
  d1 = matrix20(log1p(h$external_edge_distance_m)),
  d2 = matrix20(log1p(h$nearest_other_forest_distance_m)),
  focal_patch_area = matrix20(log1p(h$area_ha)),
  focal_patch_compactness = matrix20(h$compactness),
  nearest_nonforest_patch_area = matrix20(log1p(h$nearest_nonforest_patch_area_ha))
)

lag0 <- h[lag_year == 0L]
setorder(lag0, sample_index)
base[, mapbiomas_class_current := factor(lag0$mapbiomas_class)]
base[, tmf_class_current := factor(lag0$tmf_class)]
base[, northness := cos(aspect_deg * pi / 180)]
base[, eastness := sin(aspect_deg * pi / 180)]

climate_names <- grep("^clim_mean_", names(base), value = TRUE)
environment_names <- c(
  "longitude", "latitude", "elevation_m", "slope_deg", "northness", "eastness",
  climate_names
)

lag_lulc <- matrix(rep(0:19, times = length(radii)), nrow = nrow(base),
                   ncol = 80L, byrow = TRUE)
radius_lulc <- matrix(rep(log(radii), each = 20L), nrow = nrow(base),
                      ncol = 80L, byrow = TRUE)
lag_landscape <- matrix(rep(0:19, each = nrow(base)), nrow = nrow(base),
                        ncol = 20L, byrow = FALSE)
radius4 <- matrix(rep(log(radii), each = nrow(base)), nrow = nrow(base),
                  ncol = 4L, byrow = FALSE)

complete_matrix <- function(x) rowSums(!is.finite(x)) == 0L
complete_full <- Reduce(`&`, c(
  lapply(lulc[names(mb_groups)], complete_matrix),
  lapply(landscape, complete_matrix),
  list(complete.cases(base[, ..environment_names]),
       !is.na(base$mapbiomas_class_current), !is.na(base$tmf_class_current))
))
base[, complete_full := complete_full]

prepared <- list(
  base = base,
  lulc = lulc,
  landscape = landscape,
  axes = list(
    lulc_lag = lag_lulc,
    lulc_radius = radius_lulc,
    landscape_lag = lag_landscape,
    radius4 = radius4
  ),
  environment_names = environment_names,
  mapbiomas_groups = mb_groups,
  mapbiomas_other_is_reference = TRUE,
  history_definition = "lag 0 is observation year; lag 19 is observation year minus 19",
  transformations = c(
    d1 = "log1p", d2 = "log1p", focal_patch_area = "log1p",
    focal_patch_compactness = "identity",
    nearest_nonforest_patch_area = "log1p"
  )
)
saveRDS(prepared, out_rds, compress = "xz")

audit <- base[, .(
  n = .N,
  n_complete_full = sum(complete_full),
  proportion_complete_full = mean(complete_full)
), by = .(analysis_split, architecture_class)]
fwrite(audit, file.path(cfg$output_dir, "prepared_complete_case_audit.csv"))

levels_out <- rbindlist(list(
  data.table(variable = "mapbiomas_class_current",
             level = levels(base$mapbiomas_class_current)),
  data.table(variable = "tmf_class_current",
             level = levels(base$tmf_class_current))
))
fwrite(levels_out, file.path(cfg$output_dir, "prepared_factor_levels.csv"))
message("Prepared: ", out_rds)
