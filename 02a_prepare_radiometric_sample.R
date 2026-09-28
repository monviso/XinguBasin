#!/usr/bin/env Rscript

# STEP 02A — Prepare a spatially thinned GEDI sample for Sentinel extraction.
#
# The completed Step 01 assignment dataset is partitioned by year and class.
# This script selects at most one footprint per 5-km grid cell within each
# year x architecture-class stratum, then retains up to 1,000 observations per
# stratum. The same table is used by the direct AppEEARS and OPERA workflows.

options(stringsAsFactors = FALSE, warn = 1)

cfg <- list(
  project_dir = "E:/Xingu_rev",
  assignment_dataset =
    "E:/Xingu_rev/data/step01b_final_classification/assignments_dataset",
  output_dir = "E:/Xingu_rev/data/step02_radiometry_v3",
  # AppEEARS returns complete HLS time series for every point; 100 points per
  # class x year (up to 3,000 total) is ample and keeps the workload practical.
  sample_per_class_year = 100L,
  thinning_grid_m = 5000,
  equal_area_crs = 6933,
  seed = 42L,
  overwrite = FALSE
)

required <- c("arrow", "data.table", "dplyr", "sf")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing packages: ", paste(missing, collapse = ", "))
}
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(sf)
})

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
csv_out <- file.path(cfg$output_dir, "gedi_radiometric_sample.csv")
gpkg_out <- file.path(cfg$output_dir, "gedi_radiometric_sample.gpkg")
summary_out <- file.path(cfg$output_dir, "gedi_radiometric_sample_summary.csv")

if (file.exists(csv_out) && !cfg$overwrite) {
  stop("Output already exists: ", csv_out,
       "\nSet cfg$overwrite=TRUE only to intentionally rebuild it.")
}

set.seed(cfg$seed)
ds <- arrow::open_dataset(cfg$assignment_dataset, format = "parquet")
needed <- c("sample_id", "shot_number", "latitude", "longitude", "year",
            "architecture_class", "rh98_original")
optional <- c("mapbiomas_class", "rfld_mapbiomas_class",
              "solar_elevation", "sensitivity")
absent <- setdiff(needed, names(ds))
if (length(absent)) stop("Assignment dataset is missing: ", paste(absent, collapse = ", "))
columns_to_read <- c(needed, intersect(optional, names(ds)))

strata <- ds %>%
  distinct(year, architecture_class) %>%
  collect() %>%
  as.data.table()
setorder(strata, year, architecture_class)

samples <- vector("list", nrow(strata))
for (i in seq_len(nrow(strata))) {
  yr <- strata$year[i]
  cl <- strata$architecture_class[i]
  message("Sampling year=", yr, "; class=", cl)
  
  dt <- ds %>%
    filter(year == yr, architecture_class == cl) %>%
    select(all_of(columns_to_read)) %>%
    collect() %>%
    as.data.table()
  if (!nrow(dt)) next
  
  pts <- st_as_sf(dt, coords = c("longitude", "latitude"),
                  crs = 4326, remove = FALSE)
  xy <- st_coordinates(st_transform(pts, cfg$equal_area_crs))
  dt[, `:=`(
    grid_x = floor(xy[, 1L] / cfg$thinning_grid_m),
    grid_y = floor(xy[, 2L] / cfg$thinning_grid_m)
  )]
  
  # Randomly retain one footprint in every spatial cell, then cap the stratum.
  dt <- dt[sample.int(.N)][, .SD[1L], by = .(grid_x, grid_y)]
  if (nrow(dt) > cfg$sample_per_class_year) {
    dt <- dt[sample.int(.N, cfg$sample_per_class_year)]
  }
  dt[, spatial_block_5km := paste(grid_x, grid_y, sep = "_")]
  dt[, c("grid_x", "grid_y") := NULL]
  samples[[i]] <- dt
  rm(dt, pts, xy); gc(FALSE)
}

sample_dt <- rbindlist(samples, use.names = TRUE, fill = TRUE)
setorder(sample_dt, year, architecture_class)
setnames(sample_dt,"sample_id","source_sample_id")
sample_dt[, sample_id := sprintf("XG_%06d", seq_len(.N))]
sample_dt[, shot_number := as.character(shot_number)]
sample_dt[, `:=`(lon_lowestmode=longitude,lat_lowestmode=latitude)]
if (!"mapbiomas_class" %in% names(sample_dt)) {
  if ("rfld_mapbiomas_class" %in% names(sample_dt)) {
    setnames(sample_dt, "rfld_mapbiomas_class", "mapbiomas_class")
  } else {
    sample_dt[, mapbiomas_class := NA_real_]
  }
} else if ("rfld_mapbiomas_class" %in% names(sample_dt)) {
  # Prefer the canonical name, but use the RFLD alias to fill missing values.
  sample_dt[is.na(mapbiomas_class),
            mapbiomas_class := rfld_mapbiomas_class[is.na(mapbiomas_class)]]
  sample_dt[, rfld_mapbiomas_class := NULL]
}
if(!"sensitivity"%in%names(sample_dt)) sample_dt[,sensitivity:=NA_real_]
if(!"solar_elevation"%in%names(sample_dt)) sample_dt[,solar_elevation:=NA_real_]
setcolorder(sample_dt, c("sample_id", "shot_number", "year",
                         "architecture_class", "lon_lowestmode",
                         "lat_lowestmode"))

if (anyDuplicated(sample_dt$sample_id)) stop("sample_id is not unique")
if (any(!is.finite(sample_dt$longitude)) ||
    any(!is.finite(sample_dt$latitude))) stop("Invalid coordinates")

fwrite(sample_dt, csv_out)
sample_sf <- st_as_sf(sample_dt,
                      coords = c("lon_lowestmode", "lat_lowestmode"),
                      crs = 4326, remove = FALSE)
if (file.exists(gpkg_out) && cfg$overwrite) file.remove(gpkg_out)
st_write(sample_sf, gpkg_out, layer = "gedi_sample", quiet = TRUE)

sample_summary <- sample_dt[, .(
  n = .N,
  unique_spatial_blocks = uniqueN(spatial_block_5km),
  median_sensitivity = median(sensitivity, na.rm = TRUE),
  median_rh98 = median(rh98_original, na.rm = TRUE)
), by = .(year, architecture_class)]
fwrite(sample_summary, summary_out)

message("STEP 02A COMPLETE")
message("Radiometric sample: ", csv_out)
message("Canonical coordinates are longitude/latitude; legacy aliases are retained.")
