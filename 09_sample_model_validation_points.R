#!/usr/bin/env Rscript

# Spatially independent, class-year-balanced random sampling for the Xingu revision.
#
# Creates one shared master sample (default maximum: 500 shots per
# architecture class x GEDI year), split into approximately 80% model and 20%
# validation observations. Entire 10-km cells are assigned to only one split,
# preventing local train-validation leakage. Within each pre-assigned split,
# footprints are sampled randomly within class x year strata. This preserves
# the observed within-class spatial distribution rather than giving every
# occupied cell equal representation. Split-specific inverse sampling weights
# recover the eligible class-year distribution within each split.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(sf)
})

# -----------------------------------------------------------------------------
# 1. Configuration: edit only this block
# -----------------------------------------------------------------------------

cfg <- list(
  # A Parquet file or a directory containing the classified GEDI Parquet files.
  input_path = "E:/Xingu_rev/data/step01b_final_classification/assignments_dataset",
  output_dir = "E:/Xingu_rev/data/step09_predictor_sample_v3",
  
  # Use the complete GEDI V3 temporal coverage and every selected class.
  years = 2019:2023,
  classes = NULL,
  
  target_per_class_year = 500L,
  validation_fraction = 0.20,
  spatial_cell_m = 10000,
  projected_crs = 5880,       # SIRGAS 2000 / Brazil Polyconic
  geographic_regions_x = 2L, # 2 x 5 = ten equally sized geographic subregions
  geographic_regions_y = 5L,
  seed = 20260824L
)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  paste0(..., collapse = "")))
}

stop_if_missing <- function(x, label) {
  if (!file.exists(x) && !dir.exists(x)) stop(label, " not found: ", x)
}

stop_if_missing(cfg$input_path, "Classified GEDI input")
set.seed(cfg$seed)

# -----------------------------------------------------------------------------
# 2. Read only the columns required for sampling
# -----------------------------------------------------------------------------

column_aliases <- list(
  architecture_class = c("architecture_class", "class", "final_class",
                         "cluster", "cluster_id"),
  year = c("year", "gedi_year", "acquisition_year"),
  longitude = c("longitude", "lon", "lon_lowestmode", "x"),
  latitude = c("latitude", "lat", "lat_lowestmode", "y"),
  shot_number = c("shot_number", "shotnumber", "shot_id", "sample_id")
)

choose_column <- function(available, aliases, canonical) {
  hit <- aliases[aliases %in% available]
  if (!length(hit)) {
    stop("Could not identify '", canonical, "'. Tried: ",
         paste(aliases, collapse = ", "))
  }
  hit[[1]]
}

is_parquet_input <- dir.exists(cfg$input_path) ||
  grepl("\\.parquet$", cfg$input_path, ignore.case = TRUE)

if (is_parquet_input) {
  ds <- arrow::open_dataset(cfg$input_path, format = "parquet")
  available <- names(ds$schema)
  selected <- vapply(names(column_aliases), function(nm) {
    choose_column(available, column_aliases[[nm]], nm)
  }, character(1))
  
  log_message("Reading sampling columns from Parquet")
  dt <- as.data.table(
    ds |>
      dplyr::select(dplyr::all_of(unname(selected))) |>
      dplyr::collect()
  )
  setnames(dt, unname(selected), names(selected))
} else {
  log_message("Reading classified GEDI table")
  header <- names(data.table::fread(cfg$input_path, nrows = 0L))
  selected <- vapply(names(column_aliases), function(nm) {
    choose_column(header, column_aliases[[nm]], nm)
  }, character(1))
  dt <- data.table::fread(cfg$input_path, select = unname(selected))
  setnames(dt, unname(selected), names(selected))
}

dt[, architecture_class := as.integer(architecture_class)]
dt[, year := as.integer(year)]
dt[, longitude := as.numeric(longitude)]
dt[, latitude := as.numeric(latitude)]
# Character storage prevents loss of precision for GEDI's long identifiers.
dt[, shot_number := as.character(shot_number)]

if (is.null(cfg$classes)) {
  cfg$classes <- sort(unique(dt$architecture_class[!is.na(dt$architecture_class)]))
}

dt <- dt[
  architecture_class %in% cfg$classes &
    year %in% cfg$years &
    is.finite(longitude) & is.finite(latitude) &
    longitude >= -180 & longitude <= 180 &
    latitude >= -90 & latitude <= 90 &
    !is.na(shot_number)
]
setorder(dt, shot_number)
dt <- unique(dt, by = "shot_number")
if (!nrow(dt)) stop("No eligible observations remain after filtering")

population_counts <- dt[, .(population_n_class_year = .N),
                        by = .(architecture_class, year)]

# -----------------------------------------------------------------------------
# 3. Construct 10-km cells and ten broad geographic subregions
# -----------------------------------------------------------------------------

log_message("Constructing spatial cells in EPSG:", cfg$projected_crs)
pts <- st_as_sf(dt, coords = c("longitude", "latitude"), crs = 4326,
                remove = FALSE)
pts_projected <- st_transform(pts, cfg$projected_crs)
xy <- st_coordinates(pts_projected)

dt[, x_projected := xy[, 1]]
dt[, y_projected := xy[, 2]]
dt[, cell_x := floor(x_projected / cfg$spatial_cell_m)]
dt[, cell_y := floor(y_projected / cfg$spatial_cell_m)]
dt[, spatial_cell_id := sprintf("C%07d_%07d", cell_x, cell_y)]

x_breaks <- seq(min(dt$x_projected), max(dt$x_projected),
                length.out = cfg$geographic_regions_x + 1L)
y_breaks <- seq(min(dt$y_projected), max(dt$y_projected),
                length.out = cfg$geographic_regions_y + 1L)
dt[, region_x := pmin(cfg$geographic_regions_x,
                      findInterval(x_projected, x_breaks,
                                   all.inside = TRUE))]
dt[, region_y := pmin(cfg$geographic_regions_y,
                      findInterval(y_projected, y_breaks,
                                   all.inside = TRUE))]
dt[, geographic_region := (region_y - 1L) * cfg$geographic_regions_x + region_x]

# Assign complete spatial cells to model or validation. Assignment is stratified
# among the ten geographic subregions, but is global across class and year.
cell_table <- unique(dt[, .(spatial_cell_id, geographic_region)])
cell_table[, random_order := runif(.N)]
cell_table[, cell_rank := frank(random_order, ties.method = "first"),
           by = geographic_region]
cell_table[, n_cells_region := .N, by = geographic_region]
cell_table[, analysis_split := fifelse(
  cell_rank <= pmax(1L, round(n_cells_region * cfg$validation_fraction)),
  "validation", "model"
)]
dt[cell_table, analysis_split := i.analysis_split, on = "spatial_cell_id"]

# These are the correct denominators for sampling probabilities. The model and
# validation pools contain different complete spatial cells, so each split has
# its own eligible class-year population size.
split_population_counts <- dt[, .(
  eligible_n_split_class_year = .N,
  eligible_spatial_cells = uniqueN(spatial_cell_id),
  eligible_geographic_regions = uniqueN(geographic_region)
), by = .(analysis_split, architecture_class, year)]

# -----------------------------------------------------------------------------
# 4. Sample each class-year stratum within the pre-assigned spatial splits
# -----------------------------------------------------------------------------

target_validation <- as.integer(round(
  cfg$target_per_class_year * cfg$validation_fraction
))
target_model <- as.integer(cfg$target_per_class_year - target_validation)

sample_one_stratum <- function(x, requested_n) {
  if (!nrow(x) || requested_n <= 0L) return(x[0])
  x <- data.table::copy(x)
  take <- min(requested_n, nrow(x))
  selected <- x[sample.int(nrow(x), take, replace = FALSE)]
  selected[, selection_pass := fifelse(
    take == nrow(x), "all_available", "random_within_split_class_year"
  )]
  selected
}

sample_split <- function(split_name, requested_n) {
  source <- dt[analysis_split == split_name]
  source[, sample_one_stratum(.SD, requested_n),
         by = .(architecture_class, year)]
}

log_message("Sampling model observations: target=", target_model,
            " per class-year")
model_dt <- sample_split("model", target_model)
log_message("Sampling validation observations: target=", target_validation,
            " per class-year")
validation_dt <- sample_split("validation", target_validation)

sampled <- rbindlist(list(model_dt, validation_dt), use.names = TRUE, fill = TRUE)

# -----------------------------------------------------------------------------
# 5. Add population-representation weights and audit the split
# -----------------------------------------------------------------------------

sampled[population_counts,
        population_n_class_year := i.population_n_class_year,
        on = .(architecture_class, year)]
sampled[split_population_counts,
        `:=`(
          eligible_n_split_class_year = i.eligible_n_split_class_year,
          eligible_spatial_cells = i.eligible_spatial_cells,
          eligible_geographic_regions = i.eligible_geographic_regions
        ),
        on = .(analysis_split, architecture_class, year)]
sampled[, sampled_n_class_year := .N,
        by = .(analysis_split, architecture_class, year)]
sampled[, inclusion_probability :=
          sampled_n_class_year / eligible_n_split_class_year]
sampled[, sampling_weight_raw :=
          eligible_n_split_class_year / sampled_n_class_year]
sampled[, sampling_weight_normalized :=
          sampling_weight_raw / mean(sampling_weight_raw),
        by = analysis_split]

sampled[, sample_id := sprintf("XG_%06d", seq_len(.N))]
setcolorder(sampled, c(
  "sample_id", "shot_number", "architecture_class", "year",
  "analysis_split", "longitude", "latitude", "spatial_cell_id",
  "geographic_region", "selection_pass", "population_n_class_year",
  "eligible_n_split_class_year",
  "sampled_n_class_year", "inclusion_probability",
  "sampling_weight_raw", "sampling_weight_normalized"
))

overlap_cells <- intersect(
  sampled[analysis_split == "model", unique(spatial_cell_id)],
  sampled[analysis_split == "validation", unique(spatial_cell_id)]
)
if (length(overlap_cells)) stop("Spatial leakage detected between splits")

audit <- sampled[, .(
  population_n = unique(population_n_class_year),
  eligible_n_in_split = unique(eligible_n_split_class_year),
  sampled_n = .N,
  sampling_fraction_within_split = unique(inclusion_probability),
  unique_spatial_cells = uniqueN(spatial_cell_id),
  unique_geographic_regions = uniqueN(geographic_region),
  all_available_n = sum(selection_pass == "all_available"),
  random_sample_n = sum(selection_pass == "random_within_split_class_year"),
  raw_weight = unique(sampling_weight_raw),
  normalized_weight = unique(sampling_weight_normalized)
), by = .(analysis_split, architecture_class, year)]

split_summary <- sampled[, .(
  observations = .N,
  classes = uniqueN(architecture_class),
  years = uniqueN(year),
  spatial_cells = uniqueN(spatial_cell_id),
  geographic_regions = uniqueN(geographic_region)
), by = analysis_split]

# Compare the spatial distribution of the random sample with its eligible pool.
# Landscape predictors are not available until Step 10, but these summaries
# reveal any major geographic distortion before expensive extraction begins.
spatial_summary <- function(x, source_label) {
  x[, .(
    source = source_label,
    n = .N,
    spatial_cells = uniqueN(spatial_cell_id),
    geographic_regions = uniqueN(geographic_region),
    longitude_mean = mean(longitude),
    longitude_sd = sd(longitude),
    longitude_q05 = quantile(longitude, 0.05),
    longitude_median = median(longitude),
    longitude_q95 = quantile(longitude, 0.95),
    latitude_mean = mean(latitude),
    latitude_sd = sd(latitude),
    latitude_q05 = quantile(latitude, 0.05),
    latitude_median = median(latitude),
    latitude_q95 = quantile(latitude, 0.95)
  ), by = .(analysis_split, architecture_class, year)]
}

spatial_distribution_audit <- rbindlist(list(
  spatial_summary(dt, "eligible_population"),
  spatial_summary(sampled, "sample")
), use.names = TRUE)

population_region <- dt[, .(eligible_n_region = .N),
                        by = .(analysis_split, architecture_class, year,
                               geographic_region)]
population_region[, eligible_region_proportion :=
                    eligible_n_region / sum(eligible_n_region),
                  by = .(analysis_split, architecture_class, year)]
sample_region <- sampled[, .(sampled_n_region = .N),
                         by = .(analysis_split, architecture_class, year,
                                geographic_region)]
sample_region[, sampled_region_proportion :=
                sampled_n_region / sum(sampled_n_region),
              by = .(analysis_split, architecture_class, year)]
region_distribution_audit <- merge(
  population_region, sample_region,
  by = c("analysis_split", "architecture_class", "year", "geographic_region"),
  all = TRUE
)
region_distribution_audit[is.na(sampled_n_region), sampled_n_region := 0L]
region_distribution_audit[is.na(sampled_region_proportion),
                          sampled_region_proportion := 0]
region_distribution_audit[, proportion_difference :=
                            sampled_region_proportion - eligible_region_proportion]

# -----------------------------------------------------------------------------
# 6. Write repository-friendly tables and spatial files
# -----------------------------------------------------------------------------

columns_to_drop <- intersect(
  c("x_projected", "y_projected", "cell_x", "cell_y", "region_x", "region_y"),
  names(sampled)
)
sampled[, (columns_to_drop) := NULL]

model_out <- sampled[analysis_split == "model"]
validation_out <- sampled[analysis_split == "validation"]

fwrite(sampled, file.path(cfg$output_dir, "predictor_sample_all.csv"))
fwrite(model_out, file.path(cfg$output_dir, "predictor_sample_model.csv"))
fwrite(validation_out,
       file.path(cfg$output_dir, "predictor_sample_validation.csv"))
fwrite(audit, file.path(cfg$output_dir, "predictor_sample_audit.csv"))
fwrite(split_summary,
       file.path(cfg$output_dir, "predictor_sample_split_summary.csv"))
fwrite(population_counts,
       file.path(cfg$output_dir, "eligible_population_class_year_counts.csv"))
fwrite(split_population_counts,
       file.path(cfg$output_dir,
                 "eligible_population_split_class_year_counts.csv"))
fwrite(spatial_distribution_audit,
       file.path(cfg$output_dir,
                 "predictor_sample_spatial_distribution_audit.csv"))
fwrite(region_distribution_audit,
       file.path(cfg$output_dir,
                 "predictor_sample_region_distribution_audit.csv"))

spatial_out <- st_as_sf(sampled, coords = c("longitude", "latitude"),
                        crs = 4326, remove = FALSE)
st_write(spatial_out,
         file.path(cfg$output_dir, "predictor_sample_all.gpkg"),
         layer = "predictor_sample", delete_dsn = TRUE, quiet = TRUE)

log_message("Completed. Total sampled observations=", nrow(sampled),
            "; model=", nrow(model_out),
            "; validation=", nrow(validation_out))
print(audit)
print(split_summary)
