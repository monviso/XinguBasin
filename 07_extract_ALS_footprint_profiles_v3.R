#!/usr/bin/env Rscript

# STEP 07 — Extract independent airborne-lidar (ALS) return-height profiles
# inside nominal GEDI footprints. This script does not fit or change the GEDI
# classification. It records the ALS/GEDI acquisition-year mismatch for every
# footprint and retains it as a reported limitation rather than a filter.

options(stringsAsFactors = FALSE, warn = 1)

cfg <- list(
  project_dir = "E:/Xingu_rev",
  matches_file = "E:/Xingu_rev/data/step06_GEDI_ALS_overlap_v3/GEDI_ALS_matches.csv",
  laz_dir = "E:/Xingu_rev/data/ALS",
  output_dir = "E:/Xingu_rev/data/step07_ALS_profiles_v3",
  # Set NULL to process every matched shot. The default is a balanced external
  # validation sample and limits expensive repeated LAZ reads.
  max_shots_per_class = 150L,
  seed = 48291L,
  footprint_radius_m = 12.5,
  processing_radius_m = 50,
  coverage_cell_m = 2,
  min_point_density_m2 = 1,
  min_coverage_fraction = 0.75,
  min_ground_points = 3L,
  min_returns_in_footprint = 50L,
  maximum_normalized_height_m = 80,
  classify_ground_if_missing = TRUE,
  overwrite = FALSE
)

required <- c("data.table", "lidR", "sf")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(lidR)
  library(sf)
})
sf::sf_use_s2(FALSE)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
profile_csv <- file.path(cfg$output_dir, "ALS_footprint_profiles.csv")
qa_csv <- file.path(cfg$output_dir, "ALS_footprint_QA.csv")
sample_csv <- file.path(cfg$output_dir, "ALS_validation_sample.csv")
summary_csv <- file.path(cfg$output_dir, "ALS_extraction_summary.csv")

outputs <- c(profile_csv, qa_csv, sample_csv, summary_csv)
if (!cfg$overwrite && any(file.exists(outputs))) {
  stop("Step 07 output already exists. Set cfg$overwrite=TRUE to rebuild it.")
}

log_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

matches <- fread(cfg$matches_file, colClasses = list(character = c("shot_number", "path2")))
needed <- c("shot_number", "longitude", "latitude", "year",
            "architecture_class", "als_dataset", "als_year", "path2")
absent <- setdiff(needed, names(matches))
if (length(absent)) stop("Matches file is missing: ", paste(absent, collapse = ", "))

matches[, local_laz := file.path(cfg$laz_dir, basename(path2))]
matches[, local_laz_exists := file.exists(local_laz)]
if (!any(matches$local_laz_exists)) {
  stop("None of the matched LAZ basenames was found in ", cfg$laz_dir)
}
matches[, gedi_year := as.integer(year)]
matches[, architecture_class := as.integer(architecture_class)]
matches[, absolute_year_difference := abs(gedi_year - as.integer(als_year))]

# Select GEDI shots, not match rows: a shot can legitimately intersect more
# than one LAZ tile. Sampling is balanced by GEDI architecture class.
shot_table <- matches[, .(
  longitude = first(longitude),
  latitude = first(latitude),
  gedi_year = first(gedi_year),
  architecture_class = first(architecture_class),
  minimum_absolute_year_difference = suppressWarnings(min(absolute_year_difference, na.rm = TRUE)),
  candidate_laz_files = uniqueN(local_laz[local_laz_exists])
), by = shot_number]
shot_table[!is.finite(minimum_absolute_year_difference), minimum_absolute_year_difference := NA_real_]
# Do not allow missing local files to consume the per-class sampling quota.
shot_table <- shot_table[candidate_laz_files > 0]

set.seed(cfg$seed)
if (!is.null(cfg$max_shots_per_class)) {
  shot_table <- shot_table[, .SD[sample.int(.N, min(.N, cfg$max_shots_per_class))],
                           by = architecture_class]
}
setorder(shot_table, architecture_class, shot_number)
fwrite(shot_table, sample_csv)
selected_ids <- shot_table$shot_number
matches <- matches[shot_number %in% selected_ids & local_laz_exists == TRUE]
setorder(matches, shot_number, absolute_year_difference, als_dataset, local_laz)

header_crs <- function(path) {
  h <- lidR::readLASheader(path)
  crs <- sf::st_crs(h)
  if (is.na(crs)) stop("LAZ has no readable CRS")
  crs
}

footprint_coverage <- function(x, y, cx, cy, radius, cell_size) {
  gx <- floor((x - (cx - radius)) / cell_size)
  gy <- floor((y - (cy - radius)) / cell_size)
  occupied <- unique(paste(gx, gy, sep = "_"))
  grid_x <- seq(cx - radius + cell_size / 2, cx + radius - cell_size / 2,
                by = cell_size)
  grid_y <- seq(cy - radius + cell_size / 2, cy + radius - cell_size / 2,
                by = cell_size)
  grid <- CJ(x = grid_x, y = grid_y)
  grid <- grid[(x - cx)^2 + (y - cy)^2 <= radius^2]
  possible <- unique(paste(floor((grid$x - (cx - radius)) / cell_size),
                           floor((grid$y - (cy - radius)) / cell_size), sep = "_"))
  if (!length(possible)) return(NA_real_)
  length(intersect(occupied, possible)) / length(possible)
}

extract_candidate <- function(row) {
  path <- row$local_laz
  tryCatch({
    crs <- header_crs(path)
    centre <- st_as_sf(data.frame(lon = row$longitude,
                                  lat = row$latitude),
                       coords = c("lon", "lat"), crs = 4326)
    centre <- st_transform(centre, crs)
    xy <- st_coordinates(centre)[1L, ]
    filter <- sprintf("-inside_circle %.3f %.3f %.3f",
                      xy[1L], xy[2L], cfg$processing_radius_m)
    las <- lidR::readLAS(path, filter = filter)
    if (lidR::is.empty(las)) stop("no ALS returns in processing circle")
    
    dat <- las@data
    has_class <- "Classification" %in% names(dat)
    ground_n <- if (has_class) sum(dat$Classification %in% c(2L, 9L), na.rm = TRUE) else 0L
    ground_derived <- FALSE
    if (ground_n < cfg$min_ground_points && cfg$classify_ground_if_missing) {
      las <- lidR::classify_ground(las, lidR::csf())
      dat <- las@data
      ground_n <- sum(dat$Classification %in% c(2L, 9L), na.rm = TRUE)
      ground_derived <- TRUE
    }
    if (ground_n < cfg$min_ground_points) stop("insufficient ALS ground returns")
    
    las <- lidR::normalize_height(
      las, lidR::knnidw(k = 10L, p = 2), use_class = c(2L, 9L)
    )
    dat <- las@data
    inside <- (dat$X - xy[1L])^2 + (dat$Y - xy[2L])^2 <= cfg$footprint_radius_m^2
    z <- as.numeric(dat$Z[inside])
    x <- as.numeric(dat$X[inside])
    y <- as.numeric(dat$Y[inside])
    classification <- if ("Classification" %in% names(dat)) {
      as.integer(dat$Classification[inside])
    } else rep(NA_integer_, sum(inside))
    valid <- is.finite(z) & z >= -1 & z <= cfg$maximum_normalized_height_m
    z <- pmax(z[valid], 0)
    x <- x[valid]
    y <- y[valid]
    classification <- classification[valid]
    n <- length(z)
    if (n < cfg$min_returns_in_footprint) stop("too few ALS returns in GEDI footprint")
    
    area <- pi * cfg$footprint_radius_m^2
    density <- n / area
    coverage <- footprint_coverage(x, y, xy[1L], xy[2L],
                                   cfg$footprint_radius_m, cfg$coverage_cell_m)
    # Conventional ALS structure is calculated from vegetation returns. Ground,
    # water and noise classes are excluded; unclassified returns are retained
    # because many legacy ALS acquisitions do not classify vegetation strata.
    vegetation <- z >= 0.5 &
      (is.na(classification) | !classification %in% c(2L, 7L, 9L, 18L))
    z_veg <- z[vegetation]
    if (length(z_veg) < 20L) stop("too few vegetation returns in GEDI footprint")
    
    # ALS profile for visualization/sensitivity only: RH0 is a ground anchor,
    # while RH1-RH98 are quantiles of above-ground vegetation returns.
    rh <- c(0, as.numeric(quantile(z_veg, probs = (1:98) / 100,
                                   type = 7, na.rm = TRUE, names = FALSE)))
    names(rh) <- sprintf("als_rh%02d", 0:98)
    
    height_probabilities <- c(.50, .75, .90, .95, .98)
    height_quantiles <- as.numeric(quantile(
      z_veg, probs = height_probabilities, type = 7, names = FALSE
    ))
    names(height_quantiles) <- paste0("vegetation_height_p", height_probabilities * 100)
    cover_thresholds <- c(2, 5, 10, 20, 30)
    canopy_cover <- vapply(cover_thresholds, function(h) mean(z >= h), numeric(1))
    names(canopy_cover) <- paste0("canopy_cover_gt", cover_thresholds, "m")
    
    strata_breaks <- c(0.5, 2, 5, 10, 20, 30, Inf)
    strata_counts <- tabulate(findInterval(z_veg, strata_breaks,
                                           rightmost.closed = TRUE), nbins = 6L)
    strata_proportions <- strata_counts / sum(strata_counts)
    positive_strata <- strata_proportions[strata_proportions > 0]
    foliage_height_diversity <- -sum(positive_strata * log(positive_strata)) / log(6)
    
    passed <- density >= cfg$min_point_density_m2 &&
      coverage >= cfg$min_coverage_fraction
    
    c(list(
      status = if (passed) "passed" else "failed_QA",
      failure_reason = if (passed) NA_character_ else "density_or_coverage_below_threshold",
      local_laz = path,
      als_dataset = row$als_dataset,
      als_year = as.integer(row$als_year),
      absolute_year_difference = as.numeric(row$absolute_year_difference),
      n_returns = n,
      point_density_m2 = density,
      coverage_fraction = coverage,
      ground_returns_processing_buffer = ground_n,
      ground_classification_derived = ground_derived,
      maximum_height_m = max(z),
      vegetation_returns = length(z_veg),
      vegetation_return_fraction = mean(vegetation),
      mean_vegetation_height_m = mean(z_veg),
      sd_vegetation_height_m = sd(z_veg),
      foliage_height_diversity = foliage_height_diversity,
      gap_fraction_below_2m = 1 - canopy_cover[["canopy_cover_gt2m"]],
      qa_score = coverage + pmin(density / 20, 1)
    ), as.list(height_quantiles), as.list(canopy_cover), as.list(rh))
  }, error = function(e) list(
    status = "failed", failure_reason = conditionMessage(e),
    local_laz = path, als_dataset = row$als_dataset,
    als_year = as.integer(row$als_year),
    absolute_year_difference = as.numeric(row$absolute_year_difference),
    n_returns = NA_integer_, point_density_m2 = NA_real_,
    coverage_fraction = NA_real_, ground_returns_processing_buffer = NA_integer_,
    ground_classification_derived = NA, maximum_height_m = NA_real_,
    qa_score = -Inf
  ))
}

qa_rows <- vector("list", length(selected_ids))
profile_rows <- vector("list", length(selected_ids))
for (i in seq_along(selected_ids)) {
  id <- selected_ids[i]
  shot <- shot_table[shot_number == id][1L]
  candidates <- matches[shot_number == id]
  log_message(sprintf("[%d/%d] shot=%s; class=%d; candidate LAZ=%d",
                      i, length(selected_ids), id, shot$architecture_class,
                      nrow(candidates)))
  if (!nrow(candidates)) {
    attempts <- data.table(status = "failed", failure_reason = "local LAZ missing")
  } else {
    attempts <- rbindlist(lapply(seq_len(nrow(candidates)), function(j) {
      as.data.table(extract_candidate(candidates[j]))
    }), fill = TRUE)
  }
  attempts[, `:=`(
    shot_number = id,
    gedi_year = shot$gedi_year,
    architecture_class = shot$architecture_class,
    longitude = shot$longitude,
    latitude = shot$latitude
  )]
  attempts[, selected_candidate := FALSE]
  passed_idx <- which(attempts$status == "passed")
  if (length(passed_idx)) {
    best <- passed_idx[which.max(attempts$qa_score[passed_idx])]
    attempts[best, selected_candidate := TRUE]
    profile_rows[[i]] <- attempts[best]
  }
  qa_rows[[i]] <- attempts
  if (i %% 25L == 0L) gc(FALSE)
}

qa <- rbindlist(qa_rows, fill = TRUE)
profiles <- rbindlist(profile_rows, fill = TRUE)
fwrite(qa, qa_csv)
if (!nrow(profiles)) stop("No ALS footprint passed QA; inspect ", qa_csv)
profiles[, qa_score := NULL]
fwrite(profiles, profile_csv)

summary_dt <- profiles[, .(
  retained_footprints = .N,
  median_point_density_m2 = median(point_density_m2),
  median_coverage_fraction = median(coverage_fraction),
  median_absolute_year_difference = safe_median(absolute_year_difference),
  minimum_absolute_year_difference = suppressWarnings(min(absolute_year_difference, na.rm = TRUE)),
  maximum_absolute_year_difference = suppressWarnings(max(absolute_year_difference, na.rm = TRUE))
), by = .(architecture_class, als_dataset, als_year, gedi_year)]
summary_dt[!is.finite(minimum_absolute_year_difference), minimum_absolute_year_difference := NA_real_]
summary_dt[!is.finite(maximum_absolute_year_difference), maximum_absolute_year_difference := NA_real_]
setorder(summary_dt, architecture_class, als_dataset, als_year, gedi_year)
fwrite(summary_dt, summary_csv)

log_message("STEP 07 COMPLETE; retained footprints=", nrow(profiles),
            "; attempted shots=", length(selected_ids))
log_message("Profiles: ", profile_csv)
log_message("QA audit: ", qa_csv)

