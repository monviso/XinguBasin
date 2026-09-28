#!/usr/bin/env Rscript

# STEP 04 — Extract analysis-ready Sentinel-1 OPERA RTC-S1 at a balanced
# subset of the GEDI architecture sample. No Google Earth Engine or SNAP is
# required. OPERA values are gamma-0 power at 30 m and are supplied as COGs.
#
# Authentication is handled by ASF's official `asf_search` client through
# reticulate. Set only the Earthdata username below; the password is requested
# interactively and is never written to this script or an output file.

options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  earthdata_user = "Your_EarthData_Username",
  sample_csv = "~/data/step02_radiometry_v3/gedi_radiometric_sample.csv",
  output_dir = "~/data/step04_sentinel1_opera_v3",
  sample_per_class_year = 30L,
  max_acquisitions_per_point_year = 12L,
  extraction_buffer_m = 45,
  seed = 42L,
  collection_concept_id = "C2777436413-ASF",
  cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
  overwrite = FALSE,
  remote_cog = FALSE,
  allow_full_file_fallback = TRUE,
  retain_downloaded_cogs = FALSE,
  max_consecutive_failures = 3L
)

required <- c("data.table", "httr2", "terra", "jsonlite", "reticulate", "getPass")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(httr2); library(terra); library(jsonlite)
  library(reticulate); library(getPass)
})
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(cfg$output_dir, "cmr_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cog_cache_dir <- file.path(cfg$output_dir, "temporary_cog_cache")
dir.create(cog_cache_dir, recursive = TRUE, showWarnings = FALSE)

if (cfg$earthdata_user == "REPLACE_WITH_EARTHDATA_USERNAME") {
  stop("Set cfg$earthdata_user to the Earthdata username used successfully in ASF Vertex")
}
if (!reticulate::py_module_available("asf_search")) {
  stop("Python package 'asf-search' is unavailable. Run once in R:\n",
       "  reticulate::py_install(\"asf-search\", pip = TRUE)")
}
asf <- reticulate::import("asf_search", delay_load = FALSE)
asf_session <- asf$ASFSession()
earthdata_password <- getPass::getPass("Earthdata password for ASF: ")
tryCatch(
  asf_session$auth_with_creds(cfg$earthdata_user, earthdata_password),
  error = function(e) stop("ASF authentication failed: ", conditionMessage(e))
)
rm(earthdata_password)

sample_dt <- fread(cfg$sample_csv)
needed <- c("sample_id", "year", "architecture_class",
            "lon_lowestmode", "lat_lowestmode", "spatial_block_5km")
absent <- setdiff(needed, names(sample_dt))
if (length(absent)) stop("Sample is missing: ", paste(absent, collapse = ", "))
sample_dt[, sample_id := as.character(sample_id)]
set.seed(cfg$seed)
sar_sample <- sample_dt[, {
  take <- min(.N, cfg$sample_per_class_year)
  .SD[sample.int(.N, take)]
}, by = .(year, architecture_class)]
setorder(sar_sample, year, architecture_class, sample_id)
fwrite(sar_sample, file.path(cfg$output_dir, "sentinel1_sample.csv"))

observation_file <- file.path(cfg$output_dir, "sentinel1_observations.csv")
status_file <- file.path(cfg$output_dir, "sentinel1_extraction_status.csv")
if (cfg$overwrite) {
  if (file.exists(observation_file)) file.remove(observation_file)
  if (file.exists(status_file)) file.remove(status_file)
}
observations <- if (file.exists(observation_file)) fread(observation_file) else data.table()
status <- if (file.exists(status_file)) fread(status_file) else data.table(
  sample_id = character(), year = integer(), architecture_class = integer(),
  status = character(), acquisitions_extracted = integer(), message = character()
)
completed_ids <- if (nrow(status)) status[status == "completed", sample_id] else character()

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
safe_chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1L]])

query_cmr <- function(point_row) {
  cache_file <- file.path(cache_dir, paste0(point_row$sample_id, "_", point_row$year, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  temporal <- sprintf("%d-01-01T00:00:00Z,%d-12-31T23:59:59Z",
                      point_row$year, point_row$year)
  point <- sprintf("%.7f,%.7f", point_row$lon_lowestmode,
                   point_row$lat_lowestmode)
  response <- request(cfg$cmr_url) |>
    req_url_query(collection_concept_id = cfg$collection_concept_id,
                  point = point, temporal = temporal, page_size = 2000) |>
    req_retry(max_tries = 4) |>
    req_perform()
  payload <- resp_body_json(response, simplifyVector = FALSE)
  items <- payload$items %||% list()
  saveRDS(items, cache_file)
  items
}

attribute_value <- function(umm, possible_names) {
  attributes <- umm$AdditionalAttributes %||% list()
  for (attribute in attributes) {
    name <- toupper(safe_chr(attribute$Name))
    if (name %in% toupper(possible_names)) {
      values <- attribute$Values %||% list()
      return(safe_chr(values))
    }
  }
  NA_character_
}

items_to_candidates <- function(items) {
  rows <- lapply(items, function(item) {
    umm <- item$umm %||% list()
    granule <- safe_chr(umm$GranuleUR)
    urls <- umm$RelatedUrls %||% list()
    hrefs <- vapply(urls, function(link) safe_chr(link$URL), character(1))
    hrefs <- hrefs[!is.na(hrefs)]
    vv <- hrefs[grepl("_VV\\.tif($|\\?)", hrefs, ignore.case = TRUE)]
    vh <- hrefs[grepl("_VH\\.tif($|\\?)", hrefs, ignore.case = TRUE)]
    if (!length(vv) || !length(vh)) return(NULL)
    date_text <- sub(".*_([0-9]{8}T[0-9]{6}Z)_.*", "\\1", granule)
    acquisition_time <- as.POSIXct(date_text, format = "%Y%m%dT%H%M%SZ", tz = "UTC")
    burst <- sub(".*_(T[0-9]{3}-[0-9]{6}-IW[123])_.*", "\\1", granule)
    if (!grepl("^T[0-9]{3}-", burst)) burst <- NA_character_
    direction <- attribute_value(umm,
                                 c("ASCENDING_DESCENDING", "ORBIT_DIRECTION", "FLIGHT_DIRECTION"))
    data.table(granule = granule, acquisition_time = acquisition_time,
               burst_id = burst, orbit_direction = direction,
               vv_url = vv[1L], vh_url = vh[1L])
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

choose_consistent_timeseries <- function(candidates) {
  candidates <- unique(candidates[!is.na(acquisition_time) & !is.na(burst_id)],
                       by = c("granule", "vv_url", "vh_url"))
  if (!nrow(candidates)) return(candidates)
  dominant <- candidates[, .N, by = burst_id][order(-N), burst_id[1L]]
  candidates <- candidates[burst_id == dominant][order(acquisition_time)]
  n_keep <- min(nrow(candidates), cfg$max_acquisitions_per_point_year)
  indices <- unique(as.integer(round(seq(1, nrow(candidates), length.out = n_keep))))
  candidates[indices]
}

download_authenticated_cog <- function(url) {
  clean_name <- sub("[?].*$", "", basename(url))
  destination <- file.path(cog_cache_dir, clean_name)
  if (file.exists(destination) && file.info(destination)$size > 1024) {
    return(destination)
  }
  temporary <- paste0(destination, ".part")
  if (file.exists(temporary)) file.remove(temporary)
  response <- asf_session$get(url, timeout = 180L)
  response$raise_for_status()
  content <- response$content
  if (inherits(content, "python.builtin.object")) {
    content <- reticulate::py_to_r(content)
  }
  if (!is.raw(content)) content <- as.raw(content)
  writeBin(content, temporary)
  if (!file.exists(temporary) || file.info(temporary)$size <= 1024) {
    stop("Authenticated download did not return a valid-sized TIFF")
  }
  if (!file.rename(temporary, destination)) stop("Could not finalize cached COG")
  destination
}

read_cog_neighbourhood <- function(url, lon, lat) {
  extract_from_path <- function(path) {
    raster <- rast(path)
    point <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"),
                  crs = "EPSG:4326")
    point <- project(point, crs(raster))
    values <- extract(raster, point, buffer = cfg$extraction_buffer_m,
                      ID = FALSE)[[1L]]
    values <- values[is.finite(values) & values > 0]
    if (!length(values)) return(c(mean_power = NA_real_, sd_power = NA_real_, n_pixels = 0))
    c(mean_power = mean(values), sd_power = sd(values), n_pixels = length(values))
  }
  if (cfg$remote_cog) {
    remote <- tryCatch({extract_from_path(paste0("/vsicurl/", url))},
                       error = function(e) NULL)
    if (!is.null(remote) && is.finite(remote[["mean_power"]])) return(remote)
  }
  if (!cfg$allow_full_file_fallback) {
    return(c(mean_power = NA_real_, sd_power = NA_real_, n_pixels = 0))
  }
  destination <- download_authenticated_cog(url)
  extract_from_path(destination)
}

consecutive_failures <- 0L
for (i in seq_len(nrow(sar_sample))) {
  point_row <- sar_sample[i]
  if (point_row$sample_id %in% completed_ids) next
  message(sprintf("[%d/%d] %s; year=%d; class=%s", i, nrow(sar_sample),
                  point_row$sample_id, point_row$year,
                  point_row$architecture_class))
  outcome <- tryCatch({
    candidates <- items_to_candidates(query_cmr(point_row))
    selected <- choose_consistent_timeseries(candidates)
    if (!nrow(selected)) stop("No paired VV/VH OPERA acquisitions found")
    extracted <- rbindlist(lapply(seq_len(nrow(selected)), function(j) {
      vv <- read_cog_neighbourhood(selected$vv_url[j], point_row$lon_lowestmode,
                                   point_row$lat_lowestmode)
      vh <- read_cog_neighbourhood(selected$vh_url[j], point_row$lon_lowestmode,
                                   point_row$lat_lowestmode)
      if (!is.finite(vv[["mean_power"]]) || !is.finite(vh[["mean_power"]])) return(NULL)
      vv_db <- 10 * log10(vv[["mean_power"]])
      vh_db <- 10 * log10(vh[["mean_power"]])
      data.table(
        sample_id = point_row$sample_id, year = point_row$year,
        architecture_class = point_row$architecture_class,
        spatial_block_5km = point_row$spatial_block_5km,
        lon = point_row$lon_lowestmode, lat = point_row$lat_lowestmode,
        acquisition_time = selected$acquisition_time[j],
        month = as.integer(format(selected$acquisition_time[j], "%m")),
        burst_id = selected$burst_id[j],
        orbit_direction = selected$orbit_direction[j],
        VV_power = vv[["mean_power"]], VH_power = vh[["mean_power"]],
        VV_dB = vv_db, VH_dB = vh_db, VH_minus_VV_dB = vh_db - vv_db,
        RVI = 4 * vh[["mean_power"]] / (vv[["mean_power"]] + vh[["mean_power"]]),
        VV_spatial_sd_power = vv[["sd_power"]],
        VH_spatial_sd_power = vh[["sd_power"]],
        VV_n_pixels = vv[["n_pixels"]], VH_n_pixels = vh[["n_pixels"]]
      )
    }), use.names = TRUE, fill = TRUE)
    if (!nrow(extracted)) stop("COGs found, but no valid backscatter was extracted")
    observations <<- rbind(observations, extracted, use.names = TRUE, fill = TRUE)
    fwrite(observations, observation_file)
    list(status = "completed", n = nrow(extracted), message = "ok")
  }, error = function(e) list(status = "failed", n = 0L,
                              message = conditionMessage(e)))
  status <- status[sample_id != point_row$sample_id]
  status <- rbind(status, data.table(
    sample_id = point_row$sample_id, year = point_row$year,
    architecture_class = point_row$architecture_class,
    status = outcome$status, acquisitions_extracted = outcome$n,
    message = outcome$message
  ), fill = TRUE)
  fwrite(status, status_file)
  if (outcome$status == "failed") {
    consecutive_failures <- consecutive_failures + 1L
    warning(point_row$sample_id, ": ", outcome$message)
  } else {
    consecutive_failures <- 0L
  }
  if (!cfg$retain_downloaded_cogs) {
    cached_files <- list.files(cog_cache_dir, full.names = TRUE, all.files = FALSE)
    if (length(cached_files)) unlink(cached_files, force = TRUE)
  }
  if (consecutive_failures >= cfg$max_consecutive_failures) {
    stop(
      "Stopped after ", consecutive_failures, " consecutive failures. ",
      "When all points return HTTP 403/500, authorize the ASF application ",
      "in Earthdata Login, accept the ASF EULA, set the required study area, ",
      "then generate a new Earthdata token and restart this script."
    )
  }
}

message("STEP 04 COMPLETE")
message("Observations: ", observation_file)
message("Status: ", status_file)

