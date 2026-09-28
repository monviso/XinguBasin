#!/usr/bin/env Rscript

# Download GEDI01_B V002 for the small representative subset selected by
# Script 23, extract received waveforms, and save a tidy restartable dataset.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(rGEDI)
})

cfg <- list(
  selection_file = "E:/Xingu_rev/data/step01b_final_classification/figures/representative_waveform_shots.csv",
  output_file = "E:/Xingu_rev/data/step01b_final_classification/figures/L1B_waveforms_tidy.parquet",
  qa_file = "E:/Xingu_rev/data/step01b_final_classification/figures/L1B_waveform_retrieval_QA.csv",
  cache_dir = "E:/Xingu_rev/data/step01b_final_classification/L1B_cache",
  product = "GEDI01_B",
  version = "002",
  l2a_product = "GEDI02_A",
  # rGEDI currently resolves the V002 collection concept, not V003. V002 is
  # used only to obtain absolute lowest-mode ground elevation for aligning the
  # corresponding L1B waveform; all plotted RH metrics remain the V003 values.
  l2a_version = "002",
  date_buffer_days = 2L,
  overwrite = TRUE,
  delete_h5_after_success = TRUE
)

dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
l2_cache_dir <- file.path(cfg$cache_dir, "L2A_ground")
dir.create(l2_cache_dir, recursive = TRUE, showWarnings = FALSE)

# rGEDI searches for .netrc inside every download output directory. Create it
# once and reuse it in both caches, avoiding repeated interactive prompts.
setup_earthdata_netrc <- function(download_dirs) {
  netrc_paths <- file.path(download_dirs, ".netrc")
  existing <- netrc_paths[file.exists(netrc_paths)]
  if (length(existing)) {
    credentials <- readLines(existing[1], warn = FALSE)
  } else {
    user <- Sys.getenv("NASA_USER", unset = "")
    password <- Sys.getenv("NASA_PASSWORD", unset = "")
    if (!nzchar(user)) user <- getPass::getPass("NASA Earthdata username: ")
    if (!nzchar(password)) password <- getPass::getPass("NASA Earthdata password: ")
    if (!nzchar(user) || !nzchar(password)) stop("Earthdata credentials were not supplied")
    credentials <- c(
      "machine urs.earthdata.nasa.gov",
      sprintf("login %s", user),
      sprintf("password %s", password)
    )
  }
  for (p in netrc_paths) {
    writeLines(credentials, p, useBytes = TRUE)
    Sys.chmod(p, mode = "0600")
  }
  invisible(netrc_paths)
}

setup_earthdata_netrc(c(cfg$cache_dir, l2_cache_dir))

first_name <- function(nms, candidates, required = TRUE) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) return(hit[1])
  if (required) stop("None of these columns was found: ", paste(candidates, collapse = ", "))
  NULL
}

# Force 17-digit GEDI identifiers to character at ingestion. Allowing fread to
# infer a numeric type can lose final digits before any later conversion.
sel <- fread(
  cfg$selection_file,
  colClasses = list(character = c("sample_id", "shot_number"))
)
class_col <- first_name(names(sel), c("architecture_class", "class", "final_class"))
shot_col <- first_name(names(sel), c("shot_number", "shot_num", "shot"))
sample_col <- first_name(names(sel), c("sample_id"), required = FALSE)
lat_col <- first_name(names(sel), c("latitude", "lat_lowestmode", "lat_lm", "lat"))
lon_col <- first_name(names(sel), c("longitude", "lon_lowestmode", "lon_lm", "lon"))
date_col <- first_name(names(sel), c("date", "acquisition_date", "acquisition_datetime", "acquisition_time",
                                     "datetime", "time", "delta_time"))
ground_col <- first_name(names(sel), c(
  "elevation_lowestmode", "elev_lowestmode", "elev_lm", "elev_lm_a0"
), required = FALSE)
raw_date <- sel[[date_col]]
date_was_delta <- identical(date_col, "delta_time")
old_names <- c(class_col, shot_col, lat_col, lon_col, date_col)
new_names <- c("architecture_class", "shot_number", "latitude", "longitude", "acquisition_date")
if (!is.null(ground_col)) {
  old_names <- c(old_names, ground_col)
  new_names <- c(new_names, "ground_elevation_m")
}
setnames(sel, old_names, new_names)
if (!("ground_elevation_m" %in% names(sel))) sel[, ground_elevation_m := NA_real_]
# sample_id was created from the exact GEDI identifier before any conversion
# to double precision. Prefer it whenever present; large numeric shot_number
# values can otherwise be rounded (e.g. ...637 becomes ...600).
if (!is.null(sample_col) && sample_col %in% names(sel)) {
  sel[, shot_key := as.character(get(sample_col))]
} else {
  sel[, shot_key := as.character(shot_number)]
  warning("sample_id is absent; using shot_number, which is unsafe if stored as double")
}
sel[, shot_key := trimws(shot_key)]
if (is.numeric(raw_date)) {
  date_origin <- if (date_was_delta) "2018-01-01 00:00:00" else "1970-01-01 00:00:00"
  sel[, acquisition_date := as.POSIXct(as.numeric(raw_date), origin = date_origin, tz = "UTC")]
} else {
  sel[, acquisition_date := as.POSIXct(raw_date, tz = "UTC")]
}
if (anyNA(sel$acquisition_date)) stop("Some selected shots have missing/unparseable acquisition dates")

old <- if (file.exists(cfg$output_file) && !cfg$overwrite) {
  as.data.table(read_parquet(cfg$output_file))
} else NULL
done_keys <- if (!is.null(old)) unique(old$shot_key) else character()
qa_old <- if (file.exists(cfg$qa_file) && !cfg$overwrite) fread(cfg$qa_file) else NULL

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

match_gedi_shot <- function(x, target) {
  # Exact character/integer64 comparison first. Some rGEDI versions return
  # shot_number as double; only then use a numeric fallback within the already
  # restricted granule/date/coordinate search.
  hit <- which(as.character(x) == target)
  if (length(hit)) return(hit)
  target_num <- suppressWarnings(as.numeric(target))
  x_num <- suppressWarnings(as.numeric(x))
  which(is.finite(x_num) & is.finite(target_num) & x_num == target_num)
}

extract_one_h5 <- function(h5_path, target_shot, ground_m, class_id, sample_id = NA_character_) {
  g <- readLevel1B(level1Bpath = h5_path)
  on.exit(try(close(g), silent = TRUE), add = TRUE)
  geo <- getLevel1BGeo(level1b = g, select = c("elevation_bin0"))
  idx <- match_gedi_shot(geo$shot_number, target_shot)
  if (!length(idx)) return(NULL)
  wf <- getLevel1BWF(g, geo$shot_number[idx[1]])
  ans <- data.table(
    architecture_class = as.integer(class_id),
    sample_id = as.character(sample_id),
    shot_number = geo$shot_number[idx[1]],
    shot_key = target_shot,
    waveform_bin = seq_along(wf@dt$rxwaveform),
    elevation_m = as.numeric(wf@dt$elevation),
    height_m = as.numeric(wf@dt$elevation) - as.numeric(ground_m),
    amplitude = as.numeric(wf@dt$rxwaveform),
    source_l1b = basename(h5_path)
  )
  ans[is.finite(height_m) & is.finite(amplitude)]
}

retrieve_l2a_ground <- function(z, daterange) {
  candidates <- tryCatch(
    gedifinder(cfg$l2a_product, z$latitude, z$longitude, z$latitude, z$longitude,
               version = cfg$l2a_version, daterange = daterange),
    error = function(e) character()
  )
  candidates <- unlist(candidates, use.names = FALSE)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (!length(candidates)) return(NA_real_)

  l2_dir <- l2_cache_dir
  for (remote in candidates) {
    before <- list.files(l2_dir, pattern = "\\.h5$", full.names = TRUE)
    try(gediDownload(filepath = remote, outdir = l2_dir), silent = TRUE)
    after <- list.files(l2_dir, pattern = "\\.h5$", full.names = TRUE)
    files <- unique(c(setdiff(after, before), after))
    for (hf in files) {
      ans <- tryCatch({
        g2 <- readLevel2A(level2Apath = hf)
        on.exit(try(close(g2), silent = TRUE), add = TRUE)
        m <- as.data.table(getLevel2AM(g2))
        idx <- match_gedi_shot(m$shot_number, z$shot_key)
        row <- if (length(idx)) m[idx[1]] else m[0]
        if (!nrow(row)) NA_real_ else {
          candidates_nm <- grep("^elev_lowestmode", names(row), value = TRUE)
          if (!length(candidates_nm)) candidates_nm <- grep("lowestmode", names(row), value = TRUE)
          if (!length(candidates_nm)) NA_real_ else {
            # Prefer the selected algorithm when the field is algorithm-specific.
            alg <- if ("selected_algorithm" %in% names(z)) as.integer(z$selected_algorithm) else NA_integer_
            preferred <- if (is.finite(alg)) candidates_nm[grepl(paste0("a", alg, "$"), candidates_nm)] else character()
            nm <- if (length(preferred)) preferred[1] else candidates_nm[1]
            as.numeric(row[[nm]][1])
          }
        }
      }, error = function(e) NA_real_)
      if (is.finite(ans)) {
        if (cfg$delete_h5_after_success) unlink(hf)
        return(ans)
      }
    }
  }
  NA_real_
}

new_wave <- list()
new_qa <- list()
for (ii in seq_len(nrow(sel))) {
  z <- sel[ii]
  if (!cfg$overwrite && z$shot_key %in% done_keys) next
  message(sprintf("[%d/%d] Class %s; shot %s", ii, nrow(sel), z$architecture_class, z$shot_key))
  d0 <- as.Date(z$acquisition_date) - cfg$date_buffer_days
  d1 <- as.Date(z$acquisition_date) + cfg$date_buffer_days
  daterange <- as.character(c(d0, d1))
  ground_m <- as.numeric(z$ground_elevation_m)
  if (!is.finite(ground_m)) {
    message("  Retrieving matching L2A ground elevation")
    ground_m <- retrieve_l2a_ground(z, daterange)
  }
  if (!is.finite(ground_m)) {
    new_qa[[length(new_qa) + 1L]] <- data.table(
      architecture_class = z$architecture_class,
      sample_id = if ("sample_id" %in% names(z)) as.character(z$sample_id) else NA_character_,
      shot_key = z$shot_key, status = "missing_ground_elevation",
      candidate_granules = NA_integer_, detail = "Matching GEDI02_A ground elevation not found"
    )
    qa_now <- rbindlist(c(list(qa_old), new_qa), use.names = TRUE, fill = TRUE)
    fwrite(qa_now, cfg$qa_file)
    next
  }
  status <- "not_found"
  detail <- NA_character_
  result <- NULL
  candidates <- tryCatch(
    gedifinder(cfg$product, z$latitude, z$longitude, z$latitude, z$longitude,
               version = cfg$version, daterange = daterange),
    error = function(e) { detail <<- conditionMessage(e); character() }
  )
  candidates <- unlist(candidates, use.names = FALSE)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  if (length(candidates)) {
    for (remote in candidates) {
      before <- list.files(cfg$cache_dir, pattern = "\\.h5$", full.names = TRUE)
      dl_error <- NULL
      tryCatch(gediDownload(filepath = remote, outdir = cfg$cache_dir),
               error = function(e) dl_error <<- conditionMessage(e))
      after <- list.files(cfg$cache_dir, pattern = "\\.h5$", full.names = TRUE)
      h5 <- setdiff(after, before)
      if (!length(h5)) h5 <- after
      if (!length(h5)) {
        detail <- dl_error %||% "Download produced no HDF5 file"
        next
      }
      for (hf in h5) {
        result <- tryCatch(
          extract_one_h5(hf, z$shot_key, ground_m,
                         z$architecture_class,
                         if ("sample_id" %in% names(z)) z$sample_id else NA_character_),
          error = function(e) { detail <<- conditionMessage(e); NULL }
        )
        if (!is.null(result) && nrow(result)) {
          status <- "success"
          if (cfg$delete_h5_after_success) unlink(hf)
          break
        }
      }
      if (status == "success") break
    }
  }
  if (!is.null(result) && nrow(result)) new_wave[[length(new_wave) + 1L]] <- result
  new_qa[[length(new_qa) + 1L]] <- data.table(
    architecture_class = z$architecture_class,
    sample_id = if ("sample_id" %in% names(z)) as.character(z$sample_id) else NA_character_,
    shot_key = z$shot_key,
    status = status,
    candidate_granules = length(candidates),
    detail = detail
  )

  # Checkpoint after every shot so interrupted downloads can resume.
  if (length(new_wave)) {
    combined <- rbindlist(c(list(old), new_wave), use.names = TRUE, fill = TRUE)
    write_parquet(combined, cfg$output_file, compression = "zstd")
  }
  qa_now <- rbindlist(c(list(qa_old), new_qa), use.names = TRUE, fill = TRUE)
  fwrite(qa_now, cfg$qa_file)
}

message("Waveform retrieval complete: ", cfg$output_file)
