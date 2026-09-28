#!/usr/bin/env Rscript

# Select a small, structurally representative set of GEDI shots per final
# architecture class for L1B waveform retrieval.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

cfg <- list(
  # This contains the shots selected for NMF and their FINAL class labels, but
  # only RH0/RH98. Complete RH profiles are joined from source_files below.
  class_file = "E:/Xingu_rev/data/step01b_final_classification/nmf_training_sample.parquet",
  source_files = file.path(
    "E:/Xingu_rev/data",
    sprintf("GEDI_RFLD_sampling_frame_%d.parquet", 2019:2023)
  ),
  output_file = "E:/Xingu_rev/data/step01b_final_classification/figures/representative_waveform_shots.csv",
  # Independent waveform sample size. Use 20L-50L for normal runs or "ALL".
  n_per_class = 30L,
  seed = 42L
)

read_any <- function(path) {
  if (dir.exists(path)) return(as.data.table(collect(open_dataset(path))))
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) return(as.data.table(read_parquet(path)))
  fread(path)
}

first_name <- function(nms, candidates, required = TRUE) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) return(hit[1])
  if (required) stop("None of these columns was found: ", paste(candidates, collapse = ", "))
  NULL
}

labels <- read_any(cfg$class_file)
class_col <- first_name(names(labels), c("architecture_class", "class", "final_class"))
shot_col  <- first_name(names(labels), c("shot_number", "shot_num", "shot"))
year_col  <- first_name(names(labels), c("year", "gedi_year"), FALSE)
setnames(labels, c(class_col, shot_col), c("architecture_class", "shot_number"))
if (!is.null(year_col) && year_col != "year") setnames(labels, year_col, "year")
label_fields <- intersect(c("sample_id", "shot_number", "year", "architecture_class"), names(labels))
labels <- unique(labels[, ..label_fields])
# sample_id is the workflow's exact unique identifier and avoids possible
# precision/collision problems with large GEDI shot numbers.
join_keys <- if ("sample_id" %in% names(labels)) "sample_id" else intersect(c("shot_number", "year"), names(labels))
message("Class-label uniqueness key: ", paste(join_keys, collapse = " + "))
label_conflicts <- labels[, .(n_classes = uniqueN(architecture_class)), by = join_keys][n_classes > 1L]
if (nrow(label_conflicts)) {
  stop("Class file contains ", nrow(label_conflicts),
       " shot/year keys assigned to more than one architecture class")
}
label_duplicates <- nrow(labels) - uniqueN(labels, by = join_keys)
if (label_duplicates > 0L) {
  message("Removing ", label_duplicates, " duplicate class-label key(s)")
  labels <- unique(labels, by = join_keys)
}

if (!length(cfg$source_files) || any(!file.exists(cfg$source_files))) {
  stop("Missing source sampling-frame file(s): ",
       paste(cfg$source_files[!file.exists(cfg$source_files)], collapse = ", "))
}

# Read each annual source independently, retaining only the labelled NMF shots.
# This avoids loading all five full RH tables simultaneously.
joined <- vector("list", length(cfg$source_files))
for (fi in seq_along(cfg$source_files)) {
  f <- cfg$source_files[fi]
  src_names <- names(open_dataset(f, format = "parquet"))
  src_shot <- first_name(src_names, c("shot_number", "shot_num", "shot"))
  src_year <- first_name(src_names, c("year", "gedi_year"), FALSE)
  rh_src <- vapply(0:98, function(q) first_name(
    src_names,
    c(sprintf("rh%d", q), sprintf("rh_%d", q),
      sprintf("rh%d_original", q), sprintf("rh_%d_original", q),
      sprintf("rh_%d_a0", q), sprintf("rh_%d_ln", q))
  ), character(1))
  optional <- intersect(src_names, c(
    "longitude", "latitude", "lon_lowestmode", "lat_lowestmode", "lon_lm", "lat_lm",
    "date", "acquisition_date", "acquisition_datetime", "acquisition_time", "datetime", "delta_time",
    "beam", "beam_name",
    "granule", "granule_id", "source_granule", "source_file", "filename",
    "elevation_lowestmode", "elev_lowestmode", "elev_lm", "elev_lm_a0",
    "elevation_highestreturn", "elev_highestreturn", "elev_hr", "elev_hr_a0",
    "pai", "pai_a0", "pai_z", "plant_area_index",
    "agbd", "agbd_a0", "aboveground_biomass_density"
    ,"selected_algorithm", "selected_mode"
  ))
  src_sample <- first_name(src_names, c("sample_id"), FALSE)
  wanted <- unique(c(src_sample, src_shot, src_year, rh_src, optional))
  src <- as.data.table(read_parquet(f, col_select = all_of(wanted)))
  if (src_shot != "shot_number") setnames(src, src_shot, "shot_number")
  if (!is.null(src_year) && src_year != "year") setnames(src, src_year, "year")
  if (!("year" %in% names(src))) {
    yr <- as.integer(sub(".*(20[0-9]{2}).*", "\\1", basename(f)))
    src[, year := yr]
  }
  lab_y <- if ("year" %in% names(labels)) labels[year == unique(src$year)[1]] else labels
  keys <- if ("sample_id" %in% names(src) && "sample_id" %in% names(lab_y)) {
    "sample_id"
  } else {
    intersect(c("shot_number", "year"), names(lab_y))
  }
  if (!length(keys)) stop("No common exact join key for ", basename(f))
  # When sample_id is the key, carry only the class label from the label table.
  # Otherwise merge() would create shot_number.x/shot_number.y and year.x/year.y.
  lab_join <- lab_y[, unique(c(keys, "architecture_class")), with = FALSE]
  n_before <- nrow(src)
  src <- unique(src, by = keys)
  if (nrow(src) < n_before) {
    message("Collapsed ", n_before - nrow(src), " duplicate source key(s) in ", basename(f))
  }
  joined[[fi]] <- merge(src, lab_join, by = keys,
                        all = FALSE, sort = FALSE)
  message("Joined ", nrow(joined[[fi]]), " labelled profiles from ", basename(f))
}
d <- rbindlist(joined, use.names = TRUE, fill = TRUE)
if (!nrow(d)) stop("The class-label/source-profile join returned zero rows")

lon_col   <- first_name(names(d), c("longitude", "lon_lowestmode", "lon_lm", "lon"), FALSE)
lat_col   <- first_name(names(d), c("latitude", "lat_lowestmode", "lat_lm", "lat"), FALSE)

rh_name <- function(q) first_name(
  names(d),
  c(sprintf("rh%d", q), sprintf("rh_%d", q), sprintf("rh%d_original", q),
    sprintf("rh_%d_original", q), sprintf("rh_%d_a0", q), sprintf("rh_%d_ln", q))
)
rh_cols <- vapply(c(25L, 50L, 75L, 98L), rh_name, character(1))

setnames(d, rh_cols, c("rh25", "rh50", "rh75", "rh98"))
if (!is.null(lon_col) && lon_col != "longitude") setnames(d, lon_col, "longitude")
if (!is.null(lat_col) && lat_col != "latitude") setnames(d, lat_col, "latitude")

keep <- complete.cases(d[, .(architecture_class, shot_number, rh25, rh50, rh75, rh98)])
d <- d[keep]
set.seed(cfg$seed)

# Robust distance from the within-class median profile. This represents the
# complete RH profile better than selecting on RH98 alone.
selected <- d[, {
  z <- copy(.SD)
  for (nm in c("rh25", "rh50", "rh75", "rh98")) {
    med <- median(z[[nm]], na.rm = TRUE)
    sc <- mad(z[[nm]], center = med, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sc) || sc == 0) sc <- sd(z[[nm]], na.rm = TRUE)
    if (!is.finite(sc) || sc == 0) sc <- 1
    z[, paste0("z_", nm) := (get(nm) - med) / sc]
  }
  z[, profile_distance := sqrt(z_rh25^2 + z_rh50^2 + z_rh75^2 + z_rh98^2)]
  z[, tie_break := runif(.N)]
  setorder(z, profile_distance, tie_break)
  take_n <- if (is.character(cfg$n_per_class) && toupper(cfg$n_per_class) == "ALL") {
    .N
  } else {
    min(.N, as.integer(cfg$n_per_class))
  }
  z[seq_len(take_n)]
}, by = architecture_class]

selected[, selection_rank := seq_len(.N), by = architecture_class]
selected[, c("z_rh25", "z_rh50", "z_rh75", "z_rh98", "tie_break") := NULL]
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(selected, cfg$output_file)
message("Wrote ", nrow(selected), " representative shots to: ", cfg$output_file)
