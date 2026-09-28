#!/usr/bin/env Rscript

# STEP 02B — Download HLS Sentinel-2 point time series through AppEEARS in R.
# First-time credential setup (never put the password in this script):
#install.packages("appeears")
options(keyring_backend = "file")
appeears::rs_set_key(user = "monviso",
                       password = getPass::getPass("Earthdata password:"))

options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  earthdata_user = "monviso",
  # The encrypted file backend is more reliable than Windows Credential
  # Manager for scripted/batch access. The same backend must be used when
  # saving and retrieving the Earthdata password.
  keyring_backend = "file",
  sample_csv = "E:/Xingu_rev/data/step02_radiometry_v3/gedi_radiometric_sample.csv",
  output_dir = "E:/Xingu_rev/data/step02_radiometry_v3/appeears_downloads",
  batch_size = 200L,
  requested_layers = c("B02", "B03", "B04", "B05", "B8A", "B11", "B12",
                       "Fmask"),
  overwrite_completed_tasks = FALSE
)

options(keyring_backend = cfg$keyring_backend)

required <- c("appeears", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(appeears); library(data.table)})

# AppEEARS catalogue tables can contain nested/list columns (including NULL).
# Keep the in-memory catalogue unchanged, but flatten a copy for audit CSVs.
flatten_catalogue_for_csv <- function(x) {
  out <- as.data.table(copy(x))
  list_columns <- names(out)[vapply(out, is.list, logical(1))]
  for (column in list_columns) {
    out[, (column) := vapply(get(column), function(value) {
      if (is.null(value) || !length(value)) return(NA_character_)
      paste(as.character(unlist(value, recursive = TRUE,
                                use.names = FALSE)), collapse = " | ")
    }, character(1))]
  }
  out
}

if (cfg$earthdata_user == "REPLACE_WITH_EARTHDATA_USERNAME") {
  stop("Set cfg$earthdata_user before running")
}

# Fail once, before submitting any batches, if the credential is absent or was
# saved using another keyring backend or username.
credential_check <- tryCatch(
  appeears::rs_get_key(user = cfg$earthdata_user),
  error = function(e) e
)
if (inherits(credential_check, "error")) {
  stop(
    "No readable AppEEARS credential for user '", cfg$earthdata_user, "'.\n",
    "In a fresh interactive R session run:\n",
    "  options(keyring_backend = \"", cfg$keyring_backend, "\")\n",
    "  appeears::rs_set_key(user = \"", cfg$earthdata_user,
    "\", password = getPass::getPass(\"Earthdata password: \"))\n",
    "  appeears::rs_get_key(user = \"", cfg$earthdata_user, "\")\n",
    "Then rerun this script. Original error: ", conditionMessage(credential_check)
  )
}
if (!file.exists(cfg$sample_csv)) stop("Sample not found: ", cfg$sample_csv)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

# Discover the current HLSS30 version and layer names instead of hard-coding
# identifiers that AppEEARS may update.
products <- as.data.table(appeears::rs_products())
if (!"ProductAndVersion" %in% names(products)) stop("Unexpected product catalogue")
hls_products <- products[grepl("^HLSS30\\.", ProductAndVersion)]
fwrite(flatten_catalogue_for_csv(hls_products),
       file.path(cfg$output_dir, "available_HLSS30_products.csv"))
if (!nrow(hls_products)) stop("HLSS30 not found in the AppEEARS catalogue")
version_score <- vapply(strsplit(hls_products$ProductAndVersion, "\\."),
                        function(x) suppressWarnings(as.numeric(tail(x, 1L))),
                        numeric(1))
version_score[!is.finite(version_score)] <- 0
product <- hls_products$ProductAndVersion[which.max(version_score)]
message("Using AppEEARS product: ", product)

layers <- as.data.table(appeears::rs_layers(product))
fwrite(flatten_catalogue_for_csv(layers),
       file.path(cfg$output_dir, "available_HLSS30_layers.csv"))
layer_field <- intersect(c("Layer", "layer", "Band", "band"), names(layers))[1L]
if (is.na(layer_field)) stop("Inspect available_HLSS30_layers.csv: layer column not found")
available_layers <- as.character(layers[[layer_field]])
match_layer <- function(target) {
  hit <- available_layers[toupper(available_layers) == toupper(target)]
  if (length(hit)) return(hit[1L])
  hit <- available_layers[grepl(paste0("(^|_)", target, "$"),
                                available_layers, ignore.case = TRUE)]
  if (length(hit)) hit[1L] else NA_character_
}
selected_layers <- vapply(cfg$requested_layers, match_layer, character(1))
if (anyNA(selected_layers)) {
  stop("Missing layers: ", paste(cfg$requested_layers[is.na(selected_layers)],
                                 collapse = ", "),
       "\nInspect available_HLSS30_layers.csv")
}

sample_dt <- fread(cfg$sample_csv)
needed <- c("sample_id", "lat_lowestmode", "lon_lowestmode", "year")
absent <- setdiff(needed, names(sample_dt))
if (length(absent)) stop("Sample missing: ", paste(absent, collapse = ", "))
setorder(sample_dt, year, sample_id)

manifest_file <- file.path(cfg$output_dir, "appeears_task_manifest.csv")
manifest <- if (file.exists(manifest_file)) fread(manifest_file) else
  data.table(task = character(), year = integer(), batch = integer(),
             n_points = integer(), status = character(), message = character())

for (yr in sort(unique(sample_dt$year))) {
  year_dt <- sample_dt[year == yr]
  year_dt[, batch := ceiling(seq_len(.N) / cfg$batch_size)]
  for (batch_id in sort(unique(year_dt$batch))) {
    pts <- year_dt[batch == batch_id]
    task_name <- sprintf("Xingu_HLSS30_%d_b%03d", yr, batch_id)
    if (nrow(manifest[task == task_name & status == "completed"]) &&
        !cfg$overwrite_completed_tasks) {
      message("Skipping completed task: ", task_name); next
    }
    request_dt <- CJ(point_row = seq_len(nrow(pts)),
                     layer_row = seq_along(selected_layers))
    request_dt[, `:=`(
      task = task_name,
      subtask = pts$sample_id[point_row],
      latitude = pts$lat_lowestmode[point_row],
      longitude = pts$lon_lowestmode[point_row],
      start = sprintf("%d-01-01", yr), end = sprintf("%d-12-31", yr),
      product = product, layer = selected_layers[layer_row]
    )]
    request_df <- as.data.frame(request_dt[, .(
      task, subtask, latitude, longitude, start, end, product, layer
    )])
    message("Submitting ", task_name, "; points=", nrow(pts))
    outcome <- tryCatch({
      task_json <- appeears::rs_build_task(df = request_df)
      appeears::rs_request(request = task_json, user = cfg$earthdata_user,
                           transfer = TRUE, path = cfg$output_dir, verbose = TRUE)
      list(status = "completed", message = "downloaded")
    }, error = function(e) list(status = "failed", message = conditionMessage(e)))
    manifest <- manifest[task != task_name]
    manifest <- rbind(manifest, data.table(
      task = task_name, year = yr, batch = batch_id, n_points = nrow(pts),
      status = outcome$status, message = outcome$message
    ), fill = TRUE)
    fwrite(manifest, manifest_file)
    if (outcome$status == "failed") warning(task_name, ": ", outcome$message)
  }
}
message("STEP 02B COMPLETE: ", manifest_file)

