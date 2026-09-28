#!/usr/bin/env Rscript

# STEP 06 — Find classified GEDI footprints that overlap available airborne
# LiDAR (ALS) coverage polygons and export the exact Monsoon LAZ paths.
#
# The 13-million-shot assignment dataset is never converted to sf in full.
# ALS polygons are grouped spatially; Arrow reads only the point coordinates
# inside each small group bounding box before exact spatial intersection.

options(stringsAsFactors = FALSE, warn = 1)
cfg <- list(
  project_dir = "E:/Xingu_rev",
  assignment_dataset = "E:/Xingu_rev/data/step01b_final_classification/assignments_dataset",
  transect_dir = "E:/temp_download_toGLOBUS",
  transect_files = c(
    "FIB_2018_ref_mons.shp",
    "Ometto_2018_ref_mons.shp",
    "paisagens_2021_ref_mons.shp"
  ),
  output_dir = "E:/Xingu_rev/data/step06_GEDI_ALS_overlap_v3",
  # NULL retains every class selected by Step 01a/01b.
  classes_to_keep = NULL,
  # Centre-in-polygon is the strict/default definition. Set to 12.5 to retain
  # a GEDI shot whose nominal 25-m footprint intersects an ALS polygon.
  footprint_radius_m = 0,
  equal_area_crs = 6933,
  polygon_group_size_degrees = 0.25,
  # Keep NULL initially. Later, e.g. set 3L to retain |GEDI year - ALS year|<=3.
  maximum_year_difference = NULL,
  overwrite = FALSE
)

required <- c("arrow", "data.table", "dplyr", "sf")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(dplyr); library(sf)
})
sf::sf_use_s2(TRUE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

matched_csv <- file.path(cfg$output_dir, "GEDI_ALS_matches.csv")
matched_gpkg <- file.path(cfg$output_dir, "GEDI_ALS_matches.gpkg")
summary_csv <- file.path(cfg$output_dir, "GEDI_ALS_match_summary.csv")
manifest_csv <- file.path(cfg$output_dir, "ALS_transfer_manifest.csv")
paths_txt <- file.path(cfg$output_dir, "monsoon_laz_paths.txt")
for (path in c(matched_csv, matched_gpkg, summary_csv, manifest_csv, paths_txt)) {
  if (file.exists(path) && !cfg$overwrite) {
    stop("Output already exists: ", path,
         "\nSet cfg$overwrite=TRUE only to intentionally rebuild all outputs.")
  }
}

read_als_layer <- function(filename) {
  path <- file.path(cfg$transect_dir, filename)
  if (!file.exists(path)) stop("Missing ALS shapefile: ", path)
  layer_name <- tools::file_path_sans_ext(basename(filename))
  x <- st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
  if (is.na(st_crs(x))) stop("Missing CRS: ", path)
  x <- st_transform(x, 4326)
  x <- st_make_valid(x)
  x <- x[!st_is_empty(x), ]
  if (!"path2" %in% names(x)) stop("Column path2 is missing from ", filename)
  dataset_value <- if ("dataset" %in% names(x)) as.character(x$dataset) else layer_name
  dataset_value[is.na(dataset_value) | !nzchar(trimws(dataset_value))] <- layer_name
  als_year <- if ("year" %in% names(x)) suppressWarnings(as.integer(x$year)) else NA_integer_
  reference_id <- if ("i_ref" %in% names(x)) as.character(x$i_ref) else as.character(seq_len(nrow(x)))
  out <- st_sf(
    als_polygon_id = paste0(layer_name, "_", seq_len(nrow(x))),
    als_dataset = trimws(dataset_value),
    als_year = als_year,
    als_reference_id = reference_id,
    path2 = trimws(as.character(x$path2)),
    geometry = st_geometry(x), crs = 4326
  )
  out[!is.na(out$path2) & nzchar(out$path2), ]
}

als <- do.call(rbind, lapply(cfg$transect_files, read_als_layer))
if (!nrow(als)) stop("No usable ALS polygons were read")
if (any(!st_geometry_type(als) %in% c("POLYGON", "MULTIPOLYGON"))) {
  stop("Every ALS footprint must be a polygon or multipolygon")
}

# Assign each ALS polygon once to a compact processing group. Exact polygon
# bounding boxes are still used for every Arrow query.
centres <- suppressWarnings(st_coordinates(st_centroid(st_geometry(als))))
als$group_x <- floor(centres[, 1L] / cfg$polygon_group_size_degrees)
als$group_y <- floor(centres[, 2L] / cfg$polygon_group_size_degrees)
als$spatial_group <- paste(als$group_x, als$group_y, sep = "_")
groups <- unique(als$spatial_group)
message("ALS polygons: ", nrow(als), "; spatial processing groups: ", length(groups))

ds <- open_dataset(cfg$assignment_dataset, format = "parquet")
minimum_columns <- c("shot_number", "latitude", "longitude",
                     "year", "architecture_class")
absent <- setdiff(minimum_columns, names(ds))
if (length(absent)) stop("Assignment dataset is missing: ", paste(absent, collapse = ", "))
optional_columns <- c("sample_id", "sensitivity", "solar_elevation",
                      "rh98_original", "mapbiomas_class", "degrade_flag",
                      "quality_flag")
columns_to_read <- intersect(c(minimum_columns, optional_columns), names(ds))

matches <- vector("list", length(groups))
for (g in seq_along(groups)) {
  group_name <- groups[g]
  polygons <- als[als$spatial_group == group_name, ]
  box <- st_bbox(polygons)
  message(sprintf("[%d/%d] group=%s; polygons=%d", g, length(groups),
                  group_name, nrow(polygons)))
  
  query <- ds %>%
    filter(longitude >= !!unname(box["xmin"]),
           longitude <= !!unname(box["xmax"]),
           latitude >= !!unname(box["ymin"]),
           latitude <= !!unname(box["ymax"]))
  if (!is.null(cfg$classes_to_keep)) {
    query <- query %>% filter(architecture_class %in% cfg$classes_to_keep)
  }
  candidates <- query %>% select(all_of(columns_to_read)) %>% collect() %>%
    as.data.table()
  if (!nrow(candidates)) next
  
  points <- st_as_sf(candidates,
                     coords = c("longitude", "latitude"), crs = 4326,
                     remove = FALSE)
  if (cfg$footprint_radius_m > 0) {
    point_metric <- st_transform(points, cfg$equal_area_crs)
    polygon_metric <- st_transform(polygons, cfg$equal_area_crs)
    hit_list <- st_is_within_distance(point_metric, polygon_metric,
                                      dist = cfg$footprint_radius_m)
  } else {
    hit_list <- st_intersects(points, polygons)
  }
  point_index <- rep(seq_along(hit_list), lengths(hit_list))
  if (!length(point_index)) next
  polygon_index <- unlist(hit_list, use.names = FALSE)
  joined <- cbind(
    as.data.table(st_drop_geometry(points))[point_index],
    as.data.table(st_drop_geometry(polygons))[
      polygon_index,
      .(als_polygon_id, als_dataset, als_year, als_reference_id, path2)]
  )
  joined[, year_difference := as.integer(year) - als_year]
  joined[, absolute_year_difference := abs(year_difference)]
  if (!is.null(cfg$maximum_year_difference)) {
    joined <- joined[is.na(absolute_year_difference) |
                       absolute_year_difference <= cfg$maximum_year_difference]
  }
  matches[[g]] <- joined
  if (g %% 10L == 0L) gc(FALSE)
}

matched <- rbindlist(matches, use.names = TRUE, fill = TRUE)
if (!nrow(matched)) {
  stop("No classified GEDI shots intersect the supplied ALS polygons under the current filters")
}
matched[, shot_number := as.character(shot_number)]
matched <- unique(matched, by = c("shot_number", "als_polygon_id", "path2"))
setorder(matched, architecture_class, absolute_year_difference,
         als_dataset, path2, shot_number)
fwrite(matched, matched_csv)

matched_points <- st_as_sf(matched,
                           coords = c("longitude", "latitude"), crs = 4326,
                           remove = FALSE)
if (file.exists(matched_gpkg)) file.remove(matched_gpkg)
st_write(matched_points, matched_gpkg, layer = "GEDI_ALS_matches", quiet = TRUE)

summary_dt <- matched[, .(
  matched_pairs = .N,
  unique_gedi_shots = uniqueN(shot_number),
  unique_laz_files = uniqueN(path2),
  median_absolute_year_difference = as.numeric(median(
    absolute_year_difference,
    na.rm = TRUE
  ))
), by = .(architecture_class, als_dataset, als_year, gedi_year = year)]
setorder(summary_dt, architecture_class, als_dataset, als_year, gedi_year)
fwrite(summary_dt, summary_csv)

manifest <- unique(matched[, .(
  als_dataset, als_year, monsoon_path = path2,
  laz_filename = basename(path2)
)])
path_counts <- matched[, .(matched_gedi_shots = uniqueN(shot_number)), by = path2]
manifest <- merge(manifest, path_counts, by.x = "monsoon_path", by.y = "path2",
                  all.x = TRUE)
setorder(manifest, als_dataset, monsoon_path)
fwrite(manifest, manifest_csv)
writeLines(manifest$monsoon_path, paths_txt, useBytes = TRUE)

message("STEP 06 COMPLETE")
message("Unique GEDI shots: ", uniqueN(matched$shot_number))
message("Unique LAZ files to transfer: ", nrow(manifest))
message("Matches: ", matched_csv)
message("Transfer paths: ", paths_txt)

