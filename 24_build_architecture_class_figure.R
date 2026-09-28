#!/usr/bin/env Rscript

# Build the complete architecture-class synthesis figure:
#   mean L1B waveform + mean RH curve | spatial distribution
#   class structural summary table
# The script also exports every class panel and all plotted data.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

cfg <- list(
  # Final NMF training-shot IDs and class labels. This file contains only RH0
  # and RH98, so complete profiles are joined from source_files.
  class_file = "E:/Xingu_rev/data/step01b_final_classification/nmf_training_sample.parquet",
  source_files = file.path(
    "E:/Xingu_rev/data",
    sprintf("GEDI_RFLD_sampling_frame_%d.parquet", 2019:2023)
  ),
  
  # Points used for the spatial maps. The training sample is recommended here:
  # it is small enough for local plotting and is the dataset requested for the
  # figure. The all-shot assignment dataset may contain millions of records.
  map_file = "E:/Xingu_rev/data/step01b_final_classification/nmf_training_sample.parquet",
  
  # Basin polygon. Shapefile or GeoPackage are both accepted.
  boundary_file = "E:/Xingu_rev/data/boundaries/Xingu_basin.shp",
  
  # Tidy file produced after L1B retrieval. Required columns (aliases allowed):
  # architecture_class, shot_number, height_m, amplitude.
  # CSV or Parquet; set to NA_character_ to draw profiles without waveforms.
  waveform_file = "E:/Xingu_rev/data/step01b_final_classification/figures/L1B_waveforms_tidy.parquet",
  
  output_dir = "E:/Xingu_rev/data/step01b_final_classification/figures/architecture_synthesis",
  # These three samples are independent. Each accepts an integer or "ALL".
  profile_points_per_class = "ALL",  # mean RH curves and their ribbons
  map_points_per_class = 5000L,      # display points only; avoids solid maps
  summary_points_per_class = "ALL", # RH/PAI/AGBD table
  waveform_grid_step_m = 0.25,
  # Shared profile-axis limits. The default clips implausibly deep waveform
  # tails at -10 m while deriving the upper limit automatically. Use
  # c(NA_real_, NA_real_) for fully automatic limits or c(-15, 50), etc.
  profile_y_limits_m = c(-10, 35),
  # Identify one continuous waveform display interval from a smoothed mean
  # waveform. This removes long near-zero tails without creating gaps/dashes
  # in the original (unsmoothed) waveform curve. The threshold is a fraction
  # of the class-specific peak; increase to 0.07--0.10 for stronger trimming.
  waveform_min_display_fraction = 0.05,
  waveform_support_smooth_m = 1.0,
  waveform_support_padding_m = 0.5,
  seed = 42L,
  figure_columns = 3L,
  # Export the assembled paper figure at its actual printed dimensions. Do
  # not export a much larger figure and resize it afterwards, because that
  # also shrinks every text element.
  final_width_mm = 160,
  final_height_mm = 120,
  # Larger standalone panels are retained for inspection/supplementary use.
  single_panel_width_mm = 105,
  single_panel_height_mm = 72,
  dpi = 600,
  # Text controls. Theme text is measured in points; geom_text is measured in
  # millimetres. Keeping the two groups explicit avoids unit confusion.
  text_size = list(
    class_title_pt = 7.0,
    profile_axis_title_pt = 4.8,
    profile_axis_text_pt = 4.2,
    map_axis_text_pt = 3.4,
    rh_label_mm = 1.25,
    table_text_mm = 4.45
  ),
  geometry = list(
    map_point_mm = 0.001,
    map_point_alpha = 0.20,
    map_boundary_line_mm = 0.22,
    waveform_line_mm = 0.42,
    rh_line_mm = 0.38,
    reference_line_mm = 0.20,
    rh_segment_line_mm = 0.22,
    rh_arrow_mm = 0.85,
    waveform_ribbon_alpha = 0.28,
    rh_ribbon_alpha = 0.22
  ),
  class_colours = c(
    "#E9C411", "#EF7900", "#8A1F62", "#43205F",
    "#35A2A5", "#8DDE91", "#2374AB", "#B45F06", "#666666"
  )
)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

read_any <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path) && !dir.exists(path)) return(NULL)
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

standardize_id_columns <- function(d) {
  cc <- first_name(names(d), c("architecture_class", "class", "final_class"))
  if (cc != "architecture_class") setnames(d, cc, "architecture_class")
  sc <- first_name(names(d), c("shot_number", "shot_num", "shot"), FALSE)
  if (!is.null(sc) && sc != "shot_number") setnames(d, sc, "shot_number")
  d[, architecture_class := as.integer(gsub("[^0-9]", "", as.character(architecture_class)))]
  d
}

find_rh_columns <- function(nms) {
  out <- rep(NA_character_, 99L)
  names(out) <- as.character(0:98)
  for (q in 0:98) {
    cand <- c(sprintf("rh%d", q), sprintf("rh_%d", q),
              sprintf("rh%d_original", q), sprintf("rh_%d_original", q),
              sprintf("rh_%d_a0", q), sprintf("rh_%d_ln", q))
    hit <- cand[cand %in% nms]
    if (length(hit)) out[as.character(q)] <- hit[1]
  }
  if (anyNA(out)) stop("Missing RH columns: ", paste(names(out)[is.na(out)], collapse = ", "))
  out
}

# Join final labels back to the complete annual RH source tables.
labels <- standardize_id_columns(read_any(cfg$class_file))
if (!("shot_number" %in% names(labels))) stop("Class file has no shot_number")
label_year <- first_name(names(labels), c("year", "gedi_year"), FALSE)
if (!is.null(label_year) && label_year != "year") setnames(labels, label_year, "year")
label_fields <- intersect(c("sample_id", "shot_number", "year", "architecture_class"), names(labels))
labels <- unique(labels[, ..label_fields])
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
if (any(!file.exists(cfg$source_files))) {
  stop("Missing source sampling-frame file(s): ",
       paste(cfg$source_files[!file.exists(cfg$source_files)], collapse = ", "))
}

profile_parts <- vector("list", length(cfg$source_files))
for (fi in seq_along(cfg$source_files)) {
  f <- cfg$source_files[fi]
  src_names <- names(open_dataset(f, format = "parquet"))
  src_shot <- first_name(src_names, c("shot_number", "shot_num", "shot"))
  src_year <- first_name(src_names, c("year", "gedi_year"), FALSE)
  source_rh <- find_rh_columns(src_names)
  optional <- intersect(src_names, c(
    "longitude", "latitude", "lon_lowestmode", "lat_lowestmode", "lon_lm", "lat_lm",
    "pai", "pai_a0", "pai_z", "plant_area_index",
    "agbd", "agbd_a0", "aboveground_biomass_density"
  ))
  src_sample <- first_name(src_names, c("sample_id"), FALSE)
  wanted <- unique(c(src_sample, src_shot, src_year, source_rh, optional))
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
  lab_join <- lab_y[, unique(c(keys, "architecture_class")), with = FALSE]
  n_before <- nrow(src)
  src <- unique(src, by = keys)
  if (nrow(src) < n_before) {
    message("Collapsed ", n_before - nrow(src), " duplicate source key(s) in ", basename(f))
  }
  profile_parts[[fi]] <- merge(
    src, lab_join,
    by = keys,
    all = FALSE, sort = FALSE
  )
  message("Joined ", nrow(profile_parts[[fi]]), " labelled profiles from ", basename(f))
}
profile <- rbindlist(profile_parts, use.names = TRUE, fill = TRUE)
if (!nrow(profile)) stop("The class-label/source-profile join returned zero rows")
rh_cols <- find_rh_columns(names(profile))
class_ids <- sort(unique(profile$architecture_class))
n_class <- length(class_ids)
if (!n_class) stop("No architecture classes found")
class_cols <- setNames(rep(cfg$class_colours, length.out = n_class), class_ids)

sample_by_class <- function(x, limit, seed) {
  if (is.character(limit) && toupper(limit) == "ALL") return(copy(x))
  n_take <- as.integer(limit)
  if (!is.finite(n_take) || n_take < 1L) stop("Sampling limits must be positive integers or 'ALL'")
  set.seed(seed)
  x[, .SD[sample.int(.N, min(.N, n_take))], by = architecture_class]
}

profile_curve <- sample_by_class(profile, cfg$profile_points_per_class, cfg$seed + 1L)
profile_metrics <- sample_by_class(profile, cfg$summary_points_per_class, cfg$seed + 2L)

# Long RH data and class summaries.
rh_dt <- melt(
  profile_curve[, c("architecture_class", rh_cols), with = FALSE],
  id.vars = "architecture_class", variable.name = "rh_source", value.name = "height_m"
)
rh_lookup <- setNames(as.integer(names(rh_cols)), rh_cols)
rh_dt[, rh := rh_lookup[as.character(rh_source)]]
rh_dt <- rh_dt[is.finite(height_m)]

rh_summary <- rh_dt[, .(
  mean_height_m = mean(height_m),
  sd_height_m = sd(height_m),
  median_height_m = median(height_m),
  q25_height_m = quantile(height_m, 0.25),
  q75_height_m = quantile(height_m, 0.75),
  n = .N
), by = .(architecture_class, rh)]

# Dynamic metric summary: RH rows always; PAI and AGBD added when available.
metric_specs <- list(
  list(label = "Understory - RH25 (m)", candidates = c("rh25", "rh_25", "rh25_original", "rh_25_original", "rh_25_a0", "rh_25_ln")),
  list(label = "Midstory - RH50 (m)", candidates = c("rh50", "rh_50", "rh50_original", "rh_50_original", "rh_50_a0", "rh_50_ln")),
  list(label = "Upper canopy - RH75 (m)", candidates = c("rh75", "rh_75", "rh75_original", "rh_75_original", "rh_75_a0", "rh_75_ln")),
  list(label = "Maximum height - RH98 (m)", candidates = c("rh98", "rh_98", "rh98_original", "rh_98_original", "rh_98_a0", "rh_98_ln")),
  list(label = "PAI", candidates = c("pai", "pai_a0", "pai_z", "plant_area_index")),
  list(label = "AGBD (Mg ha⁻¹)", candidates = c("agbd", "agbd_a0", "aboveground_biomass_density"))
)

metric_rows <- list()
for (sp in metric_specs) {
  nm <- first_name(names(profile_metrics), sp$candidates, FALSE)
  if (is.null(nm)) next
  z <- profile_metrics[is.finite(get(nm)), .(
    mean = mean(get(nm)), sd = sd(get(nm)), n = .N
  ), by = architecture_class]
  z[, `:=`(metric = sp$label, display = sprintf("%.1f ± %.1f", mean, sd))]
  metric_rows[[length(metric_rows) + 1L]] <- z
}
metric_summary <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
metric_summary[, metric_order := match(metric, vapply(metric_specs, `[[`, character(1), "label"))]
setorder(metric_summary, architecture_class, metric_order)
fwrite(metric_summary, file.path(cfg$output_dir, "architecture_class_metric_summary.csv"))
fwrite(rh_summary, file.path(cfg$output_dir, "architecture_class_RH_profiles.csv"))

# Spatial sample and boundary.
set.seed(cfg$seed)
map_dt <- standardize_id_columns(read_any(cfg$map_file))
lon <- first_name(names(map_dt), c("longitude", "lon_lowestmode", "lon_lm", "lon"), FALSE)
lat <- first_name(names(map_dt), c("latitude", "lat_lowestmode", "lat_lm", "lat"), FALSE)
if (is.null(lon) || is.null(lat)) {
  message("Map file lacks coordinates; using the joined complete-profile table")
  map_dt <- copy(profile)
  lon <- first_name(names(map_dt), c("longitude", "lon_lowestmode", "lon_lm", "lon"))
  lat <- first_name(names(map_dt), c("latitude", "lat_lowestmode", "lat_lm", "lat"))
}
if (lon != "longitude") setnames(map_dt, lon, "longitude")
if (lat != "latitude") setnames(map_dt, lat, "latitude")
map_dt <- map_dt[is.finite(longitude) & is.finite(latitude)]
map_dt <- sample_by_class(map_dt, cfg$map_points_per_class, cfg$seed + 3L)
boundary <- st_read(cfg$boundary_file, quiet = TRUE)
boundary <- st_transform(boundary, 4326)

# Waveform aggregation on a common height grid.
wave <- read_any(cfg$waveform_file)
wave_summary <- NULL
if (!is.null(wave)) {
  wave <- standardize_id_columns(wave)
  hcol <- first_name(names(wave), c("height_m", "height", "relative_height_m", "elev"))
  acol <- first_name(names(wave), c("amplitude", "rxwaveform", "wf", "waveform_amplitude"))
  if (hcol != "height_m") setnames(wave, hcol, "height_m")
  if (acol != "amplitude") setnames(wave, acol, "amplitude")
  if (!("shot_number" %in% names(wave))) {
    idc <- first_name(names(wave), c("id", "waveform_id"))
    setnames(wave, idc, "shot_number")
  }
  wave <- wave[is.finite(height_m) & is.finite(amplitude)]
  # Normalize amplitude within shot: L1B amplitude units vary among shots.
  wave[, amplitude_norm := {
    a <- amplitude - min(amplitude, na.rm = TRUE)
    m <- max(a, na.rm = TRUE)
    if (is.finite(m) && m > 0) a / m else a
  }, by = .(architecture_class, shot_number)]
  # Retain valid below-ground relative heights. GEDI V3 negative RH values are
  # not automatically errors after the quality and sensitivity filtering.
  grid_min <- floor(
    min(c(-cfg$waveform_grid_step_m, wave$height_m), na.rm = TRUE) /
      cfg$waveform_grid_step_m
  ) * cfg$waveform_grid_step_m
  grid_max <- ceiling(
    max(c(45, wave$height_m), na.rm = TRUE) /
      cfg$waveform_grid_step_m
  ) * cfg$waveform_grid_step_m
  grid <- seq(grid_min, grid_max, by = cfg$waveform_grid_step_m)
  wi <- wave[, {
    ord <- order(height_m)
    hh <- height_m[ord]
    aa <- amplitude_norm[ord]
    keep <- !duplicated(hh)
    if (sum(keep) < 2) NULL else data.table(
      height_m = grid,
      amplitude_norm = approx(hh[keep], aa[keep], xout = grid,
                              rule = 1, ties = mean)$y
    )
  }, by = .(architecture_class, shot_number)]
  wave_summary <- wi[is.finite(amplitude_norm), .(
    amplitude_mean = mean(amplitude_norm),
    amplitude_sd = sd(amplitude_norm),
    n_waveforms = uniqueN(shot_number)
  ), by = .(architecture_class, height_m)]
  fwrite(wave_summary, file.path(cfg$output_dir, "architecture_class_waveform_summary.csv"))
}

# Use one height range for every class so vertical architecture remains
# comparable across panels. Automatic limits encompass the class mean +/- SD
# RH ribbons and the available mean-waveform height support, rounded outward
# to 5-m increments. Zero is always included.
rh_plot_min <- min(c(0, rh_summary$mean_height_m - rh_summary$sd_height_m), na.rm = TRUE)
rh_plot_max <- max(c(0, rh_summary$mean_height_m + rh_summary$sd_height_m), na.rm = TRUE)
wave_plot_min <- if (is.null(wave_summary)) Inf else min(wave_summary$height_m, na.rm = TRUE)
wave_plot_max <- if (is.null(wave_summary)) -Inf else max(wave_summary$height_m, na.rm = TRUE)
profile_ymin <- if (is.finite(cfg$profile_y_limits_m[[1L]])) {
  cfg$profile_y_limits_m[[1L]]
} else {
  floor(min(0, rh_plot_min, wave_plot_min, na.rm = TRUE) / 5) * 5
}
profile_ymax <- if (is.finite(cfg$profile_y_limits_m[[2L]])) {
  cfg$profile_y_limits_m[[2L]]
} else {
  ceiling(max(45, rh_plot_max, wave_plot_max, na.rm = TRUE) / 5) * 5
}
if (!is.finite(profile_ymin) || !is.finite(profile_ymax) || profile_ymin >= profile_ymax) {
  stop("Invalid shared profile height limits: ", profile_ymin, ", ", profile_ymax)
}
message("Shared architecture-profile height axis: ", profile_ymin, " to ", profile_ymax, " m")

make_table <- function(cl) {
  z <- metric_summary[architecture_class == cl]
  if (!nrow(z)) return(ggplot() + theme_void())
  z[, row := rev(seq_len(.N))]
  ggplot(z, aes(y = row)) +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = row - 0.44, ymax = row + 0.44,
                  fill = factor(row %% 2)), colour = NA) +
    geom_text(aes(x = 0.01, label = metric), hjust = 0,
              size = cfg$text_size$table_text_mm) +
    geom_text(aes(x = 0.99, label = display), hjust = 1,
              size = cfg$text_size$table_text_mm) +
    scale_fill_manual(values = c("0" = "#D2D2D2", "1" = "#F4F4F4"), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, nrow(z) + 0.5), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))
}

make_map <- function(cl) {
  pts <- map_dt[architecture_class == cl]
  ggplot() +
    geom_sf(data = boundary, fill = "#F4F4F4", colour = "black",
            linewidth = cfg$geometry$map_boundary_line_mm) +
    geom_point(data = pts, aes(longitude, latitude), colour = class_cols[as.character(cl)],
               alpha = cfg$geometry$map_point_alpha,
               size = cfg$geometry$map_point_mm) +
    coord_sf(expand = FALSE) +
    # Explicit breaks avoid a compatibility problem in some older
    # ggplot2/sf versions where a break-generating function is passed to the
    # graticule code as if it were a numeric longitude vector.
    scale_x_continuous(
      breaks = c(-55, -53, -51),
      labels = c("55°W", "53°W", "51°W")
    ) +
    scale_y_continuous(
      breaks = c(-14, -10, -6, -2),
      labels = c("14°S", "10°S", "6°S", "2°S")
    ) +
    theme_bw(base_size = cfg$text_size$map_axis_text_pt) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0),
          axis.text = element_text(size = cfg$text_size$map_axis_text_pt))
}

make_profile <- function(cl) {
  rhs <- rh_summary[architecture_class == cl]
  col <- class_cols[as.character(cl)]
  if (!is.null(wave_summary) && any(wave_summary$architecture_class == cl)) {
    ws <- copy(wave_summary[architecture_class == cl])
    setorder(ws, height_m)
    # Place waveform in an arbitrary amplitude display range and overlay the
    # cumulative-energy RH curve on the same x extent.
    amp_scale <- 100
    smooth_n <- max(
      3L,
      as.integer(round(cfg$waveform_support_smooth_m / cfg$waveform_grid_step_m))
    )
    if (smooth_n %% 2L == 0L) smooth_n <- smooth_n + 1L
    ws[, amplitude_support := frollmean(
      amplitude_mean, n = smooth_n, align = "center", fill = NA_real_
    )]
    ws[!is.finite(amplitude_support), amplitude_support := amplitude_mean]
    class_peak <- max(ws$amplitude_support, na.rm = TRUE)
    display_threshold <- cfg$waveform_min_display_fraction * class_peak
    signal_rows <- which(
      is.finite(ws$amplitude_support) &
        ws$amplitude_support >= display_threshold
    )
    if (length(signal_rows)) {
      support_min <- ws$height_m[min(signal_rows)] - cfg$waveform_support_padding_m
      support_max <- ws$height_m[max(signal_rows)] + cfg$waveform_support_padding_m
      ws[, waveform_displayed := between(height_m, support_min, support_max)]
    } else {
      ws[, waveform_displayed := TRUE]
    }
    ws[, `:=`(
      x = fifelse(waveform_displayed, amplitude_mean * amp_scale, NA_real_),
      xmin = fifelse(
        waveform_displayed,
        pmax(0, (amplitude_mean - amplitude_sd) * amp_scale),
        NA_real_
      ),
      xmax = fifelse(
        waveform_displayed,
        (amplitude_mean + amplitude_sd) * amp_scale,
        NA_real_
      )
    )]
    rhs[, x_rh := rh]
    p <- ggplot() +
      geom_ribbon(data = ws, aes(y = height_m, xmin = xmin, xmax = xmax),
                  fill = col, alpha = cfg$geometry$waveform_ribbon_alpha) +
      geom_path(data = ws, aes(x = x, y = height_m), colour = col,
                linewidth = cfg$geometry$waveform_line_mm) +
      geom_ribbon(data = rhs, aes(x = x_rh, ymin = mean_height_m - sd_height_m,
                                  ymax = mean_height_m + sd_height_m),
                  fill = "grey55", alpha = cfg$geometry$rh_ribbon_alpha) +
      geom_line(data = rhs, aes(x = x_rh, y = mean_height_m),
                colour = "grey25", linewidth = cfg$geometry$rh_line_mm)
    xlabel <- "Relative waveform amplitude / RH percentile"
  } else {
    p <- ggplot(rhs, aes(rh, mean_height_m)) +
      geom_ribbon(aes(ymin = mean_height_m - sd_height_m,
                      ymax = mean_height_m + sd_height_m),
                  fill = col, alpha = cfg$geometry$waveform_ribbon_alpha) +
      geom_line(colour = col, linewidth = cfg$geometry$waveform_line_mm)
    xlabel <- "Relative-height percentile"
  }
  
  anchors <- rhs[rh %in% c(25L, 50L, 75L, 98L)]
  p +
    geom_hline(yintercept = 0, linetype = 2,
               linewidth = cfg$geometry$reference_line_mm) +
    geom_hline(yintercept = rhs[rh == 98, mean_height_m], linetype = 2,
               linewidth = cfg$geometry$reference_line_mm) +
    geom_segment(data = anchors[rh < 98],
                 aes(x = rh, xend = rh, y = 0, yend = mean_height_m),
                 arrow = arrow(length = grid::unit(cfg$geometry$rh_arrow_mm, "mm")),
                 linewidth = cfg$geometry$rh_segment_line_mm) +
    geom_text(data = anchors, aes(x = rh, y = mean_height_m, label = paste0("RH", rh)),
              nudge_y = 1.15, size = cfg$text_size$rh_label_mm) +
    coord_cartesian(
      xlim = c(0, 105), ylim = c(profile_ymin, profile_ymax), expand = FALSE
    ) +
    labs(x = sub(" / ", " /\n", xlabel, fixed = TRUE), y = "Height (m)") +
    theme_bw(base_size = cfg$text_size$profile_axis_text_pt) +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 0, 0.5, 0.5),
      axis.text = element_text(size = cfg$text_size$profile_axis_text_pt),
      axis.title = element_text(size = cfg$text_size$profile_axis_title_pt),
      axis.title.x = element_text(margin = margin(t = 1)),
      axis.title.y = element_text(margin = margin(r = 1))
    )
}

panels <- vector("list", n_class)
for (i in seq_along(class_ids)) {
  cl <- class_ids[i]
  top <- make_profile(cl) + make_map(cl) + plot_layout(widths = c(3.45, 0.95))
  panel <- (top / make_table(cl)) +
    plot_layout(heights = c(2.55, 1.45)) +
    plot_annotation(title = paste0("Class-", cl),
                    theme = theme(plot.title = element_text(
                      size = cfg$text_size$class_title_pt,
                      hjust = 0, margin = margin(b = 0.5)
                    )))
  # Preserve nested patchwork annotations when the class panels are assembled
  # into the final multi-panel figure.
  panel_wrapped <- wrap_elements(full = panel)
  panels[[i]] <- panel_wrapped
  ggsave(file.path(cfg$output_dir, sprintf("architecture_class_%02d.png", cl)), panel_wrapped,
         width = cfg$single_panel_width_mm, height = cfg$single_panel_height_mm,
         units = "mm", dpi = cfg$dpi, bg = "white")
}

full <- wrap_plots(panels, ncol = cfg$figure_columns) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(cfg$output_dir, "architecture_classes_complete.png"), full,
       width = cfg$final_width_mm, height = cfg$final_height_mm,
       units = "mm", dpi = cfg$dpi,
       limitsize = FALSE, bg = "white")
ggsave(file.path(cfg$output_dir, "architecture_classes_complete.pdf"), full,
       width = cfg$final_width_mm, height = cfg$final_height_mm, units = "mm",
       limitsize = FALSE, bg = "white")
saveRDS(full, file.path(cfg$output_dir, "architecture_classes_complete_plot.rds"))
message("Complete figure and class panels written to: ", cfg$output_dir)
