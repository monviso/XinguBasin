#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =============================================================================
# ARGUMENTS AND PATHS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

prepared_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "E:/Xingu_rev/data/Xingu_GAM_prepared.rds"
}

output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "E:/Xingu_rev/data/history_figures/"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading prepared data: ", prepared_path)
x <- readRDS(prepared_path)

# =============================================================================
# CONFIGURATION
# =============================================================================

radii <- c(100L, 500L, 1000L, 3000L)
lags <- 0:19

class_levels <- sort(
  unique(as.integer(x$base$architecture_class))
)

# MapBiomas proportions are stored as uint8-scaled fractions:
#   0   = 0%
#   254 = 100%
#   255 = NoData
proportion_scale <- 254
proportion_nodata <- 255

# Row order: anthropogenic covers followed by natural vegetation and forest.
term_info <- data.table(
  term = c(
    "urban",
    "mining",
    "pasture",
    "agriculture",
    "silviculture_oil_palm",
    "grassland",
    "wetland",
    "water",
    "forest"
  ),
  row_label = c(
    "Urban",
    "Mining",
    "Pasture",
    "Agriculture",
    "Silviculture + oil palm",
    "Grassland",
    "Wetland",
    "Water",
    "Forest"
  )
)

# =============================================================================
# INPUT VALIDATION
# =============================================================================

required_base_columns <- c(
  "architecture_class",
  "analysis_split",
  "sampling_weight_normalized"
)

missing_base_columns <- setdiff(
  required_base_columns,
  names(x$base)
)

if (length(missing_base_columns)) {
  stop(
    "Prepared object x$base lacks columns: ",
    paste(missing_base_columns, collapse = ", ")
  )
}

if (is.null(x$lulc)) {
  stop("Prepared object does not contain x$lulc")
}

missing_terms <- setdiff(term_info$term, names(x$lulc))

if (length(missing_terms)) {
  stop(
    "Prepared object lacks LULC groups: ",
    paste(missing_terms, collapse = ", ")
  )
}

train_idx <- which(x$base$analysis_split == "model")

if (!length(train_idx)) {
  stop("No rows with analysis_split == 'model'")
}

weights <- as.numeric(
  x$base$sampling_weight_normalized[train_idx]
)

classes <- as.integer(
  x$base$architecture_class[train_idx]
)

if (length(weights) != length(train_idx)) {
  stop("Sampling-weight vector has an unexpected length")
}

if (length(classes) != length(train_idx)) {
  stop("Architecture-class vector has an unexpected length")
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

stored_proportion_to_percent <- function(value) {
  value <- as.numeric(value)
  
  # Convert explicit NoData values to NA.
  nodata <- is.finite(value) & value == proportion_nodata
  value[nodata] <- NA_real_
  
  # Stop if unexpected finite values occur.
  invalid <- is.finite(value) & (
    value < 0 |
      value > proportion_scale
  )
  
  if (any(invalid)) {
    invalid_range <- range(value[invalid], na.rm = TRUE)
    
    stop(
      "LULC proportion values outside the expected 0–255 range: ",
      paste(invalid_range, collapse = " to ")
    )
  }
  
  # Convert stored values to percentage cover.
  100 * value / proportion_scale
}

weighted_cell_mean <- function(value, weight) {
  keep <- (
    is.finite(value) &
      is.finite(weight) &
      weight > 0
  )
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(
    x = value[keep],
    w = weight[keep]
  )
}

safe_min <- function(value) {
  value <- value[is.finite(value)]
  
  if (!length(value)) {
    return(NA_real_)
  }
  
  min(value)
}

safe_max <- function(value) {
  value <- value[is.finite(value)]
  
  if (!length(value)) {
    return(NA_real_)
  }
  
  max(value)
}

# =============================================================================
# SUMMARISE ONE LULC TERM
# =============================================================================

# Matrices are assumed to be radius-major:
#
#   r100:  lag 0, ..., lag 19
#   r500:  lag 0, ..., lag 19
#   r1000: lag 0, ..., lag 19
#   r3000: lag 0, ..., lag 19

summarise_term <- function(term_value) {
  z <- x$lulc[[term_value]][
    train_idx,
    ,
    drop = FALSE
  ]
  
  expected_columns <- length(radii) * length(lags)
  
  if (ncol(z) != expected_columns) {
    stop(
      term_value,
      " has ",
      ncol(z),
      " columns; expected ",
      expected_columns
    )
  }
  
  output <- vector(
    "list",
    length(class_levels) * ncol(z)
  )
  
  q <- 0L
  
  for (class_value in class_levels) {
    class_idx <- which(classes == class_value)
    
    if (!length(class_idx)) {
      warning(
        "No model observations for architecture class ",
        class_value
      )
      
      next
    }
    
    class_weights <- weights[class_idx]
    
    for (column_index in seq_len(ncol(z))) {
      q <- q + 1L
      
      radius_index <- (
        (column_index - 1L) %/% length(lags)
      ) + 1L
      
      lag_index <- (
        (column_index - 1L) %% length(lags)
      ) + 1L
      
      stored_value <- z[class_idx, column_index]
      
      percentage_value <- stored_proportion_to_percent(
        stored_value
      )
      
      keep <- (
        is.finite(percentage_value) &
          is.finite(class_weights) &
          class_weights > 0
      )
      
      output[[q]] <- data.table(
        term = term_value,
        architecture_class = class_value,
        radius_m = radii[[radius_index]],
        temporal_lag_years = lags[[lag_index]],
        weighted_mean_percent = weighted_cell_mean(
          percentage_value,
          class_weights
        ),
        n = sum(keep),
        weighted_n = if (any(keep)) {
          sum(class_weights[keep])
        } else {
          0
        }
      )
    }
  }
  
  output <- output[!vapply(output, is.null, logical(1))]
  
  if (!length(output)) {
    return(data.table())
  }
  
  rbindlist(output)
}

# =============================================================================
# BUILD THE COMPLETE SURFACE TABLE
# =============================================================================

message("Calculating class-specific weighted LULC histories")

surface_dt <- rbindlist(
  lapply(term_info$term, summarise_term),
  use.names = TRUE,
  fill = TRUE
)

surface_dt <- merge(
  surface_dt,
  term_info,
  by = "term",
  sort = FALSE
)

# Restore the intended row order after merging.
surface_dt[
  ,
  term := factor(
    term,
    levels = term_info$term
  )
]

surface_dt[
  ,
  row_label := factor(
    row_label,
    levels = term_info$row_label
  )
]

surface_dt[
  ,
  architecture_class := factor(
    architecture_class,
    levels = class_levels,
    labels = paste0("Class-", class_levels)
  )
]

surface_dt[
  ,
  radius_factor := factor(
    radius_m,
    levels = radii,
    labels = format(
      radii,
      big.mark = ",",
      scientific = FALSE
    )
  )
]

setorder(
  surface_dt,
  term,
  architecture_class,
  radius_m,
  temporal_lag_years
)

# =============================================================================
# SAVE THE NUMERICAL HISTORIES
# =============================================================================

surface_csv <- file.path(
  output_dir,
  "Fig_LULC_space_time_weighted_means_percent.csv"
)

fwrite(surface_dt, surface_csv)

message("Written surface table: ", surface_csv)

# =============================================================================
# PLOTTING FUNCTION
# =============================================================================

make_row_plot <- function(
    term_value,
    show_x = FALSE,
    show_headers = FALSE
) {
  info <- term_info[term == term_value]
  
  dd <- surface_dt[
    as.character(term) == term_value
  ]
  
  if (!nrow(dd)) {
    stop("No summarized data available for term: ", term_value)
  }
  
  row_max <- safe_max(dd$weighted_mean_percent)
  
  if (!is.finite(row_max) || row_max <= 0) {
    row_max <- 1
  }
  
  ggplot(
    dd,
    aes(
      x = radius_factor,
      y = temporal_lag_years,
      fill = weighted_mean_percent
    )
  ) +
    geom_tile() +
    facet_grid(
      . ~ architecture_class,
      switch = "x"
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = 1,
      limits = c(0, row_max),
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::label_number(
        accuracy = 0.1,
        suffix = "%"
      ),
      na.value = "grey90",
      name = "Mean cover"
    ) +
    scale_y_continuous(
      breaks = c(0, 5, 10, 15, 19),
      limits = c(-0.5, 19.5),
      expand = c(0, 0)
    ) +
    labs(
      x = if (show_x) {
        "Spatial radius (m)"
      } else {
        NULL
      },
      y = paste0(
        as.character(info$row_label),
        "\nTemporal lag (years)"
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      
      panel.spacing.x = grid::unit(
        1.2,
        "mm"
      ),
      
      strip.text.x = if (show_headers) {
        element_text(
          face = "bold",
          size = 10
        )
      } else {
        element_blank()
      },
      
      strip.background = if (show_headers) {
        element_rect(
          fill = "grey88",
          colour = "grey45",
          linewidth = 0.3
        )
      } else {
        element_blank()
      },
      
      axis.title.y = element_text(
        size = 8.5,
        lineheight = 0.9
      ),
      
      axis.text.y = element_text(
        size = 7,
        colour = "black"
      ),
      
      axis.ticks.y = element_line(
        linewidth = 0.25
      ),
      
      axis.text.x = if (show_x) {
        element_text(
          size = 7,
          angle = 0,
          hjust = 0.5,
          colour = "black"
        )
      } else {
        element_blank()
      },
      
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.25)
      } else {
        element_blank()
      },
      
      axis.title.x = if (show_x) {
        element_text(size = 9)
      } else {
        element_blank()
      },
      
      legend.position = "right",
      
      legend.title = element_text(
        size = 7
      ),
      
      legend.text = element_text(
        size = 6
      ),
      
      legend.key.height = grid::unit(
        8,
        "mm"
      ),
      
      plot.margin = margin(
        t = 1,
        r = 2,
        b = 1,
        l = 2,
        unit = "mm"
      )
    )
}

# =============================================================================
# CREATE THE NINE ROW PLOTS
# =============================================================================

plots <- lapply(
  seq_len(nrow(term_info)),
  function(i) {
    make_row_plot(
      term_value = term_info$term[[i]],
      show_x = i == nrow(term_info),
      show_headers = i == 1L
    )
  }
)

# =============================================================================
# COMBINE ALL ROWS
# =============================================================================

draw_combined <- function(
    filename,
    device = c("png", "pdf")
) {
  device <- match.arg(device)
  
  if (device == "png") {
    png(
      filename = filename,
      width = 3900,
      height = 5700,
      res = 350,
      bg = "white"
    )
  } else {
    pdf(
      file = filename,
      width = 11.2,
      height = 16.2,
      onefile = TRUE
    )
  }
  
  on.exit(dev.off(), add = TRUE)
  
  grid::grid.newpage()
  
  layout <- grid::grid.layout(
    nrow = nrow(term_info),
    ncol = 1L
  )
  
  top_viewport <- grid::viewport(
    layout = layout
  )
  
  grid::pushViewport(top_viewport)
  
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      vp = grid::viewport(
        layout.pos.row = i,
        layout.pos.col = 1L
      )
    )
  }
  
  grid::popViewport()
}

combined_png <- file.path(
  output_dir,
  "Fig_LULC_space_time_heatmaps_percent.png"
)

combined_pdf <- file.path(
  output_dir,
  "Fig_LULC_space_time_heatmaps_percent.pdf"
)

draw_combined(
  filename = combined_png,
  device = "png"
)

draw_combined(
  filename = combined_pdf,
  device = "pdf"
)

# =============================================================================
# SAVE EACH LULC ROW INDEPENDENTLY
# =============================================================================

for (i in seq_len(nrow(term_info))) {
  term_value <- term_info$term[[i]]
  
  stem <- paste0(
    "Fig_LULC_space_time_",
    term_value,
    "_percent"
  )
  
  independent_plot <- make_row_plot(
    term_value = term_value,
    show_x = TRUE,
    show_headers = TRUE
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(stem, ".png")
    ),
    plot = independent_plot,
    width = 11.2,
    height = 2.1,
    dpi = 350,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(stem, ".pdf")
    ),
    plot = independent_plot,
    width = 11.2,
    height = 2.1
  )
}

# =============================================================================
# QUALITY-ASSURANCE OUTPUT
# =============================================================================

qa <- surface_dt[
  ,
  .(
    cells = .N,
    
    cells_with_data = sum(
      is.finite(weighted_mean_percent)
    ),
    
    minimum_cell_n = min(
      n,
      na.rm = TRUE
    ),
    
    maximum_cell_n = max(
      n,
      na.rm = TRUE
    ),
    
    minimum_weighted_mean_percent = safe_min(
      weighted_mean_percent
    ),
    
    maximum_weighted_mean_percent = safe_max(
      weighted_mean_percent
    )
  ),
  by = .(
    term,
    row_label,
    architecture_class
  )
]

qa_file <- file.path(
  output_dir,
  "Fig_LULC_space_time_QA_percent.csv"
)

fwrite(qa, qa_file)

# Additional conversion check.
conversion_qa <- surface_dt[
  ,
  .(
    minimum_percent = safe_min(
      weighted_mean_percent
    ),
    maximum_percent = safe_max(
      weighted_mean_percent
    ),
    values_below_zero = sum(
      weighted_mean_percent < 0,
      na.rm = TRUE
    ),
    values_above_100 = sum(
      weighted_mean_percent > 100,
      na.rm = TRUE
    )
  ),
  by = .(
    term,
    row_label
  )
]

conversion_qa_file <- file.path(
  output_dir,
  "Fig_LULC_proportion_conversion_QA.csv"
)

fwrite(conversion_qa, conversion_qa_file)

if (
  any(conversion_qa$values_below_zero > 0) ||
  any(conversion_qa$values_above_100 > 0)
) {
  warning(
    "Converted LULC percentages outside the expected 0–100% range. ",
    "Inspect: ",
    conversion_qa_file
  )
}

message("Completed.")
message("Combined PNG: ", combined_png)
message("Combined PDF: ", combined_pdf)
message("QA table: ", qa_file)
message("Conversion QA: ", conversion_qa_file)