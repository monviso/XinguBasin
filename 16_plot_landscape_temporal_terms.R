#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
prepared_path <- if (length(args) >= 1L) args[[1L]] else
  "E:/Xingu_rev/data/Xingu_GAM_prepared.rds"
output_dir <- if (length(args) >= 2L) args[[2L]] else
  "E:/Xingu_rev/data/history_figures/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Reading prepared data: ", prepared_path)
x <- readRDS(prepared_path)

term_info <- data.table(
  term = c(
    "d1", "d2", "focal_patch_area", "focal_patch_compactness",
    "nearest_nonforest_patch_area"
  ),
  panel = c("A", "B", "C", "D", "E"),
  label = c(
    "Distance from the forest edge",
    "Distance to the nearest other forest patch",
    "Focal forest-patch area",
    "Focal forest-patch compactness",
    "Nearest non-forest-patch area"
  ),
  y_label = c(
    "Distance (km)", "Distance (km)", "Area (km²)",
    "Compactness", "Area (km²)"
  )
)

missing_terms <- setdiff(term_info$term, names(x$landscape))
if (length(missing_terms)) {
  stop("Prepared object lacks landscape terms: ", paste(missing_terms, collapse = ", "))
}

train_idx <- which(x$base$analysis_split == "model")
if (!length(train_idx)) stop("No rows with analysis_split == 'model'")
class_levels <- sort(unique(as.integer(x$base$architecture_class)))

transform_for_plot <- function(term, z) {
  if (term %in% c("d1", "d2")) return(expm1(z) / 1000)  # m to km
  if (term %in% c("focal_patch_area", "nearest_nonforest_patch_area")) {
    return(expm1(z) / 100)                               # ha to km2
  }
  z
}

make_long <- function(term) {
  z <- x$landscape[[term]][train_idx, , drop = FALSE]
  stopifnot(ncol(z) == 20L)
  z <- transform_for_plot(term, z)
  data.table(
    sample_id = rep(x$base$sample_id[train_idx], each = 20L),
    architecture_class = factor(
      rep(x$base$architecture_class[train_idx], each = 20L),
      levels = class_levels,
      labels = paste("Class", class_levels)
    ),
    temporal_lag_years = rep(0:19, times = length(train_idx)),
    value = as.vector(t(z)),
    sampling_weight = rep(
      x$base$sampling_weight_normalized[train_idx], each = 20L
    ),
    term = term
  )
}

plot_dt <- rbindlist(lapply(term_info$term, make_long), use.names = TRUE)
plot_dt <- plot_dt[
  is.finite(value) & value >= 0 &
    is.finite(sampling_weight) & sampling_weight > 0
]
plot_dt <- merge(plot_dt, term_info, by = "term", sort = FALSE)
plot_dt[, facet_label := factor(
  paste0(panel, "  ", label),
  levels = paste0(term_info$panel, "  ", term_info$label)
)]

weighted_quantile <- function(value, weight, probability) {
  keep <- is.finite(value) & is.finite(weight) & weight > 0
  value <- value[keep]
  weight <- weight[keep]
  if (!length(value)) return(NA_real_)
  ord <- order(value)
  value <- value[ord]
  weight <- weight[ord]
  value[which(cumsum(weight) / sum(weight) >= probability)[1L]]
}

summary_dt <- plot_dt[, .(
  n = .N,
  n_samples = uniqueN(sample_id),
  weighted_mean = weighted.mean(value, sampling_weight),
  weighted_median = weighted_quantile(value, sampling_weight, 0.5),
  weighted_q25 = weighted_quantile(value, sampling_weight, 0.25),
  weighted_q75 = weighted_quantile(value, sampling_weight, 0.75)
), by = .(
  term, panel, label, y_label, facet_label,
  architecture_class, temporal_lag_years
)]

class_colours <- c(
  "Class 1" = "#EFC000",
  "Class 2" = "#E76F00",
  "Class 3" = "#8E2C62",
  "Class 4" = "#44236D",
  "Class 5" = "#3FA3A5",
  "Class 6" = "#8FD694",
  "Class 7" = "#0072B2"
)

base_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(1.25, "cm"),
    axis.title = element_text(size=18), axis.text = element_text(size=15)
  )

smooth_layer <- geom_smooth(
  method = "gam",
  formula = y ~ s(x, bs = "cs", k = 6),
  se = TRUE,
  linewidth = 0.95,
  alpha = 0.12,
  na.rm = TRUE
)

p_all <- ggplot(
  plot_dt,
  aes(
    temporal_lag_years, value,
    colour = architecture_class,
    fill = architecture_class,
    weight = sampling_weight
  )
) +
  smooth_layer +
  geom_point(
    data = summary_dt[temporal_lag_years == 0L],
    aes(temporal_lag_years, weighted_mean, colour = architecture_class),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.1,
    size = 2.5
  ) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = class_colours, drop = FALSE) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 19), limits = c(0, 19),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Temporal lag (years)",
    y = "Value (units differ among panels)",
    colour = NULL,
    fill = NULL
  ) +
  base_theme +
  guides(
    colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1.35)),
    fill = "none"
  )

p_all
ggsave(
  file.path(output_dir, "Fig_landscape_temporal_terms.png"),
  p_all, width = 10.5, height = 11, dpi = 400, bg = "white"
)
ggsave(
  file.path(output_dir, "Fig_landscape_temporal_terms.pdf"),
  p_all, width = 10.5, height = 11
)

# Also save each term as an independent panel with its correct y-axis label.
for (i in seq_len(nrow(term_info))) {
  info <- term_info[i]
  dd <- plot_dt[term == info$term]
  ss <- summary_dt[term == info$term & temporal_lag_years == 0L]
  p <- ggplot(
    dd,
    aes(
      temporal_lag_years, value,
      colour = architecture_class,
      fill = architecture_class,
      weight = sampling_weight
    )
  ) +
    geom_smooth(
      method = "gam", formula = y ~ s(x, bs = "cs", k = 6),
      se = TRUE, linewidth = 1.05, alpha = 0.13, na.rm = TRUE
    ) +
    geom_point(
      data = ss,
      aes(temporal_lag_years, weighted_mean, colour = architecture_class),
      inherit.aes = FALSE, shape = 4, stroke = 1.2, size = 3
    ) +
    scale_colour_manual(values = class_colours, drop = FALSE) +
    scale_fill_manual(values = class_colours, drop = FALSE) +
    scale_x_continuous(
      breaks = c(0, 5, 10, 15, 19), limits = c(0, 19),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      tag = info$panel,
      title = info$label,
      x = "Temporal lag (years)", y = info$y_label,
      colour = NULL, fill = NULL
    ) +
    base_theme +
    theme(
      strip.text = element_blank(),
      plot.tag = element_text(face = "bold", size = 18),
      plot.title = element_text(face = "plain", size = 15),
      aspect.ratio = 0.72
    ) +
    guides(
      colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.4)),
      fill = "none"
    )
  stem <- sprintf("Fig_landscape_%s", info$term)
  ggsave(file.path(output_dir, paste0(stem, ".png")), p,
         width = 7.2, height = 5.2, dpi = 400, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), p,
         width = 7.2, height = 5.2)
}

fwrite(
  summary_dt,
  file.path(output_dir, "Fig_landscape_temporal_terms_summary.csv")
)
fwrite(
  plot_dt,
  file.path(output_dir, "Fig_landscape_temporal_terms_plot_data.csv.gz"),
  compress = "gzip"
)

message("Training samples represented: ", uniqueN(plot_dt$sample_id))
message("Long rows after structural-NA removal: ", format(nrow(plot_dt), big.mark = ","))
message("Outputs written to: ", output_dir)
