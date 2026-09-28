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

stopifnot(
  is.data.frame(x$base),
  is.matrix(x$landscape$d1),
  nrow(x$landscape$d1) == nrow(x$base),
  ncol(x$landscape$d1) == 20L
)

# Use the spatially independent model split only. D1 was stored as log1p(m).
train_idx <- which(x$base$analysis_split == "model")
if (!length(train_idx)) stop("No rows with analysis_split == 'model'")
class_levels <- sort(unique(as.integer(x$base$architecture_class)))

d1_km <- expm1(x$landscape$d1[train_idx, , drop = FALSE]) / 1000
plot_dt <- data.table(
  sample_id = rep(x$base$sample_id[train_idx], each = 20L),
  architecture_class = factor(
    rep(x$base$architecture_class[train_idx], each = 20L),
    levels = class_levels,
    labels = paste("Class", class_levels)
  ),
  temporal_lag_years = rep(0:19, times = length(train_idx)),
  edge_distance_km = as.vector(t(d1_km)),
  sampling_weight = rep(
    x$base$sampling_weight_normalized[train_idx], each = 20L
  )
)

# Historical predictors are structurally NA before a location became forest.
plot_dt <- plot_dt[
  is.finite(edge_distance_km) & edge_distance_km >= 0 &
    is.finite(sampling_weight) & sampling_weight > 0
]
if (!nrow(plot_dt)) stop("No finite training D1 observations")

class_colours <- c(
  "Class 1" = "#EFC000",
  "Class 2" = "#E76F00",
  "Class 3" = "#8E2C62",
  "Class 4" = "#44236D",
  "Class 5" = "#3FA3A5",
  "Class 6" = "#8FD694",
  "Class 7" = "#0072B2"
)

# Weighted lag summaries are exported independently from the smoother.
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
  weighted_mean_edge_distance_km = weighted.mean(
    edge_distance_km, sampling_weight
  ),
  weighted_median_edge_distance_km = weighted_quantile(
    edge_distance_km, sampling_weight, 0.5
  ),
  weighted_q25_edge_distance_km = weighted_quantile(
    edge_distance_km, sampling_weight, 0.25
  ),
  weighted_q75_edge_distance_km = weighted_quantile(
    edge_distance_km, sampling_weight, 0.75
  )
), by = .(architecture_class, temporal_lag_years)]

p <- ggplot(
  plot_dt,
  aes(
    x = temporal_lag_years,
    y = edge_distance_km,
    colour = architecture_class,
    fill = architecture_class,
    weight = sampling_weight
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 6),
    se = TRUE,
    linewidth = 1.05,
    alpha = 0.13,
    na.rm = TRUE
  ) +
  geom_point(
    data = summary_dt[temporal_lag_years == 0L],
    aes(
      x = temporal_lag_years,
      y = weighted_mean_edge_distance_km,
      colour = architecture_class
    ),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.2,
    size = 3
  ) +
  scale_colour_manual(values = class_colours, drop = FALSE) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 19),
    limits = c(0, 19),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    tag = "A",
    title = "Distance from the forest edge",
    x = "Temporal lag (years)",
    y = "Distance (km)",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "plain", size = 13),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(1.4, "cm"),
    aspect.ratio = 0.72
  ) +
  guides(
    colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.4)),
    fill = "none"
  )
p

png_path <- file.path(output_dir, "Fig_edge_distance_panel_A.png")
pdf_path <- file.path(output_dir, "Fig_edge_distance_panel_A.pdf")
csv_path <- file.path(output_dir, "Fig_edge_distance_panel_A_summary.csv")
long_path <- file.path(output_dir, "Fig_edge_distance_panel_A_plot_data.csv.gz")

ggsave(png_path, p, width = 7.2, height = 5.2, dpi = 400, bg = "white")
ggsave(pdf_path, p, width = 7.2, height = 5.2, device = cairo_pdf)
fwrite(summary_dt, csv_path)
fwrite(plot_dt, long_path, compress = "gzip")

message("Training samples represented: ", uniqueN(plot_dt$sample_id))
message("Plot rows after structural-NA removal: ", format(nrow(plot_dt), big.mark = ","))
message("Written: ", png_path)
message("Written: ", pdf_path)
message("Written: ", csv_path)
message("Written: ", long_path)
