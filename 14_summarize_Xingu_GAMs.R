#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
fit_dir <- if (length(args) >= 1L) args[[1L]] else "/scratch/mr3882/temp/Xingu_GAM/fits_reduced_basis"
output_dir <- if (length(args) >= 2L) args[[2L]] else "/scratch/mr3882/temp/Xingu_GAM/summary"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(fit_dir, pattern = "_metrics\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No model metric files found under ", fit_dir)
m <- rbindlist(lapply(files, fread), fill = TRUE)
architecture_classes <- sort(unique(as.integer(m$focal_class)))
expected <- length(architecture_classes) * 10L * 2L
if (nrow(m) != expected) {
  warning("Expected ", expected, " metric rows (",
          length(architecture_classes) * 10L,
          " fits x 2 splits), found ", nrow(m))
}
setorder(m, focal_class, weighting, model_id, split)
m[, delta_AIC := AIC - min(AIC), by = .(focal_class, weighting)]
fwrite(m, file.path(output_dir, "all_GAM_metrics.csv"))
fwrite(m[split == "validation"], file.path(output_dir, "validation_metrics.csv"))

# Three-group commonality analysis on the seven lagged models.
lag <- unique(m[representation == "lagged" & weighting == "weighted", .(
  focal_class, weighting, model_id, explained_deviance
)])
wide <- dcast(lag, focal_class + weighting ~ model_id,
              value.var = "explained_deviance")
needed <- c("L", "S", "E", "LS", "LE", "SE", "LSE_lagged")
if (!all(needed %in% names(wide))) {
  warning("Incomplete lagged model set; commonality output not written")
} else {
  v <- wide[, .(
    focal_class, weighting,
    unique_L = LSE_lagged - SE,
    unique_S = LSE_lagged - LE,
    unique_E = LSE_lagged - LS,
    shared_LSE = L + S + E - LS - LE - SE + LSE_lagged
  )]
  v[, shared_LS := wide$L + wide$S - wide$LS - shared_LSE]
  v[, shared_LE := wide$L + wide$E - wide$LE - shared_LSE]
  v[, shared_SE := wide$S + wide$E - wide$SE - shared_LSE]
  v[, full_explained_deviance := wide$LSE_lagged]
  fwrite(v, file.path(output_dir, "three_group_commonality_deviance.csv"))
}

comparison <- m[
  split == "validation" & model_id %chin% c("LSE_lagged", "LSE_recent", "LSE_average"),
  .(focal_class, weighting, model_id, AIC, delta_AIC, explained_deviance,
    balanced_accuracy, f1, sensitivity, specificity, precision, n)
]
fwrite(comparison, file.path(output_dir, "full_model_lag_comparison.csv"))

status <- rbindlist(list(
  CJ(focal_class = architecture_classes, weighting = "weighted",
     model_id = c("L", "S", "E", "LS", "LE", "SE",
                  "LSE_lagged", "LSE_recent", "LSE_average")),
  CJ(focal_class = architecture_classes, weighting = "unweighted",
     model_id = "LSE_lagged")
))
done <- unique(m[, .(focal_class, weighting, model_id)])
status[, completed := paste(focal_class, weighting, model_id) %chin%
         paste(done$focal_class, done$weighting, done$model_id)]
fwrite(status, file.path(output_dir, "model_completion_status.csv"))
message("Summary written to ", output_dir)
