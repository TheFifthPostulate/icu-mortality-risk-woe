# paper/make_signal_auroc_meas.R -------------------------------------------------
# Per-channel out-of-fold AUROC of the measurement-only weight of evidence, for the
# F6 figure (joint vs measurement-only). The run of record exports this table for
# the joint and the conditional matrices only (signal_auroc_full/cond.csv), so we
# compute the measurement-only one here from the stored targets, with the same
# function and inputs as the `signal_auroc_full` target (_targets.R:350,
# R/09_metrics.R signal_auroc()). No refit; no change to _targets/ or out/runs/.
#
# Check first: the joint table is recomputed from the store and must equal the
# exported signal_auroc_full.csv of the run of record, otherwise the script stops.
# Output is aggregate only (one row per channel).
#
#   Rscript paper/make_signal_auroc_meas.R
#
# Writes paper/figs/coefs/signal_auroc_meas_oof.csv.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(targets))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

REF <- "out/runs/internal_20260909T112643/tables/signal_auroc_full.csv"
OUT <- "paper/figs/coefs/signal_auroc_meas_oof.csv"

lm <- tar_read(l_mats_zero)
y  <- tar_read(y_train)
stopifnot(all(c("full", "meas") %in% names(lm)), nrow(lm$full) == length(y))

# --- check: the store reproduces the exported joint table ----------------------
full_now <- signal_auroc(lm$full, y)
full_ref <- read.csv(REF, stringsAsFactors = FALSE)
m <- match(full_ref$signal, full_now$signal)
if (anyNA(m)) stop("signal sets differ between the store and ", REF)
d <- max(abs(full_now$auroc[m] - full_ref$auroc), abs(full_now$auprc[m] - full_ref$auprc))
cat(sprintf("joint check: %d channels, max |diff| = %.2e\n", nrow(full_ref), d))
if (d > 1e-9) stop("the store does not reproduce ", REF, "; not writing the measurement-only table")

# --- measurement-only ----------------------------------------------------------
M <- lm$meas
miss <- setdiff(colnames(lm$full), colnames(M))
if (length(miss)) {
  # unpaired channels: measurement-only and joint are the same model (l_long alias)
  cat("aliasing to joint (unpaired, absent from meas):", paste(miss, collapse = ", "), "\n")
  M <- cbind(M, lm$full[, miss, drop = FALSE])
}
M <- M[, colnames(lm$full), drop = FALSE]
meas <- signal_auroc(M, y)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(meas, OUT, row.names = FALSE)
cat("wrote", OUT, "with", nrow(meas), "channels\n")
