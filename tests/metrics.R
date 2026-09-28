# tests/metrics.R ------------------------------------------------------------
# Discrimination and risk ordering for the summed LLR, read off a finished
# branch_point run. FITS NOTHING — it loads `l_oof.rds` from the run directory,
# so it costs seconds and can be re-run freely.
#
#   Rscript tests/metrics.R                       # newest completed branch run
#   Rscript tests/metrics.R out/runs/branch_2026...  # a specific one
#
# Writes a new `metrics_<datetime>` run directory: tables, the two-panel risk
# curve as PNG, and a manifest naming the branch run it read.
#
# WHAT TO READ FIRST, in order:
#   1. auroc / auprc         does the score discriminate at all
#   2. cal_slope             is NAIVE SUMMING VALID. 1.000 means yes. Below 1
#                            means the L's double-count shared evidence and the
#                            sum is over-confident by roughly 1/slope.
#   3. n_sig_inversions      is the curve monotone where it matters
#   4. aggregation_gain      does summing beat the best single signal
#
# AGGREGATES ONLY (hard rule 1). The score is row-level and never printed.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
src_dir <- if (length(args)) args[1] else latest_run("branch")
if (is.null(src_dir) || !dir.exists(src_dir)) {
  stop("no completed branch run found. Run tests/branch_point.R first, or pass a run directory.",
       call. = FALSE)
}
l_path <- file.path(src_dir, "tables", "l_oof.rds")
if (!file.exists(l_path)) stop("no l_oof.rds in ", src_dir, call. = FALSE)

cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)

l_long   <- readRDS(l_path)
stay_ids <- sort(unique(l_long$stay_id))
stay_ids <- intersect(folds$stay_id[folds$split == "train"], stay_ids)

y <- tabs$cohort$mortality[match(as.character(stay_ids), as.character(tabs$cohort$stay_id))]
if (anyNA(y)) stop("metrics: a scored stay has no cohort row", call. = FALSE)
p_bar <- mean(y)

run <- new_run("metrics", cfg, note = sprintf("summed-LLR metrics read from %s",
                                              basename(src_dir)))
log_msg(run, sprintf("source run: %s | %d stays | event rate %.4f",
                     basename(src_dir), length(stay_ids), p_bar))

mats <- l_matrices(l_long, tabs, cfg, stay_ids, fill = "zero")
nb   <- cfg$metrics$n_bins %||% 20L
nboot <- cfg$metrics$n_boot %||% 200L

res <- list()
for (md in intersect(c("full", "cond", "meas", "intv"), names(mats))) {
  if (is.null(mats[[md]])) next
  res[[md]] <- metrics_report(run, mats[[md]], y, p_bar, label = md,
                              n_bins = nb, n_boot = nboot, seed = cfg$seed)
}

summ <- do.call(rbind, lapply(res, function(z) z$summary))
save_table(run, summ, "score_summary")

cat("\n=== summed LLR: discrimination ===\n\n")
print(summ[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
               "auprc", "auprc_lo", "auprc_hi", "auprc_lift")], row.names = FALSE)

cat("\n=== is naive summing valid? (slope 1 = yes) ===\n\n")
cal <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$cal)))
print(cal, row.names = FALSE)

cat("\n=== monotonicity of the 20-tile curve ===\n\n")
mono <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$mono)))
print(mono, row.names = FALSE)

cat("\n=== aggregation: does summing beat the best single signal? ===\n\n")
print(summ[, c("label", "auroc", "best_signal", "best_signal_auroc",
               "aggregation_gain")], row.names = FALSE)

cat("\n=== L_full, 20-tile observed mortality ===\n\n")
print(res$full$bins[, c("bin", "n", "deaths", "obs_rate", "lo", "hi", "mean_score")],
      row.names = FALSE)

finalize_run(run, extra = list(
  source_run = basename(src_dir),
  n_stays    = length(stay_ids),
  auroc_full = summ$auroc[summ$label == "full"],
  auprc_full = summ$auprc[summ$label == "full"],
  cal_slope_full = summ$cal_slope[summ$label == "full"]
))
cat(sprintf("\n  run directory: %s\n", run$path))
