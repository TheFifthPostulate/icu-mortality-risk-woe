# tests/branch_point.R -------------------------------------------------------
# Fit layer 1 out-of-fold, assemble the L matrices, and read the eigenspectrum.
#
# THIS IS THE BRANCH POINT. If one component dominates, layer 2 is weakly
# motivated and the paper's claims change. Nothing downstream of this should be
# built before its output has been read.
#
#   Rscript tests/branch_point.R          # all 5 folds  (~9 min, 215 fits)
#   Rscript tests/branch_point.R 1        # fold 1 only  (~1.5 min, 43 fits)
#
# The one-fold form is for checking the plumbing, not for deciding anything: it
# scores ~8,240 stays instead of ~41,200, so the correlations are noisier. The
# spectrum it prints is indicative.
#
# EVERYTHING IS SAVED TO A RUN DIRECTORY BEFORE THE ANALYSIS RUNS. 215 fits cost
# nine minutes and the analysis after them costs seconds; if anything in the
# analysis throws, the fits must survive it. `l_oof` is written the moment
# run_layer1() returns, so a failed reading can be re-read from disk without
# re-fitting anything.
#
# `meas` is fitted here because it is nearly free once the frames are built and
# it is what `full` aliases to on the 7 unpaired signals. Only `full` and `intv`
# are load-bearing; L_cond = L_full - L_intv is derived, not fitted.
#
# AGGREGATES ONLY on print (hard rule 1): eigenvalues, correlations, counts. The
# L matrix is row-level — it goes to the run directory, never to the console.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
n_folds_run <- if (length(args) && nzchar(args[1])) as.integer(args[1]) else NULL

# An OPTIONAL config path, added 2026-09-03. The ablation this script doubles as
# -- both conditional constructs switched off, which is what
# `tests/refactor_identity.R` compares against -- previously required editing
# `config/config.yml` in place. That has a cost far beyond the edit: the config
# file is a `format = "file"` target, so touching it invalidates `cfg` and every
# one of the 258 fits downstream, and the cache then has to be rebuilt twice to
# get back where it started. Passing a path instead leaves the project config
# and the `_targets/` cache untouched.
#
#   Rscript tests/branch_point.R 1 /path/to/ablation.yml
cfg_path <- if (length(args) >= 2L && nzchar(args[2])) args[2] else "config/config.yml"
cfg    <- load_config(cfg_path)
if (!identical(cfg_path, "config/config.yml")) {
  message("branch_point: config OVERRIDE in use -- ", cfg_path)
}
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)

jobs <- layer1_jobs(cfg)
jobs <- jobs[jobs$role == "oof", , drop = FALSE]          # final fits are for test/eICU
if (!is.null(n_folds_run)) jobs <- jobs[jobs$fold <= n_folds_run, , drop = FALSE]

used_folds <- sort(unique(jobs$fold))
stay_ids <- folds$stay_id[folds$split == "train" & folds$fold %in% used_folds]

run <- new_run("branch", cfg, note = sprintf(
  "out-of-fold layer 1 over fold(s) %s; %d fits; %d stays scored",
  paste(used_folds, collapse = ","), sum(jobs$fit), length(stay_ids)))

log_msg(run, sprintf("layer 1: %d fits over fold(s) %s, scoring %d stays",
                     sum(jobs$fit), paste(used_folds, collapse = ","), length(stay_ids)))

t0 <- Sys.time()
r <- run_layer1(tabs, cfg, folds, priors, jobs = jobs, keep_final = FALSE,
                verbose = FALSE)
log_msg(run, sprintf("%d fits in %.1f min", sum(jobs$fit),
                     as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# INSURANCE, written before anything can throw. `l_oof` is row-level and stays
# in the run directory; `save_object` does not print it.
save_object(run, r$l, "l_oof")
save_table(run, r$diagnostics, "layer1_diagnostics", subdir = "diagnostics")
save_object(run, jobs, "layer1_jobs")

# --- diagnostics ------------------------------------------------------------
tr <- diagnostics_summary(r$diagnostics, cfg)
if (nrow(tr)) save_table(run, tr, "layer1_triage", subdir = "diagnostics")
cat("\n  worst 5 by concurvity (observed, not the worst-case bound):\n")
d <- r$diagnostics[order(-r$diagnostics$concurvity_obs), ]
print(utils::head(d[, c("signal", "model", "n_rows", "dev_expl", "edf_total",
                        "k_index_min", "concurvity_obs")], 5), row.names = FALSE)

# --- the L matrices ---------------------------------------------------------
mats <- l_matrices(r$l, tabs, cfg, stay_ids, fill = "zero")
cat(sprintf("\n  L matrices: %d stays x %d signals; models: %s\n",
            nrow(mats$full), ncol(mats$full), paste(names(mats), collapse = ", ")))

ls_full <- l_summary(mats$full)
cat("\n  L_full, per signal (aggregates only):\n")
print(ls_full, row.names = FALSE)
save_table(run, ls_full, "l_full_summary")

# The identity that has to hold, checked rather than trusted: on the 7 unpaired
# signals L_intv is 0 by construction, so L_cond must equal L_full exactly.
unp <- cfg$signals[vapply(cfg$signals, function(sg)
  !length(interventions_of(sg, cfg)), logical(1))]
stopifnot(all(abs(mats$cond[, unp] - mats$full[, unp]) < 1e-12))
cat(sprintf("  check: L_cond == L_full on all %d unpaired signals   [ok]\n", length(unp)))

# --- the reading ------------------------------------------------------------
bp <- branch_point(r$l, tabs, cfg, stay_ids)
save_table(run, bp, "branch_point")
report_branch_point(bp)

sp <- shared_covariate_pairs(r$l, tabs, cfg, stay_ids)
save_table(run, sp, "shared_covariate_pairs")
cat("  shared-covariate pairs — how much of the correlation is the\n")
cat("  intervention block rather than physiology:\n\n")
print(sp, row.names = FALSE)

cat("\n  full spectrum, primary reading (L_full, zero-filled, 19 signals):\n\n")
es <- eigenspectrum(l_correlation(mats$full))
print(es, row.names = FALSE)
save_table(run, es, "eigenspectrum_full_zero")

# Sigma itself, for layer 2 to read later without re-fitting.
save_object(run, l_correlation(mats$full), "sigma_full_zero")

finalize_run(run, extra = list(
  folds_used   = as.list(used_folds),
  n_fits       = sum(jobs$fit),
  n_stays      = length(stay_ids),
  pc1_full_zero = bp$pc1[bp$model == "full" & bp$fill == "zero" &
                           bp$gcs == "all three"][1]
))
cat(sprintf("\n  run directory: %s\n", run$path))
