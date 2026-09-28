# tests/refactor_identity.R --------------------------------------------------
# The threading refactor must change NO NUMBER.
#
# `signal_frame()` and `fit_one()` stopped taking a bare `alpha` vector and now
# take the container from `priors_for()`. That is a pure refactor: with both
# conditional-prior constructs switched off, every formula, every model frame
# and every L value must be what the pre-refactor code produced.
#
# This is the one check that can prove it, and it can only be run while a
# pre-refactor run directory still exists. Fold 1 is enough -- 43 fits over
# ~332,000 L values -- because a threading bug would not be selective about
# which fold it corrupted.
#
# It doubles as the ablation baseline: if the constructs are ever switched off
# in config for a real run, this is the assertion that the run reproduces the
# published pre-construct numbers exactly rather than approximately.
#
# ============================================================================
# THE REFERENCE IS STALE AS OF 2026-09-03, AND FOR A DATA REASON, NOT A CODE
# REASON. READ THIS BEFORE BELIEVING A FAILURE.
# ============================================================================
# Both feature tables were re-extracted on 2026-09-02 to add the `se_` replicate
# statistics. `v2_05_features_*.sql` builds `q05`, `q95` and `value_median` with
# BigQuery's `APPROX_QUANTILES`, which is approximate, so a re-extraction can
# move those columns slightly even with no change to the query that produces
# them. Every stored `l_oof.rds` in a `branch_*` directory predates that
# re-extraction, so this test now compares current code against numbers computed
# from a DIFFERENT feature table and will report differences that are not code
# differences.
#
# WHAT DID NOT CAUSE IT. The `delta` variance-component estimator changed on the
# same day (`.fit_var_components_rep()`, see `R/04b_conditional.R`), but that
# cannot affect this test: it runs with both conditional constructs switched
# OFF, and with `magnitude_conditional` off no delta is fitted at all. The
# estimator swap was separately verified to leave every conditional-mean
# coefficient identical on 38 of 38 rows.
#
# WHAT "GIVE IT A NEW REFERENCE" MEANS, AND WHAT IT COSTS. Switch both
# constructs off in config, run the branch pipeline once against the CURRENT
# parquets, and keep that `l_oof.rds` as the reference. About seven minutes.
#
# BE CLEAR ABOUT WHAT THAT BUYS AND WHAT IT SPENDS. The original assertion was
# historical: it proved that the `priors_for()` threading refactor changed no
# number relative to the PRE-REFACTOR implementation. That proof was made and
# passed when it mattered, and regenerating the reference does not re-prove it
# -- the new reference is produced by the current code, so the test becomes a
# FORWARD regression guard rather than a proof about the past. That is still
# worth having, and it is a different thing. Do not describe a regenerated
# reference as evidence that the refactor was sound; cite the original run for
# that.
# ============================================================================
#
# AGGREGATES ONLY (hard rule 1): row counts and a maximum absolute difference.
# The L values themselves are row-level and never printed.
#
#   Rscript tests/refactor_identity.R [run_dir]     ~2 min
#     run_dir  a branch_* directory holding tables/l_oof.rds.
#              Default: the newest completed one.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) args[1] else latest_run("branch")
ref_path <- file.path(run_dir, "tables", "l_oof.rds")
if (!file.exists(ref_path)) {
  stop("no l_oof.rds under ", run_dir, ". Pass a branch_* run directory that has one.",
       call. = FALSE)
}
cat(sprintf("reference: %s\n", ref_path))

# BOTH CONSTRUCTS OFF. This is the configuration the reference run was fitted
# under, so anything that differs is the refactor and not the design.
cfg <- load_config("config/config.yml")
# RE-SCOPED 2026-09-03. The constructs-off ablation is no longer runnable and
# is no longer the useful comparison, for two independent reasons.
#
# IT DOES NOT RUN. With `intensity_conditional` off, the raw intervention
# intensity covariates enter in place of the lambda-transformed ones, and
# `config/smooth_k` was calibrated for the constructs-ON design. Since
# `inotrope` was paired to `mbp` on 2026-09-01 the ablation dies on
# `inotrope__n_agents`: k = 10 declared against 3 distinct values. Declaring a
# `smooth_k` entry for a covariate no live formula contains, purely to keep a
# vestigial path fittable, would put a fiction in the design config.
#
# IT IS NOT WHAT NEEDS GUARDING. The constructs are permanent design now, so the
# constructs-off path is exercised by nothing except this test. The threading it
# was written to protect -- `priors_for()` into `signal_frame()` and
# `fit_one()` -- is exercised by every real run, so the useful guard is against
# the PUBLISHED design rather than against an ablation of it.
#
# WHAT WAS SPENT AND WHAT WAS KEPT. The original assertion was historical: it
# proved the threading refactor changed no number relative to the PRE-REFACTOR
# implementation. That proof was made and passed when it mattered, and it is NOT
# re-proved here -- the reference below was produced by the current code, so
# this is now a FORWARD regression guard. Cite the original run for the refactor
# claim, never this one.
#
# THE REFERENCE. `out/runs/branch_20260903T085658`, fold 1, 43 fits, the
# published design against the 2026-09-02 re-extraction. Regenerate it with
# `Rscript tests/branch_point.R 1` after any deliberate change to layer 1, and
# say in the commit why the old reference stopped being valid.

tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)

jobs <- layer1_jobs(cfg)
jobs <- jobs[jobs$role == "oof" & jobs$fold == 1L, , drop = FALSE]
cat(sprintf("refitting %d jobs on fold 1 under the PUBLISHED design\n", sum(jobs$fit)))

r   <- run_layer1(tabs, cfg, folds, priors, jobs = jobs, keep_final = FALSE, verbose = FALSE)
new <- r$l
old <- readRDS(ref_path)
old <- old[old$fold == 1L, , drop = FALSE]

key <- function(d) paste(d$signal, d$model, d$stay_id)
m <- match(key(new), key(old))

pass <- TRUE
say <- function(label, ok, note = "") {
  pass <<- pass && isTRUE(ok)
  cat(sprintf("  [%s] %-46s %s\n", if (isTRUE(ok)) "ok" else "FAIL", label, note))
}

say("every new L row matches a reference row", !anyNA(m),
    sprintf("%d rows, %d unmatched", nrow(new), sum(is.na(m))))
say("row counts agree", nrow(new) == nrow(old),
    sprintf("new %d, reference %d", nrow(new), nrow(old)))

if (!anyNA(m)) {
  d <- new$l - old$l[m]
  say("L is bitwise identical", isTRUE(all.equal(new$l, old$l[m], tolerance = 0)),
      sprintf("max |difference| = %.3e", max(abs(d))))
}

cat(sprintf("\nrefactor_identity: %s\n", if (pass) "PASSED" else "FAILED"))
if (!pass) quit(status = 1L)
