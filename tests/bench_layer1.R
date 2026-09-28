# tests/bench_layer1.R -------------------------------------------------------
# Step 1 of docs/next_steps.md: time a single bam() fit before budgeting the
# full 258-fit run, and exercise fit_one() end to end on real data.
#
# Fits the largest and smallest specs plus one of each model kind, and runs one
# out-of-fold job so the prediction and L path is covered too — that path is
# never touched by a `final` job, which predicts on nothing.
#
# Aggregates only, never a row (hard rule 1). Everything printed is a count, a
# rate, or a fitted statistic.
#
# Rscript tests/bench_layer1.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

cfg <- load_config("config/config.yml")
cat(sprintf("bam: discrete=%s method=%s nthreads=%d k=%d bs=%s\n",
            cfg$bam$discrete, cfg$bam$method, cfg$bam$nthreads,
            cfg$bam$k, cfg$bam$smooth_basis))

tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)

# One of each shape. gcs_verbal/full is the widest spec in the set; the intv
# models are the cheapest and are what the +72 fits actually cost.
jobs <- list(
  list(sg = "gcs_verbal",        md = "full", role = "final", fold = NA_integer_),
  list(sg = "gcs_verbal",        md = "intv", role = "final", fold = NA_integer_),
  list(sg = "creatinine",        md = "meas", role = "final", fold = NA_integer_),
  list(sg = "mbp",               md = "full", role = "final", fold = NA_integer_),
  list(sg = "temperature",       md = "meas", role = "final", fold = NA_integer_),
  # the out-of-fold path: fits on 4 folds, predicts on the 5th, produces L
  list(sg = "mbp",               md = "full", role = "oof",   fold = 1L),
  list(sg = "mbp",               md = "intv", role = "oof",   fold = 1L)
)

rows <- list(); diags <- list()
for (j in jobs) {
  ids <- job_ids(j$role, j$fold, folds)
  pri <- priors_for(priors, j$sg, j$role, j$fold)
  t <- system.time(
    r <- fit_one(j$sg, j$md, tabs, cfg, pri,
                 fit_ids = ids$fit_ids, predict_ids = ids$predict_ids,
                 keep_model = FALSE, role = j$role, fold = j$fold)
  )
  d <- r$diagnostics
  diags[[length(diags) + 1L]] <- d
  rows[[length(rows) + 1L]] <- data.frame(
    signal = j$sg, model = j$md, role = j$role,
    secs = round(unname(t[["elapsed"]]), 2),
    n_rows = d$n_rows, n_terms = d$n_terms, n_smooth = d$n_smooth,
    edf = d$edf_total, dev_expl = d$dev_expl, converged = d$converged,
    k_index = d$k_index_min, concurv = d$concurvity_max,
    n_pred = d$n_predict, l_mean = d$l_mean, l_sd = d$l_sd,
    stringsAsFactors = FALSE)
}

bench <- do.call(rbind, rows)
cat("\n")
print(bench, row.names = FALSE)

tot <- sum(bench$secs)
cat(sprintf("\n%d fits in %.1f s (mean %.2f s/fit)\n", nrow(bench), tot, tot / nrow(bench)))
# DERIVED, NOT A LITERAL. It was `258` twice on these two lines, which was the
# budget until 2026-09-07 and is 384 now -- a projection that quietly under-
# reports by a third is worse than no projection. This is the one place in the
# project where deriving is right and asserting is wrong: `tests/smoke.R` states
# the budget as literals ON PURPOSE, because there it IS the check.
.nfit <- layer1_budget(cfg)$fits_total
cat(sprintf("serial projection for %d fits at the mean: %.1f min\n",
            .nfit, .nfit * (tot / nrow(bench)) / 60))

# Triage the sample, so a threshold that fires on every model shows up now
# rather than after the full run. A threshold nothing can pass is a threshold
# chosen wrong, and that is much cheaper to learn here than at fit 258.
cat("\n")
dg <- do.call(rbind, diags)
tr <- diagnostics_summary(dg, cfg)
if (nrow(tr)) {
  for (i in seq_len(nrow(tr))) {
    cat(sprintf("  %-18s %-5s  %s\n", tr$signal[i], tr$model[i], tr$flags[i]))
  }
}
