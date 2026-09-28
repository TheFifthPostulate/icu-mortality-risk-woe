# tests/concurvity_null.R ----------------------------------------------------
# Calibrate concurvity against a matched null, and save the table.
#
# THE PROBLEM. `mgcv::concurvity(full = FALSE)` does not return 0 for unrelated
# covariates on this model class. It largely reads the ATOM STRUCTURE of the
# columns -- above all the point mass at zero that every intervention column
# carries on unexposed stays. So docs/v2_findings_20260827.md SS5's headline
# ("205 of 215 fits flagged, observed floor 0.9146") is mostly a statement about
# how many stays were never treated, and the fixed 0.80 threshold in config was
# mis-specified for this model class.
#
# THREE SECTIONS, in increasing realism:
#
#   A  synthetic control     pure noise, one covariate, nothing else in the
#                            model. There is NOTHING to concurve with. Whatever
#                            comes back is the diagnostic reading its own input.
#   B  shared support        two INDEPENDENT noise covariates zeroed on the SAME
#                            rows. Nothing links them but the support.
#   C  permutation null      the real models. Each smooth covariate's non-zero
#                            values are permuted among the rows that hold them,
#                            preserving the zero pattern and the marginal
#                            exactly, and destroying only the association.
#
# TWO NULLS, BECAUSE ONE IS NOT VALID FOR BOTH MEASURES.
#
#   worst      -> the PERMUTATION null (C). `worst` is a property of the model
#                 matrix alone, so permuting changes exactly what the null is
#                 about and nothing else. Exact.
#   observed   -> the SYNTHETIC null (B), matched on the model's dominant atom.
#                 The permutation null is INVALID here: permuting a covariate
#                 destroys its relationship to the OUTCOME as well, so its
#                 fitted smooth collapses to flat and cannot overlap anything --
#                 MEASURED, the null mean falls to ~0.015 and the excess becomes
#                 nearly the whole observed value. No permutation preserves the
#                 outcome association while destroying the covariate one, so
#                 this cannot be patched. The synthetic null is approximate --
#                 indexed by one number, so it cannot reproduce a second atom --
#                 but it does not flatten the smooth. Every output row records
#                 which null it used.
#
# WHAT IT WRITES. A `cnull_<datetime>` run directory holding
# `concurvity_null` (one row per signal x model x measure) and
# `concurvity_null_pairs` (every pair). Downstream reporting JOINS against this
# table rather than recomputing: at one refit per permutation this is a
# calibration to be done once, never a per-fit diagnostic.
#
# AGGREGATES ONLY (hard rule 1): concurvity values, means, percentiles, counts.
#
# --- WHAT IS CALIBRATED, AND WHAT IS SKIPPED (2026-09-07) -------------------
#
# THE INTERACTION SPECS ARE IN, AND THEY NEEDED NO WIRING. The loop walks
# `models_of(sg, cfg)`, so `full_ti_trend` and `full_ti_all` entered the moment
# they entered `LAYER1_MODELS`; `signal_frame()` builds their frames from the
# same columns and `attr(d, "formula")` carries their cross terms. What did need
# fixing is the OPPOSITE problem.
#
# ALIASED SPECS ARE NOW SKIPPED, and 7 of the 10 were being calibrated twice
# BEFORE the interaction models existed. `models_of()` enumerates what a signal
# COULD have; `spec_source()` says which of those are actually fitted. An
# unpaired signal's `full` IS its `meas` -- the identical formula on the
# identical rows -- so calibrating both spent B permutations per smooth
# covariate to produce a second copy of a row already in the table. That is 29
# permuted covariates times B refits of wasted compute per run, and worse, it
# put a duplicate row into a table downstream code JOINS against, where a
# many-to-one join silently doubles rows.
#
# The interaction models added three more: `full_ti_trend` on creatinine,
# platelet and hemoglobin, whose class carries no `trend` covariate, so their
# cross-term set is empty and their formula is `full`'s.
#
# MEASURED: 74 specs enumerated, 10 aliased, 64 distinct. 293 permuted
# covariates at B refits each.
#
#   Rscript tests/concurvity_null.R              # all 64 specs, B = 20
#   Rscript tests/concurvity_null.R 10           # B = 10, half the refits
#   Rscript tests/concurvity_null.R 20 mbp,creatinine,hemoglobin   # a subset
#
# THE PLAN IS PRINTED BEFORE ANY FIT and `--plan-only` stops there. At B = 20
# this run is 5,860 refits, of which 2,480 are tensor models an order of
# magnitude more expensive than the additive ones, so "how long will this take"
# is a question the script must answer before it starts rather than after.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
PLAN_ONLY <- "--plan-only" %in% args
args    <- setdiff(args, "--plan-only")
B       <- if (length(args) >= 1L) as.integer(args[1]) else 20L
only_sg <- if (length(args) >= 2L) strsplit(args[2], ",")[[1]] else NULL

# --- config and the plan FIRST, data second ----------------------------------
# THE PLAN IS PRINTED BEFORE THE DATA IS TOUCHED, and that ordering is the point
# rather than tidiness. The first version put it after `load_tables()` and
# `layer1_priors()`, so `--plan-only` -- whose entire purpose is to answer "how
# long will this take" before committing -- spent five minutes loading a million
# rows of signal features and fitting 172 prior models in order to print an
# arithmetic projection that needs only `cfg`. That is the same defect
# `--check-specs` had in `tests/attr_replicates.R` (plan section 37): a cheap
# verification step whose only exercise is the expensive path.
cfg <- load_config("config/config.yml")
signals <- if (is.null(only_sg)) cfg$signals else intersect(only_sg, cfg$signals)

# --- the plan, before anything is fitted -------------------------------------
# `spec_source()` is the SAME function `layer1_jobs()` and `l_matrix()` read, so
# what this script calibrates is exactly what the pipeline fits. It is not a
# second list that agrees with the pipeline today.
SPECS <- do.call(rbind, lapply(signals, function(sg)
  do.call(rbind, lapply(models_of(sg, cfg), function(md) {
    if (!is.na(spec_source(sg, md, cfg))) return(NULL)   # aliased: skip
    f  <- build_formula(sg, md, cfg)
    sm <- smooth_specs(f)
    data.frame(signal = sg, model = md,
               n_terms = length(attr(stats::terms(f), "term.labels")),
               n_perm_vars = length(unique(sm$variable)),
               stringsAsFactors = FALSE)
  }))))
.enum <- sum(vapply(signals, function(sg) length(models_of(sg, cfg)), integer(1)))
cat("\n=== plan ===\n")
cat(sprintf("  signals                     : %d\n", length(signals)))
cat(sprintf("  specs enumerated            : %d\n", .enum))
cat(sprintf("  aliased, skipped            : %d\n", .enum - nrow(SPECS)))
cat(sprintf("  distinct specs calibrated   : %d\n", nrow(SPECS)))
cat(sprintf("  permuted covariates         : %d\n", sum(SPECS$n_perm_vars)))
cat(sprintf("  refits at B = %-2d            : %d\n", B, B * sum(SPECS$n_perm_vars)))
.ti <- SPECS$model %in% LAYER1_TI_MODELS
cat(sprintf("    of which tensor models    : %d (over %d spec(s))\n",
            B * sum(SPECS$n_perm_vars[.ti]), sum(.ti)))
cat("\n  A tensor refit is roughly an order of magnitude dearer than an\n")
cat("  additive one, so the tensor line dominates the wall clock.\n\n")
if (PLAN_ONLY) quit(save = "no", status = 0L)

tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
tr     <- folds$stay_id[folds$split == "train"]
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
y_all  <- tabs$cohort$mortality[match(tr, tabs$cohort$stay_id)]


run <- new_run("cnull", cfg, note = sprintf(
  "concurvity permutation null, B = %d, %d signal(s), %d distinct spec(s)",
  B, length(signals), nrow(SPECS)))
save_table(run, SPECS, "concurvity_null_plan", subdir = "diagnostics")

# ---------------------------------------------------------------------------
# A. synthetic control: one smooth, one covariate, PURE NOISE
# ---------------------------------------------------------------------------
# x is Uniform(0,1) with a share of rows set to 0. It is independent of the
# outcome and of everything else, and the model has no second covariate. The
# only other thing present is the intercept. A diagnostic measuring a
# relationship between covariates must return ~0 at every row of this table.
cat("\n=== A. synthetic control: nothing to concurve with ===\n\n")
conc1 <- function(d, form) {
  b <- mgcv::bam(form, data = d, family = stats::binomial(), discrete = TRUE,
                 method = "fREML", na.action = stats::na.fail)
  cc <- mgcv::concurvity(b, full = FALSE)
  w <- cc$worst; o <- cc$observed
  diag(w) <- NA; diag(o) <- NA
  c(worst = max(w, na.rm = TRUE), observed = max(o, na.rm = TRUE))
}
cat(sprintf("%10s %10s %10s\n", "%zero", "worst", "observed"))
synth <- list()
for (p in c(0, 0.1, 0.3, 0.5, 0.68, 0.81, 0.87, 0.95)) {
  x <- with_seed(cfg$seed, {
    v <- stats::runif(length(y_all)); v[stats::runif(length(y_all)) < p] <- 0; v
  })
  r <- conc1(data.frame(mortality = y_all, x = x),
             mortality ~ s(x, bs = "ts", k = 10))
  synth[[length(synth) + 1L]] <- data.frame(section = "A", zero_frac = p,
                                            worst = r[["worst"]], observed = r[["observed"]],
                                            stringsAsFactors = FALSE)
  cat(sprintf("%9.0f%% %10.4f %10.4f\n", 100 * p, r[["worst"]], r[["observed"]]))
}

# ---------------------------------------------------------------------------
# B. shared support: two INDEPENDENT covariates, one shared zero block
# ---------------------------------------------------------------------------
# x1 and x2 are independent draws, both zeroed on the SAME rows. Nothing links
# them except that they are constant together on the zero block. This is the
# situation in every `intv` model, where the exposure and the intensity term are
# both zero on every unexposed stay -- and it is the part no reparameterisation
# of the exposed portion can remove.
cat("\n=== B. two independent covariates sharing one zero block ===\n\n")
cat(sprintf("%10s %11s %10s %10s\n", "%zero", "cor(x1,x2)", "worst", "observed"))
for (p in c(0, 0.3, 0.5, 0.68, 0.81, 0.87)) {
  d <- with_seed(cfg$seed + 1L, {
    z  <- stats::runif(length(y_all)) < p
    x1 <- stats::runif(length(y_all)); x1[z] <- 0
    x2 <- stats::runif(length(y_all)); x2[z] <- 0
    data.frame(mortality = y_all, x1 = x1, x2 = x2)
  })
  r <- conc1(d, mortality ~ s(x1, bs = "ts", k = 10) + s(x2, bs = "ts", k = 10))
  synth[[length(synth) + 1L]] <- data.frame(section = "B", zero_frac = p,
                                            worst = r[["worst"]], observed = r[["observed"]],
                                            stringsAsFactors = FALSE)
  cat(sprintf("%9.0f%% %11.4f %10.4f %10.4f\n", 100 * p,
              stats::cor(d$x1, d$x2), r[["worst"]], r[["observed"]]))
}
SYNTH <- do.call(rbind, synth)
save_table(run, SYNTH, "concurvity_synthetic", subdir = "diagnostics")

# ---------------------------------------------------------------------------
# C. the permutation null on the real models
# ---------------------------------------------------------------------------
# Each smooth covariate's non-zero values are permuted among the rows that hold
# them. The zero pattern and the marginal distribution survive exactly; only the
# row-wise association with the other covariates is destroyed. `excess` is the
# observed value minus its own null mean, and `above` says whether the observed
# value clears the null's 97.5th percentile.
cat(sprintf("\n=== C. permutation null on the real models (B = %d per covariate) ===\n", B))
cat("    excess = observed - null mean. `above` = observed clears the null's 97.5%%.\n\n")
cat(sprintf("%-26s %-9s %6s %9s %9s %9s %-11s %6s\n",
            "model", "measure", "atom", "observed", "null", "excess",
            "null source", "above"))

rows <- list(); pairs <- list()
t0 <- Sys.time()
# ITERATES `SPECS`, NOT `models_of()`. See the header: an aliased spec is the
# identical formula on the identical rows and calibrating it produces a
# duplicate row in a table downstream code joins against.
for (.i in seq_len(nrow(SPECS))) {
  {
    sg <- SPECS$signal[.i]; md <- SPECS$model[.i]
    pri <- priors_for(priors, sg, "final")
    d <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = tr)
    f <- attr(d, "formula")
    nt <- concurvity_null(f, d, cfg, B = B, seed = cfg$seed)
    if (!nrow(nt)) next
    nt$signal <- sg; nt$model <- md
    pairs[[length(pairs) + 1L]] <- nt

    # The synthetic null replaces the permutation null wherever the latter is
    # not valid, and the source is recorded ON THE ROW rather than left to the
    # reader to remember. See the header.
    at <- dominant_atom(f, d)
    ex <- concurvity_excess(nt, sg, md)
    ex$atom_frac   <- at$frac
    ex$atom_var    <- at$variable
    ex$null_source <- ifelse(ex$null_valid, "permutation", "synthetic")
    syn <- vapply(seq_len(nrow(ex)), function(i)
      synthetic_null_at(SYNTH, at$frac, ex$measure[i]), numeric(1))
    ex$null_used   <- ifelse(ex$null_valid, ex$null_mean, syn)
    ex$excess_used <- ex$observed - ex$null_used
    rows[[length(rows) + 1L]] <- ex
    for (i in seq_len(nrow(ex))) {
      cat(sprintf("%-26s %-9s %5.1f%% %9.4f %9.4f %+9.4f %-11s %6s\n",
                  paste0(sg, "/", md), ex$measure[i], 100 * ex$atom_frac[i],
                  ex$observed[i], ex$null_used[i], ex$excess_used[i],
                  ex$null_source[i],
                  if (ex$excess_used[i] > 0) "yes" else "no"))
    }
  }
}
E <- do.call(rbind, rows)
P <- do.call(rbind, pairs)
save_table(run, E, "concurvity_null", subdir = "diagnostics")
save_table(run, P, "concurvity_null_pairs", subdir = "diagnostics")

cat(sprintf("\n%d specs calibrated in %.1f min\n", length(unique(paste(E$signal, E$model))),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (m in unique(E$measure)) {
  z <- E[E$measure == m, , drop = FALSE]
  cat(sprintf("  %-9s [%s null]  observed %.4f-%.4f (median %.4f)  ->  excess %+.4f-%+.4f (median %+.4f)\n",
              m, z$null_source[1], min(z$observed), max(z$observed),
              stats::median(z$observed), min(z$excess_used), max(z$excess_used),
              stats::median(z$excess_used)))
}
cat("\n  The `worst` row has an EXACT null and is the one to quote. The `observed`\n")
cat("  row uses the approximate synthetic null and must be labelled as such.\n")
finalize_run(run, extra = list(B = B, n_specs = length(unique(paste(E$signal, E$model)))))
cat(sprintf("\nwritten: %s\n", run$path))
