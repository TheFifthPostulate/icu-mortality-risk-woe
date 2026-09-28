# tests/conditional_review.R -------------------------------------------------
# The design measurements behind R/04b_conditional.R. Outcome-blind throughout:
# nothing here reads `mortality`, so every number is a statement about how the
# covariates are distributed, not about how they predict.
#
# Not a pass/fail test -- a measurement, to be read before either construct is
# frozen. Four sections:
#
#   A  the coupling            how strongly is intensity determined by exposure,
#                              and what shape is it? This is why a FITTED model
#                              is needed rather than a chosen transform.
#   B  the naive rate          what a plain A/D would do instead. Read beside A.
#   C  count gradation         is the deviation count really graded, or is it
#                              n_obs x 1{abnormal}? A REJECTED hypothesis -- see
#                              the note at the head of the section.
#   D  decoupling              did the constructs work? Correlation and, because
#                              correlation only sees straight lines, the R2 of a
#                              smooth of the exposure.
#
# Aggregates only, never a row (hard rule 1): correlations, fitted coefficients,
# R2, counts and percentages.
#
#   Rscript tests/conditional_review.R          ~1 min
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
tr    <- folds$stay_id[folds$split == "train"]
sf    <- tabs$signal_features
ivf   <- tabs$intervention_features

# The (exposure, accumulation) pairs, read off the design rather than listed.
iv_pairs <- local({
  ivs <- unique(unlist(lapply(cfg$signals, interventions_of, cfg = cfg)))
  out <- list()
  for (iv in ivs) {
    ls <- lambda_spec_of(iv, cfg)
    if (!is.null(ls)) out[[iv]] <- ls
  }
  out
})
if (!length(iv_pairs)) {
  message("intensity_conditional is off; sections A, B and D(intensity) are empty.")
}

# ---------------------------------------------------------------------------
# A. the coupling: how is accumulation determined by exposure?
# ---------------------------------------------------------------------------
# beta is the slope of log A on log D among EXPOSED stays:
#   beta = 0   A does not depend on D at all
#   beta = 1   A is exactly proportional to D  (steady accumulation)
#   beta > 1   A rises faster than D           (escalation)
# The point of the table is the SPREAD of beta. A single arithmetic transform
# can be correct for at most one value of beta, and the observed range is wide.
cat("\n=== A. the (exposure D, accumulation A) coupling, exposed stays ===\n\n")
cat(sprintf("%-22s %-14s %-14s %8s %8s %7s %7s\n",
            "intervention", "D", "A", "n_exp", "r(D,A)", "beta", "R2"))
for (iv in names(iv_pairs)) {
  ls <- iv_pairs[[iv]]
  z <- ivf[ivf$intervention == iv & ivf$stay_id %in% tr, , drop = FALSE]
  D <- z[[ls$exposure]] * ls$exposure_scale
  A <- z[[ls$accumulation]]
  ok <- !is.na(D) & !is.na(A) & D > 0 & A > 0
  D <- D[ok]; A <- A[ok]
  g <- stats::lm(log(A) ~ log(D))
  cat(sprintf("%-22s %-14s %-14s %8d %8.3f %7.3f %7.3f\n",
              iv, ls$exposure, ls$accumulation, length(D), stats::cor(D, A),
              stats::coef(g)[2], summary(g)$r.squared))
}

# ---------------------------------------------------------------------------
# B. why the naive rate is worse than doing nothing
# ---------------------------------------------------------------------------
# A/D removes the coupling only when beta == 1. Elsewhere it inherits a NEW
# dependence on D with the opposite sign, and where beta ~ 0 it manufactures one
# out of nothing. Same failure as the `spread = |tail - median|` composite
# already rejected in CLAUDE.md.
cat("\n=== B. what a naive rate A/D would do instead ===\n\n")
cat(sprintf("%-22s %10s %10s %10s\n", "intervention", "r(A,D)", "r(A/D,D)", "verdict"))
for (iv in names(iv_pairs)) {
  ls <- iv_pairs[[iv]]
  z <- ivf[ivf$intervention == iv & ivf$stay_id %in% tr, , drop = FALSE]
  D <- z[[ls$exposure]] * ls$exposure_scale
  A <- z[[ls$accumulation]]
  ok <- !is.na(D) & !is.na(A) & D > 0 & A > 0
  D <- D[ok]; A <- A[ok]
  r0 <- stats::cor(A, D); r1 <- stats::cor(A / D, D)
  cat(sprintf("%-22s %10.3f %10.3f %10s\n", iv, r0, r1,
              if (abs(r1) < abs(r0)) "helps" else "WORSE"))
}

# ---------------------------------------------------------------------------
# C. count gradation -- a hypothesis this table REJECTS
# ---------------------------------------------------------------------------
# The hypothesis was that `delta` should lose deviance exactly where the
# deviation count is not really graded -- where k is close to
# n_obs x 1{abnormal}, so conditioning on it strips the extensive margin out of
# the magnitude term and leaves only a weak within-abnormal contrast.
#
# It would have given an OUTCOME-BLIND rule for per-signal exclusion, which is
# the only kind of rule that may be used: selecting signals on their deviance
# change is selection on the response, and CLAUDE.md forbids exactly that move
# for diagnostic thresholds.
#
# MEASURED 2026-08-28, AND IT DOES NOT SEPARATE THE CASES. The two signals that
# lose the most deviance under `delta` are platelet (-45.0%) and spo2 (-18.3%),
# and they sit at OPPOSITE ends of every statistic here: platelet is saturated
# (sat 0.85, psat 69%, bin/full 1.06) and spo2 is not (sat 0.26, psat 0.4%,
# bin/full 0.70). Meanwhile the three GCS components have the most binary counts
# in the table (bin/full 1.26-1.32, R2grad 0.02-0.11) and are unaffected
# (-1.2%, +0.4%, -0.9%), and bilirubin_total is the most saturated signal
# present (sat 0.95, psat 89%) and GAINS 23.7%.
#
# The hypothesis is kept here, with its refutation, rather than deleted: it is
# the obvious mechanism to propose, and the table is what stops it being
# proposed again. There is currently NO validated outcome-blind criterion for
# excluding a signal from the magnitude construct.
cat("\n=== C. is k essentially n_obs x 1{abnormal}? (hypothesis: REJECTED) ===\n")
cat("  sat       mean k/n among k>0; 1.0 = every observation was deviant\n")
cat("  psat      %% of k>0 stays with k == n\n")
cat("  R2full    R2 of magnitude ~ log1p(k) + log(n), all measured stays\n")
cat("  R2bin     same with 1{k>0} in place of log1p(k)  -- the EXTENSIVE margin\n")
cat("  R2grad    R2full refitted on k>0 stays only      -- the GRADED part\n")
cat("  bin/full  share of the standardisation that is only the on/off split\n\n")
cat(sprintf("%-18s %-5s %6s %6s %6s %7s %7s %7s %8s\n",
            "signal", "tail", "%k=0", "sat", "psat", "R2full", "R2bin", "R2grad", "bin/full"))
s <- sf[sf$n_obs > 0 & sf$stay_id %in% tr, , drop = FALSE]
.r2 <- function(y, X) {
  if (stats::sd(y) == 0) return(NA_real_)
  summary(stats::lm(y ~ ., data = data.frame(y = y, X)))$r.squared
}
for (sg in cfg$signals) {
  z <- s[s$signal == sg, , drop = FALSE]
  lv <- level_vars_of(sg, cfg); side <- excursion_side_of(sg, cfg)
  for (v in if (is.na(side)) lv else if (side == "low") lv[1] else lv[2]) {
    kv <- delta_count_of(v)
    k <- z[[kv]]; n <- z$n_obs; y <- z[[v]]; ex <- k > 0
    R2f <- .r2(y, data.frame(lk = log1p(k), ln = log(n)))
    R2b <- .r2(y, data.frame(b = as.numeric(ex), ln = log(n)))
    cat(sprintf("%-18s %-5s %5.1f%% %6.3f %5.1f%% %7.3f %7.3f %7.3f %8.3f\n",
                sg, sub("k_", "", kv), 100 * mean(!ex),
                if (any(ex)) mean(k[ex] / n[ex]) else NA_real_,
                100 * (if (any(ex)) mean(k[ex] == n[ex]) else NA_real_),
                R2f, R2b,
                if (sum(ex) > 30) .r2(y[ex], data.frame(lk = log1p(k[ex]), ln = log(n[ex]))) else NA_real_,
                R2b / R2f))
  }
}

# ---------------------------------------------------------------------------
# D. did the constructs decouple what they were built to decouple?
# ---------------------------------------------------------------------------
# Correlation is reported because it is the natural summary, and the smooth R2
# beside it because correlation only sees straight lines: a covariate can be a
# perfect U-shaped function of another and still correlate at zero. The smooth
# R2 is the honest claim.
cat("\n=== D. decoupling: raw pair vs conditional coordinate ===\n")
cat("    r_new/R2 are computed on the INFORMATIVE subset (exposed, or measured),\n")
cat("    never on the full frame: a shared zero block inflates both by itself.\n\n")
pri <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
.ds <- decoupling_smooth_spec(cfg)
smooth_r2 <- function(y, x) {
  if (stats::sd(y) == 0 || stats::sd(x) == 0) return(NA_real_)
  # k must stay under the distinct-value count, exactly as in config/smooth_k:
  # some exposure columns are small counts and cannot carry the default basis.
  # The cap and the basis are DECLARED (`diagnostics.decoupling_smooth`, added
  # 2026-09-03, audit finding F11); only the reduction to the distinct-value
  # count is computed here, because that is a property of the data.
  k <- min(.ds$k_max, length(unique(x)) - 1L)
  if (k < 3L) return(NA_real_)
  f <- stats::as.formula(sprintf("y ~ s(x, bs = \"%s\", k = %d)", .ds$basis, k))
  b <- mgcv::bam(f, data = data.frame(y = y, x = x),
                 discrete = TRUE, method = "fREML")
  summary(b)$r.sq
}
cat(sprintf("%-38s %9s %9s %9s\n", "pair", "r_raw", "r_new", "smoothR2"))

for (sg in cfg$signals) {
  if (!magnitude_conditional_for(sg, cfg)) next
  p <- priors_for(pri, sg, "final")
  d <- signal_frame(sg, "meas", tabs, cfg, p, stay_ids = tr)
  z <- sf[sf$signal == sg, ]; m <- match(d$stay_id, z$stay_id)
  for (v in grep("_delta$", names(d), value = TRUE)) {
    base <- delta_base_of(v)
    pv <- if (delta_count_of(base) == "k_low") "pi_minus" else "pi_plus"
    if (!pv %in% names(d)) next
    cat(sprintf("%-38s %9.3f %9.3f %9.3f\n",
                sprintf("%s: %s ~ %s", sg, base, pv),
                stats::cor(z[[base]][m], d[[pv]]), stats::cor(d[[v]], d[[pv]]),
                smooth_r2(d[[v]], d[[pv]])))
  }
}
for (sg in cfg$signals) {
  for (iv in interventions_of(sg, cfg)) {
    ls <- lambda_spec_of(iv, cfg)
    if (is.null(ls)) next
    p <- priors_for(pri, sg, "final")
    d <- signal_frame(sg, "intv", tabs, cfg, p, stay_ids = tr)
    z <- ivf[ivf$intervention == iv, ]; m <- match(d$stay_id, z$stay_id)
    D <- z[[ls$exposure]][m]; A <- z[[ls$accumulation]][m]
    ex <- D > 0
    cat(sprintf("%-38s %9.3f %9.3f %9.3f\n",
                sprintf("%s: %s ~ %s", sg, ls$accumulation, ls$exposure),
                stats::cor(A[ex], D[ex]),
                stats::cor(d[[paste0(iv, "__lambda")]][ex], D[ex]),
                smooth_r2(d[[paste0(iv, "__lambda")]][ex], D[ex])))
    break   # one row per signal; the pair is a property of the intervention
  }
}
cat("\n")
