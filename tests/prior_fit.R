# tests/prior_fit.R ----------------------------------------------------------
# Goodness-of-fit for the three fitted, outcome-blind parameter sets: the
# Dirichlet-multinomial alpha behind `pi_hat`, the conditional means behind
# `delta`, and the intensity parameters behind `lambda`.
#
# WHY THIS IS A SEPARATE SCRIPT AND NOT PART OF branch_point.R. Three reasons,
# and the first is the one that decides it.
#
#   1. IT DOES NOT NEED THE GAMs. Everything here reads `layer1_priors()` and
#      the feature tables. Nothing depends on a fitted `bam()` object, so
#      folding it into `branch_point.R` would tie a ~30-second diagnostic to a
#      6.8-minute critical-path run and slow every re-run of the thing the whole
#      project waits on.
#   2. IT IS RE-RUN ON A DIFFERENT TRIGGER. These numbers change when alpha,
#      delta or lambda change — a config edit, a re-extraction — not when a
#      smooth changes. Coupling them to the fit would recompute 215 GAMs to
#      answer a question about a shrinkage weight.
#   3. RUN DIRECTORIES SHOULD SAY WHAT THEY ARE. `priorfit_<datetime>` is
#      inspectable on its own terms; burying these tables inside a `branch_*`
#      directory makes them findable only by someone who already knows.
#
# It is also outcome-blind end to end. `mortality` is never read here, and the
# script can be run as often as you like without touching the single test look
# or contaminating anything.
#
# WHAT THE OUTPUT IS FOR. Misspecification in these three cannot invalidate a
# result — they never see the outcome, so a wrong likelihood cannot manufacture
# discrimination. It can cost power, which the XGBoost ladder already bounds at
# 0.0072 AUROC via `xgb_feat`, and it can break transport, which is the real
# exposure because every parameter here is fitted at MIMIC and frozen into the
# bundle. These tables are the methods-appendix evidence that both were checked.
#
# AGGREGATES ONLY (hard rule 1). Every table is one row per signal, category,
# variable or intervention.
#
#   Rscript tests/prior_fit.R              # MIMIC, the default
#   Rscript tests/prior_fit.R --selftest   # no data; checks the diagnostics
#                                            themselves against known-good and
#                                            known-bad synthetic counts
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args     <- commandArgs(trailingOnly = TRUE)
selftest <- "--selftest" %in% args

# ---------------------------------------------------------------------------
# SELF-TEST. A diagnostic that cannot fail is not a diagnostic, so before the
# tables are believed on real data the checks are run against counts drawn from
# the assumed model (they must pass) and from a deliberately wrong one (they
# must fail). Needs no data and no config.
# ---------------------------------------------------------------------------
if (selftest) {
  cat("\n=== self-test: does the PPC actually detect misspecification? ===\n\n")
  set.seed(11)
  N <- 8000L
  n <- sample(1:24, N, replace = TRUE)
  alpha_true <- c(0.8, 3.0, 1.4)

  rdm <- function(n, alpha) {
    t(vapply(seq_along(n), function(i) {
      p <- stats::rgamma(length(alpha), alpha); p <- p / sum(p)
      stats::rmultinom(1, n[i], p)[, 1]
    }, numeric(length(alpha))))
  }

  # (a) counts really are Dirichlet-multinomial: dispersion ~ 1, PIT sd ~ 1
  k_ok <- rdm(n, alpha_true)
  fit  <- fit_dm_alpha(k_ok)
  cat("  (a) well-specified. true alpha", paste(alpha_true, collapse = "/"),
      "-> fitted", paste(round(fit$alpha, 3), collapse = "/"), "\n")
  print(dm_ppc(k_ok, n, fit$alpha), row.names = FALSE)

  # (b) zero-inflated: a third of stays are forced to k_low = 0. The DM has no
  #     zero-inflation component, so p0_obs must exceed p0_exp and the PIT must
  #     stop being standard normal. This is the failure mode the check exists
  #     for, because k = 0 is the modal cell for several real signals.
  k_zi <- k_ok
  hit <- sample.int(N, N %/% 3)
  k_zi[hit, 2] <- k_zi[hit, 2] + k_zi[hit, 1]
  k_zi[hit, 1] <- 0
  fit_zi <- fit_dm_alpha(k_zi)
  cat("\n  (b) zero-inflated in the LOW category (a third of stays forced to 0)\n")
  ppc_ok <- dm_ppc(k_ok, n, fit$alpha)
  ppc_zi <- dm_ppc(k_zi, n, fit_zi$alpha)
  print(ppc_zi, row.names = FALSE)

  # An explicit verdict, because "PASS if the numbers look different" is not a
  # test. The comparison is row `low` under (b) against row `low` under (a).
  a <- ppc_ok[ppc_ok$category == "low", ]; b <- ppc_zi[ppc_zi$category == "low", ]
  vg <- function(lab, x, y, better) cat(sprintf(
    "    %-26s well-specified %+.4f -> misspecified %+.4f   %s\n", lab, x, y,
    if (better) "DETECTED" else "*** NOT DETECTED ***"))
  cat("\n  verdict on the LOW category:\n")
  vg("|dispersion - 1|", abs(a$dispersion - 1), abs(b$dispersion - 1),
     abs(b$dispersion - 1) > abs(a$dispersion - 1))
  vg("PIT KS statistic",  a$pit_ks, b$pit_ks, b$pit_ks > a$pit_ks)
  vg("p0_obs - p0_exp",   a$p0_obs - a$p0_exp, b$p0_obs - b$p0_exp,
     (b$p0_obs - b$p0_exp) > (a$p0_obs - a$p0_exp))
  cat("\n  NOTE, and it matters for reading the real table: refitting alpha ABSORBS\n")
  cat("  much of the zero inflation -- alpha_low falls from ~0.80 to ~0.41 here --\n")
  cat("  so the boundary-mass gap stays small even when the model is clearly wrong.\n")
  cat("  The KS statistic is the sensitive detector (it roughly quintuples); the\n")
  cat("  p0 comparison alone would under-call the problem. Weight them that way.\n")

  # (c) delta's form check must fire on a curved truth. The fitted mean is
  #     linear in log1p(k), so a quadratic truth has to leave structure behind.
  cat("\n=== self-test: does the delta form check detect curvature? ===\n\n")
  kk <- rpois(4000, 4); nn <- pmax(rpois(4000, 12), 1L)
  lk <- log1p(kk)
  for (lab in c("linear truth", "quadratic truth")) {
    v <- if (lab == "linear truth") 2 + 1.5 * lk + 0.1 * log(nn) + rnorm(4000, sd = 0.5)
         else                       2 + 1.5 * lk - 0.6 * lk^2 + 0.1 * log(nn) + rnorm(4000, sd = 0.5)
    g <- stats::lm(v ~ lk + log(nn)); cf <- stats::coef(g)
    r <- v - (cf[1] + cf[2] * lk + cf[3] * log(nn))
    r2 <- summary(mgcv::gam(r ~ s(lk, bs = "ts", k = 6)))$r.sq
    cat(sprintf("  %-16s resid_smooth_r2 = %+.4f\n", lab, r2))
  }
  cat("\n  PASS if the linear truth is ~0 and the quadratic truth is clearly above it.\n\n")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# The real run.
# ---------------------------------------------------------------------------
cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)

run <- new_run("priorfit", cfg, note = "goodness-of-fit for pi_hat / delta / lambda")

t0 <- Sys.time()
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
log_msg(run, sprintf("layer1_priors: %.1f s | %d alpha fits, %d delta fits, %d lambda fits",
                     as.numeric(difftime(Sys.time(), t0, units = "secs")),
                     nrow(priors$signal), nrow(priors$magnitude),
                     nrow(priors$intervention)))

s <- .measured_train_rows(tabs, folds, cfg)
log_msg(run, sprintf("measured training rows: %d over %d signals",
                     nrow(s), length(unique(s$signal))))

# --- 1. shrinkage weights ---------------------------------------------------
shrink <- dm_shrinkage_table(priors$signal, s, role = "final")
save_table(run, shrink, "dm_shrinkage")

cat("\n=== pi_hat: prior strength against observed coverage ===\n\n")
print(shrink[, c("signal", "n_stays", "alpha0", "n_q25", "n_med", "n_q75",
                 "w_q25", "w_med", "w_q75")], row.names = FALSE)
cat("\n  `w_med` is the posterior weight on the OBSERVED proportion at median\n")
cat("  coverage. Near 0 means the prior dominates and pi_hat is nearly constant.\n")
cat("  Near 1 means no shrinkage, so 0-of-2 and 0-of-24 are the same evidence --\n")
cat("  the failure the Dirichlet-multinomial layer exists to prevent. There is\n")
cat("  deliberately NO threshold here; this is a table to read.\n")

# --- 2. is the DM the right model? -----------------------------------------
ppc <- dm_ppc_all(priors$signal, s, role = "final", seed = cfg$seed)
save_table(run, ppc, "dm_ppc")

cat("\n=== pi_hat: posterior-predictive check, DM against the counts ===\n\n")
print(ppc, row.names = FALSE)
cat("\n  dispersion  1.00 = correct. ABOVE 1 = data more variable than the DM\n")
cat("              allows, so alpha is pulled to a spurious consensus. BELOW 1\n")
cat("              = alpha smaller than the data support, so pi_hat under-shrinks.\n")
cat("  p0_obs vs p0_exp   the zero-inflation check. The DM has no zero-inflation\n")
cat("              component, so an excess of exact zeros is the one shape it\n")
cat("              cannot absorb -- and k = 0 is the modal cell for several signals.\n")
cat("  pit_sd      1.00 under correct specification. Uses the EXACT beta-binomial\n")
cat("              CDF (the DM's marginal), not a simulation.\n")

worst <- ppc[order(-abs(ppc$dispersion - 1)), ][1:min(8L, nrow(ppc)), ]
cat("\n  furthest from dispersion 1:\n\n")
print(worst[, c("signal", "category", "dispersion", "p0_obs", "p0_exp", "pit_sd")],
      row.names = FALSE)

# --- 3. delta and lambda ----------------------------------------------------
dd <- delta_fit_diagnostics(priors$magnitude, s, cfg, role = "final")
if (!is.null(dd)) {
  save_table(run, dd, "delta_fit")
  cat("\n=== delta: is the conditional mean's FORM adequate? ===\n\n")
  print(dd, row.names = FALSE)
  cat("\n  `resid_smooth_r2` is a smooth of the residual on log1p(k). Near 0 means\n")
  cat("  the deliberately rigid linear form captured the relationship. Materially\n")
  cat("  above 0 means it did not, and the leftover is exactly the count\n")
  cat("  dependence delta exists to remove. This is a DIFFERENT statement from\n")
  cat("  the decoupling in conditional_priors.md 7.1: that measures the finished\n")
  cat("  covariate, this measures the model that produced it.\n")
} else {
  cat("\n=== delta: no magnitude priors fitted (construct off in config) ===\n")
}

ll <- lambda_fit_diagnostics(priors$intervention, role = "final")
if (!is.null(ll)) {
  save_table(run, ll, "lambda_fit")
  cat("\n=== lambda: the fitted intensity parameters ===\n\n")
  print(ll, row.names = FALSE)
  cat("\n  `c1` (binomial) and `b1` (lognormal) are two-parameter summaries of a\n")
  cat("  population's treatment behaviour. v2_state 6 identifies them as the\n")
  cat("  cleanest cross-site quantity in the project, because they are directly\n")
  cat("  comparable between sites and independent of per-term attribution.\n")
  cat("  Run this script at eICU and difference the table: that IS the lambda\n")
  cat("  transport result.\n")
}

# --- 4. fold-to-fold stability ---------------------------------------------
stab <- prior_fold_stability(priors)
if (!is.null(stab)) {
  save_table(run, stab, "prior_fold_stability")
  cat("\n=== how much does each frozen parameter move between folds? ===\n\n")
  print(utils::head(stab, 20), row.names = FALSE)
  cat("\n  Coefficient of variation over the five out-of-fold fits. Each is fitted\n")
  cat("  on overlapping four-fifths of ONE site, so this is a LOWER BOUND on how\n")
  cat("  much the parameter would move under a genuinely different population.\n")
  cat("  A parameter that already wanders between folds will not survive a change\n")
  cat("  of site -- and that is knowable now, before eICU runs.\n")
}

finalize_run(run, extra = list(
  n_signals        = length(unique(s$signal)),
  n_measured_rows  = nrow(s),
  n_alpha_fits     = nrow(priors$signal),
  n_delta_fits     = nrow(priors$magnitude),
  n_lambda_fits    = nrow(priors$intervention),
  worst_dispersion = if (nrow(ppc)) max(abs(ppc$dispersion - 1), na.rm = TRUE) else NA_real_,
  max_delta_resid_r2 = if (!is.null(dd)) max(dd$resid_smooth_r2, na.rm = TRUE) else NA_real_))
cat(sprintf("\n  run directory: %s\n", run$path))
