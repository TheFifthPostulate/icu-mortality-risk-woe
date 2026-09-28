# tests/coupling_interaction.R -----------------------------------------------
# DOES THE MEASUREMENT-INTERVENTION COUPLING HYPOTHESIS HOLD?
#
# The design's motivating claim is that measurement and intervention form a
# JOINT stochastic block -- that treatment context changes what a physiological
# deviation MEANS, not merely that it predicts death on its own. In a GAM that
# is an INTERACTION claim, and layer 1 cannot express it: it is additive by
# frozen decision, with no `by=` smooths and no tensor terms since `o_flag` left
# on 2026-08-25.
#
# So the arm-level comparison of `llr_cond` against `llr_meas` does NOT test the
# hypothesis. Under additivity
#
#   L_full = f(M) + g(I) - logit(p_bar)
#   L_intv =        g~(I) - logit(p_bar)
#   L_cond = f(M) + [ g(I) - g~(I) ]
#
# so if M and I were independent, g and g~ would coincide, `L_cond` would equal
# `L_meas` EXACTLY, and the comparison would be a tautology returning zero. The
# near-equality that was observed therefore measures how much the intervention
# block moves when measurements are added, and says nothing about whether f
# differs across treatment strata. See
# docs/measurement_intervention_coupling_20260905.md section 1.
#
# THIS SCRIPT RUNS THE TEST THAT WAS MISSING. For each paired signal it adds a
# tensor interaction `ti(measurement_covariate, intervention_covariate)` to the
# fitted `full` model and asks whether it buys anything.
#
# WHY `ti()` AND NOT `te()` OR `by=`. `ti()` is the interaction with the main
# effects EXCLUDED and marginally centred, so the additive model is exactly
# nested inside it and the test is a clean question about the interaction alone.
# It is also why this does not reproduce the near-1 concurvity that
# `s(trend, by = o_flag)` drove in every paired dense model: a `by=` smooth
# carries the main effect again, a `ti()` does not.
#
# THIS DOES NOT REOPEN THE FROZEN NO-INTERACTION DECISION. Nothing here reaches
# an L, a score, a bundle or a target. It is diagnostic-only, which is the same
# standing exception `R/04c_prior_diagnostics.R`'s `.decoupling_gam()` holds --
# that smooth IS the diagnostic, and so is this one. It also respects the
# condition attached to the frozen decision, that an interaction brought back
# must be built on intervention INTENSITY rather than on excursion-relative
# timing: every intervention covariate crossed here is an intensity or a
# duration, and `first_hour` is refused by the formula builder regardless.
#
# WHAT THE TWO OUTCOMES MEAN, both useful and neither available now:
#
#   ti terms null across the board  the coupling hypothesis is dead in this
#                                   feature space. A clean negative that the
#                                   additive comparison could not establish.
#   ti terms not null               the frozen decision has a MEASURED cost, and
#                                   CLAUDE.md's methods caveat gets a number
#                                   instead of a hedge.
#
# Aggregates only, never a row (hard rule 1): edf, chi-square, p, deviance
# explained, AIC and counts. Nothing here prints a measurement or an id.
#
#   Rscript tests/coupling_interaction.R              all paired signals
#   Rscript tests/coupling_interaction.R mbp,spo2     just those
#   Rscript tests/coupling_interaction.R "" --no-joint   skip section B
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
only_sg <- if (length(args) >= 1L && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
do_joint <- !("--no-joint" %in% args)

cfg    <- load_config("config/config.yml")
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
tr     <- folds$stay_id[folds$split == "train"]
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)

# Paired signals only. An unpaired signal has no intervention block, so `full`
# aliases `meas` and there is no cross pair to form.
paired <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)
signals <- if (is.null(only_sg)) paired else intersect(only_sg, paired)

run <- new_run("coupling", cfg, note = sprintf(
  "measurement-intervention interaction probe, %d paired signal(s), diagnostic-only",
  length(signals)))

# --- helpers ----------------------------------------------------------------

#' Deviance explained, the quantity summary.gam() prints.
.de <- function(b) (b$null.deviance - b$deviance) / b$null.deviance

#' `.bam_fit()` with the smoothing parameters supplied.
#'
#' THIS EXISTS BECAUSE THE OBVIOUS COMPARISON IS NOT NESTED. `method = "fREML"`
#' with `bs = "ts"` estimates every smoothing parameter JOINTLY, so adding a
#' `ti()` term re-optimises the whole penalty vector and the larger model's
#' deviance is NOT bounded above by the smaller model's. MEASURED on the first
#' run of this script: `hemoglobin`'s `pi_minus x transfusion_prbc__n_hours`
#' pair reported a deviance-explained change of -0.00183, which is impossible
#' for a genuinely nested pair of maximum-likelihood fits and is entirely
#' possible for two penalised fits at different penalty optima.
#'
#' Freezing the additive model's `sp` at its own fitted values and estimating
#' only the interaction's makes the additive model exactly nested inside the
#' interaction model, so `d_dev_nested` below is monotone and readable as "what
#' the interaction adds". mgcv's convention is that a negative `sp` entry is
#' estimated, so the new terms are passed as -1.
#'
#' Every other setting still comes from `bam_settings(cfg)`, so this is the
#' pipeline's fit with one argument added and not a second fitting path.
.fit_sp <- function(f, d, cfg, sp) {
  s <- bam_settings(cfg)
  mgcv::bam(formula = f, data = d, family = s$family, method = s$method,
            discrete = s$discrete, nthreads = s$nthreads, select = s$select,
            gc.level = s$gc_level, na.action = stats::na.fail, sp = sp)
}

#' AIC with the smoothing-parameter-uncertainty correction.
#'
#' `AIC()` on a gam counts `sum(edf)`, which treats the smoothing parameters as
#' known and is anti-conservative for exactly the comparison this script makes:
#' the interaction model estimated more of them. `edf2` is mgcv's corrected
#' effective degrees of freedom, available because `method = "fREML"`, and is
#' the right count here. Both are reported so the gap between them is visible
#' rather than assumed small.
.aic2 <- function(b) {
  if (is.null(b$edf2)) return(NA_real_)
  -2 * as.numeric(stats::logLik(b)) + 2 * sum(b$edf2)
}

#' The row of summary()$s.table belonging to a term, by label prefix.
#'
#' Matched on the label mgcv itself assigns rather than on a reconstructed
#' string: `ti(a,b)` is printed without spaces and a hand-built key drifts from
#' it the moment mgcv changes its formatting.
.s_row <- function(b, prefix) {
  tb <- summary(b)$s.table
  i <- which(startsWith(rownames(tb), prefix))
  if (!length(i)) return(NULL)
  list(edf = sum(tb[i, 1]), chi2 = sum(tb[i, 3]), p = min(tb[i, 4]), n = length(i))
}

#' The interaction basis dimension for one margin.
#'
#' Capped well below the main-effect k on purpose. The question is whether an
#' interaction EXISTS, not what shape it has, and a tensor product costs the
#' product of its margins -- k = 5 on both sides is 25 coefficients per term,
#' k = 10 would be 100. The cap never exceeds config's declared `smooth_k`, so a
#' bounded discrete scale (the GCS triple, at 4 to 6 distinct values) still gets
#' a basis it can support and mgcv does not error out.
.k_ti <- function(sg, v, cfg, cap = 5L) max(3L, min(cap, smooth_k_of(sg, v, cfg)))

#' Split a fitted model's smooth covariates into the measurement and
#' intervention sets.
#'
#' Read off the FITTED object rather than re-derived from the term strings, so
#' the split cannot drift from what was actually estimated. Membership comes
#' from `measurement_terms()`, which is the same function the formula builder
#' uses, so a covariate cannot be silently assigned to the wrong block.
.blocks <- function(b, sg, cfg) {
  sv <- unique(vapply(b$smooth, function(s) s$term[1], character(1)))
  mv <- .term_vars(measurement_terms(sg, cfg))
  list(meas = intersect(sv, mv), intv = setdiff(sv, mv))
}

fmt_p <- function(p) if (is.na(p)) "     NA" else if (p < 1e-99) " <1e-99" else sprintf("%7.1e", p)

# ---------------------------------------------------------------------------
# A. one interaction at a time: every (measurement, intervention) cross pair
# ---------------------------------------------------------------------------
# The complete enumeration, so no pair is chosen after seeing a result. The
# number of tests is reported at the end and must be read with the p values --
# with roughly a hundred pairs, a handful below 0.05 is what the null predicts.
cat("\n=== A. one ti() at a time, added to the fitted `full` model ===\n\n")
cat("    %expo  share of fitted rows with a non-zero intervention value. The\n")
cat("           interaction is only identifiable where that value varies.\n")
cat("    edf/p  the ti term's effective df and mgcv's own p value. THE p VALUE\n")
cat("           IS ANTI-CONSERVATIVE HERE and is reported for triage only: it\n")
cat("           conditions on the estimated smoothing parameters as if they\n")
cat("           were known. Read d_AIC2 and d_dev_n instead.\n")
cat("    d_dev_n  deviance explained gained with the additive model's smoothing\n")
cat("           parameters FROZEN, so the two models are genuinely nested and\n")
cat("           this cannot be negative. The quantity to quote.\n")
cat("    d_dev_f  the same with all smoothing parameters re-estimated. CAN BE\n")
cat("           NEGATIVE, because the two fits sit at different penalty optima.\n")
cat("    d_AIC2 AIC change on the free fit, using mgcv's edf2 correction.\n")
cat("           Negative means the interaction is worth its parameters.\n\n")
cat(sprintf("%-16s %-19s %-28s %6s %5s %8s %9s %9s %8s\n",
            "signal", "measurement", "intervention", "%expo", "edf",
            "p", "d_dev_n", "d_dev_f", "d_AIC2"))

A <- list(); JOINT <- list()
t0 <- start_timer()

for (sg in signals) {
  pri <- priors_for(priors, sg, "final")
  d   <- signal_frame(sg, "full", tabs, cfg, pri, stay_ids = tr)
  f0  <- build_formula(sg, "full", cfg)
  b0  <- try(.bam_fit(f0, d, cfg), silent = TRUE)
  if (inherits(b0, "try-error")) {
    cat(sprintf("%-18s  BASELINE FAILED: %s\n", sg, conditionMessage(attr(b0, "condition"))))
    next
  }
  bl <- .blocks(b0, sg, cfg)
  de0 <- .de(b0); aic0 <- stats::AIC(b0); aic2_0 <- .aic2(b0)

  ti_terms <- character(0)
  for (mv in bl$meas) for (iv in bl$intv) {
    ti <- sprintf("ti(%s, %s, bs = c(\"ts\", \"ts\"), k = c(%d, %d))",
                  mv, iv, .k_ti(sg, mv, cfg), .k_ti(sg, iv, cfg))
    f1 <- stats::update(f0, stats::as.formula(paste(". ~ . +", ti)))
    b1 <- try(.bam_fit(f1, d, cfg), silent = TRUE)
    if (inherits(b1, "try-error")) {
      cat(sprintf("%-16s %-19s %-28s  FIT FAILED (free sp)\n", sg, mv, iv))
      next
    }
    # The genuinely nested comparison: the additive block's penalties frozen at
    # their additive-model values, only the interaction's estimated.
    n_new <- length(b1$sp) - length(b0$sp)
    b1n <- try(.fit_sp(f1, d, cfg, sp = c(unname(b0$sp), rep(-1, n_new))), silent = TRUE)
    de_n <- if (inherits(b1n, "try-error")) NA_real_ else .de(b1n) - de0

    ti_terms <- c(ti_terms, ti)
    r <- .s_row(b1, sprintf("ti(%s,%s)", mv, iv))
    expo <- mean(d[[iv]] != 0)
    row <- data.frame(
      signal = sg, meas_var = mv, intv_var = iv,
      frac_exposed = round(expo, 4),
      ti_edf = round(r$edf %||% NA_real_, 3),
      ti_chi2 = round(r$chi2 %||% NA_real_, 1),
      ti_p = r$p %||% NA_real_,
      dev_expl_additive = round(de0, 5),
      d_dev_nested = round(de_n, 5),
      d_dev_free = round(.de(b1) - de0, 5),
      d_aic = round(stats::AIC(b1) - aic0, 1),
      d_aic2 = round(.aic2(b1) - aic2_0, 1),
      n_sp_new = as.integer(n_new),
      n_rows = as.integer(stats::nobs(b1)),
      stringsAsFactors = FALSE)
    A[[length(A) + 1L]] <- row
    cat(sprintf("%-16s %-19s %-28s %5.1f%% %5.2f %8s %+9.5f %+9.5f %+8.1f\n",
                sg, mv, iv, 100 * expo, row$ti_edf, fmt_p(row$ti_p),
                row$d_dev_nested, row$d_dev_free, row$d_aic2))
  }

  # -------------------------------------------------------------------------
  # B. all cross interactions for this signal at once
  # -------------------------------------------------------------------------
  # The per-pair table above is many tests; this is one. It is the honest
  # headline for "does THIS signal show coupling", because it asks the question
  # once per signal rather than once per pair.
  if (do_joint && length(ti_terms)) {
    fj <- stats::update(f0, stats::as.formula(paste(". ~ . +", paste(ti_terms, collapse = " + "))))
    bj <- try(.bam_fit(fj, d, cfg), silent = TRUE)
    if (!inherits(bj, "try-error")) {
      rj <- .s_row(bj, "ti(")
      n_new <- length(bj$sp) - length(b0$sp)
      bjn <- try(.fit_sp(fj, d, cfg, sp = c(unname(b0$sp), rep(-1, n_new))), silent = TRUE)
      JOINT[[length(JOINT) + 1L]] <- data.frame(
        signal = sg, n_ti = length(ti_terms),
        ti_edf_total = round(rj$edf %||% NA_real_, 2),
        ti_chi2_total = round(rj$chi2 %||% NA_real_, 1),
        ti_p_min = rj$p %||% NA_real_,
        dev_expl_additive = round(de0, 5),
        d_dev_nested = round(if (inherits(bjn, "try-error")) NA_real_ else .de(bjn) - de0, 5),
        d_dev_free = round(.de(bj) - de0, 5),
        d_aic = round(stats::AIC(bj) - aic0, 1),
        d_aic2 = round(.aic2(bj) - aic2_0, 1),
        stringsAsFactors = FALSE)
    } else {
      cat(sprintf("%-16s  JOINT FIT FAILED\n", sg))
    }
  }
}

A <- do.call(rbind, A)
save_table(run, A, "coupling_interaction_pairs", subdir = "diagnostics")

if (length(JOINT)) {
  J <- do.call(rbind, JOINT)
  save_table(run, J, "coupling_interaction_joint", subdir = "diagnostics")
  cat("\n=== B. all cross interactions per signal, in one model ===\n\n")
  cat("    One test per signal rather than one per pair, so this is the honest\n")
  cat("    headline for `does THIS signal show coupling`.\n\n")
  cat(sprintf("%-18s %5s %8s %9s %10s %10s %9s\n",
              "signal", "n_ti", "edf", "p_min", "d_dev_n", "d_dev_f", "d_AIC2"))
  for (i in seq_len(nrow(J))) {
    cat(sprintf("%-18s %5d %8.2f %9s %+10.5f %+10.5f %+9.1f\n",
                J$signal[i], J$n_ti[i], J$ti_edf_total[i],
                fmt_p(J$ti_p_min[i]), J$d_dev_nested[i], J$d_dev_free[i], J$d_aic2[i]))
  }
}

# --- what it says -----------------------------------------------------------
cat("\n=== summary ===\n\n")
cat(sprintf("  %d cross pairs tested across %d paired signals, in %.1f minutes.\n",
            nrow(A), length(unique(A$signal)), t0()$elapsed_sec / 60))
cat(sprintf("  pairs with an interaction worth its parameters (d_AIC2 < 0): %d of %d\n",
            sum(A$d_aic2 < 0, na.rm = TRUE), nrow(A)))
cat(sprintf("  same on the uncorrected AIC, for comparison:                %d of %d\n",
            sum(A$d_aic < 0, na.rm = TRUE), nrow(A)))
cat(sprintf("  pairs at p < 0.05: %d  (the null predicts about %.0f, and the p\n",
            sum(A$ti_p < 0.05, na.rm = TRUE), 0.05 * nrow(A)))
cat("                     values are anti-conservative, so treat this as a ceiling)\n")
cat(sprintf("  largest nested deviance gain: %+.5f\n", max(A$d_dev_nested, na.rm = TRUE)))
cat(sprintf("  median nested deviance gain:  %+.5f\n",
            stats::median(A$d_dev_nested, na.rm = TRUE)))
cat(sprintf("  free-sp gains that came out NEGATIVE: %d of %d (they are not nested)\n",
            sum(A$d_dev_free < 0, na.rm = TRUE), nrow(A)))
cat("\n  READ THE GAINS AGAINST SECTION 6 OF THE COUPLING DOC. The nonlinear\n")
cat("  conditioning increment there is worth +0.0021 deviance explained and the\n")
cat("  whole propensity block +0.0101, both at the SCORE level. A per-signal\n")
cat("  interaction gain far below those is coupling that exists and does not\n")
cat("  matter; one comparable to them is coupling the additive design is losing.\n")

finalize_run(run, extra = list(n_pairs = nrow(A), n_signals = length(signals),
                               joint = do_joint))
cat(sprintf("\nwritten: %s\n", run$path))
