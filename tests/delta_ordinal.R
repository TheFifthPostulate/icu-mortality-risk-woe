# tests/delta_ordinal.R ------------------------------------------------------
# The ordinal `delta` form: does the estimator work, and does it fix what it was
# introduced to fix?
#
# WHY THIS EXISTS AS A SCRIPT RATHER THAN A ONE-OFF. The ordinal form was
# adopted on the strength of a measurement, and the measurement has to be
# re-runnable -- at MIMIC after any change to the estimator, and at eICU, where
# the frozen cut points meet a Glasgow distribution that is known to differ
# (gcs_verbal sits at its floor for 56.3% of eICU ventilated measurements
# against MIMIC's 13.0%).
#
# SIX SECTIONS. A-D are self-contained and need no data; E-G read the training
# rows.
#
#   A  parameter recovery on synthetic data with a known truth
#   B  the mid-PIT is uniform when the model is correct
#   C  the packed string round-trips exactly (the bundle path)
#   D  an off-scale value ERRORS rather than clamping (the eICU path)
#   E  linear versus ordinal on the real Glasgow rows
#   F  the proportional-odds check, and the two named escalations
#   G  calibration, uniformity and cut points per fitted row
#
# AGGREGATES ONLY (hard rule 1): coefficients, R-squareds, counts, distances.
#
#   Rscript tests/delta_ordinal.R            # everything
#   Rscript tests/delta_ordinal.R --synth    # A-D only, no data read
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args      <- commandArgs(trailingOnly = TRUE)
synth_only <- "--synth" %in% args
ok <- TRUE
say <- function(pass, fmt, ...) {
  ok <<- ok && isTRUE(pass)
  cat(sprintf(paste0("  [%s] ", fmt, "\n"), if (isTRUE(pass)) "ok" else "FAIL", ...))
}

# --- A. recovery ------------------------------------------------------------
cat("\n=== A. parameter recovery on synthetic data ===\n")
set.seed(11)
N <- 20000
n_s <- pmax(1, rpois(N, 12))
k_s <- rbinom(N, n_s, 0.12)
th_true <- c(-2.2, -1.4, -0.3, 0.9, 2.1)
a1_true <- -1.7; a2_true <- 0.4
pr <- .polr_probs(th_true, a1_true * log1p(k_s) + a2_true * log(n_s))
yc <- apply(pr, 1, function(p) sample.int(6L, 1L, prob = p))

p <- fit_delta_ordinal(yc, k_s, n_s, list(min = 1, max = 6))
say(p$converged, "converged")
say(abs(p$a1 - a1_true) < 0.1, "a1  true %+.3f  fitted %+.3f", a1_true, p$a1)
say(abs(p$a2 - a2_true) < 0.1, "a2  true %+.3f  fitted %+.3f", a2_true, p$a2)
say(max(abs(p$theta - th_true)) < 0.15,
    "max |theta error| %.4f over %d cut points", max(abs(p$theta - th_true)),
    length(th_true))

# --- B. uniformity ----------------------------------------------------------
cat("\n=== B. the mid-PIT is uniform when the model is correct ===\n")
eta <- p$a1 * log1p(k_s) + p$a2 * log(n_s) + p$a3 * as.numeric(k_s > 0)
u <- plogis(.mid_pit_logit(p$theta, eta, yc))
say(abs(mean(u) - 0.5) < 0.02, "mean(u) %.4f (want 0.500)", mean(u))
say(abs(sd(u) - sqrt(1/12)) < 0.02, "sd(u) %.4f (want %.4f)", sd(u), sqrt(1/12))
ks <- suppressWarnings(ks.test(u, "punif")$statistic)
say(ks < 0.08, "KS distance from uniform %.4f (discreteness keeps this above 0)", ks)

# --- C. round trip ----------------------------------------------------------
cat("\n=== C. the packed string form round-trips exactly ===\n")
d1 <- delta_value_ordinal(yc, k_s, n_s, p)
pp <- p; pp$theta <- .pack_num(p$theta); pp$levels <- .pack_num(p$levels)
d2 <- delta_value_ordinal(yc, k_s, n_s, pp)
say(identical(d1, d2) || isTRUE(all.equal(d1, d2, tolerance = 0)),
    "packed == unpacked, max diff %.3g", max(abs(d1 - d2)))
say(identical(.unpack_num(.pack_num(p$theta)), p$theta),
    "%%.17g encoding is lossless for theta")

# --- D. the off-scale guard -------------------------------------------------
cat("\n=== D. a value off the frozen scale ERRORS rather than clamping ===\n")
r <- tryCatch({delta_value_ordinal(c(yc[1:10], 9L), c(k_s[1:10], 1L),
                                    c(n_s[1:10], 5), p); FALSE},
              error = function(e) TRUE)
say(r, "refused a value outside the declared scale")
r2 <- tryCatch({fit_delta(1:10, 0:9, rep(5, 10), form = "ordinal", scale = NULL); FALSE},
               error = function(e) TRUE)
say(r2, "refused the ordinal form with no declared scale")

if (synth_only) {
  cat(sprintf("\n%s\n\n", if (ok) "A-D pass." else "FAILURES above."))
  quit(save = "no", status = if (ok) 0L else 1L)
}

# --- E. the real rows -------------------------------------------------------
cat("\n=== E. linear versus ordinal on the Glasgow rows ===\n")
cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
s     <- .measured_train_rows(tabs, folds, cfg)

ord_signals <- Filter(function(sg) identical(delta_form_of(sg, cfg), "ordinal"),
                      unlist(cfg$signals))
if (!length(ord_signals)) {
  cat("  no signal declares the ordinal form; nothing to compare.\n")
} else {
  cmp <- do.call(rbind, lapply(ord_signals, function(sg) {
    z  <- s[s$signal == sg, ]
    sc <- ordinal_scale_of(sg, cfg)
    lv <- .ordinal_levels(sc)
    do.call(rbind, lapply(level_vars_of(sg, cfg), function(vr) {
      cv <- delta_count_of(vr)
      v <- z[[vr]]; k <- z[[cv]]; nn <- z$n_obs
      lk <- log1p(k)
      r2 <- function(rr) if (sd(lk) == 0) NA_real_ else
        summary(mgcv::gam(rr ~ s(lk, bs = "ts", k = 6)))$r.sq
      po <- fit_delta(v, k, nn, form = "ordinal", scale = sc)
      pl <- fit_delta(v, k, nn, form = "linear")
      fv <- pl$a0 + pl$a1 * lk + pl$a2 * log(nn)
      data.frame(signal = sg, variable = vr, n = length(v),
                 lin_r2 = round(r2(.delta_residual(pl, v, k, nn)), 4),
                 ord_r2 = round(r2(.delta_residual(po, v, k, nn)), 4),
                 lin_offscale = round(mean(fv < min(lv) | fv > max(lv)), 4),
                 ord_a1 = round(po$a1, 3), ord_a3 = round(po$a3, 3),
                 ord_conv = po$converged, stringsAsFactors = FALSE)
    }))
  }))
  print(cmp, row.names = FALSE)
  # The pre-specified acceptance criterion (docs/v2_external_plan_20260901).
  tgt <- cmp[!is.na(cmp$ord_r2), ]
  say(all(tgt$ord_r2 < 0.10),
      "every ordinal row has resid_smooth_r2 < 0.10 (worst %.4f)", max(tgt$ord_r2))
  say(all(tgt$ord_r2 <= tgt$lin_r2 + 1e-9),
      "the ordinal form is no worse than the linear one on every row")
}

# --- F. proportional odds and the escalations -------------------------------
cat("\n=== F. proportional odds: does the log1p(k) slope move across cut points? ===\n")
for (sg in ord_signals) {
  z <- s[s$signal == sg, ]
  lv <- .ordinal_levels(ordinal_scale_of(sg, cfg))
  lk <- log1p(z$k_low); ln <- log(z$n_obs)
  sl <- vapply(lv[-length(lv)], function(cc) {
    yb <- as.integer(z$value_min <= cc)
    if (length(unique(yb)) < 2L) return(NA_real_)
    g <- try(suppressWarnings(glm(yb ~ lk + ln, family = binomial())), silent = TRUE)
    if (inherits(g, "try-error")) NA_real_ else unname(coef(g)[2])
  }, numeric(1))
  cat(sprintf("  %-11s cut-point slopes %s  spread %.2f\n", sg,
              paste(sprintf("%+7.2f", sl), collapse = " "),
              diff(range(sl, na.rm = TRUE))))
}
cat("\n  A large spread would say to escalate to per-cut-point slopes. MEASURED\n")
cat("  2026-09-01: that escalation LOSES to the two-regime latent predictor on\n")
cat("  both residual structure and AIC, which is why `a3` exists and partial\n")
cat("  proportional odds does not. Re-read this if the spread ever grows.\n")

# --- G. per-row calibration and uniformity ----------------------------------
cat("\n=== G. calibration, uniformity and cut points, per fitted row ===\n")
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
od <- delta_ordinal_diagnostics(priors$magnitude, s, cfg, role = "final", seed = cfg$seed)
if (is.null(od)) cat("  no ordinal rows fitted.\n") else {
  print(od[, c("signal", "variable", "n_levels", "a1", "a3", "cal_max_abs_diff",
               "pit_ks", "pit_sd", "po_slope_spread", "levels_unobserved",
               "converged")], row.names = FALSE)
  say(all(od$converged), "every ordinal fit converged")
  say(all(od$cal_max_abs_diff < 0.02),
      "worst category calibration gap %.4f of stays", max(od$cal_max_abs_diff))
  cat("\n  theta (the fitted cut points), which are what travel to eICU:\n")
  for (i in seq_len(nrow(od))) {
    cat(sprintf("    %-11s %-10s %s\n", od$signal[i], od$variable[i], od$theta[i]))
  }
}

cat(sprintf("\n%s\n\n", if (ok) "delta_ordinal: all checks pass." else
            "delta_ordinal: FAILURES above."))
quit(save = "no", status = if (ok) 0L else 1L)
