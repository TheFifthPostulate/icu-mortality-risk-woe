# tests/term_redundancy.R ----------------------------------------------------
# Is a covariate earning its place, or is it a copy of the one beside it?
#
# Not a pass/fail test — a measurement, to be read before a term set is frozen.
# Answers three questions the concurvity diagnostics raise but cannot settle,
# because concurvity says "these two overlap" and not "which one to keep":
#
#   A  Do the level terms carry the same number?      value_median vs q05/q95
#   B  Does the second level term add anything?       nested dev.expl
#   C  Extensive or intensive margin?                 ever_active vs intensity
#
# Aggregates only, never a row (hard rule 1): correlations, percentages, counts,
# and deviance explained. Nothing here prints a measurement or an id.
#
# Rscript tests/term_redundancy.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
tr    <- folds$stay_id[folds$split == "train"]
sf    <- tabs$signal_features
coh   <- tabs$cohort

dev_expl <- function(f, d) {
  b <- try(bam(f, data = d, family = binomial(), discrete = TRUE, method = "fREML"),
           silent = TRUE)
  if (inherits(b, "try-error")) return(NA_real_)
  (b$null.deviance - b$deviance) / b$null.deviance
}

# --- A. are the three level terms the same number? --------------------------
# With n_obs = 1 the median, the 5th and the 95th percentile are all the single
# observed value: three columns, one number. With n_obs = 2 the quantiles are
# the two draws and the median is between them. The percentage below is the
# direct measure of how much of the level block is duplication.
cat("\n=== A. level terms: how often is value_median LITERALLY q05? ===\n\n")
cat(sprintf("%-18s %-7s %6s %7s %7s %10s %11s %11s\n", "signal", "class",
            "p50_n", "%n==1", "%n<=2", "%med==q05", "r(med,q05)", "r(med,q95)"))
for (sg in cfg$signals) {
  keep <- sf$signal == sg & sf$n_obs > 0 & sf$stay_id %in% tr
  d <- sf[keep, c("n_obs", "value_median", "q05", "q95"), drop = FALSE]
  rr <- function(a, b) if (stats::sd(d[[a]]) == 0) NA_real_ else stats::cor(d[[a]], d[[b]])
  cat(sprintf("%-18s %-7s %6.0f %6.1f%% %6.1f%% %9.1f%% %11.3f %11.3f\n",
              sg, signal_class_of(sg, cfg), stats::median(d$n_obs),
              100 * mean(d$n_obs == 1), 100 * mean(d$n_obs <= 2),
              100 * mean(d$value_median == d$q05),
              rr("value_median", "q05"), rr("value_median", "q95")))
}

# --- B. does the second level term buy anything? ----------------------------
cat("\n=== B. nested fits: does the tail add to the median, or repeat it? ===\n")
cat("    gain = dev.expl(med + tail) - max(dev.expl(med), dev.expl(tail))\n\n")
cat(sprintf("%-18s %-7s %9s %9s %9s %9s  %s\n", "signal", "class",
            "med", "tail", "med+tail", "gain", "tail"))
for (sg in cfg$signals) {
  keep <- sf$signal == sg & sf$n_obs > 0 & sf$stay_id %in% tr
  d <- sf[keep, c("value_median", "q05", "q95", "value_min", "value_max"), drop = FALSE]
  d$y <- coh$mortality[match(sf$stay_id[keep], coh$stay_id)]
  # Follows config `level_terms`, so section B measures the term the pipeline
  # ACTUALLY fits. Section A above deliberately keeps comparing q05 to the
  # minimum — that comparison is what justifies the declaration in the first
  # place, and it has to stay independent of it.
  side <- excursion_side_of(sg, cfg)
  lv <- level_vars_of(sg, cfg)
  tv <- if (is.na(side)) lv[2] else if (side == "low") lv[1] else lv[2]
  d$tail <- d[[tv]]
  k1 <- smooth_k_of(sg, "value_median", cfg); k2 <- smooth_k_of(sg, tv, cfg)
  fm <- function(s) stats::as.formula(s)
  a  <- dev_expl(fm(sprintf("y ~ s(value_median, bs='ts', k=%d)", k1)), d)
  b  <- dev_expl(fm(sprintf("y ~ s(tail, bs='ts', k=%d)", k2)), d)
  ab <- dev_expl(fm(sprintf("y ~ s(value_median, bs='ts', k=%d) + s(tail, bs='ts', k=%d)",
                            k1, k2)), d)
  cat(sprintf("%-18s %-7s %9.5f %9.5f %9.5f %9.5f  %s\n",
              sg, signal_class_of(sg, cfg), a, b, ab, ab - max(a, b), tv))
}

# --- C. extensive vs intensive margin ---------------------------------------
# `ever_active` is exactly (intensity > 0) on 100.0000% of stays, so the two are
# not independent covariates: the question is only which one carries the signal.
#   gain_int  what the intensity smooth adds ON TOP OF the indicator
#   gain_ind  what the indicator adds ON TOP OF the smooth  <- the decisive one
cat("\n=== C. ever_active vs the intensity smooth ===\n\n")
iv <- tabs$intervention_features
shape <- unlist(cfg$intervention_shape)
cat(sprintf("%-21s %-6s %7s %9s %9s %9s %9s %9s\n", "intervention", "shape",
            "%expo", "ever", "s(int)", "both", "gain_int", "gain_ind"))
for (nm in cfg$interventions_modelled) {
  z <- iv[iv$intervention == nm, , drop = FALSE]
  z <- z[z$stay_id %in% tr, , drop = FALSE]
  if (!nrow(z)) next
  z$int <- if (shape[[nm]] == "state") as.numeric(z$exposure_frac) else as.numeric(z$n_hours)
  z$ea  <- as.numeric(z$ever_active)
  z$y   <- coh$mortality[match(z$stay_id, coh$stay_id)]
  if (anyNA(z$y) || anyNA(z$int) || mean(z$ea) < 0.005) next
  agree <- mean(z$ea == as.integer(z$int > 0))
  ku <- min(10L, max(3L, length(unique(z$int)) - 1L))
  a  <- dev_expl(y ~ ea, z)
  b  <- dev_expl(stats::as.formula(sprintf("y ~ s(int, bs='ts', k=%d)", ku)), z)
  ab <- dev_expl(stats::as.formula(sprintf("y ~ ea + s(int, bs='ts', k=%d)", ku)), z)
  if (!isTRUE(all.equal(agree, 1))) {
    cat(sprintf("  NOTE %s: ever_active != (intensity > 0) on %.4f%% of stays\n",
                nm, 100 * (1 - agree)))
  }
  cat(sprintf("%-21s %-6s %6.1f%% %9.5f %9.5f %9.5f %9.5f %9.5f\n",
              nm, shape[[nm]], 100 * mean(z$ea == 1), a, b, ab, ab - a, ab - b))
}
cat("\n  gain_ind ~ 0 means ever_active is redundant given the intensity smooth:\n")
cat("  the smooth already separates the point mass at zero from everything else.\n")
