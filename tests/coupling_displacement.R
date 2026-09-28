# tests/coupling_displacement.R ----------------------------------------------
# HOW FAR DOES A MEASUREMENT SMOOTH MOVE WHEN THE INTERVENTION BLOCK IS PRESENT?
#
# The third instrument on the coupling question, and the one that measures the
# mechanism directly rather than its consequence for a score. For each paired
# signal it evaluates every measurement smooth twice -- once from the `meas`
# model and once from the `full` model -- on a shared grid, and reports the
# displacement in NATS on the log-odds scale, split into
#
#   level   the mean of f_full(x) - f_meas(x). A constant offset, which is
#           absorbed by the intercept and carries no information about ordering.
#   shape   the standard deviation of that difference after centring. THIS is
#           the quantity the coupling hypothesis is about: it is non-zero only
#           if the measurement-risk FUNCTION changed, not merely its level.
#
# THE COMPARISON NEEDS A NULL, AND NOT AN OBVIOUS ONE. Two things move a
# measurement smooth when four intervention terms are added, and neither is
# coupling:
#
#   1. `bam(discrete = TRUE)` builds a DIFFERENT discretisation grid for the
#      same covariate depending on what else is in the model. MEASURED
#      2026-09-05: for `mbp/s(q05_delta)`, 514 unique basis values in both fits
#      and up to 0.174 apart, on model frames that are bitwise identical. See
#      docs/measurement_intervention_coupling_20260905.md section 8, which
#      identifies the branch by elimination.
#   2. `method = "fREML"` selects every smoothing parameter JOINTLY, so adding
#      any four terms reallocates the penalty budget across all of them.
#
# THE NULL IS EXACTLY THE HYPOTHESIS'S NEGATION. Within each outcome class,
# permute the ROW INDEX of the whole intervention block at once, and apply that
# one permutation to every intervention column together. This preserves
#
#   - p(I | Y) EXACTLY, as an empirical joint distribution. Every intervention
#     covariate keeps its marginal, its point mass at zero, and its correlation
#     with the other intervention covariates, separately within survivors and
#     within deaths. So the `intv` block still estimates the same contrast and
#     the intervention smooths are still as predictive as they were.
#   - the measurement block, untouched.
#
# and destroys exactly one thing: the row-wise association between a patient's
# measurements and their own treatment. That is the null hypothesis
# "M and I are conditionally independent given Y", which is precisely what the
# joint-block hypothesis denies. Under it the null model has the same term
# count, the same discretisation behaviour and the same penalty-budget
# competition, so both artefacts above are subtracted rather than argued away.
#
# WHY NOT `.permute_preserving_zeros()` FROM R/08. That helper permutes ONE
# covariate's non-zero values, which is right for the concurvity null because
# concurvity is a pairwise property of the model matrix. Here it would break the
# functional relationship between `{iv}__exposure_frac` and `{iv}__lambda`,
# which are two views of the same treatment record, and the null would then
# differ from the observed fit in a second way. A block row permutation keeps
# them consistent. It also does not need the outcome-association caveat that
# helper carries, because permuting within outcome class preserves p(I | Y).
#
# COMPARED VIA predict(), NEVER VIA COEFFICIENTS. Finding 1 above means the two
# fits do not share a basis, so their coefficients are not commensurable. Every
# number here comes from `predict(type = "terms")` on a shared grid.
#
# THE STANDARD ERRORS CARRY THE SMOOTHNESS-UNCERTAINTY CORRECTION, and getting
# that took a workaround: `unconditional = TRUE` is a NO-OP on a `bam` object.
# See `.partial()` below for the mechanism, the measurement and why it also
# bears on the attributability arm.
#
# Aggregates only, never a row (hard rule 1): displacements in nats, standard
# errors, counts. The grid is built from covariate QUANTILES and the quantiles
# themselves are never printed.
#
#   Rscript tests/coupling_displacement.R              B = 10, all paired
#   Rscript tests/coupling_displacement.R 20           B = 20
#   Rscript tests/coupling_displacement.R 10 mbp,spo2  just those
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
B       <- if (length(args) >= 1L && nzchar(args[1])) as.integer(args[1]) else 10L
only_sg <- if (length(args) >= 2L && nzchar(args[2])) strsplit(args[2], ",")[[1]] else NULL
SEED    <- 1L
N_GRID  <- 199L

cfg    <- load_config("config/config.yml")
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
tr     <- folds$stay_id[folds$split == "train"]
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)

paired  <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)
signals <- if (is.null(only_sg)) paired else intersect(only_sg, paired)

run <- new_run("coupdisp", cfg, note = sprintf(
  "measurement-smooth displacement, meas vs full, B = %d, %d signal(s)",
  B, length(signals)))

# --- helpers ----------------------------------------------------------------

#' Permute the intervention block as one row-block, within outcome class.
#'
#' One permutation applied to every intervention column at once, drawn
#' separately inside `y == 0` and `y == 1`. See the header for why this and not
#' a per-column permutation.
.permute_block_within_class <- function(d, iv_cols, y, seed) {
  with_seed(seed, {
    for (cl in c(0L, 1L)) {
      i <- which(y == cl)
      if (length(i) < 2L) next
      j <- i[sample.int(length(i))]
      d[i, iv_cols] <- d[j, iv_cols]
    }
    d
  })
}

#' A one-row newdata frame with every model variable at a neutral value.
#'
#' `type = "terms"` evaluates each smooth from its OWN covariate only, so the
#' other columns never reach the number that is reported. They are present
#' because `predict` requires the full frame and `na.fail` refuses a gap, and
#' they are set to the median (or the first level) so nothing is out of range.
.neutral_row <- function(d, vars) {
  out <- lapply(vars, function(v) {
    x <- d[[v]]
    if (is.factor(x)) factor(levels(x)[1], levels = levels(x))
    else if (is.logical(x)) FALSE
    else stats::median(x, na.rm = TRUE)
  })
  names(out) <- vars
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' The partial effect of one smooth, on a shared grid, with its standard error.
#'
#' `unconditional = TRUE` IS NOT USED, BECAUSE IT DOES NOT WORK ON A `bam`
#' OBJECT. FOUND 2026-09-05 and confirmed against the installed mgcv 1.9.4
#' source: `mgcv:::predict.bam` nulls a list of large components to save memory
#' before delegating, and `Vc` -- the smoothness-uncertainty corrected
#' covariance -- is one of them:
#'
#'   object$Sl <- object$qrx <- object$R <- object$F <- object$Ve <-
#'     object$Vc <- object$G <- ... <- NULL
#'
#' `predict.gam` then finds `Vc` missing, warns once per call, and silently
#' falls back to the conditional `Vp`. MEASURED on `hemoglobin/meas` from the
#' live bundle: `unconditional = TRUE` returns standard errors BITWISE
#' IDENTICAL to `unconditional = FALSE`, so the argument is a no-op. The
#' correction is not unavailable -- every one of the 43 bundle GAMs carries a
#' non-null `Vc` that differs from its `Vp` -- it is thrown away in transit.
#'
#' Assigning `Vp <- Vc` before predicting is exactly what `predict.gam` would
#' have done and restores the correction. It widens the standard error by a
#' median of 4.7% and by up to 29.2% on that model.
#'
#' THIS ALSO BEARS ON THE ATTRIBUTABILITY ARM.
#' `docs/v2_attributability_plan_20260902.md` section 15 proposes the analytic
#' interval from `se.fit = TRUE, unconditional = TRUE` as a headline result that
#' SHAP has no counterpart for. Every layer-1 model is a `bam`, so as written
#' that interval would be the conditional one while being reported as the
#' corrected one.
.partial <- function(b, nd, label) {
  if (!is.null(b$Vc)) b$Vp <- b$Vc
  p <- stats::predict(b, newdata = nd, type = "terms", se.fit = TRUE,
                      discrete = FALSE)
  j <- match(label, colnames(p$fit))
  if (is.na(j)) return(NULL)
  list(fit = as.numeric(p$fit[, j]), se = as.numeric(p$se.fit[, j]),
       vc = !is.null(b$Vc))
}

#' Level and shape displacement between two partial-effect curves, in nats.
.displace <- function(a, b) {
  dif <- b$fit - a$fit
  lvl <- mean(dif)
  ctr <- dif - lvl
  list(level = lvl,
       shape = stats::sd(ctr),
       max_abs = max(abs(ctr)),
       # Scaled by the size of the effect being displaced, so a shape change is
       # readable as a fraction of the curve it is a change to.
       rel = if (stats::sd(a$fit) > 0) stats::sd(ctr) / stats::sd(a$fit) else NA_real_,
       cor = if (stats::sd(a$fit) > 0 && stats::sd(b$fit) > 0)
               stats::cor(a$fit, b$fit) else NA_real_,
       # The pooled standard error of the difference at each grid point, taken
       # at its median. Both curves come from the same rows, so treating them as
       # independent OVERSTATES this -- it is a conservative yardstick and is
       # labelled as one, not a test.
       se_med = stats::median(sqrt(a$se^2 + b$se^2)))
}

# ---------------------------------------------------------------------------
cat("\n=== measurement-smooth displacement, meas -> full ===\n\n")
cat("    level    mean of f_full(x) - f_meas(x), in nats. Absorbed by the\n")
cat("             intercept; carries no ordering information.\n")
cat("    shape    sd of that difference after centring, in nats. THE quantity.\n")
cat("    null     the same shape statistic when p(I|Y) is preserved exactly and\n")
cat("             p(M,I|Y) is broken, over B permutations. This is the floor set\n")
cat("             by discretisation and by penalty reallocation.\n")
cat("    excess   shape - null_mean. Above null_hi is the finding.\n")
cat("    se_med   median pooled standard error of the difference. A yardstick,\n")
cat("             deliberately conservative; not a test.\n\n")
cat(sprintf("%-16s %-19s %8s %8s %8s %8s %9s %8s %7s\n",
            "signal", "term", "level", "shape", "null_m", "null_hi",
            "excess", "se_med", "cor"))

D <- list()
t0 <- start_timer()

for (sg in signals) {
  pri <- priors_for(priors, sg, "final")
  # ONE frame for both fits. The `meas` and `full` frames are already bitwise
  # identical on their shared columns (verified 2026-09-05), so this changes no
  # number; it removes the possibility that a future divergence between the two
  # frame builders would be read here as a displacement.
  d  <- signal_frame(sg, "full", tabs, cfg, pri, stay_ids = tr)
  yv <- as.integer(d$mortality)

  f_m <- build_formula(sg, "meas", cfg)
  f_f <- build_formula(sg, "full", cfg)
  b_m <- try(.bam_fit(f_m, d, cfg), silent = TRUE)
  b_f <- try(.bam_fit(f_f, d, cfg), silent = TRUE)
  if (inherits(b_m, "try-error") || inherits(b_f, "try-error")) {
    cat(sprintf("%-16s  FIT FAILED\n", sg)); next
  }

  meas_v <- intersect(unique(vapply(b_m$smooth, function(s) s$term[1], character(1))),
                      .term_vars(measurement_terms(sg, cfg)))
  all_v  <- unique(c(all.vars(f_f)[-1]))
  all_v  <- intersect(all_v, names(d))
  iv_cols <- setdiff(all_v, c(meas_v, "mortality", "stay_id"))

  # The null replicates, fitted once per replicate and reused for every term.
  nulls <- vector("list", B)
  for (b in seq_len(B)) {
    dd <- .permute_block_within_class(d, iv_cols, yv, seed = SEED + b)
    bb <- try(.bam_fit(f_f, dd, cfg), silent = TRUE)
    nulls[[b]] <- if (inherits(bb, "try-error")) NULL else bb
  }
  nulls <- Filter(Negate(is.null), nulls)

  for (v in meas_v) {
    lab <- paste0("s(", v, ")")
    # The grid is the covariate's own quantiles, so the displacement is weighted
    # by where patients actually are rather than uniformly across a range whose
    # tails may hold almost nobody.
    grid <- unname(stats::quantile(d[[v]], probs = seq(0.01, 0.99, length.out = N_GRID),
                                   names = FALSE, type = 7))
    nd <- .neutral_row(d, all_v)[rep(1L, N_GRID), , drop = FALSE]
    nd[[v]] <- grid
    rownames(nd) <- NULL

    pm <- .partial(b_m, nd, lab); pf <- .partial(b_f, nd, lab)
    if (is.null(pm) || is.null(pf)) next
    obs <- .displace(pm, pf)

    ns <- vapply(nulls, function(bb) {
      p <- .partial(bb, nd, lab)
      if (is.null(p)) NA_real_ else .displace(pm, p)$shape
    }, numeric(1))
    ns <- ns[!is.na(ns)]
    null_m  <- if (length(ns)) mean(ns) else NA_real_
    null_hi <- if (length(ns)) unname(stats::quantile(ns, 0.95)) else NA_real_

    row <- data.frame(
      signal = sg, term = v, n_grid = N_GRID,
      level = round(obs$level, 5), shape = round(obs$shape, 5),
      max_abs_centred = round(obs$max_abs, 5),
      shape_rel = round(obs$rel, 4), curve_cor = round(obs$cor, 5),
      se_median = round(obs$se_med, 5), vc_corrected = isTRUE(pm$vc) && isTRUE(pf$vc),
      null_mean = round(null_m, 5), null_hi95 = round(null_hi, 5),
      n_null = length(ns),
      excess = round(obs$shape - null_m, 5),
      above_null = isTRUE(obs$shape > null_hi),
      stringsAsFactors = FALSE)
    D[[length(D) + 1L]] <- row
    cat(sprintf("%-16s %-19s %+8.4f %8.4f %8.4f %8.4f %+9.4f %8.4f %7.3f%s\n",
                sg, v, row$level, row$shape, row$null_mean, row$null_hi95,
                row$excess, row$se_median, row$curve_cor,
                if (isTRUE(row$above_null)) "  *" else ""))
  }
}

D <- do.call(rbind, D)
save_table(run, D, "coupling_displacement", subdir = "diagnostics")

cat("\n=== summary ===\n\n")
cat(sprintf("  %d measurement smooths across %d paired signals, B = %d, %.1f minutes.\n",
            nrow(D), length(unique(D$signal)), B, t0()$elapsed_sec / 60))
cat(sprintf("  smooths whose shape change clears its own null (marked *): %d of %d\n",
            sum(D$above_null, na.rm = TRUE), nrow(D)))
cat(sprintf("  largest excess: %+.4f nats   median excess: %+.4f nats\n",
            max(D$excess, na.rm = TRUE), stats::median(D$excess, na.rm = TRUE)))
cat(sprintf("  smooths whose shape change exceeds its own median standard error: %d\n",
            sum(D$shape > D$se_median, na.rm = TRUE)))
cat("\n  A LEVEL SHIFT IS NOT A FINDING. It is absorbed by the intercept and\n")
cat("  changes no patient's position in the ordering. Read the shape column,\n")
cat("  against the null beside it and not against zero.\n")

finalize_run(run, extra = list(B = B, n_terms = nrow(D), n_signals = length(signals),
                               n_grid = N_GRID, seed = SEED))
cat(sprintf("\nwritten: %s\n", run$path))
