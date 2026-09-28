# R/09_metrics.R -------------------------------------------------------------
# Discrimination and risk-ordering for an aggregated score. APPLY-TIME ONLY:
# nothing here fits a layer-1 model, and nothing here may be consulted while
# choosing one.
#
# The score is the row sum of an L matrix. That is the NAIVE aggregation — it
# assumes the 19 L's are conditionally independent given the outcome, which the
# branch point already says they are not (mean |r| = 0.21, PC1 = 26.5%). Naive
# summing is therefore the BASELINE the weighted layer-2 aggregation has to
# beat, not the final answer. llr_calibration() below measures exactly how
# wrong the independence assumption is, in one number.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). Every function that writes takes a
# `run` object and goes through R/11.
#
# AGGREGATES ONLY (hard rule 1). Bins, counts, rates and areas. The score
# vector is row-level and is never printed or returned upward.
# ----------------------------------------------------------------------------

# --- discrimination ---------------------------------------------------------

#' AUROC, computed from ranks rather than from a curve.
#'
#' The Mann-Whitney identity, so it is exact rather than trapezoid-approximated
#' and `rank()`'s average-ranks tie handling is exactly the right convention for
#' tied scores. No dependency, and it cannot disagree with pROC.
.auroc <- function(score, y) {
  y <- as.integer(y)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' AUPRC as average precision, with ties collapsed.
#'
#' Average precision — sum over thresholds of (recall gain) x precision — rather
#' than trapezoidal area, because trapezoidal interpolation between PR points is
#' not achievable by any classifier and reliably over-states the area. Ties are
#' collapsed to one threshold each, which is what stops a block of equal scores
#' being credited as if it were ordered.
#'
#' The floor is the event rate, not 0.5: at 10.3% mortality an AUPRC of 0.30 is
#' a 3x lift, and reading it against 0.5 the way one reads AUROC is the standard
#' way to misjudge this number.
.auprc <- function(score, y) {
  y <- as.integer(y)
  n1 <- sum(y == 1L)
  if (n1 == 0L) return(NA_real_)
  o  <- order(score, decreasing = TRUE)
  ys <- y[o]; ss <- score[o]
  tp <- cumsum(ys); fp <- cumsum(1L - ys)
  last <- c(ss[-1L] != ss[-length(ss)], TRUE)      # end of each tie block
  tp <- tp[last]; fp <- fp[last]
  prec <- tp / (tp + fp)
  rec  <- tp / n1
  sum(diff(c(0, rec)) * prec)
}

#' The bootstrap resampling unit: patients when a grouping is given, stays
#' otherwise.
#'
#' STATISTICAL REVIEW S4 (2026-09-09). The split and folds group repeat stays
#' by patient, so the fitting protocol treats stays of one patient as
#' dependent; the reported intervals resampled individual stays and treated
#' them as independent. A cluster bootstrap draws PATIENTS with replacement
#' and keeps every stay of a drawn patient, with multiplicity, so the interval
#' respects the dependence the design declares. THE ESTIMAND IS UNCHANGED: the
#' metric is still computed over stays, one row each; only the resampling unit
#' moves. `unit` is carried into every table that uses a draw so a stay-level
#' interval can never be read as a patient-level one.
#'
#' A table produced with `group = NULL` is a conditional, stay-resampled
#' summary and says so in its `boot_unit` column. The graph and the runners
#' always pass patients.
.boot_units <- function(n, group = NULL) {
  if (is.null(group)) {
    return(list(unit = "stay", n_units = n,
                draw = function() sample.int(n, n, replace = TRUE)))
  }
  if (length(group) != n) stop(".boot_units: `group` is not aligned", call. = FALSE)
  if (anyNA(group)) stop(".boot_units: `group` has NA", call. = FALSE)
  g   <- as.character(group)
  ug  <- unique(g)
  idx <- split(seq_len(n), factor(g, levels = ug))
  m   <- length(ug)
  list(unit = "patient", n_units = m,
       draw = function() unlist(idx[sample.int(m, m, replace = TRUE)],
                                use.names = FALSE))
}

#' Discrimination for one score, with percentile bootstrap intervals.
#'
#' Bootstrap rather than DeLong so that AUROC and AUPRC get intervals from the
#' same procedure — DeLong has no AUPRC analogue, and two intervals built by
#' different machinery invite exactly the wrong comparison.
#'
#' WHAT THE INTERVAL IS. It resamples FIXED out-of-fold (or applied) scores,
#' so it is the sampling uncertainty of the metric given these fitted models.
#' It does not include refitting variability or the dependence the overlapping
#' training folds induce; that is a separate, refit-based experiment. Read it
#' as a conditional score-based interval, which is how the tables label it.
#'
#' @param score row-level; never printed
#' @param y     0/1 outcome
#' @param n_boot 0 to skip the intervals
#' @param group patient id per row; the cluster-bootstrap unit (review S4).
#'   NULL resamples stays and the row says so.
#' @return one-row data frame
score_metrics <- function(score, y, label = NA_character_, n_boot = 200L,
                          seed = 1L, group = NULL) {
  stopifnot(length(score) == length(y))
  if (!is.null(group) && length(group) != length(y)) {
    stop("score_metrics: `group` is not aligned to `score`", call. = FALSE)
  }
  ok <- !is.na(score) & !is.na(y)
  score <- score[ok]; y <- as.integer(y[ok])
  if (!is.null(group)) group <- group[ok]
  bu <- .boot_units(length(y), group)

  auroc <- .auroc(score, y)
  auprc <- .auprc(score, y)
  out <- data.frame(
    label      = label,
    n          = length(y),
    n_events   = sum(y),
    event_rate = round(mean(y), 5),
    auroc      = round(auroc, 5),
    auprc      = round(auprc, 5),
    auprc_lift = round(auprc / mean(y), 3),   # against the only honest floor
    stringsAsFactors = FALSE)

  if (n_boot > 0L) {
    b <- with_seed(seed, {
      t(vapply(seq_len(n_boot), function(i) {
        j <- bu$draw()
        c(.auroc(score[j], y[j]), .auprc(score[j], y[j]))
      }, numeric(2)))
    })
    q <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
    out$auroc_lo <- round(q(b[, 1])[1], 5); out$auroc_hi <- round(q(b[, 1])[2], 5)
    out$auprc_lo <- round(q(b[, 2])[1], 5); out$auprc_hi <- round(q(b[, 2])[2], 5)
    out$n_boot   <- n_boot
  } else {
    # A STABLE COLUMN SET (plumbing review F6, 2026-09-08). `n_boot = 0` is
    # documented above as "skip the intervals", and until now it also skipped
    # the COLUMNS, so every consumer that selects `auroc_lo` -- run/internal.R
    # prints it, the transport tables join on it -- failed on a table with no
    # such column. The interval is absent; the column that says so is not.
    out$auroc_lo <- NA_real_; out$auroc_hi <- NA_real_
    out$auprc_lo <- NA_real_; out$auprc_hi <- NA_real_
    out$n_boot   <- 0L
  }
  # The resampling unit and how many of them there were, so that a stay-level
  # interval is never mistaken for a patient-level one (review S4).
  out$boot_unit <- bu$unit
  out$n_units   <- as.integer(bu$n_units)
  out
}

# --- comparing two scores on the SAME rows ----------------------------------

#' DeLong's test for two correlated AUROCs.
#'
#' WHY THIS EXISTS. Every comparison in this project scores two methods on
#' IDENTICAL rows, so the two AUROCs are strongly positively correlated and
#' their marginal bootstrap intervals overlap long after the difference has
#' become reliable. Reading overlap as "no difference" is the standard way to
#' get this wrong, and docs/v2_state_20260828.md flags it twice: once for
#' `raw -> feat` in the XGBoost ladder and once for the severity-score arms in
#' section 5.5. A paired test is the only thing that settles either.
#'
#' This is the fast DeLong algorithm (Sun and Xu 2014): the structural
#' components V10 and V01 are obtained from midranks rather than from the
#' O(mn) pairwise comparison, so it is linear in the sample size after sorting
#' and runs on 41,000 rows in milliseconds.
#'
#' The test is on AUROC only. There is no DeLong analogue for AUPRC — for that
#' difference use the paired bootstrap below, which resamples rows once and
#' scores both methods on the same resample.
#'
#' OBSERVATION-LEVEL, AND LABELLED AS SUCH (review S4). The placement
#' covariance treats every stay as an independent observation; it does not
#' know about repeat stays of one patient. `arm_contrasts()` therefore reports
#' its AUROC interval from the patient-clustered paired bootstrap and carries
#' the DeLong z and p beside it under their own name, as the iid reference,
#' never as a cluster-adjusted test.
#'
#' @param s1,s2 two scores on the same rows, any monotone scale
#' @param y     0/1 outcome
#' @return one-row data frame: both AUROCs, the difference, its standard error,
#'   z and a two-sided p value
delong_test <- function(s1, s2, y) {
  ok <- !is.na(s1) & !is.na(s2) & !is.na(y)
  s1 <- s1[ok]; s2 <- s2[ok]; y <- as.integer(y[ok])
  m <- sum(y == 1L); n <- sum(y == 0L)
  if (m == 0L || n == 0L) stop("delong_test: one class is empty", call. = FALSE)

  # Structural components. V10 is per positive case, V01 per negative case;
  # mean(V10) == mean(V01) == AUROC, which is asserted below.
  comp <- function(s) {
    x <- s[y == 1L]; z <- s[y == 0L]
    r  <- rank(c(x, z))                       # midranks over the pooled sample
    v10 <- (r[seq_len(m)] - rank(x)) / n
    v01 <- 1 - (r[m + seq_len(n)] - rank(z)) / m
    list(v10 = v10, v01 = v01, auc = mean(v10))
  }
  a <- comp(s1); b <- comp(s2)

  V10 <- cbind(a$v10, b$v10)
  V01 <- cbind(a$v01, b$v01)
  S   <- stats::cov(V10) / m + stats::cov(V01) / n
  L   <- c(1, -1)
  var_d <- drop(L %*% S %*% L)
  d  <- a$auc - b$auc
  se <- sqrt(max(var_d, 0))
  # A zero standard error means the two scores induce the same ranking on every
  # case — comparing a score with itself, or with a monotone transformation of
  # itself. The difference is then exactly zero and p is 1, which is a real
  # answer and not a degenerate one. Only a zero SE with a non-zero difference
  # is undefined, and that cannot happen.
  z <- if (se > 0) d / se else if (d == 0) 0 else NA_real_

  data.frame(
    auroc_1 = round(a$auc, 5), auroc_2 = round(b$auc, 5),
    delta   = round(d, 5),
    se      = round(se, 5),
    ci_lo   = round(d - 1.96 * se, 5),
    ci_hi   = round(d + 1.96 * se, 5),
    z       = round(z, 4),
    p_value = signif(2 * stats::pnorm(-abs(z)), 4),
    n = length(y), n_events = m,
    stringsAsFactors = FALSE)
}

#' Paired bootstrap for the difference in ANY two-score metric.
#'
#' The AUPRC companion to `delong_test()`, and the fallback wherever an
#' analytic paired test does not exist. One row resample per replicate, both
#' scores evaluated on it, so the correlation between the two methods is
#' preserved instead of being thrown away — which is the whole point.
#'
#' THE P-VALUE IS A SHIFTED-NULL BOOTSTRAP TEST (statistical review S5,
#' 2026-09-09). It was `2 * min(mean(b <= 0), mean(b >= 0))`: zero differences
#' counted in BOTH tails, so two identical score vectors returned p = 2, and
#' undefined replicates were dropped from the quantiles but not from the means.
#' It was also not a test of anything: it counted an ordinary bootstrap
#' distribution around the observed estimate against zero. The replacement is
#' deliberate on each point the review named:
#'
#'   null        the bootstrap differences are recentred on the observed
#'               difference, `b - d0`, which is the bootstrap distribution
#'               under a null of no difference (the shift construction).
#'   statistic   |d0|, two-sided: p is the null-distribution mass at least as
#'               far from zero as the observed difference.
#'   ties        a replicate exactly as extreme counts as extreme (`>=`), so
#'               identical vectors give d0 = 0, every replicate at 0, p = 1.
#'   invalid     non-finite replicates are dropped from BOTH the interval and
#'               the test, and their count is reported.
#'   finite B    `(1 + count) / (B_ok + 1)`, so p is never 0 from B draws.
#'
#' It remains a bootstrap approximation, resampling by the same unit as the
#' interval, and is reported beside the interval rather than instead of it.
#'
#' @param metric a function(score, y) returning one number; `.auroc` or `.auprc`
#' @param group  patient id per row; the cluster-bootstrap unit (review S4)
paired_boot_diff <- function(s1, s2, y, metric = .auprc, n_boot = 200L,
                             seed = 1L, label = NA_character_, group = NULL) {
  out <- .paired_boot(s1, s2, y, metrics = list(metric = metric), n_boot = n_boot,
                      seed = seed, group = group)[[1L]]
  cbind(label = label, out, stringsAsFactors = FALSE)
}

#' The paired bootstrap, drawn ONCE for any number of metrics.
#'
#' `arm_contrasts()` needs the AUROC and the AUPRC difference on the SAME
#' resamples, so that the two intervals of one contrast are paired with each
#' other as well as across arms. Calling `paired_boot_diff()` twice under one
#' seed would achieve that only by the coincidence of identical draw sequences
#' and would draw every replicate twice; this draws each replicate once and
#' evaluates every metric on it.
#'
#' @param metrics named list of functions(score, y) -> one number
#' @return named list of one-row data frames, one per metric, each with delta,
#'   ci_lo, ci_hi, p_boot, n_boot, n_boot_ok, boot_unit
.paired_boot <- function(s1, s2, y, metrics, n_boot = 200L, seed = 1L, group = NULL) {
  stopifnot(is.list(metrics), length(metrics) > 0L, !is.null(names(metrics)))
  ok <- !is.na(s1) & !is.na(s2) & !is.na(y)
  s1 <- s1[ok]; s2 <- s2[ok]; y <- as.integer(y[ok])
  if (!is.null(group)) group <- group[ok]
  bu <- .boot_units(length(y), group)
  d0 <- vapply(metrics, function(m) m(s1, y) - m(s2, y), numeric(1))
  B <- with_seed(seed, t(vapply(seq_len(n_boot), function(i) {
    j <- bu$draw()
    vapply(metrics, function(m) m(s1[j], y[j]) - m(s2[j], y[j]), numeric(1))
  }, numeric(length(metrics)))))
  if (n_boot == 0L) B <- matrix(numeric(0), 0L, length(metrics))
  lapply(seq_along(metrics), function(k) {
    bb <- B[, k]; bb <- bb[is.finite(bb)]
    q  <- if (length(bb)) stats::quantile(bb, c(0.025, 0.975)) else c(NA_real_, NA_real_)
    p  <- if (length(bb) && is.finite(d0[k]))
      (1 + sum(abs(bb - d0[k]) >= abs(d0[k]))) / (length(bb) + 1) else NA_real_
    data.frame(delta = round(d0[k], 5),
               ci_lo = round(q[1], 5), ci_hi = round(q[2], 5),
               p_boot = signif(p, 4),
               n_boot = n_boot, n_boot_ok = length(bb), boot_unit = bu$unit,
               stringsAsFactors = FALSE, row.names = NULL)
  }) |> stats::setNames(names(metrics))
}

# --- risk ordering ----------------------------------------------------------

#' Wilson score interval. Not Wald.
#'
#' The bottom bins hold few deaths — sometimes single digits — and a Wald
#' interval there is both too narrow and capable of covering negative
#' probabilities. Wilson stays inside [0, 1] and is well behaved at small k,
#' which is exactly the regime that decides whether an apparent inversion in the
#' curve is real.
.wilson <- function(k, n, conf = 0.95) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, ctr - hw), min(1, ctr + hw))
}

#' Equal-count bins of the score, with observed mortality and its interval.
#'
#' Equal-COUNT (quantile) rather than equal-WIDTH: the score is roughly bell
#' shaped, so equal-width bins would put almost every stay in the middle three
#' and leave the tails with a handful each — and the tails are where the
#' clinically interesting behaviour is.
#'
#' THE BIN INTERVALS ARE STAY-LEVEL WILSON INTERVALS and the table says so in
#' `ci_method` (external runner review E5, 2026-09-09). The discrimination
#' intervals and paired contrasts resample PATIENTS (review S4); the per-bin
#' mortality interval still treats the stays in a bin as independent draws.
#' Repeat stays of one patient are a small fraction of both cohorts, so the
#' practical difference is small, but the two kinds of interval rest on
#' different assumptions and must not be described as one.
#'
#' @return data frame, one row per bin
risk_bins <- function(score, y, n_bins = 20L, breaks = NULL) {
  ok <- !is.na(score) & !is.na(y)
  score <- score[ok]; y <- as.integer(y[ok])

  # FROZEN CUT POINTS, when supplied. Quantiles of the score being scored are
  # right for a within-site reading and wrong for a cross-site one: eICU
  # quantiles would re-centre every bin on eICU's own distribution, so bin 20
  # would mean "the sickest 5% of eICU" rather than "a MIMIC-calibrated score
  # above x". Passing the training cut points asks the transport question
  # instead: does the same score value carry the same risk at the other site?
  # The ends are opened to +/-Inf so a score outside MIMIC's observed range
  # lands in the extreme bin rather than becoming NA and vanishing.
  if (!is.null(breaks)) {
    br <- unique(sort(as.numeric(breaks)))
    br[1L] <- -Inf; br[length(br)] <- Inf
  } else {
    br <- stats::quantile(score, probs = seq(0, 1, length.out = n_bins + 1L),
                          na.rm = TRUE, names = FALSE)
    br <- unique(br)
  }
  if (length(br) < 3L) stop("risk_bins: score has too few distinct values", call. = FALSE)
  bin <- cut(score, breaks = br, include.lowest = TRUE, labels = FALSE)

  rows <- lapply(sort(unique(bin)), function(b) {
    s <- bin == b
    k <- sum(y[s]); n <- sum(s)
    ci <- .wilson(k, n)
    data.frame(bin = b, n = n, deaths = k,
               obs_rate = k / n, lo = ci[1], hi = ci[2],
               ci_method = "wilson_stay_iid",
               mean_score = mean(score[s]),
               min_score  = min(score[s]), max_score = max(score[s]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Is the curve actually monotone, or does it just look it?
#'
#' Three numbers, because "monotone" is not one question:
#'
#'   spearman        rank correlation of bin index with observed rate. The
#'                   headline, and insensitive to the shape.
#'   n_inversions    adjacent bins where the rate goes DOWN as the score goes up.
#'   n_sig_inversions  those where the two Wilson intervals do not even overlap.
#'                   This is the one that matters: a handful of overlapping
#'                   inversions in 20 bins is sampling noise, and reporting the
#'                   raw count without this alongside overstates the problem.
#'   rate_ratio      top bin rate / bottom bin rate. The number a clinician
#'                   reads first.
monotonicity <- function(bins) {
  d <- bins[order(bins$bin), ]
  dr <- diff(d$obs_rate)
  inv <- which(dr < 0)
  sig <- inv[vapply(inv, function(i) d$lo[i] > d$hi[i + 1L], logical(1))]
  data.frame(
    n_bins           = nrow(d),
    spearman         = round(stats::cor(d$bin, d$obs_rate, method = "spearman"), 4),
    n_inversions     = length(inv),
    n_sig_inversions = length(sig),
    worst_inversion  = if (length(inv)) round(min(dr), 5) else 0,
    rate_bottom      = round(d$obs_rate[1], 5),
    rate_top         = round(d$obs_rate[nrow(d)], 5),
    rate_ratio       = round(d$obs_rate[nrow(d)] / max(d$obs_rate[1], .Machine$double.eps), 2),
    stringsAsFactors = FALSE)
}

# --- is naive summing valid? ------------------------------------------------

#' The one number that says whether the L's can simply be added.
#'
#' If each L were a correctly calibrated log-likelihood ratio and the 19 were
#' conditionally independent given the outcome, then Bayes gives EXACTLY
#'
#'     logit p(Y=1 | all L) = logit(p_bar) + sum_g L_g
#'
#' i.e. regressing the outcome on the summed score must return slope 1 and
#' intercept logit(p_bar). So fit that regression and read the slope:
#'
#'   slope ~ 1   independence holds well enough; naive summing is calibrated
#'   slope < 1   the L's are POSITIVELY correlated, so the sum double-counts
#'               shared evidence and is over-confident. The shortfall is the
#'               quantitative case for Sigma-inverse: 1/slope is roughly the
#'               factor by which the naive sum overstates its own evidence.
#'   slope > 1   under-confident, which would be surprising here
#'
#' This is a far sharper test than a calibration plot, and it is the direct
#' answer to "do I need to whiten". Given mean |r| = 0.21 across the L matrix,
#' expect a slope materially below 1.
#'
#' Fitted on the individual rows, not on the bins: binning first would throw
#' away information and make the slope depend on the bin count.
#'
#' `slope_lo`/`slope_hi` IS A PROFILE-LIKELIHOOD INTERVAL FROM AN ORDINARY GLM
#' and treats stays as independent; `ci_method` records that (external runner
#' review E5, 2026-09-09). The slope itself does not depend on `p_bar`;
#' `p_bar` only sets `expected_intercept`, and the table carries the value
#' used so a site-prevalence reference cannot be mistaken for a
#' training-prevalence one.
llr_calibration <- function(score, y, p_bar) {
  ok <- !is.na(score) & !is.na(y)
  fit <- stats::glm(as.integer(y[ok]) ~ score[ok], family = stats::binomial())
  cf <- stats::coef(fit)
  ci <- suppressMessages(try(stats::confint(fit, level = 0.95), silent = TRUE))
  slope_lo <- slope_hi <- NA_real_
  if (!inherits(ci, "try-error")) { slope_lo <- ci[2, 1]; slope_hi <- ci[2, 2] }
  data.frame(
    slope              = round(unname(cf[2]), 4),
    slope_lo           = round(slope_lo, 4),
    slope_hi           = round(slope_hi, 4),
    intercept          = round(unname(cf[1]), 4),
    expected_intercept = round(logit(p_bar), 4),
    expected_slope     = 1,
    overconfidence     = round(1 / unname(cf[2]), 3),
    p_bar_reference    = round(p_bar, 6),
    ci_method          = "glm_profile_stay_iid",
    stringsAsFactors = FALSE)
}

# --- per-signal contribution ------------------------------------------------

#' AUROC and AUPRC of each L column on its own.
#'
#' Answers "is one signal doing all the work", which is the discrimination-side
#' companion to the eigenspectrum. If the summed score barely beats the best
#' single column, aggregation is not earning its place.
#'
#' AUPRC WAS ADDED 2026-09-05 and the omission mattered. Every arm-level table
#' in the project reports AUROC and AUPRC together, for the reason
#' `score_metrics()` states: at a roughly 12% event rate AUROC is dominated by
#' the ordering of the survivors, and a column can order the cohort acceptably
#' while adding nothing where the deaths are. This table alone reported AUROC by
#' itself, so the one place the project asks "which signals carry the evidence"
#' was answered with the one metric least able to say. `auprc_lift` is AUPRC
#' divided by the event rate, which is the only honest floor: 1.0 is a column
#' that has learned nothing.
#'
#' `n_nonzero` is carried because a zero-filled column's discrimination is
#' partly a statement about how many stays have the signal at all -- an
#' unmeasured stay sits at exactly 0 by assignment (spec 5.5), which places it
#' in the middle of the ordering rather than removing it. IT IS A PROXY AND IS
#' NAMED AS ONE: the authoritative measured mask is `measured_matrix()`, which
#' reads `n_obs` and is not available from an L matrix alone. The two differ
#' only where a measured stay's L is exactly 0.0, which is possible and
#' vanishingly rare. Do not report `n_nonzero` as a coverage count.
#'
#' ORDERED BY AUROC DESCENDING, as it always was. Ordering by `auprc_lift`
#' instead was tried on 2026-09-05 and reverted: it changes the row order of a
#' table that already exists in nine internal run directories, for a cosmetic
#' gain, and it turned out that `metrics_report()` was reading `ps$auroc[1]` as
#' "the best signal" -- an implicit dependency on this sort order that nothing
#' stated. The dependency is now gone (see `metrics_report()`, which selects
#' with `which.max()`), so the order here is a presentation choice again rather
#' than a load-bearing one.
signal_auroc <- function(M, y) {
  y <- as.integer(y)
  ev <- mean(y)
  auroc <- vapply(seq_len(ncol(M)), function(j) .auroc(M[, j], y), numeric(1))
  auprc <- vapply(seq_len(ncol(M)), function(j) .auprc(M[, j], y), numeric(1))
  out <- data.frame(
    signal     = colnames(M),
    n          = nrow(M),
    n_events   = sum(y),
    event_rate = round(ev, 5),
    auroc      = round(auroc, 5),
    auprc      = round(auprc, 5),
    auprc_lift = round(auprc / ev, 3),
    n_nonzero  = as.integer(colSums(M != 0)),
    stringsAsFactors = FALSE)
  out[order(-out$auroc), , drop = FALSE]
}

# --- figures ----------------------------------------------------------------

#' The risk-ordering curve, on both scales.
#'
#' TWO panels, because they answer different questions and neither alone is
#' enough:
#'
#'   probability scale  the clinical read. Does risk rise monotonically across
#'                      the score? Wilson intervals included, because without
#'                      them a noisy bottom bin reads as a broken model.
#'   logit scale        the statistical read, with the slope-1 line the theory
#'                      predicts. Distance from that line IS the cost of naive
#'                      summing, made visible rather than asserted.
#'
#' Written through R/11's save_fig(), so this file constructs no path.
plot_risk_curve <- function(run, bins, cal, p_bar, name = "risk_curve",
                            title = "Out-of-fold risk ordering") {
  d <- bins[order(bins$bin), ]
  use_gg <- requireNamespace("ggplot2", quietly = TRUE)

  p <- save_fig(run, name, width = 11, height = 5, dpi = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (use_gg) {
    d$panel <- "observed mortality"
    e <- d
    e$panel <- "logit scale, vs the slope-1 prediction"
    lg <- function(x) log(pmax(x, 1e-4) / (1 - pmin(x, 1 - 1e-4)))
    e$obs_rate <- lg(e$obs_rate); e$lo <- lg(e$lo); e$hi <- lg(e$hi)
    dd <- rbind(d, e)
    ref <- data.frame(panel = "logit scale, vs the slope-1 prediction",
                      x = range(d$mean_score))
    ref$y <- logit(p_bar) + ref$x
    g <- ggplot2::ggplot(dd, ggplot2::aes(x = mean_score, y = obs_rate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.18) +
      ggplot2::geom_line(colour = "grey30") +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::geom_line(data = ref, ggplot2::aes(x = x, y = y),
                         linetype = "dashed", colour = "firebrick") +
      ggplot2::facet_wrap(~panel, scales = "free_y") +
      ggplot2::labs(
        title = title,
        subtitle = sprintf(
          "%d bins, equal count | calibration slope %.3f (1.000 = naive sum valid) | Spearman %.3f",
          nrow(d), cal$slope, stats::cor(d$bin, d$obs_rate, method = "spearman")),
        x = "mean summed LLR in bin", y = NULL,
        caption = "dashed: logit(p_bar) + score, the prediction under conditional independence") +
      ggplot2::theme_minimal(base_size = 11)
    print(g)
  } else {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    graphics::plot(d$mean_score, d$obs_rate, type = "b", pch = 19,
                   ylim = c(0, max(d$hi, na.rm = TRUE)),
                   xlab = "mean summed LLR in bin", ylab = "observed mortality",
                   main = title)
    graphics::arrows(d$mean_score, d$lo, d$mean_score, d$hi, code = 3,
                     angle = 90, length = 0.02, col = "grey50")
    graphics::abline(h = p_bar, lty = 3, col = "grey40")
    lg <- function(x) log(pmax(x, 1e-4) / (1 - pmin(x, 1 - 1e-4)))
    graphics::plot(d$mean_score, lg(d$obs_rate), type = "b", pch = 19,
                   xlab = "mean summed LLR in bin", ylab = "logit(observed mortality)",
                   main = sprintf("calibration slope %.3f", cal$slope))
    graphics::abline(a = logit(p_bar), b = 1, lty = 2, col = "firebrick")
  }
  invisible(p)
}

# --- the whole thing --------------------------------------------------------

#' Metrics for ANY score on the log-odds scale: discrimination, bins,
#' monotonicity, calibration, plot.
#'
#' Split out of metrics_report() so the comparison arms can be measured on
#' EXACTLY the same code path as the proposed method. A baseline scored by a
#' second implementation of AUROC is not a baseline, it is a second experiment,
#' and the 2x2 in docs/v2_analysis_tiering.md turns on the four cells being
#' commensurable.
#'
#' `score` must be on the LOG-ODDS scale for `cal` to mean anything: the
#' slope-1 reference is the Bayes statement `logit p = logit(p_bar) + score`.
#' For a probability model pass `logit(p_hat) - logit(p_bar)`, which is the same
#' quantity a summed LLR is.
#'
#' @param score  row-level, on the log-odds scale; never printed
#' @param y      0/1 outcome, aligned to `score`
#' @param p_bar  the prior the score is centred on
#' @param title  plot title; defaults to the label
#' @return list(summary, bins, mono, cal)
score_report <- function(run, score, y, p_bar, label,
                         n_bins = 20L, n_boot = 200L, seed = 1L,
                         title = NULL, breaks = NULL, group = NULL) {
  sm   <- score_metrics(score, y, label = label, n_boot = n_boot, seed = seed,
                        group = group)
  bins <- risk_bins(score, y, n_bins = n_bins, breaks = breaks)
  mono <- monotonicity(bins)
  cal  <- llr_calibration(score, y, p_bar)

  sm$spearman   <- mono$spearman
  sm$rate_ratio <- mono$rate_ratio
  sm$cal_slope  <- cal$slope

  nm <- function(x) paste0(x, "_", label)
  save_table(run, bins, nm("risk_bins"))
  save_table(run, mono, nm("monotonicity"))
  save_table(run, cal,  nm("llr_calibration"))
  plot_risk_curve(run, bins, cal, p_bar, name = nm("risk_curve"),
                  title = title %||% sprintf("Out-of-fold risk ordering — %s", label))

  list(summary = sm, bins = bins, mono = mono, cal = cal)
}

#' Metrics for one L matrix. score_report() plus the two L-specific columns.
#'
#' @param M      an L matrix from R/07 (stays x signals)
#' @param y      0/1 outcome, aligned to rownames(M)
#' @param p_bar  cohort prior, for the slope-1 reference
#' @return list(summary, bins, mono, cal, per_signal)
metrics_report <- function(run, M, y, p_bar, label = "full",
                           n_bins = 20L, n_boot = 200L, seed = 1L) {
  r  <- score_report(run, rowSums(M), y, p_bar, label = label,
                     n_bins = n_bins, n_boot = n_boot, seed = seed,
                     title = sprintf("Out-of-fold risk ordering — L_%s", label))
  ps <- signal_auroc(M, y)

  # SELECTED EXPLICITLY, NOT BY ROW POSITION. This read `ps$auroc[1]` until
  # 2026-09-05, which was correct only because `signal_auroc()` happened to sort
  # by AUROC. Re-ordering that table by AUPRC lift -- a cosmetic change made
  # while adding the AUPRC columns -- would have turned `best_signal_auroc` into
  # "the AUROC of whichever signal had the best AUPRC lift" and left
  # `aggregation_gain` overstating the benefit of aggregation, with nothing
  # erroring and no column name changing. `which.max()` removes the coupling.
  jr <- which.max(ps$auroc)
  jp <- which.max(ps$auprc_lift)
  r$summary$best_signal            <- ps$signal[jr]
  r$summary$best_signal_auroc      <- ps$auroc[jr]
  r$summary$aggregation_gain       <- round(r$summary$auroc - ps$auroc[jr], 5)
  r$summary$best_signal_auprc      <- ps$signal[jp]
  r$summary$best_signal_auprc_lift <- ps$auprc_lift[jp]

  save_table(run, ps, paste0("signal_auroc_", label))
  c(r, list(per_signal = ps))
}

# --- frozen reporting bins ---------------------------------------------------

#' The cut points a score's reporting bins are defined by, frozen from training.
#'
#' Computed once on the training score and carried in the bundle, so MIMIC-test
#' and eICU report over the SAME intervals. Without this, "bin 20" means a
#' different thing at each site and the 20-bin curves are not comparable — each
#' site's curve would be monotone-by-construction in its own quantiles while
#' saying nothing about whether the score transports.
#'
#' Returned in full (n_bins + 1 values) rather than as interior cuts, so
#' `risk_bins(breaks = )` needs no reconstruction step that could differ.
#'
#' AGGREGATES ONLY (hard rule 1): these are quantiles of a score, not rows.
#'
#' TWO MODES, BECAUSE TWO KINDS OF SCORE (statistical review S11, 2026-09-09).
#' A summed LLR is continuous and `n_bins` equal-count bins are identifiable;
#' a tie at a quantile boundary there is a defect and the strict mode stops.
#' A recalibrated SEVERITY score is an affine map of an integer, so equal
#' scores are the rule, and demanding twenty-one distinct quantiles of it
#' would let valid, well-fitted severity data abort bundle construction after
#' every expensive fit had finished. `tie_aware = TRUE` is the declared
#' severity policy: the duplicated boundaries are collapsed, the vector is
#' stamped with the bin count it actually defines, and the SAME frozen vector
#' is applied at every site, so the bins are fewer but identical everywhere.
#' Equal scores are never split to manufacture a bin.
#'
#' @param tie_aware collapse tied boundaries and record the actual bin count,
#'   instead of stopping. For discrete scores only.
#' @return numeric vector of length n_bins + 1 (strict), or shorter under
#'   `tie_aware` with attributes `n_bins_requested`, `n_bins_actual`,
#'   `tie_aware`
score_cutpoints <- function(score, n_bins = 20L, tie_aware = FALSE) {
  s <- score[!is.na(score)]
  if (length(s) < n_bins + 1L) {
    stop("score_cutpoints: fewer scores than bins", call. = FALSE)
  }
  # Tie-aware breaks are OBSERVED score values (quantile type 1, the inverse
  # ECDF), so after the collapse every bin (b_i, b_{i+1}] holds at least the
  # score b_{i+1} on the training data and `n_bins_actual` counts occupied
  # bins. Interpolated (type 7) breaks between integers can enclose no value.
  br <- unname(stats::quantile(s, probs = seq(0, 1, length.out = n_bins + 1L),
                               na.rm = TRUE, type = if (isTRUE(tie_aware)) 1L else 7L))
  if (any(duplicated(br))) {
    if (!isTRUE(tie_aware)) {
      stop("score_cutpoints: the score has ties at a bin boundary, so ", n_bins,
           " bins are not identifiable. Reduce n_bins rather than deduplicating ",
           "here — a silently shortened break vector would give the two sites ",
           "different bin counts. For a discrete severity score pass ",
           "tie_aware = TRUE, which freezes the collapsed vector and its actual ",
           "bin count.", call. = FALSE)
    }
    br <- unique(br)
    if (length(br) < 3L) {
      stop("score_cutpoints: fewer than two distinct bins survive the ties",
           call. = FALSE)
    }
  }
  attr(br, "n_bins_requested") <- as.integer(n_bins)
  attr(br, "n_bins_actual")    <- length(br) - 1L
  attr(br, "tie_aware")        <- isTRUE(tie_aware)
  br
}

# --- transport ratios --------------------------------------------------------

#' Below this excess over chance, a training AUROC cannot anchor a retention
#' ratio. Declared here, once, before any table is built (external runner
#' review E8, 2026-09-09).
RETENTION_MIN_EXCESS <- 0.05

#' Fraction of the training site's ABOVE-CHANCE discrimination that survives.
#'
#' `auroc_site / auroc_train` is a defined number but a misleading one: a model
#' at chance (0.5) against a training AUROC of 0.8 reports 62.5 percent
#' "retained", because the quotient credits the 0.5 that any score gets for
#' free. Measured from chance instead, the same case is 0 / 0.3 = 0. The raw
#' quotient is still reported, under the name `auroc_ratio`, and the absolute
#' AUROCs and their difference remain the primary transport columns; this
#' ratio is the secondary reading.
#'
#' POLICY WHEN THE TRAINING AUROC IS NEAR CHANCE: the denominator
#' `auroc_train - 0.5` below `RETENTION_MIN_EXCESS` makes the ratio
#' arbitrarily large or negative for trivial differences, so the value is NA
#' rather than a number that would be quoted. A non-finite input is NA too.
#'
#' @return numeric vector, rounded to 4 places, NA where undefined
retention_above_chance <- function(auroc_site, auroc_train,
                                   min_excess = RETENTION_MIN_EXCESS) {
  ex  <- auroc_train - 0.5
  out <- (auroc_site - 0.5) / ex
  bad <- !is.finite(ex) | !is.finite(auroc_site) | ex < min_excess
  out[bad] <- NA_real_
  round(out, 4)
}

# --- stratified discrimination ----------------------------------------------

#' AUROC AND AUPRC within each level of a grouping variable.
#'
#' Built for eICU's 166 hospitals, where the pooled value answers a different
#' question from the distribution across sites: a score can discriminate well
#' overall partly because hospitals differ in case mix, while discriminating
#' less well WITHIN any one of them. `pooled_minus_median` puts the two side by
#' side. IT IS A DESCRIPTIVE CONTRAST, NOT A DECOMPOSITION (external runner
#' review E8, 2026-09-09): the pooled AUROC weights positive-negative pairs,
#' the median weights eligible hospitals equally, and hospital selection and
#' unequal within-hospital discrimination both move the gap. Three hospitals
#' with identical prevalence and identical marginal score distributions but
#' unequal sizes and discrimination give a pooled AUROC of 0.667 against a
#' median of 0, and none of that gap is case mix. A true decomposition would
#' split positive-negative pairs into within- and between-hospital pairs on
#' one cohort; this table does not do that and must not be read as if it did.
#'
#' THREE DENOMINATORS, REPORTED SEPARATELY (review E4). `n_scored` is every
#' row offered; `n_matched` drops rows with no group (at eICU, a stay the
#' hospital table does not cover); `n_stays_kept` is the rows inside groups
#' that meet the floors. `frac_kept_of_scored` is the coverage of the scored
#' cohort, `frac_kept_of_matched` the eligibility among matched rows. The
#' pooled metrics are computed twice: `*_matched` on every matched row and
#' `*_eligible` on the rows of reported groups only, and the contrast against
#' the median uses the ELIGIBLE pooled value so both sides describe one cohort.
#' Neither is the headline pooled AUROC of the scored cohort; that lives in
#' the score summary. Empty cases return a row with NA and a `status` rather
#' than an error or an infinite extremum.
#'
#' THE SUMMARY CARRIES BOTH DISTRIBUTIONS AS OF 2026-09-07. The per-group rows
#' had always reported AUPRC; the summary quantiles reported AUROC alone, so
#' the one table in the project built to describe a SPREAD described half of it.
#'
#' AUPRC IS SUMMARISED AS A LIFT, NOT RAW, and that is not a preference. AUPRC's
#' floor is the event rate, and the event rate is exactly what differs between
#' hospitals -- eICU's reported groups run from roughly 4% to 20% mortality. A
#' median raw AUPRC over 166 hospitals is therefore a median over 166 different
#' scales, and a hospital can rank above another purely by being sicker.
#' `auprc / event_rate` is the same number divided by its own floor, so 1.0
#' means "no better than guessing here" at every event rate. The raw AUPRC stays
#' on the per-group rows, where the event rate sits beside it.
#'
#' Small groups are EXCLUDED rather than reported with a wide interval. An
#' AUROC on 30 stays with 2 deaths is not a noisy estimate of that hospital's
#' discrimination; it is dominated by which two patients died. The thresholds
#' are arguments so they are declared at the call site and land in the manifest,
#' never chosen after seeing the spread.
#'
#' AGGREGATES ONLY (hard rule 1). One row per group, counts and an AUROC. No
#' group identifier is printed by this function; the caller decides whether a
#' hospital id may be shown, and for eICU it may not leave the run directory.
#'
#' @param group      grouping value per row, aligned to `score` and `y`
#' @param min_n      minimum stays for a group to be reported
#' @param min_events minimum deaths AND minimum survivors, both required —
#'                   AUROC is undefined when either class is empty and unstable
#'                   when either is tiny
#' @return list(per_group, summary). `summary` is one row: the distribution.
group_metrics <- function(score, y, group, label = NA_character_,
                          min_n = 100L, min_events = 10L) {
  stopifnot(length(score) == length(y), length(y) == length(group))
  n_scored <- length(y)
  ok <- !is.na(score) & !is.na(y) & !is.na(group)
  n_matched <- sum(ok)
  score <- score[ok]; y <- as.integer(y[ok]); group <- as.character(group[ok])

  rows <- lapply(sort(unique(group)), function(g) {
    s <- group == g
    n <- sum(s); k <- sum(y[s])
    keep <- n >= min_n & k >= min_events & (n - k) >= min_events
    data.frame(
      group = g, n = n, deaths = k, event_rate = round(k / n, 5),
      reported = keep,
      auroc = if (keep) round(.auroc(score[s], y[s]), 5) else NA_real_,
      auprc = if (keep) round(.auprc(score[s], y[s]), 5) else NA_real_,
      # Divided by the group's OWN event rate, which is what makes two groups
      # comparable. See the header.
      auprc_lift = if (keep) round(.auprc(score[s], y[s]) / (k / n), 5) else NA_real_,
      stringsAsFactors = FALSE)
  })
  per <- if (length(rows)) do.call(rbind, rows) else
    data.frame(group = character(0), n = integer(0), deaths = integer(0),
               event_rate = numeric(0), reported = logical(0), auroc = numeric(0),
               auprc = numeric(0), auprc_lift = numeric(0), stringsAsFactors = FALSE)
  rownames(per) <- NULL

  a  <- per$auroc[per$reported]
  al <- per$auprc_lift[per$reported]
  qs <- function(v) if (length(v))
    stats::quantile(v, c(0.10, 0.25, 0.50, 0.75, 0.90), names = FALSE)
    else rep(NA_real_, 5)
  # `min()` of an empty vector is Inf with a warning; an empty reported set is
  # a legitimate outcome and must read as NA (review E4).
  mn <- function(v) if (length(v)) min(v) else NA_real_
  mx <- function(v) if (length(v)) max(v) else NA_real_
  av <- function(v) if (length(v)) mean(v) else NA_real_
  pooled <- function(s, yy) {
    two <- length(yy) > 0L && any(yy == 1L) && any(yy == 0L)
    list(auroc = if (two) .auroc(s, yy) else NA_real_,
         auprc = if (two) .auprc(s, yy) else NA_real_,
         rate  = if (length(yy)) mean(yy) else NA_real_)
  }
  elig <- group %in% per$group[per$reported]
  pm <- pooled(score, y)
  pe <- pooled(score[elig], y[elig])
  q  <- qs(a); ql <- qs(al)
  n_kept <- sum(per$n[per$reported])
  status <- if (!n_matched) "no_matched_rows"
            else if (!any(per$reported)) "no_eligible_group" else "ok"

  summ <- data.frame(
    label          = label,
    status         = status,
    n_scored       = n_scored,
    n_matched      = n_matched,
    n_unmatched    = n_scored - n_matched,
    frac_matched   = if (n_scored) round(n_matched / n_scored, 4) else NA_real_,
    n_groups       = nrow(per),
    n_reported     = sum(per$reported),
    n_stays_kept   = n_kept,
    frac_kept_of_scored  = if (n_scored) round(n_kept / n_scored, 4) else NA_real_,
    frac_kept_of_matched = if (n_matched) round(n_kept / n_matched, 4) else NA_real_,
    pooled_auroc_matched  = round(pm$auroc, 5),
    pooled_auroc_eligible = round(pe$auroc, 5),
    auroc_mean     = round(av(a), 5),
    auroc_p10      = round(q[1], 5), auroc_q1  = round(q[2], 5),
    auroc_median   = round(q[3], 5), auroc_q3  = round(q[4], 5),
    auroc_p90      = round(q[5], 5),
    auroc_iqr      = round(q[4] - q[2], 5),
    auroc_min      = round(mn(a), 5), auroc_max = round(mx(a), 5),
    # Eligible pooled minus the median over the same eligible hospitals. A
    # descriptive contrast between two weightings of one cohort; see the header
    # for why it is not a case-mix decomposition.
    pooled_minus_median = round(pe$auroc - q[3], 5),
    # The same distribution on the precision-recall side. Each pooled lift is
    # over ITS OWN pooled event rate, which is not the mean of the per-group
    # floors, so `pooled_minus_median_lift` mixes two things and is a
    # description rather than a decomposition.
    pooled_auprc_matched       = round(pm$auprc, 5),
    pooled_event_rate_matched  = round(pm$rate, 5),
    pooled_auprc_lift_matched  = round(pm$auprc / pm$rate, 5),
    pooled_auprc_eligible      = round(pe$auprc, 5),
    pooled_event_rate_eligible = round(pe$rate, 5),
    pooled_auprc_lift_eligible = round(pe$auprc / pe$rate, 5),
    auprc_lift_mean   = round(av(al), 5),
    auprc_lift_p10    = round(ql[1], 5), auprc_lift_q1     = round(ql[2], 5),
    auprc_lift_median = round(ql[3], 5), auprc_lift_q3     = round(ql[4], 5),
    auprc_lift_p90    = round(ql[5], 5),
    auprc_lift_iqr    = round(ql[4] - ql[2], 5),
    auprc_lift_min    = round(mn(al), 5), auprc_lift_max  = round(mx(al), 5),
    pooled_minus_median_lift = round(pe$auprc / pe$rate - ql[3], 5),
    stringsAsFactors = FALSE)

  per <- per[order(-per$reported, per$auroc), , drop = FALSE]
  list(per_group = per, summary = summ)
}

# --- scoring a whole site ----------------------------------------------------

#' Score every arm of a site through one loop.
#'
#' `score_report()` already guarantees that any two scores are computed
#' identically. This guarantees that any two SITES are: MIMIC-test and eICU call
#' this function with the same argument list and differ only in which rows they
#' hand it. "The same code path" then stops being a claim about two runner
#' scripts that happen to look alike and becomes a property of the code.
#'
#' @param scores  named list of row-level log-odds scores, all aligned to `y`
#' @param breaks  named list of frozen cut points, one per arm, or NULL for each
#'                site's own quantiles. Pass the bundle's `cutpoints` to ask the
#'                transport question; pass NULL to ask the within-site one. Both
#'                are worth having and they answer different things.
#' @param suffix  appended to every label, so a run can hold both binnings
#'                without one overwriting the other's tables.
#' @return list(summary, reports)
#' @param group   patient id per row of `y`; the bootstrap unit (review S4)
score_arms <- function(run, scores, y, p_bar, breaks = NULL, suffix = "",
                       n_bins = 20L, n_boot = 200L, seed = 1L, group = NULL) {
  stopifnot(is.list(scores), length(scores) > 0L)
  bad <- names(scores)[vapply(scores, function(s) length(s) != length(y), logical(1))]
  if (length(bad)) abort_values("score_arms: arm(s) not aligned to `y`", bad)
  if (!is.null(group) && length(group) != length(y)) {
    stop("score_arms: `group` is not aligned to `y`", call. = FALSE)
  }

  # A REQUESTED FROZEN BINNING MUST BE AVAILABLE FOR EVERY ARM, and an arm it
  # does not cover is an error rather than a fallback.
  #
  # THE DEFECT THIS CLOSES. The call was `breaks[[nm]]`, and `[[` on a name a
  # list does not carry returns NULL, and `risk_bins(breaks = NULL)` falls back
  # to the SCORE'S OWN QUANTILES. So an arm missing from `bundle$cutpoints` was
  # silently binned on the site's own distribution while its label still said
  # `_frozen` and it still sat in a table headed "frozen". That is the exact
  # comparison the frozen binning exists to prevent: eICU quantiles re-centre
  # every bin on eICU, so bin 20 means "the sickest 5% of eICU" rather than "a
  # MIMIC-calibrated score above x", and nothing in the output would have said
  # which of the two a row was.
  #
  # It is reachable rather than theoretical. `cutpoints` is
  # `lapply(oof_scores, score_cutpoints)`, and `oof_scores` builds `llr_meas`
  # CONDITIONALLY and then drops NULLs -- so a run without the paired `meas`
  # fits produces cut points for nine arms, an apply site requests ten, and the
  # tenth is quietly self-binned.
  #
  # THE CLASS: `[[` RETURNING NULL WHERE NULL ALREADY MEANS SOMETHING ELSE.
  # `breaks = NULL` is a legitimate instruction ("bin on this site"), so the
  # absent-key NULL is indistinguishable from the deliberate one at the point
  # of use. The fix is to distinguish them at the point of ENTRY, where the
  # caller's intent is still visible.
  if (!is.null(breaks)) {
    miss <- setdiff(names(scores), names(breaks))
    if (length(miss)) {
      abort_values(paste0("score_arms: a frozen binning was requested but the ",
                          "cut points cover no such arm(s). Falling back to ",
                          "this site's own quantiles would label a self-binned ",
                          "row as frozen"), miss)
    }
  }

  reports <- list()
  for (nm in names(scores)) {
    lab <- paste0(nm, suffix)
    reports[[nm]] <- score_report(
      run, as.numeric(scores[[nm]]), y, p_bar, label = lab,
      n_bins = n_bins, n_boot = n_boot, seed = seed,
      breaks = if (is.null(breaks)) NULL else breaks[[nm]],
      title = sprintf("Risk ordering - %s", lab), group = group)
  }
  summ <- do.call(rbind, lapply(reports, function(z) z$summary))
  rownames(summ) <- NULL
  list(summary = summ, reports = reports)
}

#' Paired contrasts between arms, on identical rows.
#'
#' Marginal confidence intervals settle nothing when every arm is scored on the
#' same stays: two intervals can overlap while the paired difference is
#' unambiguous, and they can fail to overlap while the difference is not. So
#' every contrast here is PAIRED -- DeLong for AUROC, a paired bootstrap for
#' AUPRC -- which is the same discipline `tests/metrics_severity.R` applies to
#' the severity cells.
#'
#' THE AUROC INTERVAL IS THE CLUSTERED PAIRED BOOTSTRAP, as of 2026-09-09
#' (review S4). `auroc_lo`/`auroc_hi` used to be the DeLong interval, which is
#' observation-level. Both metrics now take their interval and their p-value
#' from ONE paired bootstrap that resamples patients -- `.paired_boot()` draws
#' each replicate once and evaluates both metrics on it -- and DeLong's z and
#' p travel beside them under `delong_*` as the iid reference. `auroc_p` and
#' `auprc_p` are the shifted-null bootstrap test of `.paired_boot()`.
#'
#' @param pairs list of length-2 character vectors, c(a, b). Reported as a - b.
#' @param group patient id per row of `y`; the bootstrap unit
#' @return data frame, one row per contrast
arm_contrasts <- function(scores, y, pairs, n_boot = 200L, seed = 1L, group = NULL) {
  rows <- lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    if (is.null(scores[[a]]) || is.null(scores[[b]])) {
      abort_values("arm_contrasts: a contrast names an arm that was not scored",
                   setdiff(p, names(scores)))
    }
    s1 <- as.numeric(scores[[a]]); s2 <- as.numeric(scores[[b]])
    dl <- delong_test(s1, s2, y)
    pbb <- .paired_boot(s1, s2, y, metrics = list(auroc = .auroc, auprc = .auprc),
                        n_boot = n_boot, seed = seed, group = group)
    pa <- pbb$auroc; pb <- pbb$auprc
    data.frame(
      a = a, b = b, n = dl[["n"]], n_events = dl[["n_events"]],
      auroc_a = dl[["auroc_1"]], auroc_b = dl[["auroc_2"]],
      d_auroc = dl[["delta"]], auroc_lo = pa[["ci_lo"]], auroc_hi = pa[["ci_hi"]],
      auroc_p = pa[["p_boot"]],
      delong_se = dl[["se"]], delong_z = dl[["z"]], delong_p = dl[["p_value"]],
      d_auprc = pb[["delta"]], auprc_lo = pb[["ci_lo"]], auprc_hi = pb[["ci_hi"]],
      auprc_p = pb[["p_boot"]],
      n_boot = pb[["n_boot"]], n_boot_ok = pb[["n_boot_ok"]],
      boot_unit = pb[["boot_unit"]],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
