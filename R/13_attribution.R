# R/13_attribution.R ---------------------------------------------------------
# The attributability arm: the common currency, and the agreement metrics.
#
# Read `docs/v2_attributability_plan_20260902.md` first. This file is steps 1
# and 2 of its build order; the noise floors (step 3) and the two-by-two
# (step 4) sit on top of it and live elsewhere.
#
# WHAT THIS ARM IS FOR, IN ONE SENTENCE. It asks whether the per-signal
# decomposition our method produces is a real quantity that can be estimated,
# rather than an artefact of one particular fit, and whether it survives a
# change of hospital system better than the attributions a strong learner
# produces over the same information.
#
# THE CHOICE OF ATTRIBUTION UNIT IS THE WHOLE ARGUMENT, and it is already
# settled by measurement. `v2_findings_20260827.md` section 5 bootstrapped the
# creatinine measurement model and found the two level smooths swinging by plus
# or minus 0.53 while the total linear predictor swung by 0.063 -- the terms
# trade fit, the prediction is stable, and PER-TERM ATTRIBUTION IS NOT
# IDENTIFIED. So the unit here is the signal-level `L`, never a smooth. `L` is
# not a coefficient; it is a PREDICTION, constrained by the data directly, which
# is why section 5.2 could call L well identified in the same breath as calling
# the terms inside it unidentified.
#
# FAITHFULNESS IS FREE AND IS NOT THE RESULT. For `llr_sum` the attributions sum
# to the score exactly, by construction. TreeSHAP is also exactly additive with
# respect to the tree's margin. Neither fact separates the methods and neither
# should be reported as a finding. What is NOT free, and what this arm measures,
# is whether the decomposition is IDENTIFIED under resampling and whether it
# TRANSPORTS.
#
# EVERYTHING KEYS TO ONE OBJECT. For a fixed evaluation cohort of n stays and
# the 19 modelled signals, every method must produce an n-by-19 matrix `A` where
# `A[i, g]` is signal g's contribution to stay i's predicted log-odds. A method
# that cannot be put on that index is not in the comparison.
#
# AGGREGATES ONLY (hard rule 1). `A` is row-level and belongs in a run
# directory. Every function that RETURNS a table here returns counts,
# correlations and quantiles.
# ----------------------------------------------------------------------------

#' The methods that can produce an attribution matrix.
#'
#' `llr_meas` is included because it is the APACHE II counterpart and because
#' the intervention block is exactly what differs between it and `llr_sum` --
#' so comparing their attribution matrices isolates what interventions
#' contribute to the decomposition, which no SHAP cell can do.
#'
#' `llr_cond` was added 2026-09-05 and it is the arm the interpretability claim
#' is actually about. `L_cond = L_full - L_intv` is the measurement-deviation
#' LLR CONDITIONAL on class and intervention context -- the first term of the
#' factorisation at the end of CLAUDE.md's frozen decisions -- and the sentence
#' the whole arm exists to license ("this measurement channel contributes X
#' nats, given what was being done to the patient") is a sentence about
#' `llr_cond`, not about `llr_sum`. `llr_sum` additionally carries
#' `rowSums(L_intv)`, a pure treatment-propensity contrast that is not a
#' statement about the measurement at all.
#'
#' IT COSTS ONE SUBTRACTION PER REPLICATE. `attribution_set()` already builds
#' every L matrix in one pass, so a per-cell interval for `llr_cond` needs no
#' extra fit -- which is why it can be carried through the floors alongside the
#' others rather than being a separate expensive run.
#' THE FOUR INTERACTION ARMS JOINED 2026-09-07, when `full_ti_trend` and
#' `full_ti_all` entered `LAYER1_MODELS` and their GAMs entered the bundle.
#' Before that they could be attributed only at MIMIC, by the generator, because
#' an apply site scores off frozen models and a model that is not in the bundle
#' does not exist there. So the attribution arm could report that
#' `llr_full_ti_trend` and `llr_full` disagree about a patient's leading signal
#' at the noise floor -- and could not ask whether that holds at eICU, which is
#' the question the whole transport claim is about.
#'
#' They cost NOTHING here. `attribution_set()` builds every L matrix in one
#' `l_matrices()` call, so the four extra arms are four extra sweeps of a matrix
#' that already exists. The expense was the 126 fits, and those are the
#' pipeline's.
ATTRIBUTION_METHODS <- c(names(LLR_ARM_MATRIX), "xgb_l", "xgb_feat", "xgb_raw")

# --- the common currency -----------------------------------------------------

#' Roll a SHAP matrix up from columns to signal groups.
#'
#' The column-to-group map is the one `xgb_feature_table()` already uses --
#' everything before the `__` -- so the rollup here and the gain rollup in
#' `xgb_group_gain()` cannot drift apart.
#'
#' Summing SHAP values within a group is EXACT. Shapley values are additive over
#' features by construction, and that single property is what makes this
#' comparison possible at all: `xgb_raw` splits a signal over about ten columns
#' and `xgb_feat` over three to five, so per-column attributions are not
#' comparable between designs while per-group ones are.
#'
#' @param intervention_handling what to do with columns whose group is an
#'   INTERVENTION rather than a signal. This is a real asymmetry and it is
#'   surfaced rather than silently resolved: our `L_full` for signal g already
#'   contains g's paired interventions, so our decomposition puts intervention
#'   effects INSIDE the signal, while the tree gives them their own columns.
#'
#'   "drop"  keep only the 19 signal groups and report the discarded share of
#'           total absolute attribution. Needs no convention, and the discarded
#'           share is itself a reportable quantity. THE DEFAULT.
#'   "split" divide each intervention group's attribution equally among the
#'           signals it is paired with in `pairing.csv`, making the
#'           decomposition complete on both sides. More faithful to what our
#'           method does, but "equally" is a CONVENTION and must be declared as
#'           one wherever a number produced this way is reported.
#' @return list(A, dropped_frac)
.shap_rollup <- function(S, feature_names, signals, cfg,
                         intervention_handling = c("drop", "split")) {
  intervention_handling <- match.arg(intervention_handling)
  grp <- ifelse(grepl("__", feature_names, fixed = TRUE),
                sub("__.*$", "", feature_names), feature_names)
  A <- matrix(0, nrow = nrow(S), ncol = length(signals),
              dimnames = list(rownames(S), signals))
  for (sg in signals) {
    j <- which(grp == sg)
    if (length(j)) A[, sg] <- rowSums(S[, j, drop = FALSE])
  }
  other <- which(!grp %in% signals)
  dropped <- if (length(other)) sum(abs(S[, other, drop = FALSE])) else 0
  total   <- sum(abs(S))

  if (identical(intervention_handling, "split") && length(other)) {
    pr <- cfg$pairing
    for (iv in unique(grp[other])) {
      j <- which(grp == iv)
      tgt <- intersect(unique(pr$signal[pr$intervention == iv]), signals)
      if (!length(tgt)) next
      v <- rowSums(S[, j, drop = FALSE]) / length(tgt)
      for (sg in tgt) A[, sg] <- A[, sg] + v
    }
  }
  list(A = A, dropped_frac = if (total > 0) dropped / total else 0)
}

#' Every method's attribution matrix, computed in ONE pass over layer 1.
#'
#' Layer 1 is evaluated once and shared, because recomputing 43 smooths per
#' method would dominate the runtime of every floor replicate downstream -- and
#' the floors are the expensive part of this arm.
#'
#' THE BACKGROUND DISTRIBUTION, DECLARED. TreeSHAP here is `predcontrib`, which
#' uses the tree's own cover-weighted path-dependent expectation. The background
#' therefore TRAVELS WITH THE BOOSTER rather than being sampled at the
#' evaluation site. That is the honest and the convenient choice and it must be
#' stated wherever these numbers appear. For the `xgb_l` cell the implied
#' background is the MIMIC L distribution, which shifts at eICU; that shift is
#' part of the finding and is not a nuisance to be corrected away.
#'
#' @param w optional named weight vector over signals, the `w_g` of the design.
#'   NULL means equal weights, which is what `apply_bundle()` currently uses --
#'   `bundle$layer2` is NULL and the Sigma-inverse weights do not exist yet. It
#'   is a PARAMETER from the start so that the weighted aggregate can be run
#'   through the identical code path the day the weights land.
#' @return list(A = named list of matrices, measured, signals, w,
#'   dropped_frac, l_mats)
attribution_set <- function(bundle, tabs, cfg, stay_ids,
                            methods = ATTRIBUTION_METHODS, w = NULL,
                            intervention_handling = "drop", verbose = TRUE) {
  bad <- setdiff(methods, ATTRIBUTION_METHODS)
  if (length(bad)) abort_values("attribution_set: unknown method(s)", bad)
  stamp <- cfg$.bundle_design
  if (is.null(stamp) || !identical(stamp, .hash(bundle$cfg))) {
    stop("attribution_set: `cfg` did not come from bundle_cfg(bundle, paths). ",
         "An attribution computed under a design the bundle did not freeze is ",
         "not an attribution of that bundle (hard rule 8).", call. = FALSE)
  }
  signals <- as.character(unlist(cfg$signals))
  ids     <- as.character(stay_ids)
  ww      <- if (is.null(w)) stats::setNames(rep(1, length(signals)), signals) else w
  miss_w  <- setdiff(signals, names(ww))
  if (length(miss_w)) abort_values("attribution_set: `w` is missing signal(s)", miss_w)

  if (verbose) message("attribution_set: layer 1 over ", length(ids), " stays")
  a1 <- apply_layer1(bundle$models, tabs, cfg, bundle$priors, stay_ids, verbose = verbose)
  mats <- l_matrices(a1$l, tabs, cfg, stay_ids, fill = "zero")

  out <- list(); drop_frac <- list()
  wv <- ww[signals]

  # ONE LOOP OVER `LLR_ARM_MATRIX`, the same map `apply_bundle()` reads, so a
  # score and its attribution matrix can never be built from different L
  # matrices under one arm name. The `cond` family is derived and never fitted:
  # `l_matrices()` returns it as `full - intv`, and the weighting is applied to
  # the conditional column for the same reason it is applied to the joint one --
  # `w_g` weights a SIGNAL, and the arm differs only in which of that signal's
  # L's is being weighted.
  for (nm in intersect(names(LLR_ARM_MATRIX), methods)) {
    mt <- LLR_ARM_MATRIX[[nm]]
    if (is.null(mats[[mt]])) {
      stop("attribution_set: no `", mt, "` L matrix for arm `", nm, "`",
           call. = FALSE)
    }
    out[[nm]] <- sweep(mats[[mt]][, signals, drop = FALSE], 2L, wv, `*`)
    drop_frac[[nm]] <- 0
  }

  for (nm in intersect(XGB_DESIGNS, methods)) {
    m <- bundle$xgb[[nm]]
    if (is.null(m)) stop("attribution_set: the bundle carries no `", nm, "` booster",
                         call. = FALSE)
    if (verbose) message("attribution_set: TreeSHAP for ", nm)
    X <- switch(nm,
      xgb_l    = xgb_design_L(a1$l, tabs, cfg, stay_ids, model = "full", fill = "zero"),
      xgb_feat = xgb_design_feat(tabs, cfg, bundle$priors, stay_ids, role = "final",
                                 fold = NA_integer_, feature_names = m$feature_names),
      xgb_raw  = xgb_design_raw(tabs, cfg, stay_ids, feature_names = m$feature_names))
    Xa <- align_design(X, m$feature_names)
    S  <- stats::predict(m$booster, Xa, predcontrib = TRUE)
    # `predcontrib` appends a BIAS column. It is the model's base value, not a
    # feature's contribution, and including it in a rollup would add a constant
    # to whichever group happened to sort last.
    if (ncol(S) == length(m$feature_names) + 1L) S <- S[, seq_along(m$feature_names), drop = FALSE]
    if (ncol(S) != length(m$feature_names)) {
      stop("attribution_set: predcontrib returned ", ncol(S), " columns for ",
           length(m$feature_names), " features in `", nm, "`", call. = FALSE)
    }
    rownames(S) <- ids
    r <- .shap_rollup(S, m$feature_names, signals, cfg,
                      intervention_handling = intervention_handling)
    out[[nm]] <- r$A
    drop_frac[[nm]] <- r$dropped_frac
  }

  for (nm in names(out)) {
    if (!identical(dim(out[[nm]]), c(length(ids), length(signals)))) {
      stop("attribution_set: `", nm, "` is not n x 19", call. = FALSE)
    }
  }
  list(A = out, measured = measured_matrix(tabs, cfg, stay_ids),
       signals = signals, w = wv, dropped_frac = unlist(drop_frac),
       intervention_handling = intervention_handling, l_mats = mats)
}

# --- the metrics -------------------------------------------------------------

#' Which stay-by-signal cells a comparison may use.
#'
#' ASSIGNED ZEROS ARE EXCLUDED, and this is the single most consequential
#' bookkeeping decision in the arm. An unmeasured signal gets `L = 0` by
#' assignment at both sites, so those cells agree PERFECTLY between any two fits
#' without any estimation having happened. Counting them would inflate every
#' agreement score our method reports, by an amount that depends on the 8.7%
#' unmeasured fraction at MIMIC and on eICU's quite different coverage profile.
#'
#' SHAP has no analogue: a tree assigns a nonzero contribution to a missing
#' feature through its default direction, so its cells are never trivially
#' equal. Restricting both methods to cells measured under BOTH fits is what
#' makes the comparison like-for-like.
#'
#' The excluded fraction is returned, never assumed away.
attribution_cells <- function(m1, m2 = NULL) {
  keep <- if (is.null(m2)) m1 else (m1 & m2)
  list(keep = keep, frac_used = mean(keep), n_used = sum(keep))
}

#' Metric one and two: attribution agreement, and sign agreement.
#'
#' @param mag_floor absolute log-odds magnitude below which a cell's SIGN is
#'   noise and is excluded from the sign metric. DECLARE IT BEFORE THE RUN and
#'   state it in absolute units; choosing it as a quantile after seeing the
#'   distribution turns the metric into a description of how many near-zero
#'   cells there are.
#' @return list(overall, per_signal)
attribution_agreement <- function(A1, A2, keep, mag_floor = 0.05) {
  stopifnot(identical(dim(A1), dim(A2)), identical(dim(A1), dim(keep)))
  sg <- colnames(A1)
  sp <- function(a, b) {
    if (length(a) < 3L || stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    suppressWarnings(stats::cor(a, b, method = "spearman"))
  }
  pe <- function(a, b) {
    if (length(a) < 3L || stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    suppressWarnings(stats::cor(a, b))
  }
  v1 <- A1[keep]; v2 <- A2[keep]
  big <- pmax(abs(v1), abs(v2)) > mag_floor
  overall <- data.frame(
    n_cells = sum(keep), frac_cells_used = round(mean(keep), 5),
    spearman = round(sp(v1, v2), 5), pearson = round(pe(v1, v2), 5),
    n_sign_cells = sum(big),
    frac_above_floor = round(mean(big), 5),
    sign_agreement = round(if (any(big)) mean(sign(v1[big]) == sign(v2[big])) else NA_real_, 5),
    mag_floor = mag_floor, stringsAsFactors = FALSE)

  per <- do.call(rbind, lapply(seq_along(sg), function(j) {
    k <- keep[, j]
    a <- A1[k, j]; b <- A2[k, j]
    bg <- pmax(abs(a), abs(b)) > mag_floor
    data.frame(signal = sg[j], n_cells = sum(k),
               frac_measured = round(mean(k), 5),
               spearman = round(sp(a, b), 5), pearson = round(pe(a, b), 5),
               sign_agreement = round(if (any(bg)) mean(sign(a[bg]) == sign(b[bg])) else NA_real_, 5),
               mean_abs_1 = round(mean(abs(a)), 5), mean_abs_2 = round(mean(abs(b)), 5),
               stringsAsFactors = FALSE)
  }))
  list(overall = overall, per_signal = per[order(per$spearman), , drop = FALSE])
}

#' Metric four: within-patient bundle ranking.
#'
#' For each stay, rank the signals by `|A|` under each fit and correlate the two
#' orderings. Reported as a DISTRIBUTION over stays rather than a mean, because
#' the interesting failure is a subpopulation whose ranking scrambles, and a
#' mean hides it.
#'
#' Stays with fewer than `min_signals` usable cells are excluded and counted: a
#' Spearman over two signals is not a noisy estimate of agreement, it is a coin.
attribution_ranking <- function(A1, A2, keep, min_signals = 5L) {
  n <- nrow(A1); rho <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    j <- which(keep[i, ])
    if (length(j) < min_signals) next
    a <- abs(A1[i, j]); b <- abs(A2[i, j])
    if (stats::sd(a) == 0 || stats::sd(b) == 0) next
    rho[i] <- suppressWarnings(stats::cor(a, b, method = "spearman"))
  }
  ok <- !is.na(rho)
  q <- if (any(ok)) stats::quantile(rho[ok], c(0.1, 0.25, 0.5, 0.75, 0.9)) else rep(NA_real_, 5)
  data.frame(n_stays = n, n_scored = sum(ok),
             frac_scored = round(mean(ok), 5), min_signals = min_signals,
             p10 = round(q[1], 4), q1 = round(q[2], 4), median = round(q[3], 4),
             q3 = round(q[4], 4), p90 = round(q[5], 4),
             frac_above_0.8 = round(if (any(ok)) mean(rho[ok] > 0.8) else NA_real_, 4),
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Metric five: domain composition within a risk stratum.
#'
#' Uses the FROZEN partition in `config/domains.csv`, so the aggregation is the
#' one declared before any result existed. Within each quantile bin of the
#' score, each domain's share of total POSITIVE evidence. The most clinically
#' legible thing this arm produces, and the one where every SOFA organ has an
#' external referent.
attribution_domains <- function(A, score, domains, n_bins = 4L) {
  sg  <- colnames(A)
  dmn <- domains[domains$signal %in% sg, , drop = FALSE]
  b   <- cut(score, breaks = stats::quantile(score, seq(0, 1, length.out = n_bins + 1L)),
             include.lowest = TRUE, labels = FALSE)
  pos <- pmax(A, 0)
  rows <- lapply(sort(unique(b)), function(bb) {
    i <- which(b == bb)
    tot <- sum(pos[i, , drop = FALSE])
    do.call(rbind, lapply(split(dmn$signal, dmn$domain), function(s) {
      s <- intersect(s, sg)
      data.frame(bin = bb, n_stays = length(i),
                 domain = dmn$domain[match(s[1], dmn$signal)],
                 n_signals = length(s),
                 share_positive = round(if (tot > 0)
                   sum(pos[i, s, drop = FALSE]) / tot else NA_real_, 5),
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out[order(out$bin, -out$share_positive), , drop = FALSE]
}

#' Metric three: the marginal evidence curve.
#'
#' THE FIGURE THIS ARM SHOULD LEAD WITH, because it is the only output a
#' clinician reads directly. For signal g, the average attribution as a function
#' of the observed covariate value, on a common grid. For our method that is
#' essentially the fitted evidence curve; for SHAP it is the standard dependence
#' plot. Both are computable at both sites from both methods, which is exactly
#' what the per-smooth version in `v2_analytical_design_plan.md` could not be --
#' section 5.2 of the findings document ruled that one out because per-term
#' attribution is not identified.
#'
#' The grid is passed in rather than derived, so two sites are compared on the
#' SAME grid and neither re-centres it on its own distribution -- the identical
#' argument the frozen reporting bins make in `risk_bins()`.
#'
#' @param x per-stay covariate value for this signal, aligned to `a`
#' @param grid common bin edges. `evidence_grid()` builds one from a reference.
evidence_curve <- function(a, x, keep, grid) {
  ok <- keep & !is.na(x) & !is.na(a)
  bin <- cut(x, breaks = grid, include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(seq_len(length(grid) - 1L), function(b) {
    i <- which(ok & bin == b)
    data.frame(bin = b, lo = grid[b], hi = grid[b + 1L], n = length(i),
               mean_a = if (length(i)) round(mean(a[i]), 5) else NA_real_,
               sd_a   = if (length(i) > 1L) round(stats::sd(a[i]), 5) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' A common grid from a reference distribution. Frozen at the training site.
evidence_grid <- function(x, n_bins = 20L) {
  q <- stats::quantile(x[!is.na(x)], seq(0, 1, length.out = n_bins + 1L))
  q <- unique(as.numeric(q))
  if (length(q) < 3L) return(range(x, na.rm = TRUE))
  q[1] <- -Inf; q[length(q)] <- Inf
  q
}

#' Agreement between two evidence curves, density-weighted.
#'
#' Weighted by the evaluation density so that regions no patient occupies do not
#' dominate a correlation over an arbitrary grid.
evidence_curve_agreement <- function(c1, c2) {
  ok <- !is.na(c1$mean_a) & !is.na(c2$mean_a) & c1$n > 0 & c2$n > 0
  if (sum(ok) < 3L) return(data.frame(n_bins = sum(ok), weighted_cor = NA_real_,
                                      l2_norm = NA_real_, stringsAsFactors = FALSE))
  wt <- pmin(c1$n[ok], c2$n[ok]); wt <- wt / sum(wt)
  a <- c1$mean_a[ok]; b <- c2$mean_a[ok]
  ma <- sum(wt * a); mb <- sum(wt * b)
  va <- sum(wt * (a - ma)^2); vb <- sum(wt * (b - mb)^2)
  cv <- sum(wt * (a - ma) * (b - mb))
  sp <- stats::sd(c(a, b))
  data.frame(n_bins = sum(ok),
             weighted_cor = round(if (va > 0 && vb > 0) cv / sqrt(va * vb) else NA_real_, 5),
             # Normalised so it is comparable across signals on different scales.
             l2_norm = round(if (sp > 0) sqrt(sum(wt * (a - b)^2)) / sp else NA_real_, 5),
             stringsAsFactors = FALSE)
}

# --- the one call a floor replicate makes ------------------------------------

#' Every metric between two attribution sets, in one table each.
#'
#' This is what a bootstrap replicate, a disjoint-half comparison and a
#' cross-site comparison all call, so all three produce identically shaped
#' output and can be stacked into one table with a `perturbation` column. That
#' is what makes a floor and a finding directly comparable, which is the whole
#' point of having a floor.
compare_attributions <- function(s1, s2, method, mag_floor = 0.05,
                                 min_signals = 5L, label = NA_character_) {
  A1 <- s1$A[[method]]; A2 <- s2$A[[method]]
  if (is.null(A1) || is.null(A2)) {
    abort_values("compare_attributions: method absent from one of the sets", method)
  }
  cells <- attribution_cells(s1$measured, s2$measured)
  ag <- attribution_agreement(A1, A2, cells$keep, mag_floor = mag_floor)
  rk <- attribution_ranking(A1, A2, cells$keep, min_signals = min_signals)
  list(overall    = cbind(label = label, method = method, ag$overall),
       per_signal = cbind(label = label, method = method, ag$per_signal),
       ranking    = cbind(label = label, method = method, rk))
}
