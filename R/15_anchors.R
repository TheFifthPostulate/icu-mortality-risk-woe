# R/15_anchors.R ---------------------------------------------------------------
# Reference anchors for DISPLAYING the terms of a fitted evidence model.
#
# A fitted smooth f_j is identified only up to a constant: mgcv centers it on
# the training rows and the intercept absorbs the rest. For display, each term
# is re-anchored at a declared reference value x0_j,
#
#     f*_j(x) = f_j(x) - f_j(x0_j),
#
# and the constants move into the REFERENCE EVIDENCE of the model,
#
#     L_ref = intercept + sum_j f_j(x0_j) - logit(p_bar),
#
# the weight of evidence of a stay at every anchor. For any stay,
# L = L_ref + sum_j f*_j(x_j) (+ parametric terms, anchored at 0): every weight
# of evidence, score, band and stratum is unchanged. Presentation only; nothing
# here is fitted, and nothing here enters a score.
#
# Anchors ("a quiet day on this channel"):
#   pi_minus, pi_plus   alpha_side / (n_anchor + alpha_total): the shrunken
#                       propensity of a stay with no deviation on that side at the
#                       channel's typical charting density, n_anchor = the median
#                       number of observed hours among measured MIMIC-IV training
#                       stays (paper/make_anchor_n.R). Exactly 0 is not an
#                       achievable propensity under the Dirichlet shrinkage, and a
#                       fixed n = 24 lies outside the data of sparse channels.
#   *_delta             0: the tail as extreme as the deviation count predicts;
#                       0 for the Glasgow components too, where the ordinal
#                       residual clusters on a few values and 0 can fall between
#                       supported clusters: the curve there is interpolated.
#   trend               0: flat
#   *__exposure_frac, *__n_hours   0: not exposed
#   *__lambda           0: as intense as the exposure predicts (also the value
#                       assigned to an unexposed stay)
#   parametric terms    0: not present at admission
# ------------------------------------------------------------------------------

#' Reference value of one model variable.
#' @param var variable name as it appears in the formula (smooth `term`)
#' @param pri a priors_for() row for the model's signal
#' @param n_anchor the channel's typical number of observed hours (frozen,
#'   passed in by the caller; this file reads no paths)
term_anchor <- function(var, pri, n_anchor) {
  a <- pri$alpha   # (low, mid, high)
  if (var %in% c("pi_minus", "pi_plus")) {
    if (missing(n_anchor) || !is.finite(n_anchor) || n_anchor < 1) stop("term_anchor: n_anchor required for ", var, call. = FALSE)
    return(if (var == "pi_minus") a[1] / (n_anchor + sum(a)) else a[3] / (n_anchor + sum(a)))
  }
  0
}

#' Basis rows of one smooth at the values `x`, via mgcv::PredictMat.
.smooth_rows <- function(s, x) {
  mgcv::PredictMat(s, stats::setNames(data.frame(x), s$term[1]))
}

#' Anchored curve of one smooth on a grid: value and pointwise SE of
#' f(x) - f(x0), exact under the model's posterior covariance.
anchored_curve <- function(b, s, grid, x0) {
  idx <- s$first.para:s$last.para
  V   <- if (!is.null(b$Vc)) b$Vc else b$Vp
  D   <- sweep(.smooth_rows(s, grid), 2, as.numeric(.smooth_rows(s, x0)))
  list(fit = as.numeric(D %*% stats::coef(b)[idx]),
       se  = sqrt(pmax(rowSums((D %*% V[idx, idx, drop = FALSE]) * D), 0)))
}

#' Reference evidence of one model: the weight of evidence at every anchor,
#' with its SE. Returns the coefficient-space row `r` so a caller can evaluate
#' posterior draws of L_ref on the same draws as the channel (draws %*% r).
reference_evidence <- function(b, pri, n_anchor) {
  cf <- stats::coef(b)
  r  <- numeric(length(cf)); r[match("(Intercept)", names(cf))] <- 1
  for (s in b$smooth) {
    idx <- s$first.para:s$last.para
    r[idx] <- r[idx] + as.numeric(.smooth_rows(s, term_anchor(s$term[1], pri, n_anchor)))
  }
  V <- if (!is.null(b$Vc)) b$Vc else b$Vp
  list(value = sum(r * cf) - logit(pri$p_bar), se = sqrt(as.numeric(t(r) %*% V %*% r)), row = r)
}
