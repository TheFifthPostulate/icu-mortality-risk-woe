# R/04b_conditional.R --------------------------------------------------------
# Conditional-prior shrinkage: the same estimator as the Dirichlet-multinomial
# in R/04_features.R, with a prior mean that depends on a covariate instead of
# being a constant.
#
# THE ONE IDEA IN THIS FILE.
#
#   pi_hat_ij = (k_ij + alpha_j) / (n_i + alpha_0)
#
# is a precision-weighted blend of what THIS stay did (k_ij / n_i) with what the
# POPULATION does (alpha_j / alpha_0), weighted by n_i / (n_i + alpha_0). The
# population target is a CONSTANT: every stay is shrunk toward the same number.
#
# Both estimators here keep the blend and replace that constant with a fitted
# conditional mean. The population target becomes "what a stay like this one, in
# the one respect that mechanically drives the quantity, normally shows".
#
#   delta   the MAGNITUDE construct.  How far past its reference bound did this
#           signal go, relative to what its DEVIATION COUNT predicts?
#   lambda  the INTENSITY construct.  How hard was this stay treated, relative
#           to what its EXPOSURE DURATION predicts?
#
# WHY EACH IS NEEDED. In both cases the raw pair of covariates is two views of
# one object, and the second view is mechanically contaminated by the first:
#
#   delta   `value_min` is an order statistic. Draw 20 deviant hours and you
#           will reach further into the tail than if you draw 1, for reasons
#           that are sampling and not physiology. So the raw magnitude term
#           confounds "how sick" with "how many deviant hours" -- and the second
#           is exactly what pi_hat already carries. MEASURED 2026-08-28: 40-60%
#           of the variance in the magnitude term is predictable from the counts
#           alone (creatinine 0.512, mbp 0.604, bun 0.559, lactate 0.478,
#           temperature 0.401).
#   lambda  A stay on pressors for 20 hours has had more OPPORTUNITY to be
#           escalated to a second and third molecule than one on pressors for 4.
#           The interesting patient is the one with high intensity and LOW
#           duration, and additive `s(D) + s(A)` cannot express that.
#
# WHY NOT JUST DIVIDE. The obvious fix -- a rate, A / D -- is worse than doing
# nothing, and MEASURED 2026-08-28 in both directions. For vasopressor,
# r(D, A) = +0.398 becomes r(A/D, D) = -0.676; for inotrope, where there is no
# coupling at all (+0.012), the rate manufactures one of -0.729. Dividing by D
# only standardises correctly when A rises exactly in step with D, and the
# fitted exponent ranges from 0.004 (inotrope) to 1.686 (diuretic) across the
# five (D, A) pairs. No fixed arithmetic transform is right for more than one
# of them. Same failure as the `spread = |tail - median|` composite already
# rejected in CLAUDE.md: a raw transform swaps one dependence for another.
#
# WHY EXPOSURE ENTERS AS A PRECISION AND NOT ONLY AS A DIVISOR. With one (A, D)
# per stay the two variance components -- within-stay noise, which shrinks as
# 1/D, and genuine between-stay heterogeneity, which does not -- are separated
# ONLY by the spread of D across stays. That is the same trick that identifies
# alpha in the Dirichlet-multinomial, where the spread of n_obs separates
# multinomial sampling noise from between-stay heterogeneity. Writing A/D throws
# that structure away and the model is no longer identified.
#
# FITTED ON TRAINING ROWS ONLY, OUTCOME-BLIND. Every parameter here is estimated
# from covariates alone; mortality is never referenced. Same discipline as
# alpha, and for the same reason -- these are descriptions of how treatment and
# physiology are distributed, not of how they relate to the outcome.
#
# FROZEN AT FIT TIME, RECOMPUTED PER ROW AT APPLY TIME (hard rule 8). The
# parameters go in the bundle beside alpha; the shrunk value is then recomputed
# from each row's own counts wherever the bundle is applied. That is why the
# apply functions take parameters as an argument rather than deriving them: they
# cannot then standardise a held-out fold using that fold's own data.
#
# No paths, no clock, no data access beyond the frames handed in (hard rule 9).
# ----------------------------------------------------------------------------

# --- shared variance-component fit ------------------------------------------

#' Fit Var(r_i) = s_e / w_i + s_u by maximum likelihood.
#'
#' The heart of both constructs. `r` is a residual from a fitted conditional
#' mean and `w` is the stay's information weight -- deviant hours for `delta`,
#' treated hours for `lambda`. Two variance components:
#'
#'   s_e   within-stay noise. Scaled by 1/w, so it vanishes for a stay with a
#'         lot of data and dominates for a stay with almost none.
#'   s_u   between-stay heterogeneity. The real signal. Constant in w.
#'
#' They are separately identifiable ONLY because w varies across stays. If every
#' stay had the same w the likelihood would see s_e/w + s_u as one number.
#'
#' The posterior mean of the stay-level effect is then the precision-weighted
#' shrinkage that both apply functions use:
#'
#'   shrunk_i = r_i * s_u / (s_u + s_e / w_i)
#'
#' which is r_i when w_i is large (trust the stay) and 0 when w_i is small
#' (fall back on the population). Exactly the role n_i plays in pi_hat.
#'
#' THE SPLIT IS NOT ALWAYS IDENTIFIED, AND AN UNIDENTIFIED SPLIT IS SILENT.
#'
#' When w barely varies -- which happens whenever a signal's modelled tail is
#' nearly empty, so 1 + k is 1 for almost every stay -- the likelihood is flat
#' along s_e + s_u = Var(r). The optimiser then lands somewhere arbitrary on
#' that ridge, and if it lands near s_u = 0 the shrunk value is ZERO for every
#' stay. The column is not constant, so check_model_frame passes it; the smooth
#' is penalised away; the term silently carries nothing.
#'
#' MEASURED 2026-08-28, and this is exactly how it was caught: lactate's low
#' tail is occupied on 53 of 21,860 measured stays, and `value_min_delta` came
#' back with sd = 2.8e-06 across a range of 1e-05.
#'
#' The guard is model selection rather than a threshold on the data: fit both
#' the one-component model (s_e = 0, no shrinkage) and the two-component one,
#' and keep the two-component fit only if it earns its extra parameter on AIC.
#' Where w does not vary the likelihoods are equal, the simpler model wins on
#' the parameter count, and the fallback is NO SHRINKAGE -- which leaves the
#' conditional standardisation intact and drops only the part the data cannot
#' estimate. No tuning constant enters anywhere.
#'
#' @param r residuals from the conditional mean
#' @param w positive information weights, same length as r
#' @return list(s_e, s_u, converged, shrinks). `shrinks` is FALSE when the
#'   one-component model won, i.e. the returned value is the raw residual.
.fit_var_components <- function(r, w) {
  vr <- stats::var(r)
  if (!is.finite(vr) || vr <= 0) {
    return(list(s_e = 0, s_u = 0, converged = FALSE, shrinks = FALSE))
  }
  nll <- function(s_e, s_u) {
    v <- s_e / w + s_u
    if (any(!is.finite(v)) || any(v <= 0)) return(1e12)
    0.5 * sum(log(v) + r^2 / v)
  }
  # One component: all variance is between-stay, so nothing is shrunk away.
  aic1 <- 2 * nll(0, vr) + 2 * 1

  o <- try(stats::optim(c(log(vr * stats::median(w) / 2), log(vr / 2)),
                        function(p) nll(exp(p[1]), exp(p[2])),
                        method = "Nelder-Mead",
                        control = list(maxit = 2000, reltol = 1e-10)),
           silent = TRUE)
  if (inherits(o, "try-error") || o$convergence != 0L) {
    return(list(s_e = 0, s_u = vr, converged = FALSE, shrinks = FALSE))
  }
  s_e <- exp(o$par[1]); s_u <- exp(o$par[2])
  aic2 <- 2 * nll(s_e, s_u) + 2 * 2

  # Second guard, for the case AIC cannot see: a BOUNDARY solution at s_u -> 0.
  # Variance-component fits on the boundary are unreliable by construction, and
  # here the consequence is specific and silent -- every shrunk value becomes
  # numerically zero. Standard treatment is to fall back to the reduced model.
  # This is a numerical-zero test against the residual variance, not a
  # substantive threshold on how much heterogeneity counts as real.
  on_boundary <- s_u <= 1e-6 * vr
  if (!is.finite(aic2) || aic2 >= aic1 || on_boundary) {
    return(list(s_e = 0, s_u = vr, converged = TRUE, shrinks = FALSE))
  }
  list(s_e = s_e, s_u = s_u, converged = TRUE, shrinks = TRUE)
}

# `.shrink_resid()` LIVED HERE AND WAS DELETED 2026-09-05. It applied the
# precision weight `r * s_u / (s_u + s_e/w)` to a delta or lambda residual, and
# it is gone because that weight is gone from both constructs
# (docs/covariate_constructs_20260905.md, and the frozen decision in CLAUDE.md).
#
# THE ONE-LINE REASON: the weight is monotone INCREASING in the information
# `w`, so it pushed a high-coverage stay further from zero than a low-coverage
# stay carrying the identical residual -- which made the covariate depend MORE
# on coverage, not less. Measured over the 38 delta rows, its entire deviance
# gain was that coverage: +3.379 percentage points unconditionally, -1.591 with
# `s(log n_obs)` in both models. It reversed.
#
# `.fit_var_components()` and `.fit_var_components_rep()` below are NOT dead and
# must not be deleted with it. They still estimate `s_e` and `s_u`, which are
# reported as diagnostics: within-stay measurement noise is a real quantity
# worth publishing. We measure it; we no longer rescale with it.

#' Variance components from REPLICATES rather than from heteroscedasticity.
#'
#' ADDED 2026-09-02, and it replaces `.fit_var_components()` for `delta`
#' wherever the replicate statistics are available. The reason is a
#' measurement, not a preference.
#'
#' `.fit_var_components()` has to INFER `s_e` from the way residual variance
#' changes with `n`, using one residual per stay. That is a second-order signal
#' and it turned out not to be identified: re-running it on the real residuals
#' for all 38 delta rows, 19 returned `s_e` between 1e-9 and 1e-6 with a
#' likelihood IDENTICAL to the one-component model -- `aic2 - aic1` came back at
#' exactly +2.000, the parameter penalty and nothing else. The shrinkage weight
#' was 1.000 at every `n` and the machinery was inert, including on all three
#' Glasgow components, which is where a between-site coverage difference
#' actually bites. No selection criterion could have rescued that: the
#' likelihood is flat, so BIC, a likelihood-ratio test and cross-validation all
#' tie as well. It was never a model-selection problem.
#'
#' The mechanism is that the residual is ALREADY conditioned on `n` through
#' `a2 log(n)` in the conditional mean, which absorbs the coverage dependence in
#' LOCATION. Little second-order scale dependence is left for `s_e` to find.
#'
#' So `s_e` is estimated here from replicates the hourly lattice already
#' contains, carried into the feature table as sufficient statistics by
#' `v2_05_features_*.sql` (`se_within_ss`, `se_within_n`). Pooled across stays:
#'
#'     s_e = sum(se_within_ss) / sum(se_within_n - 1)
#'
#' which is the within-group mean square of a one-way random-effects model --
#' exactly what `s_e` was always meant to be, and estimable whether or not
#' residual variance happens to vary with `n`.
#'
#' `s_u` then follows by moments from Var(r) = s_u + s_e * E[1/n].
#'
#' EVERY DEGENERATE CASE FALLS BACK TO NO SHRINKAGE AND SAYS WHY. `se_source`
#' is returned so the diagnostics table records which rows are shrinking on
#' replicates, which are falling back, and for what reason -- the previous
#' estimator's silence about that is how the inert machinery went unnoticed.
#'
#' @param r_val   residual on the VALUE scale, the same scale the replicate
#'   statistics are on. See `fit_delta_ordinal()` for why the ordinal form
#'   builds an auxiliary value-scale residual rather than passing its own.
#' @param n       per-stay information weight (hours), aligned to `r_val`
#' @param within_ss,within_n  per-stay sum of squared deviations of the hourly
#'   values and the number of hours contributing. NULL or all-NA means the
#'   parquet predates the columns, and the caller falls back.
.fit_var_components_rep <- function(r_val, n, within_ss, within_n) {
  out <- function(s_e, s_u, conv, shr, src, df) list(
    s_e = s_e, s_u = s_u, converged = conv, shrinks = shr,
    se_source = src, se_df = df)

  ok <- is.finite(within_ss) & is.finite(within_n) & within_n >= 2
  df <- if (any(ok)) sum(within_n[ok] - 1) else 0
  vr <- stats::var(r_val)

  # No stay contributed two hours, so there is no within-stay information at
  # all. Not a failure -- a signal measured once per stay genuinely has no
  # replicate structure -- but it must not silently shrink.
  if (!any(ok) || df <= 0) {
    return(out(0, if (is.finite(vr)) vr else 0, FALSE, FALSE, "no_replicates", 0))
  }
  if (!is.finite(vr) || vr <= 0) return(out(0, 0, FALSE, FALSE, "degenerate_residual", df))

  s_e <- sum(within_ss[ok]) / df

  # s_e == 0 is a legitimate answer, not an error: the value never moves within
  # a stay, so there is no measurement noise to shrink away and a weight of 1 is
  # CORRECT. Distinguished from the old estimator, where the same number meant
  # "could not be identified".
  if (!is.finite(s_e) || s_e <= 0) return(out(0, vr, TRUE, FALSE, "no_within_variation", df))

  s_u <- vr - s_e * mean(1 / n)

  # The boundary case the old estimator also guarded: the between-stay component
  # comes out non-positive, which would drive every shrunk value to zero. It
  # means the replicate variance already exceeds the residual variance, so the
  # decomposition does not hold on these rows and the reduced model is the
  # honest fallback.
  if (!is.finite(s_u) || s_u <= 0) return(out(0, vr, TRUE, FALSE, "s_u_nonpositive", df))

  out(s_e, s_u, TRUE, TRUE, "replicated", df)
}

#' Pick the variance-component estimator: replicates when we have them.
#'
#' Backward compatible on purpose. A feature table extracted before 2026-09-02
#' carries no `se_` columns, and such a run must still reproduce rather than
#' fail -- `tests/refactor_identity.R` depends on exactly that.
.var_components <- function(r_val, n, within_ss = NULL, within_n = NULL) {
  if (is.null(within_ss) || is.null(within_n) ||
      !length(within_ss) || all(is.na(within_ss))) {
    v <- .fit_var_components(r_val, n)
    v$se_source <- "heteroscedastic_legacy"
    v$se_df <- 0
    return(v)
  }
  .fit_var_components_rep(r_val, n, within_ss, within_n)
}

#' The value-scale residual the variance decomposition is built on.
#'
#' ONE IMPLEMENTATION, BOTH FORMS, and that is the point. The weight `delta`
#' needs is the scalar `s_u / (s_u + s_e / n)`, and a RATIO of two variances
#' estimated on the SAME scale is scale-free. So the decomposition is done on
#' the value scale -- where the replicates live and where `se_within_ss` is
#' measured -- and the resulting weight is then applied to whichever residual
#' the form produces. That decouples variance estimation from residual
#' construction, which is what lets the ordinal rows use replicate statistics
#' at all: their own residual is a mid-PIT logit and is not commensurable with
#' a sum of squares on the raw Glasgow scale.
#'
#' The regression is the same three-term design both forms use, so for the
#' linear form this reproduces its own residual exactly and nothing changes.
.value_scale_resid <- function(v, k, n) {
  lk <- log1p(k); ln <- log(n); pos <- as.numeric(k > 0)
  d <- data.frame(v = v, lk = lk)
  fm <- "v ~ lk"
  if (stats::sd(ln)  > 0) { d$ln  <- ln;  fm <- paste(fm, "+ ln") }
  if (stats::sd(pos) > 0) { d$pos <- pos; fm <- paste(fm, "+ pos") }
  g <- stats::lm(stats::as.formula(fm), data = d)
  as.numeric(stats::residuals(g))
}

# --- delta: magnitude given deviation count ---------------------------------
#
# For signal g, magnitude variable v (one of q05/q95/value_min/value_max), and
# the deviation count k on v's own side:
#
#   E[v | k, n] = a0 + a1 * log(1 + k) + a2 * log(n)          fitted, 3 params
#   r_i         = v_i - E[v | k_i, n_i]                        residual
#   delta_i     = r_i * s_u / (s_u + s_e / n_i)                shrunk
#
# THE PRECISION IS n_obs, NOT k. These are two different roles and conflating
# them breaks the estimator. `k` is the CONDITIONING variable -- it is what the
# mean model standardises against, and it is the whole point of the construct.
# `n_obs` is the PRECISION -- how well determined this stay's magnitude is. The
# sampling noise in an order statistic depends on how many observations were
# drawn, not on how many of them crossed the bound: a stay with n_obs = 1 has a
# very noisy "minimum", a stay with n_obs = 24 has a well determined one. That
# is exactly the role n_obs plays in pi_hat, and using it here keeps the two
# constructs parallel.
#
# It also matters for identifiability. MEASURED 2026-08-28: with 1 + k as the
# weight, lactate's low tail is occupied on 53 of 21,860 measured stays, so the
# weight is 1 almost everywhere, the variance split is estimated from those 53
# rows, and the fit collapsed to s_u ~ 0 -- giving a `value_min_delta` column
# with sd 2.8e-06 that check_model_frame accepts as non-constant and the smooth
# then penalises to nothing. n_obs varies across every stay, so the split is
# identified wherever the signal was measured at all.
#
# `log(1 + k)` rather than `log(k)` because k = 0 IS INFORMATIVE HERE and must
# stay in the fit. This is the one place the measurement construct differs
# structurally from the intervention one, and getting it backwards would be a
# serious bug: for an intervention, D = 0 means nothing happened and the row
# carries no intensity information, so it must be excluded. For a signal, k = 0
# means "measured, never crossed the bound" -- a real observation with a real
# magnitude, whose expectation still rises with n because more draws reach
# further even inside the normal range. Excluding k = 0 rows would delete the
# majority of stays for a sparse lab and bias the mean model badly.
#
# The zero-information set for a signal is n_obs = 0, and those rows are already
# excluded upstream, in the one place the scoping is written.
#
# SIGN IS NOT FLIPPED. For a low-side magnitude, lower is worse, so delta < 0 is
# the sicker patient; for a high-side magnitude it is delta > 0. That matches
# the raw term this replaces, so a fitted smooth reads the same way round.

# --- delta, ordinal form ----------------------------------------------------
#
# WHY A SECOND FORM EXISTS AT ALL. The linear form above models
#
#   E[v | k, n] = a0 + a1 log(1+k) + a2 log(n)
#
# by OLS and subtracts. That is right for a continuous magnitude and wrong for a
# bounded ordinal one, and MEASURED 2026-09-01 on MIMIC-IV train it is wrong in
# two specific ways for the three Glasgow components.
#
#   THE SUPPORT IS VIOLATED. gcs_motor is observed on [1, 6] and the fitted
#   conditional mean runs from -1.72 to 6.05, off-scale for 7.6% of stays
#   (gcs_eyes 3.9%, gcs_verbal 2.1%). A residual taken against a negative
#   Glasgow score is a deviation from nothing.
#
#   THE RELATIONSHIP IS A STEP, NOT A SLOPE. k_low is 0 on 86.5% of measured
#   stays for gcs_motor. Among those, the fraction at the scale floor is 0.0%;
#   among stays with k_low > 0 it is 60.7% (gcs_eyes 0.0% vs 73.2%, gcs_verbal
#   0.0% vs 65.6%). So v | k is close to a two-point mixture whose MIXING
#   WEIGHT moves with k, and a line through it leaves the step behind. That
#   leftover is exactly what the form diagnostic sees: resid_smooth_r2 of 0.527,
#   0.485 and 0.399 for verbal, motor and eyes, against <= 0.082 for every other
#   (signal, variable) pair in the design.
#
# WHAT THE ORDINAL FORM DOES INSTEAD. A proportional-odds cumulative-link model
#
#   P(V <= c | k, n) = plogis(theta_c - eta),   eta = a1 log(1+k) + a2 log(n)
#
# puts the k dependence on a LATENT scale and lets the cut points place the
# categories. The step is then representable without anyone writing an
# indicator, because a jump in the DISTRIBUTION is what cut points are for. The
# support is correct by construction. There is no intercept: a0 is absorbed into
# theta, and is stored as NA so nothing downstream can read a meaningless one.
#
# THE RESIDUAL IS A MID-PIT, AND IT IS DETERMINISTIC ON PURPOSE.
#
#   u_i       = F(c_i - 1 | k_i, n_i) + 0.5 f(c_i | k_i, n_i)
#   delta_raw = qlogis(u_i)
#
# The textbook residual for a discrete response is the RANDOMISED quantile
# residual, which draws a uniform inside the probability interval. That is
# unusable here: hard rule 8 requires an apply site to recompute delta as a pure
# function of that row's counts and the frozen parameters, and an RNG draw at
# eICU would make the transported covariate irreproducible. The mid-point of the
# interval is the standard deterministic substitute, it is monotone in the
# observed category at fixed (k, n), and that monotonicity is what preserves the
# sign convention: delta > 0 is a HIGHER Glasgow score than the count and
# coverage predict, i.e. the less sick patient, exactly as the linear form's
# low-side residual already reads.
#
# THE PRECISION LAYER IS UNCHANGED. `.fit_var_components(delta_raw, n)` and
# the transformed residual is returned exactly as the linear one is (the
# linear one, so the AIC guard, the boundary guard and the `shrinks` fallback
# all behave identically. n_obs remains the precision for the same reason as
# before: a stay with n_obs = 1 has a very noisy minimum.
#
# THE SCALE IS DECLARED, NEVER DETECTED (hard rule 8). Glasgow motor is 1-6 by
# definition, not by observation. Detecting the level set would let it differ
# between sites, and it demonstrably would: eICU's gcs_verbal sits at its floor
# for 56.3% of ventilated measurements against MIMIC's 13.0% (CLAUDE.md, frozen
# decisions). A site whose data happened to miss a category would silently get a
# different model rather than a reportable difference.

#' Cut points from the unconstrained parameterisation, and back.
#'
#' theta must be strictly increasing or the model is not a distribution. It is
#' parameterised as theta_1 plus a cumulative sum of positive increments, which
#' makes that automatic and lets an unconstrained optimiser run.
#'
#' Each increment is floored at ORDINAL_EPS rather than allowed to reach zero.
#' Where a category is unobserved at the training site the likelihood pushes its
#' increment toward zero, two cut points collide, and the category gets
#' probability exactly 0 -- harmless at the fitting site, where no row is in it,
#' and NOT harmless at an apply site, where a patient in that category would get
#' an undefined mid-PIT. The floor keeps every category reachable.
ORDINAL_EPS <- 1e-4

.theta_from_par <- function(p, n_cut) {
  if (n_cut == 1L) return(p[1])
  c(p[1], p[1] + cumsum(ORDINAL_EPS + exp(p[2:n_cut])))
}

#' Category probabilities under the proportional-odds model.
#'
#' @return matrix, rows = observations, cols = categories 1..C
.polr_probs <- function(theta, eta) {
  C <- length(theta) + 1L
  cum <- vapply(theta, function(th) stats::plogis(th - eta), numeric(length(eta)))
  if (!is.matrix(cum)) cum <- matrix(cum, nrow = length(eta))
  cbind(cum, 1) - cbind(0, cum)
}

#' Cumulative probability strictly BELOW each category, i.e. F(c - 1).
.polr_cum_below <- function(theta, eta) {
  cum <- vapply(theta, function(th) stats::plogis(th - eta), numeric(length(eta)))
  if (!is.matrix(cum)) cum <- matrix(cum, nrow = length(eta))
  cbind(0, cum)
}

#' Fit the proportional-odds model by direct maximum likelihood.
#'
#' Hand-rolled with optim rather than delegated to MASS::polr, for the same
#' reason fit_dm_alpha() and fit_lambda_binom() are: the exact parameterisation
#' is what goes into the bundle and gets evaluated at eICU, so it must be ours
#' and it must be stable under a category that is thin at one site.
#'
#' @param yc integer category index in 1..C
#' @param X  design matrix, no intercept (it is absorbed into the cut points)
#' @return list(theta, beta, converged, loglik)
.polr_fit <- function(yc, X, C, maxit = 5000L) {
  n_cut <- C - 1L
  # Start from the empirical cumulative logits with zero slopes, which is the
  # exact MLE when the covariates carry nothing and a good start when they do.
  tab <- tabulate(yc, nbins = C)
  cp  <- cumsum(pmax(tab, 0.5)) / (sum(tab) + 1)
  th0 <- stats::qlogis(pmin(pmax(cp[seq_len(n_cut)], 1e-4), 1 - 1e-4))
  inc <- pmax(diff(th0), ORDINAL_EPS * 2)
  p0  <- c(th0[1], if (n_cut > 1L) log(pmax(inc - ORDINAL_EPS, 1e-8)) else numeric(0),
           rep(0, ncol(X)))

  nll <- function(p) {
    theta <- .theta_from_par(p, n_cut)
    beta  <- p[(n_cut + 1L):length(p)]
    eta   <- as.numeric(X %*% beta)
    pr    <- .polr_probs(theta, eta)
    px    <- pr[cbind(seq_along(yc), yc)]
    if (any(!is.finite(px)) || any(px <= 0)) return(1e12)
    -sum(log(px))
  }

  o <- try(stats::optim(p0, nll, method = "BFGS",
                        control = list(maxit = maxit, reltol = 1e-10)),
           silent = TRUE)
  if (inherits(o, "try-error") || !is.finite(o$value)) {
    o <- try(stats::optim(p0, nll, method = "Nelder-Mead",
                          control = list(maxit = maxit * 2L, reltol = 1e-10)),
             silent = TRUE)
  }
  if (inherits(o, "try-error") || !is.finite(o$value)) {
    return(list(theta = .theta_from_par(p0, n_cut),
                beta = rep(0, ncol(X)), converged = FALSE, loglik = NA_real_))
  }
  list(theta = .theta_from_par(o$par, n_cut),
       beta = o$par[(n_cut + 1L):length(o$par)],
       converged = o$convergence == 0L, loglik = -o$value)
}

#' The mid-probability integral transform, on the logit scale.
#'
#' Clamped away from 0 and 1 so a category the model gives near-zero mass at
#' still yields a finite residual instead of an infinite one. The clamp binds
#' only where the model already says the observation is essentially impossible,
#' which at an apply site is itself the finding.
.mid_pit_logit <- function(theta, eta, yc) {
  below <- .polr_cum_below(theta, eta)[cbind(seq_along(yc), yc)]
  pmf   <- .polr_probs(theta, eta)[cbind(seq_along(yc), yc)]
  u     <- below + 0.5 * pmf
  stats::qlogis(pmin(pmax(u, 1e-6), 1 - 1e-6))
}

#' Map declared scale bounds to the integer category vector.
#'
#' Errors rather than truncating: a value outside the DECLARED scale means the
#' declaration and the extraction disagree, and silently clamping it would turn
#' a data problem into a slightly-wrong covariate.
.ordinal_levels <- function(scale) {
  if (is.null(scale$min) || is.null(scale$max)) {
    stop("ordinal scale needs both `min` and `max`", call. = FALSE)
  }
  seq(as.integer(scale$min), as.integer(scale$max))
}

#' Fit the ordinal delta parameters for one (signal, magnitude variable).
#'
#' THE LATENT PREDICTOR HAS TWO REGIMES, AND THAT IS MEASURED RATHER THAN
#' ASSUMED:
#'
#'   eta = a1 log(1+k) + a2 log(n) + a3 I(k > 0)
#'
#' The pre-specified check was whether proportional odds holds. It was run on
#' 2026-09-01 and the answer is that PROPORTIONAL ODDS WAS NEVER THE PROBLEM.
#' Relaxing it -- one log(1+k) slope per cut point, `ppo` below -- leaves
#' resid_smooth_r2 at 0.165 / 0.063 / 0.118 for motor / eyes / verbal and loses
#' on AIC to the far cheaper alternative. Adding the single regime term instead
#' takes the same statistic to 0.0044 / 0.0029 / 0.0032 and wins AIC by
#' 10,300 / 6,814 / 7,770 points for ONE extra parameter:
#'
#'   signal      PO r2   +step r2   PPO r2    AIC PO / step / ppo
#'   gcs_motor   0.2955   0.0044    0.1647    39755 / 29455 / 34444
#'   gcs_eyes    0.0999   0.0029    0.0629    53960 / 47146 / 50598
#'   gcs_verbal  0.1976   0.0032    0.1182    45973 / 38203 / 42567
#'
#' WHY IT WORKS, and it is not the indicator doing the work by itself. Splitting
#' the PO residual by regime shows every remaining structure inside the k > 0
#' rows (r2 0.537 / 0.242 / 0.477 there, undefined at k = 0 where log(1+k) is
#' constant). A single slope forced through both regimes has to pass through the
#' k = 0 point, which pins it and leaves it wrong everywhere k > 0. The
#' indicator frees the intercept so the slope can serve the k > 0 regime alone.
#'
#' THIS IS NOT THE `ever_active` INDICATOR THAT WAS REJECTED. That one sat on a
#' model PREDICTOR, duplicating a point mass a smooth in the same variable could
#' already represent, and measured out at <= 0.0019 deviance. This one sits on
#' the CONDITIONING variable inside an outcome-blind standardisation whose
#' conditional mean is deliberately rigid -- no smooth, so it extrapolates at
#' eICU rather than bending to MIMIC -- and it removes 0.29 of leftover residual
#' structure. The boundary is also semantic rather than fitted: k = 0 is
#' "measured, never crossed the reference bound" and k > 0 is "crossed it", and
#' those are different states, not two points on one scale.
#'
#' @param scale list(min, max), declared in config. Never derived from `v`.
#' @return the same shape fit_delta() returns, plus `form`, `theta` and `a3`
fit_delta_ordinal <- function(v, k, n, scale,
                              within_ss = NULL, within_n = NULL) {
  stopifnot(length(v) == length(k), length(k) == length(n))
  if (anyNA(v) || anyNA(k) || anyNA(n)) stop("fit_delta_ordinal: NA in input", call. = FALSE)
  if (any(n <= 0)) stop("fit_delta_ordinal: n_obs <= 0 reached the magnitude fit", call. = FALSE)

  lv <- .ordinal_levels(scale)
  yc <- match(v, lv)
  if (anyNA(yc)) {
    abort_values(sprintf(paste0("fit_delta_ordinal: %d value(s) fall outside the ",
                                "DECLARED scale [%d, %d]. The declaration and the ",
                                "extraction disagree; fix one, do not clamp here"),
                         sum(is.na(yc)), min(lv), max(lv)),
                 sort(unique(v[is.na(yc)])))
  }
  if (length(unique(yc)) < 2L) {
    return(list(form = "ordinal", theta = numeric(0), a0 = NA_real_,
                a1 = 0, a2 = 0, a3 = 0, s_e = 0, s_u = 0, shrinks = FALSE,
                se_source = "degenerate_scale", se_df = 0,
                n_stays = length(v), converged = FALSE, degenerate = TRUE,
                levels = lv, loglik = NA_real_))
  }

  lk <- log1p(k); ln <- log(n); pos <- as.numeric(k > 0)
  # A constant column is dropped rather than handed to the optimiser: log(n) is
  # constant where every measured stay has identical coverage, and I(k > 0) is
  # constant where a signal's modelled tail is always or never occupied.
  use_n   <- stats::sd(ln)  > 0
  use_pos <- stats::sd(pos) > 0
  X <- cbind(lk = lk,
             if (use_n)   cbind(ln = ln)   else NULL,
             if (use_pos) cbind(pos = pos) else NULL)

  f <- .polr_fit(yc, X, C = length(lv))
  nm <- colnames(X)
  gv <- function(w) { i <- match(w, nm); if (is.na(i)) 0 else unname(f$beta[i]) }
  a1 <- gv("lk"); a2 <- gv("ln"); a3 <- gv("pos")

  eta <- as.numeric(X %*% f$beta)
  r   <- .mid_pit_logit(f$theta, eta, yc)
  # The variance decomposition runs on the VALUE scale, not on `r`. See
  # `.value_scale_resid()`: the weight is a scale-free ratio, and the replicate
  # statistics are measured on the raw ordinal scale, so decomposing there and
  # applying the weight to the mid-PIT residual is the only way these rows get
  # a replicate-based `s_e` at all.
  vc  <- .var_components(.value_scale_resid(v, k, n), n, within_ss, within_n)

  list(form = "ordinal", theta = f$theta, a0 = NA_real_,
       a1 = a1, a2 = a2, a3 = a3,
       s_e = vc$s_e, s_u = vc$s_u, shrinks = vc$shrinks,
       se_source = vc$se_source, se_df = vc$se_df,
       n_stays = length(v), converged = f$converged && vc$converged,
       degenerate = FALSE, levels = lv, loglik = f$loglik)
}

#' Apply frozen ordinal delta parameters. Recomputed per row, never re-fitted.
#'
#' RENAMED FROM `shrink_delta_ordinal()` 2026-09-05: it no longer shrinks, so a
#' `shrink_` name would have been a second declaration of something the code
#' does not do. The returned value is the mid-PIT residual against the frozen
#' cumulative link, and nothing further is applied to it.
#'
#' A value outside the frozen scale is an error here as it is at fit time. At an
#' apply site that means the two extractions disagree about what the scale IS,
#' which is a transport finding and must not be absorbed into a clamped value.
delta_value_ordinal <- function(v, k, n, par) {
  if (isTRUE(par$degenerate)) return(rep(0, length(v)))
  lv <- .par_levels(par)
  th <- .par_theta(par)
  yc <- match(v, lv)
  if (anyNA(yc)) {
    abort_values(sprintf(paste0("delta_value_ordinal: %d value(s) outside the ",
                                "FROZEN scale [%d, %d]. The site being scored ",
                                "disagrees with the training site about the scale"),
                         sum(is.na(yc)), min(lv), max(lv)),
                 sort(unique(v[is.na(yc)])))
  }
  eta <- par$a1 * log1p(k) + par$a2 * log(n) + (par$a3 %||% 0) * as.numeric(k > 0)
  .mid_pit_logit(th, eta, yc)
}

# --- serialising the variable-length pieces ---------------------------------
#
# The magnitude priors table is FLAT: one row per (signal, variable, role,
# fold), with scalar columns. theta and levels are variable-length, so they are
# stored as delimited strings rather than as list columns. A list column would
# survive the .rds path and break the .csv one, and both are written for every
# run -- so the CSV would silently become the odd one out. %.17g round-trips a
# double exactly, so nothing is lost to the encoding.

.pack_num <- function(x) if (!length(x)) "" else paste(sprintf("%.17g", x), collapse = ";")

.unpack_num <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(numeric(0))
  as.numeric(strsplit(s, ";", fixed = TRUE)[[1]])
}

.par_theta  <- function(par) if (is.character(par$theta)) .unpack_num(par$theta) else par$theta
.par_levels <- function(par) if (is.character(par$levels)) .unpack_num(par$levels) else par$levels

#' Fit the delta parameters for one (signal, magnitude variable).
#'
#' A ROUTER over the two forms. `form` is DECLARED per signal in config and
#' resolved by `delta_form_of()`; it is never inferred from the data, because a
#' form chosen by inspecting the response would differ between sites and turn a
#' transport result into a plumbing artifact (hard rule 8).
#'
#' Both forms return the same shape, so `magnitude_priors()` stores one flat
#' table and `delta_value()` routes on the stored `form` alone. That is why the
#' `.frame_magnitude()` call site needs no edit for either form: everything it
#' needs travels in `par`.
#'
#' @param v magnitude values, one per stay
#' @param k deviation count on v's side (k_low for a low-side variable)
#' @param n n_obs, covered hours
#' @param form  "linear" or "ordinal"
#' @param scale list(min, max) for the ordinal form; ignored by the linear one
#' @return list(form, a0, a1, a2, s_e, s_u, shrinks, theta, levels, n_stays,
#'   converged, degenerate)
fit_delta <- function(v, k, n, form = "linear", scale = NULL,
                      within_ss = NULL, within_n = NULL) {
  if (identical(form, "ordinal")) {
    if (is.null(scale)) {
      stop("fit_delta: the ordinal form needs a DECLARED scale. Add the signal ",
           "to `magnitude_ordinal_scale` in config -- the level set is never ",
           "detected from the data (hard rule 8).", call. = FALSE)
    }
    return(fit_delta_ordinal(v, k, n, scale, within_ss, within_n))
  }
  if (!identical(form, "linear")) {
    stop("fit_delta: unknown magnitude form '", form, "' (want linear|ordinal)",
         call. = FALSE)
  }

  stopifnot(length(v) == length(k), length(k) == length(n))
  if (anyNA(v) || anyNA(k) || anyNA(n)) stop("fit_delta: NA in input", call. = FALSE)
  if (any(n <= 0)) stop("fit_delta: n_obs <= 0 reached the magnitude fit; the ",
                        "measured-subset scoping upstream is wrong", call. = FALSE)

  # A magnitude that never varies carries nothing to standardise. Reported
  # rather than silently producing a constant-0 delta and a null smooth.
  if (stats::sd(v) == 0) {
    return(list(form = "linear", theta = numeric(0), levels = numeric(0),
                a0 = mean(v), a1 = 0, a2 = 0, a3 = 0, s_e = 0, s_u = 0, shrinks = FALSE,
                se_source = "degenerate_magnitude", se_df = 0,
                n_stays = length(v), converged = FALSE, degenerate = TRUE,
                loglik = NA_real_))
  }

  # THE k = 0 REGIME TERM. `log(1 + k)` is smooth at 0 and the model has ONE
  # intercept, so the fitted value at k = 0 is pinned to `a0 + a2 log(n)`. Since
  # 22% to 94% of measured stays sit at exactly k = 0 depending on the signal,
  # least squares drives that intercept to the k = 0 mean and the k > 0 arm is
  # then served only by `a1` rotating around a pinned point. `I(k > 0)` gives the
  # k > 0 arm its own intercept -- and because `log(1 + 0) = 0` exactly, the
  # interaction `log(1+k) : I(k>0)` IS `log(1+k)`, so the same one parameter also
  # frees the slope. Adding a step term and fitting a two-part model are the same
  # model here, and MEASURED 2026-09-02 they return identical residuals.
  #
  # THIS IS A DISCONTINUITY, NOT CURVATURE, so no reshaping of a smooth function
  # of k could absorb it. `k = 0` means "measured, never crossed the reference
  # bound" and `k >= 1` means "crossed it"; for a signal parameterised on
  # `value_min` / `value_max` the implication is exact, because a count of hours
  # beyond a bound is positive if and only if the extremum is beyond it.
  # MEASURED (tests/delta_step_evidence.R): the empirical k = 0 -> k = 1 step
  # exceeds the largest step the smooth form can express, `a1 log(2)`, on 23 of
  # 24 (signal, variable) rows, and the median |slope change| once the k = 0
  # anchor is released is 0.52 of the original slope.
  #
  # NO SEPARATION RISK HERE. `I(k > 0)` is a deterministic recode of the response
  # for an extreme-parameterised signal, which is fatal for a DISCRETE likelihood
  # and harmless for this Gaussian one: the response stays continuous inside each
  # arm and the indicator only shifts a conditional mean. The ordinal form above
  # carries the same term for the same reason and is the case where the
  # near-determinism has to be watched.
  lk <- log1p(k); ln <- log(n); pos <- as.numeric(k > 0)
  # A constant column is dropped rather than handed to lm(): log(n) is constant
  # wherever every measured stay has identical coverage, and I(k > 0) is constant
  # wherever a signal's modelled tail is always or never occupied.
  use_n   <- stats::sd(ln)  > 0
  use_pos <- stats::sd(pos) > 0
  d <- data.frame(v = v, lk = lk)
  fm <- "v ~ lk"
  if (use_n)   { d$ln  <- ln;  fm <- paste(fm, "+ ln") }
  if (use_pos) { d$pos <- pos; fm <- paste(fm, "+ pos") }
  g  <- stats::lm(stats::as.formula(fm), data = d)
  cf <- stats::coef(g)
  gv <- function(nm) { z <- unname(cf[nm]); if (length(z) != 1L || is.na(z)) 0 else z }
  a0 <- gv("(Intercept)"); a1 <- gv("lk"); a2 <- gv("ln"); a3 <- gv("pos")

  r  <- v - (a0 + a1 * lk + a2 * ln + a3 * pos)
  # `r` IS the value-scale residual for this form, so no auxiliary regression is
  # needed and the linear rows decompose on exactly the quantity they shrink.
  vc <- .var_components(r, n, within_ss, within_n)

  list(form = "linear", theta = numeric(0), levels = numeric(0),
       a0 = a0, a1 = a1, a2 = a2, a3 = a3, s_e = vc$s_e, s_u = vc$s_u,
       shrinks = vc$shrinks, se_source = vc$se_source, se_df = vc$se_df,
       n_stays = length(v), converged = vc$converged, degenerate = FALSE,
       loglik = NA_real_)
}

#' The UN-SHRUNK residual a stored delta parameter set implies, either form.
#'
#' THE SINGLE IMPLEMENTATION. `delta_value()` and every diagnostic call this,
#' so a diagnostic cannot measure a different residual from the one the model
#' frame carries. That guarantee was previously a comment rather than a
#' structure, and on 2026-09-02 it failed exactly as such comments do: `a3` was
#' added to the fit and the apply path but not to the diagnostic's private copy,
#' which then reported a step of exactly `a3` between the two regimes as if it
#' were leftover structure. The frames were correct throughout; only the
#' diagnostic was wrong. Sharing the function closes the class.
.delta_raw_resid <- function(par, v, k, n) {
  if (identical(par$form %||% "linear", "ordinal")) {
    lv <- .par_levels(par); th <- .par_theta(par)
    yc <- match(v, lv)
    if (anyNA(yc)) return(NULL)
    eta <- par$a1 * log1p(k) + par$a2 * log(n) + (par$a3 %||% 0) * as.numeric(k > 0)
    return(.mid_pit_logit(th, eta, yc))
  }
  v - (par$a0 + par$a1 * log1p(k) + par$a2 * log(n) +
         (par$a3 %||% 0) * as.numeric(k > 0))
}

#' Apply frozen delta parameters to rows. Recomputed per row, never re-fitted.
#'
#' RENAMED FROM `shrink_delta()` AND THE SHRINKAGE REMOVED, 2026-09-05. The
#' covariate is now exactly the residual against the frozen conditional mean:
#'
#'     delta_i = v_i - (a0 + a1*log(1+k_i) + a2*log(n_i))
#'
#' In words, and this is the whole of what the construct claims: how much more
#' extreme this stay's magnitude was than a stay with this deviation count and
#' this coverage typically shows. One sentence, in the original measurement
#' units, with no second dependence on coverage.
#'
#' WHAT WAS REMOVED AND WHY, because a future reader will want to add it back.
#' A precision weight `s_u / (s_u + s_e/n)` used to multiply the residual, on
#' the model that part of it is measurement noise that averages out with more
#' observation. Three measurements killed it, and the first is decisive:
#'
#'   1. It imported `n_obs`, the one variable the design excludes. Its deviance
#'      gain over the unweighted residual was +3.379 points summed over the 38
#'      delta rows, and -1.591 once `s(log n_obs)` entered both models. The gain
#'      did not survive conditioning on coverage; it REVERSED.
#'   2. Its variance law is inverted where it matters. It needs residual spread
#'      to FALL as n rises; spread RISES in 17 of 22 rows whose magnitude
#'      variable is an extreme, because `value_min`/`value_max` are MIN/MAX over
#'      raw readings and more hours means more chances at a rarer value.
#'   3. It does not transport. `s_e` and `s_u` freeze at MIMIC while `n` comes
#'      from the scored stay, so `lactate/value_min` is compressed 41% at eICU
#'      for no reason but that eICU draws one lactate where MIMIC draws two.
#'
#' What is given up is the errors-in-variables attenuation correction. That is
#' the right trade: attenuation is CONSERVATIVE, while coverage leakage
#' manufactures signal that will not transport.
#'
#' Full argument in docs/covariate_constructs_20260905.md.
#'
#' Routes on the STORED form, so a bundle fitted under one form cannot be
#' applied under another. A `par` with no `form` is treated as linear, which is
#' what every bundle written before 2026-09-01 carries.
delta_value <- function(v, k, n, par) {
  if (isTRUE(par$degenerate)) return(rep(0, length(v)))
  r <- .delta_raw_resid(par, v, k, n)
  if (is.null(r)) {
    abort_values(sprintf(paste0("delta_value [%s]: value(s) outside the FROZEN ",
                                "scale. The site being scored disagrees with the ",
                                "training site about the level set"),
                         par$form %||% "linear"),
                 sort(unique(v[is.na(match(v, .par_levels(par)))])))
  }
  r
}

# --- lambda: intensity given exposure ---------------------------------------
#
# Two families, one output scale. Both return 0 for an unexposed stay and 0 for
# a stay of exactly typical intensity, with sign meaning "more / less intense
# than the population gives a stay treated this long".
#
#   binomial    A is a bounded count of distinct molecules out of a fixed pool
#               of M, observed EXACTLY. The output is the Pearson residual
#               against the conditional mean, (A - M p(D)) / sd(A | D).
#   lognormal   A is a continuous amount accumulated over D hours, so it carries
#               genuine within-stay noise. The output is the shrunk log residual.
#
# WHY THE COUNT FAMILY DOES NOT SHRINK, and why it is not the beta-binomial.
#
# The beta-binomial was the natural choice -- it is literally the
# Dirichlet-multinomial with J = 2 categories, with the prior mean p(D) in place
# of a constant -- and the data rejected it. MEASURED 2026-08-28 on the MIMIC-IV
# training set, variance of n_agents among exposed stays relative to
# Binomial(M, p):
#
#   vasopressor   var 0.487 vs binomial 1.001   ratio 0.487
#   inotrope      var 0.026 vs binomial 0.500   ratio 0.052
#
# Both are UNDER-dispersed, and a beta-binomial can only represent OVER-
# dispersion. Its ML therefore diverges: phi ran to the optimiser's upper bound
# (162,744 ~ e^12) and the posterior mean collapsed onto the prior mean, giving
# a lambda column with sd 1.1e-05 over 99 distinct values. That is a silently
# null covariate, not a small one.
#
# The under-dispersion is structural rather than incidental: n_agents is not a
# sum of M independent Bernoullis, it is an escalation ladder that almost always
# starts at one molecule and rarely goes past two. So there is no extra-binomial
# heterogeneity to shrink toward -- and, more fundamentally, no within-stay
# sampling noise at all. "How many molecules did this patient receive" is
# observed exactly; it is not an estimate of anything.
#
# phi -> infinity IS the binomial limit, and the informative quantity there is
# the standardised residual. So the two families are the same estimator at
# different precisions: standardise against the conditional mean, then shrink by
# the observation precision. For a continuous accumulation that precision is
# finite and rises with D; for an exactly-observed count it is infinite and the
# shrinkage weight is 1. Nothing is special-cased -- one of the two cases simply
# has no noise to average away.
#
# If a site ever shows OVER-dispersed molecule counts, the beta-binomial is the
# generalisation to reach for. Declared here from MIMIC, not detected per site
# (hard rule 8): a family that switched on the observed dispersion would be a
# branch on data, and check_agent_pool() is where a site disagreement surfaces.

#' Fit the binomial intensity model for one intervention's molecule count.
#'
#' @param A distinct-molecule count, 0..M
#' @param D exposure in hours
#' @param M the molecule pool size, DECLARED in config, never detected here
fit_lambda_binom <- function(A, D, M) {
  ex <- D > 0
  if (!any(ex)) stop("fit_lambda_binom: no exposed stays", call. = FALSE)
  if (any(A[ex] > M)) {
    stop(sprintf(paste0("fit_lambda_binom: %d stay(s) exceed the declared molecule ",
                        "pool of %d. Update intervention_agent_pool in config, or the ",
                        "extraction is counting drug names rather than molecules ",
                        "(canonical_variable_spec.md SS7)."), sum(A[ex] > M), M),
         call. = FALSE)
  }
  a <- A[ex]; d <- D[ex]
  if (stats::sd(a) == 0) {
    return(list(family = "binomial", M = M, c0 = stats::qlogis(mean(a) / M),
                c1 = 0, n_stays = length(a), converged = FALSE, degenerate = TRUE))
  }

  g  <- stats::glm(cbind(a, M - a) ~ log(d), family = stats::binomial())
  c0 <- unname(stats::coef(g)[1]); c1 <- unname(stats::coef(g)[2])
  if (is.na(c1)) c1 <- 0

  list(family = "binomial", M = M, c0 = c0, c1 = c1,
       n_stays = length(a), converged = isTRUE(g$converged), degenerate = FALSE)
}

#' Fit the lognormal intensity model for one intervention.
fit_lambda_ln <- function(A, D) {
  ex <- D > 0 & A > 0
  if (!any(ex)) stop("fit_lambda_ln: no exposed stays with a positive amount", call. = FALSE)
  a <- log(A[ex]); d <- D[ex]
  if (stats::sd(a) == 0) {
    return(list(family = "lognormal", b0 = mean(a), b1 = 0, s_e = 0, s_u = 0,
                shrinks = FALSE, n_stays = length(a), converged = FALSE,
                degenerate = TRUE))
  }
  use_d <- stats::sd(log(d)) > 0
  g  <- if (use_d) stats::lm(a ~ log(d)) else stats::lm(a ~ 1)
  cf <- stats::coef(g)
  b0 <- unname(cf[1]); b1 <- unname(if (use_d && !is.na(cf[2])) cf[2] else 0)
  r  <- a - (b0 + b1 * log(d))
  vc <- .fit_var_components(r, d)
  list(family = "lognormal", b0 = b0, b1 = b1, s_e = vc$s_e, s_u = vc$s_u,
       shrinks = vc$shrinks,
       n_stays = length(a), converged = vc$converged, degenerate = FALSE)
}

#' Apply frozen lambda parameters to rows. Unexposed stays get 0 by definition.
#'
#' Assigned, never predicted -- the same treatment unmeasured stays get for L
#' (spec SS5.5) and unexposed stays get from pi_hat at n_obs = 0. An unexposed
#' stay has no intensity to be relative to, and 0 is the neutral value on this
#' scale, so it lands there by construction rather than by extrapolating a
#' conditional mean to D = 0 where log(D) is undefined.
lambda_value <- function(A, D, par) {
  out <- numeric(length(D))
  if (isTRUE(par$degenerate)) return(out)

  if (identical(par$family, "binomial")) {
    ex <- D > 0
    if (!any(ex)) return(out)
    pe <- stats::plogis(par$c0 + par$c1 * log(D[ex]))
    # The bounded count is clamped to the declared pool at apply time. A site
    # with a wider molecule list is a config change, not something to absorb
    # silently; but a single out-of-range row must not produce a non-finite
    # covariate mid-run, so it is clamped and check_agent_pool() is where the
    # discrepancy is reported.
    a  <- pmin(A[ex], par$M)
    # Pearson residual: observed minus expected, in units of the conditional SD.
    # sd is bounded away from 0 because pe comes from a logistic link.
    out[ex] <- (a - par$M * pe) / sqrt(par$M * pe * (1 - pe))
    return(out)
  }

  if (identical(par$family, "lognormal")) {
    ex <- D > 0 & A > 0
    if (!any(ex)) return(out)
    # THE SHRINKAGE WAS REMOVED HERE 2026-09-05, which makes the two families
    # consistent: the binomial arm above never shrank either. The covariate is
    # the log residual against the frozen conditional mean, i.e. how much more
    # amount this stay received than a stay treated for this long typically
    # receives.
    #
    # THE ARM'S OWN EVIDENCE, beyond the general argument in `delta_value()`.
    # The weight assumed residual variance falls with exposure hours. Measured
    # on MIMIC train, it RISES for two of the three lognormal interventions:
    # diuretic 0.337 -> 0.738 across D = 1 to 4, and transfusion_prbc's IQR
    # 0.069 -> 0.434 across D = 1 to 7. Diuretic's residual IQR at D = 1 is
    # exactly 0.000 -- single-dose stays get a standard dose -- and it jumps to
    # ~1.1 the moment titration starts. That divergence is treatment response,
    # the most patient-specific thing the covariate carries, and the weight
    # would have classified it as noise.
    #
    # AND THE SHRINKAGE WAS HURTING THE DECOUPLING IT WAS MEANT TO SERVE. On
    # transfusion_platelet, the one arm that shrank, r(lambda, log D) went from
    # exactly 0 to +0.00689 and the scale coupling r(|lambda|, log D) from
    # -0.028 to +0.142, because w(D) is monotone increasing and compresses
    # low-exposure stays. The decoupling lives in the conditional mean below,
    # never in the weight.
    out[ex] <- log(A[ex]) - (par$b0 + par$b1 * log(D[ex]))
    return(out)
  }

  stop("lambda_value: unknown family '", par$family, "'", call. = FALSE)
}

# --- which interventions carry a lambda -------------------------------------

#' The (exposure, accumulation, family) triple for an intervention, or NULL.
#'
#' Design information derived from the shape map and the agent-count whitelist,
#' both already in config. An intervention carries a lambda only where it has
#' TWO intensity covariates to reparameterise:
#'
#'   state + agent counts  ->  (exposure_frac, n_agents)      beta_binomial
#'   event                 ->  (n_hours, total_amount)        lognormal
#'   state, no agents      ->  exposure_frac alone. NO lambda: there is no
#'                             second coordinate, so nothing to condition on.
#'
#' Data-free, so R/05_formula.R can ask it while staying data-free itself.
lambda_spec_of <- function(intervention, cfg) {
  # REQUIRED, NOT DEFAULTED. This was `isTRUE(cfg$intensity_conditional %||%
  # FALSE)` until 2026-09-03, which meant an absent or misspelled key silently
  # switched the construct OFF while `config.yml` said `true` -- a construct
  # inert while the methods section claims it is active is exactly the failure
  # this project already had once with `delta` (audit finding F10).
  if (!cfg_flag(cfg, "intensity_conditional",
                what = "the lambda construct is on or off for the whole design")) {
    return(NULL)
  }
  shape <- (unlist(cfg$intervention_shape))[[intervention]]
  if (is.null(shape)) return(NULL)
  agents <- as.character(unlist(cfg$interventions_with_agent_counts) %||% character(0))

  if (identical(shape, "state") && intervention %in% agents) {
    pool <- cfg$intervention_agent_pool[[intervention]]
    if (is.null(pool)) {
      stop("lambda_spec_of: '", intervention, "' has agent counts but no ",
           "intervention_agent_pool entry in config. The pool size is DECLARED, ",
           "never detected (hard rule 8).", call. = FALSE)
    }
    return(list(exposure = "exposure_frac", accumulation = "n_agents",
                family = "binomial", M = as.integer(pool),
                exposure_scale = 24))
  }
  if (identical(shape, "event")) {
    return(list(exposure = "n_hours", accumulation = "total_amount",
                family = "lognormal", M = NA_integer_, exposure_scale = 1))
  }
  NULL
}


#' Which functional form the magnitude construct takes for a signal.
#'
#' Global default with per-signal overrides, exactly like `level_terms_override`
#' and `magnitude_conditional_override`. DECLARED, never detected: a form chosen
#' by inspecting the response would differ between sites, and the whole point of
#' freezing the design is that it cannot.
#'
#' Data-free, so R/05_formula.R can ask it while staying data-free itself.
delta_form_of <- function(signal, cfg) {
  o <- cfg$magnitude_form_override[[signal]]
  f <- if (!is.null(o)) as.character(o) else as.character(cfg$magnitude_form %||% "linear")
  if (!f %in% c("linear", "ordinal")) {
    stop("delta_form_of: illegal magnitude form '", f, "' for '", signal,
         "' (want linear|ordinal)", call. = FALSE)
  }
  if (identical(f, "ordinal") && is.null(cfg$magnitude_ordinal_scale[[signal]])) {
    stop("delta_form_of: '", signal, "' declares the ordinal form but has no ",
         "`magnitude_ordinal_scale` entry. The level set is DECLARED, never ",
         "detected (hard rule 8).", call. = FALSE)
  }
  f
}

#' The declared ordinal scale for a signal, or NULL under the linear form.
ordinal_scale_of <- function(signal, cfg) {
  if (!identical(delta_form_of(signal, cfg), "ordinal")) return(NULL)
  cfg$magnitude_ordinal_scale[[signal]]
}

#' Hold `magnitude_ordinal_scale` to the observed level set, in BOTH directions.
#'
#' Same pattern and same reason as `check_signal_tails()` and `check_smooth_k()`:
#' a static design declaration re-checked against the data at every site. Two
#' failures matter and they are not symmetric.
#'
#'   A value OUTSIDE the declared scale is fatal. The fit would error anyway;
#'   catching it here names the signal and the offending values.
#'   A declared level NEVER OBSERVED is reported but not fatal. It is legitimate
#'   -- gcs_verbal may simply have no stay at some level in a given cohort --
#'   and the fit handles it, because the cut-point increments are floored so the
#'   category stays reachable for a site that does observe it.
#'
#' AGGREGATES ONLY (hard rule 1): level sets and counts, never a row.
check_ordinal_scales <- function(tabs, cfg, strict = TRUE) {
  sigs <- Filter(function(sg) identical(delta_form_of(sg, cfg), "ordinal"),
                 unlist(cfg$signals))
  if (!length(sigs)) return(invisible(NULL))

  sf <- tabs$signal_features
  rows <- list(); bad <- character(0)
  for (sg in sigs) {
    lv <- .ordinal_levels(ordinal_scale_of(sg, cfg))
    z <- sf[sf$signal == sg & sf$n_obs > 0, , drop = FALSE]
    for (v in level_vars_of(sg, cfg)) {
      obs <- sort(unique(z[[v]]))
      obs <- obs[!is.na(obs)]
      outside <- setdiff(obs, lv)
      unused  <- setdiff(lv, obs)
      rows[[length(rows) + 1L]] <- data.frame(
        signal = sg, variable = v,
        declared = paste(range(lv), collapse = "-"),
        n_declared = length(lv), n_observed = length(obs),
        outside = paste(outside, collapse = ","),
        never_observed = paste(unused, collapse = ","),
        stringsAsFactors = FALSE)
      if (length(outside)) bad <- c(bad, sprintf("%s/%s: %s", sg, v,
                                                 paste(outside, collapse = ",")))
    }
  }
  out <- do.call(rbind, rows)
  if (length(bad)) {
    msg <- paste0("magnitude_ordinal_scale disagrees with the data; value(s) ",
                  "outside the declared scale: ", paste(bad, collapse = "; "),
                  ". Fix the declaration or the extraction -- do not clamp.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  unused <- out[nzchar(out$never_observed), , drop = FALSE]
  if (nrow(unused)) {
    message(sprintf("check_ordinal_scales: %d (signal, variable) pair(s) declare a level never observed here: %s",
                    nrow(unused), paste(sprintf("%s/%s [%s]", unused$signal,
                      unused$variable, unused$never_observed), collapse = "; ")))
  }
  out
}

#' Does this signal's magnitude term enter as a conditional (delta) one?
#'
#' Global flag with per-signal overrides, exactly like `level_terms_override`.
#'
#' OPEN DECISION, deliberately visible rather than silently defaulted: the three
#' GCS components are bounded integer scales of length 4-6, so `delta` there is
#' a continuous residual of a 4-to-6-valued quantity. That is statistically
#' well defined but semantically much coarser than the same construct on
#' creatinine, and it has NOT been measured. Set these to false in
#' `magnitude_conditional_override` if the measurement says it does not carry.
magnitude_conditional_for <- function(signal, cfg) {
  o <- cfg$magnitude_conditional_override[[signal]]
  if (!is.null(o)) return(isTRUE(o))
  # REQUIRED, NOT DEFAULTED -- see `lambda_spec_of()` above for the reasoning.
  # The per-signal OVERRIDE keeps its `%||%`-like optionality, because an
  # override that is absent genuinely means "no override"; it is the GLOBAL
  # switch that must not have a second declaration.
  cfg_flag(cfg, "magnitude_conditional",
           what = "the delta construct is on or off for the whole design")
}

#' The deviation-count column that pairs with a magnitude variable.
#' Low-side magnitudes standardise against k_low, high-side against k_high.
delta_count_of <- function(variable) {
  if (variable %in% c("q05", "value_min")) "k_low"
  else if (variable %in% c("q95", "value_max")) "k_high"
  else stop("delta_count_of: '", variable, "' is not a magnitude variable", call. = FALSE)
}

#' The formula-facing name of a conditional magnitude term.
delta_name_of <- function(variable) paste0(variable, "_delta")

#' The raw magnitude variable a `_delta` name derives from. Inverse of the above.
delta_base_of <- function(name) sub("_delta$", "", name)

# --- checks -----------------------------------------------------------------

#' Hold `intervention_agent_pool` to the observed molecule counts.
#'
#' Same pattern, and the same reason, as check_signal_tails() and
#' check_smooth_k(): a static design declaration re-checked against the data at
#' every site. A pool smaller than the observed maximum makes the beta-binomial
#' likelihood undefined; a pool much larger silently weakens the shrinkage,
#' because the prior mean p(D) is then a fraction of a pool that does not exist.
#'
#' The eICU direction is the one that matters. `canonical_variable_spec.md` SS7
#' records the trap: counting distinct `drugname` rather than distinct molecule
#' inflates eICU's n_agents by 2-4x. This check is what turns that from a
#' plausible-looking case-mix difference into a loud failure.
check_agent_pool <- function(tabs, cfg, strict = TRUE) {
  ivf <- tabs$intervention_features
  ivs <- names(cfg$intervention_agent_pool %||% list())
  if (!length(ivs)) return(invisible(NULL))
  rows <- lapply(ivs, function(iv) {
    z <- ivf[ivf$intervention == iv & !is.na(ivf$n_agents), , drop = FALSE]
    data.frame(intervention = iv,
               declared_pool = as.integer(cfg$intervention_agent_pool[[iv]]),
               max_observed  = if (nrow(z)) max(z$n_agents) else NA_integer_,
               n_exposed     = sum(z$n_agents > 0),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$ok <- !is.na(out$max_observed) & out$max_observed <= out$declared_pool
  if (any(!out$ok)) {
    b <- out[!out$ok, , drop = FALSE]
    msg <- paste0("intervention_agent_pool is smaller than the observed count for: ",
                  paste(sprintf("%s (declared %d, observed %d)", b$intervention,
                                b$declared_pool, b$max_observed), collapse = ", "),
                  ". Either the pool declaration is wrong, or the extraction is ",
                  "counting drug names rather than molecules.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  out
}
