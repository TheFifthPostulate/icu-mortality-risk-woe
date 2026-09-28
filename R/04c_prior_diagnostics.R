# R/04c_prior_diagnostics.R --------------------------------------------------
# Goodness-of-fit for the three FITTED, OUTCOME-BLIND parameter sets: the
# Dirichlet-multinomial alpha behind `pi_hat`, the conditional-mean parameters
# behind `delta`, and the intensity parameters behind `lambda`.
#
# WHAT MISSPECIFICATION CAN AND CANNOT DO HERE, because it decides what these
# diagnostics are for and it is the first thing a reviewer will ask.
#
#   IT CANNOT INVALIDATE ANYTHING. All three estimators are outcome-blind: none
#   of them sees `mortality`. A wrong likelihood therefore cannot leak label
#   information into a covariate and cannot manufacture apparent
#   discrimination. That is a structural guarantee, not an empirical claim.
#
#   IT CAN COST POWER. If the DM is a poor model for the counts, alpha is a
#   poorly estimated shrinkage strength and `pi_hat` over- or under-shrinks.
#   The covariate handed to `bam()` is then a noisier summary than it needs to
#   be. The XGBoost ladder already bounds this from the other direction:
#   `xgb_feat` loses only 0.0072 AUROC against `xgb_raw` while going from 221
#   columns to 94, so the construction demonstrably retains nearly all the
#   discriminative content of the raw columns. That is a SUFFICIENCY check, and
#   it is evidence the transformation is not destroying information.
#
#   IT CAN BREAK TRANSPORT, and this is the real exposure. Every parameter here
#   is fitted at MIMIC and frozen into the bundle (hard rule 8). If eICU's count
#   or exposure distribution is shaped differently, the frozen parameters
#   standardise it wrongly, and nothing downstream would reveal it.
#   `docs/v2_state_20260828.md` §6 names this risk explicitly and it is why the
#   conditional means are two- and three-parameter forms rather than splines.
#
# So: these are QUALITY and TRANSPORT diagnostics, not validity checks. Report
# them in a methods appendix and the reviewer question is pre-empted rather than
# invited.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). AGGREGATES ONLY (hard rule 1) —
# everything returned is a per-signal or per-intervention summary row.
# ----------------------------------------------------------------------------

# --- beta-binomial, the exact marginal of the DM ----------------------------

#' The marginal of a Dirichlet-multinomial IS a beta-binomial.
#'
#' If k ~ DM(n, alpha) then k_j ~ BetaBinomial(n, alpha_j, alpha_0 - alpha_j).
#' That identity is what makes these diagnostics exact rather than simulated:
#' every category's marginal distribution is available in closed form, so the
#' probability integral transform can be computed directly instead of
#' approximated from draws. With `n_obs` capped at 24 the CDF is at most 25
#' terms, so the cost is negligible.
.dbetabinom <- function(k, n, a, b) {
  exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
}

#' P(K <= q) for the beta-binomial, vectorised over equal-length `q` and `n`.
.pbetabinom <- function(q, n, a, b) {
  vapply(seq_along(q), function(i) {
    if (is.na(q[i]) || is.na(n[i])) return(NA_real_)
    if (q[i] < 0) return(0)
    kk <- 0:min(q[i], n[i])
    sum(.dbetabinom(kk, n[i], a, b))
  }, numeric(1))
}

# --- 1. is pi_hat's shrinkage doing its job? --------------------------------

#' Prior strength against observed coverage, per signal.
#'
#' THE CHEAPEST DIAGNOSTIC IN THE PROJECT AND THE ONE MOST WORTH HAVING.
#' `pi_hat` is `(k_j + alpha_j) / (n + alpha_0)`, so `alpha_0` is a prior sample
#' size measured in the same units as `n_obs` and the two are directly
#' comparable. The posterior weight on the observed proportion is
#' `n / (n + alpha_0)`, and reading it at the quartiles of `n` says immediately
#' whether the estimator is doing what it exists for:
#'
#'   weight near 0    the prior dominates and `pi_hat` is nearly constant across
#'                    stays. The covariate carries almost nothing and the smooth
#'                    will be penalised flat.
#'   weight near 1    no shrinkage. 0-of-2 and 0-of-24 are treated as the same
#'                    evidence, which is precisely the failure the whole
#'                    Dirichlet-multinomial layer exists to prevent
#'                    (CLAUDE.md, frozen decisions).
#'   weight in between  the estimator is separating confident from unconfident
#'                    zeros, which is the intent.
#'
#' No threshold is imposed on the WEIGHTS. This is a table to read, not a test
#' to pass: what counts as too much shrinkage is a judgement about the signal,
#' and imposing a cutoff here would be a tuning constant entering through the
#' back door.
#'
#' `shrinks` AND THE FLOOR FLAGS, added 2026-09-05 for parity with
#' `delta_fit_diagnostics()` and `lambda_fit_diagnostics()`, which have carried
#' a `shrinks` column since the conditional-prior constructs landed. READ THE
#' DEFINITION BEFORE USING IT, because it does NOT mean the same thing there and
#' here, and treating the three columns as one quantity would be a mistake.
#'
#'   In `delta` and `lambda`, `shrinks = FALSE` is a FALLBACK. The variance-
#'   component fit failed, or landed on the boundary, or lost to the one-
#'   component model on AIC, and the estimator then returns the raw residual
#'   with weight exactly 1. It is a statement that the machinery is INERT, and
#'   it is the flag that exposed 22 of 38 inert delta rows on 2026-09-02.
#'
#'   In the Dirichlet-multinomial there is no such fallback path. `alpha0` is
#'   strictly positive by construction -- `fit_dm_alpha()` floors every
#'   coordinate at 1e-10 -- so the posterior weight `n / (n + alpha0)` is ALWAYS
#'   below 1 and the estimator always shrinks something. `shrinks` here
#'   therefore asks the only question that can come back FALSE: did the fit
#'   converge, and is every one of the three coordinates off the estimator
#'   floor. It is DERIVED, from two things already in this table (`converged`
#'   and the floor flags below), so it is a convenience for reading this table
#'   beside its two siblings and it is NEVER independent evidence.
#'
#'   THE NUMBER THAT ANSWERS "IS THE SHRINKAGE DOING ANYTHING" IS STILL `w_med`,
#'   and now also `w_at_n1` -- the weight a stay measured exactly ONCE receives.
#'   That is the quantity the whole layer exists to get right, it is the same
#'   column the delta table reports, and it is where a prior too weak to
#'   separate 0-of-2 from 0-of-24 becomes visible.
#'
#' `floor_low`, `floor_mid` and `floor_high` say WHICH coordinate is pinned.
#' `degenerate` is their `any()`, collapsed at fit time inside
#' `training_priors()`, so until now the table could say that a signal had a
#' degenerate coordinate but not which one. The flags are recomputed here from
#' the stored alphas rather than by refitting, so they cost nothing and cannot
#' invalidate a single GAM.
#'
#' `tail_empty_low` and `tail_empty_high` are SITE-LOCAL, and they are the
#' reason this matters at eICU. The alphas are frozen at MIMIC and travel in the
#' bundle, but whether a tail actually carries counts is a property of the site
#' being scored. A coordinate pinned at the floor because MIMIC never observed
#' it, evaluated at a site that DOES observe it, is a frozen prior sitting where
#' the new site has data -- a transport warning rather than a MIMIC fitting
#' artefact. The four signals with a structurally empty high tail are the
#' declared ones (`spo2` and the three Glasgow components, CLAUDE.md's
#' `signal_tails` entry), so at MIMIC these flags must agree with that
#' declaration and `check_signal_tails()` is what holds them to it.
#'
#' @param sp the `signal` element of `layer1_priors()`, or any frame with
#'   `signal`, `role`, `fold`, `alpha_low/mid/high`, `alpha0`
#' @param s  measured rows for the site being described. At the training site
#'   that is `.measured_train_rows()`; at an apply site, `measured_rows()`.
#' @param role which fits to summarise; "final" is the bundle's own parameters
dm_shrinkage_table <- function(sp, s, role = "final") {
  p <- sp[sp$role == role, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  # The constant `fit_dm_alpha()` floors at, read from its ONE declaration in
  # R/04 (`DM_ALPHA_FLOOR`), times the multiple its own `degenerate` test uses.
  # This was a second literal with a comment naming where it came from until
  # 2026-09-08, which is exactly the two-declarations drift audit findings F9
  # to F11 describe.
  ALPHA_FLOOR <- DM_ALPHA_FLOOR * 10
  rows <- lapply(seq_len(nrow(p)), function(i) {
    sg <- p$signal[i]
    z  <- s[s$signal == sg, , drop = FALSE]
    n  <- z$n_obs
    if (!length(n)) return(NULL)
    qn <- stats::quantile(n, c(0.25, 0.5, 0.75), names = FALSE)
    a0 <- p$alpha0[i]
    fl <- c(low = p$alpha_low[i], mid = p$alpha_mid[i],
            high = p$alpha_high[i]) <= ALPHA_FLOOR
    data.frame(
      signal    = sg,
      n_stays   = length(n),
      alpha_low = round(p$alpha_low[i], 4),
      alpha_mid = round(p$alpha_mid[i], 4),
      alpha_high = round(p$alpha_high[i], 4),
      alpha0    = round(a0, 4),
      n_q25 = qn[1], n_med = qn[2], n_q75 = qn[3],
      # posterior weight on the OBSERVED proportion at each quartile of coverage
      w_q25 = round(qn[1] / (qn[1] + a0), 4),
      w_med = round(qn[2] / (qn[2] + a0), 4),
      w_q75 = round(qn[3] / (qn[3] + a0), 4),
      # The weight a stay measured ONCE receives. Same column the delta table
      # reports, and the sharpest single reading of whether the prior is strong
      # enough to separate a confident zero from an unconfident one.
      w_at_n1 = round(1 / (1 + a0), 4),
      # Prior strength in units of typical coverage: 1.0 means the prior is
      # worth as much as the median stay's entire observation record.
      alpha0_per_n_med = round(a0 / max(qn[2], 1), 4),
      shrinks    = isTRUE(p$converged[i]) && !any(fl),
      floor_low  = unname(fl["low"]),
      floor_mid  = unname(fl["mid"]),
      floor_high = unname(fl["high"]),
      # Site-local, from this site's own counts. See the note above on why the
      # frozen alphas and the site's own tails are different questions.
      tail_empty_low  = max(z$k_low,  na.rm = TRUE) == 0,
      tail_empty_high = max(z$k_high, na.rm = TRUE) == 0,
      converged  = p$converged[i],
      degenerate = p$degenerate[i],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out[order(out$w_med), , drop = FALSE]
}

# --- 2. is the Dirichlet-multinomial the right model for the counts? --------

#' Posterior-predictive check on the DM, per (signal, category).
#'
#' Four numbers, each answering a distinct way the model could be wrong.
#'
#'   dispersion   mean of (k - E[k|n])^2 / Var[k|n]. Exactly 1 under correct
#'                specification. ABOVE 1 means the data are MORE variable than
#'                the beta-binomial allows even after its own overdispersion —
#'                the DM is too tight and alpha is being pulled towards a
#'                spurious consensus. BELOW 1 means the fitted alpha is smaller
#'                than the data support, so `pi_hat` shrinks less than it
#'                should. This is the headline number.
#'
#'   p0_obs / p0_exp   observed against expected fraction at k = 0. The DM has
#'                no zero-inflation component, so a large excess of exact zeros
#'                is the specific failure mode it cannot absorb — and it matters
#'                here because `k = 0` is the modal outcome for several signals
#'                and is exactly the cell the shrinkage is meant to handle.
#'
#'   pfull_obs / pfull_exp   the same at k = n, the saturated tail. An excess
#'                here says a subpopulation deviates in every covered hour,
#'                which a single alpha cannot represent.
#'
#'   pit_sd / pit_ks   randomised quantile residuals, transformed to the normal
#'                scale. Under correct specification they are exactly standard
#'                normal, so `pit_sd` near 1 and a small KS statistic mean the
#'                whole shape fits, not just the first two moments. These use
#'                the exact beta-binomial CDF, not a simulation.
#'
#' Randomisation makes the PIT continuous for discrete data (Dunn and Smyth);
#' `with_seed` keeps it reproducible rather than moving between runs.
#'
#' @param k integer matrix of counts, columns (low, mid, high)
#' @param n row totals, i.e. `n_obs`
#' @param alpha length-3 fitted Dirichlet parameter
#' @return one row per category
dm_ppc <- function(k, n, alpha, seed = 1L, labels = c("low", "mid", "high")) {
  k <- as.matrix(k)
  stopifnot(ncol(k) == length(alpha), nrow(k) == length(n))
  ok <- n > 0
  k <- k[ok, , drop = FALSE]; n <- n[ok]
  a0 <- sum(alpha)

  rows <- lapply(seq_along(alpha), function(j) {
    a <- alpha[j]; b <- a0 - a
    kj <- k[, j]
    if (b <= 0 || a <= 0) {
      return(data.frame(category = labels[j], n_stays = length(kj),
                        alpha_j = a, dispersion = NA_real_,
                        p0_obs = NA_real_, p0_exp = NA_real_,
                        pfull_obs = NA_real_, pfull_exp = NA_real_,
                        pit_mean = NA_real_, pit_sd = NA_real_, pit_ks = NA_real_,
                        stringsAsFactors = FALSE))
    }
    mu <- n * a / a0
    vr <- n * a * b * (a0 + n) / (a0^2 * (a0 + 1))
    disp <- mean((kj - mu)^2 / vr)

    p0_exp    <- mean(.dbetabinom(rep(0L, length(n)), n, a, b))
    pfull_exp <- mean(.dbetabinom(n, n, a, b))

    z <- with_seed(seed, {
      lo <- .pbetabinom(kj - 1L, n, a, b)
      hi <- .pbetabinom(kj,      n, a, b)
      u  <- lo + stats::runif(length(kj)) * (hi - lo)
      stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
    })
    ks <- suppressWarnings(stats::ks.test(z, "pnorm")$statistic)

    data.frame(
      category  = labels[j],
      n_stays   = length(kj),
      alpha_j   = round(a, 4),
      dispersion = round(disp, 4),
      p0_obs    = round(mean(kj == 0), 4),
      p0_exp    = round(p0_exp, 4),
      pfull_obs = round(mean(kj == n), 4),
      pfull_exp = round(pfull_exp, 4),
      pit_mean  = round(mean(z), 4),
      pit_sd    = round(stats::sd(z), 4),
      pit_ks    = round(unname(ks), 4),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Does the DM's fit HOLD ACROSS THE COVERAGE RANGE?
#'
#' ADDED 2026-09-05, and it is the check that would have caught the `delta` and
#' `lambda` shrinkage failures years earlier than it did. Both of those rested
#' on an assumed relationship with the information weight, and both assumed
#' relationships turned out to run the wrong way in the data. This asks the same
#' question of the one construct that survives.
#'
#' A single alpha per signal is applied across a 2x to 24x range of `n_obs`.
#' If its dispersion drifts with coverage, that alpha calibrates one coverage
#' regime better than the other -- which is a specification limitation to report
#' rather than a defect to fix, since the natural remedy is an alpha that varies
#' with `n`, and that would make coverage a first-class modelled quantity in a
#' design built to keep it out.
#'
#' MEASURED at MIMIC 2026-09-05: median ratio of high-`n` to low-`n` dispersion
#' 1.244, with real spread -- `lactate` 0.10, `bicarbonate` 2.53,
#' `temperature` 2.18. Reported, not acted on. Read it at BOTH sites: eICU's
#' coverage distribution is different, so this table differenced between the
#' two is a transport reading.
#'
#' EVERY DECLARED OCCUPIABLE TAIL IS CHECKED, as of 2026-09-09 (statistical
#' review, lower-priority items). The check read `alpha_low` and `k_low` for
#' every signal, although the model carries three categories and several
#' signals enter their formula through the HIGH coordinate. The categories
#' come from `occupiable_tails_of()` -- the same accessor the formula builder
#' and the prior-fit policy read -- so a tail that is declared modelled is a
#' tail that is checked. `ratio` far from 1 in either direction is the finding.
#'
#' @return one row per (signal, declared tail), or NULL
dm_dispersion_drift <- function(sp, s, cfg, role = "final", n_bins = 4L) {
  p <- sp[sp$role == role, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  grid <- do.call(rbind, lapply(seq_len(nrow(p)), function(i) {
    tails <- intersect(c("low", "high"), occupiable_tails_of(p$signal[i], cfg))
    if (!length(tails)) return(NULL)
    data.frame(i = i, category = tails, stringsAsFactors = FALSE)
  }))
  if (is.null(grid)) return(NULL)
  rows <- lapply(seq_len(nrow(grid)), function(gi) {
    i <- grid$i[gi]; cat_ <- grid$category[gi]
    sg <- p$signal[i]
    z <- s[s$signal == sg, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    al <- c(low = p$alpha_low[i], mid = p$alpha_mid[i], high = p$alpha_high[i])
    a0 <- sum(al)
    a <- unname(al[cat_]); b <- a0 - a
    if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NULL)
    n <- z$n_obs; kj <- z[[paste0("k_", cat_)]]
    qs <- unique(stats::quantile(n, seq(0, 1, length.out = n_bins + 1L),
                                 names = FALSE))
    if (length(qs) < 3L) return(NULL)
    g <- cut(n, breaks = qs, include.lowest = TRUE, labels = FALSE)
    disp <- vapply(sort(unique(g)), function(bb) {
      sel <- g == bb
      if (sum(sel) < 50L) return(NA_real_)
      mu <- n[sel] * a / a0
      vr <- n[sel] * a * b * (a0 + n[sel]) / (a0^2 * (a0 + 1))
      mean((kj[sel] - mu)^2 / vr)
    }, numeric(1))
    disp <- disp[is.finite(disp)]
    if (length(disp) < 3L) return(NULL)
    data.frame(signal = sg, category = cat_, n_bins = length(disp),
               n_med = stats::median(n),
               disp_lo_n = round(disp[1], 4),
               disp_hi_n = round(disp[length(disp)], 4),
               ratio = round(disp[length(disp)] / max(disp[1], 1e-9), 4),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  out[order(-abs(log(pmax(out$ratio, 1e-9)))), , drop = FALSE]
}

#' `dm_ppc()` over every signal, using each signal's own final alpha.
dm_ppc_all <- function(sp, s, role = "final", seed = 1L) {
  p <- sp[sp$role == role, , drop = FALSE]
  rows <- lapply(seq_len(nrow(p)), function(i) {
    sg <- p$signal[i]
    z  <- s[s$signal == sg, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    out <- dm_ppc(cbind(z$k_low, z$k_mid, z$k_high), z$n_obs,
                  c(p$alpha_low[i], p$alpha_mid[i], p$alpha_high[i]), seed = seed)
    cbind(signal = sg, out)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out
}

# --- 3. do the conditional means have the right FORM? -----------------------

#' Reconstruct the residual a stored delta row implies, under either form.
#'
#' ONE implementation, used by every diagnostic below, so a diagnostic can never
#' measure a different residual from the one the model frame actually carries.
#' It deliberately stops short of any weighting: the form question is about
#' the CONDITIONAL MEAN, and the precision layer would mask a bad one by
#' pulling every residual toward zero.
.delta_residual <- function(par, v, k, n) .delta_raw_resid(par, v, k, n)

#' The declared decoupling-smooth specification: basis, fixed k, and the cap
#' the adaptive variant uses.
#'
#' Declared in `config/config.yml` under `diagnostics.decoupling_smooth` from
#' 2026-09-03 (audit finding F11). It was `bs = "ts", k = 6` written into this
#' file, and the R-squared it produces carries a frozen decision -- keeping
#' `delta` on for the three GCS components rests on it -- so it belongs where
#' the other bases are declared. `cfg_req()` rather than a fallback, for the
#' same reason `bam_settings()` uses it: a fallback would be a second
#' declaration of the same number.
decoupling_smooth_spec <- function(cfg) {
  what <- "declared in config/config.yml under `diagnostics.decoupling_smooth`"
  list(basis = cfg_req(cfg, "diagnostics", "decoupling_smooth", "basis", what = what),
       k     = as.integer(cfg_req(cfg, "diagnostics", "decoupling_smooth", "k", what = what)),
       k_max = as.integer(cfg_req(cfg, "diagnostics", "decoupling_smooth", "k_max", what = what)))
}

#' One decoupling smooth, built from the declared spec.
#'
#' The formula is assembled with `sprintf()` for the same reason
#' `R/05_formula.R` does it: the basis and dimension are data-free DESIGN and
#' must come from config rather than being written into a literal. `data=` is
#' supplied explicitly so the fit does not depend on what happens to be in the
#' calling frame.
.decoupling_gam <- function(rr, xx, ds) {
  f <- stats::as.formula(sprintf("rr ~ s(xx, bs = \"%s\", k = %d)", ds$basis, ds$k))
  mgcv::gam(f, data = data.frame(rr = rr, xx = xx))
}

#' Direction of a residual's spread against its information weight.
#'
#' The one-line summary of the check that killed both shrinkage constructs.
#' Strata the weight into quantile bins, take the spread within each, and
#' correlate it with the stratum's median weight. NEGATIVE is what a
#' variance-component model asserts -- more information, less spread. POSITIVE
#' means the spread grows with information, which cannot be measurement noise
#' and is therefore not something to shrink away.
#'
#' Robust (IQR) by default because a variance is dominated by the tails and the
#' question here is about the bulk.
#'
#' @return a Spearman correlation, or NA when there are too few usable strata
.spread_vs_n <- function(r, w, n_bins = 5L, min_stratum = 30L, robust = TRUE) {
  ok <- is.finite(r) & is.finite(w)
  r <- r[ok]; w <- w[ok]
  if (length(r) < 200L || length(unique(w)) < 3L) return(NA_real_)

  # STRATIFY BY DISTINCT VALUE WHERE THE WEIGHT IS A SMALL INTEGER, and by
  # quantile otherwise. This is not a refinement, it is a correctness fix made
  # on 2026-09-05 after the column came back NA for four of five interventions
  # on its first run.
  #
  # THE MECHANISM, because it will recur. Quantile stratification collapses
  # when the weight is concentrated on one value: `lambda`'s exposure `D` has a
  # median of 1 and quartiles 1/1/1/2, so `quantile(D, 0:5/5)` returns
  # 1,1,1,1,2,15, `unique()` drops the repeats, and two strata survive where
  # three are needed. The exploratory version of this check iterated distinct
  # values and never had the problem; generalising it to quantiles on the way
  # into the library silently lost the property that made it work.
  uw <- sort(unique(w))
  if (length(uw) <= 20L) {
    g <- match(w, uw)
  } else {
    qs <- unique(stats::quantile(w, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
    if (length(qs) < 3L) return(NA_real_)
    g <- cut(w, breaks = qs, include.lowest = TRUE, labels = FALSE)
  }
  tab <- do.call(rbind, lapply(sort(unique(g)), function(b) {
    sel <- g == b
    if (sum(sel) < min_stratum) return(NULL)
    data.frame(w_med = stats::median(w[sel]),
               spread = if (robust) stats::IQR(r[sel]) else stats::var(r[sel]))
  }))
  if (is.null(tab) || nrow(tab) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(tab$w_med, tab$spread, method = "spearman"))
}

#' Per-(signal, variable) fit quality for the magnitude construct, either form.
#'
#' `resid_smooth_r2` is the R-squared of a SMOOTH of the residual on `log1p(k)`.
#' Near zero means the conditional mean captured the count relationship.
#' Materially above zero means it did not, and the leftover is exactly the count
#' dependence `delta` exists to remove -- a DIFFERENT statement from the
#' decoupling in `conditional_priors.md` SS7.1, because that measures the
#' finished covariate and this measures the model that produced it.
#'
#' `resid_step_r2` is the same statistic computed WITHIN the k > 0 rows only.
#' It exists because the pooled number can be small while the model is still
#' wrong in one regime, and that is precisely how the Glasgow misspecification
#' hid: pooled 0.0999 for gcs_eyes against 0.242 inside k > 0.
#'
#' `shrink_frac` is the mean posterior shrinkage weight `s_u / (s_u + s_e / n)`.
#' At 1 the construct returns the raw residual because the variance split did
#' not earn its parameter; below 1 it is pulled toward the conditional mean.
#'
#' `offscale_frac` is the fraction of stays whose fitted conditional mean falls
#' outside the response's observed range. It is 0 by construction under the
#' ordinal form and is the single clearest symptom of the linear form applied to
#' a bounded scale.
delta_fit_diagnostics <- function(mp, s, cfg, role = "final") {
  p <- mp[mp$role == role, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  ds <- decoupling_smooth_spec(cfg)
  rows <- lapply(seq_len(nrow(p)), function(i) {
    sg <- p$signal[i]; vr <- p$variable[i]; cv <- p$count_var[i]
    z <- s[s$signal == sg, , drop = FALSE]
    if (!nrow(z) || !all(c(vr, cv) %in% names(z))) return(NULL)
    v <- z[[vr]]; k <- z[[cv]]; n <- z$n_obs
    keep <- !is.na(v) & !is.na(k) & !is.na(n) & n > 0
    v <- v[keep]; k <- k[keep]; n <- n[keep]
    if (length(v) < 50L || stats::sd(v) == 0) return(NULL)

    par <- as.list(p[i, , drop = FALSE])
    r <- .delta_residual(par, v, k, n)
    if (is.null(r)) return(NULL)

    lk <- log1p(k)
    # gam() rather than lm(), because a LINEAR residual-on-lk R-squared is
    # guaranteed ~0 by an OLS fit and would say nothing; curvature is the
    # question. Under the ordinal form the fit is not OLS at all, so even the
    # linear part is informative -- another reason to use the same smooth.
    smooth_r2 <- function(rr, xx) {
      if (length(rr) < 100L || stats::sd(xx) == 0) return(NA_real_)
      f <- try(.decoupling_gam(rr, xx, ds), silent = TRUE)
      if (inherits(f, "try-error")) NA_real_ else summary(f)$r.sq
    }
    pos <- k > 0
    sw <- if (isTRUE(p$shrinks[i]) && p$s_u[i] > 0)
      mean(p$s_u[i] / (p$s_u[i] + p$s_e[i] / n)) else if (isTRUE(p$shrinks[i])) 0 else 1

    off <- NA_real_
    if (identical(par$form %||% "linear", "linear")) {
      # THE FITTED MEAN, NOT A PRIVATE COPY OF IT (statistical review S7). `r`
      # is `v` minus the complete fitted predictor, `a3 * I(k > 0)` included,
      # so `v - r` IS the conditional mean the frame carries. The previous
      # expression omitted the step term and, on a synthetic model with a
      # step of 10, reported every stay off-scale while the true fitted mean
      # put none there.
      fv <- v - r
      off <- mean(fv < min(v) | fv > max(v))
    } else {
      off <- 0
    }

    data.frame(
      signal = sg, variable = vr, count_var = cv,
      form = par$form %||% "linear", n_stays = length(v),
      a1 = round(p$a1[i], 5), a2 = round(p$a2[i], 5),
      a3 = round(p$a3[i] %||% NA_real_, 5),
      s_e = round(p$s_e[i], 6), s_u = round(p$s_u[i], 6),
      # RENAMED FROM `shrinks` 2026-09-05. The underlying priors column keeps
      # that name -- it records what `.fit_var_components()` decided and
      # renaming it would change the `priors` target's value and rebuild all
      # 258 fits for nothing. But NOTHING SHRINKS ANY MORE, so emitting a
      # column called `shrinks` from a reporting table would be a third
      # declaration of a thing the code does not do. What the flag actually
      # says is whether within-stay noise was IDENTIFIABLE: whether the
      # two-component variance model beat the one-component model on AIC.
      # That is still worth reporting; it just no longer rescales anything.
      noise_identified = p$shrinks[i], shrink_frac_would_be = round(sw, 4),
      # WHERE s_e CAME FROM. `shrinks = FALSE` is ambiguous on its own -- it can
      # mean the value genuinely does not move within a stay (correct, weight 1)
      # or that the component was not identified (a silent failure). That
      # ambiguity is how the inert machinery went unnoticed until 2026-09-02, so
      # the reason travels with the number. "replicated" is the healthy state;
      # "heteroscedastic_legacy" means the feature table carries no `se_`
      # columns and the old estimator ran.
      se_source = p$se_source[i] %||% NA_character_,
      se_df = p$se_df[i] %||% NA_real_,
      # The weight at a single hour, which is the quantity the whole construct
      # exists to get right: it is what a stay measured once is shrunk by, and
      # it was 1.000 on 22 of 38 rows before the estimator changed.
      w_at_n1_would_be = round(if (isTRUE(p$shrinks[i]) && p$s_u[i] > 0)
                        p$s_u[i] / (p$s_u[i] + p$s_e[i]) else 1, 4),
      resid_smooth_r2 = round(smooth_r2(r, lk), 4),
      resid_step_r2 = round(smooth_r2(r[pos], lk[pos]), 4),
      frac_k_pos = round(mean(pos), 4),
      offscale_frac = round(off, 4),
      resid_sd = round(stats::sd(r), 5),
      # DOES THE RESIDUAL SPREAD FALL AS COVERAGE RISES? Added 2026-09-05 with
      # the shrinkage removal. The variance-component model that used to weight
      # this residual asserted that it must; measured, it does not, in 17 of 22
      # rows whose magnitude variable is an extreme. The weight is gone, so this
      # is no longer load-bearing -- it is kept so the assumption that was wrong
      # here is CHECKED at every site rather than assumed, and so that a future
      # site where the law genuinely holds is visible rather than inferred.
      # Spearman of stratum spread against stratum median n: negative is the
      # old model's assertion, positive is the inverted direction.
      spread_rho_iqr = round(.spread_vs_n(r, n, robust = TRUE), 3),
      spread_rho_var = round(.spread_vs_n(r, n, robust = FALSE), 3),
      degenerate = p$degenerate[i],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  if (!is.null(out)) out <- out[order(-out$resid_smooth_r2), , drop = FALSE]
  out
}

#' Ordinal-form-specific checks: calibration, uniformity, proportional odds.
#'
#' Three questions the pooled residual R-squared cannot answer, one column each.
#'
#'   CALIBRATION. Under a correct model the expected category counts match the
#'   observed ones. `cal_max_abs_diff` is the largest absolute discrepancy as a
#'   fraction of stays, over the declared levels. A large value means the cut
#'   points are in the wrong place, which the residual smooth would not see.
#'
#'   UNIFORMITY. A RANDOMISED PIT, `F(y-) + V p(y)` with V ~ Uniform(0, 1)
#'   under the caller's seed, is exactly Uniform(0, 1) when the model is right,
#'   so `pit_ks` is the Kolmogorov distance from uniform and `pit_sd` sits
#'   near 1/sqrt(12) = 0.2887 under a correct fit. CHANGED 2026-09-09
#'   (statistical review S6): the diagnostic used the MIDPOINT transform
#'   `F(y-) + 0.5 p(y)`, which is not uniform for a discrete outcome even under
#'   the true model -- a perfectly specified balanced binary outcome puts it
#'   at 0.25 and 0.75, KS distance 0.25 -- and whose discrepancy depends on the
#'   category masses, so comparing raw KS across signals mixed fit quality with
#'   discreteness. The midpoint-logit DELTA FEATURE is unchanged; only the
#'   reference distribution of the diagnostic moved.
#'
#'   PROPORTIONAL ODDS. `po_slope_spread` is the range of the `log(1+k)` slope
#'   across separate binary logits at each cut point, EACH FITTED WITH THE
#'   ORDINAL MODEL'S OWN COVARIATES (`log(1+k)`, `log(n)`, `I(k > 0)`, with the
#'   same constant-column rule `fit_delta_ordinal()` applies). The threshold
#'   fits omitted the step term until 2026-09-09 (review S7), so their slope
#'   spread was partly the step being absorbed into the slope rather than a
#'   test of proportional odds. `po_slope_se_max` is the largest standard error
#'   among those slopes: a spread that is inside it is not identified, and it
#'   is reported so the spread cannot be read on its own. Under proportional
#'   odds these are one number. A spread that is large RELATIVE to the fitted
#'   `a1` AND to its standard error is the signal to escalate to per-cut-point
#'   slopes -- which, MEASURED 2026-09-01, is NOT what the Glasgow triple
#'   needed: the two-regime latent predictor beat partial proportional odds on
#'   both this statistic and AIC.
#'
#' AGGREGATES ONLY (hard rule 1): counts, fractions and distances, never a row.
#'
#' @param seed for the randomised PIT. Required: the table must reproduce.
delta_ordinal_diagnostics <- function(mp, s, cfg, role = "final", seed) {
  if (missing(seed)) {
    stop("delta_ordinal_diagnostics: `seed` is required for the randomised PIT",
         call. = FALSE)
  }
  p <- mp[mp$role == role & (mp$form %||% "linear") == "ordinal", , drop = FALSE]
  if (!nrow(p)) return(NULL)

  rows <- lapply(seq_len(nrow(p)), function(i) {
    sg <- p$signal[i]; vr <- p$variable[i]; cv <- p$count_var[i]
    z <- s[s$signal == sg, , drop = FALSE]
    if (!nrow(z) || !all(c(vr, cv) %in% names(z))) return(NULL)
    v <- z[[vr]]; k <- z[[cv]]; n <- z$n_obs
    keep <- !is.na(v) & !is.na(k) & !is.na(n) & n > 0
    v <- v[keep]; k <- k[keep]; n <- n[keep]
    par <- as.list(p[i, , drop = FALSE])
    lv <- .par_levels(par); th <- .par_theta(par)
    yc <- match(v, lv)
    if (anyNA(yc) || !length(th)) return(NULL)

    eta <- par$a1 * log1p(k) + par$a2 * log(n) + (par$a3 %||% 0) * as.numeric(k > 0)
    pr  <- .polr_probs(th, eta)
    exp_n <- colSums(pr)
    obs_n <- tabulate(yc, nbins = length(lv))
    cal <- max(abs(obs_n - exp_n)) / length(v)

    below <- .polr_cum_below(th, eta)[cbind(seq_along(yc), yc)]
    pmf   <- .polr_probs(th, eta)[cbind(seq_along(yc), yc)]
    u  <- with_seed(seed + i, below + stats::runif(length(yc)) * pmf)
    u  <- pmin(pmax(u, 0), 1)
    ks <- suppressWarnings(stats::ks.test(u, "punif")$statistic)

    lk <- log1p(k); ln <- log(n); pos <- as.numeric(k > 0)
    dd <- data.frame(lk = lk)
    fm <- "yb ~ lk"
    if (stats::sd(ln)  > 0) { dd$ln  <- ln;  fm <- paste(fm, "+ ln") }
    if (stats::sd(pos) > 0) { dd$pos <- pos; fm <- paste(fm, "+ pos") }
    sl <- t(vapply(lv[-length(lv)], function(cc) {
      dd$yb <- as.integer(v <= cc)
      if (length(unique(dd$yb)) < 2L) return(c(NA_real_, NA_real_))
      g <- try(suppressWarnings(stats::glm(stats::as.formula(fm), data = dd,
                                           family = stats::binomial())),
               silent = TRUE)
      if (inherits(g, "try-error")) return(c(NA_real_, NA_real_))
      cf <- stats::coef(g); vc <- try(stats::vcov(g), silent = TRUE)
      c(unname(cf["lk"]),
        if (inherits(vc, "try-error")) NA_real_ else sqrt(unname(vc["lk", "lk"])))
    }, numeric(2)))
    spread <- if (all(is.na(sl[, 1]))) NA_real_ else diff(range(sl[, 1], na.rm = TRUE))
    se_max <- if (all(is.na(sl[, 2]))) NA_real_ else max(sl[, 2], na.rm = TRUE)

    data.frame(
      signal = sg, variable = vr, n_levels = length(lv), n_stays = length(v),
      theta = paste(sprintf("%.2f", th), collapse = " "),
      a1 = round(par$a1, 4), a2 = round(par$a2, 4), a3 = round(par$a3 %||% NA_real_, 4),
      cal_max_abs_diff = round(cal, 5),
      pit_type = "randomized",
      pit_ks = round(unname(ks), 4), pit_sd = round(stats::sd(u), 4),
      po_slope_spread = round(spread, 3),
      po_slope_se_max = round(se_max, 3),
      po_spread_rel_a1 = round(abs(spread / par$a1), 3),
      levels_unobserved = sum(obs_n == 0L),
      converged = p$converged[i],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out
}

#' `lambda`'s fitted parameters, per intervention, with the transport reading.
#'
#' No residual recomputation here: the binomial arm's `c1` and the lognormal
#' arm's `b1` ARE the diagnostic. They are two-parameter summaries of a
#' population's treatment behaviour — how sharply molecule count or accumulated
#' amount rises with exposure duration — and `v2_state` §6 already identifies
#' them as the cleanest cross-site quantity the project has, because they are
#' directly comparable between sites and independent of per-term attribution.
#'
#' Run this at both sites and difference the table. That IS the lambda transport
#' result, and it costs nothing beyond loading eICU.
#' @param ivf optional intervention-features rows for the site being described,
#'   and `cfg` alongside them. Supplying both adds the residual summary and the
#'   spread-versus-exposure direction, which is the check that killed this
#'   construct's shrinkage: the weight required residual spread to FALL as
#'   exposure hours rise, and it RISES for diuretic and transfusion_prbc.
#'   Optional so a caller with no feature table still gets the parameter
#'   summary.
#'
#'   THE RESIDUAL IS THE FROZEN ONE (statistical review S8, 2026-09-09). Until
#'   now the spread columns refitted `lm(log(A) ~ log(D))` on the supplied
#'   rows and read residuals from THAT fit, so at eICU they described spread
#'   after local recentring and could not see a conditional-mean shift, while
#'   the table displayed the frozen `b0`/`b1` beside them. The residual is now
#'   `log(A) - (b0 + b1 log D)` from the row's own frozen coefficients -- the
#'   quantity `lambda_value()` actually scores -- and `resid_mean` is its site
#'   mean: zero by construction on the fitting rows, and at an apply site the
#'   size of the conditional-mean transport failure. The caller decides which
#'   rows are the site; at MIMIC the graph passes TRAINING rows only.
lambda_fit_diagnostics <- function(ip, role = "final", ivf = NULL, cfg = NULL) {
  p <- ip[ip$role == role, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  out <- data.frame(
    intervention = p$intervention, family = p$family, n_stays = p$n_stays,
    M = p$M,
    c0 = round(p$c0, 5), c1 = round(p$c1, 5),
    b0 = round(p$b0, 5), b1 = round(p$b1, 5),
    s_e = round(p$s_e, 6), s_u = round(p$s_u, 6),
    # See the note in `delta_fit_diagnostics()`: renamed because nothing
    # shrinks any more, and the flag's real content is whether within-stay
    # noise was identifiable at all.
    noise_identified = p$shrinks, converged = p$converged,
    degenerate = p$degenerate,
    stringsAsFactors = FALSE)

  # The frozen residual summary and the spread-versus-exposure direction,
  # where the data to compute them exists.
  out$n_exposed_site <- NA_integer_
  out$resid_mean     <- NA_real_
  out$resid_sd       <- NA_real_
  out$spread_rho_iqr <- NA_real_
  out$spread_rho_var <- NA_real_
  if (!is.null(ivf) && !is.null(cfg)) {
    for (i in seq_len(nrow(out))) {
      iv <- out$intervention[i]
      sp <- tryCatch(lambda_spec_of(iv, cfg), error = function(e) NULL)
      if (is.null(sp) || !identical(sp$family, "lognormal")) next
      if (!is.finite(p$b0[i]) || !is.finite(p$b1[i])) next
      z <- ivf[ivf$intervention == iv, , drop = FALSE]
      if (!nrow(z) || !all(c(sp$exposure, sp$accumulation) %in% names(z))) next
      D <- z[[sp$exposure]] * sp$exposure_scale; A <- z[[sp$accumulation]]
      ex <- is.finite(D) & is.finite(A) & D > 0 & A > 0
      out$n_exposed_site[i] <- sum(ex)
      if (sum(ex) < 200L) next
      d <- D[ex]
      r <- log(A[ex]) - (p$b0[i] + p$b1[i] * log(d))
      out$resid_mean[i] <- round(mean(r), 5)
      out$resid_sd[i]   <- round(stats::sd(r), 5)
      out$spread_rho_iqr[i] <- round(.spread_vs_n(r, d, robust = TRUE), 3)
      out$spread_rho_var[i] <- round(.spread_vs_n(r, d, robust = FALSE), 3)
    }
  }
  rownames(out) <- NULL
  out
}

# --- 4. how much do the frozen parameters move between folds? ---------------

#' Fold-to-fold spread of every fitted parameter.
#'
#' THE CHEAPEST AVAILABLE PROXY FOR THE TRANSPORT QUESTION, and it needs no eICU
#' data. Each parameter is fitted five times on overlapping four-fifths of the
#' same cohort, so the spread across folds is a lower bound on how much it would
#' move under a genuinely different population. A parameter that already wanders
#' between folds of ONE site will not survive a change of site, and that is
#' knowable now rather than after the external run.
#'
#' Reported as a coefficient of variation so parameters on different scales are
#' comparable. This is a screening number, not a test: it says which parameters
#' to look at first when the eICU comparison eventually runs.
#'
#' READ `cv` BESIDE `sd`, `min` AND `max`, NEVER ON ITS OWN. This is the single
#' most important thing about this table and the three absolute columns were
#' added on 2026-09-05 because it was possible to read it wrongly.
#'
#' A coefficient of variation is `sd / mean(|x|)`, so it is scale-free -- which
#' is the whole point, and is also its failure mode. THE DENOMINATOR GOES TO
#' ZERO FOR A PARAMETER THAT IS ITSELF NEAR ZERO, and the ratio then diverges
#' while the parameter has barely moved. Every row at the top of this table at
#' MIMIC is that case rather than a wandering parameter: `gcs_motor/value_min`'s
#' `a1` has CV 1.37 and a fold-to-fold standard deviation of 0.009 on a scale
#' whose response runs 1 to 6, and `inotrope`'s `c1` has CV 0.57 and moves
#' between 0.005 and 0.023 while `vasopressor`'s sits at 0.285 to 0.290.
#'
#' The correct reading of a large CV with a tiny `sd` is NOT "this parameter is
#' unstable". It is "this parameter is indistinguishable from zero", which is a
#' statement about whether the term carries any signal at all -- often a more
#' interesting one, and a completely different one. A parameter that changes
#' SIGN across folds, which `min` and `max` now make visible, is the sharpest
#' version of it.
#'
#' The row that would genuinely worry is a large CV with an `sd` that is large
#' relative to the covariate's own scale. There is none at MIMIC.
prior_fold_stability <- function(priors) {
  cv <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L || mean(abs(x)) == 0) return(NA_real_)
    stats::sd(x) / mean(abs(x))
  }
  # One row builder, so the three blocks cannot come to report different
  # summaries of the same thing.
  row1 <- function(block, key, parameter, x) {
    x <- x[is.finite(x)]
    data.frame(block = block, key = key, parameter = parameter,
               n_folds = length(x), mean = round(mean(x), 5),
               sd = round(stats::sd(x), 5), cv = round(cv(x), 4),
               min = round(min(x), 5), max = round(max(x), 5),
               # TRUE means the parameter is not distinguishable from zero
               # across folds, which is why its CV is large. It is the flag
               # that stops a scale artefact being read as instability.
               sign_flip = min(x) < 0 && max(x) > 0,
               stringsAsFactors = FALSE)
  }
  out <- list()

  sp <- priors$signal[priors$signal$role == "oof", , drop = FALSE]
  if (nrow(sp)) {
    out$alpha <- do.call(rbind, lapply(split(sp, sp$signal), function(z)
      row1("alpha", z$signal[1], "alpha0", z$alpha0)))
  }

  mp <- priors$magnitude[priors$magnitude$role == "oof", , drop = FALSE]
  if (nrow(mp)) {
    mp$key <- paste(mp$signal, mp$variable, sep = "/")
    out$delta <- do.call(rbind, lapply(split(mp, mp$key), function(z)
      do.call(rbind, lapply(c("a1", "a2"), function(v)
        row1("delta", z$key[1], v, z[[v]])))))
  }

  ip <- priors$intervention[priors$intervention$role == "oof", , drop = FALSE]
  if (nrow(ip)) {
    out$lambda <- do.call(rbind, lapply(split(ip, ip$intervention), function(z) {
      v <- if (identical(z$family[1], "binomial")) "c1" else "b1"
      row1("lambda", z$intervention[1], v, z[[v]])
    }))
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  if (!is.null(res)) res <- res[order(-res$cv), , drop = FALSE]
  res
}
