# R/04_features.R ------------------------------------------------------------
# Two jobs, joined by one row-scoping rule:
#
#   1. Dirichlet-multinomial shrinkage of the (k_low, k_mid, k_high) triples
#      into pi_minus / pi_plus, with alpha fitted on TRAINING rows only.
#   2. The model-frame builder: one narrow frame per (signal, model), built to
#      exactly the column list R/05_formula.R's required_columns() emits.
#
# They live together because both are computed over the same set of rows — the
# measured subset (`n_obs > 0`) of a training population — and if that scoping
# were written twice the two could disagree without anything failing. See
# training_priors(), which walks it once and returns both.
#
# No paths, no clock, no data access (hard rule 9). Long format in, long format
# out: this file never builds a global wide frame (hard rule 2). The only
# widening it does is per-signal and on demand.
# ----------------------------------------------------------------------------

# --- Dirichlet-multinomial shrinkage ----------------------------------------
#
# For signal g, stay i with counts k_i = (k_low, k_mid, k_high) over n_i = n_obs
# covered hours:
#
#     k_i ~ DirichletMultinomial(n_i, alpha_g)
#
# and the shrunk propensity is the posterior mean
#
#     pi_hat_ij = (k_ij + alpha_gj) / (n_i + alpha_g0),   alpha_g0 = sum_j alpha_gj
#
# alpha_g0 is the strength of the prior in units of hours: a stay with n_obs = 2
# is pulled hard toward the signal-wide base rate, a stay with n_obs = 24 barely
# moves. That is the entire point — a single low lactate in 24 hours is not the
# same evidence as twenty of them, and raw proportions cannot tell them apart.
#
# alpha is fitted by maximum likelihood, never fixed by hand: the right amount
# of shrinkage differs by an order of magnitude between a q4h-sampled lab and a
# continuously-monitored vital, and guessing it would put a tuning constant in
# the middle of the evidence scale.

#' psi(v + a) - psi(a), tabulated for v = 0, 1, ..., vmax.
#'
#' Returned as a vector whose (v+1)-th element is the value at v.
#'
#' Computed as the partial harmonic sum sum_{t=0}^{v-1} 1/(a + t) rather than by
#' differencing two digammas. The two are equal in exact arithmetic, but when a
#' is small — which is exactly what the fit returns for a rare category, e.g. the
#' high tail of a signal that is almost never above its reference range —
#' digamma(a) diverges like -1/a and the difference of two large numbers loses
#' every significant digit. The harmonic form has no cancellation at all.
.psi_shift <- function(vmax, a) {
  if (vmax < 1L) return(0)
  c(0, cumsum(1 / (a + seq.int(0, vmax - 1L))))
}

#' Fit a Dirichlet-multinomial alpha by Minka's fixed-point iteration.
#'
#' @param k        integer matrix, one row per stay, one column per category.
#'                 Rows with rowSums == 0 contribute nothing to the likelihood
#'                 and are dropped.
#' @param tol      relative convergence tolerance on alpha
#' @param max_iter iteration cap
#' @param init     starting alpha; the default is uninformative-but-proper
#' @return list(alpha, alpha0, n, iter, converged, degenerate)
#'
#' The update is
#'
#'   alpha_j <- alpha_j * [ sum_i psi(k_ij + alpha_j) - N psi(alpha_j) ]
#'                      / [ sum_i psi(n_i + alpha_0) - N psi(alpha_0) ]
#'
#' which is a minorise-maximise step, so the likelihood is non-decreasing at
#' every iteration and no line search or gradient safeguard is needed.
#'
#' Counts are bounded (n_obs <= 24 by construction), so both sums are evaluated
#' over the tabulated support rather than over stays: each iteration costs O(25)
#' regardless of whether the fold holds 8,000 rows or 40,000. That is why
#' max_iter is generous — the whole fit is ~20 ms, and MM convergence is linear,
#' so a tight cap buys nothing and costs a spurious non-convergence flag.
#' MEASURED 2026-08-25: the 19 signals need 133-687 iterations at tol = 1e-8.
#' The floor `fit_dm_alpha()` holds every alpha coordinate above, and the ONE
#' declaration of it. `dm_shrinkage_table()` (R/04c) and `check_prior_fits()`
#' below both test coordinates against ten times this value, which is the
#' multiple `fit_dm_alpha()`'s own `degenerate` test uses; until 2026-09-08 the
#' diagnostic carried its own literal `1e-10 * 10` with a comment saying where it
#' came from, which is the second-declaration pattern audit findings F9 to F11
#' name.
DM_ALPHA_FLOOR <- 1e-10

fit_dm_alpha <- function(k, tol = 1e-8, max_iter = 10000L, init = NULL) {
  k <- as.matrix(k)
  storage.mode(k) <- "double"
  if (anyNA(k)) stop("fit_dm_alpha: count matrix has NA", call. = FALSE)
  if (any(k < 0)) stop("fit_dm_alpha: negative counts", call. = FALSE)

  n_i <- rowSums(k)
  k <- k[n_i > 0, , drop = FALSE]
  n_i <- n_i[n_i > 0]
  J <- ncol(k)
  N <- nrow(k)
  if (N == 0L) stop("fit_dm_alpha: no rows with a positive count total", call. = FALSE)

  vmax_j <- apply(k, 2L, max)
  vmax_n <- max(n_i)

  # Tabulated support: cnt[[j]][v+1] = #{i : k_ij == v}; tot[v+1] = #{i : n_i == v}
  cnt <- lapply(seq_len(J), function(j) tabulate(k[, j] + 1L, nbins = vmax_j[j] + 1L))
  tot <- tabulate(n_i + 1L, nbins = vmax_n + 1L)

  alpha <- if (is.null(init)) rep(1, J) else as.numeric(init)
  stopifnot(length(alpha) == J, all(alpha > 0))

  FLOOR <- DM_ALPHA_FLOOR   # keeps a vanishing category from becoming exactly 0
  iter <- 0L; converged <- FALSE
  while (iter < max_iter) {
    iter <- iter + 1L
    a0 <- sum(alpha)
    den <- sum(tot * .psi_shift(vmax_n, a0))
    if (!is.finite(den) || den <= 0) break
    num <- vapply(seq_len(J), function(j)
      sum(cnt[[j]] * .psi_shift(vmax_j[j], alpha[j])), numeric(1))
    new <- pmax(alpha * num / den, FLOOR)
    if (!all(is.finite(new))) break
    if (max(abs(new - alpha) / alpha) < tol) { alpha <- new; converged <- TRUE; break }
    alpha <- new
  }

  # A category never observed in this training set has MLE alpha_j -> 0, which
  # makes pi_hat_j identically 0 and the corresponding smooth a constant. It is
  # the honest estimate, but it is also a fit that will fail inside mgcv rather
  # than here, so it is reported.
  degenerate <- alpha <= FLOOR * 10 | vmax_j == 0

  list(alpha = alpha, alpha0 = sum(alpha), n = N, iter = iter,
       converged = converged, degenerate = degenerate)
}

#' Posterior-mean propensities from counts and a fitted alpha.
#'
#' @return matrix with columns pi_minus, pi_mid, pi_plus. Rows sum to 1 exactly,
#'   which is why only two of the three ever enter a formula.
shrink_pi <- function(k, alpha) {
  k <- as.matrix(k)
  storage.mode(k) <- "double"
  if (ncol(k) != length(alpha)) stop("shrink_pi: alpha length does not match categories", call. = FALSE)
  out <- (k + rep(alpha, each = nrow(k))) / (rowSums(k) + sum(alpha))
  colnames(out) <- c("pi_minus", "pi_mid", "pi_plus")
  out
}

# --- training-set priors ----------------------------------------------------

#' The measured subset of TRAIN, with fold and outcome attached. Written ONCE.
#'
#' Every quantity fitted on training rows -- alpha, p_bar, and the magnitude
#' deltas -- is defined over exactly this set. The original argument for putting
#' alpha and p_bar in one function was that "if that scoping were written twice
#' the two could disagree without anything failing"; a third estimator makes
#' that argument stronger, not weaker, so the scoping moved here rather than
#' being copied a third time.
#'
#' Test never enters. `n_obs > 0` is the measured-subset policy (spec SS5.5).
#'
#' @return signal_features rows, restricted, plus `fold` and `y`
#' The measured subset of a cohort, with the outcome attached. Site-agnostic.
#'
#' ONE COLUMN LIST, USED BY BOTH SCOPINGS. `.measured_train_rows()` below is
#' this function plus the train restriction and the fold column, and it is
#' written that way on purpose: a second copy of `need`/`opt` is exactly the
#' shape of the defect audit D exists to find. When the `se_` replicate
#' statistics were added to the extraction they reached the parquet, passed
#' validation, and never reached the estimator, because one hard-coded list
#' named them and another did not.
#'
#' @param stay_ids restrict to these stays, or NULL for the whole table
#' @return signal_features rows with `n_obs > 0`, plus `y`
measured_rows <- function(tabs, cfg, stay_ids = NULL) {
  sf <- tabs$signal_features
  need <- c("stay_id", "signal", "n_obs", "k_low", "k_mid", "k_high",
            "q05", "q95", "value_min", "value_max")
  miss <- setdiff(need, names(sf))
  if (length(miss)) abort_values("signal_features missing columns needed for priors", miss)

  # The replicate statistics for the `delta` variance components, carried
  # through only when the extraction produced them. OPTIONAL BY DESIGN: a
  # feature table built before 2026-09-02 has none, and `.var_components()`
  # falls back to the old heteroscedasticity estimator rather than failing.
  # They are NOT in `need`, because a hard requirement here would make an older
  # parquet unloadable for no benefit.
  opt <- intersect(c("se_within_n", "se_within_mean", "se_within_ss",
                     "se_hour_range_mean", "se_n_raw"), names(sf))

  keep <- sf$n_obs > 0
  if (!is.null(stay_ids)) keep <- keep & sf$stay_id %in% stay_ids
  s <- sf[keep, c(need, opt), drop = FALSE]
  s$y <- tabs$cohort$mortality[match(s$stay_id, tabs$cohort$stay_id)]
  if (anyNA(s$y)) stop("measured_rows: a signal row has no cohort match", call. = FALSE)
  s
}

.measured_train_rows <- function(tabs, folds, cfg) {
  train_ids <- folds$stay_id[folds$split == "train"]
  s <- measured_rows(tabs, cfg, stay_ids = train_ids)
  s$fold <- folds$fold[match(s$stay_id, folds$stay_id)]
  if (anyNA(s$fold)) stop(".measured_train_rows: a train row has no fold", call. = FALSE)
  # COLUMN ORDER IS PRESERVED (`..., fold, y`) exactly as it was before this
  # function was split, so that the refactor cannot move a number through some
  # positional access downstream. Verified against the stored `priors` target.
  s[, c(setdiff(names(s), c("fold", "y")), "fold", "y"), drop = FALSE]
}

#' The (signal, fold) job grid every per-signal prior is fitted over.
#'
#' `fold` names the fold the quantity is APPLIED to; it is fitted on the other
#' four. `NA` is the final fit, on all of train.
#' @param roles which roles to enumerate. The default is both, which is the
#'   primary graph. `"final"` alone is what the nested cross-fit asks for
#'   (`layer1_nested_l()`): its fold table already excludes the two held-out
#'   folds through `split`, so the fold rows would be fitted and never read.
.prior_jobs <- function(cfg, key = "signal", values = NULL,
                        roles = c("oof", "final")) {
  roles <- match.arg(roles, several.ok = TRUE)
  n_folds <- cfg$n_folds %||% 5L
  v <- values %||% unlist(cfg$signals)
  out <- rbind(
    if ("oof" %in% roles)
      expand.grid(key = v, fold = seq_len(n_folds),
                  stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE),
    if ("final" %in% roles)
      data.frame(key = v, fold = NA_integer_, stringsAsFactors = FALSE))
  names(out)[1] <- key
  out
}

#' Fit alpha and p_bar for every signal, on training rows only.
#'
#' Two roles, matching the stage table in CLAUDE.md:
#'
#'   role = "oof"    one row per (signal, fold). The quantities for fold f are
#'                   fitted on the OTHER four training folds and are used to
#'                   score fold f. `fold` names the fold they are APPLIED to,
#'                   never the folds they were fitted on.
#'   role = "final"  one row per signal, `fold` NA. Fitted on the whole training
#'                   set; these are what the bundle carries and what test and
#'                   eICU are scored with.
#'
#' p_bar is the mortality rate of the MEASURED subpopulation on the fitting
#' rows, not the cohort rate (CLAUDE.md, frozen decisions). Patients who get a
#' lactate drawn are sicker than average, and centering L on the cohort prior
#' would hand every measured patient a constant acuity offset unrelated to the
#' measured value. It is computed here, beside alpha, because both are defined
#' by the same row scoping and computing them in two places would let them
#' silently diverge.
#'
#' @param tabs  from load_tables()
#' @param folds from assign_folds()
#' @param cfg   from load_config()
#' @return data frame: signal, role, fold, n_stays, n_deaths, p_bar,
#'         alpha_low, alpha_mid, alpha_high, alpha0, iter, converged, degenerate
training_priors <- function(tabs, folds, cfg, verbose = TRUE, rows_ = NULL,
                            roles = c("oof", "final")) {
  s <- rows_ %||% .measured_train_rows(tabs, folds, cfg)
  jobs <- .prior_jobs(cfg, roles = roles)

  rows <- lapply(seq_len(nrow(jobs)), function(i) {
    sg <- jobs$signal[i]; fd <- jobs$fold[i]
    # fold f's alpha is fitted on the COMPLEMENT of fold f; the final alpha on
    # all of train.
    sel <- s$signal == sg & (if (is.na(fd)) TRUE else s$fold != fd)
    z <- s[sel, , drop = FALSE]
    if (!nrow(z)) stop("training_priors: no measured training rows for signal '", sg, "'", call. = FALSE)
    fit <- fit_dm_alpha(cbind(z$k_low, z$k_mid, z$k_high))
    data.frame(
      signal     = sg,
      role       = if (is.na(fd)) "final" else "oof",
      fold       = fd,
      n_stays    = nrow(z),
      n_deaths   = sum(z$y),
      p_bar      = mean(z$y),
      alpha_low  = fit$alpha[1], alpha_mid = fit$alpha[2], alpha_high = fit$alpha[3],
      alpha0     = fit$alpha0,
      iter       = fit$iter,
      converged  = fit$converged,
      degenerate = any(fit$degenerate),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  if (verbose) {
    bad <- !out$converged
    message(sprintf("training_priors: %d fits, %d converged%s",
                    nrow(out), sum(out$converged),
                    if (any(bad)) sprintf("; NOT converged: %s",
                                          paste(unique(out$signal[bad]), collapse = ", ")) else ""))
    if (any(out$degenerate)) {
      message("  degenerate alpha (a category never observed) for: ",
              paste(unique(out$signal[out$degenerate]), collapse = ", "))
    }
  }
  out
}

#' Fit the magnitude (delta) parameters for every (signal, magnitude variable).
#'
#' Same rows as training_priors(), by construction -- both take `rows_` from
#' .measured_train_rows(). Long rather than wide because a signal carries one
#' magnitude variable when paired and two when not, so a wide frame would be
#' half NA and the arity would be implicit.
#'
#' k = 0 STAYS IN THE FIT. See R/04b_conditional.R: for a signal, k = 0 means
#' "measured, never crossed the bound", which is a real observation with a real
#' magnitude. Only n_obs = 0 is uninformative, and that is already excluded.
#'
#' @return data frame: signal, variable, count_var, role, fold, n_stays,
#'         a0, a1, a2, s_e, s_u, converged, degenerate
magnitude_priors <- function(tabs, folds, cfg, verbose = TRUE, rows_ = NULL,
                             roles = c("oof", "final")) {
  vars <- lapply(cfg$signals, function(sg)
    if (magnitude_conditional_for(sg, cfg)) level_vars_of(sg, cfg) else character(0))
  names(vars) <- cfg$signals
  if (!length(unlist(vars))) return(.empty_magnitude_priors())

  s <- rows_ %||% .measured_train_rows(tabs, folds, cfg)
  jobs <- .prior_jobs(cfg, values = names(vars)[lengths(vars) > 0], roles = roles)

  rows <- list()
  for (i in seq_len(nrow(jobs))) {
    sg <- jobs$signal[i]; fd <- jobs$fold[i]
    sel <- s$signal == sg & (if (is.na(fd)) TRUE else s$fold != fd)
    z <- s[sel, , drop = FALSE]
    if (!nrow(z)) stop("magnitude_priors: no measured training rows for '", sg, "'", call. = FALSE)
    # DECLARED per signal, resolved once per job. The form and the scale travel
    # into the stored row so delta_value() can route without consulting config
    # at an apply site (hard rule 8).
    fm <- delta_form_of(sg, cfg)
    sc <- ordinal_scale_of(sg, cfg)
    for (v in vars[[sg]]) {
      kv <- delta_count_of(v)
      # The replicate statistics travel alongside the value. NULL when the
      # feature table predates them (2026-09-02), which `.var_components()`
      # detects and falls back on rather than failing -- an older parquet must
      # still reproduce its old numbers.
      p  <- fit_delta(z[[v]], z[[kv]], z$n_obs, form = fm, scale = sc,
                      within_ss = z[["se_within_ss"]],
                      within_n  = z[["se_within_n"]])
      rows[[length(rows) + 1L]] <- data.frame(
        signal = sg, variable = v, count_var = kv,
        role = if (is.na(fd)) "final" else "oof", fold = fd,
        form = p$form,
        n_stays = p$n_stays, a0 = p$a0, a1 = p$a1, a2 = p$a2, a3 = p$a3,
        theta = .pack_num(p$theta), levels = .pack_num(p$levels),
        loglik = p$loglik %||% NA_real_,
        s_e = p$s_e, s_u = p$s_u, shrinks = p$shrinks,
        se_source = p$se_source %||% NA_character_, se_df = p$se_df %||% 0,
        converged = p$converged, degenerate = p$degenerate,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (verbose) {
    message(sprintf("magnitude_priors: %d fits over %d signal(s), %d converged%s",
                    nrow(out), length(unique(out$signal)), sum(out$converged),
                    if (any(out$degenerate)) sprintf("; DEGENERATE: %s",
                      paste(unique(out$signal[out$degenerate]), collapse = ", ")) else ""))
  }
  out
}

.empty_magnitude_priors <- function() {
  data.frame(signal = character(0), variable = character(0),
             count_var = character(0), role = character(0), fold = integer(0),
             form = character(0), n_stays = integer(0),
             a0 = numeric(0), a1 = numeric(0), a2 = numeric(0), a3 = numeric(0),
             theta = character(0), levels = character(0), loglik = numeric(0),
             s_e = numeric(0), s_u = numeric(0), shrinks = logical(0),
             se_source = character(0), se_df = numeric(0),
             converged = logical(0), degenerate = logical(0),
             stringsAsFactors = FALSE)
}

#' Fit the intensity (lambda) parameters for every intervention that has two
#' intensity covariates.
#'
#' KEYED ON INTERVENTION, NOT ON SIGNAL, and that is a design constraint rather
#' than a convenience. R/05_formula.R's contract for the `intv` block is that
#' "every term here is a property of the treatment record ALONE", which is what
#' makes L_intv a clean estimate of log p(I|Y=1)/p(I|Y=0). If lambda were fitted
#' on each signal's measured subset, `mbp` and `heart_rate` would receive
#' different lambda columns for the same intervention and measurement
#' information would have leaked into the intervention block. The three GCS
#' models must get BYTE-IDENTICAL lambda columns even though their L_intv
#' differ, and check_lambda_invariance() asserts exactly that.
#'
#' Fitted on the EXPOSED subset. Unlike the magnitude case, an unexposed stay
#' carries no intensity information at all -- and for the lognormal family
#' log(0) is not even defined. Unexposed stays are assigned lambda = 0 at apply
#' time, never fitted.
#'
#' @return data frame: intervention, family, role, fold, exposure_var,
#'         accumulation_var, M, n_stays, c0/c1/phi (bb) or b0/b1/s_e/s_u (ln),
#'         converged, degenerate
intervention_priors <- function(tabs, folds, cfg, verbose = TRUE,
                                roles = c("oof", "final")) {
  ivs <- Filter(function(iv) !is.null(lambda_spec_of(iv, cfg)),
                unique(unlist(lapply(cfg$signals, interventions_of, cfg = cfg))))
  if (!length(ivs)) return(.empty_intervention_priors())

  ivf <- tabs$intervention_features
  train_ids <- folds$stay_id[folds$split == "train"]
  fold_of <- folds$fold[match(ivf$stay_id, folds$stay_id)]
  keep <- ivf$stay_id %in% train_ids & ivf$intervention %in% ivs
  z0 <- ivf[keep, , drop = FALSE]
  z0$fold <- fold_of[keep]
  if (anyNA(z0$fold)) stop("intervention_priors: a train row has no fold", call. = FALSE)

  jobs <- .prior_jobs(cfg, key = "intervention", values = ivs, roles = roles)
  rows <- list()
  for (i in seq_len(nrow(jobs))) {
    iv <- jobs$intervention[i]; fd <- jobs$fold[i]
    ls <- lambda_spec_of(iv, cfg)
    sel <- z0$intervention == iv & (if (is.na(fd)) TRUE else z0$fold != fd)
    z <- z0[sel, , drop = FALSE]
    if (!nrow(z)) stop("intervention_priors: no training rows for '", iv, "'", call. = FALSE)

    D <- z[[ls$exposure]] * ls$exposure_scale
    A <- z[[ls$accumulation]]
    # Shape-inapplicable columns are NA by contract (spec SS5.3); the shape map
    # says these two apply, so an NA here means the map and the data disagree.
    if (anyNA(D) || anyNA(A)) {
      stop(sprintf("intervention_priors [%s]: NA in %s/%s, which intervention_shape ",
                   "says are applicable. The shape map and the extraction disagree.",
                   iv, ls$exposure, ls$accumulation), call. = FALSE)
    }

    p <- if (identical(ls$family, "binomial")) fit_lambda_binom(A, D, ls$M)
         else fit_lambda_ln(A, D)

    rows[[length(rows) + 1L]] <- data.frame(
      intervention = iv, family = ls$family,
      role = if (is.na(fd)) "final" else "oof", fold = fd,
      exposure_var = ls$exposure, accumulation_var = ls$accumulation,
      exposure_scale = ls$exposure_scale,
      M  = if (identical(ls$family, "binomial")) as.integer(p$M) else NA_integer_,
      n_stays = p$n_stays,
      c0 = p$c0 %||% NA_real_, c1 = p$c1 %||% NA_real_,
      b0 = p$b0 %||% NA_real_, b1 = p$b1 %||% NA_real_,
      s_e = p$s_e %||% NA_real_, s_u = p$s_u %||% NA_real_,
      shrinks = p$shrinks %||% FALSE,
      converged = p$converged, degenerate = p$degenerate,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (verbose) {
    message(sprintf("intervention_priors: %d fits over %d intervention(s), %d converged%s",
                    nrow(out), length(ivs), sum(out$converged),
                    if (any(out$degenerate)) sprintf("; DEGENERATE: %s",
                      paste(unique(out$intervention[out$degenerate]), collapse = ", ")) else ""))
  }
  out
}

.empty_intervention_priors <- function() {
  data.frame(intervention = character(0), family = character(0), role = character(0),
             fold = integer(0), exposure_var = character(0),
             accumulation_var = character(0), exposure_scale = numeric(0),
             M = integer(0), n_stays = integer(0),
             c0 = numeric(0), c1 = numeric(0),
             b0 = numeric(0), b1 = numeric(0), s_e = numeric(0), s_u = numeric(0),
             shrinks = logical(0), converged = logical(0), degenerate = logical(0),
             stringsAsFactors = FALSE)
}

# --- the priors container ---------------------------------------------------

#' Every frozen quantity layer 1 needs, fitted on training rows only.
#'
#' THE SINGLE OBJECT THAT GETS THREADED. Three tables at three grains --
#' per (signal, fold), per (signal, variable, fold), per (intervention, fold) --
#' assembled once and passed whole.
#'
#' This replaced three separate arguments (`alpha`, `p_bar`, and what would have
#' been two more) for a reason specific to this pipeline: a call site that
#' passes the wrong one of four parameter sets produces a covariate standardised
#' against the wrong population, a fit that converges, and an L that looks
#' plausible. That is the failure class the project is designed against, and an
#' argument list cannot rule it out while one object can.
#'
#' The row scoping is walked ONCE and shared, so alpha, p_bar and the deltas
#' cannot silently come to disagree about which rows they were fitted on.
#'
#' THE FIT-STATUS GATE RUNS HERE, as of 2026-09-08 (plumbing review F3), so
#' that no caller -- the graph, a bootstrap replicate, a test -- can carry an
#' unsuccessful prior fit into a GAM without saying so. `strict = FALSE` is for
#' triage only. The policy itself is `check_prior_fits()` below.
#'
#' @return object of class `llr_priors`: list(signal, magnitude, intervention)
#' @param roles passed to `.prior_jobs()`. The primary graph fits both; the
#'   nested cross-fit fits `"final"` only, on a fold table whose `split` column
#'   already excludes the folds being predicted.
layer1_priors <- function(tabs, folds, cfg, verbose = TRUE, strict = TRUE,
                          roles = c("oof", "final")) {
  rows_ <- .measured_train_rows(tabs, folds, cfg)
  out <- structure(list(
    signal       = training_priors(tabs, folds, cfg, verbose = verbose, rows_ = rows_,
                                   roles = roles),
    magnitude    = magnitude_priors(tabs, folds, cfg, verbose = verbose, rows_ = rows_,
                                    roles = roles),
    intervention = intervention_priors(tabs, folds, cfg, verbose = verbose,
                                       roles = roles)
  ), class = "llr_priors")
  check_prior_fits(out, cfg, strict = strict, verbose = verbose)
  out
}

#' The prior-fit STATUS POLICY: which fit outcomes may proceed to a GAM.
#'
#' Every estimator behind `layer1_priors()` returns its last parameters with a
#' flag rather than failing -- `fit_dm_alpha()` at `max_iter`, `.polr_fit()`
#' after both optimisers fail, `fit_lambda_binom()` on a constant count -- and
#' until 2026-09-08 the flags were RECORDED and READ BY NOTHING on the fitting
#' path: `training_priors()` described a non-convergence only when `verbose`,
#' the graph calls it with `verbose = FALSE`, and `verify_bundle()` accepts any
#' finite positive alpha. A Dirichlet fit that stopped at 10,000 iterations
#' would therefore have fitted 384 GAMs and built a bundle that passed every
#' check (plumbing review F3).
#'
#' THIS IS THE ONE PLACE THE POLICY IS STATED. It is not a config key on
#' purpose: adding one to config/config.yml would change `cfg` and rebuild
#' every fit to produce identical numbers (the argument CLAUDE.md makes for
#' config/attribution.yml), and no value of it changes what is fitted -- it
#' only decides whether a fit that already happened may be USED.
#'
#' The rules, block by block. FAIL is fatal under `strict`; WARN is reported.
#'
#'   Dirichlet (`signal`)
#'     converged == FALSE                        FAIL. The MM iteration is
#'                                               monotone and 133-687 steps
#'                                               suffice at tol 1e-8; 10,000
#'                                               without convergence is a
#'                                               numerical failure, not a slow
#'                                               fit.
#'     a coordinate at the estimator floor       FAIL where the tail is declared
#'                                               OCCUPIABLE in `signal_tails`,
#'                                               PASS (expected) where it is
#'                                               declared empty. The four
#'                                               declared-empty high tails --
#'                                               spo2 and the Glasgow triple --
#'                                               sit at the floor BY DESIGN and
#'                                               enter no formula (CLAUDE.md,
#'                                               `signal_tails`). The MID
#'                                               coordinate is never allowed at
#'                                               the floor. Tested per
#'                                               coordinate from the stored
#'                                               alphas, not from the collapsed
#'                                               `degenerate` flag, because
#'                                               `any()` over three coordinates
#'                                               lets an expected high-tail floor
#'                                               mask an unexpected low-tail one.
#'   delta (`magnitude`)
#'     degenerate == TRUE                        FAIL. `delta_value()` returns a
#'                                               zero covariate for it, and no
#'                                               signal declares a degenerate
#'                                               magnitude
#'                                               (`magnitude_conditional_override`
#'                                               is empty, CLAUDE.md).
#'     form == "ordinal" & converged == FALSE    FAIL. The proportional-odds fit
#'                                               IS the design path. The stored
#'                                               flag conflates it with the
#'                                               variance components (`f$converged
#'                                               && vc$converged`) and a false
#'                                               stop is the cheaper error;
#'                                               splitting the flag adds a column
#'                                               to `priors` and costs a full
#'                                               rebuild, so it waits for one.
#'     form == "linear"  & converged == FALSE    WARN. The flag is the variance-
#'                                               component fit alone, and since
#'                                               2026-09-05 `s_e`/`s_u` reach no
#'                                               covariate (the shrinkage is
#'                                               removed); the OLS conditional
#'                                               mean has no failure mode.
#'   lambda (`intervention`)
#'     degenerate == TRUE                        FAIL. `lambda_value()` returns 0
#'                                               for every stay.
#'     family == "binomial"  & !converged        FAIL. The GLM is the design path.
#'     family == "lognormal" & !converged        WARN. Variance components only;
#'                                               `b0`/`b1` come from `lm()`.
#'
#' @return data frame, one row per (block, key, role, fold): status
#'   PASS / WARN / FAIL and a reason. Invisible under `strict`.
check_prior_fits <- function(priors, cfg, strict = TRUE, verbose = FALSE) {
  if (!inherits(priors, "llr_priors")) {
    stop("check_prior_fits: not an llr_priors container", call. = FALSE)
  }
  rows <- list()
  add <- function(block, key, role, fold, status, reason) {
    rows[[length(rows) + 1L]] <<- data.frame(
      block = block, key = key, role = role, fold = fold,
      status = status, reason = reason, stringsAsFactors = FALSE)
  }
  # The same multiple `fit_dm_alpha()`'s own `degenerate` test uses.
  floor_at <- DM_ALPHA_FLOOR * 10

  sp <- priors$signal
  for (i in seq_len(nrow(sp))) {
    sg <- sp$signal[i]
    declared <- occupiable_tails_of(sg, cfg)
    fl <- c(low = sp$alpha_low[i], mid = sp$alpha_mid[i], high = sp$alpha_high[i]) <= floor_at
    allowed <- c(low = !("low" %in% declared), mid = FALSE, high = !("high" %in% declared))
    bad <- names(fl)[fl & !allowed]
    if (!isTRUE(sp$converged[i])) {
      add("signal", sg, sp$role[i], sp$fold[i], "FAIL",
          sprintf("Dirichlet fit did not converge (%d iterations)", sp$iter[i]))
    } else if (length(bad)) {
      add("signal", sg, sp$role[i], sp$fold[i], "FAIL",
          paste0("alpha at the estimator floor on a declared-occupiable coordinate: ",
                 paste(bad, collapse = ", ")))
    } else if (any(fl)) {
      add("signal", sg, sp$role[i], sp$fold[i], "PASS",
          paste0("alpha at the floor on the declared-empty tail: ",
                 paste(names(fl)[fl], collapse = ", ")))
    } else {
      add("signal", sg, sp$role[i], sp$fold[i], "PASS", "")
    }
  }

  mg <- priors$magnitude
  for (i in seq_len(nrow(mg))) {
    key <- paste(mg$signal[i], mg$variable[i])
    fm  <- mg$form[i] %||% "linear"
    if (isTRUE(mg$degenerate[i])) {
      add("magnitude", key, mg$role[i], mg$fold[i], "FAIL",
          "degenerate delta fit; delta_value() would return a zero covariate")
    } else if (!isTRUE(mg$converged[i])) {
      if (identical(fm, "ordinal")) {
        add("magnitude", key, mg$role[i], mg$fold[i], "FAIL",
            "ordinal delta did not converge (proportional-odds fit or its variance components)")
      } else {
        add("magnitude", key, mg$role[i], mg$fold[i], "WARN",
            "variance components did not converge; s_e/s_u are diagnostic only since 2026-09-05")
      }
    } else {
      add("magnitude", key, mg$role[i], mg$fold[i], "PASS", "")
    }
  }

  ip <- priors$intervention
  for (i in seq_len(nrow(ip))) {
    iv <- ip$intervention[i]
    if (isTRUE(ip$degenerate[i])) {
      add("intervention", iv, ip$role[i], ip$fold[i], "FAIL",
          "degenerate lambda fit; lambda_value() would return 0 for every stay")
    } else if (!isTRUE(ip$converged[i])) {
      if (identical(ip$family[i], "binomial")) {
        add("intervention", iv, ip$role[i], ip$fold[i], "FAIL",
            "binomial lambda GLM did not converge")
      } else {
        add("intervention", iv, ip$role[i], ip$fold[i], "WARN",
            "variance components did not converge; s_e/s_u are diagnostic only since 2026-09-05")
      }
    } else {
      add("intervention", iv, ip$role[i], ip$fold[i], "PASS", "")
    }
  }

  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(block = character(0), key = character(0), role = character(0),
               fold = integer(0), status = character(0), reason = character(0),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL

  n_fail <- sum(out$status == "FAIL"); n_warn <- sum(out$status == "WARN")
  if (verbose || n_fail || n_warn) {
    message(sprintf("check_prior_fits: %d fit(s), %d FAIL, %d WARN", nrow(out), n_fail, n_warn))
  }
  if (n_warn) {
    w <- out[out$status == "WARN", , drop = FALSE]
    message(paste0("  [WARN] ", w$block, " ", w$key, " (", w$role,
                   ifelse(is.na(w$fold), "", paste0(" fold ", w$fold)), "): ", w$reason,
                   collapse = "\n"))
  }
  if (n_fail) {
    f <- out[out$status == "FAIL", , drop = FALSE]
    msg <- paste0("check_prior_fits: ", n_fail, " prior fit(s) failed the status policy.\n",
                  paste0("  [FAIL] ", f$block, " ", f$key, " (", f$role,
                         ifelse(is.na(f$fold), "", paste0(" fold ", f$fold)), "): ", f$reason,
                         collapse = "\n"),
                  "\nAn unsuccessful prior fit must not reach a GAM: the covariate it ",
                  "standardises would be wrong, the GAM would converge, and the L would ",
                  "look plausible. Fix the fit or declare the degeneracy; do not relax this.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  if (strict) invisible(out) else out
}

#' @export
print.llr_priors <- function(x, ...) {
  cat(sprintf("<llr_priors> signal: %d rows | magnitude: %d rows | intervention: %d rows\n",
              nrow(x$signal), nrow(x$magnitude), nrow(x$intervention)))
  invisible(x)
}

#' Accept either the container or the bare signal table.
#'
#' The low-level row accessors (priors_of, alpha_of) predate the container and
#' are still the right tool when only alpha is wanted. They take both so that a
#' diagnostic script asking one narrow question does not have to build all three
#' tables. priors_for() requires the container, because a partial answer there
#' would be exactly the silent-mismatch failure the container exists to prevent.
.sig_tbl <- function(p) if (inherits(p, "llr_priors")) p$signal else p

#' Everything one (signal, role, fold) needs, in one object.
#'
#' @return list(signal, role, fold, alpha, p_bar, magnitude, lambda)
#'   `magnitude` is named by magnitude variable, `lambda` by intervention. Both
#'   are empty lists when the corresponding construct is off, which is what the
#'   frame builders test -- they never read the config flag themselves, so the
#'   formula and the frame cannot disagree about whether a construct is active.
priors_for <- function(priors, signal, role = c("oof", "final"), fold = NA_integer_) {
  role <- match.arg(role)
  if (!inherits(priors, "llr_priors")) {
    stop("priors_for: expected the container from layer1_priors(), got a bare ",
         class(priors)[1], ". Threading the signal table alone would silently ",
         "drop the magnitude and intensity parameters.", call. = FALSE)
  }
  r <- priors_of(priors$signal, signal, role, fold)

  mg <- priors$magnitude
  sel <- mg$signal == signal & mg$role == role &
    (if (is.na(fold)) is.na(mg$fold) else !is.na(mg$fold) & mg$fold == fold)
  mag <- lapply(split(mg[sel, , drop = FALSE], mg$variable[sel]), as.list)

  ip <- priors$intervention
  sel <- ip$role == role & (if (is.na(fold)) is.na(ip$fold) else !is.na(ip$fold) & ip$fold == fold)
  lam <- lapply(split(ip[sel, , drop = FALSE], ip$intervention[sel]), as.list)

  list(signal = signal, role = role, fold = fold,
       alpha = c(r$alpha_low, r$alpha_mid, r$alpha_high),
       p_bar = r$p_bar, magnitude = mag, lambda = lam)
}

#' Hold the intervention block's purity to the data.
#'
#' `intv` terms must be properties of the treatment record ALONE. lambda is
#' keyed on intervention, so two signals sharing an intervention must receive
#' identical columns. Cheap, and it fails loudly if intervention_priors() is
#' ever re-scoped per signal.
check_lambda_invariance <- function(tabs, cfg, priors, stay_ids = NULL,
                                    role = "final", fold = NA_integer_, strict = TRUE) {
  shared <- list()
  for (sg in cfg$signals) for (iv in interventions_of(sg, cfg)) {
    if (is.null(lambda_spec_of(iv, cfg))) next
    shared[[iv]] <- c(shared[[iv]], sg)
  }
  shared <- shared[lengths(shared) > 1L]
  if (!length(shared)) return(invisible(NULL))

  bad <- character(0)
  for (iv in names(shared)) {
    cols <- lapply(shared[[iv]], function(sg) {
      pri <- priors_for(priors, sg, role, fold)
      d <- signal_frame(sg, "intv", tabs, cfg, pri, stay_ids = stay_ids)
      stats::setNames(d[[paste0(iv, "__lambda")]], d$stay_id)
    })
    ids <- Reduce(intersect, lapply(cols, names))
    ref <- cols[[1]][ids]
    for (j in seq_along(cols)[-1]) {
      if (!isTRUE(all.equal(unname(ref), unname(cols[[j]][ids])))) bad <- c(bad, iv)
    }
  }
  if (length(bad)) {
    msg <- paste0("lambda differs between signals sharing an intervention: ",
                  paste(unique(bad), collapse = ", "),
                  ". The intv block is no longer a property of the treatment ",
                  "record alone, and L_intv has stopped being log p(I|Y)/p(I|Y').")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' Is `excursion_side` pointing at a side where excursions actually happen?
#'
#' A sanity check on pairing.csv, not on the model: `excursion_side` selects the
#' tail quantile (q05 vs q95), so a side named backwards would
#' pair vasopressor to hypertension and still fit cleanly. The assertion is only
#' that the named side is non-empty — deliberately weak, because dominance is the
#' wrong test: urine output is paired on `low` although high-output hours are
#' more common, and oliguria is still the concern.
#'
#' Which pi coordinates enter is a SEPARATE question, governed by `signal_tails`
#' and checked by check_signal_tails(). Occupancy in the unnamed tail is normal
#' and is reported, not flagged.
#'
#' `n_thin` is the stay count in the named side. Occupied-but-barely is worth
#' seeing before the fit rather than inferring it from an edf of 0 afterwards:
#' lactate has 53 stays below its reference low, of 27,284 measured.
check_excursion_sides <- function(tabs, cfg, strict = TRUE) {
  sf <- tabs$signal_features
  m <- sf$n_obs > 0
  rows <- lapply(cfg$signals, function(sg) {
    z <- sf[m & sf$signal == sg, c("k_low", "k_high"), drop = FALSE]
    side <- excursion_side_of(sg, cfg)
    modelled <- if (is.na(side)) c("low", "high") else side
    stays <- c(low = sum(z$k_low > 0), high = sum(z$k_high > 0))
    data.frame(signal = sg,
               side = side %||% NA_character_,
               modelled = paste(modelled, collapse = "+"),
               stays_low = stays[["low"]], stays_high = stays[["high"]],
               n_thin = min(stays[modelled]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$ok <- out$n_thin > 0L
  if (any(!out$ok)) {
    bad <- out$signal[!out$ok]
    msg <- paste0("the modelled excursion tail is empty for: ", paste(bad, collapse = ", "),
                  ". pi would carry no information. Fix excursion_side in ",
                  "config/pairing.csv or the reference range — not here.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  out
}

#' Hold `signal_tails` to the counts, in BOTH directions.
#'
#'   declared occupiable but empty  - alpha_j falls to the estimator floor and
#'                                    pi_j becomes alpha_j/(n_obs + alpha_0), a
#'                                    monotone function of n_obs and nothing
#'                                    else. The formula would then carry
#'                                    measurement frequency under a physiology
#'                                    label — the one term deliberately excluded
#'                                    from every model.
#'   declared empty but occupied    - a free simplex coordinate is being
#'                                    discarded. Two patients with the same
#'                                    k_low and n_obs but different k_high have
#'                                    identical pi_minus, so this silently
#'                                    erases the labile-versus-one-sided
#'                                    distinction.
#'
#' Neither announces itself downstream: both produce a fit that succeeds and L's
#' that look plausible. Same pattern as the validator's `intervention_shape`
#' check — a static design map, re-checked against the data at both sites.
check_signal_tails <- function(tabs, cfg, strict = TRUE) {
  sf <- tabs$signal_features
  m <- sf$n_obs > 0
  rows <- lapply(cfg$signals, function(sg) {
    z <- sf[m & sf$signal == sg, c("k_low", "k_high"), drop = FALSE]
    declared <- occupiable_tails_of(sg, cfg)
    data.frame(signal = sg,
               declared = paste(declared, collapse = "+"),
               declared_low  = "low"  %in% declared,
               declared_high = "high" %in% declared,
               stays_low  = sum(z$k_low  > 0),
               stays_high = sum(z$k_high > 0),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$ok <- (out$declared_low  == (out$stays_low  > 0)) &
            (out$declared_high == (out$stays_high > 0))
  if (any(!out$ok)) {
    msg <- paste0("signal_tails disagrees with the counts for: ",
                  paste(out$signal[!out$ok], collapse = ", "),
                  ". Update config/signal_tails, or re-extract — do not work around it here.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  out
}

#' Hold `config/smooth_k` to the measured distinct-value counts.
#'
#' Checked in BOTH directions, as check_signal_tails() is, and the second is the
#' one that matters:
#'
#'   declared but unnecessary  the covariate comfortably supports the default k,
#'                             so the override is silently costing basis
#'                             flexibility for no reason.
#'   necessary but undeclared  bam() will fail outright at that fit. Better here,
#'                             in seconds, than at fit 173 of 258.
#'
#' Run on the population it will be fitted on. The default is full train; pass
#' a fold's fitting ids to check that subset, which is the stricter test.
#'
#' @return data frame: signal, variable, k, n_unique, declared, needed, ok
check_smooth_k <- function(tabs, cfg, priors, stay_ids = NULL, strict = TRUE) {
  kdef <- cfg$bam$k %||% 10
  rows <- list()
  for (sg in cfg$signals) {
    pri <- priors_for(priors, sg, "final", NA_integer_)
    seen <- character(0)
    for (md in models_of(sg, cfg)) {
      f <- build_formula(sg, md, cfg)
      sp <- smooth_specs(f)
      sp <- sp[!sp$variable %in% seen, , drop = FALSE]
      if (!nrow(sp)) next
      seen <- c(seen, sp$variable)
      d <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = stay_ids)
      for (i in seq_len(nrow(sp))) {
        v <- sp$variable[i]
        nu <- length(unique(d[[v]]))
        rows[[length(rows) + 1L]] <- data.frame(
          signal = sg, variable = v, k = sp$k[i], n_unique = nu,
          declared = !is.null(cfg$smooth_k[[sg]][[v]]),
          needed   = nu <= kdef,
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  out$ok <- out$declared == out$needed & out$k < out$n_unique
  rownames(out) <- NULL

  if (any(!out$ok)) {
    b <- out[!out$ok, , drop = FALSE]
    over  <- b[b$declared & !b$needed, , drop = FALSE]
    under <- b[!b$declared & b$needed, , drop = FALSE]
    infeas <- b[b$k >= b$n_unique, , drop = FALSE]
    msg <- paste0("config/smooth_k disagrees with the counts. ", paste(c(
      if (nrow(infeas)) sprintf("k >= distinct values (bam WILL fail): %s",
        paste(sprintf("%s/%s k=%d vs %d", infeas$signal, infeas$variable,
                      infeas$k, infeas$n_unique), collapse = ", ")),
      if (nrow(under)) sprintf("undeclared but needs one: %s",
        paste(sprintf("%s/%s (%d distinct)", under$signal, under$variable,
                      under$n_unique), collapse = ", ")),
      if (nrow(over)) sprintf("declared but unnecessary: %s",
        paste(sprintf("%s/%s (%d distinct)", over$signal, over$variable,
                      over$n_unique), collapse = ", "))), collapse = "; "))
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  out
}

#' Look up one row of the priors table. Errors rather than returning NULL, so a
#' wrong (role, fold) cannot silently produce an unshrunk or mis-shrunk pi.
priors_of <- function(priors, signal, role = c("oof", "final"), fold = NA_integer_) {
  role <- match.arg(role)
  # Normalised FIRST. On the container `priors$signal` is a data frame rather
  # than a column, so the comparison below would silently match nothing instead
  # of failing -- the one place this accessor could go quietly wrong.
  priors <- .sig_tbl(priors)
  sel <- priors$signal == signal & priors$role == role &
    (if (is.na(fold)) is.na(priors$fold) else !is.na(priors$fold) & priors$fold == fold)
  if (sum(sel) != 1L) {
    stop(sprintf("priors_of: expected exactly one row for (%s, %s, fold=%s), found %d",
                 signal, role, fold, sum(sel)), call. = FALSE)
  }
  priors[sel, , drop = FALSE]
}

#' The alpha triple from a priors row, in category order.
alpha_of <- function(priors, signal, role = c("oof", "final"), fold = NA_integer_) {
  r <- priors_of(priors, signal, role, fold)
  c(r$alpha_low, r$alpha_mid, r$alpha_high)
}

# --- the model frame --------------------------------------------------------

#' Build the layer-1 model frame for one (signal, model).
#'
#' The contract with R/05_formula.R runs in one direction: the formula is built
#' first, required_columns() is read off it, and this function produces exactly
#' those columns and no others. Nothing is assembled from a second list that
#' could drift — in particular the per-intervention columns are derived by
#' matching the `{intervention}__` prefix against the formula's own variables,
#' so a change to intervention_terms() propagates here with no edit.
#'
#' @param signal   one of cfg$signals
#' @param model    one of LAYER1_MODELS. `intv` exists only for paired signals
#'                 (see models_of()); on an unpaired one build_formula() refuses
#'                 it, because L_intv is 0 by construction there and assigned.
#' @param tabs     from load_tables()
#' @param cfg      from load_config()
#' @param alpha    length-3 Dirichlet alpha, from alpha_of(). Passed in, never
#'                 derived here: this function knows nothing about folds, so it
#'                 cannot accidentally shrink a held-out fold using its own data.
#' @param stay_ids restrict to these stays. NULL means every stay in the table.
#' @param measured_only keep only `n_obs > 0`. TRUE always in this pipeline; the
#'                 argument exists so the measured-subset policy is visible at
#'                 the call site rather than buried (spec §5.5).
#'
#'                 This applies to `intv` too, and that is a CONSTRAINT rather
#'                 than an inherited default. "p(mortality | interventions only)"
#'                 sounds cohort-wide and must not be: L_cond = L_full - L_intv
#'                 is only the conditional term if both models occupy the same
#'                 rows with the same p_bar. Fitted cohort-wide, `intv` would
#'                 turn the subtraction into a population contrast wearing the
#'                 conditional's name. Consequence worth stating: the three GCS
#'                 `intv` models carry identical terms but sit on three
#'                 different measured subsets, so their L_intv are
#'                 near-identical, not identical.
#' @return data frame: stay_id + exactly required_columns(formula). No NAs.
signal_frame <- function(signal, model = LAYER1_MODELS, tabs, cfg, priors,
                         stay_ids = NULL, measured_only = TRUE,
                         stage = c("fit", "predict")) {
  stage <- match.arg(stage)
  model <- match.arg(model)
  f   <- build_formula(signal, model, cfg)
  req <- required_columns(f)
  pri <- .as_row_priors(priors, signal)

  s <- .frame_rows(signal, tabs, stay_ids, measured_only)
  out <- .frame_response(s, tabs, f)
  out <- .frame_direct(out, s, req)
  out <- .frame_propensity(out, s, req, pri)
  out <- .frame_magnitude(out, s, req, pri)
  out <- .frame_intervention(out, s, req, signal, tabs, cfg, pri)

  check_model_frame(out, f, signal, model, stage = stage)
  attr(out, "signal")  <- signal
  attr(out, "model")   <- model
  attr(out, "formula") <- f
  out
}

#' Accept the container or an already-resolved row, never a bare alpha vector.
#'
#' The old signature took `alpha` as a length-3 numeric. That is rejected
#' explicitly rather than tolerated: a caller still passing one is a caller who
#' has not been told about the magnitude and intensity parameters, and silently
#' defaulting those to "absent" would build a frame missing columns the formula
#' asks for -- caught by check_model_frame, but with a confusing message.
.as_row_priors <- function(priors, signal) {
  if (is.list(priors) && !is.null(priors$alpha) && !is.null(priors$magnitude)) return(priors)
  if (is.numeric(priors)) {
    stop("signal_frame: `priors` is a bare numeric. The alpha-only signature was ",
         "replaced by priors_for(priors, signal, role, fold), which carries alpha, ",
         "p_bar, the magnitude parameters and the intensity parameters together.",
         call. = FALSE)
  }
  stop("signal_frame: `priors` must come from priors_for(); got ", class(priors)[1],
       call. = FALSE)
}

#' Row selection. The measured-subset policy lives here and nowhere else.
.frame_rows <- function(signal, tabs, stay_ids, measured_only) {
  sf <- tabs$signal_features
  sel <- sf$signal == signal
  if (!is.null(stay_ids)) sel <- sel & sf$stay_id %in% stay_ids
  if (measured_only) sel <- sel & sf$n_obs > 0
  s <- sf[sel, , drop = FALSE]
  if (!nrow(s)) stop("signal_frame: no rows for signal '", signal, "'", call. = FALSE)
  s
}

.frame_response <- function(s, tabs, f) {
  out <- data.frame(stay_id = s$stay_id)
  y <- all.vars(f)[1]
  out[[y]] <- tabs$cohort[[y]][match(s$stay_id, tabs$cohort$stay_id)]
  if (anyNA(out[[y]])) stop("signal_frame: a signal row has no cohort match", call. = FALSE)
  out
}

#' Columns carried straight through from signal_features.
#'
#' Derived names are excluded by construction rather than by a list: a `_delta`
#' column does not exist in the parquet, so `intersect(req, names(s))` never
#' matches one. The same is true of pi_minus / pi_plus.
.frame_direct <- function(out, s, req) {
  for (v in intersect(req, names(s))) out[[v]] <- s[[v]]
  out
}

#' Dirichlet-multinomial shrunk propensities.
.frame_propensity <- function(out, s, req, pri) {
  want <- intersect(c("pi_minus", "pi_plus"), req)
  if (!length(want)) return(out)
  p <- shrink_pi(cbind(s$k_low, s$k_mid, s$k_high), pri$alpha)
  for (v in want) out[[v]] <- p[, v]
  out
}

#' Conditional-prior shrunk magnitudes (R/04b_conditional.R).
#'
#' Driven off the FORMULA's variable names, like the intervention block: a
#' `{var}_delta` in `req` is what triggers the derivation, so the config flag is
#' read in exactly one place (the formula builder) and the frame cannot come to
#' disagree with it.
.frame_magnitude <- function(out, s, req, pri) {
  want <- req[grepl("_delta$", req)]
  if (!length(want)) return(out)
  for (v in want) {
    base <- delta_base_of(v)
    par  <- pri$magnitude[[base]]
    if (is.null(par)) {
      stop(sprintf("signal_frame: the formula asks for `%s` but no magnitude ",
                   "parameters were fitted for `%s`. priors_for() and the formula ",
                   "builder disagree about whether the construct is active.", v, base),
           call. = FALSE)
    }
    out[[v]] <- delta_value(s[[base]], s[[par$count_var]], s$n_obs, par)
  }
  out
}

#' Intervention blocks, one per paired intervention.
#'
#' Still name-driven: which columns are wanted comes from matching the
#' `{intervention}__` prefix against the formula's own variables, so a change to
#' intervention_terms() propagates here with no edit. The only addition is that
#' a wanted name may be DERIVED rather than carried, which is a one-entry
#' registry rather than a branch per variable.
DERIVED_INTERVENTION_VARS <- c("lambda")

.frame_intervention <- function(out, s, req, signal, tabs, cfg, pri) {
  for (iv in interventions_of(signal, cfg)) {
    pre  <- paste0(iv, "__")
    want <- req[startsWith(req, pre)]
    if (!length(want)) next
    z <- tabs$intervention_features[tabs$intervention_features$intervention == iv, , drop = FALSE]
    m <- match(s$stay_id, z$stay_id)
    if (anyNA(m)) {
      stop(sprintf("signal_frame: %d stay(s) absent from intervention_features for '%s'",
                   sum(is.na(m)), iv), call. = FALSE)
    }
    for (v in want) {
      base <- sub(pre, "", v, fixed = TRUE)
      if (base %in% DERIVED_INTERVENTION_VARS) {
        out[[v]] <- .derive_intervention(base, iv, z, m, pri)
      } else {
        out[[v]] <- z[[base]][m]
      }
    }
  }
  out
}

.derive_intervention <- function(base, iv, z, m, pri) {
  if (!identical(base, "lambda")) {
    stop("signal_frame: no derivation registered for `", iv, "__", base, "`", call. = FALSE)
  }
  par <- pri$lambda[[iv]]
  if (is.null(par)) {
    stop(sprintf("signal_frame: the formula asks for `%s__lambda` but no intensity ",
                 "parameters were fitted for `%s`. priors_for() and the formula ",
                 "builder disagree about whether the construct is active.", iv, iv),
         call. = FALSE)
  }
  D <- z[[par$exposure_var]][m] * par$exposure_scale
  A <- z[[par$accumulation_var]][m]
  if (anyNA(D) || anyNA(A)) {
    stop(sprintf("signal_frame [%s]: NA in %s/%s, which intervention_shape says ",
                 "are applicable.", iv, par$exposure_var, par$accumulation_var),
         call. = FALSE)
  }
  lambda_value(A, D, par)
}

#' Assert a frame satisfies the fit contract, before bam() sees it.
#'
#' Three things, all of which mgcv would otherwise report as something else:
#'   - exactly the required columns, no more and no fewer;
#'   - no NA anywhere. `na.action = na.fail` is in force (hard rule 3), so an NA
#'     is a hard stop inside bam() with a message that does not name the column.
#'     The NAs that exist in intervention_features are all in shape-inapplicable
#'     columns — `exposure_frac` for event shapes, `n_hours`/`total_amount` for
#'     state shapes, `n_agents` outside vasopressor and inotrope — and the
#'     formula never asks for those. If one appears here, the shape map and the
#'     data have gone out of sync, which is exactly what this catches;
#'   - no constant column feeding a smooth. mgcv errors with "A term has fewer
#'     unique covariate combinations than specified maximum degrees of freedom",
#'     which reads like a k problem and is not one.
#'   - no empty level in ANY factor. Under a smooth's `by=` mgcv cannot fit over
#'     zero rows and errors out mid-run; as a plain main effect the level's
#'     design column is all zeros and the coefficient is unidentifiable. Both
#'     matter, and the second is quieter, so both are checked here.
#'     `droplevels()` is NOT the fix: the level set has to be identical across
#'     folds and across sites or the bundle cannot predict a patient who falls
#'     in the dropped level. An empty cell is a design decision to be declared,
#'     not something to absorb silently.
check_model_frame <- function(x, f, signal = NULL, model = NULL,
                              stage = c("fit", "predict")) {
  stage <- match.arg(stage)
  tag <- sprintf("model frame [%s %s, %s]", signal %||% "?", model %||% "?", stage)
  req <- required_columns(f)
  got <- setdiff(names(x), "stay_id")
  if (length(setdiff(req, got))) abort_values(paste(tag, "is missing column(s)"), setdiff(req, got))
  if (length(setdiff(got, req))) abort_values(paste(tag, "has surplus column(s)"), setdiff(got, req))

  na_cols <- names(x)[vapply(x, anyNA, logical(1))]
  if (length(na_cols)) {
    abort_values(paste(tag, "has NA in column(s) — na.fail would stop inside bam()"), na_cols)
  }

  # --- fitting-time only ----------------------------------------------------
  # Everything above is required of ANY frame: the right columns, no NAs (which
  # na.fail would reject inside bam), and factor levels the model knows.
  #
  # The two checks below are about CONSTRUCTING a basis, which only happens when
  # fitting. A prediction frame legitimately has fewer distinct values — it is a
  # held-out fold, roughly a fifth of the rows — and requiring it to support the
  # basis its model was fitted with is simply wrong. It also fails in practice:
  # `transfusion_platelet__n_hours` has 9 distinct values across four training
  # folds and 8 within any single held-out one.
  if (identical(stage, "predict")) return(invisible(TRUE))

  smooth_vars <- .term_vars(grep("^s\\(", attr(stats::terms(f), "term.labels"), value = TRUE))
  const <- smooth_vars[vapply(smooth_vars, function(v)
    v %in% names(x) && length(unique(x[[v]])) < 2L, logical(1))]
  if (length(const)) {
    abort_values(paste(tag, "has a constant column under a smooth"), const)
  }

  # A thin-plate basis needs k <= the number of distinct covariate values.
  # mgcv's own message for this ("fewer unique covariate combinations than
  # specified maximum degrees of freedom") names neither the term nor the
  # model, which with 258 fits is close to useless. Checked here, per FITTING
  # SUBSET rather than once on full train, because an out-of-fold subset holds
  # 4/5 of the rows and can thin a covariate the full set supports.
  sp <- smooth_specs(f)
  if (nrow(sp)) {
    sp$n_unique <- vapply(sp$variable, function(v)
      if (v %in% names(x)) length(unique(x[[v]])) else NA_integer_, integer(1))
    bad <- sp[!is.na(sp$n_unique) & !is.na(sp$k) & sp$k >= sp$n_unique, , drop = FALSE]
    if (nrow(bad)) {
      stop(sprintf("%s: k is at or above the distinct-value count for %s. %s. Declare a lower k in config/smooth_k for this signal -- do NOT derive it from the data here (hard rule 8).",
                   tag, if (nrow(bad) == 1L) "a smooth" else "smooths",
                   paste(sprintf("%s: k=%d, %d distinct", bad$variable, bad$k, bad$n_unique),
                         collapse = "; ")),
           call. = FALSE)
    }
  }

  by <- .by_factors(f)
  for (v in names(x)[vapply(x, is.factor, logical(1))]) {
    cells <- table(x[[v]])
    if (any(cells == 0L)) {
      stop(sprintf("%s: `%s` has %d empty level(s) of %d, %s. Cells: %s",
                   tag, v, sum(cells == 0L), length(cells),
                   if (v %in% by) "under a by= smooth — mgcv cannot fit these"
                   else "as a main effect — the coefficient is unidentifiable",
                   paste(sprintf("%s=%d", names(cells), as.integer(cells)), collapse = ", ")),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Factors appearing as the `by=` of a smooth in a formula.
.by_factors <- function(f) {
  lab <- attr(stats::terms(f), "term.labels")
  m <- regmatches(lab, regexpr("by\\s*=\\s*[A-Za-z._][A-Za-z0-9._]*", lab))
  unique(sub("by\\s*=\\s*", "", m))
}

#' One row per (signal, model): shape of the frame the fit will see.
#'
#' Aggregates only — counts and rates, never a row (hard rule 1). This is the
#' cheap thing to look at before committing to the full run: it shows which
#' frames are small, which are near-degenerate in the response, and how the term
#' count scales with the intervention block.
frame_report <- function(tabs, cfg, priors, stay_ids = NULL,
                         role = c("final", "oof"), fold = NA_integer_) {
  role <- match.arg(role)
  rows <- list()
  for (sg in cfg$signals) {
    pri <- priors_for(priors, sg, role, fold)
    for (md in models_of(sg, cfg)) {
      # A frame that fails check_model_frame is the interesting row, so the
      # sweep records it and carries on rather than aborting at signal 5 of 19.
      d <- try(signal_frame(sg, md, tabs, cfg, pri, stay_ids = stay_ids), silent = TRUE)
      if (inherits(d, "try-error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          signal = sg, model = md, n_rows = NA_integer_, n_cols = NA_integer_,
          n_terms = length(attr(stats::terms(build_formula(sg, md, cfg)), "term.labels")),
          n_events = NA_integer_, event_rate = NA_real_,
          status = trimws(sub("^Error[^:]*:", "", conditionMessage(attr(d, "condition")))),
          stringsAsFactors = FALSE)
        next
      }
      y <- d[[all.vars(attr(d, "formula"))[1]]]
      rows[[length(rows) + 1L]] <- data.frame(
        signal    = sg,
        model     = md,
        n_rows    = nrow(d),
        n_cols    = ncol(d) - 1L,
        n_terms   = length(attr(stats::terms(attr(d, "formula")), "term.labels")),
        n_events  = sum(y),
        event_rate = round(mean(y), 4),
        status    = "ok",
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
