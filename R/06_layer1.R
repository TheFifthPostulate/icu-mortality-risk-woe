# R/06_layer1.R --------------------------------------------------------------
# Layer 1: one GAM per (signal, model), and the L value it produces.
#
#   L = logit(p_hat) - logit(p_bar_train)
#
# out-of-fold for the L's layer 2 is estimated on, and from the full-train fits
# for test and eICU.
#
# Three models per signal (R/05_formula.R): `meas`, `intv`, `full`. Their term
# sets partition exactly, so
#
#   L_cond = L_full - L_intv   is the measurement-deviation LLR conditional on
#                              class and intervention context
#
# and p_bar cancels in that subtraction, because both models are fitted on the
# same rows with the same prior. That is why fit_one() takes `p_bar` as an
# argument instead of deriving it: two models that derived it separately could
# derive it differently, and the subtraction would silently stop being the
# conditional term.
#
# HARD RULE 6. A fold fit returns predictions plus a one-row diagnostics frame
# and NOTHING ELSE. The gam object is dropped before this function returns, so
# every diagnostic has to be extracted while it is still alive — R/08 does that,
# and it is called here rather than downstream for exactly that reason. Only
# full-train fits keep the object, and R/10 strips `$model` before saving.
#
# HARD RULE 3. na.action = na.fail in every bam() call, asserted rather than
# assumed. mgcv's default silently drops NA rows.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). Nothing this file returns into the
# graph contains a wall-clock number; benchmark_fit() measures at the call site
# and keeps the result out of every value.
# ----------------------------------------------------------------------------

# --- the job table ----------------------------------------------------------

#' Every distinct layer-1 fit, enumerated. Data-free, so it can be inspected
#' before a single row is read.
#'
#' Two kinds of row that are NOT fitted, both marked `fit = FALSE`:
#'
#'   alias   For the 7 unpaired signals `meas` and `full` are the same formula
#'           by construction, so `full` aliases `meas` rather than fitting the
#'           identical model a second time. The two INTERACTION models alias
#'           there too, and for a second reason on top: an unpaired signal has
#'           no intervention to cross with. `full_ti_trend` additionally aliases
#'           `full` on the three paired signals whose class carries no `trend`
#'           covariate, where the cross-term set is empty and the two formulas
#'           are the same object.
#'   assign  `intv` does not exist for an unpaired signal: its formula would be
#'           `mortality ~ 1`, giving p_hat = p_bar and L = 0 exactly. That is
#'           the correct value, so it is ASSIGNED, never fitted — the same
#'           treatment unmeasured stays get (spec §5.5). Sending 35
#'           intercept-only models to bam() would produce the same numbers at
#'           real cost and make it look as though something was estimated.
#'
#' NEITHER RULE IS DECIDED HERE ANY MORE. Both come from `spec_source()` in
#' R/05_formula.R, which `l_matrix()` also reads, so the enumeration and the
#' pivot cannot disagree about which model a cell borrows from. They were two
#' independent copies of the same two `if` statements until 2026-09-07, which
#' survived three models and would not have survived five.
#'
#' @return data frame: signal, model, role, fold, fit, source
#'         `source` is the model whose predictions this row takes ("meas" or
#'         "full" for an alias), "zero" for an assigned row, or NA when fitted.
layer1_jobs <- function(cfg) {
  n_folds <- cfg$n_folds %||% 5L
  roles <- rbind(
    data.frame(role = "oof",   fold = seq_len(n_folds), stringsAsFactors = FALSE),
    data.frame(role = "final", fold = NA_integer_,      stringsAsFactors = FALSE)
  )

  rows <- list()
  for (sg in cfg$signals) {
    for (md in LAYER1_MODELS) {
      src <- spec_source(sg, md, cfg)
      fit <- is.na(src)
      for (i in seq_len(nrow(roles))) {
        rows[[length(rows) + 1L]] <- data.frame(
          signal = sg, model = md, role = roles$role[i], fold = roles$fold[i],
          fit = fit, source = src, stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' One-line summary of the fit budget. Cheap, and it is what to check before
#' committing to a run.
layer1_budget <- function(cfg) {
  j <- layer1_jobs(cfg)
  data.frame(
    specs_distinct = length(unique(paste(j$signal, j$model)[j$fit])),
    fits_oof       = sum(j$fit & j$role == "oof"),
    fits_final     = sum(j$fit & j$role == "final"),
    fits_total     = sum(j$fit),
    # Counted by WHAT THEY ALIAS rather than as one lump, because the two mean
    # different things: `alias_meas` is a signal with no intervention at all,
    # `alias_full` is a paired signal whose interaction set happens to be empty.
    alias_meas     = sum(!j$fit & j$source == "meas"),
    alias_full     = sum(!j$fit & j$source == "full"),
    assigned_zero  = sum(!j$fit & j$source == "zero"),
    stringsAsFactors = FALSE)
}

# --- row scoping ------------------------------------------------------------

#' Which stays a job fits on, and which it predicts on.
#'
#' The two are disjoint by construction for an out-of-fold job — that is the
#' whole point — and a final job predicts on nothing here: test is touched once,
#' by the apply path, not by the fitting path.
#'
#' @return list(fit_ids, predict_ids)
job_ids <- function(role, fold, folds) {
  tr <- folds$split == "train"
  if (identical(role, "final")) {
    return(list(fit_ids = folds$stay_id[tr], predict_ids = NULL))
  }
  if (is.na(fold)) stop("job_ids: an oof job needs a fold", call. = FALSE)
  list(fit_ids     = folds$stay_id[tr & folds$fold != fold],
       predict_ids = folds$stay_id[tr & folds$fold == fold])
}

# --- the fit ----------------------------------------------------------------

#' Fit one layer-1 GAM and return its L values plus its diagnostics row.
#'
#' @param signal,model  what to fit; `model` is one of LAYER1_MODELS
#' @param tabs   from load_tables()
#' @param cfg    from load_config()
#' @param pri    the frozen quantities for this (signal, role, fold), from
#'               priors_for(): alpha, p_bar, the magnitude parameters and the
#'               intensity parameters, in ONE object.
#'
#'               Passed in, never derived here, so a held-out fold cannot be
#'               standardised using its own data. It is one object rather than
#'               four arguments because L_full - L_intv is only the conditional
#'               term when both models centre on the identical p_bar and are
#'               built from the identical covariate construction -- and four
#'               arguments can be got individually wrong at a call site while
#'               one object cannot.
#' @param fit_ids,predict_ids  from job_ids()
#' @param keep_model  TRUE only for full-train fits (hard rule 6)
#' @param role,fold   carried through onto the diagnostics row
#' @return list(l, diagnostics, model). `l` is a data frame (stay_id, l) over
#'         the MEASURED subset of predict_ids, or NULL when nothing was
#'         predicted. Unmeasured stays are zero-filled by R/07, not here — this
#'         function never predicts on a row the model was not fitted to cover.
fit_one <- function(signal, model, tabs, cfg, pri,
                    fit_ids, predict_ids = NULL, keep_model = FALSE,
                    role = NA_character_, fold = NA_integer_) {
  f <- build_formula(signal, model, cfg)
  p_bar <- pri$p_bar

  d <- signal_frame(signal, model, tabs, cfg, pri, stay_ids = fit_ids)
  b <- .bam_fit(f, d, cfg)

  # Extracted here, while the object is alive. After this point it may be gone.
  diag <- model_diagnostics(
    b, signal = signal, model = model, role = role, fold = fold, cfg = cfg,
    extras = list(p_bar = round(p_bar, 6), n_fit_ids = length(fit_ids)))

  # ACCEPTED BEFORE IT PREDICTS. A fit that failed is an error here, never a
  # flagged row that scores anyway (review S10). Policy in `.accept_outcome_fit()`.
  .accept_outcome_fit(b, diag, cfg)

  l <- NULL
  if (length(predict_ids)) {
    # stage = "predict": the held-out fold need not support the basis the model
    # was fitted with, only supply its columns without NAs. See
    # check_model_frame() for why that distinction is load-bearing.
    nd <- signal_frame(signal, model, tabs, cfg, pri,
                       stay_ids = predict_ids, stage = "predict")
    # type = "link" IS logit(p_hat) for a binomial fit, so L is a subtraction
    # rather than a round trip through inv_logit and back.
    #
    # discrete = FALSE deliberately: the model was fitted discretely, but
    # predict.bam()'s discrete path re-uses the fitted discretisation and is the
    # fragile one on new data. Prediction is a rounding error next to the fit,
    # so the robust path is the right trade.
    eta <- as.numeric(stats::predict(b, newdata = nd, type = "link", discrete = FALSE))
    if (anyNA(eta)) {
      stop(sprintf("fit_one [%s %s]: %d prediction(s) are NA. na.fail should have ",
                   "caught this in the frame; ids are not printed (hard rule 1).",
                   signal, model, sum(is.na(eta))), call. = FALSE)
    }
    l <- data.frame(stay_id = nd$stay_id, l = eta - logit(p_bar),
                    stringsAsFactors = FALSE)
    diag$n_predict <- nrow(l)
    diag$l_mean    <- round(mean(l$l), 5)
    diag$l_sd      <- round(stats::sd(l$l), 5)
  } else {
    diag$n_predict <- 0L
    diag$l_mean <- NA_real_
    diag$l_sd   <- NA_real_
  }

  # Hard rule 6, enforced rather than trusted.
  list(l = l, diagnostics = diag, model = if (keep_model) b else NULL)
}

#' THE OUTCOME-FIT ACCEPTANCE POLICY (statistical review S10, 2026-09-09).
#'
#' Until now `fit_one()` extracted the diagnostics row and then predicted
#' regardless of what it said: `diagnostics.require_convergence: true` added a
#' triage FLAG, and nothing between a numerically unsuccessful `bam()` and an
#' exported L refused to proceed. The prior fits acquired a strict status gate
#' on 2026-09-08 (`check_prior_fits()`); this is the equivalent for the models
#' that produce the evidence itself, and it is deliberately narrow:
#'
#'   * a NON-FINITE COEFFICIENT is always fatal. No policy makes a prediction
#'     from one meaningful.
#'   * under `require_convergence`, `converged` must be TRUE. An UNKNOWN status
#'     (NA: the object reported convergence in no recognised place) is refused
#'     as well -- "not known to have failed" is not acceptance.
#'   * k-index, edf ratio, concurvity and deviance-explained thresholds stay
#'     TRIAGE. They describe a fit that may need reading; they do not say it
#'     failed, and turning them into failures would be the indiscriminate
#'     policy the review warned against.
#'
#' A quasi-separated fit ("fitted probabilities numerically 0 or 1") that still
#' converges with finite coefficients is accepted and its warning recorded by
#' the graph; the eventual L is a large finite evidence value, which is the
#' honest reading of that fit. Cached diagnostics were audited before this
#' gate was added: all 384 stored fits carry `converged = TRUE`.
.accept_outcome_fit <- function(b, diag, cfg) {
  where <- sprintf("%s/%s [%s%s]", diag$signal, diag$model, diag$role,
                   if (is.na(diag$fold)) "" else paste0(" fold ", diag$fold))
  cf <- stats::coef(b)
  if (!all(is.finite(cf))) {
    stop(sprintf(paste0("fit_one %s: %d of %d coefficient(s) are non-finite. The ",
                        "fit is not accepted and nothing is predicted from it."),
                 where, sum(!is.finite(cf)), length(cf)), call. = FALSE)
  }
  if (cfg_flag(cfg, "diagnostics", "require_convergence",
               what = "declared in config/config.yml under `diagnostics`")) {
    if (is.na(diag$converged)) {
      stop(sprintf(paste0("fit_one %s: convergence status is UNKNOWN (the object ",
                          "reports it in no recognised field). Under ",
                          "diagnostics.require_convergence an unknown status is ",
                          "not accepted."), where), call. = FALSE)
    }
    if (!isTRUE(diag$converged)) {
      stop(sprintf(paste0("fit_one %s: bam() did not converge. Under ",
                          "diagnostics.require_convergence the fit is refused ",
                          "rather than flagged; resolve the fit, do not relax ",
                          "the policy locally."), where), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' The only bam() call in the project.
#'
#' Every setting comes from config, and na.action is ASSERTED to be na.fail
#' rather than merely defaulted to it (hard rule 3). mgcv's default drops NA
#' rows silently, which would delete every unexposed patient from a fit without
#' saying so.
#' The bam settings, as ONE object, with no defaults anywhere.
#'
#' EVERY SETTING IS REQUIRED (`cfg_req`). Until 2026-09-03 each was read as
#' `bs$<key> %||% <literal>`, which wrote every one of them down twice -- in
#' `config/config.yml` and again as the fallback -- with nothing holding the
#' two together. Audit E found `nthreads` declared as 4 and defaulted to 1
#' (finding F9), and the difference is not cosmetic: MEASURED 2026-09-03 on
#' `mbp/meas` over 41,185 training stays, fitting at `nthreads = 4` against
#' `nthreads = 1` moves L by up to 8.8e-13, because parallel accumulation
#' reorders floating-point sums. At a FIXED value the fit is bitwise
#' reproducible (0.000e+00 over a re-run), so `nthreads` is a numerical setting
#' that must be pinned rather than a speed knob that may drift.
#'
#' Returning one object rather than reading eight keys inline also means the
#' settings can be reported: `bam_settings_row()` turns this into a one-line
#' record for a run log, so "which settings produced this fit" is answerable
#' from the run rather than from whatever `config.yml` happens to say today.
#'
#' `na_action` is asserted rather than defaulted (hard rule 3): mgcv's default
#' drops NA rows silently, which would delete every unexposed patient from a
#' `peak_intensity` smooth without saying so.
bam_settings <- function(cfg) {
  what <- "every bam setting is declared in config/config.yml under `bam:`"
  na_action <- cfg_req(cfg, "bam", "na_action", what = what)
  if (!identical(na_action, "na.fail")) {
    stop("bam.na_action must be na.fail (hard rule 3), got: ", na_action, call. = FALSE)
  }
  fam_name <- cfg_req(cfg, "bam", "family", what = what)
  fam <- switch(fam_name,
                binomial = stats::binomial(),
                stop("unsupported bam family: ", fam_name, call. = FALSE))
  s <- list(
    family_name = fam_name,
    family      = fam,
    method      = cfg_req(cfg, "bam", "method",       what = what),
    discrete    = cfg_flag(cfg, "bam", "discrete",    what = what),
    nthreads    = as.integer(cfg_req(cfg, "bam", "nthreads", what = what)),
    select      = cfg_flag(cfg, "bam", "select",      what = what),
    gamma       = as.numeric(cfg_req(cfg, "bam", "gamma", what = what)),
    gc_level    = as.integer(cfg_req(cfg, "bam", "gc_level", what = what)),
    na_action   = na_action)
  # gamma < 1 would select ROUGHER fits than the criterion asks for, which is
  # the opposite of why it is here. Refused rather than clamped.
  if (is.na(s$gamma) || s$gamma < 1) {
    stop("bam.gamma must be >= 1, got: ", cfg_req(cfg, "bam", "gamma"), call. = FALSE)
  }
  if (is.na(s$nthreads) || s$nthreads < 1L) {
    stop("bam.nthreads must be a positive integer, got: ",
         cfg_req(cfg, "bam", "nthreads"), call. = FALSE)
  }
  s
}

#' The effective bam settings as a one-row data frame, for a run log.
#' Aggregates only; contains no data (hard rule 1).
bam_settings_row <- function(cfg) {
  s <- bam_settings(cfg)
  data.frame(family = s$family_name, method = s$method, discrete = s$discrete,
             nthreads = s$nthreads, select = s$select, gamma = s$gamma,
             gc_level = s$gc_level,
             na_action = s$na_action, k_default = cfg_req(cfg, "bam", "k"),
             smooth_basis = cfg_req(cfg, "bam", "smooth_basis"),
             stringsAsFactors = FALSE)
}

.bam_fit <- function(f, d, cfg) {
  s <- bam_settings(cfg)
  mgcv::bam(
    formula   = f,
    data      = d,
    family    = s$family,
    method    = s$method,
    discrete  = s$discrete,
    nthreads  = s$nthreads,
    select    = s$select,
    gamma     = s$gamma,
    gc.level  = s$gc_level,
    na.action = stats::na.fail
  )
}

# --- the sweep --------------------------------------------------------------

#' Run every fitted job in the table.
#'
#' Deliberately a plain loop and not a targets pattern: `_targets.R` will map
#' over layer1_jobs() rows itself, and this exists for the serial path — smoke
#' runs, a single-signal re-fit, and the first end-to-end pass before the graph
#' is wired.
#'
#' @param jobs  subset of layer1_jobs(), or NULL for all fitted rows
#' @return list(l, diagnostics, models)
#'   l           long: signal, model, role, fold, stay_id, l
#'   diagnostics rbind of the per-fit rows
#'   models      named list of the final fits, when keep_final = TRUE
run_layer1 <- function(tabs, cfg, folds, priors, jobs = NULL,
                       keep_final = FALSE, verbose = TRUE) {
  if (is.null(jobs)) jobs <- layer1_jobs(cfg)
  jobs <- jobs[jobs$fit, , drop = FALSE]

  ls <- list(); ds <- list(); ms <- list()

  for (i in seq_len(nrow(jobs))) {
    sg <- jobs$signal[i]; md <- jobs$model[i]
    role <- jobs$role[i]; fd <- jobs$fold[i]

    ids <- job_ids(role, fd, folds)
    pri <- priors_for(priors, sg, role, fd)
    keep <- keep_final && role == "final"

    if (verbose) {
      message(sprintf("  [%3d/%3d] %-18s %-5s %-5s fold=%s",
                      i, nrow(jobs), sg, md, role, if (is.na(fd)) "-" else fd))
    }

    r <- fit_one(sg, md, tabs, cfg, pri,
                 fit_ids = ids$fit_ids, predict_ids = ids$predict_ids,
                 keep_model = keep, role = role, fold = fd)

    if (!is.null(r$l)) {
      ls[[length(ls) + 1L]] <- data.frame(
        signal = sg, model = md, role = role, fold = fd,
        r$l, stringsAsFactors = FALSE)
    }
    ds[[length(ds) + 1L]] <- r$diagnostics
    if (keep) ms[[paste(sg, md, sep = "/")]] <- r$model
  }

  out <- list(l = if (length(ls)) do.call(rbind, ls) else NULL,
              diagnostics = do.call(rbind, ds),
              models = ms)
  if (verbose) diagnostics_summary(out$diagnostics, cfg)
  out
}

# --- the nested cross-fit ---------------------------------------------------

#' Layer-1 L for the NESTED cross-fit the stacked `xgb_l` cell needs.
#'
#' STATISTICAL REVIEW S1 (2026-09-09); the design argument is in
#' `xgb_design_L_nested()` (R/09b). This function produces its ingredient: for
#' one unordered pair of folds {a, b}, every signal's `model` spec is fitted on
#' the training rows OUTSIDE both folds -- priors included, through a fold
#' table whose `split` column excludes them -- and predicted on the rows of
#' both. When the booster later holds out fold a, its training rows in fold b
#' read this L, which was fitted without a's outcomes; when it holds out b, the
#' rows in a read the same L. Ten pair fits per signal therefore serve all
#' twenty (outer, inner) cells.
#'
#' WHAT IT REUSES, DELIBERATELY. `layer1_priors(roles = "final")` on the
#' narrowed fold table is exactly the final-prior estimator on a smaller
#' training set, so the nested priors are the primary estimator and not a
#' second implementation of it; `fit_one()` is the primary fit; the alias walk
#' is `resolve_spec_source()`. Nothing here is a new model. Patient grouping is
#' inherited from the primary fold partition.
#'
#' HARD RULE 6: nothing fitted here is kept. HARD RULE 8: these fits serve the
#' training site's own comparator and never a bundle.
#'
#' @param pair  integer(2): the two folds held out together
#' @param model the design's model, `"full"` for the `xgb_l` cell
#' @return list(l, diagnostics). `l` is long: signal, model, role = "nested",
#'   fold_a, fold_b, stay_id, l, over the measured stays of both folds.
layer1_nested_l <- function(tabs, cfg, folds, pair, model = "full", verbose = FALSE) {
  pair <- sort(as.integer(pair))
  if (length(pair) != 2L || anyNA(pair) || pair[1] == pair[2]) {
    abort_values("layer1_nested_l: `pair` must be two distinct folds", pair)
  }
  tr   <- folds$split == "train"
  held <- tr & !is.na(folds$fold) & folds$fold %in% pair
  if (!any(held)) abort_values("layer1_nested_l: no training rows in folds", pair)

  # The two folds leave `train` in the table the priors are fitted from. They
  # keep their fold labels, which nothing reads under `roles = "final"`.
  f2 <- folds
  f2$split[held] <- "nested_holdout"
  pri <- layer1_priors(tabs, f2, cfg, verbose = FALSE, roles = "final")

  fit_ids     <- folds$stay_id[tr & !held]
  predict_ids <- folds$stay_id[held]

  ls <- list(); ds <- list()
  for (sg in cfg$signals) {
    src <- resolve_spec_source(sg, model, cfg)
    if (identical(src, "zero")) next
    if (verbose) message(sprintf("  nested {%d,%d} %-18s %-5s", pair[1], pair[2], sg, src))
    r <- fit_one(sg, src, tabs, cfg, priors_for(pri, sg, "final", NA_integer_),
                 fit_ids = fit_ids, predict_ids = predict_ids,
                 keep_model = FALSE, role = "nested", fold = NA_integer_)
    if (is.null(r$l) || !nrow(r$l)) {
      stop(sprintf("layer1_nested_l [%s %s]: the pair fit predicted nothing", sg, src),
           call. = FALSE)
    }
    ls[[length(ls) + 1L]] <- data.frame(
      signal = sg, model = src, role = "nested",
      fold_a = pair[1], fold_b = pair[2], r$l, stringsAsFactors = FALSE)
    d <- r$diagnostics
    d$fold_a <- pair[1]; d$fold_b <- pair[2]
    ds[[length(ds) + 1L]] <- d
  }
  list(l = do.call(rbind, ls), diagnostics = do.call(rbind, ds))
}

# --- the apply path ---------------------------------------------------------
# The half of layer 1 that FITS NOTHING. `fit_one()` fits and predicts in one
# call because a fold fit must hand back its diagnostics while the gam object is
# still alive (hard rule 6). MIMIC-test and eICU have no such fit: they receive
# the 64 final GAMs frozen into a bundle and evaluate them.
#
# Everything the evaluation needs beyond the gam object is a FITTED QUANTITY and
# arrives as an argument -- alpha, p_bar, the delta parameters, the lambda
# parameters -- through the same `priors_for()` container `fit_one()` uses. That
# is what makes hard rule 8 enforceable here rather than aspirational: there is
# no code path in this function that could re-derive any of them, because it
# never sees the tables they would be derived from.

#' Evaluate one frozen layer-1 GAM on a new set of stays.
#'
#' @param b        a fitted bam object, `$model` already stripped. Prediction
#'                 with `newdata` does not read `$model`, so a stripped object
#'                 predicts exactly as an unstripped one does.
#' @param pri      from priors_for(priors, signal, "final"). `p_bar` is the
#'                 TRAINING prior and is never re-derived from the site being
#'                 scored: doing so would silently convert the transportability
#'                 result into a plumbing artifact (hard rule 8).
#' @param stay_ids stays to score. Restricted to the measured subset inside
#'                 signal_frame(), exactly as the fitting path is -- predict()
#'                 is never called on an unmeasured row, whose L is assigned 0
#'                 at pivot time by R/07 (spec SS5.5).
#' @return data frame (stay_id, l) over the measured subset of `stay_ids`
apply_one <- function(b, signal, model, tabs, cfg, pri, stay_ids) {
  if (is.null(b)) {
    stop(sprintf("apply_one [%s %s]: no fitted model was supplied. An absent ",
                 "GAM must fail here, never fall through to a default -- a ",
                 "silently unscored signal is a transportability artifact.",
                 signal, model), call. = FALSE)
  }
  nd <- signal_frame(signal, model, tabs, cfg, pri,
                     stay_ids = stay_ids, stage = "predict")

  # discrete = FALSE for the same reason fit_one() uses it: predict.bam()'s
  # discrete path re-uses the fitted discretisation and is the fragile one on
  # new data. At an external site "new data" is the whole point.
  eta <- as.numeric(stats::predict(b, newdata = nd, type = "link", discrete = FALSE))
  if (anyNA(eta)) {
    stop(sprintf("apply_one [%s %s]: %d prediction(s) are NA. check_model_frame ",
                 "should have caught this; ids are not printed (hard rule 1).",
                 signal, model, sum(is.na(eta))), call. = FALSE)
  }
  data.frame(stay_id = nd$stay_id, l = eta - logit(pri$p_bar),
             stringsAsFactors = FALSE)
}

#' Evaluate every frozen layer-1 GAM, producing the same long L table
#' `run_layer1()` produces.
#'
#' The output shape is identical to `run_layer1()$l` on purpose: R/07's
#' `l_matrix()` is then the single pivot for both sites, so the alias rule (an
#' unpaired signal's `full` reads `meas`) and the assignment rule (its `intv` is
#' 0 on measured stays) cannot be implemented twice and drift. Only FITTED jobs
#' are emitted here, exactly as at training time.
#'
#' @param models  named list of final GAMs, keyed `"<signal>/<model>"`, from
#'                run_layer1(keep_final = TRUE) or a bundle.
#' @param priors  an `llr_priors` container carrying the FINAL rows. A container
#'                still carrying fold rows would work, but the bundle must not
#'                ship one (CLAUDE.md: never fold-level alpha, never fold GAMs).
#' @return list(l, coverage). `l` is long: signal, model, role, fold, stay_id, l.
#'         `coverage` is one row per fitted spec -- counts only (hard rule 1).
apply_layer1 <- function(models, tabs, cfg, priors, stay_ids, verbose = TRUE) {
  jobs <- layer1_jobs(cfg)
  jobs <- jobs[jobs$fit & jobs$role == "final", , drop = FALSE]

  want <- paste(jobs$signal, jobs$model, sep = "/")
  miss <- setdiff(want, names(models))
  if (length(miss)) {
    abort_values(paste0("apply_layer1: the bundle is missing final GAM(s). ",
                        "Every fitted spec must be present; there is no default"), miss)
  }

  ls <- list(); cov <- list()
  for (i in seq_len(nrow(jobs))) {
    sg <- jobs$signal[i]; md <- jobs$model[i]
    key <- paste(sg, md, sep = "/")
    pri <- priors_for(priors, sg, "final", NA_integer_)

    if (verbose) message(sprintf("  [%3d/%3d] apply %-18s %-5s", i, nrow(jobs), sg, md))

    z <- apply_one(models[[key]], sg, md, tabs, cfg, pri, stay_ids)
    ls[[i]] <- data.frame(signal = sg, model = md, role = "final",
                          fold = NA_integer_, z, stringsAsFactors = FALSE)
    cov[[i]] <- data.frame(
      signal = sg, model = md, n_scored = nrow(z),
      n_offered = length(stay_ids),
      frac_scored = round(nrow(z) / length(stay_ids), 5),
      p_bar_train = round(pri$p_bar, 6),
      l_mean = round(mean(z$l), 5), l_sd = round(stats::sd(z$l), 5),
      l_min = round(min(z$l), 5), l_max = round(max(z$l), 5),
      stringsAsFactors = FALSE)
  }
  list(l = do.call(rbind, ls), coverage = do.call(rbind, cov))
}

# --- benchmarking -----------------------------------------------------------

#' Time a single fit. Step 1 of docs/next_steps.md.
#'
#' Timing lives HERE, at the call site, and is returned to the caller rather
#' than into any value the graph stores (hard rule 7). Nothing this function
#' returns should ever become a target: `elapsed_sec` changes on every run, so a
#' target carrying it would invalidate everything downstream of it and force all
#' 258 GAMs to recompute in order to redraw a plot.
#'
#' @return list(elapsed_sec, diagnostics) — print it, do not cache it
benchmark_fit <- function(signal, model, tabs, cfg, folds, priors,
                          role = "final", fold = NA_integer_) {
  ids <- job_ids(role, fold, folds)
  pri <- priors_for(priors, signal, role, fold)

  t <- system.time(
    r <- fit_one(signal, model, tabs, cfg, pri,
                 fit_ids = ids$fit_ids, predict_ids = ids$predict_ids,
                 keep_model = FALSE, role = role, fold = fold)
  )
  el <- unname(t[["elapsed"]])
  message(sprintf("benchmark %s/%s [%s]: %.2f s, %d rows, %d terms, %d smooths",
                  signal, model, role, el, r$diagnostics$n_rows,
                  r$diagnostics$n_terms, r$diagnostics$n_smooth))
  list(elapsed_sec = el, diagnostics = r$diagnostics)
}
