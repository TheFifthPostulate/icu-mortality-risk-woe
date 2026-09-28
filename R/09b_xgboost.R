# R/09b_xgboost.R ------------------------------------------------------------
# The comparison arm. A gradient-boosted baseline scored through EXACTLY the
# code path the proposed method is scored through (R/09's score_report), so the
# cells of the 2x2 in docs/v2_analysis_tiering.md are commensurable.
#
# TWO DESIGNS, ONE LEARNER. Which input a cell uses is the whole point:
#
#   raw   the measurement and intervention columns, one wide matrix. Answers
#         "does the proposed method hold up against a strong learner given the
#         same information?"
#   L     the out-of-fold L matrix. Same learner, same attribution method, only
#         the input changes -- so a difference isolates the REPRESENTATION.
#         Against the naive sum on the same matrix it isolates the AGGREGATION.
#
# THE TRAP, and it is easy to fall into: the L matrix fed to a tree must be
# OUT-OF-FOLD, and the tree's own predictions must be out-of-fold too. Scoring an
# in-sample tree against an out-of-fold LLR sum is not a comparison. xgb_oof()
# cross-fits on the SAME folds layer 1 used, which is what makes the two scores
# commensurable rather than merely both present.
#
# MODULAR BY CONSTRUCTION, for eICU (hard rule 8):
#
#   xgb_oof()        cross-fitted, INTERNAL only. Fits.
#   xgb_fit_full()   one model on full train. Fits. Goes in the bundle.
#   xgb_apply()      applies a frozen model. Fits nothing. This is the eICU path.
#
# The design-matrix builders take `feature_names` so an external site is scored
# on the training site's columns, in the training site's order, with anything
# absent left as NA rather than silently dropped or reordered. A tree that saw
# 130 columns at MIMIC and 128 at eICU would report a transportability result
# that is a plumbing artifact.
#
# HYPERPARAMETERS ARE DECLARED IN CONFIG, NEVER TUNED HERE. Same discipline as
# the diagnostic thresholds: tuning a baseline after seeing its AUROC is how a
# comparison arm becomes a foregone conclusion. `nrounds` is the one quantity
# chosen from data, by early stopping on an INNER split of the fitting folds --
# never on the held-out fold, which would leak.
#
# HARD RULE 2, DELIBERATELY AND LOCALLY BROKEN. `xgb_design_raw()` builds a wide
# frame, which layer 1 never does. That is unavoidable for a tree baseline and
# is confined to this file: nothing in R/04-R/07 calls it, and the wide matrix
# is a model INPUT here in the same way the L matrix is a pivot of model
# OUTPUTS. Stated rather than quietly done.
#
# No paths, no clock (hard rule 9).
# ----------------------------------------------------------------------------

#' Booster parameters, read from config. Design information, not a fitted thing.
#'
#' `seed` IS REQUIRED and has no default. It used to fall back to
#' `cfg$xgboost$seed` and then to `1L`, which named a config key that does not
#' exist and never has (audit finding F9's sibling in audit E). Both call sites
#' -- `.xgb_fit1()` from `xgb_oof()`, `xgb_oof_perfold()` and
#' `xgb_fit_full()` -- thread `cfg$seed` explicitly, so the fallback was
#' unreachable as well as undeclared. An unreachable fallback to a
#' reproducibility-critical value is worth deleting rather than documenting.
xgb_params <- function(cfg, seed) {
  x <- cfg$xgboost %||% list()
  list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    max_depth        = x$max_depth %||% 4L,
    eta              = x$eta %||% 0.05,
    subsample        = x$subsample %||% 0.8,
    colsample_bytree = x$colsample_bytree %||% 0.8,
    min_child_weight = x$min_child_weight %||% 10,
    lambda           = x$lambda %||% 1,
    # REQUIRED, not defaulted. `bam.nthreads` is declared as 4 and this read
    # fell back to 1 (audit finding F9). It is not a speed knob: MEASURED
    # 2026-09-03, fitting at 4 threads against 1 moves a GAM's L values by up
    # to 8.8e-13 because parallel accumulation reorders floating-point sums,
    # and the same argument applies to a booster. `x$nthread` keeps its
    # optionality because it genuinely is an override that need not exist.
    nthread          = as.integer(x$nthread %||% cfg_req(cfg, "bam", "nthreads")),
    # WITHOUT THIS THE ARM IS NOT REPRODUCIBLE. subsample and colsample_bytree
    # are both 0.8, so xgboost draws rows and columns from its OWN generator,
    # which R's set.seed does not reach. MEASURED 2026-08-28: two runs of
    # tests/metrics_xgb.R at identical config gave best_iter 465/355/277/272/389
    # and 342/398/337/240/433, and AUROCs differing by up to 0.0005 -- small
    # against the ladder's 0.006-0.020 rungs, but a number in the paper must not
    # move when the script is re-run. `seed` is threaded from the call site so
    # each fold's booster is seeded like each fold's inner split already was.
    seed             = as.integer(seed)
  )
}

xgb_control <- function(cfg) {
  x <- cfg$xgboost %||% list()
  list(nrounds_max = x$nrounds_max %||% 2000L,
       early_stopping_rounds = x$early_stopping_rounds %||% 50L,
       valid_frac = x$valid_frac %||% 0.2)
}

# --- design matrices --------------------------------------------------------

#' The L matrix as an XGBoost design.
#'
#' A thin adapter over R/07's l_matrix() so the tree and the naive sum are fed
#' the identical object -- not two constructions of "the L matrix" that could
#' drift. `fill = "zero"` matches the primary reading: an unmeasured stay
#' contributes no evidence.
xgb_design_L <- function(l_long, tabs, cfg, stay_ids, model = "full",
                         fill = c("zero", "na")) {
  M <- l_matrix(l_long, model, tabs, cfg, stay_ids, fill = match.arg(fill))
  storage.mode(M) <- "double"
  M
}

#' The raw feature matrix: what layer 1 sees, before any of it is modelled.
#'
#' One column per (signal, measurement variable) and per (intervention,
#' intervention variable), named `{signal}__{var}` and `{intervention}__{var}`.
#' Built from the SAME whitelists the formula builder uses, so the baseline is
#' given the same information rather than a curated subset.
#'
#' TWO ASYMMETRIES THAT FAVOUR THE BASELINE, both deliberate and both to be
#' reported rather than removed:
#'
#'  1. UNMEASURED STAYS ARRIVE AS NA, and XGBoost learns a default split
#'     direction for missingness. So the tree may use "nobody ordered this test"
#'     as evidence -- exactly the mechanism the LLR design refuses on
#'     interpretability grounds (spec SS5.5, and `n_obs` never entering a
#'     formula). If the tree wins by a margin, part of that margin is a
#'     mechanism the proposed method excluded on principle. Set
#'     `missing_as_evidence = FALSE` to impute instead and measure the size of
#'     that margin directly.
#'  2. NO FORBIDDEN-VARIABLE GUARD BEYOND `dx_`/`qc_`. The baseline gets
#'     `n_obs`, which layer 1 excludes on transportability grounds. Removing it
#'     would be handicapping the baseline; keeping it means the comparison is
#'     conservative for us. `dx_` and `qc_` stay out at both sites because
#'     `qc_discharge_hospice` is unavailable at eICU and would make any
#'     cross-site result a plumbing artifact.
#'
#' @param feature_names if given, the exact columns and order to emit. Missing
#'   ones are created as NA. This is the eICU path: score an external site on
#'   the training site's design, never on its own.
#' @return numeric matrix, rownames = stay_ids
xgb_design_raw <- function(tabs, cfg, stay_ids, feature_names = NULL,
                           missing_as_evidence = TRUE) {
  ids <- as.character(stay_ids)
  sf  <- tabs$signal_features
  ivf <- tabs$intervention_features

  meas_vars <- intersect(
    c("n_obs", "k_low", "k_mid", "k_high", "q05", "q95",
      "value_min", "value_max", "value_median", "trend"),
    names(sf))
  iv_vars <- intersect(
    c("exposure_frac", "n_hours", "total_amount", "n_agents",
      "present_at_admission"),
    names(ivf))

  cols <- list()
  for (sg in cfg$signals) {
    z <- sf[sf$signal == sg, , drop = FALSE]
    m <- match(ids, as.character(z$stay_id))
    unmeasured <- is.na(m) | (!is.na(m) & z$n_obs[m] == 0)
    for (v in meas_vars) {
      x <- z[[v]][m]
      # An unmeasured stay has no value for a measurement variable. Which is
      # true, and the point of the first asymmetry above.
      x[unmeasured] <- NA_real_
      cols[[paste0(sg, "__", v)]] <- as.numeric(x)
    }
  }
  for (iv in as.character(unlist(cfg$interventions_modelled))) {
    z <- ivf[ivf$intervention == iv, , drop = FALSE]
    if (!nrow(z)) next
    m <- match(ids, as.character(z$stay_id))
    for (v in iv_vars) {
      x <- as.numeric(z[[v]][m])
      # Shape-inapplicable columns are NA by contract (spec SS5.3). Emitting a
      # column that is NA for every row would be a constant, so it is dropped.
      if (all(is.na(x))) next
      cols[[paste0(iv, "__", v)]] <- x
    }
  }

  X <- do.call(cbind, cols)
  rownames(X) <- ids

  if (!isTRUE(missing_as_evidence)) {
    # Median imputation, computed on THIS matrix. Only ever used for the
    # ablation that measures how much the missingness channel is worth; the
    # primary run leaves NA in place.
    #
    # THE FILL FOR A COLUMN WITH NO OBSERVED VALUE IS 0, STATED (plumbing review
    # F8, 2026-09-08). This read `median(...) %||% 0` until then, and `%||%`
    # tests for NULL while `median(numeric(0))` returns NA -- so the fallback
    # never ran, the column stayed all-NA, and the "no missing values" ablation
    # still handed XGBoost missing values. Reproduced on a synthetic design with
    # one wholly unmeasured signal. The class: a null-coalescing operator used
    # on a function that signals "nothing" with NA rather than NULL.
    for (j in seq_len(ncol(X))) {
      na <- is.na(X[, j])
      if (!any(na)) next
      m <- stats::median(X[!na, j])
      X[na, j] <- if (is.finite(m)) m else 0
    }
  }

  align_design(X, feature_names)
}

#' Force a design onto a fixed column set and order.
#'
#' The external-site guard. A column the training site had and this one does not
#' becomes NA (which XGBoost handles); a column this site has and the model
#' never saw is dropped, because the booster has no split for it.
align_design <- function(X, feature_names = NULL) {
  if (is.null(feature_names)) return(X)
  out <- matrix(NA_real_, nrow = nrow(X), ncol = length(feature_names),
                dimnames = list(rownames(X), feature_names))
  shared <- intersect(colnames(X), feature_names)
  if (length(shared)) out[, shared] <- X[, shared, drop = FALSE]
  out
}

# --- fitting ----------------------------------------------------------------

#' The booster fitting protocol, named so a replicate store can be held to it.
#'
#' Part of `attr_design_fingerprint()`: a change to HOW a booster is fitted,
#' with the `xgboost:` config block unchanged, must still refuse a cached SHAP
#' replicate fitted the old way. Bump it when `.xgb_fit1()` changes behaviour.
#'
#'   v1  (to 2026-09-09) early-stopped model returned as fitted on the inner
#'       training rows only; inner split by row, unstratified.
#'   v2  round count chosen on a patient-grouped, mortality-stratified inner
#'       split, then a fresh booster fitted on EVERY offered row for that count.
XGB_FIT_PROTOCOL <- "v2_grouped_refit"

#' The inner early-stopping split: grouped by patient, stratified on outcome.
#'
#' STATISTICAL REVIEW S3 (2026-09-09). The outer split and folds group patients
#' and stratify on mortality (R/03), and until now the inner split did neither:
#' it sampled row indices, so repeat stays of one patient could sit on both
#' sides of the stopping-point decision, and a small or rare-event subset could
#' hand the stopping rule a single-class validation set. This reuses the fold
#' machinery -- `.group_table()` and `.allocate()` -- so the inner split is the
#' same construction as the outer one, at the declared `valid_frac`, under the
#' caller's seed. Both sides are asserted to carry both classes; the stopping
#' rule reads logloss and a single-class side is not a policy anyone declared.
#'
#' @return list(train, valid) of row indices into `y`
.xgb_inner_split <- function(y, group, frac, seed) {
  if (length(group) != length(y)) {
    stop(".xgb_inner_split: `group` is not aligned to the rows", call. = FALSE)
  }
  if (anyNA(group)) stop(".xgb_inner_split: `group` has NA", call. = FALSE)
  g   <- as.character(group)
  pat <- .group_table(g, as.integer(y))
  lab <- with_seed(seed, .allocate(pat, k = NULL, frac = frac))
  valid <- which(g %in% pat$grp[lab == "test"])
  train <- setdiff(seq_along(y), valid)
  if (length(unique(y[valid])) < 2L || length(unique(y[train])) < 2L) {
    stop(sprintf(paste0(".xgb_inner_split: a side of the early-stopping split ",
                        "carries one class (train %d rows, valid %d rows, %d ",
                        "events offered). The stopping rule needs both."),
                 length(train), length(valid), sum(y)), call. = FALSE)
  }
  list(train = train, valid = valid)
}

#' Fit one booster: choose `nrounds` by early stopping on an inner split, then
#' refit on every offered row for that round count.
#'
#' The inner validation rows come out of the FITTING rows only. Using the
#' held-out fold would choose the stopping point on the data the score is then
#' evaluated on, which is the same leak as tuning on test.
#'
#' THE REFIT IS THE POINT (statistical review S2, 2026-09-09). Until now the
#' early-stopped model was returned as fitted, on the inner TRAINING rows only:
#' with `valid_frac: 0.2` the "full-training" booster the bundle carries had
#' seen 80% of the training set through its gradients, the cross-fitted
#' boosters about 64% of it, and `xgb_fit_full()` recorded `n_train = nrow(X)`
#' regardless. The primary GAMs are fitted on the whole training set, so the
#' comparator was being handed less data than the method under a name that said
#' otherwise. Now the early-stopped fit decides ONLY the round count; a fresh
#' booster is then fitted on all offered rows for exactly that many rounds, and
#' that is what is returned and applied. The round count is the one quantity
#' chosen from data, as before; what changed is that every offered outcome now
#' contributes a gradient. `best_iteration` is written onto the refit as an
#' attribute so `.xgb_predict()` pins the same range it always did.
#'
#' @param group patient id per row of `X`, REQUIRED (review S3)
#' @return list(booster, best_iter, nrounds, n_fit, n_select_train,
#'   n_select_valid). Sizes are the actual row counts each stage saw.
.xgb_fit1 <- function(X, y, cfg, seed, group) {
  if (missing(group)) {
    stop(".xgb_fit1: `group` (patient id per row) is required; the inner ",
         "early-stopping split is grouped by patient (review S3).", call. = FALSE)
  }
  p  <- xgb_params(cfg, seed = seed); ct <- xgb_control(cfg)
  sp <- .xgb_inner_split(y, group, ct$valid_frac, seed)

  dtr <- xgboost::xgb.DMatrix(X[sp$train, , drop = FALSE], label = y[sp$train], missing = NA)
  dva <- xgboost::xgb.DMatrix(X[sp$valid, , drop = FALSE], label = y[sp$valid], missing = NA)
  b0 <- xgboost::xgb.train(params = p, data = dtr, nrounds = ct$nrounds_max,
                           evals = list(valid = dva),
                           early_stopping_rounds = ct$early_stopping_rounds,
                           verbose = 0)
  bi <- .xgb_best_iter(b0)
  if (is.na(bi)) {
    stop(".xgb_fit1: early stopping recorded no best iteration; the round ",
         "count cannot be chosen and nothing is refitted.", call. = FALSE)
  }
  nr <- bi + 1L                       # best_iteration is 0-based; rounds are a count
  dall <- xgboost::xgb.DMatrix(X, label = y, missing = NA)
  b <- xgboost::xgb.train(params = p, data = dall, nrounds = nr, verbose = 0)
  xgboost::xgb.attr(b, "best_iteration") <- bi
  list(booster = b, best_iter = bi, nrounds = nr, n_fit = nrow(X),
       n_select_train = length(sp$train), n_select_valid = length(sp$valid))
}

#' The early-stopped round count, and prediction pinned to it.
#'
#' xgboost 3.x keeps `best_iteration` as a BOOSTER ATTRIBUTE (a character), not
#' as a list element -- `b$best_iteration` is NULL and would silently record NA.
#' VERIFIED on 3.1.2.1: predict()'s default already stops at this round, but the
#' range is passed explicitly so the score does not depend on a default that
#' could change under us.
.xgb_best_iter <- function(b) {
  v <- suppressWarnings(as.integer(xgboost::xgb.attr(b, "best_iteration")))
  if (length(v) != 1L || is.na(v)) NA_integer_ else v
}

.xgb_predict <- function(b, X) {
  d  <- xgboost::xgb.DMatrix(X, missing = NA)
  bi <- .xgb_best_iter(b)
  if (is.na(bi)) return(stats::predict(b, d))
  stats::predict(b, d, iterationrange = c(1L, bi + 1L))
}

#' Cross-fitted out-of-fold predictions on the SAME folds layer 1 used.
#'
#' This is what makes the tree score and the summed LLR commensurable. Both are
#' out-of-fold on the identical partition, so neither has seen the row it is
#' scoring.
#'
#' @param X        design matrix, rownames = stay_id
#' @param y        0/1 outcome aligned to rows of X
#' @param fold_of  fold index per row, from R/03
#' @param group    patient id per row, for the inner early-stopping split
#' @return list(score, p_hat, best_iter, n_fit). `score` is on the LOG-ODDS
#'   scale and centred on the training prior, so it is directly comparable to a
#'   summed LLR and can go straight into score_report(). `n_fit` is one row per
#'   fold with the sizes each stage of the fit actually saw.
#'
#' NOT FOR THE STACKED `L` DESIGN. A matrix of out-of-fold L values built from
#' the primary folds is itself a function of every fold's outcomes, so
#' cross-validating a booster over it with the same folds trains on features
#' that were fitted with the held-out fold's outcomes (statistical review S1).
#' The `xgb_l` cell goes through `xgb_oof_perfold()` with
#' `xgb_design_L_nested()`; this function serves designs that carry no fitted
#' quantity (`raw`).
xgb_oof <- function(X, y, fold_of, cfg, seed = 1L, group) {
  stopifnot(nrow(X) == length(y), length(y) == length(fold_of))
  if (missing(group)) stop("xgb_oof: `group` is required (review S3)", call. = FALSE)
  stopifnot(length(group) == length(y))
  p_hat <- rep(NA_real_, length(y))
  best  <- integer(0); nf <- list()
  for (f in sort(unique(fold_of))) {
    ho <- fold_of == f
    r <- .xgb_fit1(X[!ho, , drop = FALSE], y[!ho], cfg, seed = seed + f,
                   group = group[!ho])
    p_hat[ho] <- .xgb_predict(r$booster, X[ho, , drop = FALSE])
    best <- c(best, r$best_iter)
    nf[[length(nf) + 1L]] <- .xgb_fit_sizes(r, fold = f)
  }
  if (anyNA(p_hat)) stop("xgb_oof: a row was never scored; check fold_of", call. = FALSE)
  eps <- 1e-6
  p_bar <- mean(y)
  list(score = logit(pmin(pmax(p_hat, eps), 1 - eps)) - logit(p_bar),
       p_hat = p_hat, best_iter = best, n_fit = do.call(rbind, nf))
}

#' The sizes a fit saw, as one row. Counts only (hard rule 1).
.xgb_fit_sizes <- function(r, fold = NA_integer_) {
  data.frame(fold = as.integer(fold), n_fit = r$n_fit,
             n_select_train = r$n_select_train, n_select_valid = r$n_select_valid,
             best_iter = r$best_iter, nrounds = r$nrounds,
             protocol = XGB_FIT_PROTOCOL, stringsAsFactors = FALSE)
}

#' One booster on the whole training set. This is what the bundle carries and
#' what eICU is scored with (hard rule 8).
#'
#' `n_train` is now what it says (review S2): the refit inside `.xgb_fit1()`
#' sees every row of `X`. The rows the round count was chosen on are recorded
#' beside it rather than folded into a number that did not distinguish them.
#'
#' @param group patient id per row of `X`, for the inner early-stopping split
xgb_fit_full <- function(X, y, cfg, seed = 1L, group) {
  if (missing(group)) stop("xgb_fit_full: `group` is required (review S3)", call. = FALSE)
  r <- .xgb_fit1(X, y, cfg, seed = seed, group = group)
  list(booster = r$booster, feature_names = colnames(X), p_bar = mean(y),
       best_iter = r$best_iter, nrounds = r$nrounds,
       params = xgb_params(cfg, seed = seed), n_train = r$n_fit,
       n_select_train = r$n_select_train, n_select_valid = r$n_select_valid,
       protocol = XGB_FIT_PROTOCOL)
}

#' Apply a frozen booster. FITS NOTHING. The eICU and MIMIC-test path.
#'
#' The design is realigned onto the training feature set first, so a site with a
#' different column set is scored on the model's columns rather than its own.
#' `p_bar` comes from the fitted object, never from the site being scored --
#' re-deriving it here would silently convert the transportability result into a
#' plumbing artifact, exactly as re-deriving alpha would.
xgb_apply <- function(model, X) {
  Xa <- align_design(X, model$feature_names)
  p  <- .xgb_predict(model$booster, Xa)
  eps <- 1e-6
  list(score = logit(pmin(pmax(p, eps), 1 - eps)) - logit(model$p_bar),
       p_hat = p)
}

#' Gain-based importance. A cheap orientation check, NOT an attribution result.
#'
#' Gain is not SHAP and does not decompose a prediction; it is here so a run can
#' be sanity-checked ("is one column doing everything?") without pulling in the
#' attribution arm. The 2x2 attribution comparison needs TreeSHAP and is a
#' separate piece of work.
xgb_importance <- function(model, top_n = 20L) {
  imp <- xgboost::xgb.importance(model = model$booster)
  utils::head(as.data.frame(imp), top_n)
}

# --- the constructed-covariate design ---------------------------------------

#' The layer-1 covariates as an XGBoost design. The middle rung of the ladder.
#'
#' `raw` and `L` are the two ends of a chain that conflates two separate
#' changes: between them the COVARIATES change (counts -> shrunk `pi_hat`,
#' magnitude -> `delta`, intensity -> `lambda`, `n_obs` and `value_median`
#' gone) AND the REPRESENTATION collapses (a per-signal additive GAM to one
#' scalar). `xgb_l - xgb_raw` charges both to one number. This design splits it:
#'
#'   raw  -> feat   what the covariate construction costs, learner held fixed
#'   feat -> L      what collapsing 19 additive models to 19 scalars costs
#'
#' Exactly the columns `signal_frame()` hands `bam()`, for every signal, with no
#' model fitted. Nothing is re-derived here: the measurement block IS a
#' `signal_frame()` call and the intervention block IS `.frame_intervention()`,
#' so a change to the formula builder propagates with no edit, the same
#' guarantee the frame builders already give layer 1.
#'
#' FOLD-AWARENESS IS NOT OPTIONAL. `delta` and `lambda` are FITTED, so a design
#' built once from full-train priors and then cross-fitted would standardise
#' every fold's held-out rows against curves fitted on those same rows. The gap
#' would come back flattering and it would be a plumbing artifact. `role`/`fold`
#' are therefore required arguments threaded to `priors_for()`, exactly as
#' `fit_one()` threads them, and `xgb_oof_perfold()` rebuilds the design per
#' fold rather than accepting one matrix.
#'
#' LAYOUT, and why it is not one block per signal. Interventions repeat across
#' signals (the GCS triple shares nine sedation terms), so a naive union of the
#' 19 `full` frames would emit the same column up to three times. Each
#' intervention is taken ONCE instead, which is sound because `lambda` is keyed
#' on intervention alone -- asserted, not assumed, by `check_lambda_invariance()`.
#' That yields 19 `{signal}__*` groups + 13 `{intervention}__*` groups = the
#' same 32 groups `xgb_design_raw()` aggregates to (docs/v2_state §4.2), so the
#' attribution arm can compare the two cells term block for term block.
#'
#' UNMEASURED STAYS ARRIVE AS NA, and this is a real asymmetry, not an oversight.
#' The L path assigns `L = 0` -- neutral because L is a log-odds CONTRIBUTION.
#' There is no neutral value in covariate space: `pi_minus = 0` asserts "measured
#' and never low" and `q05_delta = 0` asserts "exactly at conditional
#' expectation", both confident claims about a patient nobody tested. Zero-fill
#' would be a fabricated observation. So the missingness channel survives into
#' this cell as it does into `raw`, which means `raw -> feat` isolates the
#' covariate change and NOT the missingness change. Size it the same way: the
#' `--no-missing` ablation put that channel at ~0.001 AUROC.
#'
#' The measurement frame is built with `measured_only = TRUE` and scattered, so
#' `signal_frame()` is never called on an unmeasured row -- the same discipline
#' as layer 1, where `predict()` never sees one.
#'
#' @param role,fold threaded to priors_for(). `oof`/f for the cross-fitted path,
#'   `final`/NA for the full-train booster that goes in the bundle.
#' @return numeric matrix, rownames = stay_ids
xgb_design_feat <- function(tabs, cfg, priors, stay_ids,
                            role = c("oof", "final"), fold = NA_integer_,
                            feature_names = NULL) {
  role <- match.arg(role)
  ids  <- as.character(stay_ids)
  cols <- list()

  # Measurement block. One signal at a time, on that signal's measured subset,
  # then scattered onto the full id vector with NA where unmeasured.
  for (sg in cfg$signals) {
    pri <- priors_for(priors, sg, role, fold)
    d <- signal_frame(sg, "meas", tabs, cfg, pri, stay_ids = stay_ids,
                      measured_only = TRUE, stage = "predict")
    y_name <- all.vars(attr(d, "formula"))[1]
    m <- match(ids, as.character(d$stay_id))
    for (v in setdiff(names(d), c("stay_id", y_name))) {
      cols[[paste0(sg, "__", v)]] <- .feat_numeric(d[[v]])[m]
    }
  }

  # Intervention block. Once per intervention, cohort-wide: `.frame_intervention`
  # reads only `s$stay_id`, and asserts every stay is present in
  # intervention_features, so this cannot silently under-cover.
  s_all <- data.frame(stay_id = stay_ids)
  done  <- character(0)
  for (sg in cfg$signals) {
    want_iv <- setdiff(interventions_of(sg, cfg), done)
    if (!length(want_iv)) next
    req <- required_columns(build_formula(sg, "full", cfg))
    pri <- priors_for(priors, sg, role, fold)
    d <- .frame_intervention(data.frame(stay_id = stay_ids), s_all, req,
                             sg, tabs, cfg, pri)
    for (iv in want_iv) {
      pre <- paste0(iv, "__")
      for (v in names(d)[startsWith(names(d), pre)]) cols[[v]] <- .feat_numeric(d[[v]])
    }
    done <- c(done, want_iv)
  }

  X <- do.call(cbind, cols)
  rownames(X) <- ids
  align_design(X, feature_names)
}

#' Factors and logicals into a tree's only currency, without silent surprises.
.feat_numeric <- function(x) {
  if (is.factor(x))  return(as.numeric(as.integer(x)))
  if (is.logical(x)) return(as.numeric(x))
  as.numeric(x)
}

#' Cross-fitted OOF where the DESIGN, not just the model, is refitted per fold.
#'
#' `xgb_oof()` takes one matrix, which is correct whenever the design carries no
#' fitted quantity (`raw`) or is already out-of-fold (`L`). It is wrong for any
#' design containing `delta` or `lambda`. This variant takes a builder instead
#' and calls it once per held-out fold.
#'
#' The column set is asserted identical across folds: the formula builder is
#' data-free, so a fold-varying design would mean a fold-varying FORMULA, and
#' the five boosters would not be five fits of one model.
#'
#' @param design_fn function(fold) -> design matrix for the model that holds out
#'   `fold`, covering ALL scored rows in the order `y` and `fold_of` use.
#' @param group     patient id per row, for the inner early-stopping split
xgb_oof_perfold <- function(design_fn, y, fold_of, cfg, seed = 1L, group) {
  stopifnot(length(y) == length(fold_of))
  if (missing(group)) stop("xgb_oof_perfold: `group` is required (review S3)", call. = FALSE)
  stopifnot(length(group) == length(y))
  p_hat <- rep(NA_real_, length(y))
  best  <- integer(0); nf <- list()
  ref   <- NULL
  for (f in sort(unique(fold_of))) {
    X <- design_fn(f)
    if (nrow(X) != length(y)) {
      stop("xgb_oof_perfold: design for fold ", f, " has ", nrow(X),
           " rows, expected ", length(y), call. = FALSE)
    }
    if (is.null(ref)) ref <- colnames(X)
    else if (!identical(colnames(X), ref)) {
      stop("xgb_oof_perfold: fold ", f, " produced a different column set. The ",
           "formula builder is data-free, so this means a formula changed ",
           "between folds.", call. = FALSE)
    }
    ho <- fold_of == f
    r <- .xgb_fit1(X[!ho, , drop = FALSE], y[!ho], cfg, seed = seed + f,
                   group = group[!ho])
    p_hat[ho] <- .xgb_predict(r$booster, X[ho, , drop = FALSE])
    best <- c(best, r$best_iter)
    nf[[length(nf) + 1L]] <- .xgb_fit_sizes(r, fold = f)
  }
  if (anyNA(p_hat)) stop("xgb_oof_perfold: a row was never scored", call. = FALSE)
  eps <- 1e-6
  p_bar <- mean(y)
  list(score = logit(pmin(pmax(p_hat, eps), 1 - eps)) - logit(p_bar),
       p_hat = p_hat, best_iter = best, feature_names = ref,
       n_fit = do.call(rbind, nf))
}

#' The stacked `L` design for ONE outer fold, with no outcome of that fold in
#' any feature the booster trains on.
#'
#' STATISTICAL REVIEW S1 (2026-09-09). `design_l` -- the primary out-of-fold L
#' matrix -- is out-of-fold ROW BY ROW: a stay in fold B carries an L from GAMs
#' fitted on the other four folds. But those four folds include A, so when the
#' booster holds out fold A and trains on fold B's rows, its training FEATURES
#' were fitted with fold A's outcomes. A row's own prediction being out-of-fold
#' does not make the stacked procedure out-of-fold. Reproduced synthetically in
#' the review: changing only fold A's outcomes moved the second-stage training
#' features by 0.25 while A's own features did not move.
#'
#' The honest design for outer fold A is therefore assembled from two sources:
#'
#'   rows in A       the ordinary OOF L for fold A (`l_oof`, fitted on the
#'                   other four folds). These are the VALIDATION rows and their
#'                   L never saw A.
#'   rows in B != A  the L from the NESTED pair fit {A, B} (`l_nested`): GAMs
#'                   and priors fitted on the three folds outside both A and B,
#'                   predicting B. These are the TRAINING rows and their L
#'                   never saw A either.
#'
#' `l_matrix()` then asserts that every measured training stay appears exactly
#' once, so a pair that is missing or doubled is an error, not a silent gap.
#' Patient grouping is preserved at both levels because both sources are built
#' on the primary fold partition.
#'
#' STATED, NOT HIDDEN: the two sources are centred on different `p_bar`s. A
#' pair fit's L is `logit(p_hat) - logit(p_bar)` with `p_bar` from three
#' folds; the validation rows' L uses the four-fold `p_bar`. The gap is the
#' sampling difference between two estimates of one measured-subpopulation
#' rate (of order 0.01 on the logit scale at these sizes) and is inherent to
#' nested cross-fitting of a stacked model, exactly as the primary design
#' already trains the final booster on four-fold L and applies it through
#' five-fold L. It is a property of the estimand, not a defect to correct
#' by re-centring, which would re-introduce fold outcomes into the features.
#'
#' @param l_nested  the long table from `layer1_nested_l()` over every pair
#' @param fold      the outer fold the booster will hold out
#' @return numeric matrix over ALL `stay_ids`, in their order
xgb_design_L_nested <- function(l_oof, l_nested, tabs, cfg, stay_ids, fold_of,
                                fold, model = "full") {
  stopifnot(length(stay_ids) == length(fold_of))
  need <- c("signal", "model", "stay_id", "l")
  ho_ids <- as.character(stay_ids[fold_of == fold])
  a <- l_oof[l_oof$role == "oof" & !is.na(l_oof$fold) & l_oof$fold == fold, need,
             drop = FALSE]
  if (!all(as.character(a$stay_id) %in% ho_ids)) {
    stop("xgb_design_L_nested: `l_oof` fold ", fold, " carries a stay outside ",
         "that fold. Ids not printed (hard rule 1).", call. = FALSE)
  }
  b <- l_nested[(l_nested$fold_a == fold | l_nested$fold_b == fold) &
                  !(as.character(l_nested$stay_id) %in% ho_ids), need, drop = FALSE]
  if (!nrow(b)) {
    stop("xgb_design_L_nested: no nested pair rows touch fold ", fold, call. = FALSE)
  }
  xgb_design_L(rbind(a, b), tabs, cfg, stay_ids, model = model, fill = "zero")
}

#' Every column a booster was GIVEN, with what it did with each.
#'
#' `xgb.importance()` reports only features that appear in at least one split,
#' so reading it as the design is wrong in the one direction that matters: a
#' column the booster never used is silently absent rather than reported as
#' zero. This starts from `model$feature_names` -- the design as offered -- and
#' joins importance onto it, so unused columns appear with gain 0 and
#' `used = FALSE`. "Which of the 94 covariates does the tree ignore?" is a real
#' diagnostic and it is not answerable from the importance table alone.
#'
#' `group` is the `{signal}` / `{intervention}` prefix, i.e. the term block the
#' attribution arm aggregates to (docs/v2_state SS4.2). A column carrying no `__`
#' IS its own group -- that is the `L` design, whose 19 columns are bare signal
#' names and each of which is already one term block. Without that case every
#' `xgb_l` group is NA and the rollup has nothing to aggregate.
#' attribution arm aggregates to (docs/v2_state §4.2). Aggregates only; no row
#' data ever enters this table.
#'
#' @return data frame, one row per design column, ordered by gain descending
xgb_feature_table <- function(model, design = NA_character_) {
  fn  <- model$feature_names
  imp <- as.data.frame(xgboost::xgb.importance(model = model$booster))
  m   <- match(fn, imp$Feature)
  z <- data.frame(
    design    = design,
    feature   = fn,
    group     = ifelse(grepl("__", fn, fixed = TRUE), sub("__.*$", "", fn), fn),
    used      = !is.na(m),
    gain      = ifelse(is.na(m), 0, imp$Gain[m]),
    cover     = ifelse(is.na(m), 0, imp$Cover[m]),
    frequency = ifelse(is.na(m), 0, imp$Frequency[m]),
    stringsAsFactors = FALSE)
  z[order(-z$gain, z$feature), , drop = FALSE]
}

#' Gain rolled up to the term block. The precursor to the SHAP aggregation.
#'
#' Gain is NOT an attribution result (see xgb_importance) and this is not the
#' 2x2 deliverable -- it is the cheap orientation check at the group level, and
#' the object that makes two designs with different column counts comparable at
#' all: `raw` splits a signal over ~10 columns and `feat` over ~3-5, so
#' per-column gain is not comparable between them and per-group gain is.
xgb_group_gain <- function(ft) {
  g <- ft[!is.na(ft$group), , drop = FALSE]
  z <- stats::aggregate(cbind(gain, frequency) ~ group, data = g, FUN = sum)
  n <- stats::aggregate(cbind(n_cols = used) ~ group, data = g, FUN = length)
  u <- stats::aggregate(cbind(n_used = used) ~ group, data = g, FUN = sum)
  z <- merge(merge(z, n, by = "group"), u, by = "group")
  z$design <- ft$design[1]
  z[order(-z$gain), c("design", "group", "gain", "frequency", "n_cols", "n_used")]
}
