# R/13b_attribution_floor.R --------------------------------------------------
# The noise floors: the denominator without which an agreement score is a
# number with no referent.
#
# WHAT A FLOOR IS, OPERATIONALLY. Fix an evaluation set and never move it.
# Perturb the TRAINING data. Refit everything that is fitted. Re-attribute on
# those same fixed rows. Compute the agreement metric. The distribution of that
# statistic over replicates is the floor, and a cross-site number is read
# against it in one sentence: "the MIMIC-versus-eICU sign agreement is x, the
# within-site floor is y with interquartile range a to b, so the site effect is
# or is not larger than what one site's own sampling produces."
#
# THE EVALUATION SET NEVER MOVES, AND THAT IS NOT A DETAIL. Resampling it too
# would confound sampling of the evaluation cohort with sampling of the
# estimator, and the estimator is the only thing under study.
#
# EVERYTHING THAT IS FITTED IS REFITTED. alpha, the delta conditional means, the
# lambda conditional means, the 64 final GAMs and the three boosters. Freezing
# the priors and refitting only the GAMs would understate the floor, because the
# priors are fitted parameters and their sampling variability is part of the
# estimator's. This is why a replicate costs minutes rather than seconds.
#
# TWO FLOORS, AND THEY ANSWER DIFFERENT QUESTIONS.
#
#   disjoint_half  Two halves of train sharing no patient, fitted
#                  independently. THE CORRECT NULL FOR A TRANSPORT CLAIM,
#                  because MIMIC and eICU share no patients either. If
#                  cross-site agreement falls inside this floor, the honest
#                  conclusion is that there is no detectable site effect at
#                  all -- a cleaner and stronger result than anything the
#                  bootstrap can support.
#
#   bootstrap      Resample with replacement, refit. Estimates the SAMPLING
#                  VARIABILITY of the agreement statistic. But two bootstrap
#                  resamples share about 63% of their distinct patients, so
#                  they agree MORE than two independent samples would: the
#                  bootstrap floor sits too high, and reading a cross-site
#                  number against it biases toward declaring a site effect
#                  that is really sample-to-sample variation.
#
# They are complementary rather than alternatives. The bootstrap gives the
# variance of the statistic; the disjoint half gives the location of its null.
# Both are TRAINING-SITE ONLY, so neither waits on eICU.
#
# PAIRWISE, NOT AGAINST THE FULL-SAMPLE FIT. Each comparison is between two
# replicates, because that is what a cross-site comparison is -- neither site is
# "the truth" -- and comparing every replicate against the full-sample fit would
# produce a floor that is systematically too tight.
#
# FAILURES ARE COUNTED, NEVER DROPPED. A resample changes the size of every
# signal's measured subset, and a thin signal can fall below its declared basis
# dimension in `config/smooth_k`, which `check_model_frame()` rejects outright.
# A floor computed only over the replicates that happened to succeed is a floor
# conditioned on the estimator having worked, which is exactly the conditioning
# that makes it too narrow.
#
# AGGREGATES ONLY (hard rule 1). The attribution matrices stay in memory and go
# nowhere; what leaves is correlations, quantiles and counts.
# ----------------------------------------------------------------------------

#' Refit the whole fitted stack on a subset of training stays.
#'
#' Returns a bundle-shaped object so `attribution_set()` runs against it through
#' the IDENTICAL code path it uses for the real bundle. Building a lookalike by
#' hand instead would let the replicate and the real thing drift apart, which is
#' the failure `verify_bundle()` exists to prevent one level up.
#'
#' `layer2`, `sigma`, `cutpoints` and `train_ref` are deliberately empty: a
#' replicate is never scored and never reported, it only supplies the second
#' argument to an agreement metric.
#' @param l_source WHERE THE REPLICATE'S OWN L MATRIX COMES FROM, and it is a
#'   methodological choice rather than a plumbing one.
#'
#'   The real `xgb_l` booster is trained on `l_oof` -- OUT-OF-FOLD L over train
#'   (`_targets.R`: `design_l <- xgb_design_L(l_oof, ...)`). A replicate has to
#'   supply the same thing, and until 2026-09-04 this function did not: it ran
#'   `run_layer1()` restricted to `role == "final"` and passed `r$l` to
#'   `xgb_design_L()`. But `job_ids("final", ...)` returns `predict_ids =
#'   NULL` BY DESIGN -- a final fit produces no L values, because predicting on
#'   the rows it was fitted on is what out-of-fold exists to avoid. So `r$l`
#'   was NULL on every replicate and `l_matrix()` failed on the first signal.
#'   Every one of the first five replicates failed for this reason and nothing
#'   in the loop could have worked (finding F19).
#'
#'   "out_of_fold" fits the OOF jobs as well, so the replicate's booster is
#'   trained on the same kind of object the real one is. It is FAITHFUL and it
#'   costs about five times the GAM budget per replicate.
#'
#'   "in_sample" fits only the final jobs and derives L by applying them to
#'   their own fitting rows. It is CHEAP and it changes the estimator: in-sample
#'   L is less noisy, so two replicates agree more than two genuine refits of
#'   the reported estimator would, and the floor sits TOO HIGH -- which biases
#'   toward declaring a site effect that is really sampling variation. Same
#'   direction of error as the bootstrap's 63% patient overlap. Usable for a
#'   fast first look; declare it if any number from it is reported.
#'
#'   Only `xgb_l` is affected either way. `llr_sum` and `llr_meas` come from the
#'   final GAMs directly, and `xgb_feat` and `xgb_raw` are built from covariates
#'   rather than from L.
refit_replicate <- function(tabs, cfg, folds, fit_ids, seed = 1L,
                            l_source = c("out_of_fold", "in_sample"),
                            verbose = FALSE) {
  l_source <- match.arg(l_source)
  f2 <- folds[folds$stay_id %in% fit_ids, , drop = FALSE]
  # Every fitted quantity below is scoped by `split == "train"`, so a replicate
  # is expressed by narrowing that column rather than by threading a subset
  # through six functions that would each need a new argument.
  #
  # The `fold` column is left ALONE, which is what makes "out_of_fold" work:
  # the stays in this half keep their original fold labels, so cross-fitting
  # inside the replicate partitions the half rather than re-drawing folds that
  # would differ between the two halves of a pair.
  f2$split <- "train"
  pri <- layer1_priors(tabs, f2, cfg, verbose = FALSE)
  jb  <- layer1_jobs(cfg)
  jb  <- if (identical(l_source, "out_of_fold")) jb[jb$fit, , drop = FALSE]
         else jb[jb$fit & jb$role == "final", , drop = FALSE]
  r   <- run_layer1(tabs, cfg, f2, pri, jobs = jb, keep_final = TRUE, verbose = FALSE)

  l_long <- if (identical(l_source, "out_of_fold")) r$l else
    apply_layer1(r$models, tabs, cfg, pri, fit_ids, verbose = FALSE)$l
  if (is.null(l_long) || !nrow(l_long)) {
    stop("refit_replicate: no L rows for the xgb_l design under l_source = '",
         l_source, "'. A final-only fit produces none by design.", call. = FALSE)
  }

  y <- as.integer(tabs$cohort$mortality[match(as.character(fit_ids),
                                              as.character(tabs$cohort$stay_id))])
  # The inner early-stopping split groups by patient (review S3).
  grp <- patient_group_of(tabs$cohort, cfg, fit_ids)
  dl  <- xgb_design_L(l_long, tabs, cfg, fit_ids, model = "full", fill = "zero")
  dr  <- xgb_design_raw(tabs, cfg, fit_ids)
  xgb <- list(
    xgb_l    = xgb_fit_full(dl, y, cfg, seed = seed, group = grp),
    xgb_raw  = xgb_fit_full(dr, y, cfg, seed = seed, group = grp),
    xgb_feat = xgb_fit_full(xgb_design_feat(tabs, cfg, pri, fit_ids, role = "final"),
                            y, cfg, seed = seed, group = grp))
  build_bundle(cfg = cfg, priors = pri, models = lapply(r$models, strip_gam),
               sigma = list(), eigen = list(), cutpoints = list(), xgb = xgb,
               layer2 = NULL, train_ref = list(), severity = NULL,
               domains = NULL, site = "mimic")
}

#' The stay sets each perturbation draws.
#'
#' GROUPED BY PATIENT, matching `assign_folds()`. At MIMIC the cohort is one
#' stay per patient so grouping is a no-op, but expressing it here rather than
#' assuming it means a cohort that later admits repeat stays cannot leak a
#' patient across a disjoint-half boundary without anyone noticing.
#'
#' @param group_by a cohort column, ALREADY RESOLVED by `resample_cols()`.
#'   There is no default and no fallback on purpose: the fallback this
#'   function used to carry silently downgraded patient grouping to stay
#'   grouping whenever the declared name needed resolving, which is the failure
#'   the grouping exists to prevent, arriving quietly (audit finding F13).
#' @return list of length-2 lists, each a pair of stay-id vectors to compare
.floor_pairs <- function(cohort, train_ids, perturbation, B, seed, group_by) {
  if (!group_by %in% names(cohort)) {
    abort_values(".floor_pairs: grouping column not in the cohort", group_by)
  }
  grp <- as.character(cohort[[group_by]][match(as.character(train_ids),
                                        as.character(cohort$stay_id))])
  if (anyNA(grp)) stop(".floor_pairs: a training stay has no cohort row", call. = FALSE)
  ug <- unique(grp)
  set.seed(seed)
  if (identical(perturbation, "disjoint_half")) {
    lapply(seq_len(B), function(b) {
      p <- sample(ug)
      h <- split(p, seq_along(p) %% 2L)
      list(train_ids[grp %in% h[[1]]], train_ids[grp %in% h[[2]]])
    })
  } else {
    lapply(seq_len(B), function(b) {
      a1 <- sample(ug, length(ug), replace = TRUE)
      a2 <- sample(ug, length(ug), replace = TRUE)
      # `%in%` de-duplicates, so a bootstrap here resamples PATIENTS and keeps
      # each drawn patient once rather than replicating rows. That is the
      # cluster bootstrap the grouping implies; replicating a stay would also
      # replicate it inside every prior's denominator.
      list(train_ids[grp %in% a1], train_ids[grp %in% a2])
    })
  }
}

#' Build a within-site noise floor.
#'
#' @param eval_ids the FIXED evaluation set. Never resampled.
#' @param perturbation "disjoint_half" or "bootstrap"
#' @param B replicate pairs. `v2_analysis_tiering.md` item 9 fixes the bootstrap
#'   at 50 rather than 200: a gap only visible at 200 is not a gap worth
#'   claiming. The disjoint half is cheap enough to run at 5 and is the one to
#'   run first, because it exercises the whole loop for a few minutes of compute.
#' @return list(overall, per_signal, ranking, failures, timing) -- all
#'   aggregates. `timing` is the compute cost per replicate, elapsed and CPU;
#'   it belongs in a run directory and never in a target (hard rule 7).
#' @param intervention_handling,w,min_signals THE CONVENTIONS, THREADED RATHER
#'   THAN DEFAULTED. A floor exists to be read against a finding, and that
#'   comparison is only meaningful if both were computed under the SAME
#'   conventions: the same rollup of intervention SHAP, the same layer-2
#'   weights, the same minimum signal count for a rankable stay. Until
#'   2026-09-04 this function accepted none of the three and let
#'   `attribution_set()` and `compare_attributions()` supply their own
#'   defaults, so a finding computed under `intervention_handling = "split"`
#'   or under non-unit `w` would have been read against a floor computed under
#'   neither, and nothing would have said so. Same defect class as audit
#'   findings F9 to F11, one level up (finding F17).
attribution_floor <- function(bundle, tabs, cfg, eval_ids, folds,
                              perturbation = c("disjoint_half", "bootstrap"),
                              B = 5L, methods = ATTRIBUTION_METHODS,
                              mag_floor = 0.05, seed = 1L,
                              intervention_handling = "drop", w = NULL,
                              min_signals = 5L, l_source = "out_of_fold",
                              verbose = TRUE) {
  perturbation <- match.arg(perturbation)
  train_ids <- folds$stay_id[folds$split == "train"]
  if (any(eval_ids %in% train_ids) && verbose) {
    message("attribution_floor: the evaluation set overlaps train. That is ",
            "legitimate for a floor -- the estimator is what is being ",
            "resampled -- but it is not legitimate for a performance number.")
  }
  # RESOLVED, not read raw. `cfg$split$group_by` is the declaration
  # "patient_id"; the MIMIC cohort calls that column `subject_id`, and
  # `resample_cols()` is the one place in the project that knows so. Reading
  # the key directly here meant `.floor_pairs()` looked for a column the cohort
  # does not have and fell through to its own `stay_id` fallback, so the floor
  # was grouped by STAY while claiming in its own comment to be grouped by
  # patient. It changes nothing at MIMIC, where the cohort is one stay per
  # patient, and it would leak a patient across a disjoint-half boundary in any
  # cohort that is not (audit finding F13).
  #
  # The FOLD-level grouping, because a floor partitions train, which is what
  # `assign_folds()` does at level 2. It defaults to the split-level key, so
  # the two are the same declaration until someone deliberately separates them.
  pairs <- .floor_pairs(tabs$cohort, train_ids, perturbation, B, seed,
                        group_by = resample_cols(tabs$cohort, cfg)$fold_group)

  # TIMED WITH `proc.time()`, NOT `Sys.time()`. A floor is the one arm in this
  # project whose cost is itself reportable -- 50 replicates each refitting 43
  # GAMs and three boosters is the sizing decision behind
  # `v2_analysis_tiering.md` item 9 -- so the run has to be timed. It is timed
  # off the process counters rather than the calendar, which measures the same
  # duration while keeping hard rule 9 intact; see `start_timer()` in
  # `R/00_utils.R` for why those are different things.
  timer_all <- start_timer()
  ov <- list(); ps <- list(); rk <- list(); fails <- list(); tmg <- list()
  for (b in seq_along(pairs)) {
    lab <- sprintf("%s_%02d", perturbation, b)
    timer <- start_timer()
    res <- tryCatch({
      # `cfg` is reused unchanged, and the design stamp therefore still matches:
      # a replicate is built by `build_bundle()` from the SAME cfg the real
      # bundle froze, so `.hash(replicate$cfg)` equals `cfg$.bundle_design` and
      # `attribution_set()`'s hard-rule-8 check passes for the right reason
      # rather than by being bypassed.
      r1 <- refit_replicate(tabs, cfg, folds, pairs[[b]][[1]], seed = seed + b,
                            l_source = l_source)
      r2 <- refit_replicate(tabs, cfg, folds, pairs[[b]][[2]], seed = seed + b + 1000L,
                            l_source = l_source)
      s1 <- attribution_set(r1, tabs, cfg, eval_ids, methods = methods, w = w,
                            intervention_handling = intervention_handling,
                            verbose = FALSE)
      s2 <- attribution_set(r2, tabs, cfg, eval_ids, methods = methods, w = w,
                            intervention_handling = intervention_handling,
                            verbose = FALSE)
      list(s1 = s1, s2 = s2)
    }, error = function(e) {
      fails[[length(fails) + 1L]] <<- data.frame(
        replicate = lab, message = substr(conditionMessage(e), 1, 300),
        stringsAsFactors = FALSE)
      NULL
    })
    # A FAILED REPLICATE STILL COSTS COMPUTE, and its cost is recorded for the
    # same reason its failure is: a cost summed only over the replicates that
    # succeeded understates what the arm takes to run.
    el <- timer()
    tmg[[length(tmg) + 1L]] <- data.frame(
      replicate = lab, status = if (is.null(res)) "failed" else "ok",
      elapsed_min = el$elapsed_sec / 60, cpu_min = el$cpu_sec / 60,
      stringsAsFactors = FALSE)
    if (is.null(res)) {
      if (verbose) message(sprintf("  %s: FAILED after %.1f min, recorded",
                                   lab, el$elapsed_sec / 60))
      next
    }
    for (nm in methods) {
      cp <- compare_attributions(res$s1, res$s2, nm, mag_floor = mag_floor,
                                 min_signals = min_signals, label = lab)
      ov[[length(ov) + 1L]] <- cp$overall
      ps[[length(ps) + 1L]] <- cp$per_signal
      rk[[length(rk) + 1L]] <- cp$ranking
    }
    if (verbose) {
      message(sprintf("  %s: %.1f min elapsed, %.1f min CPU", lab,
                      el$elapsed_sec / 60, el$cpu_sec / 60))
    }
  }
  all_el <- timer_all()
  if (verbose) {
    message(sprintf("%s: %d replicate(s), %.1f min elapsed, %.1f min CPU",
                    perturbation, length(pairs), all_el$elapsed_sec / 60,
                    all_el$cpu_sec / 60))
  }
  list(overall = do.call(rbind, ov), per_signal = do.call(rbind, ps),
       ranking = do.call(rbind, rk),
       # The conventions this floor was computed under, carried WITH it so a
       # finding read against it can be checked for having used the same ones.
       conventions = list(intervention_handling = intervention_handling,
                          mag_floor = mag_floor, min_signals = min_signals,
                          weighted = !is.null(w), l_source = l_source),
       failures = if (length(fails)) do.call(rbind, fails) else
         data.frame(replicate = character(0), message = character(0)),
       n_attempted = length(pairs), n_failed = length(fails),
       # WRITE IT INTO A RUN DIRECTORY, NEVER INTO A TARGET (hard rule 7). A
       # duration differs between two runs of identical code, so a target
       # carrying `timing` would rebuild everything downstream of it on every
       # run. Nothing in `_targets.R` reaches this function today and
       # `tests/audit_c_targets.R` reports it as CLOCK if anything ever does.
       timing = if (length(tmg)) do.call(rbind, tmg) else
         data.frame(replicate = character(0), status = character(0),
                    elapsed_min = numeric(0), cpu_min = numeric(0)),
       elapsed_min = all_el$elapsed_sec / 60, cpu_min = all_el$cpu_sec / 60,
       perturbation = perturbation)
}

#' Collapse a floor into the two numbers a finding is read against.
#'
#' Median and interquartile range per method, never a single number: a floor
#' reported as a point estimate invites exactly the comparison it exists to
#' prevent.
floor_summary <- function(fl) {
  # AN ALL-FAILED FLOOR IS A RESULT AND MUST NOT BE A CRASH. `fl$overall` is
  # NULL when no replicate produced a comparison, and `split(NULL, NULL)` dies
  # with "first argument must be a vector" -- which is what happened on
  # 2026-09-04 after all five replicates failed, destroying the failure table
  # that was the only useful output of a 23-minute run (finding F20). The
  # failures are the finding in that case, so the summary comes back empty and
  # correctly shaped and the caller writes everything it has.
  if (is.null(fl$overall) || !nrow(fl$overall)) {
    return(data.frame(perturbation = fl$perturbation, method = character(0),
                      n_replicates = integer(0),
                      spearman_median = numeric(0), spearman_q1 = numeric(0),
                      spearman_q3 = numeric(0), sign_median = numeric(0),
                      sign_q1 = numeric(0), sign_q3 = numeric(0),
                      stringsAsFactors = FALSE)[0, , drop = FALSE])
  }
  f <- function(d, col) {
    z <- d[[col]]; z <- z[is.finite(z)]
    if (!length(z)) return(c(NA, NA, NA, NA))
    c(stats::median(z), stats::quantile(z, 0.25), stats::quantile(z, 0.75), length(z))
  }
  do.call(rbind, lapply(split(fl$overall, fl$overall$method), function(d) {
    sp <- f(d, "spearman"); sg <- f(d, "sign_agreement")
    data.frame(perturbation = fl$perturbation, method = d$method[1],
               n_replicates = nrow(d),
               spearman_median = round(sp[1], 4),
               spearman_q1 = round(sp[2], 4), spearman_q3 = round(sp[3], 4),
               sign_median = round(sg[1], 4),
               sign_q1 = round(sg[2], 4), sign_q3 = round(sg[3], 4),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}
