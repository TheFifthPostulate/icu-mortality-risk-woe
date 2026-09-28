# _targets.R -----------------------------------------------------------------
# The dependency graph for the INTERNAL (MIMIC-IV, training) pipeline.
#
# WHAT IS IN THIS GRAPH, AND WHAT IS DELIBERATELY NOT.
#
# In:   everything computed on MIMIC TRAIN. Loading, validation, folds, the four
#       fitted parameter sets, the 320 out-of-fold GAM fits, the 64 final GAM
#       fits, the L matrices, Sigma, the eigenspectrum, the three XGBoost
#       comparison cells out-of-fold, the three full-train boosters, the frozen
#       reporting cut points, and the bundle those all go into.
#
#       384 FITS, NOT 258, AS OF 2026-09-07. `full_ti_trend` and `full_ti_all`
#       joined LAYER1_MODELS, adding 21 distinct specs and 126 fits. The budget
#       is enumerated by `layer1_budget(cfg)` and printed by every runner, so
#       the number above is a description rather than a second declaration.
#       PLUS THE NESTED CROSS-FIT AS OF 2026-09-09: ten fold-pair fits of the
#       `full` spec per signal (~190 GAMs, ten final-only prior fits) that
#       serve only the stacked `xgb_l` comparator (statistical review S1).
#       They are outside `layer1_budget()` because they enter no bundle.
#
# Out:  MIMIC TEST, and eICU. Both are APPLY sites and neither belongs in a
#       graph whose job is fitting. Test is touched ONCE (CLAUDE.md), and a
#       target that scores it would be re-run by anyone who typed `tar_make()`
#       to rebuild a plot. `run/test_look.R` and `run/external.R` load the
#       bundle this graph produces and fit nothing.
#
# HARD RULE 7, THE ONE THAT MUST NOT BE BROKEN. No target's VALUE may contain
# `Sys.time()`. The moment one does, its hash changes on every run, everything
# downstream of it rebuilds, and 384 GAMs are recomputed to redraw a figure.
# `tar_make()` computes and caches; `export_run()` -- a plain function in
# R/11_run.R, never a target -- snapshots the cached values into a dated
# directory afterwards. That split is the whole discipline.
#
# HARD RULE 1. `_targets/` holds row-level values (the L matrices, the design
# matrices). It is a LOCAL cache under data governance identical to `data/`:
# nothing in it may be printed, and no target below returns anything a console
# would render as rows. Every `tar_read()` a runner performs goes to a file.
#
# GRANULARITY: ONE BRANCH PER FOLD, NOT ONE PER FIT. 320 branches would each
# carry a dependency on `tabs` (about a million rows of signal features) and the
# per-branch overhead would exceed the fits. Five branches of 64 fits each keeps
# the useful caching -- re-run one fold, not all five -- without paying that.
#
#   Rscript -e 'targets::tar_make()'
#   Rscript -e 'targets::tar_read(eigen_full_zero)'
#   Rscript -e 'targets::tar_visnetwork()'
# ----------------------------------------------------------------------------

library(targets)

tar_option_set(
  packages = c("mgcv", "arrow", "yaml", "xgboost", "qs2", "digest"),
  format   = "rds",
  # persistent: `tabs` is the dependency of nearly every target and re-reading
  # it from the store per target would dominate the runtime.
  memory   = "persistent",
  garbage_collection = TRUE,
  error    = "stop"
)

tar_source("R")

list(

  # --- inputs, tracked as files ---------------------------------------------
  # format = "file" hashes CONTENTS, so editing a whitelist invalidates exactly
  # what depends on it. The parquets are tracked too: a re-extraction must
  # rebuild the fits, and discovering that by hand is how a stale result gets
  # published.
  tar_target(config_file,  "config/config.yml",   format = "file"),
  tar_target(pairing_file, "config/pairing.csv",  format = "file"),
  tar_target(domains_file, "config/domains.csv",  format = "file"),

  # THE CODE ITSELF, tracked as an input. Added 2026-09-03 (audit finding F3).
  # `build_bundle()` used to call `.source_hashes()`, which listed and hashed
  # `R/` from inside the bundle target -- so the contents of `R/` decided part
  # of that target's value while `targets` could not see it, and the value was
  # not a function of its declared dependencies. Tracking `R/` as a file target
  # converts that hidden read into a declared edge.
  #
  # WHAT IT COSTS: editing ANY file under `R/`, comment-only edits included,
  # now rebuilds `code_hashes`, `bundle`, `bundle_checks` and
  # `bundle_contents`. Those four are assembly and verification, not fitting --
  # no GAM recomputes -- so the price is seconds and the provenance is honest.
  tar_target(source_files,
             sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
             format = "file"),
  tar_target(code_hashes, source_hashes(source_files)),

  tar_target(cfg, load_config(config_file, pairing_file)),

  tar_target(mimic_files, unlist(cfg$paths$mimiciv, use.names = FALSE),
             format = "file"),

  # --- load and validate ----------------------------------------------------
  # `mimic_files` is referenced rather than used, and it MUST be. Without the
  # reference the file target is tracked and nothing depends on it, so a
  # re-extraction changes the hash, `mimic_files` rebuilds, and `tabs` is
  # skipped -- which is exactly the stale-result failure the comment above
  # claims the file tracking prevents. FOUND 2026-09-03, when the `se_` columns
  # were added: the pipeline validated a cached table that predated them and
  # errored a hundred lines away.
  tar_target(tabs, {
    mimic_files
    load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
  }),

  # The nine checks. A target rather than a side effect, so a failed validation
  # blocks every fit downstream of it instead of printing a warning nobody reads.
  tar_target(validation, validate_tables(tabs, cfg, strict = TRUE)),

  tar_target(domains, load_domains(domains_file)),

  # --- folds ----------------------------------------------------------------
  tar_target(folds, {
    validation
    assign_folds(tabs$cohort, cfg)
  }),
  # `cohort` is passed so the check reads the SAME four resampling
  # declarations `assign_folds()` read, through the same `resample_cols()`.
  # Without it the check can only assert `subject_id` and `mortality`, which
  # would be checking something other than what was done the moment those keys
  # name anything else (audit finding F5).
  tar_target(folds_ok, check_folds(folds, cfg, cohort = tabs$cohort)),
  tar_target(fold_table, fold_summary(folds)),

  tar_target(train_ids, folds$stay_id[folds$split == "train"]),
  tar_target(y_train, {
    v <- tabs$cohort$mortality[match(as.character(train_ids),
                                     as.character(tabs$cohort$stay_id))]
    if (anyNA(v)) stop("y_train: a training stay has no cohort row", call. = FALSE)
    as.integer(v)
  }),
  tar_target(fold_of_train, {
    v <- folds$fold[match(train_ids, folds$stay_id)]
    if (anyNA(v)) stop("fold_of_train: a training stay has no fold", call. = FALSE)
    v
  }),
  tar_target(p_bar_cohort_train, mean(y_train)),

  # THE BOOTSTRAP AND INNER-SPLIT UNIT: the patient of every training stay
  # (statistical review S3/S4, 2026-09-09). Resolved through `resample_cols()`
  # -- the one place `patient_id` becomes `subject_id` -- so the reporting
  # intervals and the boosters' early-stopping split group by exactly the
  # column the folds group by. Character, never printed (hard rule 1).
  tar_target(group_of_train, patient_group_of(tabs$cohort, cfg, train_ids)),

  # --- the four fitted parameter sets ---------------------------------------
  # One object at three grains (alpha/p_bar per signal-fold, delta per
  # signal-variable-fold, lambda per intervention-fold). Built once so the row
  # scoping cannot silently diverge between them.
  tar_target(priors, {
    folds_ok
    layer1_priors(tabs, folds, cfg, verbose = FALSE)
  }),

  # Static design declarations re-checked against the data, in BOTH directions.
  # These are the checks that must also run at eICU, which is why each of them
  # takes `tabs` rather than reading a file.
  tar_target(design_checks, {
    list(
      resampling_cols = check_resampling_cols(tabs, cfg, strict = TRUE),
      excursion_sides = check_excursion_sides(tabs, cfg, strict = TRUE),
      signal_tails    = check_signal_tails(tabs, cfg, strict = TRUE),
      ordinal_scales  = check_ordinal_scales(tabs, cfg, strict = TRUE),
      smooth_k        = check_smooth_k(tabs, cfg, priors, stay_ids = train_ids,
                                       strict = TRUE),
      agent_pool      = check_agent_pool(tabs, cfg, strict = TRUE),
      lambda_invariant = check_lambda_invariance(tabs, cfg, priors,
                                                 stay_ids = train_ids,
                                                 role = "final", strict = TRUE)
    )
  }),

  # Prior-fit diagnostics (R/04c). Free once the priors exist, and the
  # `lambda_fit` table is the cleanest cross-site quantity the project has --
  # differencing it between MIMIC and eICU is the whole of the prior-transport
  # reading.
  tar_target(dm_shrinkage, dm_shrinkage_table(priors$signal,
                                              .measured_train_rows(tabs, folds, cfg),
                                              role = "final")),
  tar_target(delta_fit, delta_fit_diagnostics(priors$magnitude,
                                              .measured_train_rows(tabs, folds, cfg),
                                              cfg, role = "final")),
  tar_target(delta_ordinal_fit, delta_ordinal_diagnostics(priors$magnitude,
                                              .measured_train_rows(tabs, folds, cfg),
                                              cfg, role = "final", seed = cfg$seed)),
  # TRAINING ROWS ONLY (statistical review S8, 2026-09-09). This passed every
  # intervention row, so the frozen-fit table's residual columns were computed
  # on a population that included MIMIC TEST covariates. The residual is now
  # the frozen one and the rows are the ones it was fitted on.
  tar_target(lambda_fit, lambda_fit_diagnostics(
    priors$intervention, role = "final",
    ivf = tabs$intervention_features[
      tabs$intervention_features$stay_id %in% train_ids, , drop = FALSE],
    cfg = cfg)),
  tar_target(prior_stability, prior_fold_stability(priors)),

  # THE PRIOR-FIT STATUS POLICY, AS A TABLE. `layer1_priors()` already applied
  # it (strict) before the `priors` target could exist, so every row here is
  # PASS or WARN by construction; the target exists so the policy's verdict on
  # every (block, key, fold) is EXPORTED beside the diagnostics rather than
  # living only as the absence of an error. Added 2026-09-08 (plumbing review
  # F3): until then a Dirichlet fit that hit `max_iter` was recorded as
  # `converged = FALSE` in a column nothing read, and fitted 384 GAMs.
  tar_target(prior_fit_checks, check_prior_fits(priors, cfg, strict = FALSE)),

  # THE DIRICHLET-MULTINOMIAL'S OWN GOODNESS OF FIT. `dm_ppc_all()` has existed
  # since the prior diagnostics were written and was called by NOBODY in either
  # runner -- only by `tests/prior_fit.R`, which nobody has to run -- so the one
  # construct whose shrinkage survives had never had its fit checked at eICU.
  # Same class of gap as audit findings F1 and F4. Wired 2026-09-05.
  #
  # MEASURED at MIMIC: median dispersion 0.977 against a target of 1.000, median
  # PIT sd 0.999, median KS 0.0098. Five of 57 cells exceed |dispersion - 1| >
  # 0.5 and four of those are the structurally empty high tails, which sit at
  # the estimator floor by construction and enter no formula.
  # NAMED `dm_fit_ppc` AND NOT `dm_ppc`: `tar_source("R")` loads functions into
  # the same namespace `targets` resolves target names in, so a target named
  # `dm_ppc` SHADOWS the function `dm_ppc()` -- which `dm_ppc_all()` calls
  # internally. targets warns ("Ignoring global objects that conflict with
  # target names") rather than erroring, so the failure would have surfaced
  # inside the fit as a confusing type error. Found on the first
  # `tar_outdated()` after wiring it, 2026-09-05.
  tar_target(dm_fit_ppc, dm_ppc_all(priors$signal,
                                .measured_train_rows(tabs, folds, cfg),
                                role = "final", seed = cfg$seed)),
  # And whether that fit HOLDS across the coverage range, which is the question
  # `delta` and `lambda` both failed.
  tar_target(dm_drift, dm_dispersion_drift(priors$signal,
                                           .measured_train_rows(tabs, folds, cfg),
                                           cfg, role = "final")),

  # --- layer 1, out of fold -------------------------------------------------
  tar_target(jobs, layer1_jobs(cfg)),
  tar_target(budget, layer1_budget(cfg)),
  tar_target(fold_id, seq_len(cfg$n_folds %||% 5L)),

  tar_target(
    layer1_oof_branch,
    {
      design_checks
      jb <- jobs[jobs$fit & jobs$role == "oof" & jobs$fold == fold_id, , drop = FALSE]
      run_layer1(tabs, cfg, folds, priors, jobs = jb,
                 keep_final = FALSE, verbose = FALSE)
    },
    pattern   = map(fold_id),
    iteration = "list"
  ),

  tar_target(l_oof, do.call(rbind, lapply(layer1_oof_branch, function(z) z$l))),
  tar_target(diag_oof, do.call(rbind, lapply(layer1_oof_branch,
                                             function(z) z$diagnostics))),
  tar_target(triage_oof, diagnostics_summary(diag_oof, cfg)),

  # --- layer 1, NESTED, for the stacked xgb_l cell --------------------------
  # STATISTICAL REVIEW S1 (2026-09-09). `design_l` is out-of-fold row by row
  # but every row's L was fitted with the other four folds' outcomes, so a
  # booster cross-validated over it on the same folds trains on features that
  # saw the held-out fold. For every unordered pair of folds {a, b} the `full`
  # spec is refitted -- priors included -- on the three folds outside both and
  # predicted on both; `xgb_design_L_nested()` then assembles, for outer fold
  # a, training features from the pair fits and validation features from the
  # ordinary fold-a OOF fit. Ten pair fits per signal, 19 signals: ~190 GAMs
  # and ten final-only prior fits, kept out of the bundle (hard rule 8) and
  # out of every other arm. Branched per pair so a failed pair reruns alone.
  tar_target(nested_pairs, t(utils::combn(seq_len(cfg$n_folds %||% 5L), 2L))),
  tar_target(pair_id, seq_len(nrow(nested_pairs))),
  tar_target(
    layer1_nested_branch,
    {
      design_checks
      layer1_nested_l(tabs, cfg, folds, pair = nested_pairs[pair_id, ],
                      model = "full", verbose = FALSE)
    },
    pattern   = map(pair_id),
    iteration = "list"
  ),
  tar_target(l_nested, do.call(rbind, lapply(layer1_nested_branch, function(z) z$l))),
  tar_target(diag_nested, do.call(rbind, lapply(layer1_nested_branch,
                                                function(z) z$diagnostics))),

  # --- layer 1, final -------------------------------------------------------
  # The 43 fits the bundle carries. `$model` is stripped HERE rather than in
  # build_bundle(), so the fitted frames never enter the store either: every
  # diagnostic was already extracted by R/08 while the object was alive
  # (hard rule 6), and nothing downstream of this point reads them.
  tar_target(
    layer1_final_raw,
    {
      design_checks
      jb <- jobs[jobs$fit & jobs$role == "final", , drop = FALSE]
      run_layer1(tabs, cfg, folds, priors, jobs = jb,
                 keep_final = TRUE, verbose = FALSE)
    }
  ),
  tar_target(final_models, lapply(layer1_final_raw$models, strip_gam)),
  tar_target(diag_final, layer1_final_raw$diagnostics),

  # --- layer 2 --------------------------------------------------------------
  tar_target(l_mats_zero, l_matrices(l_oof, tabs, cfg, train_ids, fill = "zero")),
  tar_target(l_mats_na,   l_matrices(l_oof, tabs, cfg, train_ids, fill = "na")),

  # Sigma for EVERY L matrix, including the four interaction ones. It is a
  # correlation matrix over 19 columns and costs nothing next to a fit, and the
  # question it answers is a real one: the branch point asks whether the 19 L's
  # carry more than one direction of information, and an interaction model is
  # the most plausible way that answer could move. If PC1 rises under
  # `full_ti_all`, the cross terms are re-introducing a shared factor the
  # additive design had kept out.
  tar_target(sigma, {
    mats <- c(list(full = l_mats_zero$full, cond = l_mats_zero$cond,
                   intv = l_mats_zero$intv),
              l_mats_zero[c(LAYER1_TI_MODELS, sub("^full", "cond", LAYER1_TI_MODELS))])
    lapply(mats, l_correlation)
  }),
  tar_target(sigma_na, list(full = l_correlation(l_mats_na$full))),

  tar_target(eigen_full_zero, eigenspectrum(sigma$full)),
  tar_target(eigen_cond_zero, eigenspectrum(sigma$cond)),
  tar_target(eigen_full_na,   eigenspectrum(sigma_na$full)),
  tar_target(eigen_ti, do.call(rbind, lapply(
    c("full_ti_trend", "full_ti_all", "cond_ti_trend", "cond_ti_all"),
    function(nm) cbind(matrix_name = nm, eigenspectrum(sigma[[nm]]))))),

  # The branch point, both fills and both GCS treatments, in one table. The
  # interaction models are IN it as of 2026-09-07 rather than in a table of
  # their own, because the reading is the SPREAD across rows and a row kept
  # somewhere else is a row nobody compares.
  # `mats =` PASSES THE MATRICES THE GRAPH ALREADY HOLDS. Without it
  # `branch_point()` pivots `l_oof` twice more, reproducing `l_mats_zero` and
  # `l_mats_na` -- twelve pivots per run since the set went from four matrices
  # to six -- and, worse, computing "the L matrices" in a second place that
  # nothing held to the first.
  tar_target(branch_pt, branch_point(l_oof, tabs, cfg, train_ids,
    models = c("full", "cond", "intv", LAYER1_TI_MODELS,
               sub("^full", "cond", LAYER1_TI_MODELS)),
    mats = list(zero = l_mats_zero, na = l_mats_na))),
  tar_target(shared_pairs, shared_covariate_pairs(l_oof, tabs, cfg, train_ids)),
  tar_target(l_full_summary, l_summary(l_mats_zero$full)),
  tar_target(l_cond_summary, l_summary(l_mats_zero$cond)),
  tar_target(l_ti_summary, do.call(rbind, lapply(
    c("full_ti_trend", "full_ti_all", "cond_ti_trend", "cond_ti_all"),
    function(nm) cbind(matrix_name = nm, l_summary(l_mats_zero[[nm]]))))),
  # Per-signal discrimination, now on BOTH the joint and the conditional matrix
  # and now carrying AUPRC beside AUROC. At a 12% event rate AUROC alone is the
  # wrong instrument for "is one signal doing all the work": a column can look
  # respectable on AUROC and add nothing where the deaths are. `auprc_lift` is
  # AUPRC over the event rate, which is the only honest floor.
  tar_target(signal_auroc_full, signal_auroc(l_mats_zero$full, y_train)),
  tar_target(signal_auroc_cond, signal_auroc(l_mats_zero$cond, y_train)),
  # Per-signal, per L matrix, in ONE long table. The comparison worth making is
  # not "which signal leads" -- `signal_auroc_full` already says that -- but
  # whether adding the cross terms moves which signals carry the evidence, and
  # that is a difference between two rows of one table rather than two tables.
  tar_target(signal_auroc_ti, do.call(rbind, lapply(
    c("full_ti_trend", "full_ti_all", "cond_ti_trend", "cond_ti_all"),
    function(nm) cbind(matrix_name = nm,
                       signal_auroc(l_mats_zero[[nm]], y_train))))),

  # --- the comparison arms, out of fold -------------------------------------
  # All three cross-fitted on the SAME folds layer 1 used. `xgb_feat` goes
  # through xgb_oof_perfold() because its design carries fitted quantities
  # (delta, lambda) and one matrix cannot serve five folds without leaking.
  tar_target(design_l,   xgb_design_L(l_oof, tabs, cfg, train_ids,
                                      model = "full", fill = "zero")),
  tar_target(design_raw, xgb_design_raw(tabs, cfg, train_ids,
                                        missing_as_evidence = TRUE)),

  # `xgb_l` IS NESTED (review S1): one design per outer fold, built from the
  # pair fits above, through the same per-fold path `xgb_feat` uses. `design_l`
  # itself now serves only the FINAL booster, which is the documented stacking
  # design -- trained on OOF features, applied through full-train GAMs -- and
  # not an internal estimate.
  tar_target(oof_xgb_l, xgb_oof_perfold(
    function(f) xgb_design_L_nested(l_oof, l_nested, tabs, cfg, train_ids,
                                    fold_of_train, fold = f, model = "full"),
    y_train, fold_of_train, cfg, seed = cfg$seed, group = group_of_train)),
  tar_target(oof_xgb_raw, xgb_oof(design_raw, y_train, fold_of_train, cfg,
                                  seed = cfg$seed, group = group_of_train)),
  tar_target(oof_xgb_feat, xgb_oof_perfold(
    function(f) xgb_design_feat(tabs, cfg, priors, train_ids, role = "oof", fold = f),
    y_train, fold_of_train, cfg, seed = cfg$seed, group = group_of_train)),
  # The sizes every booster stage saw, per cell and fold (review S2): the
  # rows the round count was chosen on and the rows the refit was fitted on.
  tar_target(xgb_oof_sizes, do.call(rbind, list(
    cbind(design = "xgb_l",    oof_xgb_l$n_fit),
    cbind(design = "xgb_feat", oof_xgb_feat$n_fit),
    cbind(design = "xgb_raw",  oof_xgb_raw$n_fit)))),

  # Every arm's out-of-fold score, on the log-odds scale and centred on a
  # training prior, in one aligned matrix. This is what the reporting cut points
  # are frozen from and what `train_ref` summarises.
  #
  # EVERY ARM IN `BUNDLE_ARMS` IS REQUIRED, as of 2026-09-08 (plumbing review
  # F2). Until then `llr_meas` was built conditionally and a `Filter(Negate(
  # is.null))` dropped it when absent, so a run whose `meas` assembly had
  # failed upstream produced nine arms, nine sets of cut points and a bundle
  # that passed its checks -- and the apply sites then requested ten. The LLR
  # arms are read through `LLR_ARM_MATRIX`, the SAME map `apply_bundle()`
  # reads, so the two sites cannot disagree about which matrix an arm sums.
  # `llr_cond` (2026-09-05) and the four interaction arms (2026-09-07) are in
  # that map; the order below is the order the arms have always had here, so
  # the stored value is unchanged.
  tar_target(oof_scores, {
    s <- lapply(LLR_ARM_MATRIX, function(mt) {
      M <- l_mats_zero[[mt]]
      if (is.null(M)) {
        stop("oof_scores: no `", mt, "` L matrix. Every arm in BUNDLE_ARMS is ",
             "required; an absent one must fail here rather than shorten the ",
             "cut-point list and the training reference.", call. = FALSE)
      }
      rowSums(M)
    })
    s$xgb_l    <- oof_xgb_l$score
    s$xgb_feat <- oof_xgb_feat$score
    s$xgb_raw  <- oof_xgb_raw$score
    if (!identical(names(s), BUNDLE_ARMS)) {
      abort_values("oof_scores: arm set is not BUNDLE_ARMS, in order", names(s))
    }
    for (nm in names(s)) {
      if (is.null(s[[nm]]) || length(s[[nm]]) != length(train_ids)) {
        stop("oof_scores: arm `", nm, "` is absent or not aligned to train_ids", call. = FALSE)
      }
    }
    s
  }),

  # --- frozen reporting bins ------------------------------------------------
  # Computed on the OUT-OF-FOLD scores, not the final-fit ones. An apply site is
  # scored by models that never saw it, so the training distribution it is being
  # compared against must be the one that also never saw its own rows.
  tar_target(cutpoints, lapply(oof_scores, score_cutpoints,
                               n_bins = cfg$metrics$n_bins %||% 20L)),

  # --- training reference numbers -------------------------------------------
  # AUROC/AUPRC per arm on MIMIC train, out of fold. Carried in the bundle so a
  # transport table can name its training-side comparator without re-reading the
  # internal run. No files written here: score_report() needs a run object and
  # runs are the export layer, not the graph.
  tar_target(train_ref, {
    nb <- cfg$metrics$n_boot %||% 200L
    nbins <- cfg$metrics$n_bins %||% 20L
    summ <- do.call(rbind, lapply(names(oof_scores), function(nm) {
      s  <- oof_scores[[nm]]
      m  <- score_metrics(s, y_train, label = nm, n_boot = nb, seed = cfg$seed,
                          group = group_of_train)
      # score_report() adds these by way of a run object, and the graph has
      # none (hard rule 9). They are computed directly here because the
      # transport table compares them across sites: AUROC is rank-based and
      # cannot see a calibration shift, so `cal_slope` is the only column in
      # this table that can.
      mo <- monotonicity(risk_bins(s, y_train, n_bins = nbins))
      ca <- llr_calibration(s, y_train, p_bar_cohort_train)
      m$spearman   <- mo$spearman
      m$rate_ratio <- mo$rate_ratio
      m$cal_slope  <- ca$slope
      m
    }))
    list(site = "mimic", stage = "train_oof",
         n = length(train_ids), n_events = sum(y_train),
         p_bar_cohort = p_bar_cohort_train,
         arms = summ)
  }),

  # The same table, flat, so it can be EXPORTED AS A CSV. `train_ref` is a
  # nested list and `export_run()` routes it to `objects`, which writes an .rds
  # and no CSV -- so the internal run's arm-level AUROC and AUPRC existed only
  # inside an .rds, and the only CSV in the run directory with an AUROC column
  # was the per-signal one, which carried no AUPRC at all. Added 2026-09-05.
  tar_target(train_ref_arms, {
    d <- train_ref$arms
    cbind(site = train_ref$site, stage = train_ref$stage,
          p_bar_cohort = round(train_ref$p_bar_cohort, 6), d,
          stringsAsFactors = FALSE)
  }),

  # --- the frozen boosters --------------------------------------------------
  # One booster per design on ALL training rows. This is what an apply site is
  # scored with (hard rule 8); the out-of-fold cells above exist only to give
  # the training site its own honest comparator.
  tar_target(xgb_final, list(
    xgb_l    = xgb_fit_full(design_l,   y_train, cfg, seed = cfg$seed,
                            group = group_of_train),
    xgb_raw  = xgb_fit_full(design_raw, y_train, cfg, seed = cfg$seed,
                            group = group_of_train),
    xgb_feat = xgb_fit_full(
      xgb_design_feat(tabs, cfg, priors, train_ids, role = "final"),
      y_train, cfg, seed = cfg$seed, group = group_of_train)
  )),

  tar_target(feature_manifest, do.call(rbind, lapply(names(xgb_final), function(nm)
    xgb_feature_table(xgb_final[[nm]], nm)))),
  tar_target(group_gain, do.call(rbind, lapply(split(feature_manifest,
                                                     feature_manifest$design),
                                               xgb_group_gain))),


  # --- the severity comparators, frozen -------------------------------------
  # APACHE II and SOFA are recomputed from frozen point tables (R/09c, R/09d)
  # and nothing about that construction is fitted, so it needs no bundle. ONE
  # THING IS FITTED: the logistic recalibration that puts an integer score on
  # the log-odds scale, which `recalibrate_oof()` does fold-wise here and which
  # an apply site -- having no folds -- cannot re-derive without fitting on its
  # own outcomes. The transportable object is the FULL-TRAIN intercept and
  # slope, and that is what the bundle carries (hard rule 8).
  #
  # The coverage restriction is part of the frozen object too, and it is applied
  # to EVERY cell rather than to the baselines alone: comparing a restricted
  # baseline against an unrestricted proposed method is a different and much
  # weaker experiment.
  tar_target(severity_train, severity_raw(tabs, cfg, train_ids)),
  tar_target(severity_diag,  severity_diagnostics(severity_train)),

  # What travels. Full train, restricted rows.
  tar_target(severity_recal, fit_severity_recal(severity_train, y_train)),

  # What the TRAINING SITE reports for itself. Out-of-fold, so it is the honest
  # comparator and is built the same way every other `train_ref` arm is. This is
  # never what an apply site uses -- that is `severity_recal` above -- and the
  # two differing slightly is expected rather than a defect.
  tar_target(severity_oof, {
    k  <- severity_train$keep
    yk <- y_train[k]
    fk <- fold_of_train[k]
    pb <- mean(yk)
    stats::setNames(lapply(SEVERITY_CELLS, function(nm)
      recalibrate_oof(severity_train$raw[[nm]][k], yk, fk, pb)), SEVERITY_CELLS)
  }),

  # TIE-AWARE (review S11): a recalibrated integer score can tie at a quantile
  # boundary on valid data, and the strict mode would abort the bundle after
  # every fit had finished. The frozen vector carries its actual bin count.
  tar_target(severity_cutpoints, lapply(severity_oof, score_cutpoints,
                                        n_bins = cfg$metrics$n_bins %||% 20L,
                                        tie_aware = TRUE)),

  # Every arm on the RESTRICTED rows, so the severity comparison and the bundle
  # comparison are never quoted off different row sets. `llr_sum` shifting
  # between this table and `train_ref` is the size of the selection the coverage
  # floors introduce, and it is a number to report rather than to assume away.
  tar_target(severity_train_ref, {
    k  <- severity_train$keep
    yk <- y_train[k]
    pb <- mean(yk)
    nb <- cfg$metrics$n_boot %||% 200L
    cells <- c(lapply(oof_scores, function(v) v[k]), severity_oof)
    gk <- group_of_train[k]
    summ <- do.call(rbind, lapply(names(cells), function(nm)
      score_metrics(cells[[nm]], yk, label = nm, n_boot = nb, seed = cfg$seed,
                    group = gk)))
    list(site = "mimic", stage = "train_oof_restricted",
         n = sum(k), n_events = sum(yk), p_bar = pb,
         frac_kept = mean(k), arms = summ)
  }),

  tar_target(severity_slot, make_severity_bundle(
    settings  = severity_train$settings,
    recal     = severity_recal,
    p_bar     = attr(severity_recal, "p_bar"),
    n_train   = sum(severity_train$keep),
    frac_kept = mean(severity_train$keep),
    cutpoints = severity_cutpoints,
    train_ref = severity_train_ref,
    # Frozen with the rest of the severity design, so an apply site cannot
    # choose a different comparison population than the training site did.
    symmetric = cfg$severity_symmetric)),

  # The post-loader schema of the training tables, frozen so an apply site
  # checks its columns, classes and factor levels against what the models were
  # FITTED ON rather than against the MIMIC extraction on disk at apply time
  # (external runner review E7, 2026-09-09). Names, classes and level sets
  # only; no row value enters it.
  tar_target(schema_sig, schema_signature(tabs)),

  # --- the bundle -----------------------------------------------------------
  # Everything above that a non-fitting site cannot re-derive, in one object,
  # verified before it leaves the graph. `layer2` is NULL: the Sigma-inverse
  # weights do not exist yet, and a placeholder vector of ones would let a later
  # apply site produce an unweighted aggregate under a weighted name.
  tar_target(bundle, build_bundle(
    cfg        = cfg,
    priors     = priors,
    models     = final_models,
    sigma      = sigma,
    eigen      = list(full_zero = eigen_full_zero,
                      cond_zero = eigen_cond_zero,
                      full_na   = eigen_full_na),
    cutpoints  = cutpoints,
    xgb        = xgb_final,
    layer2     = NULL,
    train_ref  = train_ref,
    severity   = severity_slot,
    schema     = schema_sig,
    domains    = domains,
    site       = "mimic",
    # Provenance as a DECLARED dependency (audit finding F3), not a read from
    # inside the target.
    source_hashes = code_hashes)),

  tar_target(bundle_checks, verify_bundle(bundle, cfg = cfg, strict = TRUE)),
  tar_target(bundle_contents, bundle_summary(bundle))
)
