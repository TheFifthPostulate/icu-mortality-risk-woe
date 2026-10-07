# tests/attr_external_common.R -------------------------------------------------
# THE IDENTITY OF AN eICU REPRODUCIBILITY STORE, defined once and sourced by
# both its generator (`tests/attr_external_bags.R`) and its consumer
# (`tests/attr_metrics.R --site eicu`). It lives under tests/ and not under R/
# on purpose: a new file under R/ changes the `source_files` target and
# rebuilds the bundle hash chain (`_targets.R`), and nothing here belongs to
# the pipeline. Added 2026-10-06 (paper/plans/plan_eicu_reproducibility_arm.md).
#
# WHAT THE STORE IS. The internal attribution run retrained every method on 38
# shared bags with the out-of-fold cross-validation protocol and compared the
# replicates on the MIMIC-IV training stays. This store retrains the paper's
# three methods WHOLE on the same 38 bags (no folds), with the random-variable
# parameters held at the FINAL bundle values, and applies each retrained model
# to the eICU cohort. Hard rule 9: these are perturbation fits on training
# data, for reproducibility only; none of them becomes a scoring bundle, and
# nothing is fitted on eICU outcomes.
# ------------------------------------------------------------------------------

EICU_GEN_VERSION <- "eicu_wholebag_v1"
EICU_METHODS     <- c("llr_meas", "llr_full", "shap_xgb_feat")

#' The design key and the fingerprint of an eICU whole-bag store.
#'
#' The design key names the files: it carries the internal store's design key
#' (so the bags and the specifications are those of the internal run), the
#' frozen bundle's design, the eICU stay ids and the method set. The
#' fingerprint adds what the key does not cover: the frozen priors, the frozen
#' model coefficients and booster, the booster settings and protocol, and the
#' evaluation seeds and counts. Both are recomputed by the consumer.
eicu_store_identity <- function(int_design_key, bundle, cfg, ecfg, stay_ids,
                                methods = EICU_METHODS) {
  design <- list(
    site            = "eicu",
    fit             = EICU_GEN_VERSION,
    internal_design = attr_key_hash(int_design_key),
    bundle_design   = attr_key_hash(bundle$cfg),
    eicu_ids        = attr_key_hash(as.character(stay_ids)),
    methods         = sort(as.character(methods)))
  mk <- sort(names(bundle$models))
  fp <- list(
    design        = attr_key_hash(design),
    bundle_cfg    = attr_key_hash(bundle$cfg),
    priors        = attr_key_hash(bundle$priors),
    model_coefs   = attr_key_hash(lapply(bundle$models[mk], stats::coef)),
    booster       = attr_key_hash(xgboost::xgb.save.raw(bundle$xgb$xgb_feat$booster)),
    xgboost       = cfg$xgboost,
    xgb_fit_protocol = XGB_FIT_PROTOCOL,
    bootstrap_seed_base = as.integer(cfg_req(ecfg, "bootstrap", "seed_base")),
    bootstrap_b   = as.integer(cfg_req(ecfg, "bootstrap", "b")),
    seed_base     = as.integer(cfg_req(ecfg, "levels", "seed", "seed_base")),
    seed_b        = as.integer(cfg_req(ecfg, "levels", "seed", "b")),
    seed_b_per_bag = as.integer(cfg_req(ecfg, "levels", "seed", "b_per_bag")),
    draw_seed_base = as.integer(cfg_req(ecfg, "levels", "sample", "draw_seed_base")),
    draw_b        = as.integer(cfg_req(ecfg, "levels", "sample", "b")),
    generator     = EICU_GEN_VERSION)
  list(design_key = design, fingerprint = fp)
}

#' The replicate plan of an eICU store: the internal plan restricted to the
#' store's methods. The routes, counts and coordinates are those of
#' `attr_replicate_plan()`, so the consumer forms the same contrasts.
eicu_store_plan <- function(ecfg, n_folds, methods = EICU_METHODS) {
  e <- ecfg
  e$methods <- as.character(methods)
  attr_replicate_plan(e, n_folds = n_folds, n_pass_fits = 0L)
}
