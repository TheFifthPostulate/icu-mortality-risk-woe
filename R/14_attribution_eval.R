# R/14_attribution_eval.R -----------------------------------------------------
# THE ATTRIBUTION-EVALUATION LIBRARY. One vocabulary, one design fingerprint,
# one metric family.
#
# Read `docs/attribution_analysis_plan_20260906.md` section 10 first. This file
# is step 3 of its order of work, and it is CONSOLIDATION rather than new
# statistics: every metric here was already written, once each, inside
# `tests/attribution_ties.R`, `tests/coupling_attribution.R` and
# `tests/shap_noise_floor.R`. Three copies of `prep()` and `agree_k()` existed
# and were byte-identical by luck rather than by construction, which is exactly
# the arrangement in which one of them silently acquires a different convention.
# `tests/attr_metrics.R` asserts that this file reproduces the scattered tables
# on the identical cached inputs, and that assertion is the correctness gate for
# the whole sub-pipeline.
#
# WHY GENERATION AND MEASUREMENT ARE SPLIT. Every script named above builds its
# own inputs and then measures them, so a change to a metric costs a refit and a
# change to a specification costs a rewrite. Here the unit is an ATTRIBUTION
# REPLICATE
#
#     (method, boot_id, seed_id, draw_id)  ->  an N x G matrix in nats, oof
#
# and every metric is a function of a PAIR of replicates. The pair's place in
# the reproducibility hierarchy is decided by WHICH COORDINATE DIFFERS, by
# `attr_pair_contrast()` and nowhere else, so no caller can label a
# specification comparison as a noise measurement.
#
# THE COORDINATE REPLACED A `(method, level, index)` TRIPLE ON 2026-09-06, and
# it had to. `level` was the generator's ROUTE and `index` was one integer
# standing in for a bag id, a seed id or a posterior draw id depending on which
# loop emitted the row, so a replicate that differs from another in TWO ways --
# a re-seeded fit on a resampled bag -- could not be described at all. That
# replicate is the one that carries the level-2-into-level-3 propagation.
#
# AND EVERY LEVEL IS NOW A DISTRIBUTION, NOT A NUMBER. Until 2026-09-06 each
# level was reported from replicate 1 against replicate 2 and the remaining
# pairs were generated and discarded. `attr_stat_summary()` fixes the column
# set -- min, p05, p10, q1, median, mean, q3, p90, p95, max, IQR, sd and two
# skew indicators -- so every level is reported the same way and two can be
# read side by side. `attr_dominance()` compares two of those distributions
# without dividing one by the other.
#
# HARD RULE 9. No path is built here, no clock is read, and nothing touches the
# filesystem. `tests/attr_replicates.R` owns every path; this file receives
# matrices and returns tables.
#
# HARD RULE 1. Two functions return ROW-LEVEL vectors -- `attr_cosine()` returns
# one similarity per patient and `attr_displacement()` returns a sample of
# per-cell displacements. Both are inputs to a summary, never to a print;
# `attr_dist_summary()` is the intended consumer and returns quantiles.
# ----------------------------------------------------------------------------

# --- the vocabulary ----------------------------------------------------------
#
# THE ARM NAMES ARE THE LADDER'S, NOT `R/13_attribution.R`'s. The two files
# answer different questions and their method vocabularies genuinely differ.
# `ATTRIBUTION_METHODS` in R/13 names the arms an APPLY SITE can produce from a
# bundle -- `llr_sum`, `llr_cond`, `llr_meas` and the three boosters -- because
# that arm is about transport. The ladder here additionally carries the two
# INTERACTION specifications, which no bundle contains and which exist only
# out-of-fold at the training site. Collapsing the two vocabularies would
# require either inventing bundle entries for arms that were never frozen, or
# dropping the interaction arms, and both are worse than two named sets with a
# stated relationship: `llr_full` here IS `llr_sum` there (before layer-2
# weights), and `llr_cond` and `llr_meas` are the same object in both.

#' The five arms that must be GENERATED. Everything else is a subtraction.
ATTR_BASE_ARMS <- c("meas", "intv", "full", "full_ti_trend", "full_ti_all")

#' The three that are derived. `cond = full - intv` is the attribution unit:
#' `full` bundles the measurement evidence with a pure treatment-propensity
#' contrast, so "signal g contributed X nats" is, on `full`, partly a claim
#' about who got treated.
ATTR_DERIVED_ARMS <- c("cond", "cond_ti_trend", "cond_ti_all")

ATTR_LLR_ARMS  <- c(ATTR_BASE_ARMS, ATTR_DERIVED_ARMS)
ATTR_SHAP_ARMS <- c("shap_xgb_feat", "shap_xgb_raw")

#' Every method this library will accept.
ATTR_EVAL_METHODS <- c(paste0("llr_", ATTR_LLR_ARMS), ATTR_SHAP_ARMS)

#' What is being varied between two replicates of the same method.
#'
#' `spec` is the degenerate level at which a method has exactly one replicate:
#' the arm as fitted. It exists so that a level-4 comparison (different method)
#' has a well-formed level to sit at, rather than being a special case.
ATTR_LEVELS <- c("spec", "seed", "sample")

#' Method family. `llr` arms are deterministic under re-seeding; `shap` arms are
#' not, and that asymmetry is a reportable property rather than a gap.
attr_method_family <- function(method) {
  ifelse(method %in% ATTR_SHAP_ARMS, "shap",
         ifelse(method %in% paste0("llr_", ATTR_LLR_ARMS), "llr", NA_character_))
}

#' The ladder arm name behind a method name. `llr_cond_ti_all` -> `cond_ti_all`.
attr_arm_of <- function(method) {
  ifelse(method %in% ATTR_SHAP_ARMS, method, sub("^llr_", "", method))
}

attr_check_methods <- function(methods) {
  bad <- setdiff(methods, ATTR_EVAL_METHODS)
  if (length(bad)) abort_values("unknown attribution method(s)", bad)
  invisible(methods)
}

#' Does this method have any seed-level variation at all?
#'
#' `bam` with fREML is deterministic at fixed settings -- audit finding F9
#' recorded a bitwise-identical refit -- and the only stochastic element in the
#' GAM path is fold assignment, which is seeded and frozen. So LEVEL 2 IS
#' EXACTLY ZERO FOR EVERY LLR ARM, and the generator must refuse to spend
#' compute pretending otherwise. It is also why level 3 is the only common
#' currency between the LLR arms and SHAP: for SHAP most of the instability is
#' already at level 2, for the LLR arms it is all at level 4, and the two are
#' not comparable until both are expressed against sampling noise.
attr_has_seed_noise <- function(method) {
  identical(unname(attr_method_family(method)), "shap")
}

# --- deriving the eight LLR arms from the five that are fitted ---------------

#' `cond` arms from the base arms, by subtraction. Costs no fit.
#'
#' @param base named list of matrices over ATTR_BASE_ARMS, identically shaped.
#' @return the same list plus the three derived arms.
attr_derive_arms <- function(base) {
  miss <- setdiff(ATTR_BASE_ARMS, names(base))
  if (length(miss)) abort_values("attr_derive_arms: missing base arm(s)", miss)
  d <- dim(base[[ATTR_BASE_ARMS[1]]])
  for (a in ATTR_BASE_ARMS) {
    if (!identical(dim(base[[a]]), d)) {
      stop("attr_derive_arms: arm `", a, "` is not shaped like `",
           ATTR_BASE_ARMS[1], "`", call. = FALSE)
    }
  }
  base$cond          <- base$full          - base$intv
  base$cond_ti_trend <- base$full_ti_trend - base$intv
  base$cond_ti_all   <- base$full_ti_all   - base$intv
  base[ATTR_LLR_ARMS]
}

# --- the interaction term set -------------------------------------------------
#
# `attr_ti_terms()` AND `attr_smooth_vars()` ARE GONE FROM THIS FILE as of
# 2026-09-07. They are `interaction_terms()` and `smooth_term_vars()` in
# R/05_formula.R, which emit byte-identical strings.
#
# THE MOVE IS THE SAME MOVE, MADE TWICE, FOR THE SAME REASON. On 2026-09-06 the
# term builder moved out of `tests/coupling_attribution.R` and into this file,
# because the ladder's design fingerprint had acquired a second reader and a
# fingerprint computed by two copies of a function is not a fingerprint. On
# 2026-09-07 `full_ti_trend` and `full_ti_all` joined `LAYER1_MODELS`, so the
# cross terms became something the PIPELINE fits, and a term set defined in the
# attribution library and re-derived by `build_formula()` would have been two
# definitions of one design -- the identical failure, one level up.
#
# `attr_design_key()` below still hashes those strings, and they have not moved,
# which is what lets the 371 replicates on disk keep their `ti_all` / `ti_trend`
# fields through the spec change.

# --- the design fingerprint --------------------------------------------------

#' Short content hash. Pure computation; builds no path and reads no clock.
attr_key_hash <- function(x) substr(digest::digest(x, algo = "md5"), 1L, 12L)

#' Everything that decides what an L ladder CONTAINS.
#'
#' WHY THIS EXISTS. On 2026-09-06 `tests/attribution_vs_shap.R` was run against
#' a ladder built before `bam.gamma` changed from 1 to 1.5. Every SHAP row in
#' that run is current -- `xgb_feat`'s design is built from `pi_hat`, `delta`
#' and `lambda`, none of which `bam` fits -- and every L row is stale. Nothing
#' errored, nothing looked wrong, and the two halves of the table were simply
#' from different pipelines. `tests/coupling_attribution.R` already refused a
#' stale cache through a fingerprint of its own; the consumer had no such check.
#'
#' The fields are chosen so that ANY edit which would change a fitted L changes
#' the key:
#'   `bam`       carries `gamma`, `method`, `discrete`, `select`, `nthreads` and
#'               the default basis. `nthreads` is in there because audit finding
#'               F10 measured it moving fitted L values at the 1e-13 level.
#'   `formulas`  every fitted layer-1 spec, deparsed. Catches a whitelist edit,
#'               a `smooth_k` override, a level-term rule and a pairing change.
#'
#'               IT IS OVER-BROAD IN ONE DIRECTION, and 2026-09-07 is when that
#'               mattered. Adding the two interaction models to `LAYER1_MODELS`
#'               put 21 new specs into `layer1_jobs()`, so this field changed --
#'               while every EXISTING spec's formula, rows, folds and `bam`
#'               settings were untouched and every stored L was bitwise
#'               unchanged. An independent new spec cannot move another spec's
#'               L, so the field reports a difference that is real about the
#'               DESIGN and false about the DATA. The remedy is a migration with
#'               evidence attached, not a regeneration: verify the cheap
#'               replicates bitwise, then rewrite the stored key. Narrowing the
#'               field to "the formulas of the specs this store actually holds"
#'               was considered and rejected -- it would stop catching a spec
#'               being REMOVED, which is a change that does invalidate a store.
#'   `ti`        both interaction scopes' cross terms, and `k_ti`.
#'   `folds`     hashed rather than stored: fold membership is per stay and
#'               therefore row-level (hard rule 1), while its hash is not.
#'   `signals`   the modelled vocabulary, so a signal entering or leaving is
#'               visible even if no surviving formula changed.
#'
#' What is deliberately NOT in it: paths, run ids, package versions, the `R/`
#' source hashes. Those belong in the run manifest, which records them already.
#' A key that changes when a comment changes would refuse every cache and teach
#' the reader to pass `--force`.
#'
#' @param fold_vec per-stay fold assignment for the training rows, in the order
#'   the ladder's rows are stored. Hashed here, never returned.
#' @param k_ti the interaction basis cap. Defaults to the config value, which is
#'   what the pipeline now fits at; passed explicitly by the two staleness
#'   guards, which have to key a STORED object to the value it was built under
#'   rather than to the value in force today.
#' @param folds_hash an already-computed `folds` field, for a caller that has
#'   no fold vector -- an apply site checking a training-site store against
#'   the bundle's frozen design (review finding A10). Every OTHER field is then
#'   recomputed from `cfg` and compared; the folds field is carried through
#'   unchanged and the caller must say so in whatever it reports.
attr_design_key <- function(cfg, fold_vec, k_ti = cfg_req(cfg, "k_ti"),
                            folds_hash = NULL) {
  sigs   <- as.character(unlist(cfg$signals))
  paired <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, sigs)
  jobs   <- layer1_jobs(cfg)
  specs  <- unique(paste(jobs$signal, jobs$model)[jobs$fit])
  if (is.null(folds_hash)) {
    if (is.null(fold_vec)) {
      stop("attr_design_key: either `fold_vec` or `folds_hash` is required",
           call. = FALSE)
    }
    folds_hash <- attr_key_hash(as.character(fold_vec))
  }
  list(
    bam = bam_settings_row(cfg),
    signals = paste(sigs, collapse = ","),
    paired  = paste(paired, collapse = ","),
    formulas = vapply(specs, function(s) {
      p <- strsplit(s, " ", fixed = TRUE)[[1]]
      paste(deparse(build_formula(p[1], p[2], cfg)), collapse = " ")
    }, character(1)),
    ti_all   = vapply(paired, function(sg)
      paste(interaction_terms(sg, cfg, "all", k_ti), collapse = " + "), character(1)),
    ti_trend = vapply(paired, function(sg)
      paste(interaction_terms(sg, cfg, "trend", k_ti), collapse = " + "), character(1)),
    k_ti = as.integer(k_ti),
    folds = as.character(folds_hash))
}

# --- THE STORE FINGERPRINT: WHAT THE DESIGN KEY DOES NOT COVER --------------
#
# ADDED 2026-09-09, REVIEW FINDINGS A1 AND A2. `attr_design_key()` decides the
# FILE NAME a replicate is stored under, and it was also the only thing the
# generator's resume and the consumer's staleness guard compared. It covers the
# GAM design and nothing else: a change to `xgboost.max_depth`, to the
# vasopressor agent pool that shapes intensity covariates, to the frozen priors
# the covariates are built from, to the posterior draw seed, or to the number
# of posterior draws (which changes EVERY earlier draw -- see
# `ATTR_GENERATOR_VERSION`) left the key identical, so a resume would keep old
# replicates under a design that no longer produces them and the consumer
# would read them without complaint. Reproduced synthetically in
# `tests/attr_eval_unit.R`.
#
# THE FINGERPRINT IS A SECOND, WIDER IDENTITY CHECKED BESIDE THE KEY, NOT A
# NEW KEY. Widening the key itself would rename every file in the store and
# turn a validation fix into a five-hour regeneration; the 2026-09-07
# migration already showed that a key change with no data change is a
# migration with evidence attached, not a rebuild. So the file names stay, and
# `design.qs2` additionally carries this fingerprint. The generator refuses to
# resume into a store whose fingerprint differs, the consumer refuses to read
# one, and `attr_fingerprint_diff()` names the field AND the routes that
# field can move, so the operator knows which replicates are actually stale.
#
# A store written before the fingerprint existed has none. It is stamped by
# the generator on its next resume, AFTER the LLR ladder replicates on disk
# are verified bitwise against the targets cache (free) and the SHAP ladder is
# refitted and verified (about a minute and a half). The evidence is written
# beside the stamp; a stamp with no evidence is a `--force` wearing a hat.

#' Identity of the generator's own algorithms, for the fingerprint.
#'
#' `posterior_draw_layout` names the fact that `rmvn_clamped()` fills its
#' standard-normal matrix column-major over `nd * p` values, so draw i depends
#' on the TOTAL number of draws requested: extending B = 40 to B = 60 changes
#' the first 40 draws while their coordinates and file names stay the same
#' (finding A2). A stable-under-extension layout is the right long-term fix;
#' switching to it changes every posterior replicate on disk, so it is a
#' declared version bump that forces posterior regeneration, not an edit.
ATTR_GENERATOR_VERSION <- list(
  posterior_draw_layout = "rnorm_colmajor_nd_dependent_v1",
  replicate_key         = "coord_key_v1_20260906")

#' Which routes a fingerprint field can invalidate. Reporting only.
ATTR_FINGERPRINT_ROUTES <- c(
  design     = "all",
  cfg_design = "all",
  xgboost    = "ladder(shap),seed,bootstrap(shap),bootstrap_seeded",
  priors     = "all",
  eval_bootstrap_seed_base = "bootstrap,bootstrap_seeded",
  eval_seed_base           = "seed,bootstrap_seeded",
  eval_draw_seed_base      = "posterior",
  eval_sample_b            = "posterior",
  eval_llr_route           = "posterior",
  generator  = "posterior",
  # The booster fitting PROTOCOL, not its hyperparameters: `xgboost` above is
  # the config block, which did not change when `.xgb_fit1()` moved to a
  # grouped inner split and a full-row refit (statistical review S2/S3,
  # 2026-09-09). Every SHAP replicate depends on the protocol.
  xgb_fit_protocol = "ladder(shap),seed(shap),bootstrap(shap),bootstrap_seeded(shap)")

#' Which replicates a fingerprint mismatch actually invalidates.
#'
#' ROUTE-SELECTIVE, AS OF 2026-09-09. A fingerprint difference used to stop
#' the generator and the consumer outright, whatever field moved: after the
#' booster protocol changed, the five-hour LLR bootstrap replicates -- valid by
#' estimand -- were unreadable until the SHAP routes had been regenerated into
#' a new store. `ATTR_FINGERPRINT_ROUTES` already names, per field, the routes
#' that field can move and the method family it moves them for; this parses
#' those declarations so the two scripts can quarantine exactly the stale
#' replicates and keep the rest.
#'
#' A declaration is `route` or `route(family)`; `all` means every replicate.
#' Parsing is strict: an unknown route or family is an error, because a
#' declaration that fails to parse must not silently invalidate nothing.
#'
#' @param fd the diff table from `attr_fingerprint_diff()`
#' @return list(all = logical, rules = data.frame(route, family)) where
#'   `family` is "shap", "llr" or "any"
attr_stale_routes <- function(fd) {
  routes <- c("ladder", "seed", "bootstrap", "bootstrap_seeded", "posterior")
  if (!nrow(fd)) return(list(all = FALSE, rules = data.frame(route = character(0),
                                                             family = character(0),
                                                             stringsAsFactors = FALSE)))
  toks <- unique(trimws(unlist(strsplit(as.character(fd$routes_affected), ","))))
  if (any(toks %in% c("all", "unknown_field"))) {
    return(list(all = TRUE, rules = data.frame(route = routes, family = "any",
                                               stringsAsFactors = FALSE)))
  }
  m <- regmatches(toks, regexec("^([a-z_]+)(?:\\(([a-z]+)\\))?$", toks, perl = TRUE))
  bad <- toks[vapply(m, length, integer(1)) == 0L]
  if (length(bad)) abort_values("attr_stale_routes: unparseable route declaration", bad)
  rules <- data.frame(
    route  = vapply(m, function(z) z[2], character(1)),
    family = vapply(m, function(z) if (nzchar(z[3])) z[3] else "any", character(1)),
    stringsAsFactors = FALSE)
  if (!all(rules$route %in% routes)) {
    abort_values("attr_stale_routes: unknown route", setdiff(rules$route, routes))
  }
  if (!all(rules$family %in% c("any", "shap", "llr"))) {
    abort_values("attr_stale_routes: unknown method family", setdiff(rules$family, c("any", "shap", "llr")))
  }
  list(all = FALSE, rules = unique(rules))
}

#' The manifest rows a set of stale-route rules names. Logical over `man`.
attr_stale_mask <- function(man, rules) {
  if (!nrow(man)) return(logical(0))
  fam <- attr_method_family(man$method)
  out <- rep(FALSE, nrow(man))
  for (i in seq_len(nrow(rules))) {
    hit <- man$route == rules$route[i] &
      (rules$family[i] == "any" | (!is.na(fam) & fam == rules$family[i]))
    out <- out | hit
  }
  out
}

#' The wider identity of a replicate store.
#'
#' @param cfg the fitted-design config (targets `cfg`).
#' @param fold_vec,k_ti as `attr_design_key()`.
#' @param priors the `llr_priors` container the covariates are built from. Only
#'   its hash is kept: the object is per (signal, fold) parameters, never rows.
#' @param eval_cfg parsed `config/attribution_eval.yml`.
#' @param design_key an already-computed `attr_design_key()`, to avoid
#'   rebuilding 64 formulas twice.
attr_design_fingerprint <- function(cfg, fold_vec, k_ti, priors, eval_cfg,
                                    design_key = NULL) {
  dk <- if (is.null(design_key)) attr_design_key(cfg, fold_vec, k_ti) else design_key
  # THE SAME SUBSET THE BUNDLE FREEZES. `BUNDLE_DESIGN_KEYS` is documented as
  # "exactly the set of cfg references in R/", so it is the sanctioned
  # statement of what counts as design; `paths` is local by construction and
  # is the one top-level key it excludes.
  cd <- cfg[intersect(BUNDLE_DESIGN_KEYS, names(cfg))]
  list(
    design     = attr_key_hash(dk),
    cfg_design = attr_key_hash(cd),
    xgboost    = cfg$xgboost,
    priors     = attr_key_hash(priors),
    eval_bootstrap_seed_base = as.integer(cfg_req(eval_cfg, "bootstrap", "seed_base")),
    eval_seed_base           = as.integer(cfg_req(eval_cfg, "levels", "seed", "seed_base")),
    eval_draw_seed_base      = as.integer(cfg_req(eval_cfg, "levels", "sample", "draw_seed_base")),
    eval_sample_b            = as.integer(cfg_req(eval_cfg, "levels", "sample", "b")),
    eval_llr_route           = as.character(cfg_req(eval_cfg, "levels", "sample", "llr_route")),
    generator  = ATTR_GENERATOR_VERSION,
    xgb_fit_protocol = XGB_FIT_PROTOCOL)
}

#' Which fingerprint fields differ, with the routes each can invalidate.
attr_fingerprint_diff <- function(f1, f2) {
  d <- attr_design_diff(f1, f2)
  d$routes_affected <- unname(ATTR_FINGERPRINT_ROUTES[d$field])
  d$routes_affected[is.na(d$routes_affected)] <- "unknown_field"
  d
}

#' What the level-3 replicates actually vary, stated once and emitted by both
#' scripts (review finding A5). The generator's header explains each line at
#' length; this is the version that travels inside a run directory, so a table
#' is never read without the estimand it was computed under.
ATTR_ESTIMAND_NOTES <- data.frame(
  field = c("resample_kind", "resample_unit", "in_bag_fraction",
            "held_fixed_under_resampling", "level3_estimand",
            "level3_posterior_estimand", "shap_matrix"),
  value = c(
    "distinct members of a with-replacement patient draw (an m-out-of-n subsample without multiplicity), NOT a multiplicity-preserving bootstrap",
    "patient (fold_group_by), so repeat admissions stay together",
    "about 0.632 of patients per bag; see diagnostics/bootstrap_bags",
    "out-of-fold priors (alpha, delta, lambda), fold assignment, covariate construction, smooth_k; layer-1 GAMs and boosters are refitted, transformations are not",
    "sampling variability of the layer-1 fit CONDITIONAL on frozen covariate construction; understates whole-pipeline variability by the prior-estimation share for both families",
    "posterior draws from one frozen out-of-fold fit (contrast L3P): estimation uncertainty conditional on the observed sample, a different estimand from L3 and never pooled with it",
    "19 signal groups of TreeSHAP; the BIAS and 12 intervention groups are dropped, so rowSums() is NOT the booster's margin (finding A4); the discrimination axis reads the stored margin instead"),
  stringsAsFactors = FALSE)

#' Which fields of two design keys differ, as a table a human can act on.
#'
#' An assertion that says only "these disagree" sends the reader back to
#' CLAUDE.md's working-style rule with no way to obey it: a failing check does
#' not say WHICH BRANCH moved. This names the field, which is usually enough to
#' identify the branch outright -- a `bam` row means the config changed, a
#' `formulas` row means a whitelist or pairing change, a `folds` row means the
#' split or seed moved.
attr_design_diff <- function(k1, k2) {
  nm <- union(names(k1), names(k2))
  one <- function(v) paste(utils::capture.output(utils::str(v)), collapse = " | ")
  rows <- lapply(nm, function(f) {
    a <- k1[[f]]; b <- k2[[f]]
    if (isTRUE(all.equal(a, b))) return(NULL)
    data.frame(field = f,
               n_elements_a = length(a), n_elements_b = length(b),
               hash_a = attr_key_hash(a), hash_b = attr_key_hash(b),
               detail = substr(one(a), 1L, 200L), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(field = character(0), n_elements_a = integer(0),
                      n_elements_b = integer(0), hash_a = character(0),
                      hash_b = character(0), detail = character(0),
                      stringsAsFactors = FALSE)
  }
  out
}

# --- the replicate plan ------------------------------------------------------

#' Every replicate the evaluation will need, enumerated before a row is read.
#'
#' The counterpart of `layer1_jobs()`, and for the same reason: a fit budget
#' that can only be discovered by starting the run is a budget nobody checks.
#' `fits` is the number of `bam` or `xgboost` calls a replicate costs, so
#' `sum(plan$fits, na.rm = TRUE)` is the whole compute bill and is inspectable
#' in one line.
#'
#' EVERY ROW CARRIES THE REPLICATE COORDINATE `(boot_id, seed_id, draw_id)`, so
#' the plan table is also the statement of which contrasts the run will be able
#' to measure. A level that no replicate supports is visible here rather than as
#' an empty table six hours later.
#'
#' THE SIX GENERATION ROUTES, and why a replicate belongs to exactly one:
#'
#'   `ladder`            boot 0, seed 0, draw 0. The arm as fitted out of fold.
#'   `seed`              boot 0, seed s. SHAP only, because `attr_has_seed_noise()`
#'                       is FALSE for every LLR arm and a replicate provably
#'                       identical to another is not an observation.
#'   `bootstrap`         boot b, seed 0. A refit on SHARED BAG b. Both families.
#'   `bootstrap_seeded`  boot b, seed s. SHAP only. THE LEVEL-2-INTO-LEVEL-3
#'                       PROPAGATION: bag b carries more than one algorithmic
#'                       draw, so seed noise and sampling noise can be varied
#'                       together (contrast L3T) as well as separately.
#'   `posterior`         boot 0, draw d. LLR only. Estimation uncertainty
#'                       conditional on the observed sample -- contrast L3P,
#'                       which is a DIFFERENT ESTIMAND from L3 and is never
#'                       pooled with it.
#'   `none`              a marker for a level a method provably cannot have.
#'
#' @param eval_cfg the parsed `config/attribution_eval.yml`.
#' @param n_folds folds per SHAP replicate; the cost of one `shap_oof()`.
#' @param n_pass_fits the fit cost of ONE full pass over every base arm and
#'   every fold. Charged to synthetic rows rather than spread over the
#'   replicates, because one pass serves every LLR method at once and dividing
#'   it seven ways would make no single row's cost mean anything. Two synthetic
#'   rows use it: the posterior basis (one pass) and the LLR bootstrap
#'   (`llr_bootstrap_b` passes). `sum(plan$fits)` therefore stays the real bill,
#'   which is what `layer1_budget()` exists to protect.
attr_replicate_plan <- function(eval_cfg, n_folds, n_pass_fits) {
  methods <- as.character(cfg_req(eval_cfg, "methods"))
  attr_check_methods(methods)
  b_seed    <- as.integer(cfg_req(eval_cfg, "levels", "seed", "b"))
  b_per_bag <- as.integer(cfg_req(eval_cfg, "levels", "seed", "b_per_bag"))
  b_boot    <- as.integer(cfg_req(eval_cfg, "bootstrap", "b"))
  b_draw    <- as.integer(cfg_req(eval_cfg, "levels", "sample", "b"))
  b_lboot   <- as.integer(cfg_req(eval_cfg, "levels", "sample", "llr_bootstrap_b"))
  route     <- as.character(cfg_req(eval_cfg, "levels", "sample", "llr_route"))
  if (!route %in% c("posterior", "bootstrap", "none")) {
    abort_values("attribution_eval.levels.sample.llr_route must be posterior/bootstrap/none",
                 route)
  }
  # THE INDEX OFFSETS ARE PART OF THE KEY AND MUST NOT COLLIDE. `index` exists
  # only to keep two replicates of one method at one level distinguishable;
  # the coordinate columns are what anything downstream reads. Refused loudly
  # rather than left to produce a silent overwrite.
  if (b_boot > 400L || b_lboot > 400L) {
    abort_values("attribution_eval: bootstrap.b and levels.sample.llr_bootstrap_b must be <= 400 (index offsets)",
                 c(b_boot, b_lboot))
  }
  if (b_lboot > b_boot) {
    abort_values(paste0("attribution_eval: levels.sample.llr_bootstrap_b exceeds ",
                        "bootstrap.b, so LLR would be asked for bags that were ",
                        "never declared. The shared manifest is bootstrap.b bags ",
                        "long and every method draws from it"),
                 c(llr_bootstrap_b = b_lboot, bootstrap_b = b_boot))
  }

  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
  R <- function(method, route, index, boot_id, seed_id, draw_id, fits, note) {
    add(method = method, route = route, index = as.integer(index),
        boot_id = as.integer(boot_id), seed_id = as.integer(seed_id),
        draw_id = as.integer(draw_id), fits = as.integer(fits), note = note)
  }

  for (m in methods) {
    fam <- unname(attr_method_family(m))

    # The arm itself. Coordinate (0, 0, 0), always index 1.
    R(m, "ladder", 1L, 0L, 0L, 0L, if (fam == "shap") n_folds else NA_integer_,
      "the arm as fitted, out of fold")

    # Seed noise at the original sample. SHAP only, and the reason is measured.
    if (attr_has_seed_noise(m)) {
      for (s in seq_len(b_seed)) {
        R(m, "seed", s, 0L, s, 0L, n_folds, "re-seeded booster, identical design")
      }
    } else if (b_seed > 0L) {
      R(m, "none", NA_integer_, 0L, NA_integer_, 0L, 0L,
        "bam/fREML is deterministic: level 2 is exactly zero, not unmeasured")
    }

    # Sampling: the shared bags.
    if (fam == "shap") {
      for (b in seq_len(b_boot)) {
        R(m, "bootstrap", b, b, 0L, 0L, n_folds, "refit on shared bag b")
        # The propagation. Each bag additionally carries `b_per_bag` re-seeded
        # fits, which is what turns L2 and L3 into a factorial rather than two
        # unrelated one-way experiments.
        for (s in seq_len(b_per_bag)) {
          R(m, "bootstrap_seeded", 1000L * s + b, b, s, 0L, n_folds,
            "shared bag b, re-seeded: supports contrast L3T")
        }
      }
    } else {
      if (b_draw > 0L && identical(route, "posterior")) {
        for (d in seq_len(b_draw)) {
          R(m, "posterior", d, 0L, 0L, d, 0L,
            "posterior draw from the frozen out-of-fold beta and Vc")
        }
      }
      if (b_lboot > 0L) {
        for (b in seq_len(b_lboot)) {
          R(m, "bootstrap", 500L + b, b, 0L, 0L, NA_integer_,
            "refit on shared bag b; cost is charged to the shared pass row")
        }
      }
      if (b_draw == 0L && b_lboot == 0L) {
        R(m, "none", NA_integer_, NA_integer_, 0L, 0L, 0L,
          "sampling level disabled for this family in config")
      }
    }
  }

  # The two shared-cost rows. Neither is a replicate, both are real compute.
  llr_present <- any(attr_method_family(methods) == "llr")
  if (llr_present && b_draw > 0L && identical(route, "posterior")) {
    R("(posterior basis)", "posterior_basis", NA_integer_, 0L, 0L, 0L,
      n_pass_fits,
      "one pass over every out-of-fold spec-fold, to obtain beta and Vc")
  }
  if (llr_present && b_lboot > 0L) {
    R("(llr bootstrap passes)", "bootstrap_basis", NA_integer_, NA_integer_,
      0L, 0L, b_lboot * n_pass_fits,
      "one full refit pass per shared bag, serving every LLR arm at once")
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' The content fingerprint a cached replicate is stored under.
#'
#' Same discipline as `tests/coupling_attribution.R`'s `.fingerprint()`: a
#' replicate may be reused only when the design that would produce it is the one
#' that did. `design_key` carries the config, the formulas and the folds;
#' `method`, `level` and `index` carry the replicate's identity; `extra` carries
#' anything route-specific -- the actual RNG seed for a booster, the bootstrap
#' seed for a resample -- so two replicates that differ only in their draw get
#' different keys.
attr_replicate_key <- function(method, level, index, design_key, extra = NULL) {
  if (!level %in% ATTR_LEVELS) abort_values("unknown replicate level", level)
  attr_check_methods(method)
  h <- attr_key_hash(list(method = method, level = level,
                          index = as.integer(index), design = design_key,
                          extra = extra))
  sprintf("%s__%s__%02d__%s", method, level, as.integer(index), h)
}

# `attr_pair_level()` AND `ATTR_LEVEL_LABELS` WERE REMOVED ON 2026-09-06, and
# the removal is the point rather than tidying. They decided a pair's hierarchy
# level from `(method, level, index)`, where `level` was the generator's ROUTE
# and `index` was one integer standing in for a bag id, a seed id or a posterior
# draw id depending on which loop produced the row. That triple cannot express a
# pair which differs in two coordinates on purpose -- a re-seeded fit on a
# resampled bag -- so the design below could not be labelled by it at all.
#
# They are not kept alongside `attr_pair_contrast()`. This file's own header
# says the pair's place in the hierarchy is decided "by `attr_pair_level()` and
# nowhere else", and the whole force of that sentence is the "nowhere else":
# two functions that both decide which level a comparison sits at is precisely
# the arrangement in which one of them silently acquires a different convention.
# `tests/attr_metrics.R` was the only caller.

# --- THE REPLICATE COORDINATE, AND THE CONTRASTS DERIVED FROM IT -------------
#
# ADDED 2026-09-06, AND IT REPLACES `level` AS THE FUNDAMENTAL DIMENSION.
#
# The generator's `level` field ("spec", "seed", "sample") is a ROUTE: it says
# which loop produced a replicate. It was also being used as the reproducibility
# hierarchy, and those two are not the same thing. A SHAP replicate fitted on
# bootstrap bag 7 with a re-seeded booster differs from the base fit in TWO
# ways, and `level = "sample"` records only one of them.
#
# So a replicate now carries three integer coordinates beside its method:
#
#   boot_id   0 = the original training sample; b >= 1 = SHARED BAG b
#   seed_id   0 = the base algorithmic seed;    s >= 1 = seed offset s
#   draw_id   0 = the fit itself;               d >= 1 = posterior draw d
#
# and the hierarchy level of a PAIR is derived from WHICH COORDINATES DIFFER,
# by `attr_pair_contrast()` and nowhere else. That is the same guarantee
# `attr_pair_level()` gave for the (method, level, index) triple, extended to a
# design that can vary two things at once on purpose.
#
# WHY `boot_id` IS SHARED ACROSS METHODS AND THIS IS THE WHOLE POINT. Until
# 2026-09-06 the SHAP bootstrap used `bag_of(seed_base + i)` and the LLR
# bootstrap used `bag_of(seed_base + 90000 + i)`, so SHAP replicate 3 and LLR
# replicate 3 were fitted on DIFFERENT resamples and nothing could be paired.
# One bag per `boot_id`, fed to every method, makes the whole experiment a
# BLOCKED design: a level-4 comparison holds the sample fixed, and a
# cross-method level-3 comparison can be paired on the bag pair.

#' Which routes are mutually comparable, and which must stay isolated.
#'
#' ADDED 2026-09-06 TO FIX A DEFECT THAT MADE `L3T` UNCONSTRUCTIBLE. The first
#' version of `attr_contrast_pairs()` grouped the manifest by `(method, route)`
#' and formed pairs only within a group. That is the right instinct applied with
#' the wrong key: it kept the posterior draws away from the bootstrap refits,
#' which is essential, but it ALSO kept `bootstrap` away from
#' `bootstrap_seeded`, and those two routes exist precisely so that their CROSS
#' pairs can be formed.
#'
#' The consequence was silent and in the worse direction. `L3T` -- the whole
#' point of generating a re-seeded fit per bag -- was never constructed, so its
#' table was empty. And `bootstrap_seeded` paired against itself changes only
#' `boot_id` (every one of those replicates is at the same `seed_id`), so it
#' produced a SECOND `L3` DISTRIBUTION that was correctly labelled `L3` and
#' looked like a bonus rather than like the missing propagation.
#'
#' THE RULE IS NOT "never pair different routes". It is:
#'
#'   `refit`      `ladder`, `seed`, `bootstrap`, `bootstrap_seeded`. All four
#'                are genuine refits of the same estimator, differing only in
#'                their coordinate, so any coordinate contrast between them is
#'                meaningful. This is what makes L2-within-a-bag and L3T
#'                expressible at all.
#'   `posterior`  draws from ONE frozen fit. Isolated, because a posterior draw
#'                and a refit estimate different quantities -- measured
#'                2026-09-06 at 1.17 to 1.57 times apart -- and a pair drawn
#'                across the two would be neither.
ATTR_ROUTE_FAMILY <- c(ladder = "refit", seed = "refit", bootstrap = "refit",
                       bootstrap_seeded = "refit", posterior = "posterior")

attr_route_family <- function(route) {
  f <- unname(ATTR_ROUTE_FAMILY[as.character(route)])
  if (anyNA(f)) {
    abort_values("attr_route_family: unknown route(s)",
                 unique(route[is.na(f)]))
  }
  f
}

#' The contrast a pair of replicate coordinates measures.
#'
#' @param a,b lists with `method`, `boot_id`, `seed_id`, `draw_id`.
#' @return a list with `code`, `hierarchy_level` and `held`.
#'
#' THE SIX CONTRASTS.
#'
#'   L1   identity        nothing differs. The nats budget, which sets scale.
#'   L2   seed            algorithmic noise: same spec, same sample, new RNG.
#'   L3   sample          sampling noise: same spec, same RNG, new bag.
#'   L3T  sample_total    OPERATIONAL reproducibility: new bag AND new RNG.
#'                        This is level-2 uncertainty PROPAGATED INTO level 3,
#'                        obtained by design rather than by a variance model:
#'                        when you retrain on another sample you also get
#'                        another random draw, and L3T is that joint quantity.
#'                        L3 and L2 are its two isolated components.
#'   L3P  posterior       estimation noise conditional on the observed sample:
#'                        two posterior draws from one frozen fit. A DIFFERENT
#'                        ESTIMAND from L3 and it must never be pooled with it
#'                        -- measured 2026-09-06 at 1.17 to 1.57 times the
#'                        bootstrap's spread, so pooling would report the wider
#'                        of two quantities as the narrower one's value.
#'   L4   spec            specification: different method, EVERY coordinate
#'                        held. Refused otherwise.
#'
#' A PAIR THAT DIFFERS IN `method` AND IN A COORDINATE IS REFUSED, unchanged
#' from `attr_pair_level()`: it measures specification confounded with
#' estimation noise and there is no level for that.
attr_pair_contrast <- function(a, b) {
  need <- c("method", "boot_id", "seed_id", "draw_id")
  for (f in need) {
    if (is.null(a[[f]]) || is.null(b[[f]])) {
      stop("attr_pair_contrast: a replicate coordinate needs ",
           paste(need, collapse = ", "), call. = FALSE)
    }
  }
  dm <- !identical(as.character(a$method), as.character(b$method))
  db <- !identical(as.integer(a$boot_id),  as.integer(b$boot_id))
  ds <- !identical(as.integer(a$seed_id),  as.integer(b$seed_id))
  dd <- !identical(as.integer(a$draw_id),  as.integer(b$draw_id))

  if (dm) {
    if (db || ds || dd) {
      stop("attr_pair_contrast: methods differ AND a coordinate differs ",
           "(boot ", db, ", seed ", ds, ", draw ", dd, "). That pair measures ",
           "specification confounded with estimation noise, which is not one ",
           "of the six contrasts. Compare at a common coordinate.", call. = FALSE)
    }
    return(list(code = "L4", hierarchy_level = 4L,
                held = sprintf("boot=%d seed=%d draw=%d", as.integer(a$boot_id),
                               as.integer(a$seed_id), as.integer(a$draw_id))))
  }
  if (dd && (db || ds)) {
    stop("attr_pair_contrast: a posterior draw index differs at the same time ",
         "as a bag or a seed. Posterior draws exist only at boot_id 0 and ",
         "seed_id 0, so this pair cannot have been generated and is a ",
         "bookkeeping error rather than a measurement.", call. = FALSE)
  }

  # A ZERO ON THE BAG OR DRAW AXIS MEANS "NOT PERTURBED", NOT "PERTURBATION
  # NUMBER ZERO", and the two refusals below are what stops that being read the
  # wrong way. They were not needed while `attr_contrast_pairs()` grouped by
  # route, because no group then contained both an unperturbed and a perturbed
  # replicate; merging `ladder` into the `refit` family with `bootstrap` makes
  # them load-bearing, and adding one without the other would have traded a
  # missing contrast for a wrong one.
  #
  # BAG. `boot_id = 0` is the FULL training set; `boot_id = b` is a 63.2%
  # patient subsample, because `bag_of()` keeps the distinct draws. Pairing them
  # would confound "a different sample" with "a smaller sample", and the
  # difference in sample SIZE is the larger of the two effects. Level 3 asks
  # what happens between two plausible training samples OF THE SAME KIND.
  #
  # DRAW. `draw_id = 0` is the fit itself and `draw_id = d` is a draw from its
  # posterior. Pairing them is a deviation-from-centre measurement, which is
  # smaller than a pairwise one by a factor that depends on B -- the same
  # argument `attr_displacement()` makes for being pairwise rather than
  # deviation-from-median.
  #
  # SEED HAS NO SUCH ASYMMETRY and is deliberately exempt: `seed_id = 0` is
  # `cfg$seed + fold` and `seed_id = s` is that plus an offset. Both are
  # arbitrary RNG seeds, neither is an absence of perturbation, so `ladder`
  # against a `seed` replicate is a legitimate L2 pair and is now formed.
  if (db && min(as.integer(a$boot_id), as.integer(b$boot_id)) == 0L) {
    stop("attr_pair_contrast: one replicate is at boot_id 0 (the FULL training ",
         "sample) and the other is on a resampled bag (a 63% patient ",
         "subsample). That pair confounds a different sample with a smaller ",
         "one and is not level 3, which compares two plausible training samples ",
         "of the same kind. boot_id 0 pairs only at level 4, where the ",
         "coordinate is held.", call. = FALSE)
  }
  if (dd && min(as.integer(a$draw_id), as.integer(b$draw_id)) == 0L) {
    stop("attr_pair_contrast: one replicate is the fit itself (draw_id 0) and ",
         "the other is a posterior draw from it. That is a ",
         "deviation-from-centre measurement, not a pairwise one, and it is ",
         "smaller by a factor that depends on B.", call. = FALSE)
  }

  if (dd) return(list(code = "L3P", hierarchy_level = 3L, held = "sample"))
  if (db && ds) return(list(code = "L3T", hierarchy_level = 3L, held = "nothing"))
  if (db) return(list(code = "L3", hierarchy_level = 3L,
                      held = sprintf("seed=%d", as.integer(a$seed_id))))
  # THE L2 STRATUM IS COARSE ON PURPOSE. Reporting "boot=7" would give one pair
  # per bag and forty cells of one observation each, which is not a
  # distribution. "original_sample" against "resampled" is the split the design
  # was built to make answerable: is seed noise worse on a resample than on the
  # full training set?
  if (ds) return(list(code = "L2", hierarchy_level = 2L,
                      held = if (as.integer(a$boot_id) == 0L) "original_sample"
                             else "resampled"))
  list(code = "L1", hierarchy_level = 1L, held = "everything")
}

#' Human-readable name for each contrast. Emitted on every table so a number
#' can never be quoted without the thing it varied.
ATTR_CONTRAST_LABELS <- c(
  L1  = "identity (nothing differs)",
  L2  = "seed noise (same spec, same sample, re-seeded fit)",
  L3  = "sampling noise (same spec, same seed, resampled training set)",
  L3T = "total repeated-fit variability (resampled AND re-seeded)",
  L3P = "posterior estimation noise (same sample, two draws from one fit)",
  L4  = "specification (different spec, same sample and same draw)")

#' The full set of within-method comparable pairs in a replicate manifest.
#'
#' ONE FUNCTION SO NO SCRIPT DECIDES FOR ITSELF which pairs are which contrast.
#' Returns row-index pairs into `man` together with the contrast each measures,
#' so the consumer's only job is to load the two matrices and apply a metric.
#'
#' @param man a manifest with `method`, `route`, `boot_id`, `seed_id`, `draw_id`.
#' @return data frame with `i`, `j`, `method`, `route_family`, `route_i`,
#'   `route_j`, `code`, `hierarchy_level`, `held` and the two coordinates.
attr_contrast_pairs <- function(man) {
  attr_require_coords(man, "attr_contrast_pairs")
  empty <- data.frame(i = integer(0), j = integer(0), method = character(0),
                      route_family = character(0), route_i = character(0),
                      route_j = character(0), code = character(0),
                      hierarchy_level = integer(0), held = character(0),
                      boot_i = integer(0), boot_j = integer(0),
                      seed_i = integer(0), seed_j = integer(0),
                      stringsAsFactors = FALSE)
  # GROUPED BY ROUTE **FAMILY**, NOT BY ROUTE, AND THE DIFFERENCE IS A FIXED
  # DEFECT. Grouping by route kept the posterior draws away from the bootstrap
  # refits, which is essential, but it also kept `bootstrap` away from
  # `bootstrap_seeded` -- and those two exist precisely so their CROSS pairs can
  # be formed. `L3T` was therefore never constructed, and `bootstrap_seeded`
  # against itself produced a second `L3` that looked like a bonus rather than
  # like the missing propagation. See `attr_route_family()`.
  grp <- paste(man$method, attr_route_family(man$route), sep = "|")
  ii <- jj <- integer(0)
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) < 2L) next
    cb <- utils::combn(length(idx), 2L)
    ii <- c(ii, idx[cb[1, ]]); jj <- c(jj, idx[cb[2, ]])
  }
  if (!length(ii)) return(empty)

  # VECTORISED, BUT THE DECISION IS STILL MADE IN ONE PLACE. A store with 89
  # SHAP replicates in the refit family has 3,916 pairs, and calling
  # `attr_pair_contrast()` once per pair and `rbind`-ing a data frame each time
  # took longer than every metric in this file put together. But a pair's
  # contrast is a function of a SMALL PATTERN -- which coordinates differ, and
  # whether the bag or draw axis touches its unperturbed zero -- so the function
  # is called once per pattern and its answer broadcast. `attr_pair_contrast()`
  # remains the only thing that decides which contrast a pair sits at, which is
  # the property that matters; what is dropped is per-pair call overhead.
  #
  # THE ZERO FLAGS ARE PART OF THE PATTERN AND WERE NOT, WHICH WAS A LATENT
  # DEFECT IN THE BROADCAST. The first version built its representative pair as
  # (0,0,0) against (db, ds, dd), so a `db` pattern was always probed as
  # boot 0 against boot 1 -- which is now a REFUSED pair. Every genuine L3 pair
  # between two resampled bags would have inherited that refusal and the L3
  # table would have emptied. The representative has to have the same zero
  # structure as the rows it stands for.
  bi <- man$boot_id[ii]; bj <- man$boot_id[jj]
  si <- man$seed_id[ii]; sj <- man$seed_id[jj]
  di <- man$draw_id[ii]; dj <- man$draw_id[jj]
  db <- bi != bj; ds <- si != sj; dd <- di != dj
  bz <- pmin(bi, bj) == 0L
  dz <- pmin(di, dj) == 0L
  pat <- db + 2L * ds + 4L * dd + 8L * bz + 16L * dz
  code <- rep(NA_character_, length(ii)); hl <- rep(NA_integer_, length(ii))
  rep2 <- function(differs, zero) if (differs) {
    if (zero) c(0L, 1L) else c(1L, 2L)
  } else if (zero) c(0L, 0L) else c(1L, 1L)
  for (u in sort(unique(pat))) {
    w <- which(pat == u)[1L]
    B <- rep2(db[w], bz[w]); D <- rep2(dd[w], dz[w])
    S <- if (ds[w]) c(0L, 1L) else c(0L, 0L)
    co <- try(attr_pair_contrast(
      list(method = "x", boot_id = B[1], seed_id = S[1], draw_id = D[1]),
      list(method = "x", boot_id = B[2], seed_id = S[2], draw_id = D[2])),
      silent = TRUE)
    if (inherits(co, "try-error")) next     # refused pattern: dropped, not guessed
    code[pat == u] <- co$code
    hl[pat == u]   <- co$hierarchy_level
  }
  ok <- !is.na(code) & code != "L1"         # L1 is two names for one matrix
  if (!any(ok)) return(empty)
  # `held` is the STRATUM the pair belongs to, recomputed per row because the
  # pattern carries which coordinates differ but not the value being held. It
  # is what a distribution is grouped by, so it has to be coarse enough to
  # leave more than one pair in a cell: see the L2 note in
  # `attr_pair_contrast()`.
  held <- ifelse(code == "L3", sprintf("seed=%d", si),
          ifelse(code == "L2", ifelse(bz, "original_sample", "resampled"),
          ifelse(code == "L3P", "sample", "nothing")))
  data.frame(i = ii[ok], j = jj[ok], method = man$method[ii][ok],
             route_family = attr_route_family(man$route[ii][ok]),
             route_i = man$route[ii][ok], route_j = man$route[jj][ok],
             code = code[ok], hierarchy_level = hl[ok], held = held[ok],
             boot_i = bi[ok], boot_j = bj[ok], seed_i = si[ok], seed_j = sj[ok],
             stringsAsFactors = FALSE)
}

#' Refuse a manifest that predates the replicate coordinates.
#'
#' A STALENESS GUARD ON THE PRODUCER'S OUTPUT, and it is the same argument
#' `tests/attribution_vs_shap.R` makes for the ladder: a guard the producer
#' applies to its own cache protects the producer and nothing else. A generator
#' run from before 2026-09-06 has a manifest whose `index` column mixes a bag
#' id, a seed id and a posterior draw id, and reading it under the new scheme
#' would silently pair a posterior draw with a bootstrap refit.
attr_require_coords <- function(man, who) {
  need <- c("method", "route", "boot_id", "seed_id", "draw_id")
  miss <- setdiff(need, names(man))
  if (length(miss)) {
    stop(who, ": the replicate manifest is missing column(s) ",
         paste(miss, collapse = ", "), ".\nA manifest without the replicate ",
         "coordinates was written by tests/attr_replicates.R from before ",
         "2026-09-06, when `level` was still the fundamental dimension and a ",
         "single `index` column carried a bag id, a seed id and a posterior ",
         "draw id at once. Reading it under the current scheme would pair a ",
         "posterior draw with a bootstrap refit and call the result sampling ",
         "noise.\nRe-run the generator against that directory:\n",
         "  Rscript tests/attr_replicates.R --levels spec,seed,sample --resume <dir>\n",
         "It is content-addressed and resumable, so every replicate already on ",
         "disk is skipped and only the manifest is rewritten.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Refuse a manifest whose rows cannot all be replicates.
#'
#' ADDED 2026-09-09 (review finding A9). `attr_require_coords()` checked that
#' the coordinate COLUMNS exist and nothing about their contents: a duplicated
#' coordinate, a missing one, a route the library has no family for, or two
#' rows claiming different shapes all passed, and `attr_contrast_pairs()` would
#' then pair a replicate with a copy of itself and call it L1 -- dropped, but
#' silently -- or pair two rows for one coordinate and count the contrast
#' twice. Reproduced synthetically in `tests/attr_eval_unit.R`.
attr_validate_manifest <- function(man, who) {
  attr_require_coords(man, who)
  if (!nrow(man)) return(invisible(TRUE))
  attr_check_methods(unique(man$method))
  attr_route_family(unique(man$route))
  for (f in c("boot_id", "seed_id", "draw_id")) {
    v <- man[[f]]
    if (anyNA(v) || any(v != as.integer(v)) || any(v < 0)) {
      stop(who, ": manifest column `", f, "` has a missing, negative or ",
           "non-integer coordinate in ", sum(is.na(v) | v < 0), " row(s).",
           call. = FALSE)
    }
  }
  co <- paste(man$method, man$boot_id, man$seed_id, man$draw_id)
  if (anyDuplicated(co)) {
    abort_values(paste0(who, ": the manifest lists the same replicate ",
                        "coordinate more than once. Every metric would count ",
                        "that coordinate twice. Duplicated (method boot seed ",
                        "draw)"), unique(co[duplicated(co)]))
  }
  if (anyDuplicated(man$key)) {
    abort_values(paste0(who, ": duplicated replicate key(s)"),
                 unique(man$key[duplicated(man$key)]))
  }
  if (all(c("n_rows", "n_cols") %in% names(man))) {
    if (length(unique(man$n_rows)) != 1L || length(unique(man$n_cols)) != 1L) {
      stop(who, ": the manifest describes matrices of more than one shape (",
           paste(unique(paste0(man$n_rows, "x", man$n_cols)), collapse = ", "),
           "). Every metric compares row i with row i.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Refuse a replicate matrix that is not the object every metric assumes.
#'
#' Every function in this file compares ROW i OF A AGAINST ROW i OF B, so a
#' matrix in a different patient order, a different signal order, or with a
#' non-finite cell produces a complete table of confident numbers about the
#' wrong thing. `load_rep()` in the consumer reads a matrix and handed it
#' straight on; a reordered synthetic matrix was compared without error
#' (review finding A9). Row identity is compared as a whole character vector,
#' never printed (hard rule 1).
attr_check_replicate <- function(M, stay_ids, signals, what = "replicate") {
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("attr_check_replicate: ", what, " is not a numeric matrix", call. = FALSE)
  }
  if (nrow(M) != length(stay_ids) || ncol(M) != length(signals)) {
    stop("attr_check_replicate: ", what, " is ", nrow(M), "x", ncol(M),
         ", expected ", length(stay_ids), "x", length(signals), call. = FALSE)
  }
  if (!identical(rownames(M), as.character(stay_ids))) {
    stop("attr_check_replicate: ", what, " is not in the store's stay_id row ",
         "order (or carries no rownames). Row i would be compared against a ",
         "different patient's row i.", call. = FALSE)
  }
  if (!identical(colnames(M), as.character(signals))) {
    stop("attr_check_replicate: ", what, " does not carry the signals in the ",
         "declared order.", call. = FALSE)
  }
  if (!all(is.finite(M))) {
    stop("attr_check_replicate: ", what, " has ", sum(!is.finite(M)),
         " non-finite cell(s). A failed fit must exclude its bag, never leave ",
         "NA or Inf behind.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Planned against present, per (method, route): is this store complete?
#'
#' ADDED 2026-09-09 (review finding A9). `reconcile_manifest()` lists what is
#' on disk and says nothing about what is NOT: a store missing half its
#' bootstrap replicates produced a complete-looking manifest and every
#' downstream table silently used the subset. This names every planned
#' coordinate as present, tombstoned (excluded with a recorded cause) or
#' missing, and both scripts write it.
#'
#' @param plan `attr_replicate_plan()` output.
#' @param man the reconciled manifest.
#' @param tombstoned integer vector of bag ids excluded for the LLR bootstrap.
attr_replicate_coverage <- function(plan, man, tombstoned = integer(0)) {
  routes <- c("ladder", "seed", "bootstrap", "bootstrap_seeded", "posterior")
  pr <- plan[!is.na(plan$index) & plan$route %in% routes, , drop = FALSE]
  pk <- paste(pr$method, pr$boot_id, pr$seed_id, pr$draw_id)
  mk <- paste(man$method, man$boot_id, man$seed_id, man$draw_id)
  pr$present <- pk %in% mk
  pr$tombstoned <- !pr$present & pr$route == "bootstrap" &
    attr_method_family(pr$method) == "llr" & pr$boot_id %in% tombstoned
  pr$missing <- !pr$present & !pr$tombstoned
  ag <- aggregate(cbind(planned = 1L, present = pr$present,
                        tombstoned = pr$tombstoned, missing = pr$missing),
                  by = list(method = pr$method, route = pr$route), FUN = sum)
  ag$complete <- ag$missing == 0L
  # Replicates on disk the plan does not name: a stale config that shrank a
  # count, or a file from another design. Reported, never used.
  extra <- sum(!(mk %in% pk))
  attr(ag, "n_unplanned_on_disk") <- extra
  ag[order(ag$method, ag$route), , drop = FALSE]
}

# --- THE DISTRIBUTION SUMMARY ------------------------------------------------

#' Every summary statistic a disagreement distribution is reported with.
#'
#' WHY THIS EXISTS AS ONE FUNCTION. Until 2026-09-06 every level of the
#' hierarchy was reported as a SINGLE NUMBER computed from the first pair of
#' replicates -- `mats[[1]]` against `mats[[2]]` -- and the remaining 778 pairs
#' were generated and then thrown away. A single number cannot say how wide the
#' disagreement gets, what the worst observed case is, or whether the
#' distribution is skewed, and those are the questions a reproducibility claim
#' has to answer. Fixing the column set in one function means every level is
#' reported the same way, so two levels can be read side by side without the
#' reader first checking which statistics each one happened to carry.
#'
#' THE MEAN AND THE MEDIAN ARE BOTH HERE ON PURPOSE. Their difference is the
#' cheapest honest skew indicator: a disagreement distribution whose mean sits
#' well above its median has a heavy right tail, which means the median
#' understates how bad a bad resample is. `skew_g1` is the standard sample
#' skewness beside it.
#'
#' THE MAXIMUM IS A DIAGNOSTIC AND NOT A CRITERION, and the distinction is
#' declared here rather than left to the reader. A maximum mechanically grows
#' as more pairs are generated and is set by one pathological fit, so a
#' threshold on it would tighten every time B increased. `p95` is the
#' robustness criterion; `max` is the worst observed case.
attr_stat_summary <- function(v) {
  v <- as.numeric(v); x <- v[is.finite(v)]
  if (!length(x)) {
    return(data.frame(n = length(v), n_scored = 0L, min = NA_real_,
                      p05 = NA_real_, p10 = NA_real_, q1 = NA_real_,
                      median = NA_real_, mean = NA_real_, q3 = NA_real_,
                      p90 = NA_real_, p95 = NA_real_, max = NA_real_,
                      iqr = NA_real_, sd = NA_real_,
                      mean_minus_median = NA_real_, skew_g1 = NA_real_,
                      stringsAsFactors = FALSE))
  }
  q <- unname(stats::quantile(x, c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)))
  m <- mean(x)
  # `skew_g1` = m3 / m2^1.5, the moment estimator. NA rather than 0/0 for a
  # degenerate distribution, which is the honest answer when every pair agrees.
  m2 <- mean((x - m)^2); m3 <- mean((x - m)^3)
  data.frame(n = length(v), n_scored = length(x),
             min = round(min(x), 6), p05 = round(q[1], 6), p10 = round(q[2], 6),
             q1 = round(q[3], 6), median = round(q[4], 6), mean = round(m, 6),
             q3 = round(q[5], 6), p90 = round(q[6], 6), p95 = round(q[7], 6),
             max = round(max(x), 6), iqr = round(q[5] - q[3], 6),
             sd = round(if (length(x) > 1L) stats::sd(x) else NA_real_, 6),
             mean_minus_median = round(m - q[4], 6),
             skew_g1 = round(if (m2 > 0) m3 / m2^1.5 else NA_real_, 4),
             stringsAsFactors = FALSE)
}

# --- FAST TOP-K AGREEMENT, SO A DISTRIBUTION IS AFFORDABLE -------------------

#' Bitmask encoding of each patient's top-k set, one integer per patient.
#'
#' WHY A BITMASK. The distribution tables need top-k agreement over hundreds of
#' PAIRS, and the obvious route -- `attr_prep()` then `attr_agree_k(delta = 0)`
#' -- allocates two n x p logical matrices per call and needs one `order()` per
#' patient per replicate. Encoding the top-k set as an integer whose set bits
#' are the chosen columns makes strict agreement a single integer comparison,
#' `mean(mask_a == mask_b)`, and makes the per-replicate cost k calls to
#' `max.col` at C level.
#'
#' THE SEMANTICS ARE `attr_agree_k(delta = 0)`'s, EXACTLY, INCLUDING TIES AT
#' THE k-th BOUNDARY -- and until 2026-09-09 they were not (review finding
#' A3). The first version compared the two top-k SETS for equality, with ties
#' broken to the lowest column index. `attr_agree_k()` is tie-tolerant even at
#' delta = 0: it forgives a member of the symmetric difference whenever the
#' other arm holds it EXACTLY tied with its own k-th value. Hand-checkable:
#' A = (3, 2, 1, 1) picks {1, 2, 3}, B = (1, 1, 3, 2) picks {1, 3, 4}; set
#' equality says disagree, the tie-aware rule says agree because column 4 is
#' tied at A's third-place boundary and column 2 at B's. The gate compared the
#' two functions only on the anchor ladder, where no boundary tie happened to
#' occur, and passed. That is a test of the inputs, not of the definition.
#'
#' The tie-tolerant reading is the documented one ("two arms agree when
#' NEITHER strongly prefers its own choice"), it is what the delta grid
#' reduces to at zero, and it is what `tests/attr_external.R` uses, so it is
#' the definition both sites now share. The fast path keeps its speed by also
#' encoding, per k, the BOUNDARY-TIE SET: every column whose |value| equals the
#' patient's k-th largest. Agreement is then two bitwise tests per patient:
#' every column b has and a lacks must lie in a's tie set, and vice versa. The
#' counterexample above is a regression test in the gate and in
#' `tests/attr_eval_unit.R`.
#'
#' Requires at most 30 columns, which the design guarantees: 19 signals or 11
#' domains. Refused rather than silently overflowing an integer.
#'
#' @return per k, a list with `top` (bitmask of the top-k set) and `tie`
#'   (bitmask of the columns tied with the k-th largest value).
attr_topk <- function(M, ks) {
  p <- ncol(M)
  if (p > 30L) {
    stop("attr_topk: bitmask encoding needs at most 30 columns, got ", p,
         call. = FALSE)
  }
  ks <- sort(unique(as.integer(ks)))
  if (max(ks) > p) abort_values("attr_topk: k exceeds the column count", max(ks))
  ab <- abs(M); work <- ab; n <- nrow(work)
  pow <- bitwShiftL(1L, seq_len(p) - 1L)
  cur <- integer(n); out <- list()
  for (k in seq_len(max(ks))) {
    j <- max.col(work, ties.method = "first")
    cur <- cur + pow[j]
    if (k %in% ks) {
      kth <- ab[cbind(seq_len(n), j)]
      tie <- numeric(n)
      for (jj in seq_len(p)) tie <- tie + pow[jj] * (ab[, jj] == kth)
      out[[as.character(k)]] <- list(top = cur, tie = as.integer(tie))
    }
    if (k < max(ks)) work[cbind(seq_len(n), j)] <- -Inf
  }
  out
}

#' Tie-tolerant top-k agreement at delta = 0 between two encodings.
#'
#' A violation by arm a is a column in b's set, absent from a's, that is NOT
#' tied at a's k-th boundary -- which is precisely `|A[g]| < kth_a`, because a
#' column outside the top-k set can never exceed the k-th value. Symmetric in
#' the other direction. Identical to `attr_agree_k(pa, pb, k, 0)`.
attr_topk_agree <- function(ta, tb, k) {
  kk <- as.character(as.integer(k))
  if (is.null(ta[[kk]]) || is.null(tb[[kk]])) {
    abort_values("attr_topk_agree: no encoding cached for k", k)
  }
  Sa <- ta[[kk]]$top; Sb <- tb[[kk]]$top
  Ta <- ta[[kk]]$tie; Tb <- tb[[kk]]$tie
  viol_a <- bitwAnd(bitwAnd(Sb, bitwNot(Sa)), bitwNot(Ta)) != 0L
  viol_b <- bitwAnd(bitwAnd(Sa, bitwNot(Sb)), bitwNot(Tb)) != 0L
  mean(!(viol_a | viol_b))
}

#' A bounded, deterministic sample of index pairs from `n` replicates.
#'
#' B = 40 gives 780 pairs and every metric that needs the full matrices costs a
#' pass over 783,750 cells, so an exhaustive sweep is neither affordable nor
#' needed for a quantile. Same discipline as `attr_displacement()`: the bound
#' and the seed are declared in config before generation, so the pair set is a
#' property of the plan rather than of the run.
attr_pair_index <- function(n, max_pairs, seed) {
  if (n < 2L) return(matrix(integer(0), nrow = 2L))
  cb <- utils::combn(n, 2L)
  if (ncol(cb) <= max_pairs) return(cb)
  with_seed(seed, cb[, sort(sample.int(ncol(cb), max_pairs)), drop = FALSE])
}

# --- COMPARING TWO DISTRIBUTIONS WITHOUT DIVIDING ---------------------------

#' Probabilistic dominance of one disagreement distribution over another.
#'
#' THE REPLACEMENT FOR THE LEVEL-4-OVER-LEVEL-3 RATIO, and the reasons are
#' arithmetic rather than taste.
#'
#' A ratio has a denominator that can approach zero, so it explodes for exactly
#' the arms whose sampling noise is smallest -- which is to say, the best ones.
#' It also collapses two distributions to one number and then invites a
#' threshold on it to be read as a decision, when the two distributions are not
#' on one conceptual axis at all: a large level 3 means THIS METHOD IS NOT
#' REPRODUCIBLE, a large level 4 means TWO SPECIFICATIONS TELL DIFFERENT
#' STORIES, and the second is not a defect in either of them.
#'
#' `p_dominates` is P(X > Y) + 0.5 P(X = Y) over all cross pairs, the
#' common-language effect size. It is bounded in [0, 1], needs no denominator,
#' and answers what the ratio was reaching for: how often does a specification
#' difference exceed the sampling noise it has to clear? 0.5 is
#' indistinguishable; 1.0 is always larger.
#'
#' `median_gap` is on the disagreement scale, so it reads directly as "this
#' many more patients per hundred change their top signal because of the
#' specification than because of the resample".
#'
#' EXACT, NOT SAMPLED. Both inputs are already per-pair summaries of the order
#' of hundreds of values, so the full outer comparison is affordable and there
#' is no seed to declare.
#' Which way a per-pair metric points, so `attr_dominance()` compares like
#' with like.
#'
#' ADDED 2026-09-09 (review finding A8). Every per-pair scalar was fed to the
#' same "is level 4 larger than level 3" comparison, and for the SHARE metrics
#' -- `share_ab_median`, `share_ab_p05` and their mirrors, the share of the
#' other arm's evidence the leader still carries -- larger means LESS
#' disagreement. Synthetic shares of 0.8-0.9 against 0.1-0.2 reported
#' dominance 1 in the direction the table labels "more disagreement". The
#' consumer now maps each metric through `attr_orient()` before comparing.
#'
#' @return "disagreement" (larger is more disagreement) or "agreement".
attr_metric_orientation <- function(metric) {
  ifelse(grepl("^share_", metric), "agreement", "disagreement")
}

#' A metric's values oriented so that larger is more disagreement.
#'
#' Shares live in [0, 1], so `1 - share` is the share LOST, on the same scale.
#' Returns the oriented values and the name they should be reported under.
attr_orient <- function(values, metric) {
  if (identical(attr_metric_orientation(metric), "agreement")) {
    list(value = 1 - values, metric = paste0("one_minus_", metric),
         orientation = "larger_is_more_disagreement (transformed)")
  } else {
    list(value = values, metric = metric,
         orientation = "larger_is_more_disagreement")
  }
}

attr_dominance <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) {
    return(data.frame(n_x = length(x), n_y = length(y), p_dominates = NA_real_,
                      median_x = NA_real_, median_y = NA_real_,
                      median_gap = NA_real_, p95_x = NA_real_, p95_y = NA_real_,
                      p95_gap = NA_real_, stringsAsFactors = FALSE))
  }
  mx <- stats::median(x); my <- stats::median(y)
  qx <- unname(stats::quantile(x, 0.95)); qy <- unname(stats::quantile(y, 0.95))
  data.frame(n_x = length(x), n_y = length(y),
             p_dominates = round(mean(outer(x, y, ">")) +
                                 0.5 * mean(outer(x, y, "==")), 4),
             median_x = round(mx, 6), median_y = round(my, 6),
             median_gap = round(mx - my, 6),
             p95_x = round(qx, 6), p95_y = round(qy, 6),
             p95_gap = round(qx - qy, 6), stringsAsFactors = FALSE)
}

#' The difference between two methods measured on THE SAME bag pairs.
#'
#' WHAT THE SHARED BOOTSTRAP MANIFEST BUYS, and the only thing in this file
#' that would be impossible without it. Two methods' level-3 distributions can
#' always be compared marginally, but a marginal comparison still carries the
#' variation caused by WHICH PATIENTS HAPPENED TO ENTER each bag. When both
#' methods are fitted on the identical bags the comparison is made bag-pair by
#' bag-pair and that variation differences out. It is the same paired-versus-
#' marginal argument `tests/metrics_severity.R` already makes for the severity
#' arms: every cell scored on identical rows, so marginal intervals settle
#' nothing.
#'
#' @param a,b data frames with `boot_i`, `boot_j` and `value`, one row per pair.
#' @return the distribution of `value_a - value_b` over the bag pairs they share.
attr_paired_delta <- function(a, b) {
  key <- function(z) paste(pmin(z$boot_i, z$boot_j), pmax(z$boot_i, z$boot_j),
                           sep = "-")
  ka <- key(a); kb <- key(b)
  # REFUSED, NOT RESOLVED BY TAKING THE FIRST. `match()` below returns the first
  # hit, so a duplicated bag pair would silently pick one of two values and the
  # caller would never learn which. That is reachable: once each bag carries
  # more than one seed, a method has an L3 distribution at `seed=0` AND at
  # `seed=1`, and the same (b, b') appears in both. Those are two valid
  # measurements of one contrast on one bag pair, and averaging or taking either
  # is a decision the CALLER has to make explicitly by filtering to a stratum --
  # not one this function may take on its behalf.
  for (nm in c("a", "b")) {
    kk <- if (nm == "a") ka else kb
    if (anyDuplicated(kk)) {
      d <- unique(kk[duplicated(kk)])
      abort_values(paste0("attr_paired_delta: argument `", nm, "` has more than ",
                          "one row for the same bag pair, so the pairing would ",
                          "be ambiguous. Filter to one stratum (one route and ",
                          "one held seed) before calling. Duplicated bag ",
                          "pair(s)"), utils::head(d, 5))
    }
  }
  sh <- intersect(ka, kb)
  if (!length(sh)) {
    return(cbind(data.frame(n_shared_bag_pairs = 0L),
                 attr_stat_summary(numeric(0))))
  }
  d <- a$value[match(sh, ka)] - b$value[match(sh, kb)]
  # `frac_positive` is the sign test's statistic and is the reading a median
  # cannot give: a median gap of +0.02 could come from every bag pair moving
  # +0.02 or from three quarters moving +0.05 and a quarter moving -0.05, and
  # those support very different sentences.
  cbind(data.frame(n_shared_bag_pairs = length(sh),
                   frac_positive = round(mean(d > 0, na.rm = TRUE), 5),
                   stringsAsFactors = FALSE),
        attr_stat_summary(d))
}

# --- aggregation to the reporting partition ----------------------------------

#' Roll a signal-level attribution matrix up to the frozen 11 domains.
#'
#' `config/domains.csv` is an aggregation guide for layer 2 only and changes no
#' model. Equal weight within a domain, which is what `D_k` reduces to and the
#' same convention `R/09d_sofa.R` uses.
attr_to_domain <- function(M, domains) {
  sg  <- colnames(M)
  dm  <- domains$domain[match(sg, domains$signal)]
  if (anyNA(dm)) abort_values("attr_to_domain: signal(s) absent from domains.csv",
                              sg[is.na(dm)])
  dnm <- sort(unique(dm))
  D <- matrix(0, nrow(M), length(dnm), dimnames = list(rownames(M), dnm))
  for (k in dnm) {
    j <- which(dm == k)
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE])
  }
  D
}

# --- the metric family -------------------------------------------------------
#
# Every function below takes matrices and returns a data frame or a numeric
# vector. None of them knows what a run directory is.

#' Precomputed orderings for one attribution matrix.
#'
#' `ord` is the column order by descending |A|, so the k-th largest value and
#' top-k membership both fall out by indexing rather than by a second pass over
#' the rows. One `apply` per matrix instead of one per comparison per k per
#' delta, which is the difference between seconds and an hour.
attr_prep <- function(M) {
  ab <- abs(M); n <- nrow(ab)
  ord <- t(apply(-ab, 1, order))
  # MEMOISED ON 2026-09-06, and it is a speedup with no numeric consequence.
  # The distribution tables call `attr_agree_k()` once per (pair, k, delta), so
  # a prep that recomputes `inS(k)` -- an n x p logical, 783,750 cells -- on
  # every call rebuilt the same matrix six times per pair per k. The cache is
  # keyed on k, lives in this prep's own environment, and dies with it. The
  # gate re-runs every legacy table through this function and asserts the
  # numbers are unchanged.
  .cS <- new.env(parent = emptyenv()); .cK <- new.env(parent = emptyenv())
  list(ab = ab, ord = ord,
       kth = function(k) {
         kk <- as.character(k)
         if (is.null(.cK[[kk]])) .cK[[kk]] <- ab[cbind(seq_len(n), ord[, k])]
         .cK[[kk]]
       },
       inS = function(k) {
         kk <- as.character(k)
         if (is.null(.cS[[kk]])) {
           S <- matrix(FALSE, n, ncol(ab))
           S[cbind(rep(seq_len(n), times = k), as.vector(ord[, seq_len(k)]))] <- TRUE
           .cS[[kk]] <- S
         }
         .cS[[kk]]
       },
       total = rowSums(ab), n = n, p = ncol(ab))
}

#' Tie-tolerant agreement at rank k.
#'
#' THE DEFINITION, AND WHY THIS ONE. Two arms agree at rank k with tolerance
#' `delta` when NEITHER strongly prefers its own choice to the other's. A real
#' disagreement exists if some g in S_b but not S_a has |A_a[g]| < kth_a - delta
#' -- arm a says b's pick is materially worse -- or symmetrically. At delta = 0
#' this is exact set equality, so the grid contains the strict number.
#'
#' SYMMETRIC ON PURPOSE. A one-sided version would call it agreement whenever
#' arm a happens to be flat, even if arm b is emphatic that a's pick is wrong.
#' Requiring both directions means a disagreement is forgiven only when both
#' arms are genuinely undecided.
#'
#' @param delta scalar, or a per-row vector for a relative tolerance.
#' @return logical, one per patient.
attr_agree_k <- function(pa, pb, k, delta) {
  Sa <- pa$inS(k); Sb <- pb$inS(k)
  ka <- pa$kth(k); kb <- pb$kth(k)
  viol_a <- rowSums(Sb & !Sa & (pa$ab < ka - delta)) > 0
  viol_b <- rowSums(Sa & !Sb & (pb$ab < kb - delta)) > 0
  !(viol_a | viol_b)
}

#' Row minimum of `ab` over the cells where `mask` is TRUE; Inf where none.
.masked_row_min <- function(ab, mask) {
  out <- rep(Inf, nrow(ab))
  for (j in seq_len(ncol(ab))) {
    v <- ab[, j]; v[!mask[, j]] <- Inf
    out <- pmin(out, v)
  }
  out
}

#' Agreement rate at rank k over a WHOLE GRID of tolerances, in one pass.
#'
#' ADDED 2026-09-09 (review finding A7). The magnitude tables evaluated one
#' noise-calibrated delta only, although config declares three tolerance
#' families -- noise-calibrated at each declared quantile, the common absolute
#' grid, and the relative grid -- and the file header promised all three. The
#' cost of `attr_agree_k()` is two n x p logical passes PER DELTA, which is why
#' the grid was never evaluated on the distribution tables. But the delta
#' enters only through `min_{g in Sb \ Sa} |A[g]| < kth_a - delta`, so the
#' masked row minimum can be taken ONCE per (pair, k) and every delta is then
#' an O(n) comparison. Same answer as `attr_agree_k()` at every scalar delta;
#' asserted in `tests/attr_eval_unit.R`.
#'
#' `relative = TRUE` scales the tolerance by each arm's OWN row total --
#' `delta * total_a` on a's side and `delta * total_b` on b's -- so the
#' tolerance is symmetric between the arms. The legacy relative grid of
#' `tests/attribution_ties.R`, reproduced by the gate through
#' `attr_agree_k(pa, pb, k, d * pa$total)`, applied ARM a's total to both
#' sides; that is order-dependent and is kept only where a legacy table has to
#' be reproduced. Tables built here carry `delta_kind = relative_own_total`.
#'
#' @return numeric vector of agreement rates, one per delta.
attr_agree_k_grid <- function(pa, pb, k, deltas, relative = FALSE) {
  Sa <- pa$inS(k); Sb <- pb$inS(k)
  ka <- pa$kth(k); kb <- pb$kth(k)
  mina <- .masked_row_min(pa$ab, Sb & !Sa)
  minb <- .masked_row_min(pb$ab, Sa & !Sb)
  sa <- if (relative) pa$total else 1
  sb <- if (relative) pb$total else 1
  vapply(deltas, function(d) {
    mean(!((mina < ka - d * sa) | (minb < kb - d * sb)))
  }, numeric(1))
}

#' Magnitude-gated sign flip.
#'
#' A cell flips when the two arms disagree about the DIRECTION of that patient's
#' evidence from that signal. GATED: the flip counts only when max(|a|,|b|)
#' exceeds tau nats, so a reversal of a quantity that was never distinguishable
#' from zero is not scored as a disagreement. tau = 0 is the ungated number and
#' is an UPPER BOUND that should not be quoted alone.
#'
#' @param keep optional logical matrix of cells the comparison may use. For the
#'   LLR arms this is the measured mask: an unmeasured signal gets L = 0 by
#'   ASSIGNMENT, so those cells agree perfectly without any estimation having
#'   happened, and counting them inflates every agreement score. SHAP has no
#'   analogue -- a tree assigns a nonzero contribution to a missing feature
#'   through its default direction -- so its `keep` is legitimately all cells,
#'   and the asymmetry is surfaced by the `n_cells` column rather than hidden.
attr_sign_flip <- function(A, B, keep = NULL, taus) {
  stopifnot(identical(dim(A), dim(B)))
  if (is.null(keep)) keep <- matrix(TRUE, nrow(A), ncol(A))
  fl <- sign(A) != sign(B)
  mx <- pmax(abs(A), abs(B))
  data.frame(tau = taus, n_cells = sum(keep),
             flip = round(vapply(taus, function(t) mean(fl[keep] & mx[keep] > t),
                                 numeric(1)), 5),
             stringsAsFactors = FALSE)
}

#' Resolution: how many contributions are in contention at all.
#'
#' NOT a comparison. A property of ONE arm, and the ceiling on what any
#' per-patient explanation can claim: if the typical patient has four signals
#' within a tenth of a nat of the leader, then "the top signal for this patient"
#' is not a supportable statement however stable it happens to be across
#' specifications. This is what catches DEGENERATE STABILITY -- a model that
#' assigns everything to one signal is perfectly stable and useless.
attr_tieset <- function(pa, deltas) {
  lead <- pa$kth(1L)
  do.call(rbind, lapply(deltas, function(d) {
    sz <- rowSums(pa$ab >= lead - d)
    data.frame(delta = d, n_columns = pa$p,
               tieset_median = stats::median(sz), tieset_mean = round(mean(sz), 3),
               tieset_p90 = unname(stats::quantile(sz, 0.90)),
               frac_unique_leader = round(mean(sz == 1L), 5),
               stringsAsFactors = FALSE)
  }))
}

#' Signed cosine similarity of two attribution vectors, per patient.
#'
#' THE ONE MAGNITUDE-SENSITIVE METRIC IN THE FAMILY, and it fills a real hole:
#' every other metric here is rank-based or sign-based and therefore discards
#' magnitude entirely. Cosine is scale-invariant but magnitude-sensitive -- it
#' notices when two arms agree on the ordering and disagree about HOW MUCH.
#'
#' SIGNED, not absolute. Using |A| would make a sign flip in a large component
#' invisible, which is the failure the gated sign-flip metric exists to catch;
#' on the signed vector a flip in a dominant component drives the similarity
#' toward -1 and is penalised properly.
#'
#' Returns one value per patient, NA where either vector is the zero vector (a
#' patient with no measured signal has no direction to compare). ROW-LEVEL:
#' summarise with `attr_dist_summary()`, never print.
attr_cosine <- function(A, B, na = NULL, nb = NULL,
                        min_norm = .Machine$double.eps) {
  stopifnot(identical(dim(A), dim(B)))
  # THE NORMS ARE A PROPERTY OF ONE MATRIX, NOT OF THE PAIR, so they are
  # accepted precomputed. Over 703 pairs of one method's replicates the naive
  # version recomputes each replicate's norm 37 times and spends two thirds of
  # its work there; passing `attr_row_norm()` once per replicate turns three
  # passes over 783,750 cells per pair into one. Identical arithmetic either
  # way -- the default path still computes them.
  if (is.null(na)) na <- attr_row_norm(A)
  if (is.null(nb)) nb <- attr_row_norm(B)
  ok <- na > min_norm & nb > min_norm
  out <- rep(NA_real_, nrow(A))
  out[ok] <- rowSums(A[ok, , drop = FALSE] * B[ok, , drop = FALSE]) / (na[ok] * nb[ok])
  out
}

#' Per-row Euclidean norm. Cache one per replicate; see `attr_cosine()`.
attr_row_norm <- function(A) sqrt(rowSums(A * A))

#' Cosine DISSIMILARITY, so that larger means more disagreement.
#'
#' WHY THE FLIP IS NOT COSMETIC. `attr_dominance()` answers "how often does a
#' level-4 pair disagree MORE than a level-3 pair", and that question only has
#' one meaning if every metric fed to it is oriented the same way. Top-k
#' disagreement is `1 - agreement` and rises with disagreement; cosine
#' SIMILARITY falls with it. Putting the two in one table without flipping one
#' would give the same column name opposite meanings on adjacent rows, which is
#' the kind of thing a reader has no way to detect.
#'
#' Bounded in [0, 2] on signed vectors: 0 identical direction, 1 orthogonal,
#' 2 exactly opposed.
attr_cosine_dissim <- function(A, B, na = NULL, nb = NULL) 1 - attr_cosine(A, B, na, nb)

#' Per-row Pearson correlation of two matrices, vectorised.
#'
#' On rank matrices this IS the per-patient Spearman correlation. Written out
#' rather than looped: `cor()` once per patient is 41,250 calls per comparison,
#' which dominated the ladder script; this is three matrix passes. NA where a
#' row has no variance, which happens when a patient has one measured signal.
attr_row_cor <- function(X, Y) {
  Xc <- X - rowMeans(X); Yc <- Y - rowMeans(Y)
  den <- sqrt(rowSums(Xc^2) * rowSums(Yc^2))
  ifelse(den > 0, rowSums(Xc * Yc) / den, NA_real_)
}

#' Within-patient ranking agreement: top-1, top-3 and the rank correlation.
#'
#' `rho` is CONSERVATIVE for the LLR arms: unmeasured signals sit at exactly 0
#' in both arms, so they are tied and agreeing, which inflates it. `top1` and
#' `top3` do not have that problem, because a zero cell can never be a top
#' contributor for a patient with any measured signal. Read those two.
attr_rank_agreement <- function(A, B) {
  ra <- t(apply(abs(A), 1, rank, ties.method = "average"))
  rb <- t(apply(abs(B), 1, rank, ties.method = "average"))
  rho <- attr_row_cor(ra, rb)
  t3 <- function(M) t(apply(abs(M), 1, function(v) sort(order(-v)[seq_len(3L)])))
  data.frame(top1_agree = round(mean(max.col(abs(A), ties.method = "first") ==
                                     max.col(abs(B), ties.method = "first")), 5),
             top3_agree = round(mean(rowSums(t3(A) == t3(B)) == 3L), 5),
             rho_median = round(stats::median(rho, na.rm = TRUE), 5),
             rho_p10 = round(unname(stats::quantile(rho, 0.10, na.rm = TRUE)), 5),
             n_patients = nrow(A), stringsAsFactors = FALSE)
}

#' The nats budget and its concentration. HIERARCHY LEVEL 1.
#'
#' Level 1 is not a comparison and is not a defect being measured: it is the
#' SCALE against which every tolerance at levels 2 to 4 has to be read. An
#' absolute tolerance of 0.25 nats forgives twice as much of a SHAP
#' disagreement as of an L disagreement, because the median patient carries
#' 3.24 nats of |SHAP| against 6.84 of |L_cond|. Reporting an absolute-delta
#' curve without this table beside it is the scale confound.
#'
#' `hhi` is the Herfindahl index of the |A| shares: 1/p when every column
#' contributes equally, 1 when one column carries everything.
attr_budget <- function(M) {
  ab <- abs(M); tot <- rowSums(ab)
  sh <- ab / pmax(tot, .Machine$double.eps)
  data.frame(n_columns = ncol(M), n_patients = nrow(M),
             budget_median = round(stats::median(tot), 5),
             budget_p90 = round(unname(stats::quantile(tot, 0.90)), 5),
             hhi_median = round(stats::median(rowSums(sh^2)), 5),
             max_share_median = round(stats::median(apply(sh, 1, max)), 5),
             stringsAsFactors = FALSE)
}

#' Quantile summary of a ROW-LEVEL distribution -- one value per patient.
#'
#' A THIN WRAPPER OVER `attr_stat_summary()` SINCE 2026-09-06, and the two are
#' kept apart because they summarise different populations, not because they
#' compute different things. `attr_stat_summary()` summarises a distribution
#' over PAIRS OF REPLICATES; this one summarises a distribution over PATIENTS,
#' where NA is meaningful -- `attr_cosine()` returns NA for a patient whose
#' attribution vector is exactly zero, and `frac_scored` is how many patients
#' the metric could be computed for at all. `frac_above` is the tail share
#' against a declared threshold, which a pair-level table has no use for.
#'
#' The quantile arithmetic was duplicated here until 2026-09-06 and the two
#' copies had already drifted -- this one rounded to five decimals and reported
#' neither p05 nor the mean, so a cosine table could not be read beside a
#' disagreement table without the reader noticing which statistics each carried.
attr_dist_summary <- function(v, above = 0.9) {
  ok <- !is.na(v) & is.finite(v)
  cbind(attr_stat_summary(v),
        data.frame(frac_scored = round(mean(ok), 5), threshold = above,
                   frac_above = round(if (any(ok)) mean(v[ok] > above) else NA_real_, 5),
                   stringsAsFactors = FALSE))
}

# --- the delta policy --------------------------------------------------------

#' A sample of per-cell absolute displacements between replicates of one method.
#'
#' The input to the noise-calibrated tolerance. `delta` stops being a declared
#' grid and becomes a measured quantity: two contributions are tied if they
#' differ by less than this method's own noise moves them.
#'
#' PAIRWISE, NOT DEVIATION-FROM-CENTRE. The tolerance is used to decide whether
#' TWO replicates disagree, so the distribution it is calibrated on must be the
#' distribution of |a - b| between two replicates. A deviation-from-median
#' summary answers a different question and is smaller by a factor that depends
#' on B, which would silently tighten the tolerance as B grows.
#'
#' BOUNDED AND SEEDED, because B = 50 is 1,225 pairs of 783,750-cell matrices
#' and the exact distribution is neither affordable nor needed for a quantile.
#' `max_pairs` pairs are drawn without replacement from the full set and
#' `cell_frac` of the kept cells is sampled from each, under `seed`. Both are
#' declared in config before generation.
#'
#' @param mats list of matrices, all identically shaped: the replicates of ONE
#'   (method, level).
#' @param keep optional logical cell mask, as in `attr_sign_flip()`.
#' @return numeric vector of sampled |a - b|. ROW-LEVEL; summarise, never print.
attr_displacement <- function(mats, keep = NULL, max_pairs, cell_frac, seed) {
  B <- length(mats)
  if (B < 2L) return(numeric(0))
  d <- dim(mats[[1]])
  for (m in mats) if (!identical(dim(m), d)) {
    stop("attr_displacement: replicates are not identically shaped", call. = FALSE)
  }
  if (is.null(keep)) keep <- matrix(TRUE, d[1], d[2])
  idx <- which(keep)
  if (!length(idx)) return(numeric(0))

  pairs <- utils::combn(B, 2L)
  with_seed(seed, {
    if (ncol(pairs) > max_pairs) pairs <- pairs[, sample.int(ncol(pairs), max_pairs), drop = FALSE]
    n_take <- max(1L, as.integer(round(cell_frac * length(idx))))
    unlist(lapply(seq_len(ncol(pairs)), function(j) {
      i <- if (n_take < length(idx)) sample(idx, n_take) else idx
      abs(mats[[pairs[1, j]]][i] - mats[[pairs[2, j]]][i])
    }), use.names = FALSE)
  })
}

#' Turn a measured displacement distribution into a tolerance.
#'
#' THE HAZARD, STATED SO IT CANNOT BE FORGOTTEN: a noisier method earns a bigger
#' cushion under a noise-calibrated delta, which flatters it. That is why the
#' COMMON absolute delta is not optional and why every agreement table is
#' reported three ways -- at delta = 0 (strict, and comparable with every number
#' reported before this policy existed), at `delta_noise(method)` (is this
#' disagreement bigger than the method's own noise?), and at a common absolute
#' delta (do the two agree in absolute evidence terms?). The three answer
#' different questions and all three belong in the paper.
#'
#' @param disp numeric displacements from `attr_displacement()`.
#' @param quantiles the declared quantiles, from config. The median and the 90th
#'   percentile are both defensible and both should be reported.
attr_delta_policy <- function(disp, quantiles) {
  if (!length(disp)) {
    return(data.frame(q = quantiles, delta = NA_real_, n_displacements = 0L,
                      stringsAsFactors = FALSE))
  }
  data.frame(q = quantiles,
             delta = round(unname(stats::quantile(disp, quantiles)), 6),
             n_displacements = length(disp), stringsAsFactors = FALSE)
}

# --- the selection rule ------------------------------------------------------

#' Score one method on the three declared axes, and apply the rule.
#'
#' THE STABILITY AXIS WAS REDEFINED ON 2026-09-06 AND THE LEVEL-4-OVER-LEVEL-3
#' RATIO IS DEMOTED. It is still computed and still reported; it no longer
#' decides anything. Four reasons, and the fourth is the one that settles it.
#'
#'   1. THE DENOMINATOR CAN BE SMALL, so the ratio explodes for exactly the arms
#'      whose sampling noise is lowest -- which is to say, the best ones.
#'   2. IT COLLAPSES TWO DISTRIBUTIONS TO ONE NUMBER, and the spread of each was
#'      the thing worth knowing. Both are now reported through
#'      `attr_stat_summary()` and compared through `attr_dominance()`.
#'   3. IT WAS COMPUTED FROM ONE PAIR OF REPLICATES. Replicate 1 against
#'      replicate 2 is one draw from a distribution with a p95 and a maximum,
#'      and a selection rule resting on one draw is a rule resting on an
#'      arbitrary index.
#'   4. LEVEL 3 AND LEVEL 4 ARE NOT ON ONE CONCEPTUAL AXIS. A large level 3
#'      means THIS METHOD IS NOT REPRODUCIBLE, which is a defect. A large level
#'      4 means TWO SPECIFICATIONS TELL DIFFERENT EXPLANATORY STORIES, which is
#'      a finding and is not a defect in either of them. Dividing one by the
#'      other produces a quantity whose numerator and denominator want opposite
#'      readings, and the frozen rule could not even say which pairs counted as
#'      level 4 (plan section 21.5). A rule that cannot name its own pair set is
#'      not made sound by choosing one after the numbers are visible.
#'
#' SO THE STABILITY AXIS IS NOW INTRINSIC: a method's OWN level-3 distribution,
#' scored at two points. `max_l3_median` stops an arm whose typical resample
#' already disagrees with itself; `max_l3_p95` stops an arm that is usually fine
#' and occasionally terrible. Both are needed -- a median alone lets one
#' horrible tail hide, and a p95 alone lets one extreme pair condemn an
#' otherwise reproducible method. The maximum is reported and is deliberately
#' NOT a threshold, because a maximum grows mechanically with B and a criterion
#' on it would tighten every time more replicates were generated.
#'
#' RESOLUTION AND DISCRIMINATION ARE UNCHANGED. Stability alone is still a trap:
#' a model that assigns everything to one signal is perfectly stable and
#' useless, which is what the resolution axis catches, and the score still has
#' to work, which is what discrimination catches.
#'
#' @param tab data frame, one row per method, with `method`, `l3_median`,
#'   `l3_p95`, `frac_unique_leader`, `auroc`.
#' @param rule the `selection` block of `config/attribution_eval.yml`.
attr_selection_score <- function(tab, rule) {
  need <- c("method", "l3_median", "l3_p95", "frac_unique_leader", "auroc")
  miss <- setdiff(need, names(tab))
  if (length(miss)) abort_values("attr_selection_score: missing column(s)", miss)
  max_med  <- as.numeric(cfg_req(rule, "max_l3_median"))
  max_p95  <- as.numeric(cfg_req(rule, "max_l3_p95"))
  min_res  <- as.numeric(cfg_req(rule, "min_resolution"))
  max_drop <- as.numeric(cfg_req(rule, "max_auroc_drop"))
  ref      <- as.character(cfg_req(rule, "auroc_reference"))
  if (!ref %in% tab$method) abort_values("attr_selection_score: auroc_reference absent", ref)

  a_ref <- tab$auroc[match(ref, tab$method)]
  # NA IS NOT A PASS, AND IT IS NOT A VERDICT EITHER. Two distinct things had to
  # be separated here. An unmeasured axis must never count as a cleared one --
  # that is how a method with no evidence behind it ends up looking admissible
  # -- so every `pass_*` is FALSE when its input is missing. But an arm scored
  # on two axes out of three has not been through the rule at all, and the first
  # version of this function let `rowSums` propagate the NA into `n_pass` and
  # then into `verdict`, so a partially generated replicate store printed a
  # column of `<NA>` verdicts that read as a failure to compute rather than as a
  # refusal to judge. `measured` is now explicit and the verdict is `unscored`.
  measured <- !is.na(tab$l3_median) & !is.na(tab$l3_p95) &
              !is.na(tab$frac_unique_leader) & !is.na(tab$auroc)
  z <- function(v) ifelse(is.na(v), FALSE, v)
  out <- data.frame(method = tab$method,
    l3_median = round(tab$l3_median, 5), l3_p95 = round(tab$l3_p95, 5),
    pass_stability = z(tab$l3_median <= max_med & tab$l3_p95 <= max_p95),
    resolution = round(tab$frac_unique_leader, 5),
    pass_resolution = z(tab$frac_unique_leader >= min_res),
    auroc = round(tab$auroc, 5), auroc_drop = round(a_ref - tab$auroc, 5),
    pass_discrimination = z((a_ref - tab$auroc) <= max_drop),
    all_axes_measured = measured,
    stringsAsFactors = FALSE)
  out$n_pass <- rowSums(out[, c("pass_stability", "pass_resolution",
                                "pass_discrimination")])
  # "Dominated on all three is out" -- the rule is stated as an exclusion, not
  # as a ranking, because a ranking would invite reading the top row as a
  # winner when the three axes are not commensurable.
  out$verdict <- ifelse(!measured, "unscored",
                 ifelse(out$n_pass == 0L, "excluded",
                 ifelse(out$n_pass == 3L, "admissible", "trade-off")))
  out[order(!out$all_axes_measured, -out$n_pass, -out$auroc), , drop = FALSE]
}

#' The demoted ratio, computed for the record and marked as reporting-only.
#'
#' KEPT RATHER THAN DELETED because every number reported before 2026-09-06 was
#' on this scale and a reader comparing the two documents needs the bridge. It
#' is computed from the two DISTRIBUTIONS now rather than from one pair, so the
#' `_median` and `_p95` columns say which point of each is being divided, and
#' `attr_dominance()` beside it is the quantity to actually read.
#'
#' @param l4,l3 numeric vectors of disagreements.
attr_stability_ratio <- function(l4, l3) {
  l4 <- l4[is.finite(l4)]; l3 <- l3[is.finite(l3)]
  q <- function(v, p) if (length(v)) unname(stats::quantile(v, p)) else NA_real_
  m4 <- q(l4, 0.5); m3 <- q(l3, 0.5); p4 <- q(l4, 0.95); p3 <- q(l3, 0.95)
  data.frame(n_l4 = length(l4), n_l3 = length(l3),
             ratio_median = round(m4 / pmax(m3, 1e-9), 4),
             ratio_p95 = round(p4 / pmax(p3, 1e-9), 4),
             status = "reporting_only_demoted_20260906",
             stringsAsFactors = FALSE)
}

#' Leader collapse: where does ONE arm's most important signal land in the other?
#'
#' WHAT THE REST OF THE FAMILY CANNOT SEE. Every other metric here is either
#' TOP-FOCUSED (`attr_agree_k`, `attr_tieset`) or a WHOLE-VECTOR AGGREGATE
#' (`attr_cosine`, `attr_sign_flip`). Between them they report HOW MUCH two arms
#' disagree and nothing about WHERE the disagreement lives. Those are different
#' questions and the difference is consequential.
#'
#' MEASURED 2026-09-07, and it is the reason this exists. `llr_full` against
#' `llr_cond` has a top-1 disagreement of 0.470 -- 47% of patients get a
#' different leading signal. That single number is equally consistent with two
#' stories: every leader nudged aside by its near-tie neighbour, or a minority
#' of leaders collapsing outright. It is the second. Median rank displacement is
#' ZERO (53% keep the same leader) while the 95th percentile is 12 ranks of
#' about 18, and 3.3% of patients see their `full` leader land in `cond`'s
#' bottom three. `llr_full` against `llr_full_ti_trend` has a top-1 disagreement
#' of 0.136 and a displacement p95 of ONE rank, with a bottom-three rate BELOW
#' the `cond` arm's own sampling noise. Two "level-4 specification differences"
#' of the same kind by every existing metric, and completely different in
#' character.
#'
#' Cosine partly sees this (0.106 against 0.0035) but cannot separate "3% of
#' patients collapse catastrophically" from "everyone shifted mildly", because
#' it is an average over the whole vector. This metric LOCALISES: it says
#' whether the disagreement is concentrated or diffuse.
#'
#' --- THE MAGNITUDE FORM IS PRIMARY, AND THE RANK FORM IS THE DIAGNOSTIC ---
#'
#' `share_*` is the share of the OTHER arm's own kept evidence carried by this
#' arm's leader. It needs no ranking, no tie-breaking, no bottom-boundary
#' convention and no minimum-measured-count filter, which is why it leads.
#'
#' A LITERAL "LAST RANK" METRIC WAS CONSIDERED AND REJECTED. For an LLR arm the
#' bottom of the ranking is mostly ASSIGNED ZEROS -- an unmeasured signal gets
#' `L = 0` by construction (CLAUDE.md frozen decisions) and among tied zeros the
#' order falls out of the column index, identically in both arms. So "the
#' last-ranked signal" would frequently be a fact about the signal vocabulary's
#' column order rather than about the model. SHAP has no assigned zeros -- a
#' tree routes a missing feature through its default direction and gives it a
#' real contribution -- so its last rank is a genuine but tiny-magnitude entry,
#' where signal-to-noise is worst by construction. The two families would
#' differ for reasons that have nothing to do with attribution quality. Ranking
#' here is therefore restricted to KEPT cells and reported as displacement
#' rather than as an absolute position.
#'
#' ASYMMETRIC ON PURPOSE. `A`'s leader may collapse in `B` while `B`'s leader
#' survives in `A`; the two directions are different facts and both are
#' returned. For a within-method contrast they should agree by symmetry, which
#' makes their difference a free consistency check.
#'
#' @param keep optional logical cell mask, as in `attr_sign_flip()`. For the LLR
#'   arms this is the measured mask, so the denominator is the patient's own
#'   MEASURED evidence and an assigned zero can never be a leader. SHAP has no
#'   analogue and legitimately passes NULL, which the `n_cells` column surfaces.
#' @param shares thresholds on the carried share. A signal carrying less than
#'   `1/p` of the evidence is already below the flat-vector average, so the
#'   declared grid should sit well under that.
#' --- A PERFORMANCE DEFECT, AND WHY THE FIX WAS TO COMPUTE LESS -----------
#'
#' Everything here is `max.col`, `rowSums` and one vectorised comparison, so it
#' affords the full pair set the way `attr_cosine()` does. IT DID NOT START
#' THAT WAY. The rank displacement was computed by building two complete
#' n x p rank matrices with `t(apply(M, 1, rank))` -- 19 ranks per patient, to
#' read one of them -- which is `attr_prep()`'s cost class, precisely what the
#' full pair-set passes exist to avoid.
#'
#' MEASURED, twice, because the first fix was the wrong one. Unguarded, it was
#' 12,817 within-method pairs times two aggregation levels times two rank
#' sweeps of 41,250 rows: 2.11 BILLION `rank()` calls, which took the consumer
#' from 17 minutes to over 24 and still running. Guarding it behind a `ranks`
#' switch cut the hot path but left 7.4 minutes in a 56-pair anchor table whose
#' own comment called the cost "irrelevant". Only then was the actual question
#' asked -- what is the rank FOR? -- and the answer was one number per patient,
#' obtainable by counting entries that outrank it.
#'
#' THE CLASS: OPTIMISING THE GUARD RATHER THAN THE COMPUTATION. A switch that
#' turns an expensive path off is a way of not asking whether the path needed
#' to be expensive. The switch is gone; the diagnostic is now free and is always
#' computed.
#'
#' SPLIT INTO CELLS AND SUMMARY ON 2026-09-07. `attr_leader_cells()` computes
#' the per-patient quantities once; `attr_leader_collapse()` summarises them to
#' one row per threshold exactly as before, and
#' `attr_leader_patient_summary()` reports the same quantities as DISTRIBUTIONS
#' OVER PATIENTS with the fixed `attr_stat_summary()` column set. The split
#' exists because `leader_collapse.csv` carried a median and a p95 of the rank
#' displacement and nothing else, so the shape of that distribution -- which is
#' the whole point of the metric, concentrated collapse against diffuse drift --
#' was not reported, and an apply site with no replicates had no distributional
#' form of it at all.
#'
#' @return `attr_leader_cells()`: a list of ROW-LEVEL vectors (hard rule 1:
#'   summarise, never print). `attr_leader_collapse()`: one row per threshold;
#'   the rank columns repeat across rows.
attr_leader_cells <- function(A, B, keep = NULL) {
  stopifnot(identical(dim(A), dim(B)))
  n <- nrow(A); p <- ncol(A)
  if (is.null(keep)) keep <- matrix(TRUE, n, p)
  # `abs()` ONCE EACH. The first version called `abs(A)` and `abs(B)` four times
  # between them and reached the row totals through `ifelse(kept, abs(A), 0)`,
  # which allocates two 783,750-cell matrices per call and is the slowest way in
  # R to zero a masked matrix. Multiplying by the logical mask coerces it to 0/1
  # and does the same job in one pass. Measured together with the `ranks` guard:
  # 0.228 s per pair down to well under a tenth of that.
  aA <- abs(A); aB <- abs(B)
  aa <- aA; bb <- aB
  aa[!keep] <- -Inf; bb[!keep] <- -Inf          # a dropped cell can never lead

  # THE `-Inf` IS LOAD-BEARING AND WAS `Inf` IN THE FIRST DRAFT OF THE PROBE
  # THAT MOTIVATED THIS FUNCTION. `max.col` returns the position of the LARGEST
  # value, so substituting `Inf` for a masked cell makes every masked cell win
  # and the leader is drawn from precisely the cells that were meant to be
  # excluded. It failed loudly there only because the downstream rank lookup
  # then returned NA and the scored count collapsed; with a different downstream
  # it would have returned a confident wrong answer.
  la <- max.col(aa, ties.method = "first")      # A's leader, per patient
  lb <- max.col(bb, ties.method = "first")      # B's leader, per patient
  ia <- cbind(seq_len(n), la); ib <- cbind(seq_len(n), lb)

  kept  <- keep
  tot_a <- rowSums(aA * kept)
  tot_b <- rowSums(aB * kept)
  # Share of the OTHER arm's evidence carried by THIS arm's leader.
  sh_a_in_b <- aB[ia] / pmax(tot_b, .Machine$double.eps)
  sh_b_in_a <- aA[ib] / pmax(tot_a, .Machine$double.eps)
  ok <- rowSums(kept) >= 2L & tot_a > 0 & tot_b > 0

  # Rank displacement, among KEPT cells only, WITHOUT BUILDING A RANK MATRIX.
  #
  # ONLY ONE RANK PER PATIENT IS WANTED -- the position of the other arm's
  # leader -- and that is simply a count of how many entries outrank it. One
  # vectorised comparison per direction. The first version called
  # `t(apply(M, 1, rank))` twice, computing all 19 ranks for every patient in
  # order to read one of them: 8.1 seconds per pair, which put 7.4 minutes into
  # a 56-pair table described in its own comment as one "where the cost is
  # irrelevant". It is under a tenth of a second now, and the `ranks` switch
  # that existed to hide the cost is gone with it.
  #
  # `bb > bb[ia]` recycles the length-n leader vector down each column, so
  # element (i, j) is compared against patient i's own leader value -- masked
  # cells hold -Inf and can never outrank a kept one. The convention is MIN
  # RANK (count of strictly greater) rather than `ties.method = "first"`; on
  # continuous L values exact ties are vanishing, and for a displacement
  # diagnostic the min-rank reading is the more natural one anyway.
  nk <- rowSums(kept)
  disp_a <- rowSums(bb > bb[ia])                 # signals outranking A's leader in B
  disp_b <- rowSums(aa > aa[ib])

  # `same_leader` IS THE CONDITIONING VARIABLE FOR THE SECOND POPULATION
  # (review finding A8): the summaries below are reported over all scored
  # patients AND over the patients whose leader actually differs, and a
  # first-column tie policy decides equality here exactly as it decides `la`
  # and `lb`.
  list(ok = ok, nk = nk, n_cells = sum(keep),
       sh_a_in_b = sh_a_in_b, sh_b_in_a = sh_b_in_a,
       disp_a = disp_a, disp_b = disp_b, same_leader = la == lb,
       n_patients = n, n_pairs = 1L)
}

attr_leader_collapse <- function(A, B, keep = NULL, shares) {
  cl <- attr_leader_cells(A, B, keep)
  ok <- cl$ok; nk <- cl$nk
  sh_a_in_b <- cl$sh_a_in_b; sh_b_in_a <- cl$sh_b_in_a
  disp_a <- cl$disp_a; disp_b <- cl$disp_b

  q <- function(v, pr) if (any(ok)) unname(stats::quantile(v[ok], pr)) else NA_real_
  rt <- function(v) if (any(ok)) round(mean(v[ok] >= 2 * nk[ok] / 3), 5) else NA_real_
  data.frame(
    share_threshold = shares,
    n_scored = sum(ok), n_cells = cl$n_cells,
    collapse_a_in_b = round(vapply(shares, function(s) mean(sh_a_in_b[ok] < s), numeric(1)), 5),
    collapse_b_in_a = round(vapply(shares, function(s) mean(sh_b_in_a[ok] < s), numeric(1)), 5),
    share_a_in_b_median = round(q(sh_a_in_b, 0.50), 5),
    share_a_in_b_p05 = round(q(sh_a_in_b, 0.05), 5),
    share_b_in_a_median = round(q(sh_b_in_a, 0.50), 5),
    share_b_in_a_p05 = round(q(sh_b_in_a, 0.05), 5),
    disp_a_median = q(disp_a, 0.50), disp_a_p95 = q(disp_a, 0.95),
    disp_b_median = q(disp_b, 0.50), disp_b_p95 = q(disp_b, 0.95),
    disp_a_p90 = q(disp_a, 0.90), disp_b_p90 = q(disp_b, 0.90),
    # Rate rather than a boundary count, so it means the same thing for a
    # patient with 8 kept signals and one with 18.
    frac_a_into_bottom_third = rt(disp_a),
    frac_b_into_bottom_third = rt(disp_b),
    stringsAsFactors = FALSE)
}

#' Concatenate the scored entries of several `attr_leader_cells()` results.
#'
#' THIS IS HOW A MODEL-AGAINST-ITSELF DISTRIBUTION IS BUILT. One bag pair gives
#' one displacement per patient; pooling the scored entries over every sampled
#' bag pair gives the distribution of "how far does this patient's leader move
#' when the model is refitted on another plausible sample", over patients AND
#' pairs at once. Its median, quartiles and tail are on the rank scale and read
#' directly, which the per-pair median cannot: that is identically zero for
#' every LLR pair because most patients keep their leader, so a distribution
#' over pairs of it says nothing, and a per-pair MEAN of a rank is not a
#' quantity anyone can report.
.pool_leader_cells <- function(cls) {
  pool <- .pool_leader_new()
  for (cl in cls) pool <- .pool_leader_add(pool, cl)
  pool
}

#' An empty pool, and the append that grows it one pair at a time.
#'
#' INCREMENTAL AS OF 2026-09-09 (review finding A12). The consumer kept every
#' pair's FULL cells list -- seven vectors over all 41,250 patients -- for
#' every pair of every cell until the end of the run, and only then pooled the
#' scored entries. Keeping only the scored entries, appended as each pair is
#' finished and summarised as soon as its cell is complete, bounds the memory
#' to one cell's pooled vectors. `.pool_leader_extract()` is the compact form
#' a caller can write to local scratch when the loop order forces a cell to be
#' finished much later than it is started (the bag-outer level-4 loop).
.POOL_FIELDS <- c("ok", "nk", "sh_a_in_b", "sh_b_in_a", "disp_a", "disp_b",
                  "same_leader")

.pool_leader_new <- function() {
  list(chunks = list(), n_patients = 0L, n_pairs = 0L, n_cells = NA_integer_,
       compact = TRUE)
}

#' Concatenate a pool's chunks into flat vectors, once.
#'
#' CHUNKED, NOT GROWN, because the first incremental version appended with
#' `c()` on every pair and a 150-pair cell became a quadratic copy: the
#' cosine-and-leader pass took 69 minutes where the list-pooled version took
#' a few. One concatenation per field at the end is linear.
.pool_leader_flatten <- function(pool) {
  if (is.null(pool$chunks)) return(pool)
  out <- pool[c("n_patients", "n_pairs", "n_cells")]
  for (f in .POOL_FIELDS) {
    out[[f]] <- unlist(lapply(pool$chunks, `[[`, f), use.names = FALSE)
  }
  out$compact <- TRUE
  out
}

.pool_leader_extract <- function(cl) {
  ok <- cl$ok
  list(ok = rep(TRUE, sum(ok)), nk = cl$nk[ok], n_patients = length(ok),
       n_pairs = 1L, n_cells = NA_integer_,
       sh_a_in_b = cl$sh_a_in_b[ok], sh_b_in_a = cl$sh_b_in_a[ok],
       disp_a = cl$disp_a[ok], disp_b = cl$disp_b[ok],
       same_leader = cl$same_leader[ok], compact = TRUE)
}

.pool_leader_add <- function(pool, cl) {
  x <- if (isTRUE(cl$compact)) cl else .pool_leader_extract(cl)
  if (is.null(pool$chunks)) pool <- .pool_leader_new()
  if (!is.null(x$chunks)) x <- .pool_leader_flatten(x)
  pool$chunks[[length(pool$chunks) + 1L]] <- x[.POOL_FIELDS]
  pool$n_patients <- pool$n_patients + x$n_patients
  pool$n_pairs    <- pool$n_pairs + x$n_pairs
  pool
}

#' The leader-collapse quantities as DISTRIBUTIONS OVER PATIENTS, from a cells
#' list -- a single pair's, or a pooled one from `.pool_leader_cells()`.
#'
#' Four metrics: the rank displacement of A's leader inside B and of B's leader
#' inside A, and the share of the other arm's kept evidence each leader
#' carries. Each row carries the fixed `attr_stat_summary()` column set over
#' the scored entries, `frac_scored` (a patient is scored when both arms have
#' at least two kept cells and non-zero evidence), and for the displacement
#' rows `frac_bottom_third`, the rate at which the leader lands in the bottom
#' third of the other arm's kept ranking. The share rows carry the collapse
#' rate at each declared threshold in `collapse_<share>` columns.
#'
#' TWO POPULATIONS PER METRIC AS OF 2026-09-09 (review finding A8). The
#' scored set includes every patient whose leader did NOT move, so a
#' displacement distribution over it has a median of zero by construction
#' whenever most leaders survive, and it answers "how far does the typical
#' leader move" rather than "where does a leader go WHEN it moves". Both are
#' worth having and they were not distinguishable before: `population` is
#' now `all_scored` or `leader_differs`, the second conditioned on the two
#' arms naming different leaders under the `tie_policy` column's rule
#' (first column among exact ties, the same rule `max.col` applies). A
#' pooled cell counts PATIENT-PAIR observations rather than distinct patients,
#' and says so in `unit`.
#'
#' ONE SCHEMA FOR EVERY SITE AND EVERY CONTRAST. An apply site has one fit per
#' arm and no distribution over replicates, but it has 95,507 patients; the
#' training site has 38 shared bags and pools over them. Both are distributions
#' over patients and carry the same columns, so the tables read side by side.
#' What differs -- one fit or many -- is recorded in the caller's `held`
#' column, never inferred from the numbers.
attr_leader_distribution <- function(cl, shares) {
  cl <- .pool_leader_flatten(cl)
  ok_all <- cl$ok; nk <- cl$nk
  n_pat <- if (is.null(cl$n_patients)) length(ok_all) else cl$n_patients
  n_pairs <- if (is.null(cl$n_pairs)) 1L else cl$n_pairs
  unit <- if (n_pairs > 1L) "patient_pairs" else "patients"
  same <- cl$same_leader
  if (is.null(same)) same <- rep(FALSE, length(ok_all))
  pops <- list(all_scored = ok_all, leader_differs = ok_all & !same)
  one <- function(metric, direction, v, is_disp, pop, ok) {
    x <- v; x[!ok] <- NA_real_
    out <- cbind(
      data.frame(metric = metric, direction = direction, population = pop,
                 unit = unit, tie_policy = "first_column",
                 n_patients = n_pat, n_pairs = n_pairs, n_scored = sum(ok),
                 frac_scored = round(sum(ok) / n_pat, 5),
                 frac_leader_differs = round(if (any(ok_all))
                   mean(!same[ok_all]) else NA_real_, 5),
                 stringsAsFactors = FALSE),
      attr_stat_summary(x),
      data.frame(frac_bottom_third = if (is_disp && any(ok))
                   round(mean(v[ok] >= 2 * nk[ok] / 3), 5) else NA_real_,
                 stringsAsFactors = FALSE))
    for (s in shares) {
      out[[sprintf("collapse_%g", s)]] <- if (!is_disp && any(ok))
        round(mean(v[ok] < s), 5) else NA_real_
    }
    out
  }
  rows <- list()
  for (pop in names(pops)) {
    ok <- pops[[pop]]
    rows[[length(rows) + 1L]] <- rbind(
      one("rank_displacement", "a_leader_in_b", cl$disp_a, TRUE, pop, ok),
      one("rank_displacement", "b_leader_in_a", cl$disp_b, TRUE, pop, ok),
      one("share_of_other_evidence", "a_leader_in_b", cl$sh_a_in_b, FALSE, pop, ok),
      one("share_of_other_evidence", "b_leader_in_a", cl$sh_b_in_a, FALSE, pop, ok))
  }
  do.call(rbind, rows)
}

#' Convenience: the distribution table for one pair of matrices.
attr_leader_patient_summary <- function(A, B, keep = NULL, shares) {
  attr_leader_distribution(attr_leader_cells(A, B, keep), shares)
}

#' The per-pair scalars a pair-level distribution is built from. No mean of a
#' rank anywhere: the displacement is summarised by its median, p90 and p95.
attr_leader_pair_scalars <- function(cl, shares) {
  ok <- cl$ok
  q <- function(v, pr) if (any(ok)) unname(stats::quantile(v[ok], pr)) else NA_real_
  rt <- function(v) if (any(ok)) mean(v[ok] >= 2 * cl$nk[ok] / 3) else NA_real_
  cr <- function(v) vapply(shares, function(s) if (any(ok)) mean(v[ok] < s) else NA_real_, numeric(1))
  c(stats::setNames(cr(cl$sh_a_in_b), sprintf("collapse_ab_%g", shares)),
    stats::setNames(cr(cl$sh_b_in_a), sprintf("collapse_ba_%g", shares)),
    disp_ab_median = q(cl$disp_a, 0.5), disp_ab_p90 = q(cl$disp_a, 0.9),
    disp_ab_p95 = q(cl$disp_a, 0.95),
    disp_ba_median = q(cl$disp_b, 0.5), disp_ba_p90 = q(cl$disp_b, 0.9),
    disp_ba_p95 = q(cl$disp_b, 0.95),
    bottom3rd_ab = rt(cl$disp_a), bottom3rd_ba = rt(cl$disp_b),
    share_ab_median = q(cl$sh_a_in_b, 0.5), share_ab_p05 = q(cl$sh_a_in_b, 0.05),
    share_ba_median = q(cl$sh_b_in_a, 0.5), share_ba_p05 = q(cl$sh_b_in_a, 0.05))
}
