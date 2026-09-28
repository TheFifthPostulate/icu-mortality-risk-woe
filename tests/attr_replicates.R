# tests/attr_replicates.R -----------------------------------------------------
# THE EXPENSIVE GENERATOR, AND THE ONLY SCRIPT IN THIS SUB-PIPELINE THAT FITS.
#
# Read `docs/attribution_analysis_plan_20260906.md` section 10, then PART FOUR
# (section 25 onward), which records the 2026-09-06 redesign. It builds
# ATTRIBUTION REPLICATES
#
#     (method, boot_id, seed_id, draw_id)  ->  an N x 19 matrix in nats, oof
#
# writes one content-addressed file per replicate plus a manifest, and skips any
# replicate whose key already exists. `tests/attr_metrics.R` reads what this
# writes and computes every metric in minutes.
#
# --- THE SHARED BOOTSTRAP MANIFEST, WHICH IS THE POINT OF THE REDESIGN -------
#
# Bag `b` is `bag_of(bootstrap.seed_base + b)` and EVERY METHOD THAT RESAMPLES
# USES THAT BAG for `boot_id = b`. Until 2026-09-06 the SHAP bootstrap used
# `bag_of(seed_base + i)` and the LLR bootstrap used
# `bag_of(seed_base + 90000 + i)`, so SHAP replicate 3 and LLR replicate 3 were
# fitted on DIFFERENT resamples. Nothing wrong followed from that as long as
# each method's level 3 was only ever read on its own; but it made two things
# impossible, and both are things the arm needs.
#
#   A LEVEL-4 DISTRIBUTION. Level 4 is specification: arm A against arm B. With
#   independent resamples the only honest level-4 number is the one from the two
#   ORIGINAL fits -- a single value. Fit A and B on the SAME bag b and the
#   comparison is paired, the variation caused by which patients entered bag b
#   is shared and differences out, and across bags you get {D(A_b, B_b)}: a
#   distribution with a median, a p95 and a worst case.
#
#   A PAIRED CROSS-METHOD LEVEL 3. For one pair of bags (b, b'), SHAP's sampling
#   disagreement and an LLR arm's sampling disagreement are now measured on
#   identical resamples, so they can be differenced bag-pair by bag-pair rather
#   than compared as two marginal distributions. That is the comparison the
#   paper wants: not "does LLR agree with SHAP", which there is every reason it
#   should not, but "under identical training-sample perturbations, whose
#   patient-level explanation moves less".
#
# `bootstrap.seed_base` IS UNCHANGED from the value the 40 SHAP bootstrap
# replicates on disk were generated under, deliberately, so bag b is the bag it
# always was and those replicates are still valid. The five OLD LLR bootstrap
# replicates are NOT preserved: they were on bags `seed_base + 90000 + i`, which
# are different bags, and a replicate on a different bag is a different
# replicate however convenient it would be to keep.
#
# --- THE SIX ROUTES, AND WHY A REPLICATE BELONGS TO EXACTLY ONE -------------
#
#   ladder            (0,0,0). The arm as fitted out of fold. The five base L
#                     matrices come from the targets cache and a `coupattr`
#                     run whose design key matches; the three `cond` arms are
#                     subtractions. Cheap.
#   seed              (0,s,0). SHAP only. `attr_has_seed_noise()` is FALSE for
#                     every LLR arm because `bam` with fREML is deterministic at
#                     fixed settings (audit finding F9 recorded a bitwise
#                     identical refit), so a re-seeded GAM replicate is not an
#                     observation and is not generated. That zero is a
#                     REPORTABLE PROPERTY and the plan table carries a
#                     `route: none` row for it rather than omitting the method.
#   bootstrap         (b,0,0). A refit on shared bag b. Both families.
#   bootstrap_seeded  (b,s,0). SHAP only. THE LEVEL-2-INTO-LEVEL-3 PROPAGATION:
#                     bag b carries a second algorithmic draw, so seed noise and
#                     sampling noise can be varied TOGETHER as well as
#                     separately. See the next block.
#   posterior         (0,0,d). LLR only. Estimation uncertainty conditional on
#                     the observed sample.
#   none              a marker for a level a method provably cannot have.
#
# --- HOW LEVEL 2 IS PROPAGATED INTO LEVEL 3, WITHOUT A VARIANCE MODEL -------
#
# The question is what happens to a SHAP attribution when you retrain on another
# plausible sample. In real use that changes two things at once: the training
# rows AND the booster's random draw. Until 2026-09-06 the two were measured
# separately and never together -- level 2 varied the seed at a fixed sample,
# level 3 varied the sample at a fixed seed -- so neither number was the
# operational quantity and there was no way to combine them without assuming a
# variance decomposition that the response variable (a per-patient categorical
# top-1 choice) does not support.
#
# `levels.seed.b_per_bag` fixes it BY DESIGN rather than by a model. Each shared
# bag gets one additional re-seeded fit, and the same store then yields all
# three:
#
#   L2   (b,0) vs (b,1)     seed alone, and now conditional on a bag, so the
#                           question "is seed noise bigger on some resamples
#                           than others" is answerable rather than assumed away
#   L3   (b,s) vs (b',s)    sample alone, seed held
#   L3T  (b,0) vs (b',1)    BOTH -- the operational quantity
#
# A Bayesian hierarchical model over pairwise disagreement indicators would need
# careful treatment of the dependence that pairwise comparison induces, and it
# would be a model used to rescue an ambiguous design. A nested design measures
# all three directly, and a variance decomposition, if wanted later, becomes a
# second-stage summary of a good experiment rather than a substitute for one.
#
# THE COST IS `bootstrap.b` EXTRA SHAP FITS -- 40 boosters, about 28 minutes.
#
# --- WHY POSTERIOR SIMULATION IS KEPT, AND WHAT IT IS NOT --------------------
#
# A single LLR bootstrap replicate refits 320 out-of-fold spec-folds and
# measured about 7.4 minutes; the 320-fit basis pass for the posterior route
# cost under ten minutes ONCE, after which every replicate is matrix algebra. So
# 40 posterior replicates cost less than two bootstrap ones.
#
# THEY ARE NOT THE SAME QUANTITY, AND THE 2026-09-06 MEASUREMENT REVERSED THE
# PREDICTED DIRECTION. Posterior simulation conditions on the observed sample
# and propagates ESTIMATION uncertainty; a bootstrap resamples and propagates
# SAMPLING uncertainty. The posterior spread was predicted to be the narrower of
# the two and MEASURED AT 1.17 TO 1.57 TIMES THE BOOTSTRAP'S. The mechanism is
# the penalty: every layer-1 smooth is `bs = "ts"` under `gamma = 1.5`, so `Vc`
# describes how far the coefficients COULD be while a bootstrap describes how
# far the penalised PROCEDURE actually moves, and under heavy shrinkage the
# second is smaller. The nominal posterior width is not wrong; it answers a
# different question, and the question level 3 asks is the bootstrap's.
#
# SO THE BOOTSTRAP ROUTE IS PRIMARY FOR THE LLR ARMS AS OF 2026-09-06, and the
# posterior route is retained as the cheap high-B companion. `attr_pair_contrast()`
# gives them different contrast codes (L3 against L3P) and the consumer groups
# by route, so the two can never be pooled into one noise scale.
#
# --- THE ONE APPROXIMATION IN THE DERIVED ARMS, STATED ----------------------
#
# `cond = full - intv` is a difference of TWO SEPARATELY FITTED MODELS, and a
# POSTERIOR draw perturbs each independently. Their estimation errors are in
# fact positively correlated -- `intv`'s term set is a subset of `full`'s and
# both are fitted on the same rows -- so an independent draw OVERSTATES the
# `cond` arms' estimation noise by roughly the covariance it discards. Measured:
# the `cond` arms' posterior-to-bootstrap ratio is 1.42 to 1.57 against 1.17 to
# 1.21 for `meas` and the `full` arms, and on the bootstrap route `llr_cond` is
# QUIETER than `llr_full` while on the posterior route it is noisier, which is
# the reverse of what removing a source of variation should do.
#
# THE BOOTSTRAP REPLICATES DO NOT HAVE THIS PROBLEM -- a bootstrap refits both
# models on the same resample, so their errors covary correctly -- which is the
# third reason the bootstrap route is now primary.
#
# --- THE RESAMPLE IS 63.2% OF PATIENTS, NOT A FULL BOOTSTRAP ----------------
#
# `bag_of()` draws patients WITH replacement and then keeps the DISTINCT ones,
# so about 63.2% of patients are in bag. It cannot do better: `.frame_rows()`
# in `R/04_features.R` selects with `%in%`, so a stay_id repeated in the fit set
# is silently deduplicated, and expressing a true with-replacement bootstrap
# would mean changing a load-bearing pipeline function to serve a diagnostic.
#
# THE CONSEQUENCE, WITH ITS DIRECTION. An m-out-of-n subsample has roughly
# sqrt(n/m) = 1.26 times the sampling spread of a full-size resample, so every
# level-3 number produced this way is INFLATED by about that factor.
#
# IT IS NOW SYMMETRIC BETWEEN THE FAMILIES, WHICH IT WAS NOT BEFORE. The
# subsample applies to every SHAP level-3 replicate and to every LLR BOOTSTRAP
# replicate, and both now draw from the SAME bags, so the cross-family level-3
# comparison is like for like and the inflation is common to both sides of it.
# The posterior replicates resample nothing and are therefore NOT comparable
# with SHAP; they are the within-method high-B companion and the consumer keeps
# them in their own contrast.
#
# --- WHAT A BOOTSTRAP REPLICATE HERE DOES NOT RESAMPLE ----------------------
#
# EVERY REPLICATE HOLDS THE FROZEN PRIOR PARAMETERS FIXED. `alpha`, the `delta`
# coefficients and the `lambda` coefficients are taken from the pipeline's own
# per-fold `priors` object and are NOT re-estimated on the resample, for either
# family. So a level-3 replicate measures the sampling variability of LAYER 1
# CONDITIONAL ON THE COVARIATE CONSTRUCTION, not of the whole stack.
#
# THAT IS DELIBERATE AND IS NOT THE SAME AS `attribution_floor()`, which refits
# the whole stack (`R/13b_attribution_floor.R`). Holding priors fixed for BOTH
# families keeps the SHAP and LLR level-3 numbers on the same footing, which is
# the entire reason level 3 exists. THE DIRECTION OF THE OMISSION: it
# understates level-3 spread for both families by whatever the prior estimation
# contributes.
#
# --- HARD RULE 6 -------------------------------------------------------------
#
# NOT BREACHED, AND IT DID NOT NEED TO BE RELAXED. Rule 6 forbids a FOLD FIT'S
# GAM OBJECT being kept: `fit_one()` returns predictions and diagnostics and
# never the model. Nothing here keeps one either. The fold fit lives inside one
# iteration of a loop, its `beta` and `Vc` are consumed into posterior draws in
# that same iteration, and what reaches the disk is an L matrix -- the same kind
# of object `l_oof.rds` already is. Nothing enters `_targets/` and nothing
# enters a bundle.
#
# --- RESUMABILITY AND IMMUTABILITY ------------------------------------------
#
# `--resume <dir>` appends to an existing generator run, which is the one place
# in this project where a run directory is written to twice. It is safe for a
# specific reason rather than by exception: every replicate file is named by
# `attr_replicate_key()`, a content hash of the method, stage, index and the
# full design. A key that already exists cannot be produced with different
# content, so appending is MONOTONE -- it can add files, never change one. A
# design change produces different keys and therefore a different set of files.
#
# THE MANIFEST IS AN UPSERT AS OF 2026-09-06, AND THAT IS A DELIBERATE
# EXCEPTION TO THE MONOTONICITY ABOVE. The replicate coordinates
# `boot_id`/`seed_id`/`draw_id` are new columns, so a resumed run must be able
# to write them onto rows whose FILES are already correct and already skipped.
# The upsert rewrites metadata about a replicate and never the replicate, and
# the key -- which is a hash of what produced the matrix -- is what it keys on.
#
# Aggregates only on print (hard rule 1). The replicate matrices and the bag
# memberships are row-level, are saved the way `l_oof.rds` already is, and are
# never printed; the bag table reports counts and a membership hash.
#
#   Rscript tests/attr_replicates.R --plan-only
#   Rscript tests/attr_replicates.R --check-specs     # seconds, fits nothing
#   Rscript tests/attr_replicates.R --levels spec,seed
#   Rscript tests/attr_replicates.R --levels sample --resume out/runs/attrgen_...
#   Rscript tests/attr_replicates.R --levels sample --resume <dir> --retry-excluded
#
# RUN `--check-specs` BEFORE `--levels sample`. It walks every (arm, signal)
# through the same `arm_formula()` the fitting loops use, asserts the enumerated
# spec-fold count matches the plan, and takes seconds. The five-hour route is
# the only thing that exercises those loops otherwise, which is how a guard
# left one statement too late got found by a five-hour run rather than by a
# five-second one.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(targets); library(mgcv); library(xgboost); library(arrow); library(yaml)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

# --- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
LEVELS_WANTED <- strsplit(.opt("--levels", "spec,seed"), ",", fixed = TRUE)[[1]]
RESUME    <- .opt("--resume", NULL)
PLAN_ONLY <- "--plan-only" %in% args
# Exercises every (arm, signal) through `arm_formula()` and fits nothing. It
# lives behind a flag rather than running unconditionally because it needs
# `spec_of()`, which is defined inside the `sample` block; see `check_specs()`.
CHECK_SPECS <- "--check-specs" %in% args
# Clears the excluded-bag tombstone so bags recorded as unsupported are
# attempted again. For use after something that would change the outcome; a bag
# excluded for `basis_exceeds_distinct_values` under an unchanged design will
# simply be excluded again, at about seven minutes each.
RETRY_EXCL <- "--retry-excluded" %in% args
# Refit the SHAP ladder booster (about a minute and a half) and assert the
# stored ladder replicate is bitwise what the live design produces, even when
# its margin file already exists. The refit happens without the flag whenever
# the margin file is absent, which is every store written before 2026-09-09.
VERIFY_SHAP <- "--verify-shap" %in% args

bad <- setdiff(LEVELS_WANTED, ATTR_LEVELS)
if (length(bad)) abort_values("--levels must be a subset of spec,seed,sample", bad)

# --- the design --------------------------------------------------------------
ecfg   <- yaml::read_yaml("config/attribution_eval.yml")
cfg    <- tar_read(cfg)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
y      <- as.integer(tar_read(y_train))
priors <- tar_read(priors)
Lz     <- tar_read(l_mats_zero)
tabs   <- tar_read(tabs)

# `k_ti` AND THE LADDER PATH COME FROM CONFIG AS OF 2026-09-06. Both were
# literals in this file -- `K_TI <- 5L` and a mandatory `--ladder` argument the
# operator had to remember -- and `k_ti` was ALSO a literal in
# `tests/coupling_attribution.R`. Two literals for one quantity is the F9-to-F11
# pattern, and it is worse than usual here because `k_ti` is a FIELD OF THE
# DESIGN KEY: the two scripts disagreeing about it produces two design keys and
# fires the staleness guard on a difference nobody declared.
# `k_ti` MOVED TO config/config.yml ON 2026-09-07, when the interaction
# models entered LAYER1_MODELS and the cap became something the PIPELINE
# fits with. `config/attribution_eval.yml` promises in its own header that
# no key in it changes what is fitted, and that promise was about to become
# false. Read through `cfg` with no fallback, exactly as before.
K_TI      <- as.integer(cfg_req(cfg, "k_ti"))
.ladder_cfg <- as.character(cfg_req(ecfg, "ladder", "run"))
LADDER_D  <- .opt("--ladder", if (nzchar(.ladder_cfg)) .ladder_cfg else NA_character_)
if (is.na(LADDER_D) || !nzchar(LADDER_D)) LADDER_D <- NULL

sigs    <- as.character(unlist(cfg$signals))
paired  <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, sigs)
n_fold  <- cfg$n_folds %||% 5L
ids_ch  <- as.character(tr)
fold_k  <- folds$fold[match(tr, folds$stay_id)]
# The boosters' inner early-stopping split groups by patient (statistical
# review S3, 2026-09-09), resolved through the same accessor the folds use.
grp_k   <- patient_group_of(tabs$cohort, cfg, tr)
meas_ok <- measured_matrix(tabs, cfg, tr)

DESIGN  <- attr_design_key(cfg, fold_k, K_TI)
# THE WIDER IDENTITY, CHECKED BESIDE THE KEY (review findings A1, A2). The key
# names the files; the fingerprint says whether the design that would produce
# them today is the one that did. See `attr_design_fingerprint()`.
FP      <- attr_design_fingerprint(cfg, fold_k, K_TI, priors, ecfg, design_key = DESIGN)
METHODS <- as.character(cfg_req(ecfg, "methods"))
attr_check_methods(METHODS)

# `deferred_methods` IS ENFORCED, NOT A NOTE. It was documentation only until
# audit B reported it dead on 2026-09-06, and a deferred list nothing reads is a
# list that stops being true the first time someone edits `methods:` without
# scrolling down. `shap_xgb_raw` is on it because its columns carry `n_obs` and
# a missingness channel that map onto the 19 signals only loosely, so its groups
# are not commensurable with the L matrix's -- which is a reason it must not be
# generated, not a reason to remember not to.
.deferred <- as.character(cfg_req(ecfg, "deferred_methods"))
if (length(intersect(METHODS, .deferred))) {
  abort_values(paste0("config/attribution_eval.yml lists method(s) in BOTH ",
                      "`methods` and `deferred_methods`. A deferred method has ",
                      "no defensible column-to-signal map; remove it from one ",
                      "list or the other"),
               intersect(METHODS, .deferred))
}

# ============================================================================
# THE SPEC VOCABULARY, AND THE CHECK THAT COSTS SECONDS
# ============================================================================
#
# HOISTED ABOVE THE RUN DIRECTORY ON 2026-09-06, which is the point rather than
# tidiness. `check_specs()` has to be able to answer `--check-specs` and QUIT
# before `new_run()` is called: the first version of it sat inside the `sample`
# block, so invoking the check created a run directory, walked the ladder route,
# refitted the SHAP spec and seed replicates into it, and had started on the
# forty SHAP bootstrap bags before it reached the thing it was meant to check.
# A verification step that costs ten minutes and leaves a stray run directory
# behind is one nobody runs, which defeats the purpose of having it.
#
# Nothing here reads data or fits anything. It needs `cfg`, `sigs`, `paired` and
# `K_TI`, all of which are resolved above.

.ti_of <- function(sg, scope) interaction_terms(sg, cfg, scope, K_TI)

#' Which formula, and whether this (arm, signal) is fitted at all.
#'
#' Mirrors `layer1_jobs()`'s alias and assign rules exactly. An unpaired
#' signal's `full` IS its `meas` and its `intv` is zero by assignment; a
#' signal with no `trend` covariate has no trend cross term, so its
#' `full_ti_trend` IS its `full`. Getting any of these wrong would fit a model
#' the pipeline never fits and call the result the same arm.
spec_of <- function(arm, sg) {
  is_p <- sg %in% paired
  if (arm == "meas") return(list(model = "meas", ti = character(0), fit = TRUE))
  if (arm == "intv") return(list(model = "intv", ti = character(0), fit = is_p))
  if (arm == "full") return(list(model = if (is_p) "full" else "meas",
                                 ti = character(0), fit = TRUE))
  scope <- sub("^full_ti_", "", arm)
  ti <- if (is_p) .ti_of(sg, scope) else character(0)
  list(model = if (is_p) "full" else "meas", ti = ti, fit = TRUE)
}

#' The formula for one (arm, signal), or NULL where that pair is never fitted.
#'
#' THE `fit` GUARD COMES FIRST AND THAT IS NOT AN OPTIMISATION. `build_formula()`
#' REFUSES an `intv` spec on an unpaired signal, by design and loudly: its
#' formula would be `mortality ~ 1`, so `L_intv` is 0 there by construction and
#' is ASSIGNED rather than fitted (CLAUDE.md, frozen decisions, "Three models
#' per signal"). Seven of the nineteen signals are unpaired, so the very first
#' call `boot_arm("intv", bag)` makes reaches `temperature` and stops the run.
#'
#' THIS FUNCTION IS WHERE THE 2026-09-06 REFACTOR PUT THE DEFECT. `boot_arm()`
#' and `posterior_arm()` each had the loop body inline and each had the guard
#' in the right place:
#'
#'     sp <- spec_of(arm, sg)
#'     if (!sp$fit) next
#'     f  <- build_formula(sg, sp$model, cfg)
#'
#' Hoisting the two duplicated lines into a shared helper is right -- one
#' definition of what an arm's formula is -- but the hoist took the
#' `build_formula()` call with it and LEFT THE GUARD BEHIND at the call site,
#' where it now runs one statement too late. The callers read
#' `af <- arm_formula(...)` then `if (!af$spec$fit) next`, which looks correct
#' and is exactly backwards: the refusal has already fired inside the call
#' whose result the guard was going to inspect.
#'
#' THE CLASS, because it will recur: EXTRACTING A HELPER FROM A GUARDED BLOCK
#' MOVES THE GUARDED CODE AND NOT THE GUARD. A guard and the expression it
#' protects have to travel together, so the extracted function owns both or
#' neither. `--check-specs` below exercises every (arm, signal) pair in
#' seconds and would have caught it before a five-hour run did.
arm_formula <- function(arm, sg) {
  sp <- spec_of(arm, sg)
  if (!sp$fit) return(list(spec = sp, formula = NULL))
  f <- build_formula(sg, sp$model, cfg)
  if (length(sp$ti)) {
    f <- stats::update(f, stats::as.formula(
      paste(". ~ . +", paste(sp$ti, collapse = " + "))))
  }
  list(spec = sp, formula = f)
}

# --- THE CHEAP CHECK THAT WOULD HAVE CAUGHT IT ----------------------------
#
# `--check-specs` walks every (arm, signal) pair, calls `arm_formula()`, and
# reports the fit budget without fitting anything. It runs in seconds.
#
# WHY IT IS NEEDED AT ALL, given that the generator prints a plan first. The
# plan is enumerated from the CONFIG -- method names and replicate counts --
# and never touches `spec_of()` or `build_formula()`, so it cannot see a
# formula that refuses to be built. And the two loops that do are entered only
# under `--levels sample`, which is the five-hour path. Every cheaper
# invocation (`--plan-only`, `--levels spec,seed`) leaves them unexecuted, and
# so did every verification run on 2026-09-06: the posterior arms were fully
# cached, so `posterior_arm()` was skipped, and the LLR bootstrap route did
# not exist yet. A code path whose only exercise is the expensive run is a
# code path that gets debugged by the expensive run.
check_specs <- function() {
  cat("\n=== --check-specs: every (arm, signal), no fit ===\n\n")
  arms <- c("intv", "meas", "full", "full_ti_trend", "full_ti_all")
  rows <- list()
  for (arm in arms) for (sg in sigs) {
    af <- try(arm_formula(arm, sg), silent = TRUE)
    if (inherits(af, "try-error")) {
      cat(sprintf("  *** %-16s %-18s FAILED: %s", arm, sg,
                  conditionMessage(attr(af, "condition"))))
      rows[[length(rows) + 1L]] <- data.frame(arm = arm, signal = sg,
        paired = sg %in% paired, fit = NA, model = NA_character_,
        n_ti = NA_integer_, n_terms = NA_integer_, ok = FALSE,
        stringsAsFactors = FALSE)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(arm = arm, signal = sg,
      paired = sg %in% paired, fit = af$spec$fit, model = af$spec$model,
      n_ti = length(af$spec$ti),
      n_terms = if (is.null(af$formula)) NA_integer_
                else length(attr(stats::terms(af$formula), "term.labels")),
      ok = TRUE, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (!all(out$ok)) {
    print(out[!out$ok, ], row.names = FALSE)
    stop("check-specs: ", sum(!out$ok), " (arm, signal) pair(s) could not be ",
         "resolved. The `sample` routes would stop on the first of them.",
         call. = FALSE)
  }
  ag <- aggregate(list(fitted = out$fit, assigned = !out$fit),
                  by = list(arm = out$arm), FUN = sum)
  ag$spec_folds <- ag$fitted * n_fold
  print(ag, row.names = FALSE)
  cat(sprintf("\n  %d (arm, signal, fold) cells per full pass; %d distinct fits\n",
              sum(ag$spec_folds), PASS_FITS))
  cat(sprintf("  after alias memoisation (%d shared specifications).\n",
              length(SHARED_SPECS)))

  # THE INVARIANT THE POSTERIOR ROUTE DEPENDS ON, asserted rather than assumed.
  # `posterior_arm()` draws with `with_seed(DRAW0 + ctr)` where `ctr` counts
  # (signal, fold) WITHIN an arm, so two arms give a column the SAME draw only
  # while their counters stay aligned. That alignment is what makes an ALIASED
  # column -- an unpaired signal, whose `full_ti_all` IS its `full` -- come out
  # bitwise identical between the two arms rather than differing by independent
  # draw noise, which would put noise into `cond_ti - cond` for a signal where
  # nothing changed.
  #
  # It holds only if every non-`intv` arm fits EVERY signal. Skipping the 7
  # unpaired signals in `full` would look like a saving of 35 fits and would
  # silently desynchronise `ctr` against `meas` and both `full_ti_*` arms. This
  # is the check that catches that, and it costs no fit.
  #
  # THE COMPARISON THAT USED TO BE HERE WAS CIRCULAR ONCE `PASS_FITS` WAS FIXED.
  # It checked the enumeration against `PASS_FITS`, which is now derived from
  # `spec_of()` by the same rule -- so it could only ever agree. It earned its
  # keep exactly once, on 2026-09-06, by catching the 320-against-440
  # disagreement that fixing `PASS_FITS` removed; keeping it afterwards would be
  # a test of arithmetic against itself.
  bad <- out[out$arm != "intv" & !out$fit, , drop = FALSE]
  if (nrow(bad)) {
    print(bad, row.names = FALSE)
    stop("check-specs: ", nrow(bad), " (arm, signal) pair(s) outside `intv` are ",
         "not fitted. The posterior route's common-random-numbers alignment ",
         "requires every non-`intv` arm to fit every signal, or aliased columns ",
         "get independent draws.", call. = FALSE)
  }
  wrong <- out[out$arm == "intv" & out$fit != out$paired, , drop = FALSE]
  if (nrow(wrong)) {
    print(wrong, row.names = FALSE)
    stop("check-specs: `intv` is fitted for an unpaired signal or skipped for a ",
         "paired one. `L_intv` is 0 by construction on the unpaired and is ",
         "assigned, never fitted.", call. = FALSE)
  }
  cat(sprintf("  every (arm, signal) resolves; every non-`intv` arm fits all %d\n",
              length(sigs)))
  cat(sprintf("  signals (the ctr-alignment invariant) and `intv` fits exactly the %d paired.\n",
              length(paired)))
  invisible(out)
}
#' The `cond` counterpart of a base arm, or NA where there is none.
cond_of <- function(arm) if (arm == "meas") NA_character_ else
  paste0("llr_", sub("^full", "cond", arm))

# The fit cost of ONE full pass over every base arm and every fold, DERIVED
# FROM `arm_formula()` rather than restated from the alias rules.
#
# ONE FIT PER DISTINCT (signal, formula, model) AS OF 2026-09-09 (review
# finding A11). `spec_of()` enumerates 88 (arm, signal) cells per fold, but
# only 64 distinct specifications exist among them: an unpaired signal's
# `full`, `full_ti_trend` and `full_ti_all` ARE its `meas`, and a paired
# signal with no `trend` covariate has `full_ti_trend` equal to `full`.
# Until today both fitting loops refitted each alias -- 440 spec-folds per
# pass against the pipeline's 320, and 4,800 redundant `bam()` calls across
# the 40 bootstrap bags -- on the argument that the posterior route's
# common-random-numbers alignment needed the counter `ctr` to advance
# identically in every arm. It does; but advancing a counter does not require
# a fit. The loops now memoise by `.spec_key()` within one bag or one
# posterior pass, the counter advances on a memo hit exactly as on a fit, and
# the hit is REFUSED unless its recorded counter equals the current one, so
# the aliased column is bitwise what the refit would have produced rather
# than assumed to be. `tests/attr_eval_unit.R` asserts the bam determinism
# that identity rests on.
#
# So the budget follows the code: `PASS_FITS` counts the distinct keys the
# loops will actually fit, and `--check-specs` reports the enumerated cells
# beside them.
.spec_key <- function(sg, af) {
  paste(sg, af$spec$model, paste(deparse(af$formula), collapse = ""), sep = "|")
}
.pass_keys <- function() {
  ks <- character(0)
  for (a in c("meas", "intv", "full", "full_ti_trend", "full_ti_all")) {
    for (s in sigs) {
      af <- arm_formula(a, s)
      if (af$spec$fit) ks <- c(ks, .spec_key(s, af))
    }
  }
  ks
}
.PK        <- .pass_keys()
PASS_CELLS <- n_fold * length(.PK)
PASS_FITS  <- n_fold * length(unique(.PK))
# Keys that more than one arm resolves to: the ones the memo must hold.
SHARED_SPECS <- unique(.PK[duplicated(.PK)])
# `--check-specs` runs here rather than directly after `check_specs()` is
# defined, because it reports the distinct-fit budget computed just above.
if (CHECK_SPECS) {
  check_specs()
  quit(save = "no", status = 0L)
}
PLAN <- attr_replicate_plan(ecfg, n_folds = n_fold, n_pass_fits = PASS_FITS)

cat("\n=== the replicate plan, before a row is read ===\n\n")
cat(sprintf("  design key %s   (gamma = %s, %d fitted specs, %d folds, k_ti = %d)\n\n",
            attr_key_hash(DESIGN), DESIGN$bam$gamma, length(DESIGN$formulas),
            n_fold, K_TI))
# `route: none` and the two shared-cost rows are MARKERS, not replicates, and
# are counted as zero replicates so the column means what it says. They are
# still printed: a level that is provably empty for a method is a result, and
# omitting the row would make it look unmeasured.
agg <- aggregate(cbind(replicates = as.integer(!is.na(PLAN$index)), fits = PLAN$fits),
                 by = list(method = PLAN$method, route = PLAN$route),
                 FUN = sum, na.rm = TRUE)
print(agg[order(agg$route, agg$method), ], row.names = FALSE)
cat(sprintf("\n  total model fits implied: %d\n", sum(PLAN$fits, na.rm = TRUE)))
cat(sprintf("  one full LLR pass = %d distinct fits (%d cells); the posterior basis is one\n",
            PASS_FITS, PASS_CELLS))
cat(sprintf("  such pass and the LLR bootstrap is %d of them\n",
            as.integer(cfg_req(ecfg, "levels", "sample", "llr_bootstrap_b"))))
cat(sprintf("  levels requested this invocation: %s\n",
            paste(LEVELS_WANTED, collapse = ", ")))
if (PLAN_ONLY) quit(save = "no", status = 0L)

# --- the run directory, the only path builder here ---------------------------
if (!is.null(RESUME)) {
  if (!dir.exists(RESUME)) stop("--resume: no such directory: ", RESUME, call. = FALSE)
  run <- structure(list(prefix = "attrgen", id = basename(RESUME), path = RESUME,
                        started = Sys.time(), config = cfg,
                        log_file = file.path(RESUME, "log.txt")),
                   class = "llr_run")
  cat(sprintf("\n  resuming %s\n", basename(RESUME)))
} else {
  run <- new_run("attrgen", cfg, note = sprintf(
    "attribution replicates, levels %s, design %s",
    paste(LEVELS_WANTED, collapse = "+"), attr_key_hash(DESIGN)))
}
SUB <- cfg_req(ecfg, "storage", "subdir")
dir.create(file.path(run$path, SUB), showWarnings = FALSE, recursive = TRUE)
save_table(run, PLAN, "replicate_plan", subdir = "diagnostics")

# --- WHAT A RESUME MUST AGREE WITH BEFORE IT WRITES A BYTE -------------------
#
# ADDED 2026-09-09 (review findings A1, A2, A9). A resume checked NOTHING about
# the store it was resuming into. `have()` computed keys from the LIVE design,
# so a resume after a design change would have written a second design's
# files beside the first's and reconciled a manifest naming only the live
# ones; and the key itself omits the booster settings, the covariate
# transformations, the priors, the draw seed and the draw count, so a change
# to any of those would have kept every old replicate. Three things are
# checked here, and the third is stamped rather than asserted when the store
# predates it.
#
#   the design key      must be identical. A differing key is a migration
#                       (`tests/attr_store_migrate.R`) or a regeneration,
#                       never a resume.
#   the training rows   `stay_id` and the measured mask must be identical.
#   the fingerprint     must be identical where the store carries one. Where
#                       it carries none, the ladder replicates on disk are
#                       verified bitwise against what the live design produces
#                       -- free for the LLR arms, one booster refit for SHAP --
#                       and the fingerprint is stamped WITH that evidence at
#                       the close. Until both halves are verified the store
#                       stays unstamped and the consumer refuses it.
DESIGN_P <- file.path(run$path, "design.qs2")
MEAS_P   <- file.path(run$path, "measured.qs2")
OLD_DG   <- if (!is.null(RESUME) && file.exists(DESIGN_P)) qs2::qs_read(DESIGN_P) else NULL
NEED_STAMP <- FALSE
STALE_RULES <- NULL
if (!is.null(OLD_DG)) {
  if (is.null(OLD_DG$design_key)) {
    stop("the store carries no design_key and predates 2026-09-06; it cannot ",
         "be resumed, only regenerated.", call. = FALSE)
  }
  d <- attr_design_diff(OLD_DG$design_key, DESIGN)
  if (nrow(d)) {
    print(d[, c("field", "hash_a", "hash_b")], row.names = FALSE)
    stop("the store in ", basename(run$path), " was built under a different ",
         "design key (", nrow(d), " field(s) differ, above). A resume would ",
         "write a second design's replicates beside the first's. If only ",
         "`formulas` differs by addition, run tests/attr_store_migrate.R; ",
         "otherwise regenerate into a new run.", call. = FALSE)
  }
  if (!identical(as.character(OLD_DG$stay_id), ids_ch)) {
    stop("the store's stay_id vector is not the current training set. That is ",
         "a cohort change and nothing here can resume across it.", call. = FALSE)
  }
  if (file.exists(MEAS_P)) {
    old_meas <- qs2::qs_read(MEAS_P)
    if (!identical(dim(old_meas), dim(meas_ok)) || !identical(dimnames(old_meas), dimnames(meas_ok)) ||
        any(old_meas != meas_ok)) {
      stop("the store's measured mask differs from the live one in ",
           sum(old_meas != meas_ok), " cell(s). The extraction has changed ",
           "under the store; regenerate.", call. = FALSE)
    }
    rm(old_meas)
  }
  if (!is.null(OLD_DG$fingerprint)) {
    d <- attr_fingerprint_diff(OLD_DG$fingerprint, FP)
    if (nrow(d)) {
      cat("\n!!! STORE FINGERPRINT MISMATCH !!!\n\n")
      print(d[, c("field", "hash_a", "hash_b", "routes_affected")], row.names = FALSE)
      # ROUTE-SELECTIVE (2026-09-09). A field that moves every route still
      # stops the resume; a field that moves only some routes quarantines
      # exactly those replicates below (`STALE_RULES`, once the manifest is
      # loaded) and lets the rest of the store stand. Files are moved aside,
      # never deleted, and the fingerprint is re-stamped at the close.
      .st <- attr_stale_routes(d)
      if (.st$all) {
        stop("the store's fingerprint differs from the live design in ", nrow(d),
             " field(s) (above) and at least one of them invalidates EVERY ",
             "route. Regenerate into a new run, or restore the setting.",
             call. = FALSE)
      }
      STALE_RULES <- .st$rules
      cat("\n  the differing field(s) invalidate only the route(s) below; those\n")
      cat("  replicates will be quarantined and regenerated, the rest kept:\n\n")
      print(STALE_RULES, row.names = FALSE)
    } else {
      cat(sprintf("  store fingerprint %s matches the live design.\n", attr_key_hash(FP)))
    }
  } else {
    NEED_STAMP <- TRUE
    cat("  the store carries no fingerprint (written before 2026-09-09). The\n")
    cat("  ladder replicates will be verified against the live design and the\n")
    cat("  fingerprint stamped with that evidence at the close.\n")
  }
} else {
  qs2::qs_save(meas_ok, MEAS_P)
  qs2::qs_save(list(design_key = DESIGN, stay_id = tr, k_ti = K_TI,
                    naming_design_key = DESIGN,
                    fingerprint = FP,
                    fingerprint_evidence = data.frame(
                      method = "(new store)", route = "-", status = "new_store",
                      max_abs_diff = NA_real_, stringsAsFactors = FALSE)),
               DESIGN_P)
}

# --- THE NAMING KEY: WHICH DESIGN THE FILES ON DISK ARE NAMED UNDER ----------
#
# FOUND 2026-09-09 BY THE FIRST RESUME AFTER THE 2026-09-07 MIGRATION, and it
# is a defect in that migration. `attr_replicate_key()` hashes the WHOLE design
# key into every file name. `tests/attr_store_migrate.R` rewrote `design.qs2`
# to the post-spec-change key with bitwise evidence -- correctly -- but left
# the 642 files named under the pre-change key, and bridged the two only
# through `manifest_replicates.csv`, which the consumer reads and the
# generator does not. So a resume computed every key under the LIVE design,
# found no file, reported every ladder replicate as newly written, and
# `reconcile_manifest()` rewrote the manifest down to the eight files it could
# see. Nothing on disk was lost, but a `--levels sample` resume would have
# refitted five hours of replicates that were already there.
#
# The store therefore records TWO keys. `design_key` is the design the store
# describes and is what the resume check and the consumer compare against the
# live design. `naming_design_key` is the design its files are named under,
# and is the one `.coord_key()` hashes. For a store migrated before this field
# existed, the naming key is recovered from the `design_pre*.qs2` the
# migration kept beside it, and must hash to the store's `migrated_from`. A
# file-rename migration to the current key is the cleaner permanent form; it
# touches every file name and is left as an authorised step, not a resume.
NAMING_DESIGN <- DESIGN
if (!is.null(OLD_DG)) {
  if (!is.null(OLD_DG$naming_design_key)) {
    NAMING_DESIGN <- OLD_DG$naming_design_key
  } else if (!is.null(OLD_DG$migrated_from)) {
    pres <- list.files(run$path, pattern = "^design_pre.*[.]qs2$", full.names = TRUE)
    found <- NULL
    for (p in pres) {
      pd <- qs2::qs_read(p)
      if (!is.null(pd$design_key) &&
          identical(attr_key_hash(pd$design_key), OLD_DG$migrated_from)) {
        found <- pd$design_key; break
      }
    }
    if (is.null(found)) {
      stop("the store was migrated from design ", OLD_DG$migrated_from,
           " but no design_pre*.qs2 beside it carries that key, so the files' ",
           "naming key cannot be recovered. Refusing to resume: every existing ",
           "replicate would be invisible and refitted.", call. = FALSE)
    }
    NAMING_DESIGN <- found
  }
}
if (!identical(attr_key_hash(NAMING_DESIGN), attr_key_hash(DESIGN))) {
  cat(sprintf("  files are named under design %s (pre-migration); the store describes %s.\n",
              attr_key_hash(NAMING_DESIGN), attr_key_hash(DESIGN)))
  cat("  Keys are computed under the naming key; the design check above used the live one.\n")
}

# Verification evidence gathered by `emit()` and the ladder checks, written
# beside the fingerprint at the close.
EVID <- list()
.evid <- function(method, route, status, max_abs_diff = NA_real_) {
  EVID[[length(EVID) + 1L]] <<- data.frame(
    method = method, route = route, status = status,
    max_abs_diff = max_abs_diff, stringsAsFactors = FALSE)
}

# ============================================================================
# THE SHARED BOOTSTRAP MANIFEST
# ============================================================================

#' Bootstrap the TRAINING rows, grouped by patient.
#'
#' BY PATIENT AND NOT BY STAY, for the same reason the CV folds are grouped by
#' `patient_id`: a stay-level resample splits a patient's repeat admissions
#' across the in-bag and out-of-bag sets and leaks. Returns a logical over `tr`
#' rather than an index vector, because a duplicated row would break the
#' one-row-per-stay invariant every L matrix depends on -- so this is an
#' m-out-of-n bootstrap WITHOUT replacement at the patient level, which is a
#' subsample and is described as one throughout.
bag_of <- function(seed) {
  # RESOLVED, NEVER NAMED. `resample_cols()` is the one place the grouping key
  # is decided (R/03_folds.R), so the bootstrap groups by whatever
  # `fold_group_by` resolves to rather than by a literal that would drift from
  # it -- the same reason `check_resampling_cols()` exists.
  gcol <- resample_cols(tabs$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  if (anyNA(gid)) stop("bag_of: `", gcol, "` is missing for some training stay",
                       call. = FALSE)
  ug   <- unique(gid)
  with_seed(seed, {
    keep <- sample(ug, length(ug), replace = TRUE)
    gid %in% unique(keep)
  })
}

BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
# Hoisted above the manifest section on 2026-09-09: `.coord_key()` reads
# `SEED_B` for the seed routes, and the ladder verification below now calls
# it before the SHAP block that used to define these.
SEED_B    <- as.integer(cfg_req(ecfg, "levels", "seed", "seed_base"))
B_PER_BAG <- as.integer(cfg_req(ecfg, "levels", "seed", "b_per_bag"))
B_SEED    <- as.integer(cfg_req(ecfg, "levels", "seed", "b"))

#' Bag `b` of the shared manifest. THE ONE PLACE A BAG IS DEFINED.
#'
#' Memoised because every method asks for the same forty bags and `bag_of()` is
#' a sort and a match over 41,250 stays. The memo is what makes "the same bag"
#' a fact about this process rather than a claim about two seeds agreeing.
.BAGS <- new.env(parent = emptyenv())
bag_for <- function(b) {
  b <- as.integer(b)
  if (b < 1L || b > B_BOOT) {
    abort_values(paste0("bag_for: boot_id outside the declared shared manifest ",
                        "(1..", B_BOOT, "). Every method draws from the same ",
                        "bags, so a bag outside the manifest is not a bag any ",
                        "other method could have used"), b)
  }
  k <- paste0("b", b)
  if (is.null(.BAGS[[k]])) .BAGS[[k]] <- bag_of(BAG_SEED0 + b)
  .BAGS[[k]]
}

# The bag table: COUNTS AND A HASH, never membership (hard rule 1). The hash is
# what lets a later run prove it used the same bags without either run printing
# who was in them.
BAGT <- do.call(rbind, lapply(seq_len(B_BOOT), function(b) {
  g <- bag_for(b)
  data.frame(boot_id = b, seed = BAG_SEED0 + b, n_in_bag = sum(g),
             n_train = length(g), frac_in_bag = round(mean(g), 5),
             membership_hash = attr_key_hash(which(g)), stringsAsFactors = FALSE)
}))
save_table(run, BAGT, "bootstrap_bags", subdir = "diagnostics")
cat(sprintf("\n  shared bootstrap manifest: %d bags, in-bag fraction %.4f to %.4f\n",
            B_BOOT, min(BAGT$frac_in_bag), max(BAGT$frac_in_bag)))
cat(sprintf("  manifest hash %s  (every method draws from these bags)\n",
            attr_key_hash(BAGT$membership_hash)))

# --- THE TOMBSTONE, hoisted here on 2026-09-09 so the close-of-run coverage
# table can read it whatever `--levels` was asked for --------------------------
#
# CONTENT ADDRESSING MAKES A SUCCESS SKIPPABLE AND SAYS NOTHING ABOUT A
# FAILURE. `have()` asks whether a replicate file exists, so an excluded bag
# -- which by construction writes nothing -- is indistinguishable from a bag
# that was never attempted. The first resume after the exclusions refitted
# bags 20 and 26 in full, about 14 minutes, to exclude them again for the
# same reason. So an exclusion is RECORDED rather than merely reported.
#
# KEYED BY DESIGN AND BY BAG IDENTITY (review finding A9). The record carried
# the design key and the bag INDEX, so a change to `bootstrap.seed_base` --
# which moves every bag while leaving the design key untouched -- would have
# carried an exclusion onto a bag that was never tried. `bag_seed` and the
# bag's `membership_hash` now travel with it and both must match the live
# bag. A record written before these columns existed is backfilled from the
# live bags, which is correct exactly because `seed_base` has not moved since
# 2026-09-06 (the SHAP replicate keys, which DO hash the bag seed, would have
# been orphaned otherwise), and the backfill is announced.
#
# CARRIED ACROSS A PROVEN MIGRATION. `tests/attr_store_migrate.R` rewrote the
# design key on 2026-09-07 with bitwise evidence and left the tombstone under
# the OLD key, so the filter below silently dropped both records and the next
# `--levels sample` resume would have spent 14 minutes re-excluding bags 20
# and 26. A record whose design is the store's `migrated_from` key is the
# same exclusion under the migrated design and is carried forward.
#
# `--retry-excluded` clears this design's records, for the case where the
# exclusion was caused by something since fixed.
TOMB_P    <- file.path(run$path, "llr_bootstrap_excluded.csv")
TOMB_COLS <- c("boot_id", "design", "bag_seed", "membership_hash",
               "n_spec_folds", "cause")
.empty_tomb <- function() data.frame(
  boot_id = integer(0), design = character(0), bag_seed = integer(0),
  membership_hash = character(0), n_spec_folds = integer(0),
  cause = character(0), stringsAsFactors = FALSE)
.DK <- attr_key_hash(DESIGN)
TOMB_ALL <- if (file.exists(TOMB_P)) utils::read.csv(TOMB_P, stringsAsFactors = FALSE) else .empty_tomb()
if (nrow(TOMB_ALL)) {
  if (!all(c("bag_seed", "membership_hash") %in% names(TOMB_ALL))) {
    file.copy(TOMB_P, sub("[.]csv$", "_pre20260909.csv", TOMB_P), overwrite = FALSE)
    TOMB_ALL$bag_seed <- BAG_SEED0 + as.integer(TOMB_ALL$boot_id)
    TOMB_ALL$membership_hash <- BAGT$membership_hash[match(TOMB_ALL$boot_id, BAGT$boot_id)]
    cat(sprintf("  %d legacy tombstone row(s) backfilled with the live bag identity\n",
                nrow(TOMB_ALL)))
    cat("  (bootstrap.seed_base unchanged since 2026-09-06; original kept as *_pre20260909.csv)\n")
  }
  mig_from <- if (!is.null(OLD_DG) && !is.null(OLD_DG$migrated_from)) OLD_DG$migrated_from else NA_character_
  carry <- !is.na(mig_from) & TOMB_ALL$design == mig_from
  if (any(carry)) {
    TOMB_ALL$design[carry] <- .DK
    cat(sprintf("  %d tombstone row(s) carried from migrated design %s to %s\n",
                sum(carry), mig_from, .DK))
  }
  TOMB_ALL <- TOMB_ALL[, TOMB_COLS, drop = FALSE]
  utils::write.csv(TOMB_ALL, TOMB_P, row.names = FALSE)
}
cur_hash <- BAGT$membership_hash[match(TOMB_ALL$boot_id, BAGT$boot_id)]
live <- TOMB_ALL$design == .DK & !is.na(cur_hash) &
  TOMB_ALL$membership_hash == cur_hash &
  TOMB_ALL$bag_seed == BAG_SEED0 + as.integer(TOMB_ALL$boot_id)
n_other <- sum(TOMB_ALL$design == .DK & !live)
if (n_other) {
  cat(sprintf("  %d tombstone row(s) for this design refer to a DIFFERENT bag and are ignored\n",
              n_other))
}
TOMB <- if (RETRY_EXCL) .empty_tomb() else TOMB_ALL[live, , drop = FALSE]
if (RETRY_EXCL && any(live)) {
  TOMB_ALL <- TOMB_ALL[!live, , drop = FALSE]
  utils::write.csv(TOMB_ALL, TOMB_P, row.names = FALSE)
  cat("  --retry-excluded: this design's tombstones cleared\n")
}
if (nrow(TOMB)) {
  cat(sprintf("  %d bag(s) tombstoned as unsupported under design %s: %s\n",
              nrow(TOMB), .DK, paste(sort(TOMB$boot_id), collapse = ", ")))
  cat("  (pass --retry-excluded to attempt them again)\n")
}
tombstone_bag <- function(b, n_spec_folds, cause) {
  row <- data.frame(boot_id = b, design = .DK, bag_seed = BAG_SEED0 + b,
                    membership_hash = BAGT$membership_hash[match(b, BAGT$boot_id)],
                    n_spec_folds = n_spec_folds, cause = cause,
                    stringsAsFactors = FALSE)
  TOMB     <<- rbind(TOMB, row)
  TOMB_ALL <<- rbind(TOMB_ALL, row)
  # Written IMMEDIATELY, not at the close, so an interrupted run does not
  # lose the knowledge and repeat the work.
  utils::write.csv(TOMB_ALL, TOMB_P, row.names = FALSE)
}

# ============================================================================
# THE MANIFEST, AND THE UPSERT
# ============================================================================

MAN_COLS <- c("method", "route", "stage", "index", "boot_id", "seed_id",
              "draw_id", "key", "n_rows", "n_cols")
MAN_P <- file.path(run$path, "manifest_replicates.csv")
MAN <- data.frame(method = character(0), route = character(0), stage = character(0),
                  index = integer(0), boot_id = integer(0), seed_id = integer(0),
                  draw_id = integer(0), key = character(0), n_rows = integer(0),
                  n_cols = integer(0), stringsAsFactors = FALSE)
if (file.exists(MAN_P)) {
  .old <- utils::read.csv(MAN_P, stringsAsFactors = FALSE)
  if (all(MAN_COLS %in% names(.old))) {
    MAN <- .old[, MAN_COLS, drop = FALSE]
  } else {
    # A manifest from before 2026-09-06 has `level` where `stage` now is and
    # carries no coordinates. It is DISCARDED rather than migrated, and that is
    # not laziness: an old `index` is a bag id, a seed id or a posterior draw
    # id depending on which loop wrote it, which is exactly the ambiguity the
    # coordinates exist to remove, and the old LLR bootstrap rows are on bags
    # that no longer exist. `reconcile_manifest()` below rebuilds every row
    # from the plan and from what is on disk, so nothing is lost that is still
    # valid, and the old file is kept beside it.
    file.rename(MAN_P, file.path(run$path, "manifest_replicates_pre20260906.csv"))
    cat("\n  the manifest on disk predates the replicate coordinates. It has been\n")
    cat("  renamed to manifest_replicates_pre20260906.csv and will be rebuilt\n")
    cat("  from the plan. NO REPLICATE FILE IS TOUCHED -- a key that exists\n")
    cat("  cannot hold different content.\n")
  }
}

#' THE ONE PLACE A COORDINATE BECOMES A STORAGE KEY.
#'
#' `stage` is the storage stage ("spec", "seed", "sample") that
#' `attr_replicate_key()` validates, `index` keeps two replicates of one method
#' at one stage distinguishable, and `extra` is what actually distinguishes
#' them in the content hash. Everything downstream reads the COORDINATES; this
#' mapping exists so the key can stay byte-compatible with replicates already
#' on disk while the coordinates are new.
#'
#' WHY IT IS A FUNCTION RATHER THAN SIX INLINE `emit()` ARGUMENTS. The mapping
#' was written out at each of the six call sites in the first draft, and
#' `reconcile_manifest()` below needs it a seventh time. Seven copies of an
#' encoding is the F9-to-F11 pattern with the copies inside one file: a
#' `bootstrap_seed` that drifts at one site produces a key nothing can find,
#' and the symptom is a replicate silently refitted rather than an error.
#'
#' THE INDEX OFFSETS. `bootstrap` for an LLR arm is 500 + b so it cannot
#' collide with a posterior draw index; `bootstrap_seeded` is 1000*s + b so a
#' bag's re-seeded fit cannot collide with its seed-0 fit or with another bag's.
#' `attr_replicate_plan()` asserts `bootstrap.b <= 400`, which is what makes
#' both offsets safe.
.coord_key <- function(method, route, boot_id = 0L, seed_id = 0L, draw_id = 0L) {
  fam <- unname(attr_method_family(method))
  b <- as.integer(boot_id); sd <- as.integer(seed_id); dw <- as.integer(draw_id)
  z <- switch(route,
    ladder = if (fam == "shap") list("spec", 1L, list(seed_offset = 0L))
             else list("spec", 1L, NULL),
    seed = list("seed", sd, list(seed_offset = SEED_B * sd)),
    bootstrap = if (fam == "shap")
        list("sample", b, list(bootstrap_seed = BAG_SEED0 + b))
      else
        list("sample", 500L + b,
             list(route = "bootstrap", bootstrap_seed = BAG_SEED0 + b)),
    bootstrap_seeded = list("sample", 1000L * sd + b,
        list(bootstrap_seed = BAG_SEED0 + b, seed_offset = SEED_B * sd)),
    posterior = list("sample", dw, list(route = "posterior")),
    abort_values(".coord_key: no key encoding for route", route))
  # `NAMING_DESIGN`, not `DESIGN`: the key the files on disk were named under.
  # The two differ only for a store migrated across the 2026-09-07 spec change.
  list(stage = z[[1]], index = z[[2]], extra = z[[3]],
       key = attr_replicate_key(method, z[[1]], z[[2]], NAMING_DESIGN, z[[3]]))
}

.path_of <- function(k) file.path(run$path, SUB, paste0(k$key, ".qs2"))

# --- QUARANTINE OF FINGERPRINT-STALE REPLICATES ------------------------------
# Runs once the manifest and the key functions exist. Every replicate the
# stale rules name is MOVED into a subdirectory named by the old fingerprint,
# with the SHAP ladder margin file when that ladder is stale, and dropped from
# the manifest, so `have()` is FALSE and the route regenerates it under the
# live design. Nothing is deleted. File names are hashes (hard rule 1).
STALE_MOVED <- 0L
if (!is.null(STALE_RULES)) {
  stale <- attr_stale_mask(MAN, STALE_RULES)
  qdir  <- file.path(run$path, SUB,
                     paste0("stale_", substr(attr_key_hash(OLD_DG$fingerprint), 1, 12)))
  mv <- function(from) {
    if (!file.exists(from)) return(invisible(FALSE))
    dir.create(qdir, showWarnings = FALSE, recursive = TRUE)
    ok <- file.rename(from, file.path(qdir, basename(from)))
    if (!isTRUE(ok)) stop("quarantine: could not move ", basename(from), call. = FALSE)
    STALE_MOVED <<- STALE_MOVED + 1L
    invisible(TRUE)
  }
  for (i in which(stale)) mv(file.path(run$path, SUB, paste0(MAN$key[i], ".qs2")))
  if (any(stale & MAN$route == "ladder" & attr_method_family(MAN$method) == "shap")) {
    mv(file.path(run$path, SUB, paste0(.coord_key("shap_xgb_feat", "ladder")$key,
                                       "__margin.qs2")))
  }
  if (any(stale)) {
    print(table(method = MAN$method[stale], route = MAN$route[stale]))
  }
  MAN <- MAN[!stale, , drop = FALSE]
  utils::write.csv(MAN, MAN_P, row.names = FALSE)
  cat(sprintf("  quarantined %d file(s) into %s; %d replicate(s) kept.\n",
              STALE_MOVED, basename(qdir), nrow(MAN)))
}

#' Write the matrix if it is not there, VERIFY it if it is, and upsert its
#' manifest row either way.
#'
#' VERIFIED, NOT SKIPPED, AS OF 2026-09-09 (review finding A1). "A key that
#' already exists cannot be produced with different content" was true only
#' of the fields the key hashes. When a caller has the matrix in hand anyway
#' -- the ladder route always does -- an existing file is compared bitwise
#' and a difference stops the run, because it means the store no longer
#' describes the design whatever the key says. Routes that fit only when
#' `have()` is FALSE never reach the comparison and cost nothing extra.
#'
#' @return "written" or "verified", invisibly.
emit <- function(M, method, route, boot_id = 0L, seed_id = 0L, draw_id = 0L) {
  k <- .coord_key(method, route, boot_id, seed_id, draw_id)
  status <- "written"
  if (file.exists(.path_of(k))) {
    old <- qs2::qs_read(.path_of(k))
    dm <- if (identical(dim(old), dim(M)) && identical(dimnames(old), dimnames(M)))
      max(abs(old - M)) else Inf
    if (!is.finite(dm) || dm > 0) {
      stop("emit: replicate `", k$key, "` (", method, ", ", route, ", boot ",
           boot_id, ", seed ", seed_id, ", draw ", draw_id, ") already exists ",
           "on disk with DIFFERENT content: max |stored - live| = ", format(dm),
           ". The design key is unchanged, so the store was built under a ",
           "setting the key does not cover. Do not resume into it.",
           call. = FALSE)
    }
    status <- "verified"
    .evid(method, route, "verified_bitwise", dm)
  } else {
    qs2::qs_save(M, .path_of(k))
    .evid(method, route, "written")
  }
  row <- data.frame(method = method, route = route, stage = k$stage,
                    index = as.integer(k$index), boot_id = as.integer(boot_id),
                    seed_id = as.integer(seed_id), draw_id = as.integer(draw_id),
                    key = k$key, n_rows = nrow(M), n_cols = ncol(M),
                    stringsAsFactors = FALSE)
  MAN <<- rbind(MAN[MAN$key != k$key, , drop = FALSE], row)
  utils::write.csv(MAN, MAN_P, row.names = FALSE)
  invisible(status)
}

have <- function(method, route, boot_id = 0L, seed_id = 0L, draw_id = 0L) {
  file.exists(.path_of(.coord_key(method, route, boot_id, seed_id, draw_id)))
}

#' Rebuild the manifest from the PLAN and from what is on disk.
#'
#' THE MANIFEST IS DERIVED, NOT APPENDED, AND THIS IS THE FIX FOR A REAL
#' DEFECT IN THE FIRST DRAFT. There, a resumed run rewrote the manifest from
#' the replicates it VISITED, so `--levels spec,seed` on a store that also held
#' 320 `sample` replicates would have written a manifest describing 9 of them.
#' Nothing would have errored: the consumer would have read a complete-looking
#' manifest, found no level-3 pairs, and reported that the sampling contrast had
#' no replicates -- which is indistinguishable from the truthful version of that
#' message. The class is worth naming because it is the same one as audit
#' finding F3's gap (section 24 of the plan): A DERIVED INDEX THAT IS UPDATED
#' INCREMENTALLY REPORTS THE LAST UPDATE, NOT THE STATE.
#'
#' So the manifest is a function of the plan and the filesystem, computed at
#' the end of every invocation whatever `--levels` was asked for. A replicate
#' whose file is present gets a row; one whose file is absent does not.
reconcile_manifest <- function(plan) {
  rows <- list()
  pr <- plan[!is.na(plan$index) & plan$route %in%
               c("ladder", "seed", "bootstrap", "bootstrap_seeded", "posterior"), ]
  for (i in seq_len(nrow(pr))) {
    k <- .coord_key(pr$method[i], pr$route[i], pr$boot_id[i], pr$seed_id[i],
                    pr$draw_id[i])
    if (!file.exists(.path_of(k))) next
    old <- MAN[MAN$key == k$key, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      method = pr$method[i], route = pr$route[i], stage = k$stage,
      index = as.integer(k$index), boot_id = as.integer(pr$boot_id[i]),
      seed_id = as.integer(pr$seed_id[i]), draw_id = as.integer(pr$draw_id[i]),
      key = k$key,
      n_rows = if (nrow(old)) old$n_rows[1] else NA_integer_,
      n_cols = if (nrow(old)) old$n_cols[1] else NA_integer_,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(MAN[0, , drop = FALSE])
  out <- do.call(rbind, rows)
  # A replicate whose row was rebuilt from the plan has no stored shape. One
  # read per such row, and only for rows the plan did not emit this invocation.
  need <- which(is.na(out$n_rows))
  for (i in need) {
    M <- qs2::qs_read(file.path(run$path, SUB, paste0(out$key[i], ".qs2")))
    out$n_rows[i] <- nrow(M); out$n_cols[i] <- ncol(M)
  }
  out
}

# ============================================================================
# ROUTE `ladder` -- stage `spec`, coordinate (0, 0, 0)
# ============================================================================


#' The five base out-of-fold L matrices, unperturbed. ALL FIVE NOW COME FROM
#' THE TARGETS CACHE.
#'
#' UNTIL 2026-09-07 THIS FUNCTION REFITTED THE TWO INTERACTION ARMS, or read
#' them out of a `coupattr` run named in `attribution_eval.ladder.run`, because
#' the pipeline did not fit them: `LAYER1_MODELS` was `meas`, `full`, `intv`,
#' and the interaction models existed only inside analysis scripts. Now they are
#' pipeline specs, so `l_mats_zero` carries `full_ti_trend` and `full_ti_all`
#' beside the other three and there is nothing left to refit. Six minutes saved
#' is the small part; the real change is that the attribution arm's ANCHOR is
#' now literally the pipeline's own out-of-fold L rather than a second fit of
#' the same formula, so "the level-4 contrast is between two arms the paper
#' reports" stops being a claim and becomes an identity.
#'
#' `ladder.run` IS NOW A CROSS-CHECK RATHER THAN A SOURCE. If it names a run,
#' every arm in it is compared against the targets store and any difference
#' stops the generator. That is the evidence the spec change is meant to
#' produce: the 371 replicates on disk were built from the cached ladder, so a
#' bitwise match between that ladder and the pipeline's new fits is what says
#' the stored replicates still describe the design. A refit that reproduces its
#' reference is worth more as an assertion than as a saving.
base_arms <- function() {
  need <- c("meas", "intv", "full", "full_ti_trend", "full_ti_all")
  miss <- setdiff(need, names(Lz))
  if (length(miss)) {
    abort_values(paste0("base_arms: `l_mats_zero` is missing L matri(ces). The ",
                        "interaction models must be in LAYER1_MODELS and the ",
                        "graph must be up to date -- run targets::tar_make()"),
                 miss)
  }
  B <- Lz[need]

  if (!is.null(LADDER_D)) {
    if (!dir.exists(LADDER_D)) {
      stop("attribution_eval.ladder.run names a directory that does not exist: ",
           LADDER_D, "\nSet the key to \"\" to skip the cross-check.", call. = FALSE)
    }
    lad <- readRDS(file.path(LADDER_D, "tables", "l_oof_ladder.rds"))
    if (is.null(lad$design_key)) {
      stop("the ladder in ", basename(LADDER_D), " carries no design_key, so it ",
           "predates 2026-09-06 and cannot be checked. Re-run ",
           "tests/coupling_attribution.R.", call. = FALSE)
    }
    # KEYED TO THE LADDER'S OWN `k_ti`, not to today's. The stored object is
    # being checked against the pipeline, and the question is whether the L
    # values agree -- keying it to the current config would refuse the
    # comparison on a field the comparison is supposed to test.
    d <- attr_design_diff(lad$design_key,
                          attr_design_key(cfg, fold_k, lad$k_ti %||% K_TI))
    # `formulas` differs BY CONSTRUCTION after the spec change: the pipeline now
    # enumerates 64 specs where the ladder was built under 43. Every other field
    # must still match, and the L comparison below is what decides whether the
    # differing field mattered.
    d <- d[d$field != "formulas", , drop = FALSE]
    if (nrow(d)) {
      print(d[, c("field", "hash_a", "hash_b")], row.names = FALSE)
      stop("the ladder in ", basename(LADDER_D), " was built under a different ",
           "design (", nrow(d), " field(s) differ besides `formulas`, above).",
           call. = FALSE)
    }
    stopifnot(identical(as.character(lad$stay_id), ids_ch))
    cat(sprintf("  cross-checking every base arm against %s\n", basename(LADDER_D)))
    bad <- character(0)
    for (a in need) {
      if (is.null(lad$arms[[a]])) next
      dm <- max(abs(B[[a]] - lad$arms[[a]]))
      cat(sprintf("    %-16s max |pipeline - ladder| = %.3e\n", a, dm))
      if (dm > 0) bad <- c(bad, a)
    }
    if (length(bad)) {
      abort_values(paste0("base_arms: the pipeline's L differs from the cached ",
                          "ladder for arm(s) below. The 371 stored replicates ",
                          "were built from the ladder, so a difference here ",
                          "means the store no longer describes the design and ",
                          "the migration must not proceed"), bad)
    }
  }
  B
}

#' The LLR ladder: written where absent, VERIFIED bitwise where present.
#'
#' `write = FALSE` is the close-of-run form used when `spec` was not among
#' the requested levels but the store still needs its fingerprint stamped:
#' every ladder replicate on disk is compared against the targets cache and
#' nothing is written. The comparison is free -- `l_mats_zero` is already in
#' memory -- which is why it is not optional.
llr_ladder <- function(write) {
  A <- attr_derive_arms(base_arms())
  for (m in intersect(METHODS, paste0("llr_", ATTR_LLR_ARMS))) {
    L <- A[[attr_arm_of(m)]][, sigs, drop = FALSE]
    if (write || have(m, "ladder")) {
      st <- emit(L, m, "ladder")
      cat(sprintf("  %-22s spec  %s\n", m, st))
    } else {
      .evid(m, "ladder", "absent")
      cat(sprintf("  %-22s spec  absent (not requested)\n", m))
    }
  }
  invisible(TRUE)
}

if ("spec" %in% LEVELS_WANTED) {
  cat("\n=== ROUTE `ladder`: stage `spec`, coordinate (0,0,0) ===\n\n")
  t0 <- start_timer()
  llr_ladder(write = TRUE)
  cat(sprintf("\n  ladder route (LLR arms) complete in %.1f minutes.\n",
              t0()$elapsed_sec / 60))
}

# ============================================================================
# THE SHAP ROUTES -- `ladder`, `seed`, `bootstrap`, `bootstrap_seeded`
# ============================================================================

DESIGN_CACHE <- new.env(parent = emptyenv())
#' The per-fold `xgb_feat` design, built once and reused across every replicate.
#'
#' The design depends on the FOLD, because `pi_hat`, `delta` and `lambda` are
#' evaluated from out-of-fold priors, and on nothing else -- not on the seed and
#' not on which rows a replicate has in bag. Rebuilding it per replicate was
#' costing about as much as the boosters it feeds. 5 folds x 99 columns x 41,250
#' rows is roughly 163 MB held.
design_for <- function(f) {
  k <- paste0("f", f)
  if (is.null(DESIGN_CACHE[[k]])) {
    DESIGN_CACHE[[k]] <- xgb_design_feat(tabs, cfg, priors, tr, role = "oof", fold = f)
  }
  DESIGN_CACHE[[k]]
}

#' One out-of-fold SHAP matrix over the 19 signal groups.
#'
#' `xgb_feat`'s 99 columns ARE 31 groups -- the 19 signals and the 12
#' interventions -- but this function keeps only the 19 signal groups. THE 29
#' INTERVENTION COLUMNS ARE DROPPED, NOT ROLLED UP: `grp <- sub("__.*$", "",
#' ref)` labels every column by its group and `for (g in sigs)` iterates the 19
#' signal names only, so an intervention group is never selected and never
#' summed. Measured 2026-09-06 on one fold: 10.1% of |SHAP| mass sits in the
#' dropped columns, and reconstructing the model's log-odds from the kept 19
#' misses by a median of 0.293 nats -- TreeSHAP's additivity, the property that
#' makes summing within a group legitimate at all, holds for the true 31-group
#' partition and NOT for this one.
#'
#' THIS IS A DELIBERATE TRUNCATION WITH A DEFENSIBLE READING, KEPT BECAUSE
#' TESTING IT AGAINST THE ALTERNATIVE CHANGED NOTHING THAT MATTERS. What
#' remains after dropping the intervention columns is the measurement
#' contribution CONDITIONAL ON the intervention features already being in the
#' model, which makes the returned matrix structurally the `llr_cond` analogue
#' rather than `llr_full`'s. Two other rollups were measured against it on 12
#' shared bags, same boosters: redistributing each intervention's mass onto its
#' cognate signal(s) (additive, same column count) moved the level-3 top-1
#' disagreement by at most 0.015 against a SHAP-vs-LLR gap of 0.085 to 0.13, so
#' the headline stability comparison does not depend on this choice. It DOES
#' change which signal is named as a patient's top contributor, at an 8.8%
#' rate, so any per-patient claim of the form "SHAP named signal g" carries that
#' as a caveat the noise-calibrated tie tolerance does not cover -- it is not
#' noise, it is this construction choice. Full derivation in
#' `docs/attribution_analysis_plan_20260906.md` section 41.
#'
#' @param seed_offset added to `cfg$seed + fold`. 0 is `seed_id = 0`.
#' @param bag optional logical over `tr`: which training stays are in bag. THE
#'   EVALUATION COHORT NEVER MOVES -- SHAP is computed on the entire held-out
#'   fold either way -- because resampling it too would confound sampling of the
#'   evaluation set with sampling of the estimator, and the estimator is the
#'   only thing under study.
#' @return a list: `G`, the 19-group matrix that is the replicate; `eta`, the
#'   out-of-fold logit prediction of the same boosters (review finding A4);
#'   `best_iter` per fold; `additivity_max_abs`, the largest gap between the
#'   full 100-column TreeSHAP sum and `eta`, which is the check that the
#'   contribution path and the prediction path describe the same trees.
shap_oof <- function(seed_offset = 0L, bag = NULL) {
  S <- ref <- NULL
  eta <- rep(NA_real_, length(tr)); best <- integer(0); addit <- numeric(0)
  eps <- 1e-6
  for (f in sort(unique(fold_k))) {
    X <- design_for(f)
    if (is.null(ref)) {
      ref <- colnames(X)
      S <- matrix(NA_real_, length(tr), length(ref),
                  dimnames = list(ids_ch, ref))
    } else if (!identical(colnames(X), ref)) {
      stop("fold ", f, " produced a different column set", call. = FALSE)
    }
    ho  <- fold_k == f
    use <- !ho & (if (is.null(bag)) TRUE else bag)
    b <- .xgb_fit1(X[use, , drop = FALSE], y[use], cfg, seed = cfg$seed + f + seed_offset,
                   group = grp_k[use])$booster
    ctr <- stats::predict(b, xgboost::xgb.DMatrix(X[ho, , drop = FALSE], missing = NA),
                          predcontrib = TRUE)
    stopifnot(ncol(ctr) == length(ref) + 1L)   # the trailing BIAS column
    S[ho, ] <- ctr[, seq_along(ref), drop = FALSE]
    ph <- .xgb_predict(b, X[ho, , drop = FALSE])
    eta[ho] <- logit(pmin(pmax(ph, eps), 1 - eps))
    best  <- c(best, .xgb_best_iter(b))
    addit <- c(addit, max(abs(rowSums(ctr) - eta[ho])))
  }
  stopifnot(!anyNA(S), !anyNA(eta))
  grp <- sub("__.*$", "", ref)
  G <- matrix(0, nrow(S), length(sigs), dimnames = list(ids_ch, sigs))
  for (g in sigs) {
    j <- which(grp == g)
    if (length(j)) G[, g] <- rowSums(S[, j, drop = FALSE])
  }
  list(G = G, eta = stats::setNames(eta, ids_ch), p_bar = mean(y),
       best_iter = best, additivity_max_abs = max(addit))
}

#' Where the SHAP ladder's out-of-fold MARGIN lives, beside its replicate.
#'
#' ADDED 2026-09-09 (review finding A4). The stored SHAP matrix is 19 signal
#' groups with the BIAS and the 12 intervention groups dropped, so its row
#' sums are NOT the booster's prediction and an AUROC of them is not the
#' model's discrimination. The consumer's discrimination axis reads this file
#' instead: `eta` is `logit(p_hat)` from the same five out-of-fold boosters
#' that produced the ladder replicate, one value per training stay, clamped
#' exactly as `xgb_oof_perfold()` clamps. Row-level, saved like every
#' replicate, never printed.
.shap_margin_path <- function() {
  file.path(run$path, SUB, paste0(.coord_key("shap_xgb_feat", "ladder")$key,
                                  "__margin.qs2"))
}

if ("shap_xgb_feat" %in% METHODS) {
  if ("spec" %in% LEVELS_WANTED) {
    mp <- .shap_margin_path()
    if (have("shap_xgb_feat", "ladder") && file.exists(mp) && !VERIFY_SHAP) {
      cat("  shap_xgb_feat spec  cached (margin present)\n")
    } else {
      t0 <- start_timer()
      # A refit when the replicate exists is the SHAP half of the fingerprint
      # evidence: `emit()` compares the group matrix bitwise and stops on a
      # difference (finding A1). It costs one booster per fold.
      r  <- shap_oof(0L)
      st <- emit(r$G, "shap_xgb_feat", "ladder")
      qs2::qs_save(list(eta = r$eta, p_bar = r$p_bar, best_iter = r$best_iter,
                        additivity_max_abs = r$additivity_max_abs,
                        replicate_key = .coord_key("shap_xgb_feat", "ladder")$key,
                        clamp_eps = 1e-6), mp)
      cat(sprintf("  shap_xgb_feat spec %s in %.1f min; margin saved (TreeSHAP additivity max |sum - eta| = %.2e)\n",
                  st, t0()$elapsed_sec / 60, r$additivity_max_abs))
    }
  }
  if ("seed" %in% LEVELS_WANTED && B_SEED > 0L) {
    cat("\n=== ROUTE `seed`: contrast L2 at the original sample, SHAP only ===\n\n")
    cat("    Nothing changes but the random draw. This is estimation noise in a\n")
    cat("    method with NO SPECIFICATION CHOICE AT ALL, and it is the reference\n")
    cat("    point every level-4 number has to be read against.\n\n")
    for (s in seq_len(B_SEED)) {
      if (have("shap_xgb_feat", "seed", seed_id = s)) {
        cat(sprintf("  seed %2d  cached\n", s)); next
      }
      t0 <- start_timer()
      emit(shap_oof(SEED_B * s)$G, "shap_xgb_feat", "seed", seed_id = s)
      cat(sprintf("  seed %2d  offset %5d  %.1f min\n", s, SEED_B * s,
                  t0()$elapsed_sec / 60))
    }
  }
}

# ============================================================================
# STAGE `sample` -- the shared bags
# ============================================================================

if ("sample" %in% LEVELS_WANTED) {
  B_DRAW  <- as.integer(cfg_req(ecfg, "levels", "sample", "b"))
  ROUTE   <- as.character(cfg_req(ecfg, "levels", "sample", "llr_route"))
  B_LBOOT <- as.integer(cfg_req(ecfg, "levels", "sample", "llr_bootstrap_b"))
  DRAW0   <- as.integer(cfg_req(ecfg, "levels", "sample", "draw_seed_base"))

  # --- SHAP: the shared bags, at seed 0 and at each extra seed ---------------
  if ("shap_xgb_feat" %in% METHODS) {
    cat(sprintf("\n=== ROUTE `bootstrap`: SHAP on the %d shared bags ===\n\n", B_BOOT))
    cat(sprintf("    Plus %d re-seeded fit(s) per bag, which is the whole of the\n",
                B_PER_BAG))
    cat("    level-2-into-level-3 propagation: with two seeds on one bag, seed\n")
    cat("    noise and sampling noise can be varied together (L3T) as well as\n")
    cat("    separately (L2, L3), without assuming a variance decomposition.\n\n")
    for (b in seq_len(B_BOOT)) {
      # seed_id 0. `.coord_key()` reproduces the PRE-2026-09-06 `extra` for this
      # route exactly, so the 40 replicates already on disk keep their keys and
      # are skipped rather than refitted to reproduce themselves bitwise.
      if (have("shap_xgb_feat", "bootstrap", boot_id = b)) {
        cat(sprintf("  bag %2d  seed 0  cached\n", b))
      } else {
        t0 <- start_timer()
        emit(shap_oof(0L, bag = bag_for(b))$G, "shap_xgb_feat", "bootstrap",
             boot_id = b)
        cat(sprintf("  bag %2d  seed 0  %.1f min\n", b, t0()$elapsed_sec / 60))
      }
      for (s in seq_len(B_PER_BAG)) {
        if (have("shap_xgb_feat", "bootstrap_seeded", boot_id = b, seed_id = s)) {
          cat(sprintf("  bag %2d  seed %d  cached\n", b, s)); next
        }
        t0 <- start_timer()
        emit(shap_oof(SEED_B * s, bag = bag_for(b))$G, "shap_xgb_feat",
             "bootstrap_seeded", boot_id = b, seed_id = s)
        cat(sprintf("  bag %2d  seed %d  %.1f min\n", b, s, t0()$elapsed_sec / 60))
      }
    }
  }

  LLR_M <- intersect(METHODS, paste0("llr_", ATTR_LLR_ARMS))


  # --- LLR bootstrap on the shared bags: the PRIMARY level-3 route ----------
  if (length(LLR_M) && B_LBOOT > 0L) {
    cat(sprintf("\n=== ROUTE `bootstrap`: LLR on the same %d shared bags ===\n\n",
                B_LBOOT))
    cat("    THIS IS NOW THE PRIMARY LEVEL-3 QUANTITY FOR THE LLR ARMS. It is\n")
    cat("    the only route on which a level-4 comparison can be paired -- arm\n")
    cat("    A and arm B fitted on the identical resample -- and the only one\n")
    cat("    directly comparable with SHAP, which draws from these same bags.\n")
    cat("    It is also the route on which the `cond` arms' two models covary\n")
    cat("    correctly, which the posterior route cannot reproduce.\n\n")

    #' One arm on one bag, or the reason it could not be produced.
    #'
    #' THE WHOLE SPEC-FOLD IS INSIDE THE `try()`, AND IT WAS NOT. `signal_frame()`
    #' sat outside it, so `check_model_frame()`'s refusal propagated as an error
    #' and halted the run rather than being recorded. That refusal is not
    #' hypothetical and is not a defect: `smooth_k` is calibrated, per CLAUDE.md's
    #' frozen decision, from THE THINNEST OF THE FIVE FOLD FITTING SUBSETS, and a
    #' bootstrap bag is a 63% subsample OF one of those, thinner than anything the
    #' declaration was ever measured against. A count-valued duration covariate
    #' can therefore lose distinct values and fall below its declared basis --
    #' measured on bag 20, where `urine_output_rate`/`intv` had `diuretic__n_hours`
    #' at exactly 10 distinct values against a declared k of 10.
    #'
    #' LOWERING `smooth_k` IS NOT THE REMEDY. It is a field of the design key, so
    #' it would invalidate all 258 cached fits and every bundle, and it would be
    #' changing what the pipeline models in order to serve a diagnostic. The
    #' declaration is correct for the pipeline; the bootstrap route is simply
    #' fitting on data the pipeline never fits on.
    #'
    #' The old code's `next` on a fit failure is also gone, for the reason the
    #' caller documents: a skipped spec-fold leaves a column of zeros, and L = 0
    #' means UNMEASURED in this design.
    #' MEMOISED WITHIN A BAG AS OF 2026-09-09 (review finding A11). `memo` is
    #' an environment keyed by `.spec_key()` and fold: the first arm to need
    #' a (signal, formula, model, fold) fits it, every later arm that resolves
    #' to the same specification reads the result. An unpaired signal's
    #' `full`, `full_ti_trend` and `full_ti_all` are its `meas`, so this
    #' removes 24 of the 88 cells per fold. A failed fit is memoised too, so
    #' every arm that would have used it records its own exclusion row, as
    #' before, and the bag is still excluded whole.
    boot_arm <- function(arm, bag, memo) {
      M <- matrix(0, length(tr), length(sigs), dimnames = list(ids_ch, sigs))
      bad <- list()
      for (sg in sigs) {
        af <- arm_formula(arm, sg)
        if (!af$spec$fit) next
        key <- .spec_key(sg, af)
        for (k in seq_len(n_fold)) {
          mk  <- paste(key, k, sep = "|")
          r   <- memo[[mk]]
          hit <- !is.null(r)
          if (!hit) {
            r <- try({
              pri <- priors_for(priors, sg, "oof", fold = k)
              ji  <- job_ids("oof", k, folds)
              fit_ids <- intersect(ji$fit_ids, tr[bag])
              d_fit <- signal_frame(sg, af$spec$model, tabs, cfg, pri, stay_ids = fit_ids)
              d_prd <- signal_frame(sg, af$spec$model, tabs, cfg, pri,
                                    stay_ids = ji$predict_ids, stage = "predict")
              bm  <- .bam_fit(af$formula, d_fit, cfg)
              eta <- as.numeric(stats::predict(bm, newdata = d_prd, type = "link",
                                               discrete = FALSE))
              list(rows = match(as.character(d_prd$stay_id), ids_ch),
                   val = eta - logit(pri$p_bar))
            }, silent = TRUE)
            memo[[mk]] <- r
            memo$.n_fit <- (memo$.n_fit %||% 0L) + 1L
          } else {
            memo$.n_hit <- (memo$.n_hit %||% 0L) + 1L
          }
          if (inherits(r, "try-error")) {
            msg <- conditionMessage(attr(r, "condition"))
            bad[[length(bad) + 1L]] <- data.frame(
              arm = arm, signal = sg, fold = k,
              cause = if (grepl("distinct-value count", msg, fixed = TRUE))
                        "basis_exceeds_distinct_values" else "fit_error",
              detail = substr(gsub("\\s+", " ", msg), 1L, 200L),
              memo_hit = hit, stringsAsFactors = FALSE)
            next
          }
          M[r$rows, sg] <- r$val
        }
      }
      list(M = M, bad = if (length(bad)) do.call(rbind, bad) else NULL)
    }

    # The tombstone is read once, above the manifest section; see the block
    # that defines `TOMB`, `TOMB_ALL` and `tombstone_bag()`.
    if (nrow(TOMB)) cat("\n")

    BFAIL <- list()
    for (b in seq_len(B_LBOOT)) {
      if (b %in% TOMB$boot_id) {
        cat(sprintf("  bag %2d  excluded (tombstoned)\n", b)); next
      }
      if (all(vapply(LLR_M, function(m) have(m, "bootstrap", boot_id = b),
                     logical(1)))) {
        cat(sprintf("  bag %2d  cached\n", b)); next
      }
      t0 <- start_timer()
      bg <- bag_for(b)
      memo <- new.env(parent = emptyenv())
      # ATOMIC PER BAG: every arm is computed into memory and NOTHING is written
      # until the bag is known to be complete. Five matrices of 41,250 x 19 is
      # 31 MB, which is nothing, and the alternative is a bag that is partially
      # on disk. `intv` happens to be computed first, so bag 20's failure
      # aborted before any write -- but a failure in `full_ti_all` would have
      # left `meas` and `full` on disk for that bag, and `reconcile_manifest()`
      # would then have described a partially-populated bag as a bag.
      # `intv` IS FITTED ONLY WHEN SOMETHING SUBTRACTS IT (re-review, A11): a
      # configuration with no `cond` arm and no `llr_intv` has no use for the
      # 60 intervention spec-folds per bag. Every current configuration asks
      # for the `cond` arms, so this changes no fit today.
      need_intv <- any(c("llr_intv", paste0("llr_", ATTR_DERIVED_ARMS)) %in% LLR_M)
      bad <- list(); Iv <- NULL
      if (need_intv) {
        rv <- boot_arm("intv", bg, memo); Iv <- rv$M
        bad <- list(rv$bad)
      }
      keep <- list()
      if (need_intv && "llr_intv" %in% LLR_M) keep[["llr_intv"]] <- Iv
      for (arm in c("meas", "full", "full_ti_trend", "full_ti_all")) {
        m_direct <- paste0("llr_", arm); m_cond <- cond_of(arm)
        if (!length(intersect(c(m_direct, m_cond), LLR_M))) next
        ra <- boot_arm(arm, bg, memo); bad[[length(bad) + 1L]] <- ra$bad
        if (m_direct %in% LLR_M) keep[[m_direct]] <- ra$M
        if (!is.na(m_cond) && m_cond %in% LLR_M) keep[[m_cond]] <- ra$M - Iv
      }
      BAD <- do.call(rbind, Filter(Negate(is.null), bad))
      n_fit <- memo$.n_fit %||% 0L; n_hit <- memo$.n_hit %||% 0L
      rm(memo)
      mins <- round(t0()$elapsed_sec / 60, 2)

      # THE BAG IS EXCLUDED WHOLE, OR KEPT WHOLE. A spec-fold that produced no
      # prediction would leave a COLUMN OF ZEROS, and L = 0 means UNMEASURED in
      # this design (frozen decision: "Unmeasured stays get L = 0, assigned
      # directly -- never an indicator"). No metric here can tell an assigned
      # zero from a failed fit, so a zeroed column is a corruption that reads as
      # data. It is worse still for a `cond` arm: `cond = full - intv`, so a
      # missing `intv` fold silently turns `cond` into `full` for those rows.
      #
      # A replicate that differs from its siblings in WHAT IT CONTAINS rather
      # than in its draw is not a draw from the same distribution, so it cannot
      # be pooled with them. Excluding and COUNTING is the pattern this project
      # already uses for a unit that cannot support the analysis -- see the
      # hospital inclusion floors in `config/external.yml`, where a hospital
      # below the floor is excluded and counted rather than reported with a wide
      # interval.
      if (!is.null(BAD)) {
        BFAIL[[length(BFAIL) + 1L]] <- cbind(
          data.frame(boot_id = b, frac_in_bag = round(mean(bg), 5),
                     minutes = mins, stringsAsFactors = FALSE), BAD)
        tombstone_bag(b, nrow(BAD), paste(sort(unique(BAD$cause)), collapse = "+"))
        cat(sprintf("  bag %2d  EXCLUDED: %d spec-fold(s) unsupported (%s)  %.1f min\n",
                    b, nrow(BAD), paste(unique(BAD$cause), collapse = ", "), mins))
        for (q in seq_len(min(3L, nrow(BAD)))) {
          cat(sprintf("            %s / %s / fold %d: %s\n", BAD$arm[q],
                      BAD$signal[q], BAD$fold[q], BAD$cause[q]))
        }
        rm(keep, Iv); gc(verbose = FALSE)
        next
      }
      for (m in names(keep)) emit(keep[[m]], m, "bootstrap", boot_id = b)
      rm(keep, Iv); gc(verbose = FALSE)
      cat(sprintf("  bag %2d  in bag %s of patients  ok  %.1f min  (%d fits, %d alias hits)\n",
                  b, pct(mean(bg)), mins, n_fit, n_hit))
    }
    # THE FIT FAILURES ARE COUNTED AND SAVED, NOT SWALLOWED. `mgcv::bgam.fitd`
    # emitted `fitted probabilities numerically 0 or 1 occurred` on the 63%
    # bags during the 2026-09-06 validation refits and it did NOT occur on the
    # full folds, so separation is a property of the smaller fit set. A silent
    # `next` on a failed fit leaves a column of zeros, which is an ASSIGNED
    # value the metrics cannot distinguish from an unmeasured signal -- so the
    # count has to be visible per bag rather than inferred from a warning that
    # scrolled past. Plan section 21.7 flagged this exact risk if the bootstrap
    # route was ever promoted from validation to primary, which it now is.
    if (length(BFAIL)) {
      BF <- do.call(rbind, BFAIL)
      save_table(run, BF, "llr_bootstrap_excluded_bags", subdir = "diagnostics")
      nb <- length(unique(BF$boot_id))
      cat(sprintf("\n  *** %d BAG(S) EXCLUDED of %d: %s\n", nb, B_LBOOT,
                  paste(sort(unique(BF$boot_id)), collapse = ", ")))
      cat("      A bag is excluded WHOLE when any spec-fold cannot be produced,\n")
      cat("      because a skipped spec-fold leaves a column of zeros and L = 0\n")
      cat("      means UNMEASURED in this design. See\n")
      cat("      llr_bootstrap_excluded_bags.csv for the arm, signal and cause.\n")
      cat(sprintf("      THE EXCLUSION RATE IS ITSELF REPORTABLE: %.1f%% of bags\n",
                  100 * nb / B_LBOOT))
      cat("      cannot support the declared design. If it is large, that is a\n")
      cat("      finding about the bootstrap route and not a chore.\n")
      print(table(BF$cause))
    }
  }

  # --- LLR posterior draws: the cheap high-B companion ----------------------
  if (length(LLR_M) && B_DRAW > 0L && identical(ROUTE, "posterior")) {
    cat(sprintf("\n=== ROUTE `posterior`: LLR estimation noise, B = %d ===\n\n", B_DRAW))
    cat("    Contrast L3P, NOT L3. It conditions on the observed sample. Kept\n")
    cat("    because 40 draws cost less than two bootstrap replicates and its\n")
    cat("    within-method spread is worth having at high B; never pooled with\n")
    cat("    the bootstrap route and never compared with SHAP.\n\n")

    #' Draw from N(mu, V) with a clamped spectrum.
    #'
    #' `Vc` is symmetric but numerically indefinite often enough that a Cholesky
    #' fails on a real fit. Eigen with negatives clamped to zero is the standard
    #' repair; the number clamped is RETURNED rather than swallowed, so a fit
    #' whose covariance was materially repaired is visible in the diagnostics
    #' table instead of silently producing narrow draws.
    rmvn_clamped <- function(nd, mu, V) {
      e <- eigen(V, symmetric = TRUE)
      neg <- sum(e$values < 0)
      d <- sqrt(pmax(e$values, 0))
      Z <- matrix(stats::rnorm(nd * length(mu)), nd, length(mu))
      list(draws = sweep(Z %*% (t(e$vectors) * d), 2L, mu, `+`), n_clamped = neg)
    }

    #' B posterior L matrices for one base arm, out of fold.
    #'
    #' The fold fit lives inside this loop and dies with it (hard rule 6). What
    #' leaves is B columns of held-out log-odds per spec-fold.
    #' MEMOISED ACROSS ARMS AS OF 2026-09-09 (review finding A11), for the
    #' specifications in `SHARED_SPECS` only. The alias identity this route
    #' always claimed -- an unpaired signal's `full_ti_all` column equals
    #' its `full` column bitwise -- rested on refitting the same formula and
    #' drawing under the same seed `DRAW0 + ctr`. The memo keeps the first
    #' arm's draws and hands them to the next; the counter still advances,
    #' and a hit whose recorded counter differs from the live one is refused
    #' rather than substituted, so the memo can only ever return what the
    #' refit would have. `tests/attr_eval_unit.R` asserts the bam
    #' determinism the identity rests on.
    posterior_arm <- function(arm, nd, memo) {
      L <- replicate(nd, matrix(0, length(tr), length(sigs),
                                dimnames = list(ids_ch, sigs)), simplify = FALSE)
      diag <- list(); ctr <- 0L
      for (sg in sigs) {
        af <- arm_formula(arm, sg)
        if (!af$spec$fit) next     # unpaired `intv`: assigned zero, never fitted
        key <- .spec_key(sg, af); shared <- key %in% SHARED_SPECS
        for (k in seq_len(n_fold)) {
          ctr <- ctr + 1L
          mk <- paste(key, k, sep = "|")
          h  <- if (shared) memo[[mk]] else NULL
          if (!is.null(h)) {
            if (!identical(h$ctr, ctr)) {
              stop("posterior_arm: memo hit for `", sg, "` fold ", k, " at ",
                   "draw counter ", ctr, " but it was drawn at counter ", h$ctr,
                   ". The common-random-numbers alignment the alias identity ",
                   "depends on has broken; refusing to substitute.", call. = FALSE)
            }
            rows <- h$rows; Eta <- h$Eta
            dg <- h$diag; dg$arm <- arm; dg$memo_hit <- TRUE
          } else {
            pri <- priors_for(priors, sg, "oof", fold = k)
            ji  <- job_ids("oof", k, folds)
            d_fit <- signal_frame(sg, af$spec$model, tabs, cfg, pri, stay_ids = ji$fit_ids)
            d_prd <- signal_frame(sg, af$spec$model, tabs, cfg, pri,
                                  stay_ids = ji$predict_ids, stage = "predict")
            b  <- .bam_fit(af$formula, d_fit, cfg)
            # `Vc` and NOT `Vp`. `Vc` carries the smoothing-parameter uncertainty
            # correction, which is the whole reason a posterior draw is a
            # defensible stand-in for a refit. It must be read OFF THE OBJECT
            # here: `mgcv:::predict.bam` sets `object$Vc <- NULL` before
            # delegating, so `unconditional = TRUE` is a silent no-op on a `bam`
            # (see tests/coupling_displacement.R's `.partial()`), and anything
            # that goes through predict() to get it comes back uncorrected.
            V  <- if (!is.null(b$Vc)) b$Vc else b$Vp
            Xp <- stats::predict(b, newdata = d_prd, type = "lpmatrix", discrete = FALSE)
            dr <- with_seed(DRAW0 + ctr, rmvn_clamped(nd, stats::coef(b), V))
            rows <- match(as.character(d_prd$stay_id), ids_ch)
            Eta  <- Xp %*% t(dr$draws) - logit(pri$p_bar)
            dg <- data.frame(arm = arm, signal = sg, fold = k,
              n_coef = length(stats::coef(b)), used_vc = !is.null(b$Vc),
              n_eigen_clamped = dr$n_clamped,
              se_median = round(stats::median(apply(Eta, 1, stats::sd)), 5),
              memo_hit = FALSE, stringsAsFactors = FALSE)
            if (shared) memo[[mk]] <- list(ctr = ctr, rows = rows, Eta = Eta, diag = dg)
          }
          for (i in seq_len(nd)) L[[i]][rows, sg] <- Eta[, i]
          diag[[length(diag) + 1L]] <- dg
        }
        cat(sprintf("    %-18s %-14s %d folds\n", sg, arm, n_fold))
      }
      list(L = L, diag = do.call(rbind, diag))
    }

    DG <- list()
    PMEMO <- new.env(parent = emptyenv())
    # `intv` first and cached to disk, so `cond` can be formed later without
    # holding two full sets of B matrices in memory at once.
    #
    # KEYED AND COMPLETE AS OF 2026-09-09 (review finding A2). The cache lived
    # at `_intv/intv_XX.qs2` under no design key, and only the LAST file's
    # existence decided whether all of it was reusable, so a resume across a
    # design or draw-seed change could pair new `full` draws with old `intv`
    # draws inside `cond`. The directory is now named by the design key, the
    # draw seed, the draw count and the draw layout, every file must exist,
    # and each is checked against the store's row and column identity before
    # it is subtracted. An old unkeyed cache is left in place and ignored.
    intv_dir <- file.path(run$path, SUB, "_intv", attr_key_hash(list(
      design = DESIGN, draw_seed_base = DRAW0, b = B_DRAW,
      layout = ATTR_GENERATOR_VERSION$posterior_draw_layout)))
    intv_files <- file.path(intv_dir, sprintf("intv_%02d.qs2", seq_len(B_DRAW)))
    fully_cached <- function(m) all(vapply(seq_len(B_DRAW),
      function(i) have(m, "posterior", draw_id = i), logical(1)))
    cond_needed <- Filter(function(m) !fully_cached(m),
                          intersect(LLR_M, paste0("llr_", ATTR_DERIVED_ARMS)))
    if (length(cond_needed)) {
      dir.create(intv_dir, showWarnings = FALSE, recursive = TRUE)
      if (!all(file.exists(intv_files))) {
        cat("\n  base arm `intv` (needed for every cond arm)\n")
        t0 <- start_timer()
        r <- posterior_arm("intv", B_DRAW, PMEMO); DG[[length(DG) + 1L]] <- r$diag
        for (i in seq_len(B_DRAW)) qs2::qs_save(r$L[[i]], intv_files[i])
        rm(r); gc()
        cat(sprintf("  intv done in %.1f min\n", t0()$elapsed_sec / 60))
      } else {
        cat("\n  base arm `intv` cached (keyed)\n")
      }
    }

    for (arm in c("meas", "full", "full_ti_trend", "full_ti_all")) {
      m_direct <- paste0("llr_", arm); m_cond <- cond_of(arm)
      want <- intersect(c(m_direct, m_cond), LLR_M)
      if (!length(want)) next
      if (all(vapply(want, fully_cached, logical(1)))) {
        cat(sprintf("\n  base arm `%s` fully cached\n", arm)); next
      }
      cat(sprintf("\n  base arm `%s`\n", arm))
      t0 <- start_timer()
      r <- posterior_arm(arm, B_DRAW, PMEMO); DG[[length(DG) + 1L]] <- r$diag
      for (i in seq_len(B_DRAW)) {
        if (m_direct %in% want) emit(r$L[[i]], m_direct, "posterior", draw_id = i)
        if (!is.na(m_cond) && m_cond %in% want) {
          Iv <- qs2::qs_read(intv_files[i])
          attr_check_replicate(Iv, tr, sigs, sprintf("intv posterior draw %d", i))
          emit(r$L[[i]] - Iv, m_cond, "posterior", draw_id = i)
        }
      }
      rm(r); gc()
      cat(sprintf("  %s done in %.1f min\n", arm, t0()$elapsed_sec / 60))
    }
    rm(PMEMO); gc(verbose = FALSE)
    if (length(DG)) save_table(run, do.call(rbind, DG), "posterior_fit_diagnostics",
                               subdir = "diagnostics")
  }
}

# --- close -------------------------------------------------------------------
#
# THE MANIFEST IS DERIVED FROM THE PLAN AND THE FILESYSTEM, whatever `--levels`
# this invocation asked for. See `reconcile_manifest()` for the defect that
# makes this necessary rather than tidy.
MAN <- reconcile_manifest(PLAN)
MAN <- MAN[order(MAN$method, MAN$route, MAN$index), , drop = FALSE]
utils::write.csv(MAN, MAN_P, row.names = FALSE)
attr_validate_manifest(MAN, "tests/attr_replicates.R")
cat(sprintf("\n=== %d replicates on disk in %s ===\n\n", nrow(MAN), basename(run$path)))
print(table(MAN$method, MAN$route))

# --- COVERAGE: planned against present, and the fingerprint stamp -----------
#
# ADDED 2026-09-09 (review finding A9). A manifest lists what exists and says
# nothing about what was planned and is not there, so a store missing half
# its bootstrap replicates looked complete. Every planned coordinate is now
# present, tombstoned with a cause, or MISSING, and the run manifest says
# which.
COV <- attr_replicate_coverage(PLAN, MAN, tombstoned = TOMB$boot_id)
save_table(run, COV, "replicate_coverage", subdir = "diagnostics")
save_table(run, ATTR_ESTIMAND_NOTES, "design_notes", subdir = "diagnostics")
cat("\n  coverage (planned / present / tombstoned / missing):\n\n")
print(COV, row.names = FALSE)
STORE_COMPLETE <- all(COV$complete)
if (!STORE_COMPLETE) {
  cat(sprintf("\n  *** %d planned replicate(s) MISSING. The store is PARTIAL; the\n",
              sum(COV$missing)))
  cat("      consumer refuses it without --allow-partial.\n")
}
if (attr(COV, "n_unplanned_on_disk") > 0L) {
  cat(sprintf("  note: %d replicate(s) on disk are not in the current plan (ignored).\n",
              attr(COV, "n_unplanned_on_disk")))
}

# The fingerprint stamp, with its evidence, for a store that predates it.
if (NEED_STAMP) {
  cat("\n=== fingerprint stamp ===\n\n")
  llr_m <- intersect(METHODS, paste0("llr_", ATTR_LLR_ARMS))
  if (!("spec" %in% LEVELS_WANTED)) {
    cat("  verifying the LLR ladder replicates on disk (no write):\n")
    llr_ladder(write = FALSE)
  }
  EV <- if (length(EVID)) do.call(rbind, EVID) else NULL
  ok_llr <- !is.null(EV) && all(vapply(llr_m, function(m)
    any(EV$method == m & EV$route == "ladder" & EV$status == "verified_bitwise"),
    logical(1)))
  ok_shap <- !("shap_xgb_feat" %in% METHODS) ||
    (file.exists(.shap_margin_path()) && !is.null(EV) &&
     any(EV$method == "shap_xgb_feat" & EV$route == "ladder" &
         EV$status == "verified_bitwise"))
  if (ok_llr && ok_shap) {
    file.copy(DESIGN_P, file.path(run$path, "design_pre_fingerprint.qs2"), overwrite = FALSE)
    stamped <- OLD_DG
    stamped$fingerprint <- FP
    stamped$fingerprint_evidence <- EV[EV$route == "ladder", , drop = FALSE]
    stamped$fingerprint_stamped <- "2026-09-09 ladder verified bitwise (LLR from targets cache, SHAP by refit)"
    stamped$naming_design_key <- NAMING_DESIGN
    qs2::qs_save(stamped, DESIGN_P)
    cat(sprintf("  design.qs2 stamped with fingerprint %s; previous file kept as design_pre_fingerprint.qs2\n",
                attr_key_hash(FP)))
  } else {
    cat("  NOT stamped. Verification evidence is incomplete:\n")
    cat(sprintf("    LLR ladder verified : %s\n", ok_llr))
    cat(sprintf("    SHAP ladder verified: %s\n", ok_shap))
    cat("  Re-run with `--levels spec --resume <dir>` (about two minutes) to verify\n")
    cat("  both halves; the consumer refuses the store until then.\n")
  }
}
# The re-stamp after a route-selective quarantine. The files the old
# fingerprint described are gone from the store, so the store now describes
# the live design; whether every regenerated coordinate is back is what the
# coverage table above says, and a partial store is refused by the consumer
# exactly as before.
if (!is.null(STALE_RULES)) {
  old_h <- substr(attr_key_hash(OLD_DG$fingerprint), 1, 12)
  file.copy(DESIGN_P, file.path(run$path, paste0("design_pre_fingerprint_", old_h, ".qs2")),
            overwrite = FALSE)
  stamped <- qs2::qs_read(DESIGN_P)
  stamped$fingerprint <- FP
  stamped$fingerprint_restamped <- sprintf(
    "stale routes quarantined (%d file(s) into stale_%s) and regenerated: %s",
    STALE_MOVED, old_h,
    paste(sprintf("%s(%s)", STALE_RULES$route, STALE_RULES$family), collapse = ","))
  stamped$fingerprint_evidence <- rbind(
    stamped$fingerprint_evidence,
    if (length(EVID)) do.call(rbind, EVID)[, c("method", "route", "status", "max_abs_diff")])
  qs2::qs_save(stamped, DESIGN_P)
  cat(sprintf("\n  design.qs2 re-stamped with fingerprint %s (previous kept as design_pre_fingerprint_%s.qs2)\n",
              attr_key_hash(FP), old_h))
  if (!STORE_COMPLETE) {
    cat("  the store is PARTIAL until the quarantined routes are regenerated:\n")
    cat("  resume with the levels those routes belong to.\n")
  }
}
if (length(EVID)) {
  save_table(run, do.call(rbind, EVID), "emit_evidence", subdir = "diagnostics")
}
cat("\n  contrasts this store supports (within-method pairs):\n\n")
CP <- attr_contrast_pairs(MAN)
if (nrow(CP)) {
  # `route_family` AND `held`, NOT `route`. `attr_contrast_pairs()` stopped
  # returning a `route` column when route grouping became route-FAMILY grouping
  # (the L3T fix); the consumer was updated and this print was not, so
  # `CP$route` was NULL and `aggregate()` stopped with "arguments must have
  # same length" -- AFTER every replicate and the manifest were written but
  # BEFORE `finalize_run()`. Same rename class as the two before it: a field
  # renamed in a library and followed in one caller of two.
  print(aggregate(list(n_pairs = CP$i),
                  by = list(method = CP$method, family = CP$route_family,
                            contrast = CP$code, stratum = CP$held),
                  FUN = length), row.names = FALSE)
} else {
  cat("    none yet -- generate stage `sample`.\n")
}
finalize_run(run, extra = list(
  design_key = attr_key_hash(DESIGN),
  fingerprint = attr_key_hash(FP),
  levels_generated = paste(LEVELS_WANTED, collapse = ","),
  n_replicates = nrow(MAN),
  n_planned = sum(COV$planned),
  n_missing = sum(COV$missing),
  n_tombstoned = sum(COV$tombstoned),
  store_complete = STORE_COMPLETE,
  tombstoned_bags = paste(sort(TOMB$boot_id), collapse = ","),
  n_shared_bags = B_BOOT,
  bag_manifest_hash = attr_key_hash(BAGT$membership_hash),
  k_ti = K_TI,
  pass_fits_distinct = PASS_FITS,
  pass_cells_enumerated = PASS_CELLS,
  ladder_reused = if (is.null(LADDER_D)) "" else basename(LADDER_D),
  resample_kind = ATTR_ESTIMAND_NOTES$value[ATTR_ESTIMAND_NOTES$field == "resample_kind"],
  held_fixed_under_resampling = ATTR_ESTIMAND_NOTES$value[ATTR_ESTIMAND_NOTES$field == "held_fixed_under_resampling"],
  resumed = !is.null(RESUME)))
cat(sprintf("\nwritten: %s\n", run$path))
