# tests/attr_external.R -------------------------------------------------------
# THE ATTRIBUTION ARM AT eICU. Applies the frozen bundle, builds every arm the
# bundle carries, and measures level 1 and level 4 with the same library the
# MIMIC analysis uses. FITS NOTHING.
#
# Read `docs/attribution_analysis_plan_20260906.md` PART FOUR first, then
# `docs/interaction_spec_change_20260907.md`.
#
# --- WHAT TRANSPORTS, AND WHAT DOES NOT ------------------------------------
#
# ALL EIGHT ARMS EXIST HERE AS OF 2026-09-07. This header said the opposite
# until that date, and correctly:
#
#     "FOUR OF THE EIGHT ARMS EXIST HERE. The bundle carries 43 frozen specs
#      and every one of them is `meas`, `full` or `intv`; NO INTERACTION SPEC
#      IS IN IT [...] THE INTERACTION QUESTION CANNOT BE ASKED AT eICU AT ALL.
#      Getting it would need four new specs fitted at full train and frozen
#      into a new bundle -- a `_targets.R` change and a new design key -- or a
#      fit at eICU, which hard rule 8 forbids."
#
# That is exactly what was done. `full_ti_trend` and `full_ti_all` joined
# `LAYER1_MODELS`, the bundle now carries 64 frozen specs, and `apply_bundle()`
# builds all six L matrices at any apply site. So the seven LLR arms and
# `shap_xgb_feat` are all available and there are 28 level-4 pairs here, the
# same 28 as at MIMIC -- which is what makes the two sites' level-4 tables
# comparable row for row rather than as a subset against a whole.
#
# NOTHING IS FITTED HERE AND THE HARD-RULE-8 POSITION IS UNCHANGED. The
# interaction smooths were fitted at MIMIC on full train and frozen; eICU
# evaluates them per patient exactly as it evaluates every other smooth.
#
# LEVELS 2, 3 AND 3T DO NOT EXIST AT AN APPLY SITE, and not merely because they
# would be expensive. THE ATTRIBUTION OF A GIVEN PATIENT IS DETERMINISTIC HERE.
# Every parameter is frozen, so there is no estimation noise per patient to
# measure; resampling the eICU EVALUATION cohort would change only which
# patients enter a summary, not any individual attribution, and that is ordinary
# sampling error of a mean rather than the estimator variability level 3 is
# about. Level 3 is a property of the FITTING PROCEDURE and is fully answered at
# the training site. It is not a site-transportable quantity and this script
# does not pretend to produce one.
#
# LEVEL 4 TRANSPORTS FOR FREE, because it holds the data fixed and varies the
# specification: apply two frozen arms to the same eICU rows and compare. That
# is the result this script exists for.
#
# --- WHY THE CROSS-SITE COMPARISON IS DISTRIBUTIONAL AND NOT PAIRED --------
#
# Every metric in `R/14_attribution_eval.R` compares two matrices ROW BY ROW,
# same patient. Arm A against arm B WITHIN eICU is therefore exact. Arm A at
# MIMIC against arm A at eICU is not a comparison this family can make -- there
# are no shared patients. So the transport reading is a comparison of SUMMARIES:
# does the level-4 ladder have the same shape and the same ordering at both
# sites, and does resolution survive? Those are stated as such and never as a
# per-patient agreement.
#
# Aggregates only (hard rule 1).
#
#   Rscript tests/attr_external.R
#   Rscript tests/attr_external.R --mimic out/runs/attrmetrics_...
# ----------------------------------------------------------------------------

# `mgcv` IS NOT OPTIONAL HERE, and omitting it does not fail where you would
# expect. The bundle's models are class c("bam","gam","glm","lm"); without mgcv
# ATTACHED, `predict()` cannot find `predict.bam` or `predict.gam` and falls
# through to `predict.glm`, which demands the `qr` component that a stripped
# `bam` does not carry. The error surfaces deep inside `apply_bundle()` as
# "lm object does not have a proper 'qr' component", which reads like a corrupt
# bundle rather than a missing import.
suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
# COMPLETE RUNS ONLY, AND PROVENANCE CHECKED BELOW (review finding A10). The
# default was `require_complete = FALSE`, so an interrupted metrics run that
# happened to be newest would have been the comparator, and the join was on
# method NAMES alone: nothing said the MIMIC store and this bundle described
# the same design.
MIMIC_RUN <- .opt("--mimic", latest_run("attrmetrics", require_complete = TRUE))

ecfg  <- yaml::read_yaml("config/attribution_eval.yml")
xcfg  <- yaml::read_yaml("config/external.yml")
COL_SH   <- as.numeric(cfg_req(ecfg, "delta", "collapse_shares"))
DELTA_ABS <- as.numeric(cfg_req(ecfg, "delta", "common_absolute"))
TAUS     <- as.numeric(cfg_req(ecfg, "delta", "sign_taus"))
KS       <- as.integer(cfg_req(ecfg, "aggregation", "top_k"))
LVLS     <- as.character(cfg_req(ecfg, "aggregation", "levels"))

# --- the bundle, and the design it froze -------------------------------------
#
# THE DESIGN COMES FROM THE BUNDLE AND NOT FROM config/config.yml. That is the
# frozen decision `bundle_cfg()` exists to enforce: without it, editing a
# whitelist or a `smooth_k` entry between the MIMIC fit and this run would
# change what is being transported while changing nothing visible.
bpath  <- as.character(cfg_req(xcfg, "bundle"))
bundle <- load_bundle(bpath)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
cat(sprintf("\n=== eICU attribution, bundle %s ===\n\n", basename(dirname(bpath))))

# `load_tables(paths, cfg, site, verbose)` -- four arguments, as
# `run/external.R` calls it. Passing `cfg` alone put the config's top-level
# names where the path list belonged, and the missing-entry guard caught it.
tabs      <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
stay_ids  <- tabs$cohort$stay_id
ids_ch    <- as.character(stay_ids)
sigs      <- as.character(unlist(cfg$signals))
domains   <- bundle$domains
# THERE IS NO OUTCOME VECTOR HERE, AND ITS ABSENCE IS THE POINT. `y_eicu` was
# assigned on this line until 2026-09-07 and read by nothing. That is worth a
# sentence rather than a silent deletion: every metric in this script compares
# two attribution matrices with each other, `mortality` appears in none of them,
# and that is exactly why the arm spends no single-look budget
# (config/attribution.yml's argument, applied at the external site). A dangling
# outcome vector invites the next reader to use it and quietly turn an agreement
# statistic into a performance one.
cat(sprintf("  eICU stays scored : %d\n", length(stay_ids)))
cat(sprintf("  signals           : %d\n", length(sigs)))

# --- the seven LLR arms, from frozen parameters ------------------------------
#
# ARM NAMES ARE `R/14`'s, NOT `R/10`'s, and the difference is deliberate rather
# than an oversight. `apply_bundle()` calls the joint arm `llr_sum` because it
# is a row SUM of a score; the attribution library calls it `llr_full` because
# it is the `full` MODEL's matrix, and every contrast in `R/14` is keyed on the
# model. `LLR_ARM_MATRIX` is the map between the two and is read here rather
# than restated, so a new arm cannot be added to the bundle and forgotten here.
t0 <- start_timer()
ap <- apply_bundle(bundle, tabs, cfg, stay_ids,
                   arms = names(LLR_ARM_MATRIX), verbose = FALSE)
cat(sprintf("  apply_bundle      : %.1f min\n", t0()$elapsed_sec / 60))
# THE ARM SET COMES FROM `config/attribution_eval.yml`, not from
# `ATTR_LLR_ARMS`. The library constant includes `intv`, which the analysis
# deliberately excludes -- a pure intervention-propensity contrast is not an
# attribution of measurement evidence -- and the point of this script is a
# level-4 table over the SAME arm set MIMIC uses. Reading the config is what
# makes the two sites' pair sets identical by construction rather than by two
# lists that happen to agree today.
.mw   <- as.character(cfg_req(ecfg, "methods"))
attr_check_methods(.mw)
.llrw <- vapply(Filter(function(m) attr_method_family(m) == "llr", .mw),
                attr_arm_of, character(1))
A <- stats::setNames(
  lapply(.llrw, function(a) {
    M <- ap$l_mats[[a]]
    if (is.null(M)) {
      stop("attr_external: the bundle produced no `", a, "` L matrix. Every ",
           "arm the attribution library names must be present, or the two ",
           "sites' level-4 tables are over different pair sets and the ",
           "transport reading is a subset against a whole.", call. = FALSE)
    }
    M[, sigs, drop = FALSE]
  }),
  paste0("llr_", .llrw))

# --- SHAP at eICU, from the frozen full-train booster ------------------------
#
# The same rollup `shap_oof()` uses at MIMIC, with the FROZEN booster in place
# of the five per-fold ones and `role = "final"` priors from the bundle. The
# 10.1% intervention-mass truncation of plan section 41 applies here unchanged
# and is deliberate. It makes this object STRUCTURALLY analogous to `llr_cond`
# -- measurement groups with the intervention block set aside -- and no more
# than that: dropping SHAP groups is not the subtraction of two fitted models,
# and the two are not the same estimand (review finding A4). Nothing here
# reads the row sums of this matrix as the booster's score.
t0 <- start_timer()
m  <- bundle$xgb$xgb_feat
X  <- xgb_design_feat(tabs, cfg, bundle$priors, stay_ids, role = "final",
                      fold = NA_integer_, feature_names = m$feature_names)
X  <- align_design(X, m$feature_names)   # as `xgb_apply()` does, defensively
# `m$booster`, NOT `m$model` -- the bundle's xgb entry is
# (booster, feature_names, p_bar, best_iter, params, n_train).
#
# NO `iterationrange`, MATCHING `shap_oof()` AT MIMIC. The booster carries
# best_iteration = 449 and `.xgb_predict()` passes the range explicitly, which
# looked like a divergence: SHAP would describe all trees while the score
# described 450. MEASURED 2026-09-07 -- `predict()` in this xgboost version
# honours `best_iteration` by itself, for scores AND for `predcontrib`, to a
# difference of exactly zero. So the two paths already agree and adding the
# range would change nothing.
ctr <- stats::predict(m$booster, xgboost::xgb.DMatrix(X, missing = NA),
                      predcontrib = TRUE)
ref <- colnames(X)
stopifnot(ncol(ctr) == length(ref) + 1L)      # the trailing BIAS column
S   <- ctr[, seq_along(ref), drop = FALSE]
grp <- sub("__.*$", "", ref)
G   <- matrix(0, nrow(S), length(sigs), dimnames = list(ids_ch, sigs))
for (g in sigs) { j <- which(grp == g); if (length(j)) G[, g] <- rowSums(S[, j, drop = FALSE]) }
A$shap_xgb_feat <- G

# THE DISCARDED MASS IS MEASURED, NOT ASSUMED. Every column whose group is not
# one of the 19 signals -- the intervention block -- is dropped by the loop
# above, and that truncation is deliberate: it is what makes SHAP-19 the
# structural counterpart of `llr_cond` rather than of `llr_sum`.
#
# It was SILENT here until 2026-09-07, and silence is the problem. Plan section
# 41 measures the discarded share at MIMIC as 10.1%, and every reading that
# compares the two sites' SHAP arms assumes eICU's share is comparable. eICU's
# intervention coverage is not MIMIC's -- `invasive_vent__exposure_frac` alone
# behaves completely differently -- so the share is a quantity that can move,
# and one that must be reported beside any cross-site SHAP claim rather than
# carried over from the other site's document.
.kept  <- grp %in% sigs
.absS  <- abs(S)
.trunc <- data.frame(
  site = "eicu",
  n_cols_total = length(ref),
  n_cols_kept  = sum(.kept),
  n_cols_dropped = sum(!.kept),
  # Share of TOTAL |contribution| mass that the rollup discards, per patient,
  # summarised over patients. A mean over patients would be dominated by the
  # few with large budgets.
  drop_frac_median = round(stats::median(
    rowSums(.absS[, !.kept, drop = FALSE]) /
      pmax(rowSums(.absS), .Machine$double.eps)), 5),
  drop_frac_p90 = round(stats::quantile(
    rowSums(.absS[, !.kept, drop = FALSE]) /
      pmax(rowSums(.absS), .Machine$double.eps), 0.90, names = FALSE), 5),
  drop_frac_pooled = round(sum(.absS[, !.kept, drop = FALSE]) /
                             sum(.absS), 5),
  stringsAsFactors = FALSE)
rm(.absS)
cat(sprintf("  shap rollup       : %d of %d columns kept; discarded |mass| ",
            .trunc$n_cols_kept, .trunc$n_cols_total))
cat(sprintf("median %.4f, p90 %.4f, pooled %.4f\n",
            .trunc$drop_frac_median, .trunc$drop_frac_p90, .trunc$drop_frac_pooled))
cat(sprintf("  shap contributions: %.1f min\n\n", t0()$elapsed_sec / 60))

meas <- measured_matrix(tabs, cfg, stay_ids)

# ROW AND COLUMN IDENTITY, NOT JUST SHAPE. The check here compared `dim()` only,
# which is the weakest form of the property every metric in this file depends
# on: `attr_agree_k()`, `attr_cosine()` and `attr_leader_collapse()` all compare
# ROW i OF A AGAINST ROW i OF B, so two matrices of the same shape in different
# patient orders produce a complete table of confident, meaningless numbers.
#
# The risk is real and one-sided. The seven LLR arms come from one
# `l_matrices()` call and cannot disagree with each other. The SHAP arm is built
# on a different path entirely -- `xgb_design_feat()` scatters onto `ids` and
# `align_design()` reorders columns -- so it is the one that could drift, and it
# is exactly the arm every cross-family comparison in this script involves.
#
# TRACED 2026-09-07 and currently correct: `xgb_design_feat()` sets
# `rownames(X) <- ids` after a `match(ids, ...)` scatter, and `align_design()`
# preserves rownames. Asserted anyway, because "I traced it once" is not a
# property of the code.
for (nm in names(A)) {
  if (!identical(dim(A[[nm]]), dim(A$llr_meas))) {
    stop("attr_external: arm `", nm, "` is not shaped like `llr_meas`",
         call. = FALSE)
  }
  if (!identical(rownames(A[[nm]]), ids_ch)) {
    stop("attr_external: arm `", nm, "` is not in `stay_ids` row order. Every ",
         "metric here compares row i against row i, so a reordered arm ",
         "produces a full table of numbers about the wrong patients.",
         call. = FALSE)
  }
  if (!identical(colnames(A[[nm]]), sigs)) {
    stop("attr_external: arm `", nm, "` does not carry the signals in ",
         "`cfg$signals` order.", call. = FALSE)
  }
}

run <- new_run("attrext", cfg, note = sprintf(
  "eICU attribution: level 1 and level 4 for the %d arms the bundle carries",
  length(A)))
save_table(run, data.frame(arm = names(A), n_rows = nrow(A[[1]]),
                           n_cols = ncol(A[[1]]), stringsAsFactors = FALSE),
           "arms_scored", subdir = "diagnostics")
save_table(run, .trunc, "shap_rollup_truncation", subdir = "diagnostics")

to_level <- function(M, lv) if (lv == "signal") M else attr_to_domain(M, domains)
keep_at  <- function(lv) if (lv == "signal") meas else attr_to_domain(meas * 1, domains) > 0
# `keep_for()` WAS HERE AND WAS CALLED BY NOTHING. Removed 2026-09-07 with the
# rest of this pass: the mask a pair needs depends on BOTH arms (LLR-vs-SHAP
# takes no mask at all), so a per-arm helper could never have been the thing the
# pair loop wanted, and the loop had always computed the pair mask inline.

# --- EVERY LEVELLED MATRIX AND EVERY PREP, BUILT ONCE ------------------------
#
# `attr_prep()` MEMOISES WITHIN ONE PREP OBJECT AND NOT ACROSS CALLS. Its cache
# (`.cS` / `.cK`) lives in that prep's own environment and dies with it, so
# calling `attr_prep(A[[a]])` again rebuilds `ord <- t(apply(-ab, 1, order))`
# from scratch -- an `order()` per patient over 95,507 rows, which is the single
# most expensive operation in this script.
#
# The pair loop below called it TWICE PER PAIR PER LEVEL: 28 pairs x 2 levels x
# 2 sides is 112 preps where only 8 arms x 2 levels = 16 distinct ones exist.
# Seven eighths of the dominant cost, spent reproducing matrices already in
# memory. `to_level()` was recomputed the same way.
#
# THE CLASS: A MEMOISED FUNCTION WHOSE CACHE IS NARROWER THAN ITS CALL PATTERN.
# The memoisation note inside `attr_prep()` is about repeated `inS(k)` lookups
# on ONE prep, which is a different axis from repeated preps of one matrix, and
# reading "memoised" at the call site is what makes the redundancy invisible.
# Hoisting is the fix; enlarging the cache would key it on a 1.8 million-cell
# matrix.
AL <- PP <- KP <- list()
for (lv in LVLS) {
  KP[[lv]] <- keep_at(lv)
  AL[[lv]] <- lapply(A, to_level, lv = lv)
  PP[[lv]] <- lapply(AL[[lv]], attr_prep)
}

# ============================================================================
# LEVEL 1 at eICU: the budget and the resolution ceiling
# ============================================================================
cat("=== L1 at eICU: evidence budget and resolution ===\n\n")
B1 <- R1 <- list()
for (nm in names(A)) for (lv in LVLS) {
  B1[[length(B1) + 1L]] <- cbind(site = "eicu", arm = nm, agg = lv,
                                 attr_budget(AL[[lv]][[nm]]))
  R1[[length(R1) + 1L]] <- cbind(site = "eicu", arm = nm, agg = lv,
                                 attr_tieset(PP[[lv]][[nm]], DELTA_ABS))
}
B1 <- do.call(rbind, B1); R1 <- do.call(rbind, R1)
save_table(run, B1, "l1_budget_eicu", subdir = "diagnostics")
save_table(run, R1, "l1_resolution_eicu", subdir = "diagnostics")
z <- B1[B1$agg == "signal", ]
cat(sprintf("  %-16s %10s %10s %10s\n", "arm", "budget", "budget_p90", "hhi"))
for (i in seq_len(nrow(z))) cat(sprintf("  %-16s %10.4f %10.4f %10.5f\n",
  z$arm[i], z$budget_median[i], z$budget_p90[i], z$hhi_median[i]))

# RESOLUTION IS THE ONE TO WATCH. eICU's coverage is thinner -- nine lab signals
# fall from a median of two measured hours to one -- and thinner coverage means
# smaller |L|, which means more near-ties, which means fewer patients with a
# unique leader. If resolution collapses here then PER-PATIENT ATTRIBUTION
# CLAIMS DO NOT TRANSPORT even where the level-4 ladder does, and that is a
# limitation of the method rather than a property of one site.
cat("\n  resolution (unique leader), signal level:\n")
z <- R1[R1$agg == "signal" & R1$delta == 0.10, ]
for (i in seq_len(nrow(z))) cat(sprintf("  %-16s delta 0.10  unique leader %.4f  median tie set %.0f\n",
  z$arm[i], z$frac_unique_leader[i], z$tieset_median[i]))

# ============================================================================
# LEVEL 4 at eICU: specification, the whole metric family
# ============================================================================
cat("\n=== L4 at eICU: specification, on identical eICU rows ===\n\n")
CMP <- utils::combn(names(A), 2L, simplify = FALSE)
L4 <- COS <- FLIP <- COLL <- list()
for (p in CMP) for (lv in LVLS) {
  Aa <- AL[[lv]][[p[1]]]; Bb <- AL[[lv]][[p[2]]]
  pa <- PP[[lv]][[p[1]]]; pb <- PP[[lv]][[p[2]]]
  kp <- if (attr_method_family(p[1]) == "llr" &&
            attr_method_family(p[2]) == "llr") KP[[lv]] else NULL
  for (k in KS) {
    r <- vapply(DELTA_ABS, function(d) mean(attr_agree_k(pa, pb, k, d)), numeric(1))
    L4[[length(L4) + 1L]] <- data.frame(
      site = "eicu", method_a = p[1], method_b = p[2], agg = lv, k = k,
      delta = DELTA_ABS, disagree = round(1 - r, 5), stringsAsFactors = FALSE)
  }
  # THIS COSINE ROW IS A DISTRIBUTION OVER PATIENTS. MIMIC's cosine rows, in
  # `pairwise_metric_distributions.csv`, are a distribution over BAG PAIRS of
  # per-patient MEDIANS. The two share column names and summarise different
  # populations, so a reader lining them up would be comparing a between-patient
  # spread against a between-replicate spread. `attr_dist_summary()`'s own
  # header states the split; this note is here because the two tables are one
  # directory apart and nothing else says so. The transport section below joins
  # only the top-k ladder, which IS the same quantity at both sites.
  COS[[length(COS) + 1L]] <- cbind(
    data.frame(site = "eicu", method_a = p[1], method_b = p[2], agg = lv,
               population = "patients", stringsAsFactors = FALSE),
    attr_dist_summary(attr_cosine_dissim(Aa, Bb), above = 0.05))
  FLIP[[length(FLIP) + 1L]] <- cbind(
    data.frame(site = "eicu", method_a = p[1], method_b = p[2], agg = lv,
               stringsAsFactors = FALSE),
    attr_sign_flip(Aa, Bb, keep = kp, taus = TAUS))
  # THE SAME SCHEMA AS MIMIC's `leader_collapse_distributions`. An apply site
  # has no replicates, so the only distribution it can report for the rank
  # displacement and the evidence share is the one over PATIENTS on the single
  # frozen fits; that is the MIMIC table's `held = "boot=0"` stratum, and the
  # two sites read side by side on exactly that row set. Neither is comparable
  # with the between-replicate rows in `pairwise_metric_distributions`.
  COLL[[length(COLL) + 1L]] <- cbind(
    site = "eicu",
    data.frame(contrast = "L4", held = "boot=0", method_a = p[1], method_b = p[2],
               agg = lv, route = "frozen", n_pairs = 1L, masked = !is.null(kp),
               stringsAsFactors = FALSE),
    attr_leader_patient_summary(Aa, Bb, keep = kp, shares = COL_SH))
}
L4 <- do.call(rbind, L4); COS <- do.call(rbind, COS)
FLIP <- do.call(rbind, FLIP); COLL <- do.call(rbind, COLL)
save_table(run, L4, "l4_agreement_eicu", subdir = "diagnostics")
save_table(run, COS, "l4_cosine_eicu", subdir = "diagnostics")
save_table(run, FLIP, "l4_signflip_eicu", subdir = "diagnostics")
save_table(run, COLL, "leader_collapse_distributions", subdir = "diagnostics")

z <- L4[L4$agg == "signal" & L4$k == 1L & L4$delta == 0, ]
cat(sprintf("  %-16s %-16s %10s\n", "method a", "method b", "top1 disag"))
for (i in seq_len(nrow(z))) cat(sprintf("  %-16s %-16s %10.4f\n",
  z$method_a[i], z$method_b[i], z$disagree[i]))

# ============================================================================
# THE TRANSPORT READING: does the level-4 ladder have the same shape?
# ============================================================================
#
# DISTRIBUTIONAL, NOT PAIRED, and the header says why. What is compared is the
# LADDER: the same six pairs, the same metric, at two sites. A pair's value may
# move; what matters for the claim is whether the ORDERING survives, because the
# ordering is what the specification argument rests on.
cat("\n=== transport: the eICU ladder against the MIMIC one ===\n\n")
#' Establish that the MIMIC comparator describes THIS bundle's design.
#'
#' ADDED 2026-09-09 (review finding A10). The metrics run names its generator
#' store; the store carries the design key it was built under. Every field of
#' that key except `folds` is recomputed from the BUNDLE's frozen config and
#' compared -- `bam`, the signal vocabulary, the pairing, every fitted
#' formula, the cross terms and `k_ti`. The fold assignment is a training-site
#' object with no counterpart at an apply site, so it is carried through and
#' reported as `not_checkable_at_apply_site` rather than silently passed.
#' Matching METHOD NAMES never established a matching experiment.
transport_provenance <- function(mimic_run) {
  mm <- read_manifest(mimic_run)
  if (!identical(mm$status, "complete")) {
    stop("attr_external: the MIMIC comparator ", basename(mimic_run),
         " is not marked complete.", call. = FALSE)
  }
  gen <- mm$generator
  if (is.null(gen) || !nzchar(gen)) {
    stop("attr_external: ", basename(mimic_run), " names no generator store.",
         call. = FALSE)
  }
  gen_d <- file.path(dirname(mimic_run), gen)
  dgp <- file.path(gen_d, "design.qs2")
  if (!file.exists(dgp)) {
    stop("attr_external: the generator store ", gen, " behind ",
         basename(mimic_run), " has no design.qs2.", call. = FALSE)
  }
  dg <- qs2::qs_read(dgp)
  bk <- attr_design_key(cfg, fold_vec = NULL, k_ti = dg$k_ti,
                        folds_hash = dg$design_key$folds)
  d  <- attr_design_diff(dg$design_key, bk)
  if (nrow(d)) {
    print(d[, c("field", "hash_a", "hash_b")], row.names = FALSE)
    stop("attr_external: the MIMIC comparator ", basename(mimic_run),
         " was built under a different design from this bundle (", nrow(d),
         " field(s) differ, above). Its level-4 ladder is not the ladder of ",
         "the specifications this bundle froze.", call. = FALSE)
  }
  mm_methods <- strsplit(mm$methods %||% "", ",", fixed = TRUE)[[1]]
  data.frame(
    mimic_run = basename(mimic_run), generator = gen,
    generator_design_key = attr_key_hash(dg$design_key),
    bundle_design_key_excluding_folds = attr_key_hash(bk),
    fields_checked = paste(setdiff(names(bk), "folds"), collapse = ","),
    folds = "not_checkable_at_apply_site",
    store_fingerprint = if (is.null(dg$fingerprint)) "absent" else attr_key_hash(dg$fingerprint),
    mimic_store_status = mm$store_status %||% "unrecorded",
    mimic_bag_population = mm$common_bags %||% "unrecorded",
    methods_match = setequal(mm_methods, .mw),
    mimic_basis = "out-of-fold 5-fold fits on MIMIC train (OOF priors, OOF GAMs, OOF boosters)",
    eicu_basis  = "single frozen full-train bundle applied once; final priors, final GAMs, final booster",
    levels_present_eicu = "L1,L4",
    levels_absent_eicu  = "L2,L3,L3T,L3P",
    stringsAsFactors = FALSE)
}
if (is.null(MIMIC_RUN) || !dir.exists(MIMIC_RUN)) {
  cat("  no complete MIMIC attrmetrics run given; eICU tables written, comparison skipped.\n")
} else {
  PROV <- transport_provenance(MIMIC_RUN)
  save_table(run, PROV, "transport_provenance", subdir = "diagnostics")
  cat(sprintf("  comparator %s (generator %s) matches the bundle design on every\n",
              PROV$mimic_run, PROV$generator))
  cat("  checkable field; the fold assignment cannot be checked at an apply site.\n")
  if (!PROV$methods_match) {
    cat("  *** the comparator's method set differs from this run's; the join is\n")
    cat("      over the intersection only.\n")
  }
  mp <- file.path(MIMIC_RUN, "diagnostics", "disagreement_distributions.rds")
  if (!file.exists(mp)) {
    cat("  ", basename(MIMIC_RUN), " carries no disagreement_distributions.\n", sep = "")
  } else {
    md <- readRDS(mp)
    m4 <- md[md$contrast == "L4" & md$agg == "signal" & md$k == 1L, ]
    if ("bag_population" %in% names(m4)) m4 <- m4[m4$bag_population == "common", ]
    e4 <- L4[L4$agg == "signal" & L4$k == 1L & L4$delta == 0, ]
    key <- function(a, b) paste(pmin(a, b), pmax(a, b))
    m4$kk <- key(m4$method_a, m4$method_b); e4$kk <- key(e4$method_a, e4$method_b)
    j <- merge(e4[, c("kk", "method_a", "method_b", "disagree")],
               m4[, c("kk", "anchor_boot0", "median", "p95")], by = "kk")
    if (!nrow(j)) {
      cat("  no shared pairs between the two runs.\n")
    } else {
      names(j)[names(j) == "disagree"] <- "eicu"
      names(j)[names(j) == "anchor_boot0"] <- "mimic_anchor"
      j$shift <- round(j$eicu - j$mimic_anchor, 5)
      # IS THE eICU VALUE INSIDE MIMIC'S OWN LEVEL-4 SPREAD? That spread is the
      # bootstrap distribution over shared bags, so a shift inside it is not
      # distinguishable from MIMIC's own resampling variability.
      j$within_mimic_p95 <- j$eicu <= j$p95
      # THE TWO BASES DIFFER IN TRAINING REGIME AS WELL AS IN SITE (finding
      # A10): MIMIC's numbers are out-of-fold five-fold fits, eICU's are one
      # frozen full-train fit. A shift between them carries both.
      j$mimic_basis <- "oof_5fold_train"; j$eicu_basis <- "final_full_train_frozen"
      j$mimic_run <- basename(MIMIC_RUN)
      j <- j[order(j$mimic_anchor), ]
      save_table(run, j, "l4_transport", subdir = "diagnostics")
      cat(sprintf("  %-16s %-16s %8s %8s %8s %s\n", "method a", "method b",
                  "MIMIC", "eICU", "shift", "in MIMIC p95"))
      for (i in seq_len(nrow(j))) cat(sprintf("  %-16s %-16s %8.4f %8.4f %+8.4f %s\n",
        j$method_a[i], j$method_b[i], j$mimic_anchor[i], j$eicu[i], j$shift[i],
        if (isTRUE(j$within_mimic_p95[i])) "yes" else "NO"))
      rho <- if (nrow(j) > 2) stats::cor(j$mimic_anchor, j$eicu, method = "spearman") else NA_real_
      cat(sprintf("\n  rank correlation of the ladder across sites: %.4f (n = %d pairs)\n",
                  rho, nrow(j)))
      cat("  The ORDERING is what the specification argument rests on; a shift in\n")
      cat("  level that preserves the ordering is a weaker finding than a reordering.\n")
      cat("  MIMIC values are out-of-fold five-fold fits; eICU values are one frozen\n")
      cat("  full-train fit, so a shift carries a training-regime change as well as a site change.\n")
    }
  }
}

save_table(run, ATTR_ESTIMAND_NOTES, "design_notes", subdir = "diagnostics")
finalize_run(run, extra = list(
  bundle = basename(dirname(bpath)),
  arms = paste(names(A), collapse = ","),
  n_stays = length(stay_ids),
  levels_present = "L1,L4",
  topk_definition = "tie_tolerant_at_delta_0 (attr_agree_k semantics; identical to the MIMIC fast path as of 2026-09-09)",
  mimic_run = if (is.null(MIMIC_RUN)) "" else basename(MIMIC_RUN),
  # `arms_absent` IS GONE. It named the four interaction arms until 2026-09-07,
  # when they entered the bundle; an "absent arms" field that is empty is a
  # field that will be read as "nothing was checked" the next time an arm goes
  # missing, and the guard in the arm loop above is the real check now.
  n_arms = length(A),
  levels_absent = "L2,L3,L3T (an apply site fits nothing)"))
cat(sprintf("\nwritten: %s\n", run$path))
