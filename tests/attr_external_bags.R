# tests/attr_external_bags.R ---------------------------------------------------
# THE eICU REPRODUCIBILITY ARM: GENERATOR. Added 2026-10-06
# (paper/plans/plan_eicu_reproducibility_arm.md).
#
# Retrains the paper's three attribution methods -- the measurement and paired
# weights of evidence (`llr_meas`, `llr_full`) and SHAP of the constructed
# booster (`shap_xgb_feat`) -- WHOLE on the 38 shared bags of the internal run,
# and applies every retrained model to the eICU cohort. The stored object is
# the same n x 19 contribution matrix the internal store holds, one row per
# eICU stay, at the same replicate coordinates, so `tests/attr_metrics.R
# --site eicu` forms the same contrasts at eICU: L1, L2, L3, L3T, L3P and L4.
#
# WHY THIS EXISTS. `tests/attr_external.R` attributes eICU once, with the frozen
# bundle, and its header argued that levels 2, 3 and 3T "do not exist at an
# apply site". That is true of the frozen bundle only. Retraining on the
# training bags and APPLYING each retrained model at eICU measures directly
# whether a retrained model names the same leader for an eICU patient.
#
# THE DESIGN, AND WHAT IT HOLDS FIXED.
#   bags        bag b = bag_of(bootstrap.seed_base + b), b = 1..40, the shared
#               manifest of the internal run; each bag's membership hash is
#               checked against the internal store's `bootstrap_bags` table.
#               The LLR arms skip the bags the internal run tombstoned (20 and
#               26), so the LLR bag set is the internal one; SHAP is refitted on
#               all 40, as internally. The common population is 38 bags.
#   fit         WHOLE BAG, no folds: every GAM and the booster train on the
#               bag's training stays, about 63% of them. An internal bag refit
#               trained on the bag intersected with four folds, about 50%. Both
#               families are refitted on the same whole bags, so a comparison
#               BETWEEN them at eICU is like for like; but every eICU floor is
#               the floor of whole-bag fits, and its absolute level is not the
#               level fold-intersection fits would give. Set an eICU floor
#               beside a MIMIC-IV floor with that difference stated.
#   acceptance  each GAM is fitted through `fit_one()`, the path of the final
#               fits being perturbed, which applies `.accept_outcome_fit()`
#               (finite coefficients, convergence). The internal bag refits
#               called `.bam_fit()` directly, without acceptance. So a fit that
#               fails acceptance here excludes its bag, with the cause in the
#               tombstone, where internally it would have entered.
#   priors      the FINAL bundle priors (alpha, delta and lambda coefficients,
#               p_bar), held fixed, as the internal bags hold the out-of-fold
#               priors fixed.
#   seeds       a bag booster uses cfg$seed (+ levels.seed.seed_base * s when
#               re-seeded); the seed route refits the booster on the WHOLE
#               training set with cfg$seed + seed_base * s, s = 1..8, beside
#               the anchor, the bundle's own booster (cfg$seed).
#   posterior   40 coefficient draws from each FINAL bundle GAM, smoothing-
#               corrected covariance, seeded draw_seed_base + the signal's
#               index (internally draw_seed_base + a (signal, fold) counter:
#               a different stream of the same construction); hard rule 9
#               calls this application, not refitting.
#   anchor      coordinate (0, 0, 0): the frozen bundle applied, as
#               `tests/attr_external.R` does.
#
# HARD RULES. Rule 1: output is counts, hashes, timings and summaries only.
# Rule 7: no GAM or booster is persisted; each is applied in memory and only
# the contribution matrix is written. Rule 9: training-data perturbation fits
# for reproducibility only. No file under R/ is added or changed, so no target
# is invalidated.
#
# MEMORY. Measured 2026-10-06: the tables and the bundle hold about 1 GB, a
# bag-size GAM 4 MB. R's maximum is logged after every step, and the run stops
# itself above MEM_LIMIT_GB (resumable).
#
#   Rscript tests/attr_external_bags.R --plan-only
#   Rscript tests/attr_external_bags.R --verify-anchor --bags 1     # benchmark
#   Rscript tests/attr_external_bags.R --verify-anchor --resume out/runs/attrextgen_...
#   options: --routes anchor,seed,bootstrap,posterior   --bags 1,2,3
#            --internal out/runs/attrgen_20260909T115047
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(targets); library(mgcv); library(xgboost); library(arrow); library(yaml)
  library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/attr_external_common.R")

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
PLAN_ONLY     <- "--plan-only" %in% args
VERIFY_ANCHOR <- "--verify-anchor" %in% args
RESUME        <- .opt("--resume", NULL)
INTERNAL      <- .opt("--internal", "out/runs/attrgen_20260909T115047")
ROUTES_ALL    <- c("anchor", "seed", "bootstrap", "posterior")
ROUTES        <- strsplit(.opt("--routes", paste(ROUTES_ALL, collapse = ",")), ",",
                          fixed = TRUE)[[1]]
if (length(setdiff(ROUTES, ROUTES_ALL))) {
  abort_values("--routes must be a subset of anchor,seed,bootstrap,posterior",
               setdiff(ROUTES, ROUTES_ALL))
}
BAGS_ARG <- .opt("--bags", NULL)

# Tolerances, declared here and used below: the refit-against-bundle check is
# exact in practice (measured 0) and allowed 1e-8; SHAP additivity is float32
# arithmetic inside xgboost (measured about 1e-5) and allowed 1e-3.
TOL_REFIT    <- 1e-8
TOL_ADDITIVE <- 1e-3
MEM_LIMIT_GB <- 12

invisible(gc(reset = TRUE))
#' R's maximum memory since the start, in GB (gc() "max used", both cell types).
mem_max_gb <- function() { g <- gc(verbose = FALSE); sum(g[, ncol(g)]) / 1024 }
mem_check <- function(where) {
  mx <- mem_max_gb()
  cat(sprintf("  [memory] %-28s R max used %.2f GB\n", where, mx))
  if (mx > MEM_LIMIT_GB) {
    stop("R's maximum memory passed ", MEM_LIMIT_GB, " GB at ", where,
         ". Stopped before going further; the store is resumable.", call. = FALSE)
  }
  invisible(mx)
}

t_all <- start_timer()
ecfg <- yaml::read_yaml("config/attribution_eval.yml")
xcfg <- yaml::read_yaml("config/external.yml")
METHODS <- EICU_METHODS
attr_check_methods(METHODS)

# --- the frozen bundle, and the config it carries ----------------------------
bpath  <- as.character(cfg_req(xcfg, "bundle"))
bundle <- load_bundle(bpath, verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
sigs   <- as.character(unlist(cfg$signals))
cat(sprintf("\n=== eICU reproducibility arm: whole-bag refits, bundle %s ===\n\n",
            basename(dirname(bpath))))

# --- the training side, from the targets store (read only) -------------------
cfg_t    <- tar_read(cfg)
folds    <- tar_read(folds)
tr       <- tar_read(train_ids)
y        <- as.integer(tar_read(y_train))
priors_t <- tar_read(priors)
tabs_m   <- tar_read(tabs)
stopifnot(length(y) == length(tr))

# GUARD 1: the targets store and the bundle describe one design and one set of
# final priors. A refit that used the store's config and the bundle's priors
# under two different designs would be a third model.
.dk <- intersect(BUNDLE_DESIGN_KEYS, names(cfg_t))
if (!identical(attr_key_hash(cfg_t[.dk]), attr_key_hash(bundle$cfg))) {
  stop("the targets store's config and the bundle's frozen design differ; the ",
       "bundle was not built from the store this script reads.", call. = FALSE)
}
if (!identical(attr_key_hash(priors_final(priors_t)), attr_key_hash(bundle$priors))) {
  stop("the targets store's final priors differ from the bundle's.", call. = FALSE)
}
cat("  guard 1: targets store and bundle share the design and the final priors\n")

# GUARD 2: the internal store is the one the bags and the tombstones come from,
# and it describes this bundle's design (the check `attr_external.R`'s
# `transport_provenance()` makes, fields other than the folds).
.im <- read_manifest(INTERNAL)
if (!identical(.im$status, "complete")) {
  stop("the internal store ", basename(INTERNAL), " is not marked complete.", call. = FALSE)
}
int_dg <- qs2::qs_read(file.path(INTERNAL, "design.qs2"))
.bk <- attr_design_key(cfg, fold_vec = NULL, k_ti = int_dg$k_ti,
                       folds_hash = int_dg$design_key$folds)
.dd <- attr_design_diff(int_dg$design_key, .bk)
if (nrow(.dd)) {
  print(.dd[, c("field", "hash_a", "hash_b")], row.names = FALSE)
  stop("the internal store was built under a different design from this bundle.",
       call. = FALSE)
}
if (!identical(as.character(int_dg$stay_id), as.character(tr))) {
  stop("the internal store's training rows are not the targets store's.", call. = FALSE)
}
cat(sprintf("  guard 2: internal store %s matches the bundle design\n", basename(INTERNAL)))

# --- the eICU side ------------------------------------------------------------
tabs_e   <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
stay_ids <- tabs_e$cohort$stay_id
ids_ch   <- as.character(stay_ids)
n_e      <- length(stay_ids)
meas_e   <- measured_matrix(tabs_e, cfg, stay_ids)
cat(sprintf("  eICU stays: %d; training stays: %d\n", n_e, length(tr)))

# --- the bags: the internal run's shared manifest -----------------------------
BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
SEED_B    <- as.integer(cfg_req(ecfg, "levels", "seed", "seed_base"))
B_SEED    <- as.integer(cfg_req(ecfg, "levels", "seed", "b"))
B_PER_BAG <- as.integer(cfg_req(ecfg, "levels", "seed", "b_per_bag"))
B_DRAW    <- as.integer(cfg_req(ecfg, "levels", "sample", "b"))
DRAW0     <- as.integer(cfg_req(ecfg, "levels", "sample", "draw_seed_base"))
if (B_PER_BAG != 1L) abort_values("this generator implements levels.seed.b_per_bag = 1", B_PER_BAG)

#' Copied verbatim from `tests/attr_replicates.R` (L742-757), which defines it
#' inside its own top-level code and cannot be sourced. Identity is not
#' assumed: every bag's membership hash is checked against the internal
#' store's table below.
bag_of <- function(seed) {
  gcol <- resample_cols(tabs_m$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs_m$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  if (anyNA(gid)) stop("bag_of: `", gcol, "` is missing for some training stay",
                       call. = FALSE)
  ug   <- unique(gid)
  with_seed(seed, {
    keep <- sample(ug, length(ug), replace = TRUE)
    gid %in% unique(keep)
  })
}
BAG <- lapply(seq_len(B_BOOT), function(b) bag_of(BAG_SEED0 + b))
BAGT <- do.call(rbind, lapply(seq_len(B_BOOT), function(b) data.frame(
  boot_id = b, seed = BAG_SEED0 + b, n_in_bag = sum(BAG[[b]]), n_train = length(tr),
  frac_in_bag = round(mean(BAG[[b]]), 5),
  membership_hash = attr_key_hash(which(BAG[[b]])), stringsAsFactors = FALSE)))
.ibt <- readRDS(file.path(INTERNAL, "diagnostics", "bootstrap_bags.rds"))
.ibt <- .ibt[match(BAGT$boot_id, .ibt$boot_id), ]
if (!identical(.ibt$membership_hash, BAGT$membership_hash)) {
  stop("the bags drawn here differ from the internal store's bags (membership ",
       "hashes differ for boot_id ",
       paste(BAGT$boot_id[.ibt$membership_hash != BAGT$membership_hash], collapse = ","),
       ").", call. = FALSE)
}
cat(sprintf("  bags: %d drawn, every membership hash equals the internal store's\n", B_BOOT))

# The bags the internal run excluded for the LLR arms, read from its tombstone
# rather than written down here.
.it <- utils::read.csv(file.path(INTERNAL, "llr_bootstrap_excluded.csv"), stringsAsFactors = FALSE)
.it <- .it[.it$design %in% c(attr_key_hash(int_dg$design_key), int_dg$migrated_from), , drop = FALSE]
.it <- .it[.it$membership_hash == BAGT$membership_hash[match(.it$boot_id, BAGT$boot_id)], , drop = FALSE]
INT_EXCL <- sort(unique(as.integer(.it$boot_id)))
cat(sprintf("  LLR bags excluded in the internal run, excluded here too: %s\n",
            paste(INT_EXCL, collapse = ", ")))

BAGS_RUN <- if (is.null(BAGS_ARG)) seq_len(B_BOOT) else
  as.integer(strsplit(BAGS_ARG, ",", fixed = TRUE)[[1]])
if (any(!BAGS_RUN %in% seq_len(B_BOOT))) abort_values("--bags outside 1..B", BAGS_RUN)

# --- identity, plan ------------------------------------------------------------
ID   <- eicu_store_identity(int_dg$design_key, bundle, cfg, ecfg, stay_ids, METHODS)
EDES <- ID$design_key
.DK  <- attr_key_hash(EDES)
PLAN <- eicu_store_plan(ecfg, n_folds = as.integer(cfg_req(cfg, "n_folds")), METHODS)
PLAN <- PLAN[!is.na(PLAN$index) &
               PLAN$route %in% c("ladder", "seed", "bootstrap", "bootstrap_seeded", "posterior"), ]
cat(sprintf("\n  replicate plan (planned; LLR bags %s will be tombstoned, so %d LLR bags are fitted):\n",
            paste(INT_EXCL, collapse = ", "), B_BOOT - length(INT_EXCL)))
print(table(PLAN$method, PLAN$route))
if (PLAN_ONLY) {
  cat(sprintf("\n  design key %s; plan-only, nothing fitted or written.\n", .DK))
  quit(save = "no", status = 0L)
}

# --- the run directory -----------------------------------------------------------
if (!is.null(RESUME)) {
  if (!dir.exists(RESUME)) stop("--resume: no such directory: ", RESUME, call. = FALSE)
  run <- structure(list(prefix = "attrextgen", id = basename(RESUME), path = RESUME,
                        started = Sys.time(), config = cfg,
                        log_file = file.path(RESUME, "log.txt")), class = "llr_run")
  cat(sprintf("\n  resuming %s\n", basename(RESUME)))
} else {
  run <- new_run("attrextgen", cfg, note = sprintf(
    "eICU reproducibility: whole-bag refits on the internal bags, applied to eICU; design %s",
    .DK))
}
SUB <- cfg_req(ecfg, "storage", "subdir")
dir.create(file.path(run$path, SUB), showWarnings = FALSE, recursive = TRUE)
DESIGN_P <- file.path(run$path, "design.qs2")
MEAS_P   <- file.path(run$path, "measured.qs2")
if (file.exists(DESIGN_P)) {
  old <- qs2::qs_read(DESIGN_P)
  if (!identical(old$design_key, EDES)) {
    stop("--resume: the store's design key differs from the live one; this is a ",
         "different experiment, not a resume.", call. = FALSE)
  }
  .fd <- attr_fingerprint_diff(old$fingerprint, ID$fingerprint)
  if (nrow(.fd)) {
    print(.fd[, intersect(c("field", "hash_a", "hash_b"), names(.fd))], row.names = FALSE)
    stop("--resume: the store's fingerprint differs from the live design.", call. = FALSE)
  }
  if (!identical(old$stay_id, stay_ids)) stop("--resume: eICU rows differ.", call. = FALSE)
  if (!identical(qs2::qs_read(MEAS_P), meas_e)) stop("--resume: measured mask differs.", call. = FALSE)
  cat("  resume: design key, fingerprint, rows and mask all match the store\n")
} else {
  qs2::qs_save(list(site = "eicu", design_key = EDES, fingerprint = ID$fingerprint,
                    stay_id = stay_ids, methods = METHODS, generator = EICU_GEN_VERSION,
                    bundle = bpath, internal_store = basename(INTERNAL),
                    internal_design_key = int_dg$design_key, k_ti = int_dg$k_ti,
                    internal_excluded_bags = INT_EXCL),
               DESIGN_P)
  qs2::qs_save(meas_e, MEAS_P)
}
save_table(run, BAGT, "bootstrap_bags", subdir = "diagnostics")
save_table(run, PLAN, "replicate_plan", subdir = "diagnostics")

# The tombstone, in the internal store's format, so the consumer reads it the
# same way. The internal exclusions are carried over with their cause.
TOMB_P <- file.path(run$path, "llr_bootstrap_excluded.csv")
TOMB <- if (file.exists(TOMB_P)) utils::read.csv(TOMB_P, stringsAsFactors = FALSE) else
  data.frame(boot_id = integer(0), design = character(0), bag_seed = integer(0),
             membership_hash = character(0), n_spec_folds = integer(0),
             cause = character(0), stringsAsFactors = FALSE)
tombstone_bag <- function(b, cause) {
  if (b %in% TOMB$boot_id) return(invisible(FALSE))
  TOMB <<- rbind(TOMB, data.frame(boot_id = b, design = .DK, bag_seed = BAG_SEED0 + b,
                                  membership_hash = BAGT$membership_hash[b],
                                  n_spec_folds = NA_integer_, cause = cause,
                                  stringsAsFactors = FALSE))
  utils::write.csv(TOMB, TOMB_P, row.names = FALSE)
  invisible(TRUE)
}
for (b in INT_EXCL) tombstone_bag(b, paste0("excluded_in_internal_run:",
                                            .it$cause[match(b, .it$boot_id)]))

# --- keys and the manifest (the encoding of tests/attr_replicates.R) ------------
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
        list("sample", 500L + b, list(route = "bootstrap", bootstrap_seed = BAG_SEED0 + b)),
    bootstrap_seeded = list("sample", 1000L * sd + b,
        list(bootstrap_seed = BAG_SEED0 + b, seed_offset = SEED_B * sd)),
    posterior = list("sample", dw, list(route = "posterior")),
    abort_values(".coord_key: no key encoding for route", route))
  list(stage = z[[1]], index = z[[2]],
       key = attr_replicate_key(method, z[[1]], z[[2]], EDES, z[[3]]))
}
.path_of <- function(k) file.path(run$path, SUB, paste0(k$key, ".qs2"))
have <- function(method, route, boot_id = 0L, seed_id = 0L, draw_id = 0L)
  file.exists(.path_of(.coord_key(method, route, boot_id, seed_id, draw_id)))

MAN_COLS <- c("method", "route", "stage", "index", "boot_id", "seed_id", "draw_id",
              "key", "n_rows", "n_cols")
MAN_P <- file.path(run$path, "manifest_replicates.csv")
MAN <- if (file.exists(MAN_P)) utils::read.csv(MAN_P, stringsAsFactors = FALSE)[, MAN_COLS] else
  data.frame(method = character(0), route = character(0), stage = character(0),
             index = integer(0), boot_id = integer(0), seed_id = integer(0),
             draw_id = integer(0), key = character(0), n_rows = integer(0),
             n_cols = integer(0), stringsAsFactors = FALSE)

#' The manifest as a function of the PLAN and the files on disk (the rule of
#' `reconcile_manifest()` in tests/attr_replicates.R L1064-1092). A run
#' interrupted between writing a replicate and writing its manifest row leaves
#' a file the manifest does not name; `have()` then skips it on resume and an
#' incremental manifest would never list it. Rebuilt at the start of every
#' invocation and at the close.
reconcile_manifest <- function() {
  rows <- list()
  for (i in seq_len(nrow(PLAN))) {
    k <- .coord_key(PLAN$method[i], PLAN$route[i], PLAN$boot_id[i], PLAN$seed_id[i],
                    PLAN$draw_id[i])
    if (!file.exists(.path_of(k))) next
    old <- MAN[MAN$key == k$key, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      method = PLAN$method[i], route = PLAN$route[i], stage = k$stage,
      index = as.integer(k$index), boot_id = as.integer(PLAN$boot_id[i]),
      seed_id = as.integer(PLAN$seed_id[i]), draw_id = as.integer(PLAN$draw_id[i]),
      key = k$key,
      n_rows = if (nrow(old)) old$n_rows[1] else NA_integer_,
      n_cols = if (nrow(old)) old$n_cols[1] else NA_integer_,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(MAN[0, , drop = FALSE])
  out <- do.call(rbind, rows)
  for (i in which(is.na(out$n_rows))) {
    M <- qs2::qs_read(file.path(run$path, SUB, paste0(out$key[i], ".qs2")))
    attr_check_replicate(M, ids_ch, sigs, out$key[i])
    out$n_rows[i] <- nrow(M); out$n_cols[i] <- ncol(M)
  }
  out
}
.n_before <- nrow(MAN)
MAN <- reconcile_manifest()
utils::write.csv(MAN, MAN_P, row.names = FALSE)
if (nrow(MAN) != .n_before) {
  cat(sprintf("  manifest reconciled from the plan and the disk: %d -> %d rows\n",
              .n_before, nrow(MAN)))
}

#' Write a replicate, or verify it bitwise if it is already on disk.
emit <- function(M, method, route, boot_id = 0L, seed_id = 0L, draw_id = 0L) {
  attr_check_replicate(M, ids_ch, sigs, paste(method, route, boot_id, seed_id, draw_id))
  k <- .coord_key(method, route, boot_id, seed_id, draw_id)
  if (file.exists(.path_of(k))) {
    old <- qs2::qs_read(.path_of(k))
    if (!identical(dimnames(old), dimnames(M)) || max(abs(old - M)) > 0) {
      stop("emit: replicate ", k$key, " exists with DIFFERENT content.", call. = FALSE)
    }
  } else {
    qs2::qs_save(M, .path_of(k))
  }
  row <- data.frame(method = method, route = route, stage = k$stage,
                    index = as.integer(k$index), boot_id = as.integer(boot_id),
                    seed_id = as.integer(seed_id), draw_id = as.integer(draw_id),
                    key = k$key, n_rows = nrow(M), n_cols = ncol(M),
                    stringsAsFactors = FALSE)
  MAN <<- rbind(MAN[MAN$key != k$key, , drop = FALSE], row)
  utils::write.csv(MAN, MAN_P, row.names = FALSE)
  invisible(k$key)
}

# --- shared ingredients, built once ----------------------------------------------
PRI <- lapply(stats::setNames(sigs, sigs), function(sg) priors_for(bundle$priors, sg, "final", NA_integer_))

# The evidence-model specifications: the measurement model of every signal, and
# the paired model where `spec_source()` says it is fitted (an unpaired signal's
# paired column is its measurement column, the alias rule). The anchor check
# below proves this set and this alias rule reproduce `apply_bundle()` bitwise.
SPECS <- do.call(rbind, lapply(sigs, function(sg) {
  src_full <- spec_source(sg, "full", cfg)
  data.frame(signal = sg,
             model = if (is.na(src_full)) c("meas", "full") else "meas",
             stringsAsFactors = FALSE)
}))
SPECS$key <- paste0(SPECS$signal, "/", SPECS$model)
FULL_ALIAS <- setdiff(sigs, SPECS$signal[SPECS$model == "full"])
cat(sprintf("\n  evidence-model specs per bag: %d (paired models %d; paired column aliased to measurement for %d signals)\n",
            nrow(SPECS), sum(SPECS$model == "full"), length(FULL_ALIAS)))

# eICU frames, one per spec, built once: the priors are frozen, so the frame a
# refitted model predicts on is the same for every bag. Mirrors `apply_one()`.
ND <- ROWS <- list()
for (i in seq_len(nrow(SPECS))) {
  sg <- SPECS$signal[i]; md <- SPECS$model[i]
  nd <- signal_frame(sg, md, tabs_e, cfg, PRI[[sg]], stay_ids = stay_ids, stage = "predict")
  ND[[SPECS$key[i]]]   <- nd
  ROWS[[SPECS$key[i]]] <- match(as.character(nd$stay_id), ids_ch)
}

#' `apply_one()` on a cached frame: the same predict call and the same
#' subtraction, without rebuilding the frame per bag.
apply_cached <- function(b, key, sg) {
  eta <- as.numeric(stats::predict(b, newdata = ND[[key]], type = "link", discrete = FALSE))
  if (anyNA(eta)) stop("apply_cached [", key, "]: NA predictions (ids not printed).", call. = FALSE)
  eta - logit(PRI[[sg]]$p_bar)
}

#' One empty eICU contribution matrix.
.empty_M <- function() matrix(0, n_e, length(sigs), dimnames = list(ids_ch, sigs))

# The booster designs. The training design is the one the pipeline's final
# booster was fitted on; the eICU design is aligned to its feature names.
Xtr <- xgb_design_feat(tabs_m, cfg, bundle$priors, tr, role = "final")
FN  <- bundle$xgb$xgb_feat$feature_names
if (!identical(colnames(Xtr), FN)) {
  stop("the training booster design's columns differ from the bundle booster's.", call. = FALSE)
}
X_e <- align_design(xgb_design_feat(tabs_e, cfg, bundle$priors, stay_ids, role = "final",
                                    fold = NA_integer_, feature_names = FN), FN)
DM_e <- xgboost::xgb.DMatrix(X_e, missing = NA)
grp  <- patient_group_of(tabs_m$cohort, cfg, tr)
GRP_OF <- sub("__.*$", "", FN)

#' SHAP of a booster on eICU, grouped to the 19 signals (bias and intervention
#' groups dropped), exactly as `tests/attr_external.R` L186-194. Additivity --
#' the full contribution row summing to the margin -- is checked, since a
#' failure would mean the contribution path and the prediction path disagree.
shap_eicu <- function(booster, what) {
  ctr <- stats::predict(booster, DM_e, predcontrib = TRUE)
  stopifnot(ncol(ctr) == length(FN) + 1L)
  eta <- stats::predict(booster, DM_e, outputmargin = TRUE)
  addit <- max(abs(rowSums(ctr) - eta))
  if (!is.finite(addit) || addit > TOL_ADDITIVE) {
    stop("SHAP additivity failed for ", what, ": max |rowSums(contrib) - margin| = ",
         format(addit), " > ", TOL_ADDITIVE, ".", call. = FALSE)
  }
  S <- ctr[, seq_along(FN), drop = FALSE]
  G <- .empty_M()
  for (g in sigs) { j <- which(GRP_OF == g); if (length(j)) G[, g] <- rowSums(S[, j, drop = FALSE]) }
  list(G = G, additivity_max_abs = addit)
}
fit_booster <- function(use, seed) {
  xgb_fit_full(Xtr[use, , drop = FALSE], y[use], cfg, seed = seed, group = grp[use])
}

DIAG <- list()
.diag <- function(route, method, boot_id = 0L, seed_id = 0L, additivity_max_abs = NA_real_,
                  max_abs_diff_vs_bundle = NA_real_, best_iter = NA_integer_,
                  elapsed_sec = NA_real_) {
  DIAG[[length(DIAG) + 1L]] <<- data.frame(
    route = route, method = method, boot_id = as.integer(boot_id), seed_id = as.integer(seed_id),
    additivity_max_abs = additivity_max_abs, max_abs_diff_vs_bundle = max_abs_diff_vs_bundle,
    best_iter = as.integer(best_iter), elapsed_sec = round(elapsed_sec, 2),
    r_max_used_gb = round(mem_max_gb(), 3), stringsAsFactors = FALSE)
}
mem_check("inputs loaded")

# =================================================================================
# ROUTE anchor: the frozen bundle applied, coordinate (0, 0, 0)
# =================================================================================
if ("anchor" %in% ROUTES) {
  cat("\n=== anchor: the frozen bundle applied at eICU ===\n")
  t0 <- start_timer()
  ap <- apply_bundle(bundle, tabs_e, cfg, stay_ids, arms = names(LLR_ARM_MATRIX), verbose = FALSE)
  A0 <- list(llr_meas = ap$l_mats$meas[, sigs, drop = FALSE],
             llr_full = ap$l_mats$full[, sigs, drop = FALSE])
  rm(ap)
  # The assembly every bag replicate goes through must reproduce `apply_bundle()`
  # on the bundle's own models, or a bag replicate would differ from the anchor
  # for a reason that is not the bag.
  Mm <- .empty_M(); Mf <- .empty_M()
  for (i in seq_len(nrow(SPECS))) {
    sg <- SPECS$signal[i]; key <- SPECS$key[i]
    v <- apply_cached(bundle$models[[key]], key, sg)
    if (SPECS$model[i] == "meas") Mm[ROWS[[key]], sg] <- v else Mf[ROWS[[key]], sg] <- v
  }
  for (sg in FULL_ALIAS) Mf[, sg] <- Mm[, sg]
  A0c <- list(llr_meas = Mm, llr_full = Mf); rm(Mm, Mf)
  for (m in names(A0)) {
    d <- max(abs(A0c[[m]] - A0[[m]]))
    cat(sprintf("  assembly check %-9s max|cached assembly - apply_bundle| = %.3e\n", m, d))
    if (!identical(rownames(A0[[m]]), ids_ch) || d > 0) {
      stop("the cached assembly does not reproduce apply_bundle() for ", m, ".", call. = FALSE)
    }
    emit(A0[[m]], m, "ladder")
    .diag("ladder", m, max_abs_diff_vs_bundle = d)
  }
  rm(A0, A0c)
  sh0 <- shap_eicu(bundle$xgb$xgb_feat$booster, "the bundle booster")
  emit(sh0$G, "shap_xgb_feat", "ladder")
  .diag("ladder", "shap_xgb_feat", additivity_max_abs = sh0$additivity_max_abs,
        best_iter = bundle$xgb$xgb_feat$best_iter, elapsed_sec = t0()$elapsed_sec)
  rm(sh0)
  cat(sprintf("  anchor written (%.1f min)\n", t0()$elapsed_sec / 60))
  mem_check("anchor")
}

# The anchor reproduction check: refit EVERY spec, and the booster, on the WHOLE
# training set through the path every bag uses, and compare with the bundle's
# own models on eICU.
if (VERIFY_ANCHOR) {
  cat(sprintf("\n=== anchor reproduction: %d whole-training-set refits and the booster against the bundle ===\n",
              nrow(SPECS)))
  for (i in seq_len(nrow(SPECS))) {
    sg <- SPECS$signal[i]; md <- SPECS$model[i]; key <- SPECS$key[i]
    t0 <- start_timer()
    r <- fit_one(sg, md, tabs_m, cfg, PRI[[sg]], fit_ids = tr, keep_model = TRUE, role = "final")
    d <- max(abs(apply_cached(r$model, key, sg) - apply_cached(bundle$models[[key]], key, sg)))
    rm(r)
    cat(sprintf("  %-26s max|refit - bundle| on eICU = %.3e  (%.1f s)\n", key, d, t0()$elapsed_sec))
    .diag("verify_anchor", key, max_abs_diff_vs_bundle = d, elapsed_sec = t0()$elapsed_sec)
    if (!is.finite(d) || d > TOL_REFIT) {
      stop("the whole-set refit does not reproduce the bundle model ", key,
           " (max|d| = ", format(d), "). Identify the branch before going on.", call. = FALSE)
    }
  }
  t0 <- start_timer()
  bst <- fit_booster(rep(TRUE, length(tr)), cfg$seed)
  d <- max(abs(shap_eicu(bst$booster, "the whole-set booster refit")$G -
                 shap_eicu(bundle$xgb$xgb_feat$booster, "the bundle booster")$G))
  cat(sprintf("  booster refit (seed %d) vs bundle booster, eICU SHAP: max|d| = %.3e  best_iter %d vs %d  (%.1f s)\n",
              as.integer(cfg$seed), d, bst$best_iter, bundle$xgb$xgb_feat$best_iter, t0()$elapsed_sec))
  .diag("verify_anchor", "shap_xgb_feat", max_abs_diff_vs_bundle = d, best_iter = bst$best_iter,
        elapsed_sec = t0()$elapsed_sec)
  if (!is.finite(d) || d > TOL_REFIT) {
    stop("the whole-set booster refit does not reproduce the bundle booster ",
         "(max|d| = ", format(d), ").", call. = FALSE)
  }
  rm(bst)
  mem_check("anchor reproduction")
}

# =================================================================================
# ROUTE seed: the booster re-seeded on the WHOLE training set, (0, s, 0)
# =================================================================================
if ("seed" %in% ROUTES) {
  cat("\n=== seed: booster refitted on the whole training set, re-seeded ===\n")
  for (s in seq_len(B_SEED)) {
    if (have("shap_xgb_feat", "seed", seed_id = s)) next
    t0 <- start_timer()
    bst <- fit_booster(rep(TRUE, length(tr)), cfg$seed + SEED_B * s)
    sh <- shap_eicu(bst$booster, sprintf("seed %d", s))
    emit(sh$G, "shap_xgb_feat", "seed", seed_id = s)
    .diag("seed", "shap_xgb_feat", seed_id = s, additivity_max_abs = sh$additivity_max_abs,
          best_iter = bst$best_iter, elapsed_sec = t0()$elapsed_sec)
    cat(sprintf("  seed %d written (%.1f min)\n", s, t0()$elapsed_sec / 60))
    rm(bst, sh)
  }
  mem_check("seed route")
}

# =================================================================================
# ROUTE bootstrap (and bootstrap_seeded): whole-bag refits, (b, 0, 0) and (b, 1, 0)
# =================================================================================
if ("bootstrap" %in% ROUTES) {
  cat("\n=== bootstrap: whole-bag refits applied at eICU ===\n")
  for (b in BAGS_RUN) {
    bag <- BAG[[b]]
    # SHAP: the bag booster and its re-seeded refit, on every bag, as internally.
    for (s in c(0L, seq_len(B_PER_BAG))) {
      rt <- if (s == 0L) "bootstrap" else "bootstrap_seeded"
      if (have("shap_xgb_feat", rt, boot_id = b, seed_id = s)) next
      t0 <- start_timer()
      bst <- fit_booster(bag, cfg$seed + SEED_B * s)
      sh <- shap_eicu(bst$booster, sprintf("bag %d seed %d", b, s))
      emit(sh$G, "shap_xgb_feat", rt, boot_id = b, seed_id = s)
      .diag(rt, "shap_xgb_feat", boot_id = b, seed_id = s,
            additivity_max_abs = sh$additivity_max_abs, best_iter = bst$best_iter,
            elapsed_sec = t0()$elapsed_sec)
      cat(sprintf("  bag %2d shap %-16s (%.1f min)\n", b, rt, t0()$elapsed_sec / 60))
      rm(bst, sh)
    }
    # LLR: both arms from one pass of whole-bag fits; the internal exclusions
    # skipped. Each model is applied at eICU as soon as it is fitted and then
    # dropped, so one GAM is held at a time; the matrices are emitted only when
    # every spec of the bag has succeeded (a bag is excluded whole).
    if (b %in% TOMB$boot_id) { cat(sprintf("  bag %2d llr excluded (tombstoned)\n", b)); next }
    if (have("llr_meas", "bootstrap", boot_id = b) && have("llr_full", "bootstrap", boot_id = b)) next
    t0 <- start_timer()
    Mm <- .empty_M(); Mf <- .empty_M(); failed <- NULL
    for (i in seq_len(nrow(SPECS))) {
      sg <- SPECS$signal[i]; md <- SPECS$model[i]; key <- SPECS$key[i]
      r <- try(fit_one(sg, md, tabs_m, cfg, PRI[[sg]], fit_ids = tr[bag],
                       keep_model = TRUE, role = "final"), silent = TRUE)
      if (inherits(r, "try-error")) {
        msg <- conditionMessage(attr(r, "condition"))
        failed <- paste0(key, ": ", substr(gsub("\\s+", " ", msg), 1L, 160L))
        break
      }
      v <- apply_cached(r$model, key, sg); rm(r)
      if (md == "meas") Mm[ROWS[[key]], sg] <- v else Mf[ROWS[[key]], sg] <- v
    }
    if (!is.null(failed)) {
      tombstone_bag(b, paste0(if (grepl("distinct-value count", failed, fixed = TRUE))
                                "basis_exceeds_distinct_values" else "fit_error", ": ", failed))
      cat(sprintf("  bag %2d llr EXCLUDED whole: %s\n", b, failed))
      rm(Mm, Mf); next
    }
    for (sg in FULL_ALIAS) Mf[, sg] <- Mm[, sg]
    emit(Mm, "llr_meas", "bootstrap", boot_id = b)
    emit(Mf, "llr_full", "bootstrap", boot_id = b)
    rm(Mm, Mf); invisible(gc(verbose = FALSE))
    .diag("bootstrap", "llr_meas+llr_full", boot_id = b, elapsed_sec = t0()$elapsed_sec)
    cat(sprintf("  bag %2d llr %d fits (%.1f min)\n", b, nrow(SPECS), t0()$elapsed_sec / 60))
    mem_check(sprintf("bag %d", b))
  }
}

# =================================================================================
# ROUTE posterior: draws from the FINAL bundle GAMs, (0, 0, d), one arm at a time
# =================================================================================
rmvn_clamped <- function(nd, mu, V) {      # verbatim, tests/attr_replicates.R L1613
  e <- eigen(V, symmetric = TRUE)
  neg <- sum(e$values < 0)
  d <- sqrt(pmax(e$values, 0))
  Z <- matrix(stats::rnorm(nd * length(mu)), nd, length(mu))
  list(draws = sweep(Z %*% (t(e$vectors) * d), 2L, mu, `+`), n_clamped = neg)
}
#' The B_DRAW x n eICU log-odds draws of one final bundle model, minus logit(p_bar).
#' The seed is DRAW0 + the signal's index, the same in both arms, so an
#' unpaired signal's paired column is its measurement column draw for draw.
draw_spec <- function(sg, key) {
  b  <- bundle$models[[key]]
  V  <- if (!is.null(b$Vc)) b$Vc else b$Vp
  Xp <- stats::predict(b, newdata = ND[[key]], type = "lpmatrix", discrete = FALSE)
  dr <- with_seed(DRAW0 + match(sg, sigs), rmvn_clamped(B_DRAW, stats::coef(b), V))
  list(Eta = Xp %*% t(dr$draws) - logit(PRI[[sg]]$p_bar),
       diag = data.frame(spec = key, n_coef = length(stats::coef(b)), used_vc = !is.null(b$Vc),
                         n_eigen_clamped = dr$n_clamped, stringsAsFactors = FALSE))
}
if ("posterior" %in% ROUTES) {
  PD <- list()
  for (arm in c("llr_meas", "llr_full")) {
    if (all(vapply(seq_len(B_DRAW), function(d) have(arm, "posterior", draw_id = d), logical(1)))) next
    cat(sprintf("\n=== posterior, %s: %d draws from each final bundle GAM, applied at eICU ===\n",
                arm, B_DRAW))
    t0 <- start_timer()
    L <- replicate(B_DRAW, .empty_M(), simplify = FALSE)
    # The arm's model for each signal: the paired model where it is fitted, the
    # measurement model otherwise (and always, for the measurement arm).
    for (sg in sigs) {
      md  <- if (arm == "llr_full" && !(sg %in% FULL_ALIAS)) "full" else "meas"
      key <- paste0(sg, "/", md)
      ds  <- draw_spec(sg, key)
      for (d in seq_len(B_DRAW)) L[[d]][ROWS[[key]], sg] <- ds$Eta[, d]
      PD[[length(PD) + 1L]] <- cbind(arm = arm, ds$diag)
      rm(ds)
    }
    for (d in seq_len(B_DRAW)) emit(L[[d]], arm, "posterior", draw_id = d)
    rm(L); invisible(gc(verbose = FALSE))
    cat(sprintf("  %d posterior replicates written (%.1f min)\n", B_DRAW, t0()$elapsed_sec / 60))
    mem_check(sprintf("posterior %s", arm))
  }
  if (length(PD)) save_table(run, do.call(rbind, PD), "posterior_draw_diagnostics", subdir = "diagnostics")
}

# =================================================================================
# CLOSE: diagnostics, manifest, coverage, status
# =================================================================================
if (length(DIAG)) {
  DG <- do.call(rbind, DIAG)
  .dp <- file.path(run$path, "diagnostics", "fit_diagnostics.rds")
  if (file.exists(.dp)) {
    old <- readRDS(.dp)
    # A table written before 2026-10-06's review carried the refit-against-bundle
    # difference in `additivity_max_abs` and named the LLR bag rows `llr`; it is
    # converted to the present columns rather than dropped.
    if (!"max_abs_diff_vs_bundle" %in% names(old)) {
      va <- old$route == "verify_anchor"
      old$max_abs_diff_vs_bundle <- ifelse(va, old$additivity_max_abs, NA_real_)
      old$additivity_max_abs[va] <- NA_real_
      old$method[old$method == "llr"] <- "llr_meas+llr_full"
      old$r_max_used_gb <- NA_real_
      old <- old[, names(DG)]
    }
    DG <- rbind(old, DG)
  }
  # One row per (route, method, boot_id, seed_id): a resume that re-verifies the
  # anchor replaces its row rather than adding a second one. The benchmark's
  # four-spec verification rows are superseded by the 31-spec ones by this rule.
  DG <- DG[!duplicated(DG[, c("route", "method", "boot_id", "seed_id")], fromLast = TRUE), ]
  save_table(run, DG, "fit_diagnostics", subdir = "diagnostics")
}
MAN <- reconcile_manifest()
MAN <- MAN[order(MAN$method, MAN$route, MAN$boot_id, MAN$seed_id, MAN$draw_id), ]
utils::write.csv(MAN, MAN_P, row.names = FALSE)
attr_validate_manifest(MAN, "tests/attr_external_bags.R")
TB <- sort(unique(as.integer(TOMB$boot_id)))
COV <- attr_replicate_coverage(PLAN, MAN, tombstoned = TB)
save_table(run, COV, "replicate_coverage", subdir = "diagnostics")
cat("\n=== coverage: planned / present / tombstoned / missing ===\n\n")
print(COV, row.names = FALSE)
complete <- all(COV$missing == 0L)
.mx <- mem_max_gb()
.extra <- list(generator = EICU_GEN_VERSION, site = "eicu", internal_store = basename(INTERNAL),
               bundle = bpath, design_key = .DK, n_replicates = nrow(MAN),
               n_planned = sum(COV$planned), n_missing = sum(COV$missing),
               tombstoned_bags = paste(TB, collapse = ","),
               store_complete = if (complete) "yes" else "no",
               r_max_used_gb = round(.mx, 2),
               methods = paste(METHODS, collapse = ","))
if (complete) {
  finalize_run(run, extra = .extra)
} else {
  write_manifest(run, status = "running", extra = .extra)
}
cat(sprintf("\nstore %s: %d replicates, %s (%.1f min, R max used %.2f GB)\n", run$path, nrow(MAN),
            if (complete) "COMPLETE" else "incomplete; resume with --resume",
            t_all()$elapsed_sec / 60, .mx))
