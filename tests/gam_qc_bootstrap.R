# tests/gam_qc_bootstrap.R -----------------------------------------------------
# STABILITY QC: BOOTSTRAP EVIDENCE CURVES AND SURFACES FOR THE FINAL GAMs.
#
# Read docs/gam_qc_plan_20260908.md sections 5 and 6 first.
#
# --- WHY THIS SCRIPT FITS, WHEN A BOOTSTRAP STORE ALREADY EXISTS -------------
#
# out/runs/attrgen_20260906T130150 holds 38 LLR bootstrap replicates per method.
# EVERY ONE OF THEM IS AN L MATRIX (41,250 stays x 19 signals) AND NOT A FITTED
# MODEL. The generator's fold fits live inside one loop iteration and die with
# it (hard rule 6, and tests/attr_replicates.R's own header), so no smooth from
# any bag exists on disk. A per-patient L spread can be read from that store; a
# per-term CURVE spread cannot, because a curve is a function of a covariate and
# the store keeps only its value at each patient's own covariate.
#
# So the curves have to be refitted. THE ONE DESIGN DECISION THAT MATTERS is
# WHICH resamples to refit on, and the answer is the SHARED BAGS of
# config/attribution_eval.yml: bag b here is bag b of the L replicates, derived
# by the same patient-grouped draw from the same seed and VERIFIED by membership
# hash against attrgen's bootstrap_bags.csv before a single fit. That makes the
# curve spread and the already-measured per-patient L spread two readings of one
# resample rather than two experiments, and it is what lets a reviewer ask "when
# this patient's evidence moved under bag 7, which smooth moved?".
#
# --- WHAT A REPLICATE IS HERE ------------------------------------------------
#
# A FINAL-SHAPED fit -- the spec's formula on every in-bag training stay, with
# the bundle's FROZEN final priors -- evaluated with predict(type = "terms") on
# the exact grids tests/gam_qc.R wrote. Priors are held fixed for the same
# reason the L generator holds them fixed: the question is the sampling
# variability of the smooth conditional on the covariate construction, on the
# same footing as the L replicates. The direction of that omission is stated
# in the plan: it understates the spread by whatever prior estimation adds.
#
# WHAT REACHES DISK: one .qs2 per bag holding, for every (spec, term), the term's
# value at each grid point. Never a gam object (hard rule 6 in spirit: nothing
# here can be applied to a patient), never a row. Content-keyed by the bundle's
# design hash, the bag seed and the grid hash, so a bag fitted under a different
# design or a different grid cannot be summarised as this one.
#
# THE FITS ARE M-OUT-OF-N SUBSAMPLES (about 63% of patients), exactly as the L
# replicates are and for the same reason (`%in%` deduplicates), so the spread is
# inflated by about sqrt(n/m) = 1.26 relative to a full-size resample. Common to
# every curve and to the L replicates; stated, not corrected.
#
# --- WHAT IS REPORTED, PER SMOOTH, INSIDE THE SUPPORTED REGION ONLY ---------
#
#   sd_ratio        median over supported grid points of bootstrap SD over the
#                   analytic Vc-corrected SE. The analytic band is nearly free
#                   and SHAP has no counterpart; this says whether it can be
#                   trusted. Expected below 1 under heavy `ts` shrinkage (the
#                   penalised procedure moves less than Vc says it could) and
#                   above 1 where the smooth is genuinely unstable.
#   shape_cor       per bag, the correlation of the bag curve with the final
#                   curve over supported points; median and minimum over bags.
#   dir_agree       share of bags whose monotone direction inside support
#                   matches the final fit's.
#   extrema_agree   share of bags with the same number of interior extrema.
#   sign_stable     share of supported grid points where at least
#                   `sign_stable_frac` of bags agree with the final fit on the
#                   sign of the CENTRED curve.
#   fit_in_band     share of supported points where the final fit lies inside
#                   the bags' 2.5-97.5% interval.
#   sd_in / sd_out  bootstrap SD inside against outside support -- the
#                   measured cost of reading a curve where nobody is.
#
# AGGREGATES ONLY (hard rule 1). Bags are reported as counts and hashes.
#
#   Rscript tests/gam_qc_bootstrap.R --plan-only
#   Rscript tests/gam_qc_bootstrap.R                          # latest gamqc, n_bags from config
#   Rscript tests/gam_qc_bootstrap.R --qc out/runs/gamqc_... --bags 20
#   Rscript tests/gam_qc_bootstrap.R --resume out/runs/gamboot_...   # continue
#   Rscript tests/gam_qc_bootstrap.R --resume <dir> --summarise-only  # no fit
#   Rscript tests/gam_qc_bootstrap.R --signals mbp,spo2 --bags 3      # a look
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/gam_qc_common.R")   # helpers shared with tests/gam_qc.R, defined once

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
QC_D      <- .opt("--qc")
RESUME    <- .opt("--resume")
ONLY_SG   <- .opt("--signals"); if (!is.null(ONLY_SG)) ONLY_SG <- strsplit(ONLY_SG, ",")[[1]]
ONLY_MD   <- .opt("--models");  if (!is.null(ONLY_MD)) ONLY_MD <- strsplit(ONLY_MD, ",")[[1]]
PLAN_ONLY <- "--plan-only" %in% args
SUMM_ONLY <- "--summarise-only" %in% args
NO_FIGS   <- "--no-figs" %in% args

# --- config -------------------------------------------------------------------
qcfg <- yaml::read_yaml("config/gam_qc.yml")
ecfg <- yaml::read_yaml("config/attribution_eval.yml")
BUNDLE_P  <- as.character(cfg_req(qcfg, "bundle"))
ATTRGEN_D <- as.character(cfg_req(qcfg, "attrgen"))
BC        <- cfg_req(qcfg, "bootstrap")
for (k in c("n_bags", "sign_stable_frac", "shape_cor_min", "sec_per_fit", "sec_per_fit_ti")) cfg_req(BC, k)
Z         <- as.numeric(cfg_req(qcfg, "band", "z"))
STABLE    <- as.numeric(BC$sign_stable_frac)
SHAPE_MIN <- as.numeric(BC$shape_cor_min)
# The SAME direction and flatness rules the QC script reads, so "reproduces
# the final fit's direction" means the direction the QC table reported.
RHO_MONO  <- as.numeric(cfg_req(qcfg, "interp", "monotone_rho_min"))
FLAT_NATS <- as.numeric(cfg_req(qcfg, "interp", "flat_range_nats"))
N_BAGS    <- as.integer(.opt("--bags", BC$n_bags))
# THE SHARED BAGS: declared once, in the attribution plan, and read from there.
BAG_SEED0 <- as.integer(cfg_req(ecfg, "bootstrap", "seed_base"))
B_BOOT    <- as.integer(cfg_req(ecfg, "bootstrap", "b"))
if (N_BAGS < 1L || N_BAGS > B_BOOT) {
  stop("--bags must be within the shared manifest 1..", B_BOOT, " (config/attribution_eval.yml ",
       "bootstrap.b); a bag outside it is not a bag any L replicate used.", call. = FALSE)
}

cfg_local <- load_config("config/config.yml")
bundle    <- load_bundle(BUNDLE_P, cfg = cfg_local, strict = TRUE)
cfg       <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
DESIGN_H  <- .hash(bundle$cfg)

# --- the QC run whose grids this reproduces -----------------------------------
if (is.null(QC_D)) QC_D <- latest_run("gamqc")
if (is.null(QC_D) || !dir.exists(QC_D)) {
  stop("no gamqc run. Run tests/gam_qc.R first, or pass --qc <dir>.", call. = FALSE)
}
PROV <- readRDS(file.path(QC_D, "diagnostics", "provenance.rds"))
qc_design <- PROV$value[PROV$item == "bundle_design_hash"]
if (!identical(qc_design, DESIGN_H)) {
  stop("gam_qc_bootstrap: the QC run ", basename(QC_D), " was made under design ", qc_design,
       " and this bundle is ", DESIGN_H, ". Its grids describe other smooths.", call. = FALSE)
}
CUR1 <- readRDS(file.path(QC_D, "tables", "curves_1d.rds"))
SUR2 <- if (file.exists(file.path(QC_D, "tables", "surfaces_2d.rds")))
  readRDS(file.path(QC_D, "tables", "surfaces_2d.rds")) else NULL
# Hashed as bare columns, not as data frames: a data frame's row names ride
# along in the digest and would make two identical grids hash apart.
GRID_H <- attr_key_hash(list(
  c1 = list(CUR1$key, CUR1$term, CUR1$grid_x),
  s2 = if (is.null(SUR2)) NULL else list(SUR2$key, SUR2$term, SUR2$x1, SUR2$x2)))

SPECS <- unique(CUR1[, c("key", "signal", "model")])
if (!is.null(ONLY_SG)) SPECS <- SPECS[SPECS$signal %in% ONLY_SG, , drop = FALSE]
if (!is.null(ONLY_MD)) SPECS <- SPECS[SPECS$model  %in% ONLY_MD, , drop = FALSE]
rownames(SPECS) <- NULL
if (!nrow(SPECS)) stop("no spec selected", call. = FALSE)
n_ti <- sum(SPECS$model %in% LAYER1_TI_MODELS)

cat("\n=== plan ===\n")
cat(sprintf("  bundle          : %s  (design %s)\n", BUNDLE_P, DESIGN_H))
cat(sprintf("  grids from      : %s  (grid hash %s)\n", basename(QC_D), GRID_H))
cat(sprintf("  specs           : %d (%d interaction)\n", nrow(SPECS), n_ti))
cat(sprintf("  bags            : %d of the %d shared bags (seed_base %d)\n", N_BAGS, B_BOOT, BAG_SEED0))
cat(sprintf("  fits            : %d final-shaped fits, none kept\n", N_BAGS * nrow(SPECS)))
cat(sprintf("  rough wall clock: %.0f min at %s s per additive fit and %s s per tensor fit\n",
            N_BAGS * ((nrow(SPECS) - n_ti) * as.numeric(BC$sec_per_fit) + n_ti * as.numeric(BC$sec_per_fit_ti)) / 60,
            BC$sec_per_fit, BC$sec_per_fit_ti))
if (PLAN_ONLY) quit(save = "no", status = 0L)

# --- data and the bags --------------------------------------------------------
tabs  <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
tr    <- folds$stay_id[folds$split == "train"]       # the order tar_read(train_ids) has

#' Bag `seed`, as tests/attr_replicates.R draws it. A CHECKED DUPLICATE: the
#' definition of record lives in that script, and every bag drawn here is
#' compared by membership hash to the one recorded there before use. A hash
#' that differs stops the run -- the two definitions have drifted, and the fix
#' is to hoist bag_of() into R/ so there is one, not to trust either copy.
bag_of <- function(seed) {
  gcol <- resample_cols(tabs$cohort, cfg)$fold_group
  src  <- if (gcol %in% names(folds)) folds else tabs$cohort
  gid  <- as.character(src[[gcol]])[match(tr, src$stay_id)]
  if (anyNA(gid)) stop("bag_of: `", gcol, "` is missing for some training stay", call. = FALSE)
  ug <- unique(gid)
  with_seed(seed, {
    keep <- sample(ug, length(ug), replace = TRUE)
    gid %in% unique(keep)
  })
}
REF_BAGS <- NULL
rb <- file.path(ATTRGEN_D, "diagnostics", "bootstrap_bags.rds")
if (file.exists(rb)) REF_BAGS <- readRDS(rb)
BAGS <- lapply(seq_len(N_BAGS), function(b) bag_of(BAG_SEED0 + b))
BAGT <- do.call(rbind, lapply(seq_len(N_BAGS), function(b) {
  g <- BAGS[[b]]; h <- attr_key_hash(which(g))
  ref <- if (is.null(REF_BAGS)) NA_character_ else REF_BAGS$membership_hash[match(b, REF_BAGS$boot_id)]
  data.frame(boot_id = b, seed = BAG_SEED0 + b, n_in_bag = sum(g), n_train = length(g),
             frac_in_bag = round(mean(g), 5), membership_hash = h, reference_hash = ref,
             verified = identical(h, ref), stringsAsFactors = FALSE)
}))
if (is.null(REF_BAGS)) {
  cat(sprintf("\n  NO REFERENCE BAG TABLE at %s: bags are drawn but UNVERIFIED against the L replicates.\n", rb))
} else if (!all(BAGT$verified)) {
  print(BAGT[!BAGT$verified, c("boot_id", "membership_hash", "reference_hash")], row.names = FALSE)
  stop("gam_qc_bootstrap: ", sum(!BAGT$verified), " bag(s) differ from the shared manifest in ",
       basename(ATTRGEN_D), ". The bag definition here and in tests/attr_replicates.R have ",
       "drifted, or the training ids moved. Name the branch; do not refit.", call. = FALSE)
} else {
  cat(sprintf("\n  %d bag(s) verified by membership hash against %s (in-bag %.3f to %.3f)\n",
              N_BAGS, basename(ATTRGEN_D), min(BAGT$frac_in_bag), max(BAGT$frac_in_bag)))
}

# --- the run ------------------------------------------------------------------
if (!is.null(RESUME)) {
  if (!dir.exists(RESUME)) stop("--resume: no such directory ", RESUME, call. = FALSE)
  run <- structure(list(prefix = "gamboot", id = basename(RESUME), path = RESUME,
                        started = Sys.time(), config = cfg_local,
                        log_file = file.path(RESUME, "log.txt")), class = "llr_run")
  log_msg(run, "resumed")
} else {
  run <- new_run("gamboot", cfg_local, note = sprintf(
    "bootstrap evidence curves: %d spec(s) x %d shared bag(s); grids from %s",
    nrow(SPECS), N_BAGS, basename(QC_D)))
}
SUB <- file.path(run$path, "replicates"); dir.create(SUB, showWarnings = FALSE)
save_table(run, BAGT, "bootstrap_bags", subdir = "diagnostics")
bag_file <- function(b) file.path(SUB, sprintf("bag_%02d__%s.qs2", b,
  attr_key_hash(list(design = DESIGN_H, seed = BAG_SEED0 + b, grid = GRID_H))))
FAIL_P <- file.path(run$path, "replicates_failed.csv")

# ==============================================================================
# THE FITS
# ==============================================================================
if (!SUMM_ONLY) {
  cat("\n=== fitting ===\n")
  for (b in seq_len(N_BAGS)) {
    if (file.exists(bag_file(b))) { cat(sprintf("  bag %2d  cached\n", b)); next }
    t0 <- start_timer()
    ids <- tr[BAGS[[b]]]
    curves <- list(); surfaces <- list(); failed <- list(); n_fit <- 0L
    # One frame per (signal, column set) per bag: the two interaction models
    # share `full`'s columns exactly, so this is three builds per paired
    # signal rather than five. See qc_frame_for().
    fcache <- new.env(parent = emptyenv())
    for (.i in seq_len(nrow(SPECS))) {
      key <- SPECS$key[.i]; sg <- SPECS$signal[.i]; md <- SPECS$model[.i]
      pri <- priors_for(bundle$priors, sg, "final")
      r <- try({
        d  <- qc_frame_for(fcache, sg, md, tabs, cfg, pri, stay_ids = ids)
        bm <- .bam_fit(attr(d, "formula"), d, cfg)
        vars <- setdiff(names(d), c("stay_id", "mortality"))
        nr0 <- qc_neutral_row(d, vars)
        c1 <- CUR1[CUR1$key == key, ]
        for (t in unique(c1$term)) {
          z <- c1[c1$term == t, ]; v <- z$variable[1]
          nd <- nr0[rep(1L, nrow(z)), , drop = FALSE]; nd[[v]] <- z$grid_x
          p <- stats::predict(bm, newdata = nd, type = "terms", discrete = FALSE)
          curves[[paste(key, t)]] <- as.numeric(p[, match(t, colnames(p))])
        }
        if (!is.null(SUR2)) {
          s2 <- SUR2[SUR2$key == key, ]
          for (t in unique(s2$term)) {
            z <- s2[s2$term == t, ]
            nd <- nr0[rep(1L, nrow(z)), , drop = FALSE]; nd[[z$var1[1]]] <- z$x1; nd[[z$var2[1]]] <- z$x2
            p <- stats::predict(bm, newdata = nd, type = "terms", discrete = FALSE)
            surfaces[[paste(key, t)]] <- as.numeric(p[, match(t, colnames(p))])
          }
        }
        rm(bm, d); TRUE
      }, silent = TRUE)
      if (inherits(r, "try-error")) {
        msg <- conditionMessage(attr(r, "condition"))
        failed[[length(failed) + 1L]] <- data.frame(
          boot_id = b, key = key, signal = sg, model = md,
          cause = if (grepl("distinct-value count", msg, fixed = TRUE)) "basis_exceeds_distinct_values"
                  else if (grepl("constant column", msg, fixed = TRUE)) "constant_column" else "fit_error",
          detail = substr(gsub("\\s+", " ", msg), 1L, 200L), stringsAsFactors = FALSE)
        next
      }
      n_fit <- n_fit + 1L
    }
    FAILED <- if (length(failed)) do.call(rbind, failed) else NULL
    qs2::qs_save(list(boot_id = b, seed = BAG_SEED0 + b, design = DESIGN_H, grid = GRID_H,
                      n_in_bag = length(ids), n_fit = n_fit, curves = curves, surfaces = surfaces,
                      failed = FAILED), bag_file(b))
    if (!is.null(FAILED)) {
      utils::write.csv(rbind(if (file.exists(FAIL_P)) utils::read.csv(FAIL_P, stringsAsFactors = FALSE) else NULL,
                             FAILED), FAIL_P, row.names = FALSE)
    }
    el <- t0()
    log_msg(run, sprintf("bag %2d: %d fit(s), %d failed, %.1f min", b, n_fit,
                         if (is.null(FAILED)) 0L else nrow(FAILED), el$elapsed_sec / 60))
  }
}

# ==============================================================================
# THE SUMMARY
# ==============================================================================
files <- vapply(seq_len(N_BAGS), bag_file, character(1))
files <- files[file.exists(files)]
if (!length(files)) stop("no bag file on disk for this design and grid", call. = FALSE)
REPS <- lapply(files, qs2::qs_read)
cat(sprintf("\n=== summarising %d bag(s) ===\n\n", length(REPS)))
FAILED <- do.call(rbind, Filter(Negate(is.null), lapply(REPS, `[[`, "failed")))
if (!is.null(FAILED)) {
  save_table(run, FAILED, "replicates_failed", subdir = "diagnostics")
  cat(sprintf("  %d (bag, spec) fit(s) failed and are counted, not silently dropped:\n", nrow(FAILED)))
  print(table(FAILED$key, FAILED$cause))
}

# The extremum counter and the direction rule are the QC script's own
# (tests/gam_qc_common.R), bound here to the same config values it read.
n_extrema <- function(f) qc_n_extrema(f)
dir_of    <- function(x, f) qc_dir_of(x, f, RHO_MONO)

BC1 <- list(); TSUM <- list()
for (.i in seq_len(nrow(SPECS))) {
  key <- SPECS$key[.i]
  c1 <- CUR1[CUR1$key == key, ]
  for (t in unique(c1$term)) {
    z <- c1[c1$term == t, ]; nm <- paste(key, t)
    M <- do.call(rbind, Filter(Negate(is.null), lapply(REPS, function(r) r$curves[[nm]])))
    if (is.null(M) || nrow(M) < 2L) next
    sup <- z$supported; i <- which(sup)
    bm <- colMeans(M); bs <- apply(M, 2, stats::sd)
    q  <- apply(M, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
    BC1[[length(BC1) + 1L]] <- data.frame(
      key = key, signal = z$signal[1], model = z$model[1], term = t, variable = z$variable[1],
      grid_x = z$grid_x, fit = z$fit, se = z$se, supported = sup, n_bags = nrow(M),
      boot_mean = round(bm, 6), boot_sd = round(bs, 6), boot_q025 = round(q[1, ], 6), boot_q975 = round(q[2, ], 6),
      stringsAsFactors = FALSE)
    if (length(i) < 3L) next
    fi <- z$fit[i]; xi <- z$grid_x[i]; Mi <- M[, i, drop = FALSE]
    fc <- fi - mean(fi); Mc <- sweep(Mi, 1L, rowMeans(Mi))
    # A final curve shrunk to a flat line has no shape to reproduce: its
    # correlation, direction and sign pattern are round-off and are reported
    # as NA rather than as instability.
    flat <- (max(fi) - min(fi)) < FLAT_NATS
    shape_cor <- if (flat) NA_real_ else
      apply(Mi, 1, function(r) if (stats::sd(r) > 0 && stats::sd(fi) > 0) stats::cor(r, fi) else NA_real_)
    d0 <- if (flat) 0 else dir_of(xi, fi)
    agree_sign <- colMeans(sweep(sign(Mc), 2L, sign(fc), `==`))
    TSUM[[length(TSUM) + 1L]] <- data.frame(
      key = key, signal = z$signal[1], model = z$model[1], term = t, variable = z$variable[1],
      n_bags = nrow(M), n_supported = length(i),
      boot_sd_med_in = round(stats::median(bs[i]), 5),
      boot_sd_med_out = if (any(!sup)) round(stats::median(bs[!sup]), 5) else NA_real_,
      se_med_in = round(stats::median(z$se[i]), 5),
      sd_ratio_med = round(stats::median(bs[i] / z$se[i]), 4),
      sd_ratio_max = round(max(bs[i] / z$se[i]), 4),
      flat_final = flat,
      fit_in_band = round(mean(fi >= q[1, i] & fi <= q[2, i]), 4),
      shape_cor_med = if (flat) NA_real_ else round(stats::median(shape_cor, na.rm = TRUE), 4),
      shape_cor_min = if (flat) NA_real_ else round(min(shape_cor, na.rm = TRUE), 4),
      dir_final = d0,
      dir_agree = if (flat) NA_real_ else round(mean(apply(Mi, 1, function(r) dir_of(xi, r)) == d0), 4),
      n_extrema_final = if (flat) NA_integer_ else n_extrema(fi),
      extrema_agree = if (flat) NA_real_ else round(mean(apply(Mi, 1, n_extrema) == n_extrema(fi)), 4),
      extrema_within1 = if (flat) NA_real_ else round(mean(abs(apply(Mi, 1, n_extrema) - n_extrema(fi)) <= 1L), 4),
      sign_stable = if (flat) NA_real_ else round(mean(agree_sign >= STABLE), 4),
      level_sd = round(stats::sd(rowMeans(Mi)), 5),
      range_in_final = round(max(fi) - min(fi), 5),
      range_in_boot_q025 = round(unname(stats::quantile(apply(Mi, 1, function(r) max(r) - min(r)), 0.025)), 5),
      stringsAsFactors = FALSE)
  }
}
BC1 <- do.call(rbind, BC1); TSUM <- do.call(rbind, TSUM)

BS2 <- NULL; SSUM <- NULL
if (!is.null(SUR2)) {
  bs2 <- list(); ss <- list()
  for (.i in seq_len(nrow(SPECS))) {
    key <- SPECS$key[.i]; s2 <- SUR2[SUR2$key == key, ]
    for (t in unique(s2$term)) {
      z <- s2[s2$term == t, ]; nm <- paste(key, t)
      M <- do.call(rbind, Filter(Negate(is.null), lapply(REPS, function(r) r$surfaces[[nm]])))
      if (is.null(M) || nrow(M) < 2L) next
      sup <- z$supported; i <- which(sup)
      bs <- apply(M, 2, stats::sd); bm <- colMeans(M)
      bs2[[length(bs2) + 1L]] <- data.frame(
        key = key, term = t, x1 = z$x1, x2 = z$x2, fit = z$fit, se = z$se, supported = sup,
        n_bags = nrow(M), boot_mean = round(bm, 6), boot_sd = round(bs, 6), stringsAsFactors = FALSE)
      if (length(i) < 3L) next
      fi <- z$fit[i]; Mi <- M[, i, drop = FALSE]
      sc <- apply(Mi, 1, function(r) if (stats::sd(r) > 0 && stats::sd(fi) > 0) stats::cor(r, fi) else NA_real_)
      ss[[length(ss) + 1L]] <- data.frame(
        key = key, signal = z$signal[1], model = z$model[1], term = t, var1 = z$var1[1], var2 = z$var2[1],
        n_bags = nrow(M), n_supported = length(i),
        boot_sd_med_in = round(stats::median(bs[i]), 5),
        boot_sd_med_out = if (any(!sup)) round(stats::median(bs[!sup]), 5) else NA_real_,
        se_med_in = round(stats::median(z$se[i]), 5),
        sd_ratio_med = round(stats::median(bs[i] / z$se[i]), 4),
        surface_cor_med = round(stats::median(sc, na.rm = TRUE), 4),
        surface_cor_min = round(min(sc, na.rm = TRUE), 4),
        range_in_final = round(max(fi) - min(fi), 5),
        range_in_boot_q025 = round(unname(stats::quantile(apply(Mi, 1, function(r) max(r) - min(r)), 0.025)), 5),
        range_in_boot_med = round(stats::median(apply(Mi, 1, function(r) max(r) - min(r))), 5),
        stringsAsFactors = FALSE)
    }
  }
  BS2 <- do.call(rbind, bs2); SSUM <- do.call(rbind, ss)
}

# spec-level roll-up
.min_na <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
SP <- do.call(rbind, lapply(split(TSUM, TSUM$key), function(z) data.frame(
  key = z$key[1], signal = z$signal[1], model = z$model[1], n_terms_1d = nrow(z),
  n_terms_flat = sum(z$flat_final),
  n_bags = min(z$n_bags), sd_ratio_med = round(stats::median(z$sd_ratio_med), 4),
  sd_ratio_max = round(max(z$sd_ratio_max), 4),
  shape_cor_min = round(.min_na(z$shape_cor_med), 4),
  n_terms_dir_unstable = sum(z$dir_agree < STABLE & z$dir_final != 0, na.rm = TRUE),
  n_terms_shape_unstable = sum(z$shape_cor_med < SHAPE_MIN, na.rm = TRUE),
  sign_stable_min = round(.min_na(z$sign_stable), 4),
  stringsAsFactors = FALSE)))
if (!is.null(SSUM)) {
  m <- match(SP$key, SSUM$key)
  agg <- do.call(rbind, lapply(split(SSUM, SSUM$key), function(z) data.frame(
    key = z$key[1], n_terms_2d = nrow(z), surf_sd_ratio_med = round(stats::median(z$sd_ratio_med), 4),
    surf_cor_min = round(min(z$surface_cor_med), 4), stringsAsFactors = FALSE)))
  SP <- merge(SP, agg, by = "key", all.x = TRUE)
}
SP <- SP[match(intersect(SPECS$key, SP$key), SP$key), ]; rownames(SP) <- NULL

save_table(run, BC1,  "boot_curves_1d", csv = FALSE)
save_table(run, TSUM, "boot_term_summary")
if (!is.null(BS2)) { save_table(run, BS2, "boot_surfaces_2d", csv = FALSE); save_table(run, SSUM, "boot_surface_summary") }
save_table(run, SP, "boot_spec_summary")

# --- figures: spaghetti of bag curves against the analytic band -------------
if (!NO_FIGS) {
  for (key in unique(BC1$key)) {
    z0 <- BC1[BC1$key == key, ]; terms <- unique(z0$term); n <- length(terms)
    nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
    save_fig(run, paste0("boot_", gsub("/", "_", key)), width = 4 * nc, height = 3.2 * nr, dpi = 120)
    graphics::par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2.5, 0.8), mgp = c(2.1, 0.7, 0))
    for (t in terms) {
      # `o` is taken BEFORE z is reordered, because the bag matrices are in
      # the stored grid order and must be permuted by the same index.
      z <- z0[z0$term == t, ]; o <- order(z$grid_x); z <- z[o, ]
      M <- do.call(rbind, Filter(Negate(is.null), lapply(REPS, function(r) r$curves[[paste(key, t)]])))
      ts <- TSUM[TSUM$key == key & TSUM$term == t, ]
      yl <- range(c(z$fit - Z * z$se, z$fit + Z * z$se, M), na.rm = TRUE)
      graphics::plot(z$grid_x, z$fit, type = "n", ylim = yl, xlab = z$variable[1], ylab = "partial effect (nats)",
                     main = if (nrow(ts)) sprintf("%s  sd ratio %.2f  shape cor %.2f", t, ts$sd_ratio_med, ts$shape_cor_med) else t,
                     cex.main = 0.85)
      r <- rle(!z$supported); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
      for (k in which(r$values)) {
        x0 <- if (starts[k] > 1L) (z$grid_x[starts[k] - 1L] + z$grid_x[starts[k]]) / 2 else graphics::par("usr")[1]
        x1 <- if (ends[k] < nrow(z)) (z$grid_x[ends[k]] + z$grid_x[ends[k] + 1L]) / 2 else graphics::par("usr")[2]
        graphics::rect(x0, yl[1] - diff(yl), x1, yl[2] + diff(yl), col = "#00000014", border = NA)
      }
      for (j in seq_len(nrow(M))) graphics::lines(z$grid_x, M[j, o], col = "#80808040", lwd = 0.7)
      graphics::lines(z$grid_x, z$fit - Z * z$se, lty = 2, col = "#1F3F7A")
      graphics::lines(z$grid_x, z$fit + Z * z$se, lty = 2, col = "#1F3F7A")
      graphics::lines(z$grid_x, z$fit, lwd = 2, col = "#1F3F7A")
      graphics::abline(h = 0, lty = 3, col = "grey40")
    }
    grDevices::dev.off()
  }
}

# --- the reading ---------------------------------------------------------------
cat(sprintf("%-28s %5s %8s %8s %7s %7s %7s\n", "spec", "bags", "sd_ratio", "sd_rmax", "shp_min", "dir-", "sign"))
for (i in seq_len(nrow(SP))) cat(sprintf("%-28s %5d %8.3f %8.3f %7.3f %7d %7.2f\n", SP$key[i], SP$n_bags[i],
  SP$sd_ratio_med[i], SP$sd_ratio_max[i], SP$shape_cor_min[i], SP$n_terms_dir_unstable[i], SP$sign_stable_min[i]))
cat(sprintf("\n  %d smooth(s) over %d spec(s), %d bag(s)\n", nrow(TSUM), nrow(SP), length(REPS)))
cat(sprintf("  bootstrap SD / analytic Vc SE inside support: median %.3f, IQR %.3f-%.3f, max %.3f (%s %s)\n",
            stats::median(TSUM$sd_ratio_med), stats::quantile(TSUM$sd_ratio_med, 0.25),
            stats::quantile(TSUM$sd_ratio_med, 0.75), max(TSUM$sd_ratio_max),
            TSUM$key[which.max(TSUM$sd_ratio_max)], TSUM$variable[which.max(TSUM$sd_ratio_max)]))
cat("    (below 1: the penalised procedure moves less than Vc says it could, so the analytic band is conservative;\n")
cat("     above 1: the band understates sampling variability. Inflated by ~1.26 by the 63% subsample.)\n")
cat(sprintf("  bootstrap SD outside support / inside: median ratio %.2f\n",
            stats::median(TSUM$boot_sd_med_out / TSUM$boot_sd_med_in, na.rm = TRUE)))
cat(sprintf("  smooths shrunk to a flat line in the final fit (no shape to reproduce): %d of %d\n",
            sum(TSUM$flat_final), nrow(TSUM)))
cat(sprintf("  shape correlation with the final curve: median %.3f; smooths with median below %.2f: %d; below 0.5: %d\n",
            stats::median(TSUM$shape_cor_med, na.rm = TRUE), SHAPE_MIN,
            sum(TSUM$shape_cor_med < SHAPE_MIN, na.rm = TRUE), sum(TSUM$shape_cor_med < 0.5, na.rm = TRUE)))
cat(sprintf("  monotone direction reproduced in >= %.0f%% of bags: %d of %d smooths with a direction\n",
            100 * STABLE, sum(TSUM$dir_agree >= STABLE & TSUM$dir_final != 0, na.rm = TRUE),
            sum(TSUM$dir_final != 0, na.rm = TRUE)))
cat(sprintf("  extrema count reproduced exactly: median %.2f of bags; within one: median %.2f\n",
            stats::median(TSUM$extrema_agree, na.rm = TRUE), stats::median(TSUM$extrema_within1, na.rm = TRUE)))
cat(sprintf("  final fit inside the bags' 95%% band: median %.2f of the supported grid\n",
            stats::median(TSUM$fit_in_band, na.rm = TRUE)))
if (!is.null(SSUM)) {
  cat(sprintf("  ti surfaces: %d; sd ratio median %.3f; surface correlation median %.3f, min %.3f;\n",
              nrow(SSUM), stats::median(SSUM$sd_ratio_med), stats::median(SSUM$surface_cor_med), min(SSUM$surface_cor_med)))
  cat(sprintf("    supported range: final median %.3f nats, bootstrap 2.5%% quantile median %.3f -- an interaction whose\n",
              stats::median(SSUM$range_in_final), stats::median(SSUM$range_in_boot_q025)))
  cat("    lower quantile is near zero is one the resample does not reliably see\n")
}
finalize_run(run, extra = list(n_bags = length(REPS), n_specs = nrow(SP), qc_run = basename(QC_D),
                               grid_hash = GRID_H, design_hash = DESIGN_H))
cat(sprintf("\nwritten: %s\n", run$path))
