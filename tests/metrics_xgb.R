# tests/metrics_xgb.R --------------------------------------------------------
# The XGBoost comparison arm. Reads a finished branch_point run, fits nothing in
# layer 1, and scores three cells through the SAME code path as the proposed
# method (R/09's score_report) so the numbers are commensurable:
#
#   llr_sum    the naive summed LLR                    <- the proposed method
#   xgb_l      XGBoost on the out-of-fold L matrix     <- isolates AGGREGATION
#   xgb_feat   XGBoost on the layer-1 covariates      <- isolates the COVARIATE
#                                                        construction (--feat)
#   xgb_raw    XGBoost on the raw feature matrix       <- isolates the whole
#                                                         approach; this is the
#                                                         headline baseline
#
# Reading the three:
#   llr_sum vs xgb_raw   does the method hold up against a strong learner given
#                        the same information?
#   xgb_l vs llr_sum     same input, different aggregator. Where a nonlinear
#                        aggregator would beat summing, and therefore where
#                        Sigma-inverse has to earn its place.
#   xgb_l vs xgb_raw     same learner and same attribution method, only the
#                        input changes. Isolates the REPRESENTATION.
#
# TWO ASYMMETRIES THAT FAVOUR THE BASELINE, both deliberate:
#   - `xgb_raw` gets NA for unmeasured stays and XGBoost learns a default split
#     direction, so it may use missingness as evidence -- the exact mechanism the
#     LLR design refuses (spec SS5.5). `--no-missing` re-runs it with imputation
#     to measure what that channel is worth.
#   - `xgb_raw` gets `n_obs`, which no formula may carry. Removing it would be
#     handicapping the baseline.
# Both make the comparison conservative for the proposed method. Say so.
#
# EICU. Everything here routes through R/09b, which separates fitting from
# applying: `xgb_oof()` cross-fits internally, `xgb_fit_full()` + `xgb_apply()`
# are the external path. This script also writes the full-train booster to the
# run directory so an eICU run can load and apply it without refitting
# (hard rule 8).
#
# AGGREGATES ONLY (hard rule 1). Scores are row-level and never printed.
#
#   Rscript tests/metrics_xgb.R                        # newest branch run
#   Rscript tests/metrics_xgb.R out/runs/branch_...    # a specific one
#   Rscript tests/metrics_xgb.R - --no-missing         # the missingness ablation
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
src_dir <- if (length(args) >= 1L && args[1] != "-") args[1] else latest_run("branch")
no_miss <- "--no-missing" %in% args
do_feat <- "--feat" %in% args
if (is.null(src_dir) || !dir.exists(src_dir)) {
  stop("no completed branch run found. Run tests/branch_point.R first, or pass a run directory.",
       call. = FALSE)
}
l_path <- file.path(src_dir, "tables", "l_oof.rds")
if (!file.exists(l_path)) stop("no l_oof.rds in ", src_dir, call. = FALSE)

cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)

l_long   <- readRDS(l_path)
stay_ids <- sort(unique(l_long$stay_id))
stay_ids <- intersect(folds$stay_id[folds$split == "train"], stay_ids)

y <- tabs$cohort$mortality[match(as.character(stay_ids), as.character(tabs$cohort$stay_id))]
if (anyNA(y)) stop("metrics_xgb: a scored stay has no cohort row", call. = FALSE)
fold_of <- folds$fold[match(stay_ids, folds$stay_id)]
if (anyNA(fold_of)) stop("metrics_xgb: a scored stay has no fold", call. = FALSE)
# Patient per stay: the boosters' inner split and the bootstrap unit (review
# S3/S4, 2026-09-09), resolved through the accessor the folds use.
grp <- patient_group_of(tabs$cohort, cfg, stay_ids)
p_bar <- mean(y)

run <- new_run("xgb", cfg, note = sprintf(
  "XGBoost comparison arm from %s; %d stays; missing_as_evidence = %s",
  basename(src_dir), length(stay_ids), !no_miss))
log_msg(run, sprintf("source run: %s | %d stays | event rate %.4f | folds %s",
                     basename(src_dir), length(stay_ids), p_bar,
                     paste(sort(unique(fold_of)), collapse = ",")))

nb    <- cfg$metrics$n_bins %||% 20L
nboot <- cfg$metrics$n_boot %||% 200L

# --- designs ----------------------------------------------------------------
M_full <- xgb_design_L(l_long, tabs, cfg, stay_ids, model = "full", fill = "zero")
X_raw  <- xgb_design_raw(tabs, cfg, stay_ids, missing_as_evidence = !no_miss)
log_msg(run, sprintf("designs: L %d x %d, raw %d x %d (%.1f%% NA)",
                     nrow(M_full), ncol(M_full), nrow(X_raw), ncol(X_raw),
                     100 * mean(is.na(X_raw))))

# The constructed-covariate design is FOLD-DEPENDENT: `delta` and `lambda` are
# fitted, so one matrix cannot serve five folds without leaking. Built per fold
# by xgb_oof_perfold() through this closure.
if (do_feat) {
  t0 <- Sys.time()
  priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
  log_msg(run, sprintf("layer1_priors: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  feat_fn <- function(f) xgb_design_feat(tabs, cfg, priors, stay_ids, role = "oof", fold = f)
  X_feat1 <- feat_fn(1L)
  log_msg(run, sprintf("designs: feat %d x %d (%.1f%% NA), fold-1 instance",
                       nrow(X_feat1), ncol(X_feat1), 100 * mean(is.na(X_feat1))))
}

# --- the three cells --------------------------------------------------------
res <- list()

res$llr_sum <- score_report(run, rowSums(M_full), y, p_bar, label = "llr_sum",
                            n_bins = nb, n_boot = nboot, seed = cfg$seed,
                            title = "Out-of-fold risk ordering — summed LLR")

t0 <- Sys.time()
# NOT THE HONEST INTERNAL ESTIMATE (statistical review S1, 2026-09-09). This
# cross-validates a booster over the primary out-of-fold L matrix on the same
# folds, so its training features were fitted with the held-out fold's
# outcomes. The graph's `oof_xgb_l` target is the nested version
# (`xgb_design_L_nested()`); this script keeps the single-level cell only as a
# fast ladder comparison and its `xgb_l` row must not be quoted as the
# internal result.
log_msg(run, "xgb_l here is the SINGLE-LEVEL stacked estimate (review S1); quote the graph's nested oof_xgb_l")
oL <- xgb_oof(M_full, y, fold_of, cfg, seed = cfg$seed, group = grp)
log_msg(run, sprintf("xgb_l: %d folds, best_iter %s, %.1f s",
                     length(oL$best_iter), paste(oL$best_iter, collapse = "/"),
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
res$xgb_l <- score_report(run, oL$score, y, p_bar, label = "xgb_l",
                          n_bins = nb, n_boot = nboot, seed = cfg$seed,
                          title = "Out-of-fold risk ordering — XGBoost on L")

t0 <- Sys.time()
oR <- xgb_oof(X_raw, y, fold_of, cfg, seed = cfg$seed, group = grp)
log_msg(run, sprintf("xgb_raw: %d folds, best_iter %s, %.1f s",
                     length(oR$best_iter), paste(oR$best_iter, collapse = "/"),
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
res$xgb_raw <- score_report(run, oR$score, y, p_bar, label = "xgb_raw",
                            n_bins = nb, n_boot = nboot, seed = cfg$seed,
                            title = "Out-of-fold risk ordering — XGBoost on raw features")

if (do_feat) {
  t0 <- Sys.time()
  oF <- xgb_oof_perfold(feat_fn, y, fold_of, cfg, seed = cfg$seed, group = grp)
  log_msg(run, sprintf("xgb_feat: %d folds, best_iter %s, %.1f s",
                       length(oF$best_iter), paste(oF$best_iter, collapse = "/"),
                       as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res$xgb_feat <- score_report(run, oF$score, y, p_bar, label = "xgb_feat",
                               n_bins = nb, n_boot = nboot, seed = cfg$seed,
                               title = "Out-of-fold risk ordering - XGBoost on layer-1 covariates")
  res <- res[c("llr_sum", "xgb_l", "xgb_feat", "xgb_raw")]
}

summ <- do.call(rbind, lapply(res, function(z) z$summary))
save_table(run, summ, "score_summary")

# --- the frozen boosters, for eICU -----------------------------------------
# Fitted on ALL training rows, which is what an external site must be scored
# with. Saved here so run/external.R can load and apply without refitting.
mL <- xgb_fit_full(M_full, y, cfg, seed = cfg$seed, group = grp)
mR <- xgb_fit_full(X_raw,  y, cfg, seed = cfg$seed, group = grp)
mF <- if (do_feat) xgb_fit_full(group = grp,
  xgb_design_feat(tabs, cfg, priors, stay_ids, role = "final"), y, cfg, seed = cfg$seed) else NULL
save_object(run, list(xgb_l = mL, xgb_raw = mR, xgb_feat = mF), "xgb_models")
save_table(run, xgb_importance(mR, 25L), "importance_raw", subdir = "diagnostics")
save_table(run, xgb_importance(mL, 25L), "importance_l",   subdir = "diagnostics")
if (do_feat) save_table(run, xgb_importance(mF, 25L), "importance_feat", subdir = "diagnostics")

# The DESIGN of every cell, one row per column offered, with what the booster
# did with it. This is the authoritative feature list -- `xgb.importance()`
# alone omits any column that was never split on, which is exactly the column
# you want to know about.
mods <- list(xgb_l = mL, xgb_raw = mR)
if (do_feat) mods$xgb_feat <- mF
ft <- do.call(rbind, lapply(names(mods), function(k) xgb_feature_table(mods[[k]], k)))
save_table(run, ft, "feature_manifest", subdir = "diagnostics")
gg <- do.call(rbind, lapply(split(ft, ft$design), xgb_group_gain))
save_table(run, gg, "group_gain", subdir = "diagnostics")

# --- report -----------------------------------------------------------------
cat("\n=== discrimination: the three cells ===\n\n")
print(summ[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
               "auprc", "auprc_lo", "auprc_hi", "auprc_lift")], row.names = FALSE)

cat("\n=== monotonicity of the 20-tile curve ===\n\n")
mono <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$mono)))
print(mono, row.names = FALSE)

cat("\n=== calibration slope (1.000 = the score is already in log-odds units) ===\n\n")
cal <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$cal)))
print(cal, row.names = FALSE)

cat("\n=== the three readings ===\n\n")
g <- function(a, b, f) summ[[f]][summ$label == a] - summ[[f]][summ$label == b]
cat(sprintf("  llr_sum vs xgb_raw   dAUROC %+.4f  dAUPRC %+.4f   does the method hold up?\n",
            g("llr_sum", "xgb_raw", "auroc"), g("llr_sum", "xgb_raw", "auprc")))
cat(sprintf("  xgb_l   vs llr_sum   dAUROC %+.4f  dAUPRC %+.4f   what a nonlinear aggregator adds\n",
            g("xgb_l", "llr_sum", "auroc"), g("xgb_l", "llr_sum", "auprc")))
cat(sprintf("  xgb_l   vs xgb_raw   dAUROC %+.4f  dAUPRC %+.4f   what the L representation costs or buys\n",
            g("xgb_l", "xgb_raw", "auroc"), g("xgb_l", "xgb_raw", "auprc")))

if (do_feat) {
  cat("\n=== the ladder: where the gap to a strong learner actually goes ===\n\n")
  cat(sprintf("  raw  -> feat   dAUROC %+.4f  dAUPRC %+.4f   the COVARIATE construction\n",
              g("xgb_feat", "xgb_raw", "auroc"), g("xgb_feat", "xgb_raw", "auprc")))
  cat(sprintf("  feat -> L      dAUROC %+.4f  dAUPRC %+.4f   the per-signal COLLAPSE to 19 scalars\n",
              g("xgb_l", "xgb_feat", "auroc"), g("xgb_l", "xgb_feat", "auprc")))
  cat(sprintf("  L    -> sum    dAUROC %+.4f  dAUPRC %+.4f   linear AGGREGATION (the Sigma-inverse case)\n",
              g("llr_sum", "xgb_l", "auroc"), g("llr_sum", "xgb_l", "auprc")))
  cat(sprintf("\n  These three sum to the headline llr_sum - xgb_raw = %+.4f AUROC.\n",
              g("llr_sum", "xgb_raw", "auroc")))
  cat(sprintf("  feat design: %d columns, %.1f%% NA. Unmeasured stays are NA, not 0:\n",
              ncol(X_feat1), 100 * mean(is.na(X_feat1))))
  cat("  0 is neutral for L (a log-odds contribution) but asserts a physiological\n")
  cat("  fact in covariate space. So raw -> feat isolates the covariates and NOT\n")
  cat("  the missingness channel, which --no-missing sizes at ~0.001 AUROC.\n")
}
cat(sprintf("\n  raw design: %d columns, %.1f%% NA, missing_as_evidence = %s\n",
            ncol(X_raw), 100 * mean(is.na(X_raw)), !no_miss))
cat("  A margin for xgb_raw is partly the missingness channel the LLR design\n")
cat("  refuses on principle. Re-run with --no-missing to size it.\n")

cat("\n=== 20-tile observed mortality, xgb_raw ===\n\n")
print(res$xgb_raw$bins[, c("bin", "n", "deaths", "obs_rate", "lo", "hi", "mean_score")],
      row.names = FALSE)

cat("\n=== design coverage: what each booster was given, and what it used ===\n\n")
cov <- do.call(rbind, lapply(split(ft, ft$design), function(z) data.frame(
  design = z$design[1], n_cols = nrow(z), n_used = sum(z$used),
  n_unused = sum(!z$used), n_groups = length(unique(stats::na.omit(z$group))),
  gain_top1 = round(max(z$gain), 4), gain_top5 = round(sum(sort(z$gain, decreasing = TRUE)[1:5]), 4))))
print(cov, row.names = FALSE)
cat("\n  A column with used = FALSE was offered and never split on. Gain is an\n")
cat("  orientation check, NOT an attribution result -- under correlated columns\n")
cat("  it is arbitrary WITHIN a group, which is why the rollup below exists.\n")

for (d in intersect(c("xgb_feat", "xgb_raw", "xgb_l"), names(mods))) {
  z <- gg[gg$design == d, , drop = FALSE]
  cat(sprintf("\n=== %s: gain by term block (%d groups, top 12) ===\n\n", d, nrow(z)))
  print(utils::head(transform(z, gain = round(gain, 4))[, -1], 12L), row.names = FALSE)
  u <- ft[ft$design == d & !ft$used, "feature"]
  if (length(u)) cat(sprintf("\n  never split on (%d): %s\n", length(u), paste(u, collapse = ", ")))
}


# Does the ATTRIBUTION ordering survive each rung, and not just the AUROC?
# Gain is an orientation check and not the 2x2 deliverable (that needs TreeSHAP),
# but it is free here and it is a second, independent reading of the same ladder.
if (do_feat) {
  w <- reshape(gg[, c("design", "group", "gain")], idvar = "group",
               timevar = "design", direction = "wide")
  names(w) <- sub("^gain[.]", "", names(w))
  cc <- do.call(rbind, lapply(list(c("xgb_raw", "xgb_feat"), c("xgb_feat", "xgb_l"),
                                   c("xgb_raw", "xgb_l")), function(p) {
    z <- w[stats::complete.cases(w[, p]), , drop = FALSE]
    data.frame(a = p[1], b = p[2], n_shared = nrow(z),
               spearman = round(stats::cor(z[[p[1]]], z[[p[2]]], method = "spearman"), 3),
               pearson  = round(stats::cor(z[[p[1]]], z[[p[2]]], method = "pearson"), 3))
  }))
  save_table(run, cc, "group_gain_concordance", subdir = "diagnostics")
  cat("\n=== does the term-block ORDERING survive each rung? ===\n\n")
  print(cc, row.names = FALSE)
  cat("\n  Gain, not SHAP: orientation only, and an agreement number is not\n")
  cat("  interpretable without the within-design bootstrap floor of v2_state SS4.3(a).\n")
}

finalize_run(run, extra = list(
  source_run = basename(src_dir),
  n_stays = length(stay_ids),
  missing_as_evidence = !no_miss,
  auroc_llr_sum = summ$auroc[summ$label == "llr_sum"],
  auroc_xgb_l   = summ$auroc[summ$label == "xgb_l"],
  auroc_xgb_raw = summ$auroc[summ$label == "xgb_raw"],
  auroc_xgb_feat = if (do_feat) summ$auroc[summ$label == "xgb_feat"] else NA_real_,
  feat_cell = do_feat))
cat(sprintf("\n  run directory: %s\n", run$path))
