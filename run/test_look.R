# run/test_look.R ------------------------------------------------------------
# THE SINGLE TEST LOOK. Loads a bundle, scores MIMIC-IV TEST, and fits nothing.
#
# WHAT THIS SPENDS, AND WHY IT NOW COVERS EVERY ARM. CLAUDE.md reserves one look
# at the 20% test split, re-designated 2026-08-31 to the APACHE II / SOFA
# comparison. Decided 2026-09-01: that one look scores every arm in a single
# pass rather than the severity cells alone, because the eICU transport claim
# needs a MIMIC comparator computed on held-out data through the identical code
# path. Comparing eICU against MIMIC OUT-OF-FOLD TRAIN would confound transport
# with the difference between two kinds of held-out, and the confound runs in
# the direction that flatters the method.
#
# The APACHE II / SOFA cells remain the PRE-REGISTERED CONFIRMATORY result. The
# other arms are descriptive reference points for the external comparison, and
# must be reported as such.
#
# IT IS A SEPARATE SCRIPT FROM run/internal.R ON PURPOSE. The look must be spent
# deliberately. If scoring test were a target, anyone typing `tar_make()` to
# rebuild a figure would spend it.
#
# NOTHING IS FITTED AND NOTHING IS RE-DERIVED (hard rule 8). alpha, p_bar, the
# delta and lambda parameters, the 43 GAMs, the three boosters and the reporting
# cut points all come out of the bundle. `assign_folds()` is called, and that is
# NOT a fit: it is a deterministic partition of the cohort from `seed` and
# `split`, both of which are read from the BUNDLE'S frozen design rather than
# from config/config.yml, so the test set here is provably the same test set the
# training run held out.
#
# AGGREGATES ONLY (hard rule 1). Scores and L's go to the run directory as .rds.
#
#   Rscript run/test_look.R                      # bundle from config/internal.yml
#   Rscript run/test_look.R out/runs/internal_.../bundle.qs2
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
rc   <- yaml::read_yaml("config/internal.yml")
tl   <- rc$test_look %||% list()

bundle_path <- if (length(args) && nzchar(args[1])) args[1] else tl$bundle
if (is.null(bundle_path) || !nzchar(bundle_path)) {
  stop("no bundle. Set test_look.bundle in config/internal.yml, or pass a path.\n",
       "It is an EXPLICIT path on purpose: a published result must not change ",
       "because a newer internal run appeared.", call. = FALSE)
}

# The local config is loaded ONLY to supply data paths and run-level settings,
# and to be compared against the bundle's frozen design by verify_bundle().
cfg_local <- load_config(rc$config)
bundle    <- load_bundle(bundle_path, cfg = cfg_local, strict = TRUE)

# THE DESIGN COMES FROM THE BUNDLE. Everything below reads `cfg`, which is the
# frozen design with MIMIC's paths grafted on -- never `cfg_local`.
cfg  <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
tabs <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)

folds    <- assign_folds(tabs$cohort, cfg)
test_ids <- folds$stay_id[folds$split == "test"]
y_test   <- as.integer(tabs$cohort$mortality[match(as.character(test_ids),
                                                   as.character(tabs$cohort$stay_id))])
if (anyNA(y_test)) stop("test_look: a test stay has no cohort row", call. = FALSE)
p_bar_test <- mean(y_test)
# The bootstrap unit (statistical review S4): the patient of every test stay,
# resolved through the same accessor the folds use. Never printed.
grp_test <- patient_group_of(tabs$cohort, cfg, test_ids)

arms <- as.character(unlist(tl$arms %||% BUNDLE_ARMS))
nb   <- cfg_local$metrics$n_bins %||% 20L
nbt  <- cfg_local$metrics$n_boot %||% 200L

run <- new_run("test", cfg_local, note = sprintf(
  "THE SINGLE TEST LOOK. bundle %s; %d test stays; %d arms.",
  basename(dirname(bundle_path)), length(test_ids), length(arms)))
log_msg(run, "bundle: ", bundle_path)
log_msg(run, sprintf("test set: %d stays, %d deaths, event rate %.4f",
                     length(test_ids), sum(y_test), p_bar_test))
save_table(run, verify_bundle(bundle, cfg = cfg_local, strict = FALSE),
           "bundle_checks", subdir = "diagnostics")
save_table(run, bundle_summary(bundle), "bundle_contents", subdir = "diagnostics")

# --- apply ------------------------------------------------------------------
# The identical call run/external.R makes. Every number below is produced by
# apply_bundle(); the two runners differ only in which rows they hand it.
t0 <- Sys.time()
ap <- apply_bundle(bundle, tabs, cfg, test_ids, arms = arms, verbose = FALSE)
log_msg(run, sprintf("apply_bundle: %.1f min",
                     as.numeric(difftime(Sys.time(), t0, units = "mins"))))

save_object(run, ap$l_long, "l_test")
save_object(run, ap$scores, "scores_test")
save_table(run, ap$coverage, "layer1_coverage", subdir = "diagnostics")
if (!is.null(ap$arm_table)) save_table(run, ap$arm_table, "arm_design", subdir = "diagnostics")

# --- score ------------------------------------------------------------------
# Twice, on purpose, because the two binnings answer different questions:
#
#   frozen  the cut points MIMIC TRAIN was binned on. What the external run will
#           also use, so the three 20-bin curves are directly comparable.
#   self    test's own quantiles. The within-site reading, and the one that is
#           monotone by construction if the score orders anything at all.
#
# Reporting only the second would hide a calibration shift; only the first would
# make a shift look like a discrimination failure.
frozen <- if (isTRUE(tl$frozen_bins %||% TRUE)) bundle$cutpoints else NULL
sc_f <- score_arms(run, ap$scores, y_test, p_bar_test, breaks = frozen,
                   suffix = "_frozen", n_bins = nb, n_boot = nbt, seed = cfg$seed,
                   group = grp_test)
sc_s <- score_arms(run, ap$scores, y_test, p_bar_test, breaks = NULL,
                   suffix = "_self",   n_bins = nb, n_boot = nbt, seed = cfg$seed,
                   group = grp_test)

summ <- rbind(cbind(binning = "frozen", sc_f$summary),
              cbind(binning = "self",   sc_s$summary))
save_table(run, summ, "score_summary")

ct <- arm_contrasts(ap$scores, y_test, LADDER_CONTRASTS, n_boot = nbt, seed = cfg$seed,
                    group = grp_test)
save_table(run, ct, "arm_contrasts")

# --- the transport comparator ----------------------------------------------
# Test against the training site's own out-of-fold numbers, which the bundle
# carries. This is the MIMIC half of the eICU comparison, and computing it here
# rather than at eICU is what keeps the external run from having to re-read an
# internal run directory.
tr <- bundle$train_ref
delta <- merge(
  data.frame(label = sc_s$summary$label, auroc_test = sc_s$summary$auroc,
             auprc_test = sc_s$summary$auprc, stringsAsFactors = FALSE),
  data.frame(label = paste0(tr$arms$label, "_self"),
             auroc_train_oof = tr$arms$auroc, auprc_train_oof = tr$arms$auprc,
             stringsAsFactors = FALSE),
  by = "label", all = TRUE)
delta$d_auroc <- round(delta$auroc_test - delta$auroc_train_oof, 5)
delta$d_auprc <- round(delta$auprc_test - delta$auprc_train_oof, 5)
save_table(run, delta, "train_oof_vs_test")

# --- the severity arm -------------------------------------------------------
# APACHE II and SOFA, through `severity_arm()` -- the SAME function
# run/external.R calls, so the two sites differ only in which rows they hand it.
# Nothing is fitted: the point tables are frozen in code and the recalibration
# intercept and slope come out of the bundle (hard rule 8).
#
# The coverage floors restrict every cell, ours included, so `llr_sum` appears
# twice in this run on two different row sets. The gap between them is the size
# of the selection the floors introduce and is reported rather than assumed away.
#
# FOUR STATES, NOT TWO (external runner review E1, 2026-09-09, applied here for
# parity). The arm is `disabled` by config, `unavailable` when the bundle has
# no severity slot or the site configures no severity table -- both declared
# conditions, recorded in the manifest -- or it runs. If it runs and raises,
# the error is no longer swallowed: `run_stage()` stamps the manifest `failed`
# with the stage and reason and re-raises, so a manifest that says `complete`
# once again means every requested arm finished.
sa <- NULL
sev_status <- {
  if (!isTRUE(tl$severity$enabled %||% TRUE)) "disabled"
  else if (is.null(bundle$severity)) "unavailable: bundle carries no severity slot"
  else if (is.null(tabs$severity)) "unavailable: no severity table configured for this site"
  else "requested"
}
if (sev_status == "requested") {
  sa <- run_stage(run, "severity_arm",
    severity_arm(run, bundle, tabs, cfg, test_ids, y_test, ap$scores, group = grp_test,
                 l_full = ap$l_mats$full, domains = bundle$domains,
                 n_bins = nb, n_boot = nbt, seed = cfg$seed, verbose = FALSE),
    extra = list(bundle_path = bundle_path, n_test = length(test_ids),
                 n_events = sum(y_test), arms = as.list(arms)))
  sev_status <- "complete"
} else {
  log_msg(run, "severity arm: ", sev_status)
}


# --- report -----------------------------------------------------------------
cat("\n=== THE SINGLE TEST LOOK: MIMIC-IV held-out test ===\n\n")
cat(sprintf("  bundle      %s\n", bundle_path))
cat(sprintf("  test set    %d stays, %d deaths (%.2f%%)\n",
            length(test_ids), sum(y_test), 100 * p_bar_test))
cat(sprintf("  arms        %s\n\n", paste(arms, collapse = ", ")))

cat("=== discrimination (self-binned; AUROC and AUPRC do not depend on bins) ===\n\n")
print(sc_s$summary[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
                       "auprc", "auprc_lift", "spearman", "cal_slope")],
      row.names = FALSE)

cat("\n=== the same scores binned on MIMIC-train cut points ===\n")
cat("  A calibration shift shows here and nowhere else: AUROC is rank-based and\n")
cat("  cannot see one. `cal_slope` away from 1.000 is the number to read.\n\n")
print(sc_f$summary[, c("label", "n", "spearman", "rate_ratio", "cal_slope")],
      row.names = FALSE)

cat("\n=== paired contrasts, identical rows ===\n")
cat("  Interval and `auroc_p` are the paired bootstrap resampled by `boot_unit`\n")
cat("  (patients). DeLong is the observation-level iid reference and is printed\n")
cat("  separately below; the two rest on different assumptions.\n\n")
print(ct[, c("a", "b", "d_auroc", "auroc_lo", "auroc_hi", "auroc_p", "boot_unit",
             "d_auprc", "auprc_lo", "auprc_hi", "auprc_p")], row.names = FALSE)
cat("\n  DeLong, iid sensitivity (not cluster-adjusted):\n\n")
print(ct[, c("a", "b", "d_auroc", "delong_se", "delong_z", "delong_p")], row.names = FALSE)

cat("\n=== the ladder, and whether it sums ===\n\n")
g <- function(a, b) ct$d_auroc[ct$a == a & ct$b == b]
rungs <- c(covariates = g("xgb_feat", "xgb_raw"),
           collapse   = g("xgb_l", "xgb_feat"),
           aggregation = g("llr_sum", "xgb_l"))
cat(sprintf("  covariate construction   %+.4f\n", rungs[1]))
cat(sprintf("  per-signal collapse      %+.4f\n", rungs[2]))
cat(sprintf("  linear aggregation       %+.4f\n", rungs[3]))
cat(sprintf("  --------------------------------\n  sum                      %+.4f\n", sum(rungs)))
cat(sprintf("  llr_sum - xgb_raw        %+.4f   (these must agree)\n",
            g("llr_sum", "xgb_raw")))

cat("\n=== train out-of-fold versus test ===\n")
cat("  A large drop is overfitting to the training folds; a small one is what\n")
cat("  the eICU comparison is read against.\n\n")
print(delta, row.names = FALSE)

if (!is.null(sa)) report_severity_arm(sa, "MIMIC-IV held-out test")

cat("\n=== layer-1 coverage at test ===\n\n")
print(ap$coverage[, c("signal", "model", "n_scored", "frac_scored",
                      "p_bar_train", "l_mean", "l_sd")], row.names = FALSE)

finalize_run(run, extra = list(
  bundle_path = bundle_path,
  n_test = length(test_ids), n_events = sum(y_test),
  arms = as.list(arms),
  auroc = as.list(stats::setNames(sc_s$summary$auroc, sc_s$summary$label)),
  severity_scored = !is.null(sa),
  severity_status = sev_status,
  severity_n_kept = if (is.null(sa)) NA_integer_ else sum(sa$keep),
  note_test_look = paste(
    "The single test look, spent on all arms. APACHE II / SOFA are the",
    "pre-registered confirmatory cells; the rest are descriptive reference",
    "points for the eICU comparison.")))

cat(sprintf("\n  run directory: %s\n\n", run$path))
