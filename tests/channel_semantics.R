# tests/channel_semantics.R ---------------------------------------------------
# IS +0.2 OF EVIDENCE FROM ONE SIGNAL THE SAME RISK AS +0.2 FROM ANOTHER?
#
# Three readings of that question, on the joint-model evidence matrix
# (`L_full`) at MIMIC-IV train (five-fold out-of-fold) and at eICU (the frozen
# bundle applied once), plus the layer-2 constructions that the question led
# to. Written 2026-09-13/14 while planning the Methods section; the numbers
# are transcribed in docs/channel_calibration_20260913.md.
#
#   A. MARGINAL CALIBRATION PER CHANNEL. Among stays where a signal was
#      measured, logit(mortality) regressed on that signal's L. Slope 1 and
#      intercept equal to the measured-subpopulation prior mean the channel's
#      evidence is on the common log-odds unit. This is the empirical basis
#      for summing the nineteen values at all.
#   B. RARITY WITHIN CHANNEL. Quantiles of L per signal, and the percentile at
#      which +0.2 and +0.5 sit. Same risk increment, different abnormality; a
#      display quantity for the patient card, never summed.
#   C. INCREMENTAL EVIDENCE: LAYER-2 CONSTRUCTIONS. (i) The Gaussian
#      shared-covariance weights w = S^{-1} d with S the pooled WITHIN-CLASS
#      covariance on measured entries; (ii) the pattern-dependent version that
#      uses the sub-block S_MM for the signals a patient actually has, which is
#      the Gaussian marginalised over unmeasured coordinates and never treats
#      an unmeasured signal as an observation of zero; (iii) a logistic stack
#      on the nineteen OOF values, cross-fitted over the training folds for
#      its internal estimate and frozen from full train for eICU. Each is
#      scored against the unweighted sum at both sites, paired on patients.
#
# NOTHING HERE IS A SCORING ARM. S, d and the stack coefficients are estimated
# on MIMIC train and applied to eICU as frozen constants, which is the
# fit/apply direction the bundle enforces; the eICU-own covariance row is a
# diagnostic ceiling and is labelled as such. No bundle is written.
#
# THE MIMIC ROWS FOR (i) AND (ii) ARE IN-SAMPLE for S and d (19 x 19 second
# moments on 41,000 stays; optimism is small but non-zero and is stated in the
# table). The logistic stack's MIMIC row is cross-fitted and honest.
#
# THE L MATRICES COME FROM `l_matrix()`, not from a hand pivot of the long
# table. The saved long tables hold only fitted specs, so `full` has no rows
# for the seven unpaired signals -- it is an alias of `meas` there -- and a
# pivot on `model == "full"` silently zeroes lactate, bilirubin, BUN, sodium,
# bicarbonate, WBC and temperature. That is exactly what happened in the
# scratch version of this check on 2026-09-13 and it turned a null transport
# result into a spurious 0.047 loss. `l_matrix()` walks the alias through
# `resolve_spec_source()`, and the rebuilt sum is asserted equal to the run's
# stored `llr_sum` before any number is computed.
#
# Aggregates only (hard rule 1): every table is per signal, per pattern
# class, or per arm. Scores and L are never printed or saved.
#
#   Rscript tests/channel_semantics.R
#   Rscript tests/channel_semantics.R --internal out/runs/internal_<id> \
#                                     --external out/runs/external_<id> --n-boot 200
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(qs2); library(Matrix)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
INTERNAL <- .opt("--internal", latest_run("internal", require_complete = TRUE))
EXTERNAL <- .opt("--external", latest_run("external", require_complete = TRUE))
N_BOOT   <- as.integer(.opt("--n-boot", "200"))
SEED     <- 5333L
if (is.null(INTERNAL) || is.null(EXTERNAL)) stop("no complete internal/external run found")

# --- provenance: the external run must have applied the internal run's bundle --
xcfg   <- yaml::read_yaml("config/external.yml")
bpath  <- as.character(cfg_req(xcfg, "bundle"))
xman   <- read_manifest(EXTERNAL)
if (!identical(normalizePath(dirname(bpath)), normalizePath(INTERNAL))) {
  stop("config/external.yml names bundle ", bpath, " but --internal is ", INTERNAL,
       ": the eICU L must come from the same bundle as the MIMIC OOF L", call. = FALSE)
}
cat(sprintf("\n=== channel semantics: MIMIC %s, eICU %s ===\n\n",
            basename(INTERNAL), basename(EXTERNAL)))

# --- MIMIC train, out of fold ------------------------------------------------
cfg_m  <- load_config("config/config.yml")
tabs_m <- load_tables(cfg_m$paths$mimiciv, cfg_m, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs_m$cohort, cfg_m)
tr     <- folds$stay_id[folds$split == "train"]
fold_m <- folds$fold[folds$split == "train"]
y_m    <- as.integer(tabs_m$cohort$mortality[match(tr, tabs_m$cohort$stay_id)])
grp_m  <- patient_group_of(tabs_m$cohort, cfg_m, tr)
l_oof  <- readRDS(file.path(INTERNAL, "tables", "l_oof.rds"))
Mn     <- l_matrix(l_oof, "full", tabs_m, cfg_m, tr, fill = "na")
sig    <- colnames(Mn)
M0     <- Mn; M0[is.na(M0)] <- 0
oof    <- readRDS(file.path(INTERNAL, "tables", "oof_scores.rds"))$llr_sum
if (!is.null(names(oof))) oof <- oof[as.character(tr)]
if (length(oof) != nrow(M0) || max(abs(rowSums(M0) - as.numeric(oof))) > 1e-9) {
  stop("MIMIC: rebuilt L_full does not reproduce the run's stored llr_sum", call. = FALSE)
}

# --- eICU, frozen bundle -----------------------------------------------------
bundle <- load_bundle(bpath, verbose = FALSE)
cfg_e  <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
tabs_e <- load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE)
ids_e  <- sort(unique(tabs_e$cohort$stay_id))
y_e    <- as.integer(tabs_e$cohort$mortality[match(ids_e, tabs_e$cohort$stay_id)])
grp_e  <- patient_group_of(tabs_e$cohort, cfg_e, ids_e)
l_e    <- readRDS(file.path(EXTERNAL, "tables", "l_eicu.rds"))
En     <- l_matrix(l_e, "full", tabs_e, cfg_e, ids_e, fill = "na")
stopifnot(identical(colnames(En), sig))
E0     <- En; E0[is.na(E0)] <- 0
sc_e   <- readRDS(file.path(EXTERNAL, "tables", "scores_eicu.rds"))$llr_sum
sc_e   <- sc_e[as.character(ids_e)]
if (anyNA(sc_e) || max(abs(rowSums(E0) - as.numeric(sc_e))) > 1e-9) {
  stop("eICU: rebuilt L_full does not reproduce the run's stored llr_sum", call. = FALSE)
}
cat(sprintf("MIMIC train %d stays, %d deaths; eICU %d stays, %d deaths; both L matrices reproduce the stored sums.\n",
            length(tr), sum(y_m), length(ids_e), sum(y_e)))

run <- new_run("chansem", cfg_m, note = sprintf(
  "per-channel calibration, rarity and layer-2 constructions; MIMIC %s, eICU %s",
  basename(INTERNAL), basename(EXTERNAL)))
prov <- data.frame(item = c("internal_run", "external_run", "bundle", "n_boot", "seed"),
                   value = c(basename(INTERNAL), basename(EXTERNAL), bpath, N_BOOT, SEED))
save_table(run, prov, "provenance", subdir = "diagnostics")

# ============ A. marginal calibration per channel ============================
calib_site <- function(M, yy, site) {
  do.call(rbind, lapply(sig, function(g) {
    m <- !is.na(M[, g])
    cal <- llr_calibration(M[m, g], yy[m], p_bar = mean(yy[m]))
    data.frame(site = site, signal = g, n_measured = sum(m),
               slope = cal$slope, slope_lo = cal$slope_lo, slope_hi = cal$slope_hi,
               intercept_vs_site_prior = round(cal$intercept - cal$expected_intercept, 4),
               p_bar_site_measured = cal$p_bar_reference, stringsAsFactors = FALSE)
  }))
}
cal <- rbind(calib_site(Mn, y_m, "mimic_oof"), calib_site(En, y_e, "eicu"))
save_table(run, cal, "channel_calibration", subdir = "diagnostics")
cat("\nA. per-channel calibration slope of logit(mortality) on L (1 = common unit):\n")
print(cal[, c("site", "signal", "n_measured", "slope", "slope_lo", "slope_hi", "intercept_vs_site_prior")],
      row.names = FALSE)

binned <- function(M, yy) {
  do.call(rbind, lapply(sig, function(g) {
    m <- !is.na(M[, g]); L <- M[m, g]; v <- yy[m]; pbar <- mean(v)
    q <- unique(stats::quantile(L, seq(0, 1, length.out = 21)))
    b <- cut(L, q, include.lowest = TRUE)
    xm <- tapply(L, b, mean); pm <- tapply(v, b, mean); n <- tapply(v, b, length)
    ok <- pm > 0 & pm < 1
    data.frame(signal = g, bin = seq_along(xm)[ok], l_mean = round(xm[ok], 4),
               n = as.integer(n[ok]), logit_obs_minus_prior = round(logit(pm[ok]) - logit(pbar), 4),
               row.names = NULL)
  }))
}
bm <- binned(Mn, y_m); be <- binned(En, y_e)
save_table(run, cbind(site = "mimic_oof", bm), "channel_calibration_bins_mimic", subdir = "diagnostics")
save_table(run, cbind(site = "eicu", be), "channel_calibration_bins_eicu", subdir = "diagnostics")
save_fig(run, "channel_calibration", width = 12, height = 6, dpi = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
cols <- grDevices::hcl.colors(length(sig), "Dark 3")
for (nm in c("MIMIC-IV train, out of fold", "eICU, frozen bundle")) {
  d <- if (grepl("MIMIC", nm)) bm else be
  rng <- range(c(d$l_mean, d$logit_obs_minus_prior))
  plot(NA, xlim = rng, ylim = rng, xlab = "evidence L (bin mean)",
       ylab = "logit(observed mortality) - logit(measured prior)", main = nm)
  abline(0, 1, col = "grey40", lwd = 2)
  for (i in seq_along(sig)) { dd <- d[d$signal == sig[i], ]
    lines(dd$l_mean, dd$logit_obs_minus_prior, col = cols[i], lwd = 1.2)
    points(dd$l_mean, dd$logit_obs_minus_prior, col = cols[i], pch = 16, cex = 0.5) }
}
legend("bottomright", legend = sig, col = cols, lwd = 1.5, cex = 0.55, ncol = 2, bg = "white")
dev.off()

# ============ B. rarity within channel =======================================
rar <- do.call(rbind, lapply(sig, function(g) {
  L <- Mn[!is.na(Mn[, g]), g]
  q <- stats::quantile(L, c(.01, .10, .50, .90, .99))
  data.frame(signal = g, n_measured = length(L),
             p01 = round(q[1], 3), p10 = round(q[2], 3), p50 = round(q[3], 3),
             p90 = round(q[4], 3), p99 = round(q[5], 3),
             pct_below_0.2 = round(100 * mean(L <= 0.2), 1),
             pct_below_0.5 = round(100 * mean(L <= 0.5), 1),
             max_abs = round(max(abs(L)), 3), row.names = NULL)
}))
save_table(run, rar, "channel_rarity_mimic", subdir = "diagnostics")
cat("\nB. rarity: quantiles of L per channel (MIMIC train, measured stays):\n")
print(rar, row.names = FALSE)

# ============ C. layer-2 constructions =======================================
within_cov <- function(M, yy) {
  n0 <- sum(yy == 0); n1 <- sum(yy == 1)
  S <- (n0 * stats::cov(M[yy == 0, ], use = "pairwise.complete.obs") +
        n1 * stats::cov(M[yy == 1, ], use = "pairwise.complete.obs")) / (n0 + n1)
  pd <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > 1e-8
  if (!pd) S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  d <- colMeans(M[yy == 1, ], na.rm = TRUE) - colMeans(M[yy == 0, ], na.rm = TRUE)
  list(S = S, d = d, was_pd = pd)
}
mm <- within_cov(Mn, y_m); ee <- within_cov(En, y_e)

# C1: does the within-class covariance itself transport?
Wm <- stats::cov2cor(mm$S); We <- stats::cov2cor(ee$S); ut <- upper.tri(Wm)
ij <- which(ut, arr.ind = TRUE)
pairs <- data.frame(a = sig[ij[, 1]], b = sig[ij[, 2]],
                    r_within_mimic = round(Wm[ut], 4), r_within_eicu = round(We[ut], 4))
pairs$abs_diff <- round(abs(pairs$r_within_mimic - pairs$r_within_eicu), 4)
pairs <- pairs[order(-pairs$abs_diff), ]
save_table(run, pairs, "within_class_correlation_pairs", subdir = "diagnostics")
cov_sum <- data.frame(
  quantity = c("cor_offdiag_across_sites", "mean_abs_diff", "max_abs_diff",
               "sd_ratio_eicu_over_mimic_min", "sd_ratio_median", "sd_ratio_max",
               "d_ratio_eicu_over_mimic_min", "d_ratio_median", "d_ratio_max",
               "S_mimic_was_pd", "S_eicu_was_pd"),
  value = c(round(stats::cor(Wm[ut], We[ut]), 4), round(mean(pairs$abs_diff), 4), max(pairs$abs_diff),
            round(range(sqrt(diag(ee$S)) / sqrt(diag(mm$S)))[1], 3),
            round(stats::median(sqrt(diag(ee$S)) / sqrt(diag(mm$S))), 3),
            round(range(sqrt(diag(ee$S)) / sqrt(diag(mm$S)))[2], 3),
            round(range(ee$d / mm$d)[1], 3), round(stats::median(ee$d / mm$d), 3),
            round(range(ee$d / mm$d)[2], 3), as.integer(mm$was_pd), as.integer(ee$was_pd)))
save_table(run, cov_sum, "within_class_covariance_transport", subdir = "diagnostics")
cat("\nC1. within-class correlation transport:\n"); print(cov_sum, row.names = FALSE)

# C2: weights. Fixed (all-19) Gaussian, pattern-dependent, and the logistic stack.
w_of <- function(m, S, d) {
  w <- rep(NA_real_, length(sig)); names(w) <- sig
  if (sum(m) >= 1L) w[m] <- solve(S[m, m, drop = FALSE], d[m])
  w
}
score_pat <- function(M, S, d) {
  mask <- !is.na(M)
  pat  <- apply(mask, 1L, function(r) paste(as.integer(r), collapse = ""))
  up   <- unique(pat)
  W    <- vapply(up, function(p) w_of(as.logical(as.integer(strsplit(p, "")[[1]])), S, d),
                 numeric(length(sig)))
  colnames(W) <- up
  X <- M; X[!mask] <- 0
  list(score = colSums(W[, pat, drop = FALSE] * t(X), na.rm = TRUE), pat = pat, W = W)
}
pm        <- score_pat(Mn, mm$S, mm$d)        # MIMIC, own S,d (in-sample second moments)
pe_frozen <- score_pat(En, mm$S, mm$d)        # eICU, frozen MIMIC S,d: the transport
pe_own    <- score_pat(En, ee$S, ee$d)        # eICU-own S,d: diagnostic ceiling, not a score
w_fix_m   <- solve(mm$S, mm$d)
w_fix_e   <- solve(ee$S, ee$d)

df_m <- data.frame(y = y_m, M0)
fit_full <- stats::glm(y ~ ., data = df_m, family = stats::binomial())
nn_fit <- function(d) {
  keep <- sig
  repeat {
    g  <- stats::glm(stats::reformulate(keep, "y"), data = d, family = stats::binomial())
    cf <- stats::coef(g)[-1]
    if (all(cf >= 0)) return(g)
    keep <- keep[cf >= 0]
  }
}
fit_nn <- nn_fit(df_m)
lin <- function(g, X) as.vector(cbind(1, X[, names(stats::coef(g))[-1], drop = FALSE]) %*% stats::coef(g))
cross_fit <- function(fitter) {
  s <- numeric(length(y_m))
  for (f in sort(unique(fold_m))) {
    trn <- fold_m != f
    g <- fitter(df_m[trn, , drop = FALSE])
    s[!trn] <- lin(g, M0[!trn, , drop = FALSE])
  }
  s
}
s_stack_cf <- cross_fit(function(d) stats::glm(y ~ ., data = d, family = stats::binomial()))
s_nn_cf    <- cross_fit(nn_fit)

co <- stats::coef(fit_full)[-1]; se <- sqrt(diag(stats::vcov(fit_full)))[-1]
nn_co <- stats::coef(fit_nn)[-1]
weights <- data.frame(
  signal = sig,
  d_mimic = round(mm$d, 4), d_eicu = round(ee$d, 4),
  w_gauss_mimic = round(w_fix_m, 4), w_gauss_mimic_scaled = round(w_fix_m / mean(w_fix_m), 3),
  w_gauss_eicu_own = round(w_fix_e, 4), w_gauss_eicu_own_scaled = round(w_fix_e / mean(w_fix_e), 3),
  stack_coef = round(co, 4), stack_se = round(se, 4), stack_z = round(co / se, 2),
  nn_stack_coef = round(ifelse(sig %in% names(nn_co), nn_co[sig], NA_real_), 4),
  row.names = NULL)
save_table(run, weights, "layer2_weights", subdir = "diagnostics")
cat("\nC2. weights per signal (Gaussian, both sites; logistic stack, MIMIC full train):\n")
print(weights[, c("signal", "w_gauss_mimic_scaled", "w_gauss_eicu_own_scaled", "stack_coef", "stack_z", "nn_stack_coef")],
      row.names = FALSE)

# pattern bookkeeping and how far frozen weights sit from eICU-own weights per pattern
tab_m <- sort(table(pm$pat), decreasing = TRUE); tab_e <- sort(table(pe_frozen$pat), decreasing = TRUE)
common <- intersect(colnames(pe_frozen$W), colnames(pe_own$W))
cw <- vapply(common, function(p) {
  a <- pe_frozen$W[, p]; b <- pe_own$W[, p]; ok <- !is.na(a)
  c(cor = if (sum(ok) >= 3L) stats::cor(a[ok], b[ok]) else NA_real_,
    sign_flips = sum(sign(a[ok]) != sign(b[ok])))
}, numeric(2))
cnt_e <- as.numeric(table(pe_frozen$pat)[common])
sw <- function(v, q) unname(stats::quantile(rep(v, cnt_e), q, na.rm = TRUE))
pat_sum <- data.frame(
  quantity = c("n_patterns_mimic", "n_patterns_eicu", "share_mimic_in_top_pattern",
               "share_mimic_top10", "share_eicu_top10", "share_eicu_in_pattern_unseen_at_mimic",
               "weight_cor_frozen_vs_own_p05", "weight_cor_frozen_vs_own_median",
               "weight_cor_frozen_vs_own_p95", "sign_flips_median_over_eicu_stays"),
  value = c(length(tab_m), length(tab_e), round(tab_m[1] / length(y_m), 4),
            round(sum(tab_m[seq_len(min(10, length(tab_m)))]) / length(y_m), 4),
            round(sum(tab_e[seq_len(min(10, length(tab_e)))]) / length(y_e), 4),
            round(mean(!(pe_frozen$pat %in% names(tab_m))), 4),
            round(sw(cw["cor", ], .05), 4), round(sw(cw["cor", ], .5), 4), round(sw(cw["cor", ], .95), 4),
            sw(cw["sign_flips", ], .5)))
save_table(run, pat_sum, "missingness_patterns", subdir = "diagnostics")
cat("\nC3. missingness patterns and per-pattern weight agreement:\n"); print(pat_sum, row.names = FALSE)

# C4: the arms, both sites
arms_m <- list(unweighted_sum = rowSums(M0), gaussian_fixed = as.vector(M0 %*% w_fix_m),
               gaussian_pattern = pm$score, logistic_stack = s_stack_cf, nonneg_stack = s_nn_cf)
arms_e <- list(unweighted_sum = rowSums(E0), gaussian_fixed = as.vector(E0 %*% w_fix_m),
               gaussian_pattern = pe_frozen$score, logistic_stack = lin(fit_full, E0),
               nonneg_stack = lin(fit_nn, E0), gaussian_pattern_eicu_own_DIAGNOSTIC = pe_own$score)
basis_m <- c(unweighted_sum = "none", gaussian_fixed = "S,d in-sample", gaussian_pattern = "S,d in-sample",
             logistic_stack = "cross-fitted over folds", nonneg_stack = "cross-fitted over folds")
row_of <- function(nm, s, yy, site, basis) {
  cal <- llr_calibration(s, yy, p_bar = mean(yy))
  data.frame(site = site, arm = nm, basis = basis, auroc = round(.auroc(s, yy), 4),
             auprc = round(.auprc(s, yy), 4), cal_slope = cal$slope,
             cal_intercept_shift = round(cal$intercept - cal$expected_intercept, 4),
             stringsAsFactors = FALSE)
}
arms <- rbind(
  do.call(rbind, lapply(names(arms_m), function(nm) row_of(nm, arms_m[[nm]], y_m, "mimic_oof", basis_m[[nm]]))),
  do.call(rbind, lapply(names(arms_e), function(nm) row_of(nm, arms_e[[nm]], y_e, "eicu",
    if (grepl("DIAGNOSTIC", nm)) "S,d estimated on eICU outcomes; NOT a transportable score" else "frozen from MIMIC full train"))))
arms$retained <- NA_real_
for (nm in names(arms_m)) {
  a_m <- arms$auroc[arms$site == "mimic_oof" & arms$arm == nm]
  a_e <- arms$auroc[arms$site == "eicu" & arms$arm == nm]
  arms$retained[arms$site == "eicu" & arms$arm == nm] <- round((a_e - 0.5) / (a_m - 0.5), 4)
}
save_table(run, arms, "layer2_arms", subdir = "diagnostics")
cat("\nC4. layer-2 arms at both sites:\n"); print(arms, row.names = FALSE)

# C5: paired contrasts against the unweighted sum, patient bootstrap, both sites
contrast <- function(nm, site) {
  A <- if (site == "mimic_oof") arms_m else arms_e
  yy <- if (site == "mimic_oof") y_m else y_e
  g  <- if (site == "mimic_oof") grp_m else grp_e
  r <- .paired_boot(A[[nm]], A$unweighted_sum, yy, metrics = list(auroc = .auroc, auprc = .auprc),
                    n_boot = N_BOOT, seed = SEED, group = g)
  data.frame(site = site, contrast = paste(nm, "minus unweighted_sum"),
             d_auroc = r$auroc$delta, auroc_lo = r$auroc$ci_lo, auroc_hi = r$auroc$ci_hi,
             d_auprc = r$auprc$delta, auprc_lo = r$auprc$ci_lo, auprc_hi = r$auprc$ci_hi,
             p_boot = r$auroc$p_boot, n_boot = r$auroc$n_boot_ok, boot_unit = r$auroc$boot_unit,
             stringsAsFactors = FALSE)
}
con <- rbind(
  do.call(rbind, lapply(c("gaussian_fixed", "gaussian_pattern", "logistic_stack", "nonneg_stack"),
                        contrast, site = "mimic_oof")),
  do.call(rbind, lapply(c("gaussian_fixed", "gaussian_pattern", "logistic_stack", "nonneg_stack"),
                        contrast, site = "eicu")))
save_table(run, con, "layer2_contrasts", subdir = "diagnostics")
cat("\nC5. paired contrasts vs the unweighted sum (patient bootstrap):\n"); print(con, row.names = FALSE)

finalize_run(run, note = "channel semantics complete; no bundle written, no scoring arm changed")
cat("\nrun:", run$path, "\n")
