# tests/coupling_attribution.R -----------------------------------------------
# WHAT DOES MODEL SPECIFICATION DO TO THE PATIENT EVIDENCE VECTOR?
#
# tests/coupling_oof.R established that the measurement x intervention
# interaction survives out of fold and is worth about +0.0038 AUROC over
# `llr_meas`. It also produced, almost incidentally, the finding that matters
# more: the SUMMED SCORE barely moves while the INDIVIDUAL ATTRIBUTIONS move a
# great deal -- more than half of some signals' cells shift by over 0.10 nats,
# and single-digit percentages of cells change the SIGN of their evidence.
#
# This script is that observation done properly. It builds every L matrix in the
# ladder, out of fold, and asks what each specification does to the distribution
# of evidence across a patient's 19 signals.
#
#   L_meas          measurement block alone
#   L_intv          intervention block alone            (identical in every arm)
#   L_full          additive joint                       <- the pipeline
#   L_full_ti_trend additive + trend x intervention only
#   L_full_ti_all   additive + every measurement x intervention cross term
#
#   L_cond          = L_full          - L_intv           <- the attribution unit
#   L_cond_ti_trend = L_full_ti_trend - L_intv
#   L_cond_ti_all   = L_full_ti_all   - L_intv
#
# WHY `L_cond` IS THE ATTRIBUTION UNIT AND `L_full` IS NOT. `L_full` bundles the
# measurement evidence with a pure treatment-propensity contrast, so a claim of
# the form "signal g contributed X nats for this patient" is, on `L_full`, partly
# a claim about who got treated. Both spaces are reported, but the cond-space
# rows are the ones a paper's attribution claim rests on.
#
# THE SIGN-FLIP METRIC IS NOW MAGNITUDE-GATED, AND THE OLD VERSION WAS WRONG TO
# REPORT UNGATED. `tests/coupling_oof.R` counted `sign(a) != sign(b)` over all
# measured cells, which scores a cell moving from +0.001 to -0.001 exactly as it
# scores one moving from +0.8 to -0.8. The first is arithmetic noise about a
# quantity that was never distinguishable from zero. Every flip rate below is
# therefore reported at a grid of thresholds tau, where a flip counts only if
# max(|a|, |b|) > tau. tau = 0 reproduces the old number and is kept so the two
# can be compared; it is an UPPER BOUND and should not be quoted alone.
#
# COMPUTE. Only the two interaction arms need fitting, at 5 folds each.
# `L_meas`, `L_intv` and `L_full` come from the targets cache and are the
# pipeline's own output. `L_intv` is untouched by construction: p(I | Y) depends
# on the treatment record alone, so no cross term can enter the `intv` model.
# The trend arm is fitted for the 9 paired signals that HAVE a trend covariate;
# `creatinine`, `platelet` and `hemoglobin` are class-gated out of `trend`
# (config `trend_classes`), so their trend arm equals the additive arm exactly
# and is assigned, never fitted -- the same treatment `layer1_jobs()` gives an
# unpaired signal's `intv`.
#
# The term lists are built from `measurement_terms()` / `intervention_terms()`
# rather than by reading a fitted object, so no throwaway base fit is needed.
# That halves the fit count against tests/coupling_oof.R.
#
# CACHING. Every L matrix is saved to the run directory, and `--reuse <dir>`
# loads them back instead of refitting. A stored arm is only reused when its
# design fingerprint -- the formula strings, the bam settings and the fold
# assignment -- matches what this run would produce, so a config edit forces a
# refit rather than silently reusing a stale matrix.
#
# Aggregates only on print (hard rule 1). The L matrices are row-level, are
# saved the way `l_oof.rds` already is, and are never printed.
#
#   Rscript tests/coupling_attribution.R
#   Rscript tests/coupling_attribution.R --reuse out/runs/coupattr_...
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args   <- commandArgs(trailingOnly = TRUE)
reuse  <- if ("--reuse" %in% args) args[which(args == "--reuse") + 1L] else NULL
K_TI   <- 5L
TAUS   <- c(0, 0.05, 0.10, 0.25, 0.50)   # nats. Declared before the run.

cfg    <- tar_read(cfg)
tabs   <- tar_read(tabs)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
y      <- as.integer(tar_read(y_train))
priors <- tar_read(priors)
Lz     <- tar_read(l_mats_zero)
domains <- load_domains("config/domains.csv")

sigs    <- unlist(cfg$signals)
paired  <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)
n_fold  <- cfg$n_folds %||% 5L
ids_ch  <- as.character(tr)
meas_ok <- measured_matrix(tabs, cfg, tr)

run <- new_run("coupattr", cfg, note =
  "attribution under specification change: additive vs trend-only vs full interaction")

# --- term construction, data-free -------------------------------------------
#
# `.ti_terms()` AND `.smooth_vars()` MOVED TO `R/14_attribution_eval.R` ON
# 2026-09-06, verbatim, as `attr_ti_terms()` and `attr_smooth_vars()`. They had
# to move because the ladder's design fingerprint is now checked by a SECOND
# script (`tests/attribution_vs_shap.R`), and a fingerprint computed by two
# copies of a function is not a fingerprint -- the copies drift and the check
# starts passing on a design mismatch. `tests/attr_metrics.R` asserts the moved
# functions reproduce the pre-move term sets exactly.

.ti_terms <- function(sg, cfg, scope = c("all", "trend")) {
  interaction_terms(sg, cfg, match.arg(scope), K_TI)
}

#' What a stored ARM must match before it may be reused.
#'
#' Narrower than the ladder key below on purpose: an arm is one interaction
#' scope's out-of-fold matrix, so the only formulas that can invalidate it are
#' its own.
.fingerprint <- function(scope) {
  list(scope = scope, k_ti = K_TI,
       formulas = vapply(paired, function(sg) {
         ti <- .ti_terms(sg, cfg, scope)
         paste(c(deparse(build_formula(sg, "full", cfg)), ti), collapse = " + ")
       }, character(1)),
       bam = bam_settings_row(cfg),
       folds = paste0(folds$fold[match(tr, folds$stay_id)], collapse = ""))
}

#' What the WHOLE LADDER must match before a consumer may read it.
#'
#' THE ARM FINGERPRINT ABOVE IS NOT ENOUGH, and the gap is what
#' `attrshap_20260906T101626` fell into. Five of the ladder's eight arms come
#' from the targets cache (`meas`, `intv`, `full` and the two `cond` arms built
#' off them), and NOTHING in `.fingerprint()` covers those: it hashes the
#' interaction formulas and the `full` formula, so a change to `bam.gamma`
#' shows up, but a change to a `meas` whitelist or a `smooth_k` override on an
#' unpaired signal would not. `attr_design_key()` covers all 43 fitted specs.
#' Stamped onto `l_oof_ladder.rds` so a consumer can refuse a stale ladder
#' rather than silently mixing it with a current SHAP floor.
.ladder_key <- function() {
  attr_design_key(cfg, folds$fold[match(tr, folds$stay_id)], K_TI)
}

# --- the out-of-fold interaction fits ---------------------------------------

fit_arm <- function(scope) {
  M <- Lz$full                       # unpaired signals carry over unchanged
  for (sg in paired) {
    ti <- .ti_terms(sg, cfg, scope)
    if (!length(ti)) {               # assigned, never fitted
      cat(sprintf("  %-18s %-6s  no cross term exists; additive column assigned\n",
                  sg, scope))
      next
    }
    f0 <- build_formula(sg, "full", cfg)
    f1 <- stats::update(f0, stats::as.formula(paste(". ~ . +", paste(ti, collapse = " + "))))
    col <- rep(NA_real_, length(tr))
    for (k in seq_len(n_fold)) {
      pri <- priors_for(priors, sg, "oof", fold = k)
      ji  <- job_ids("oof", k, folds)
      d_fit <- signal_frame(sg, "full", tabs, cfg, pri, stay_ids = ji$fit_ids)
      d_prd <- signal_frame(sg, "full", tabs, cfg, pri,
                            stay_ids = ji$predict_ids, stage = "predict")
      b <- try(.bam_fit(f1, d_fit, cfg), silent = TRUE)
      if (inherits(b, "try-error")) {
        cat(sprintf("  %-18s %-6s fold %d FAILED\n", sg, scope, k)); next
      }
      eta <- as.numeric(stats::predict(b, newdata = d_prd, type = "link", discrete = FALSE))
      col[match(as.character(d_prd$stay_id), ids_ch)] <- eta - logit(pri$p_bar)
    }
    col[is.na(col)] <- Lz$full[is.na(col), sg]
    M[, sg] <- col
    cat(sprintf("  %-18s %-6s %3d ti term(s)   cor with additive = %.4f\n",
                sg, scope, length(ti),
                stats::cor(Lz$full[meas_ok[, sg], sg], col[meas_ok[, sg]])))
  }
  attr(M, "fingerprint") <- .fingerprint(scope)
  M
}

load_or_fit <- function(scope) {
  nm <- paste0("l_oof_ti_", scope)
  if (!is.null(reuse)) {
    p <- file.path(reuse, "tables", paste0(nm, ".rds"))
    if (file.exists(p)) {
      M <- readRDS(p)
      if (identical(attr(M, "fingerprint"), .fingerprint(scope)) &&
          identical(dimnames(M), dimnames(Lz$full))) {
        cat(sprintf("  reused %s from %s\n", nm, basename(reuse)))
        return(M)
      }
      cat(sprintf("  %s in %s does NOT match this design; refitting\n",
                  nm, basename(reuse)))
    }
  }
  fit_arm(scope)
}

cat("\n=== building the interaction arms, out of fold ===\n\n")
t0 <- start_timer()
L_ti_all   <- load_or_fit("all")
L_ti_trend <- load_or_fit("trend")
save_object(run, L_ti_all,   "l_oof_ti_all")
save_object(run, L_ti_trend, "l_oof_ti_trend")

# --- the eight matrices -----------------------------------------------------
A <- list(
  meas          = Lz$meas,
  intv          = Lz$intv,
  full          = Lz$full,
  full_ti_trend = L_ti_trend,
  full_ti_all   = L_ti_all,
  cond          = Lz$full          - Lz$intv,
  cond_ti_trend = L_ti_trend       - Lz$intv,
  cond_ti_all   = L_ti_all         - Lz$intv)
# `design_key` travels WITH the ladder so a consumer can check what produced it.
# Storing it beside the arms rather than only in the manifest matters: the
# manifest is a separate file that a consumer would have to know to read, and
# the failure this guards against is a consumer that reads only the .rds.
save_object(run, list(arms = A, stay_id = tr, design_key = .ladder_key(),
                      k_ti = K_TI), "l_oof_ladder")

# --- 1. what each arm is worth, so the trade is visible ---------------------
cat("\n=== 1. discrimination of every arm, out of fold ===\n\n")
cat(sprintf("%-16s %9s %9s %9s %9s\n", "arm", "AUROC", "AUPRC", "cal_slope", "top1%"))
tp1 <- function(x) { k <- ceiling(0.01 * length(x)); mean(y[order(-x)][seq_len(k)]) }
S <- list()
for (a in names(A)) {
  s <- rowSums(A[[a]])
  cal <- llr_calibration(s, y, mean(y))
  S[[a]] <- data.frame(arm = a, auroc = round(.auroc(s, y), 5),
                       auprc = round(.auprc(s, y), 5), cal_slope = cal$slope,
                       top1 = round(tp1(s), 5), stringsAsFactors = FALSE)
  cat(sprintf("%-16s %9.5f %9.5f %9.4f %9.4f\n", a, S[[a]]$auroc, S[[a]]$auprc,
              S[[a]]$cal_slope, S[[a]]$top1))
}
save_table(run, do.call(rbind, S), "ladder_arms", subdir = "diagnostics")

# --- the comparisons this script is about -----------------------------------
CMP <- list(
  c("cond",          "cond_ti_all"),    # specification uncertainty, full
  c("cond",          "cond_ti_trend"),  # specification uncertainty, trend only
  c("cond_ti_trend", "cond_ti_all"),    # how much is NOT trend
  c("meas",          "cond"),           # what conditioning itself does
  c("meas",          "cond_ti_all"),    # conditioning with interaction
  c("full",          "full_ti_all"),    # the joint-space version
  c("full",          "full_ti_trend"))

# --- 2. magnitude-gated sign flip -------------------------------------------
cat("\n=== 2. sign flip, gated on magnitude ===\n\n")
cat("    A cell flips when the two arms disagree about the DIRECTION of that\n")
cat("    patient's evidence from that signal. Gated: the flip counts only when\n")
cat("    max(|a|,|b|) exceeds tau nats, so a reversal of a quantity that was\n")
cat("    never distinguishable from zero is not scored as a disagreement.\n")
cat("    tau = 0 is the ungated number tests/coupling_oof.R reported. It is an\n")
cat("    UPPER BOUND and should not be quoted alone.\n\n")
cat(sprintf("%-16s %-16s %8s", "arm a", "arm b", "cells"))
for (tt in TAUS) cat(sprintf("  tau=%.2f", tt))
cat("\n")
F1 <- list()
for (p in CMP) {
  a <- A[[p[1]]]; b <- A[[p[2]]]
  ok <- meas_ok
  fl <- sign(a) != sign(b)
  mx <- pmax(abs(a), abs(b))
  r <- vapply(TAUS, function(tt) mean(fl[ok] & mx[ok] > tt), numeric(1))
  F1[[length(F1) + 1L]] <- data.frame(a = p[1], b = p[2], n_cells = sum(ok),
                                      tau = TAUS, flip = round(r, 5),
                                      stringsAsFactors = FALSE)
  cat(sprintf("%-16s %-16s %8d", p[1], p[2], sum(ok)))
  for (v in r) cat(sprintf("  %8.4f", v))
  cat("\n")
}
save_table(run, do.call(rbind, F1), "sign_flip_gated", subdir = "diagnostics")

cat("\n--- per signal, at tau = 0 and tau = 0.10, for cond vs cond_ti_all ---\n\n")
a <- A$cond; b <- A$cond_ti_all
fl <- sign(a) != sign(b); mx <- pmax(abs(a), abs(b))
F2 <- list()
cat(sprintf("%-18s %9s %9s %9s %10s %10s\n", "signal", "flip_t0", "flip_t10",
            "flip_t25", "med|d|", "med|L| flip"))
for (sg in sigs) {
  ok <- meas_ok[, sg]
  if (!any(ok)) next
  d <- b[ok, sg] - a[ok, sg]
  fs <- fl[ok, sg]; ms <- mx[ok, sg]
  F2[[sg]] <- data.frame(arm_a = "cond", arm_b = "cond_ti_all",
    signal = sg, n_cells = sum(ok),
    flip_t0 = round(mean(fs), 5), flip_t10 = round(mean(fs & ms > 0.10), 5),
    flip_t25 = round(mean(fs & ms > 0.25), 5),
    med_abs_delta = round(stats::median(abs(d)), 5),
    med_mag_at_flip = round(if (any(fs)) stats::median(ms[fs]) else NA_real_, 5),
    stringsAsFactors = FALSE)
  cat(sprintf("%-18s %9.4f %9.4f %9.4f %10.4f %10.4f\n", sg, F2[[sg]]$flip_t0,
              F2[[sg]]$flip_t10, F2[[sg]]$flip_t25, F2[[sg]]$med_abs_delta,
              F2[[sg]]$med_mag_at_flip))
}
save_table(run, do.call(rbind, F2), "sign_flip_by_signal", subdir = "diagnostics")

# --- 3. within-patient attribution ranking ----------------------------------
#
# THE METRIC AN ATTRIBUTION CLAIM ACTUALLY RESTS ON. A reader of a per-patient
# explanation does not read 19 nats values; they read an ORDER -- which signal
# mattered most for this patient, and which three. So the question is whether
# two specifications agree about that order.
#
# `top1` is the share of patients for whom both arms name the SAME signal as the
# largest contributor by |L|. `top3` is the share whose top-three SETS agree
# exactly. `rho` is the per-patient Spearman correlation of the |L| ranking over
# all 19 signals.
#
# `rho` is CONSERVATIVE: unmeasured signals sit at exactly 0 in both arms, so
# they are tied and agreeing, which inflates it. `top1` and `top3` do not have
# that problem because a zero cell can never be a top contributor for a patient
# with any measured signal. Read those two.
cat("\n=== 3. within-patient attribution ranking ===\n\n")
cat("    top1  share of patients whose LARGEST-contributing signal is the same\n")
cat("    top3  share whose top-three set is identical\n")
cat("    rho   median per-patient Spearman of the |L| ranking (conservative:\n")
cat("          unmeasured signals are tied at 0 in both arms)\n\n")
rank_rows <- function(M) t(apply(abs(M), 1, rank, ties.method = "average"))
top_k_set <- function(M, k) t(apply(abs(M), 1, function(v) sort(order(-v)[seq_len(k)])))

#' Per-row Pearson correlation of two matrices, vectorised.
#'
#' On rank matrices this IS the per-patient Spearman correlation. Written out
#' rather than looped: `cor()` once per patient is 41,250 calls per comparison
#' and seven comparisons, which dominates the whole script; this is three matrix
#' passes. NA where a row has no variance, which happens when a patient has one
#' measured signal or none.
row_cor <- function(X, Y) {
  Xc <- X - rowMeans(X); Yc <- Y - rowMeans(Y)
  den <- sqrt(rowSums(Xc^2) * rowSums(Yc^2))
  ifelse(den > 0, rowSums(Xc * Yc) / den, NA_real_)
}

# Precomputed ONCE PER ARM rather than once per comparison. Seven comparisons
# over eight arms would otherwise rank the same matrix up to three times.
PRE <- lapply(A, function(M) list(
  rank = rank_rows(M),
  top1 = max.col(abs(M), ties.method = "first"),
  top3 = top_k_set(M, 3L)))

cat(sprintf("%-16s %-16s %9s %9s %9s %12s\n", "arm a", "arm b",
            "top1", "top3", "med rho", "n_patients"))
R3 <- list()
for (p in CMP) {
  pa <- PRE[[p[1]]]; pb <- PRE[[p[2]]]
  rho <- row_cor(pa$rank, pb$rank)
  t1 <- mean(pa$top1 == pb$top1)
  t3 <- mean(rowSums(pa$top3 == pb$top3) == 3L)
  R3[[length(R3) + 1L]] <- data.frame(a = p[1], b = p[2],
    top1_agree = round(t1, 5), top3_agree = round(t3, 5),
    rho_median = round(stats::median(rho, na.rm = TRUE), 5),
    rho_p10 = round(unname(stats::quantile(rho, 0.10, na.rm = TRUE)), 5),
    n_patients = nrow(pa$rank), stringsAsFactors = FALSE)
  cat(sprintf("%-16s %-16s %9.4f %9.4f %9.4f %12d\n", p[1], p[2], t1, t3,
              stats::median(rho, na.rm = TRUE), nrow(pa$rank)))
}
save_table(run, do.call(rbind, R3), "attribution_ranking", subdir = "diagnostics")

# --- 4. how the evidence is distributed across a patient's vector -----------
#
# Two different things can happen when a specification changes. The patient's
# TOTAL evidence can change, or the same total can be redistributed across their
# signals. These separate them.
#
# `hhi` is the Herfindahl index of the |L| shares over the 19 signals: 1/19 when
# every signal contributes equally, 1 when one signal carries everything. It is
# the concentration of a patient's explanation.
cat("\n=== 4. evidence budget and its concentration ===\n\n")
cat("    budget  median per-patient sum of |L| over the 19 signals, in nats\n")
cat("    hhi     median Herfindahl concentration of that patient's shares\n")
cat("            (1/19 = 0.053 is perfectly spread, 1 is one signal only)\n\n")
cat(sprintf("%-16s %10s %10s %10s\n", "arm", "budget", "hhi", "max share"))
B4 <- list()
for (a in names(A)) {
  ab <- abs(A[[a]]); tot <- rowSums(ab)
  sh <- ab / pmax(tot, .Machine$double.eps)
  B4[[a]] <- data.frame(arm = a, budget_median = round(stats::median(tot), 4),
    budget_p90 = round(unname(stats::quantile(tot, 0.90)), 4),
    hhi_median = round(stats::median(rowSums(sh^2)), 5),
    max_share_median = round(stats::median(apply(sh, 1, max)), 5),
    stringsAsFactors = FALSE)
  cat(sprintf("%-16s %10.4f %10.5f %10.5f\n", a, B4[[a]]$budget_median,
              B4[[a]]$hhi_median, B4[[a]]$max_share_median))
}
save_table(run, do.call(rbind, B4), "evidence_budget", subdir = "diagnostics")

# --- 5. domain level, the unit the paper reports -----------------------------
#
# `config/domains.csv` is the frozen reporting partition: 19 signals into 11
# disjoint domains. Equal weight within a domain, which is what `D_k` reduces to
# and the same convention `R/09d_sofa.R` uses. An attribution claim in the paper
# is made at THIS level, so the flip and ranking metrics are repeated here.
cat("\n=== 5. domain level (11 domains, equal weight within domain) ===\n\n")
dm <- domains$domain[match(sigs, domains$signal)]
dom_names <- sort(unique(dm))
to_dom <- function(M) {
  D <- matrix(0, nrow(M), length(dom_names), dimnames = list(rownames(M), dom_names))
  for (k in dom_names) {
    j <- which(dm == k)
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE])
  }
  D
}
dmeas <- to_dom(meas_ok * 1) > 0
DOM <- lapply(A, function(M) { D <- to_dom(M); list(
  D = D, rank = rank_rows(D), top1 = max.col(abs(D), ties.method = "first")) })
cat(sprintf("%-16s %-16s %9s %9s %9s %9s\n", "arm a", "arm b",
            "flip_t0", "flip_t10", "top1", "med rho"))
D5 <- list()
for (p in CMP) {
  qa <- DOM[[p[1]]]; qb <- DOM[[p[2]]]
  fl <- sign(qa$D) != sign(qb$D); mx <- pmax(abs(qa$D), abs(qb$D))
  rho <- row_cor(qa$rank, qb$rank)
  t1 <- mean(qa$top1 == qb$top1)
  D5[[length(D5) + 1L]] <- data.frame(a = p[1], b = p[2], n_domains = ncol(qa$D),
    flip_t0 = round(mean(fl[dmeas]), 5),
    flip_t10 = round(mean(fl[dmeas] & mx[dmeas] > 0.10), 5),
    top1_agree = round(t1, 5),
    rho_median = round(stats::median(rho, na.rm = TRUE), 5),
    stringsAsFactors = FALSE)
  cat(sprintf("%-16s %-16s %9.4f %9.4f %9.4f %9.4f\n", p[1], p[2],
              D5[[length(D5)]]$flip_t0, D5[[length(D5)]]$flip_t10, t1,
              D5[[length(D5)]]$rho_median))
}
save_table(run, do.call(rbind, D5), "domain_attribution", subdir = "diagnostics")

# --- 6. does the disagreement live where the coupling lives? ----------------
#
# THE FALSIFIABLE PREDICTION. If specification disagreement is really about
# treatment context, it must concentrate among stays that RECEIVED the paired
# intervention. If instead it is spread evenly, it is estimation noise wearing a
# coupling label. Among unexposed stays every intervention covariate is constant
# (verified in tests/coupling_strata.R), so a `ti` surface there is a function of
# the measurement alone and can still move the fit -- the prediction is a
# DIFFERENCE in flip rate, not a zero.
cat("\n=== 6. flip rate by exposure to the signal's own paired intervention ===\n\n")
iv <- tabs$intervention_features
iv <- iv[iv$stay_id %in% tr, , drop = FALSE]
cat(sprintf("%-18s %10s %10s %10s %10s %9s\n", "signal", "n_expo", "flip_expo",
            "n_unexp", "flip_unexp", "ratio"))
E6 <- list()
for (sg in paired) {
  z <- iv[iv$intervention %in% interventions_of(sg, cfg), , drop = FALSE]
  hit <- tapply(as.integer(z$ever_active),
                factor(z$stay_id, levels = ids_ch), max, default = 0L)
  ex <- as.integer(hit[ids_ch]); ex[is.na(ex)] <- 0L; ex <- ex > 0L
  ok <- meas_ok[, sg]
  f <- (sign(A$cond[, sg]) != sign(A$cond_ti_all[, sg])) &
       pmax(abs(A$cond[, sg]), abs(A$cond_ti_all[, sg])) > 0.10
  fe <- mean(f[ok & ex]); fu <- mean(f[ok & !ex])
  E6[[sg]] <- data.frame(arm_a = "cond", arm_b = "cond_ti_all", tau = 0.10,
    interventions = paste(interventions_of(sg, cfg), collapse = "+"),
    signal = sg, n_exposed = sum(ok & ex),
    flip_exposed = round(fe, 5), n_unexposed = sum(ok & !ex),
    flip_unexposed = round(fu, 5),
    ratio = round(fe / pmax(fu, 1e-9), 3), stringsAsFactors = FALSE)
  cat(sprintf("%-18s %10d %10.4f %10d %10.4f %9.2f\n", sg, sum(ok & ex), fe,
              sum(ok & !ex), fu, fe / pmax(fu, 1e-9)))
}
save_table(run, do.call(rbind, E6), "flip_by_exposure", subdir = "diagnostics")

# --- 7. how much of the interaction effect is trend? ------------------------
cat("\n=== 7. the trend share of the interaction effect ===\n\n")
d_all   <- A$cond_ti_all   - A$cond
d_trend <- A$cond_ti_trend - A$cond
cat(sprintf("%-18s %12s %12s %10s %10s\n", "signal", "med|d| all",
            "med|d| trend", "share", "cor(d,d)"))
T7 <- list()
for (sg in paired) {
  ok <- meas_ok[, sg]
  da <- d_all[ok, sg]; dt <- d_trend[ok, sg]
  ma <- stats::median(abs(da)); mt <- stats::median(abs(dt))
  T7[[sg]] <- data.frame(baseline = "cond", arm_all = "cond_ti_all",
    arm_trend = "cond_ti_trend", signal = sg, med_abs_all = round(ma, 5),
    med_abs_trend = round(mt, 5), trend_share = round(mt / pmax(ma, 1e-12), 4),
    cor_deltas = round(if (stats::sd(dt) > 0) stats::cor(da, dt) else NA_real_, 5),
    has_trend = trend_enabled_for(sg, cfg), stringsAsFactors = FALSE)
  cat(sprintf("%-18s %12.5f %12.5f %10.4f %10.4f\n", sg, ma, mt,
              mt / pmax(ma, 1e-12), T7[[sg]]$cor_deltas))
}
save_table(run, do.call(rbind, T7), "trend_share", subdir = "diagnostics")

cat(sprintf("\n  complete in %.1f minutes.\n", t0()$elapsed_sec / 60))
finalize_run(run, extra = list(k_ti = K_TI, taus = paste(TAUS, collapse = ","),
                               reused = !is.null(reuse),
                               design_key_hash = attr_key_hash(.ladder_key())))
cat(sprintf("\nwritten: %s\n", run$path))
