# tests/coupling_oof.R -------------------------------------------------------
# THE OUT-OF-FOLD INTERACTION ANALYSIS. The number that can be quoted.
#
# tests/coupling_score.R measured the same thing in-sample and every figure it
# produced is an upper bound, because the interaction arm carries far more
# parameters and is flattered more by being scored on its own training rows.
# This script removes that: `full + ti` is fitted on four folds and predicted on
# the fifth, exactly as `fit_one()` does for the additive arm, so both arms are
# out-of-fold and the difference between them is honest.
#
# ALL CROSS PAIRS, NO PRE-SELECTION. The interaction set is the same 98 pairs
# tests/coupling_interaction.R enumerated -- every measurement smooth covariate
# crossed with every intervention smooth covariate, per paired signal. Choosing
# a subset by in-sample gain and then measuring it out-of-fold would inflate the
# result, which is the one mistake this run exists to avoid. Scope AFTER seeing
# an unbiased number, not before.
#
# ONLY `full + ti` IS FITTED. `L_intv` is untouched by construction: `p(I | Y)`
# depends on the treatment record alone, so no measurement x intervention term
# can enter the `intv` model, and its out-of-fold L values are read from the
# targets cache. Same for the additive `L_full` and `L_meas`. So the additive
# arm here IS the pipeline's own output, not a reimplementation of it, and the
# 7 unpaired signals are bitwise identical in both arms.
#
# 12 signals x 5 folds = 60 fits.
#
# WHY THE SUBTRACTION IS STILL THE CONDITIONAL LLR. `p(M, I | Y) =
# p(M | I, Y) p(I | Y)` is an identity and does not depend on the model's
# functional form, so `L_full_ti - L_intv` estimates
# `log[ p(M|I,Y=1) / p(M|I,Y=0) ]` for the same reason `L_full - L_intv` does.
# What the tensor changes is the form available to that estimate: additively the
# M-dependence is one function shared across every treatment context, and with
# the interaction it can differ by context. Every M x I term belongs to the
# conditional-measurement factor BY THE PROBABILITY STRUCTURE rather than by
# declaration, which is exactly what `o_flag` could not claim -- it was built
# from the measurement stream under an intervention name and would have
# contaminated `L_intv` itself. A `ti` term never enters `intv`.
#
# THE PER-CELL OBJECT. For each stay and each paired signal this also records
#
#   Delta = L_cond_ti - L_cond = L_full_ti - L_full
#
# -- `L_intv` cancels exactly -- which is how many nats the interaction moves
# THIS patient's evidence from THIS signal. That is the unit the attributability
# arm reports in (docs/v2_attributability_plan_20260902.md PART TWO), so the
# matrix is saved for it rather than recomputed later.
#
# Aggregates only on print (hard rule 1). The cell matrix is row-level, is saved
# to the run directory the way `l_oof.rds` already is, and is never printed or
# summarised per stay.
#
#   Rscript tests/coupling_oof.R
#   Rscript tests/coupling_oof.R mbp,spo2      just those signals
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
only_sg <- if (length(args) >= 1L && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
K_TI    <- 5L
BIG     <- 0.10   # "a materially displaced cell", in nats. Declared before the run.

cfg    <- tar_read(cfg)
tabs   <- tar_read(tabs)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
y      <- as.integer(tar_read(y_train))
priors <- tar_read(priors)
Lz     <- tar_read(l_mats_zero)          # additive, out-of-fold, the pipeline's own

paired  <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)
signals <- if (is.null(only_sg)) paired else intersect(only_sg, paired)

run <- new_run("coupoof", cfg, note = sprintf(
  "out-of-fold measurement x intervention interaction, all cross pairs, %d signal(s)",
  length(signals)))

# --- helpers (mirrors of tests/coupling_interaction.R) -----------------------

.k_ti <- function(sg, v, cfg, cap = K_TI) max(3L, min(cap, smooth_k_of(sg, v, cfg)))

#' Every measurement x intervention cross term for a signal, as `ti()` strings.
#'
#' The block split is read off a FITTED object so it cannot drift from what was
#' estimated; membership comes from `measurement_terms()`, the same function the
#' formula builder uses. `present_at_admission` is parametric and has no smooth,
#' so it is not crossed.
.ti_terms <- function(b, sg, cfg) {
  sv <- unique(vapply(b$smooth, function(s) s$term[1], character(1)))
  mv <- intersect(sv, .term_vars(measurement_terms(sg, cfg)))
  iv <- setdiff(sv, mv)
  unlist(lapply(mv, function(m) vapply(iv, function(i)
    sprintf("ti(%s, %s, bs = c(\"ts\", \"ts\"), k = c(%d, %d))",
            m, i, .k_ti(sg, m, cfg), .k_ti(sg, i, cfg)), character(1))))
}

# --- the out-of-fold interaction fits ---------------------------------------

sigs   <- unlist(cfg$signals)
L_plus <- Lz$full                        # unpaired columns carry over unchanged
n_fold <- cfg$n_folds %||% 5L
ids_ch <- as.character(tr)

cat("\n=== out-of-fold fits: full + every cross ti(), 4 folds -> the 5th ===\n\n")
cat(sprintf("%-18s %6s %8s %10s %10s %9s\n",
            "signal", "n_ti", "n_pred", "sd L_full", "sd L_plus", "cor"))
t0 <- start_timer()
FAIL <- character(0)

for (sg in signals) {
  f0 <- build_formula(sg, "full", cfg)
  col <- rep(NA_real_, length(tr))

  for (k in seq_len(n_fold)) {
    pri <- priors_for(priors, sg, "oof", fold = k)
    ji  <- job_ids("oof", k, folds)
    d_fit <- signal_frame(sg, "full", tabs, cfg, pri, stay_ids = ji$fit_ids)
    d_prd <- signal_frame(sg, "full", tabs, cfg, pri,
                          stay_ids = ji$predict_ids, stage = "predict")

    # The term list is built from a fit on THIS fold's rows, so a covariate that
    # is degenerate in one fold cannot silently import another fold's basis.
    b0 <- try(.bam_fit(f0, d_fit, cfg), silent = TRUE)
    if (inherits(b0, "try-error")) { FAIL <- c(FAIL, sprintf("%s/f%d base", sg, k)); next }
    ti <- .ti_terms(b0, sg, cfg)
    f1 <- stats::update(f0, stats::as.formula(paste(". ~ . +", paste(ti, collapse = " + "))))
    b1 <- try(.bam_fit(f1, d_fit, cfg), silent = TRUE)
    if (inherits(b1, "try-error")) { FAIL <- c(FAIL, sprintf("%s/f%d ti", sg, k)); next }

    eta <- as.numeric(stats::predict(b1, newdata = d_prd, type = "link", discrete = FALSE))
    if (anyNA(eta)) stop(sprintf("coupling_oof [%s fold %d]: NA prediction(s)", sg, k),
                         call. = FALSE)
    j <- match(as.character(d_prd$stay_id), ids_ch)
    col[j] <- eta - logit(pri$p_bar)
  }

  # A stay measured for this signal but never predicted is a fold-scoping
  # failure, not a fill case -- the same check `l_matrix()` makes.
  meas <- Lz$full[, sg] != 0 | Lz$meas[, sg] != 0
  gap  <- is.na(col) & meas
  if (any(gap)) {
    cat(sprintf("%-18s  WARNING: %d measured stay(s) have no interaction L; ",
                sg, sum(gap)))
    cat("additive column kept for those\n")
  }
  col[is.na(col)] <- Lz$full[is.na(col), sg]
  L_plus[, sg] <- col

  ok <- meas
  cat(sprintf("%-18s %6d %8d %10.4f %10.4f %9.5f\n", sg, length(ti), sum(ok),
              stats::sd(Lz$full[ok, sg]), stats::sd(L_plus[ok, sg]),
              stats::cor(Lz$full[ok, sg], L_plus[ok, sg])))
}
if (length(FAIL)) cat(sprintf("\n  %d fold fit(s) failed: %s\n",
                              length(FAIL), paste(FAIL, collapse = ", ")))

# --- A. the five arms, all out-of-fold --------------------------------------
scores <- list(
  llr_meas    = rowSums(Lz$meas),
  llr_sum     = rowSums(Lz$full),
  llr_cond    = rowSums(Lz$full - Lz$intv),
  llr_sum_ti  = rowSums(L_plus),
  llr_cond_ti = rowSums(L_plus - Lz$intv))

cat("\n=== A. out-of-fold discrimination, five arms, identical rows ===\n\n")
cat(sprintf("%-14s %9s %9s %9s %9s\n", "arm", "AUROC", "AUPRC", "cal_slope", "top1%"))
tp1 <- function(x) { k <- ceiling(0.01 * length(x)); mean(y[order(-x)][seq_len(k)]) }
S <- list()
for (a in names(scores)) {
  cal <- llr_calibration(scores[[a]], y, mean(y))
  S[[a]] <- data.frame(arm = a, auroc = round(.auroc(scores[[a]], y), 5),
                       auprc = round(.auprc(scores[[a]], y), 5),
                       cal_slope = cal$slope, top1 = round(tp1(scores[[a]]), 5),
                       out_of_fold = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("%-14s %9.5f %9.5f %9.4f %9.4f\n", a, S[[a]]$auroc, S[[a]]$auprc,
              S[[a]]$cal_slope, S[[a]]$top1))
}
save_table(run, do.call(rbind, S), "oof_score_arms", subdir = "diagnostics")

cat("\n=== the contrasts ===\n\n")
prs <- list(c("llr_sum_ti", "llr_sum"), c("llr_cond_ti", "llr_cond"),
            c("llr_cond_ti", "llr_meas"), c("llr_cond", "llr_meas"),
            c("llr_sum", "llr_meas"))
C <- list()
for (p in prs) {
  dt <- delong_test(scores[[p[1]]], scores[[p[2]]], y)
  pb <- paired_boot_diff(scores[[p[1]]], scores[[p[2]]], y, metric = .auprc, n_boot = 200L)
  C[[length(C) + 1L]] <- data.frame(a = p[1], b = p[2], d_auroc = dt$delta,
                                    z = dt$z, p_value = dt$p_value,
                                    d_auprc = pb$delta %||% NA_real_,
                                    stringsAsFactors = FALSE)
  cat(sprintf("  %-13s vs %-13s  d_AUROC = %+.5f (z=%+6.2f, p=%.3g)   d_AUPRC = %+.5f\n",
              p[1], p[2], dt$delta, dt$z, dt$p_value, pb$delta %||% NA_real_))
}
save_table(run, do.call(rbind, C), "oof_score_contrasts", subdir = "diagnostics")

# --- B. per-cell displacement: the attributability object -------------------
#
# Delta = L_cond_ti - L_cond = L_full_ti - L_full, because L_intv cancels
# exactly. Reported over MEASURED cells only: an unmeasured stay has L = 0 by
# assignment in both arms, so including it would dilute every statistic with
# structural zeros.
cat("\n=== B. per-cell displacement, in nats, over measured cells ===\n\n")
cat("    How many nats the interaction moves ONE patient's evidence from ONE\n")
cat("    signal. The unit the attributability arm reports in.\n")
cat(sprintf("    frac_big is the share of cells moved by more than %.2f nats.\n", BIG))
cat("    sign_flip is the share whose evidence changes DIRECTION -- the signal\n")
cat("    counted for survival in one arm and against it in the other.\n\n")
D <- L_plus - Lz$full
cat(sprintf("%-18s %9s %9s %9s %9s %10s %10s\n",
            "signal", "n_cells", "med|d|", "p90|d|", "max|d|", "frac_big", "sign_flip"))
P <- list()
for (sg in signals) {
  ok <- Lz$full[, sg] != 0 | Lz$meas[, sg] != 0
  d  <- D[ok, sg]; a <- Lz$full[ok, sg]; b <- L_plus[ok, sg]
  P[[sg]] <- data.frame(
    signal = sg, n_cells = sum(ok),
    med_abs = round(stats::median(abs(d)), 5),
    p90_abs = round(unname(stats::quantile(abs(d), 0.90)), 5),
    max_abs = round(max(abs(d)), 5),
    frac_big = round(mean(abs(d) > BIG), 5),
    sign_flip = round(mean(sign(a) != sign(b)), 5),
    cor_cols = round(stats::cor(a, b), 5),
    sd_additive = round(stats::sd(a), 5), sd_ti = round(stats::sd(b), 5),
    stringsAsFactors = FALSE)
  cat(sprintf("%-18s %9d %9.4f %9.4f %9.4f %10.4f %10.4f\n", sg, P[[sg]]$n_cells,
              P[[sg]]$med_abs, P[[sg]]$p90_abs, P[[sg]]$max_abs,
              P[[sg]]$frac_big, P[[sg]]$sign_flip))
}
save_table(run, do.call(rbind, P), "oof_cell_displacement", subdir = "diagnostics")
# Row-level, so it is SAVED and never printed -- the same treatment `l_oof.rds`
# gets. This is the object the attributability arm consumes.
save_object(run, list(L_plus = L_plus, delta = D, stay_id = tr), "l_oof_ti")

# --- C. per-signal discrimination and the evidence geometry -----------------
cat("\n=== C. per-signal L column, additive vs interaction ===\n\n")
sa <- signal_auroc(Lz$full, y); sb <- signal_auroc(L_plus, y)
m <- match(sa$signal, sb$signal)
G <- data.frame(signal = sa$signal, auroc_additive = sa$auroc, auroc_ti = sb$auroc[m],
                d_auroc = round(sb$auroc[m] - sa$auroc, 5),
                lift_additive = sa$auprc_lift, lift_ti = sb$auprc_lift[m],
                stringsAsFactors = FALSE)
G <- G[order(-G$d_auroc), ]
save_table(run, G, "oof_signal_auroc", subdir = "diagnostics")
cat(sprintf("%-18s %10s %10s %10s %10s %10s\n", "signal", "auroc_add", "auroc_ti",
            "d_auroc", "lift_add", "lift_ti"))
for (i in seq_len(nrow(G))) {
  if (!G$signal[i] %in% signals) next
  cat(sprintf("%-18s %10.5f %10.5f %+10.5f %10.3f %10.3f\n", G$signal[i],
              G$auroc_additive[i], G$auroc_ti[i], G$d_auroc[i],
              G$lift_additive[i], G$lift_ti[i]))
}

cat("\n=== D. eigenspectrum of the conditional L matrix ===\n\n")
E <- list()
for (nm in c("cond", "cond_ti")) {
  M <- if (nm == "cond") Lz$full - Lz$intv else L_plus - Lz$intv
  sp <- spectrum_summary(l_correlation(M), label = nm); E[[nm]] <- sp
  cat(sprintf("  %-9s PC1 = %.5f   PC1+2 = %.5f   pr = %.3f   mean|r| = %.4f\n",
              nm, sp$pc1, sp$pc12, sp$pr, sp$mean_abs_r))
}
save_table(run, do.call(rbind, E), "oof_spectrum", subdir = "diagnostics")

cat(sprintf("\n  %d out-of-fold fits in %.1f minutes.\n",
            2 * length(signals) * n_fold, t0()$elapsed_sec / 60))
cat("  THESE ARE THE QUOTABLE NUMBERS. Compare them against\n")
cat("  out/runs/coupscore_* , which measured the same contrasts in sample and\n")
cat("  is an upper bound on every one of them.\n")

finalize_run(run, extra = list(k_ti = K_TI, big_nats = BIG,
                               n_signals = length(signals), n_failed = length(FAIL)))
cat(sprintf("\nwritten: %s\n", run$path))
