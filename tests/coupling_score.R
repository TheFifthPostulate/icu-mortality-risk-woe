# tests/coupling_score.R -----------------------------------------------------
# DOES THE INTERACTION CHANGE THE PREDICTION, OR ONLY THE FIT?
#
# tests/coupling_interaction.R measures what a `ti(measurement, intervention)`
# term buys INSIDE one layer-1 model, in deviance explained. That is not the
# same question as whether it changes the SCORE, because an L is a difference of
# log-odds and nineteen of them are summed before anything is ranked. A term can
# improve a model's fit and move no patient's position in the ordering.
#
# This script closes that gap. It builds the L matrix twice -- once from the
# additive `full` fits the pipeline actually uses, once from `full` plus every
# measurement x intervention cross term -- sums each into a score, and compares
# them on identical rows through R/09's `score_report()` machinery.
#
# IN-SAMPLE, AND SAYS SO EVERYWHERE. Both arms are fitted on full train and
# predicted on full train, so every absolute number here is optimistic. The
# DIFFERENCE between two arms measured the same way is the quantity, and it is
# an UPPER BOUND on what an out-of-fold version would find: the interaction has
# many more parameters than the additive model, so in-sample it is flattered by
# more. If the difference is small here it cannot be large out-of-fold.
#
# WHAT IS FITTED AND WHAT IS REUSED. Only the 12 `full+` models are fitted. The
# additive `full`, `meas` and `intv` predictions come from the live bundle's
# final GAMs, which are the pipeline's own objects -- so the additive arm is not
# a reimplementation that might differ, it is the thing itself.
#
# THE UNPAIRED SEVEN ARE UNTOUCHED. They have no intervention, so no cross term
# exists and their L columns are identical in both arms by construction. Any
# difference in the summed score comes from the 12 paired signals alone, which
# is the design's own control.
#
# Aggregates only, never a row (hard rule 1). Scores and L columns are
# row-level and are never printed.
#
#   Rscript tests/coupling_score.R
#   Rscript tests/coupling_score.R out/runs/internal_20260905T121017
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml); library(qs2)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args     <- commandArgs(trailingOnly = TRUE)
bundle_d <- if (length(args) >= 1L && nzchar(args[1])) args[1] else
              latest_run("internal", require_complete = FALSE)
K_TI     <- 5L

bundle <- qs2::qs_read(file.path(bundle_d, "bundle.qs2"))
cfg    <- load_config("config/config.yml")
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
tr     <- folds$stay_id[folds$split == "train"]
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
y      <- tabs$cohort$mortality[match(tr, tabs$cohort$stay_id)]

run <- new_run("coupscore", cfg, note = sprintf(
  "in-sample score effect of measurement x intervention interactions, bundle %s",
  basename(bundle_d)))

paired <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)

# --- helpers (mirrors of tests/coupling_interaction.R) -----------------------

.k_ti <- function(sg, v, cfg, cap = K_TI) max(3L, min(cap, smooth_k_of(sg, v, cfg)))

.blocks <- function(b, sg, cfg) {
  sv <- unique(vapply(b$smooth, function(s) s$term[1], character(1)))
  mv <- .term_vars(measurement_terms(sg, cfg))
  list(meas = intersect(sv, mv), intv = setdiff(sv, mv))
}

#' L for one fitted model on one set of rows: logit(p_hat) - logit(p_bar).
#'
#' Identical arithmetic to `fit_one()`, and `p_bar` comes from the SAME frozen
#' prior object both arms use, so `L_full+ - L_intv` is the conditional term for
#' exactly the reason `L_full - L_intv` is.
.L <- function(b, nd, p_bar) {
  as.numeric(stats::predict(b, newdata = nd, type = "link", discrete = FALSE)) -
    logit(p_bar)
}

# --- build both L matrices --------------------------------------------------

sigs <- unlist(cfg$signals)
mk <- function() matrix(0, nrow = length(tr), ncol = length(sigs),
                        dimnames = list(as.character(tr), sigs))
L_full <- mk(); L_plus <- mk(); L_intv <- mk(); L_meas <- mk()
meas_ok <- measured_matrix(tabs, cfg, tr)

cat("\n=== fitting the interaction arm (in-sample) ===\n\n")
cat(sprintf("%-18s %6s %8s %10s %10s %10s\n",
            "signal", "n_ti", "rows", "sd L_full", "sd L_plus", "cor"))
t0 <- start_timer()

for (sg in sigs) {
  pri     <- priors_for(priors, sg, "final")
  is_pair <- sg %in% paired
  ids     <- tr[meas_ok[, sg]]

  d_full <- signal_frame(sg, "full", tabs, cfg, pri, stay_ids = ids)
  d_meas <- signal_frame(sg, "meas", tabs, cfg, pri, stay_ids = ids)
  j <- match(as.character(d_full$stay_id), as.character(tr))

  b_full <- bundle$models[[paste0(sg, "/", if (is_pair) "full" else "meas")]]
  b_meas <- bundle$models[[paste0(sg, "/meas")]]
  L_full[j, sg] <- .L(b_full, d_full, pri$p_bar)
  L_meas[j, sg] <- .L(b_meas, d_meas, pri$p_bar)

  if (!is_pair) {
    # No intervention: L_intv is 0 by assignment on measured stays (spec 5.5),
    # and no cross term exists, so the interaction arm equals the additive arm.
    L_plus[j, sg] <- L_full[j, sg]
    next
  }

  b_intv <- bundle$models[[paste0(sg, "/intv")]]
  L_intv[j, sg] <- .L(b_intv, signal_frame(sg, "intv", tabs, cfg, pri, stay_ids = ids),
                      pri$p_bar)

  f0 <- build_formula(sg, "full", cfg)
  bl <- .blocks(b_full, sg, cfg)
  ti <- unlist(lapply(bl$meas, function(mv) vapply(bl$intv, function(iv)
    sprintf("ti(%s, %s, bs = c(\"ts\", \"ts\"), k = c(%d, %d))",
            mv, iv, .k_ti(sg, mv, cfg), .k_ti(sg, iv, cfg)), character(1))))
  f1 <- stats::update(f0, stats::as.formula(paste(". ~ . +", paste(ti, collapse = " + "))))
  b1 <- try(.bam_fit(f1, d_full, cfg), silent = TRUE)
  if (inherits(b1, "try-error")) {
    cat(sprintf("%-18s  FIT FAILED, falling back to the additive column\n", sg))
    L_plus[j, sg] <- L_full[j, sg]
    next
  }
  L_plus[j, sg] <- .L(b1, d_full, pri$p_bar)

  cat(sprintf("%-18s %6d %8d %10.4f %10.4f %10.5f\n", sg, length(ti), length(ids),
              stats::sd(L_full[j, sg]), stats::sd(L_plus[j, sg]),
              stats::cor(L_full[j, sg], L_plus[j, sg])))
}

# --- the five arms ----------------------------------------------------------
scores <- list(
  llr_meas      = rowSums(L_meas),
  llr_sum       = rowSums(L_full),
  llr_cond      = rowSums(L_full - L_intv),
  llr_sum_ti    = rowSums(L_plus),
  llr_cond_ti   = rowSums(L_plus - L_intv))

cat("\n=== in-sample discrimination, all five arms on identical rows ===\n\n")
cat("    EVERY NUMBER IS IN-SAMPLE and therefore optimistic. Read the\n")
cat("    DIFFERENCES between arms, and read them as UPPER BOUNDS.\n\n")
cat(sprintf("%-14s %9s %9s %9s %9s\n", "arm", "AUROC", "AUPRC", "cal_slope", "top1%"))
S <- list()
tp1 <- function(x) { k <- ceiling(0.01 * length(x)); mean(y[order(-x)][seq_len(k)]) }
for (a in names(scores)) {
  cal <- llr_calibration(scores[[a]], y, mean(y))
  S[[a]] <- data.frame(arm = a, auroc = round(.auroc(scores[[a]], y), 5),
                       auprc = round(.auprc(scores[[a]], y), 5),
                       cal_slope = cal$slope, top1 = round(tp1(scores[[a]]), 5),
                       in_sample = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("%-14s %9.5f %9.5f %9.4f %9.4f\n", a, S[[a]]$auroc, S[[a]]$auprc,
              S[[a]]$cal_slope, S[[a]]$top1))
}
save_table(run, do.call(rbind, S), "coupling_score_arms", subdir = "diagnostics")

cat("\n=== paired tests, the contrasts that answer the question ===\n\n")
prs <- list(c("llr_sum_ti", "llr_sum"), c("llr_cond_ti", "llr_cond"),
            c("llr_cond_ti", "llr_meas"), c("llr_sum", "llr_meas"))
C <- list()
for (p in prs) {
  dt <- delong_test(scores[[p[1]]], scores[[p[2]]], y)
  C[[length(C) + 1L]] <- data.frame(a = p[1], b = p[2], auroc_a = dt$auroc_1,
                                    auroc_b = dt$auroc_2, d_auroc = dt$delta,
                                    z = dt$z, p_value = dt$p_value,
                                    stringsAsFactors = FALSE)
  cat(sprintf("  %-13s vs %-13s  d_AUROC = %+.5f   z = %+7.2f   p = %.3g\n",
              p[1], p[2], dt$delta, dt$z, dt$p_value))
}
save_table(run, do.call(rbind, C), "coupling_score_contrasts", subdir = "diagnostics")

# --- does the interaction change the evidence GEOMETRY? ---------------------
cat("\n=== eigenspectrum of the conditional L matrix, with and without ===\n\n")
E <- list()
for (nm in c("cond", "cond_ti")) {
  M <- if (nm == "cond") L_full - L_intv else L_plus - L_intv
  sp <- spectrum_summary(l_correlation(M), label = nm)
  E[[nm]] <- sp
  cat(sprintf("  %-9s PC1 = %.5f   PC1+2 = %.5f   participation ratio = %.3f\n",
              nm, sp$pc1, sp$pc12, sp$pr))
}
save_table(run, do.call(rbind, E), "coupling_score_spectrum", subdir = "diagnostics")

cat(sprintf("\n  fitted 12 interaction models in %.1f minutes.\n", t0()$elapsed_sec / 60))
cat("  IF d_AUROC for llr_sum_ti vs llr_sum is small IN SAMPLE, it cannot be\n")
cat("  large out of fold, and the interaction is a fit-quality finding rather\n")
cat("  than a prediction finding.\n")

finalize_run(run, extra = list(bundle = basename(bundle_d), k_ti = K_TI,
                               n_paired = length(paired)))
cat(sprintf("\nwritten: %s\n", run$path))
