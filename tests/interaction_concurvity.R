# tests/interaction_concurvity.R ----------------------------------------------
# THE GATE THE REOPENED DECISION NAMED. Reads what the graph already computed
# and fits nothing. Runs in seconds.
#
# --- WHY THIS EXISTS ---------------------------------------------------------
#
# CLAUDE.md's frozen decision "Layer 1 carries no interaction terms at all" gave
# two grounds. The conceptual one -- that `s(trend, by = o_flag)` was built on
# excursion-relative timing -- is answered by construction, because a
# `ti(measurement, intervention intensity)` term refers to no excursion and
# `o_flag` no longer exists. The EMPIRICAL one was never answered:
#
#     "the `by=` version also drove worst-case concurvity to ~1 in every paired
#      dense model."
#
# The counter-argument -- that `ti()` excludes the main effects a `by=` smooth
# re-carries, so it should not reproduce that number -- lived in prose in
# `tests/coupling_interaction.R` and had never been checked against the
# `concurvity()` values the pipeline itself records for every fit.
#
# It is checkable for free, because `R/08_diagnostics.R` extracts
# `concurvity_max`, `concurvity_obs` and `concurvity_est` from EVERY fit. So the
# 21 interaction specs are already measured on the same instrument, on the same
# rows, against the same config threshold as the 43 additive ones. This script
# does the reading.
#
# --- WHAT IT REPORTS, AND WHAT WOULD COUNT AS A FAILURE ----------------------
#
# The comparison is WITHIN A SIGNAL: `full_ti_trend` and `full_ti_all` against
# that same signal's `full`. A cross-signal comparison would be dominated by
# which signals have three paired interventions rather than by what the
# interaction did.
#
#   d_worst   the ti model's `concurvity_max` minus its own `full`'s. This is
#             the number the frozen decision is about. `by=` drove it to ~1;
#             the question is whether `ti()` does.
#   d_obs     the same on `concurvity_obs`. `worst` is a deliberately
#             pessimistic upper bound over the whole span, `observed` is what
#             the fitted functions actually do, and R/08's header says a flagged
#             row is unreadable without both.
#   excess    `observed` minus the SYNTHETIC null's upper bound, joined from a
#             `cnull_*` run. A concurvity of 0.85 is not interpretable on its
#             own -- a matched null with no shared structure already produces a
#             substantial value at this design -- so the calibrated excess is
#             the honest quantity. `worst` is calibrated against the PERMUTATION
#             null and `observed` against the SYNTHETIC one; see
#             `tests/concurvity_null.R`'s header for why they cannot be swapped.
#
# THERE IS NO PASS THRESHOLD HERE ON PURPOSE. `diagnostics.concurvity_max` in
# config is the pipeline's declared flag and every fit is already checked
# against it by `diagnostics_summary()`. What this adds is the DIFFERENCE, which
# no threshold can express: an interaction model at 0.87 where its own additive
# model is at 0.85 is a different fact from one at 0.87 where its additive model
# is at 0.30, and only the first is evidence for keeping the arm.
#
# The verdict is a judgement and is stated as one. If the ti specs reproduce
# ~1 where their additive counterparts do not, that is a reason to report the
# interaction arms as a bounded limitation rather than as scored arms -- and
# the decision to reopen was taken knowing that.
#
# Aggregates only (hard rule 1): one row per (signal, model). No stay is read.
#
#   Rscript tests/interaction_concurvity.R
#   Rscript tests/interaction_concurvity.R --null out/runs/cnull_20260828T144343
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(targets); library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default = NA_character_) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[i + 1L]
}
NULL_D <- .opt("--null", latest_run("cnull", require_complete = FALSE))

cfg <- tar_read(cfg)
THR <- as.numeric(cfg_req(cfg, "diagnostics", "concurvity_max"))

cat("\n=== interaction concurvity: the gate the reopened decision named ===\n\n")
cat(sprintf("  config flag `diagnostics.concurvity_max` : %.2f\n", THR))

# --- the two tables the graph already holds ----------------------------------
DF <- tar_read(diag_final)
DO <- tar_read(diag_oof)
keep <- c("signal", "model", "concurvity_max", "concurvity_obs",
          "concurvity_est", "concurvity_term", "n_terms", "n_smooth")
miss <- setdiff(keep, names(DF))
if (length(miss)) abort_values("diag_final is missing column(s)", miss)

# FINAL FITS ARE THE PRIMARY READING and out-of-fold is the stability check.
# The bundle carries the final fits, so a concurvity claim about what
# transports is a claim about those; the five fold fits say whether the number
# is a property of the design or of one subset.
fin <- DF[, keep, drop = FALSE]

wide <- function(d, col) {
  z <- stats::reshape(d[, c("signal", "model", col)], idvar = "signal",
                      timevar = "model", direction = "wide")
  names(z) <- sub(paste0("^", col, "\\."), "", names(z))
  z
}
W  <- wide(fin, "concurvity_max")
WO <- wide(fin, "concurvity_obs")

paired <- Filter(function(s) length(interventions_of(s, cfg)) > 0L,
                 as.character(unlist(cfg$signals)))
rows <- list()
for (sg in paired) {
  i <- match(sg, W$signal)
  if (is.na(i)) next
  for (md in LAYER1_TI_MODELS) {
    # A ti model that ALIASES `full` was never fitted, so it has no diagnostics
    # row and there is nothing to compare. Reported as such rather than skipped
    # silently: three of the twelve paired signals are in this state and a
    # reader counting rows should see why.
    src <- spec_source(sg, md, cfg)
    if (!is.na(src)) {
      rows[[length(rows) + 1L]] <- data.frame(
        signal = sg, model = md, fitted = FALSE, aliases = src,
        full_worst = W$full[i], ti_worst = NA_real_, d_worst = NA_real_,
        full_obs = WO$full[i], ti_obs = NA_real_, d_obs = NA_real_,
        n_ti = 0L, over_threshold = NA, stringsAsFactors = FALSE)
      next
    }
    if (!md %in% names(W)) next
    rows[[length(rows) + 1L]] <- data.frame(
      signal = sg, model = md, fitted = TRUE, aliases = NA_character_,
      full_worst = W$full[i], ti_worst = W[[md]][i],
      d_worst = round(W[[md]][i] - W$full[i], 4),
      full_obs = WO$full[i], ti_obs = WO[[md]][i],
      d_obs = round(WO[[md]][i] - WO$full[i], 4),
      n_ti = length(interaction_terms(sg, cfg, .ti_scope_of(md))),
      over_threshold = W[[md]][i] > THR, stringsAsFactors = FALSE)
  }
}
CMP <- do.call(rbind, rows)
CMP <- CMP[order(-CMP$d_worst, na.last = TRUE), , drop = FALSE]

cat("\n=== per signal: each ti model against its OWN `full` (final fits) ===\n")
cat("  `d_worst` is the number the frozen decision is about. The `by=` version\n")
cat("  it rejected drove worst-case concurvity to about 1 in every paired dense\n")
cat("  model; a `ti()` excludes the main effects a `by=` re-carries, and this\n")
cat("  table is whether that argument survives measurement.\n\n")
print(CMP[CMP$fitted, c("signal", "model", "n_ti", "full_worst", "ti_worst",
                        "d_worst", "full_obs", "ti_obs", "d_obs",
                        "over_threshold")], row.names = FALSE)
if (any(!CMP$fitted)) {
  cat("\n  not fitted (empty cross-term set, aliases another model):\n\n")
  print(CMP[!CMP$fitted, c("signal", "model", "aliases")], row.names = FALSE)
}

f <- CMP[CMP$fitted, , drop = FALSE]
cat(sprintf("\n  fitted interaction specs      : %d\n", nrow(f)))
cat(sprintf("  worst-case concurvity, range  : %.4f to %.4f\n",
            min(f$ti_worst), max(f$ti_worst)))
cat(sprintf("  their `full` counterparts     : %.4f to %.4f\n",
            min(f$full_worst), max(f$full_worst)))
cat(sprintf("  d_worst  median %.4f   max %.4f\n",
            stats::median(f$d_worst), max(f$d_worst)))
cat(sprintf("  d_obs    median %.4f   max %.4f\n",
            stats::median(f$d_obs), max(f$d_obs)))
cat(sprintf("  over the config flag of %.2f  : %d of %d ti specs, %d of %d `full` specs\n",
            THR, sum(f$ti_worst > THR), nrow(f),
            sum(W$full > THR, na.rm = TRUE), sum(!is.na(W$full))))
cat(sprintf("  at or above 0.99 (the `by=` regime): %d of %d\n",
            sum(f$ti_worst >= 0.99), nrow(f)))

# --- against the calibrated null ---------------------------------------------
NUL <- NULL
if (!is.null(NULL_D) && !is.na(NULL_D) && dir.exists(NULL_D)) {
  np <- file.path(NULL_D, "diagnostics", "concurvity_null.rds")
  if (file.exists(np)) NUL <- readRDS(np)
}
if (is.null(NUL)) {
  cat("\n  NO CALIBRATED NULL JOINED. `tests/concurvity_null.R` has not been run,\n")
  cat("  or --null names no such run, so the numbers above are raw. A raw\n")
  cat("  concurvity is not interpretable on its own: a matched null with no\n")
  cat("  shared structure already produces a substantial value at this design.\n")
} else {
  cat(sprintf("\n=== against the calibrated null in %s ===\n", basename(NULL_D)))
  cat("  THE NULL PREDATES THE INTERACTION SPECS, so it has no row for a ti\n")
  cat("  model. What it can calibrate is each ti model's own `full`, which is\n")
  cat("  the baseline `d_worst` is measured from -- so the reading is: is the\n")
  cat("  additive model's concurvity already excess over its null, and does the\n")
  cat("  interaction add to that or sit inside it?\n\n")
  nz <- NUL[NUL$measure == "observed" & NUL$model == "full",
            c("signal", "observed", "null_hi", "excess_used", "null_used"),
            drop = FALSE]
  J <- merge(f[, c("signal", "model", "full_obs", "ti_obs", "d_obs")],
             nz, by = "signal", all.x = TRUE)
  J$ti_excess_vs_full_null <- round(J$ti_obs - J$null_hi, 4)
  J$full_excess_vs_null    <- round(J$full_obs - J$null_hi, 4)
  J <- J[order(-J$ti_excess_vs_full_null), , drop = FALSE]
  print(J[, c("signal", "model", "null_used", "null_hi", "full_obs",
              "full_excess_vs_null", "ti_obs", "ti_excess_vs_full_null",
              "d_obs")], row.names = FALSE)
  cat("\n  The ti model's null is NOT this null -- a tensor block changes the\n")
  cat("  model matrix and therefore what a matched null is -- so\n")
  cat("  `ti_excess_vs_full_null` is an approximation and is labelled as one.\n")
  cat("  Re-run tests/concurvity_null.R to calibrate the interaction specs\n")
  cat("  properly if this table is close to a decision.\n")
}

# --- out-of-fold stability ----------------------------------------------------
oof <- DO[DO$model %in% LAYER1_TI_MODELS,
          c("signal", "model", "fold", "concurvity_max"), drop = FALSE]
if (nrow(oof)) {
  S <- do.call(rbind, lapply(split(oof, paste(oof$signal, oof$model)), function(z)
    data.frame(signal = z$signal[1], model = z$model[1], n_folds = nrow(z),
               min = round(min(z$concurvity_max), 4),
               median = round(stats::median(z$concurvity_max), 4),
               max = round(max(z$concurvity_max), 4),
               range = round(max(z$concurvity_max) - min(z$concurvity_max), 4),
               stringsAsFactors = FALSE)))
  S <- S[order(-S$max), , drop = FALSE]
  cat("\n=== out of fold: is the number a property of the design or of a subset? ===\n\n")
  print(utils::head(S, 12), row.names = FALSE)
  cat(sprintf("\n  fold-to-fold range: median %.4f, max %.4f over %d spec(s)\n",
              stats::median(S$range), max(S$range), nrow(S)))
}

run <- new_run("ticonc", cfg, note = sprintf(
  "interaction concurvity: %d fitted ti specs against their own `full`", nrow(f)))
save_table(run, CMP, "interaction_concurvity")
if (exists("J")) save_table(run, J, "interaction_concurvity_vs_null", subdir = "diagnostics")
if (nrow(oof)) save_table(run, S, "interaction_concurvity_oof", subdir = "diagnostics")
finalize_run(run, extra = list(threshold = THR,
                               null_run = if (is.null(NUL)) "" else basename(NULL_D)))

cat("\n  NO VERDICT IS EMITTED. `d_worst` near zero means the `ti()` argument\n")
cat("  survived measurement; `ti_worst` near 1 where `full_worst` is not means\n")
cat("  it did not, and the interaction arms should then be reported as a\n")
cat("  bounded limitation rather than scored. That is a judgement about the\n")
cat("  design and no threshold decides it.\n\n")
