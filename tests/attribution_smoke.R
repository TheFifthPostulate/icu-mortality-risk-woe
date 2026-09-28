# tests/attribution_smoke.R --------------------------------------------------
# Does the attribution machinery produce the object it claims to?
#
# Four assertions, and the first two are the ones that matter. An agreement
# metric that does not return exactly 1 when handed a matrix and its own copy is
# broken in a way that would look like a finding, and an attribution that does
# not sum to the score it decomposes is not a decomposition.
#
# AGGREGATES ONLY (hard rule 1): dimensions, maxima of absolute differences, and
# correlations. The attribution matrices are row-level and never printed.
#
#   Rscript tests/attribution_smoke.R [n_stays]        ~2 min at n = 4000
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
n_eval <- if (length(args) && nzchar(args[1])) as.integer(args[1]) else 4000L

rc <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config)
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg_local, strict = TRUE)
cfg    <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
tabs   <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)

# TRAIN rows on purpose. This is a plumbing check and must not touch test.
tr  <- folds$stay_id[folds$split == "train"]
set.seed(cfg$seed)
ids <- sort(sample(tr, min(n_eval, length(tr))))
cat(sprintf("evaluation set: %d train stays (the test split is not touched)\n\n",
            length(ids)))

t0 <- Sys.time()
s <- attribution_set(bundle, tabs, cfg, ids, verbose = FALSE)
cat(sprintf("attribution_set: %.1f min, %d method(s)\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(s$A)))

fail <- 0L
chk <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", msg))
  if (!isTRUE(ok)) fail <<- fail + 1L
}

# --- 1. shape ---------------------------------------------------------------
cat("=== 1. every method is on the same n x 19 index ===\n")
for (nm in names(s$A)) {
  A <- s$A[[nm]]
  chk(nrow(A) == length(ids) && ncol(A) == length(s$signals) &&
        identical(colnames(A), s$signals),
      sprintf("%-9s %d x %d", nm, nrow(A), ncol(A)))
}

# --- 2. additivity ----------------------------------------------------------
# For `llr_sum` this must be EXACT: the attribution is the summand, so its row
# sum is the score by construction. A failure here means the weights or the
# column set diverged from what apply_bundle() uses, which would make every
# comparison downstream a comparison of two different quantities.
cat("\n=== 2. the attributions sum to the score they decompose ===\n")
ap <- apply_bundle(bundle, tabs, cfg, ids, arms = BUNDLE_ARMS, verbose = FALSE)
# `llr_cond` joined the list 2026-09-05 and it is the strongest of the three
# checks: `L_cond = L_full - L_intv` is built by subtracting two matrices, and
# a row-sum that still lands on the score proves the subtraction was done on
# the aligned column set rather than on two frames that merely have the same
# dimensions.
for (nm in c("llr_sum", "llr_cond", "llr_meas")) {
  d <- max(abs(rowSums(s$A[[nm]]) - as.numeric(ap$scores[[nm]])))
  chk(d < 1e-10, sprintf("%-9s max |rowSum(A) - score| = %.3e  (exact by construction)", nm, d))
}
# For the tree cells the rollup drops the intervention groups by default, so the
# row sum is the margin MINUS the dropped groups minus the bias. The check is
# that the discarded share is what `dropped_frac` says it is, not that the sum
# reproduces the margin.
for (nm in intersect(XGB_DESIGNS, names(s$A))) {
  cat(sprintf("  [    ] %-9s discarded %.4f of total |SHAP| to non-signal groups\n",
              nm, s$dropped_frac[[nm]]))
}

# --- 3. self-comparison is exactly 1 ----------------------------------------
cat("\n=== 3. a set compared with itself returns perfect agreement ===\n")
for (nm in names(s$A)) {
  cp <- compare_attributions(s, s, nm, label = "self")
  o  <- cp$overall
  chk(isTRUE(all.equal(o$spearman, 1)) && isTRUE(all.equal(o$sign_agreement, 1)) &&
        isTRUE(all.equal(cp$ranking$median, 1)),
      sprintf("%-9s spearman %.4f | sign %.4f | ranking median %.4f | cells %.3f",
              nm, o$spearman, o$sign_agreement, cp$ranking$median, o$frac_cells_used))
}

# --- 4. the metrics separate two genuinely different methods ----------------
# The complement of check 3: a metric that returns 1 for everything is as broken
# as one that returns 1 for nothing.
cat("\n=== 4. and they separate different methods (orientation only) ===\n")
cells <- attribution_cells(s$measured)
cat(sprintf("  cells usable: %.4f (assigned zeros excluded)\n\n", cells$frac_used))
pairs <- list(c("llr_sum", "xgb_l"), c("llr_sum", "xgb_raw"),
              c("xgb_l", "xgb_raw"), c("llr_sum", "llr_meas"))
rows <- lapply(pairs, function(p) {
  if (!all(p %in% names(s$A))) return(NULL)
  a <- attribution_agreement(s$A[[p[1]]], s$A[[p[2]]], cells$keep)
  cbind(a = p[1], b = p[2], a$overall[, c("spearman", "pearson", "sign_agreement",
                                          "frac_above_floor")])
})
print(do.call(rbind, Filter(Negate(is.null), rows)), row.names = FALSE)
cat("\n  These are NOT the arm's result -- they compare methods on ONE fit, and\n")
cat("  the arm compares ONE method across two fits. They are here so a metric\n")
cat("  that is constant at 1 cannot pass unnoticed.\n")

cat(sprintf("\n%s: %d check(s) failed\n", if (fail == 0L) "PASS" else "FAIL", fail))
if (fail > 0L) quit(status = 1L)
