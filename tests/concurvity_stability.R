# tests/concurvity_stability.R -----------------------------------------------
# Does concurvity actually destabilise anything, and if so what?
#
# Concurvity is a statement about the DESIGN MATRIX: it says the space spanned
# by one smooth's basis overlaps the space spanned by the others. What it does
# not say is which quantities that overlap makes unreliable. This measures it.
#
# Refit the same model on bootstrap resamples of the training rows and track,
# at a FIXED set of evaluation points:
#
#   per-term    the contribution of each individual smooth
#   block       the summed contribution of a group of concurve smooths
#   eta         the total linear predictor, which is what L is built from
#
# If concurvity is doing what the theory says, the per-term spread should be
# large while the block and eta spread stay small. That is the difference
# between "this model is unreliable" and "this model's DECOMPOSITION is
# unreliable", and the whole interpretability claim turns on which one it is.
#
# Evaluation points are a fixed subsample of REAL rows, held constant across
# resamples: a synthetic grid would put the smooths off the data manifold,
# where any two functions can disagree freely and the answer would be rigged.
#
# AGGREGATES ONLY (hard rule 1). Only standard deviations ACROSS resamples,
# averaged over evaluation points, are printed. No row is ever emitted.
#
#   Rscript tests/concurvity_stability.R [n_boot]
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args   <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args)) as.integer(args[1]) else 25L
N_EVAL <- 400L

cfg    <- load_config("config/config.yml")
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
tr     <- folds$stay_id[folds$split == "train"]

# The two ends of the concurvity range in the real run: creatinine/meas carried
# observed 1.00 on s(value_median) ~ s(value_max); mbp/full carried 0.98 on the
# vasopressor pair and is the widest paired model.
cases <- list(
  list(sg = "creatinine", md = "meas"),
  list(sg = "mbp",        md = "full")
)

for (cs in cases) {
  sg <- cs$sg; md <- cs$md
  pri <- priors_for(priors, sg, "final", NA_integer_)
  d  <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = tr)
  f  <- attr(d, "formula")

  set.seed(cfg$seed)
  eval_rows <- sample.int(nrow(d), min(N_EVAL, nrow(d)))
  nd <- d[eval_rows, , drop = FALSE]

  ctr <- matrix(NA_real_, nrow = B, ncol = 0)
  terms_mat <- NULL
  eta_mat   <- matrix(NA_real_, nrow = B, ncol = nrow(nd))

  for (b in seq_len(B)) {
    set.seed(cfg$seed + b)
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    fit <- try(mgcv::bam(f, data = d[idx, , drop = FALSE], family = stats::binomial(),
                         method = "fREML", discrete = TRUE,
                         nthreads = cfg$bam$nthreads %||% 1L,
                         na.action = stats::na.fail), silent = TRUE)
    if (inherits(fit, "try-error")) next
    tm <- stats::predict(fit, newdata = nd, type = "terms", discrete = FALSE)
    if (is.null(terms_mat)) {
      terms_mat <- array(NA_real_, dim = c(B, nrow(nd), ncol(tm)),
                         dimnames = list(NULL, NULL, colnames(tm)))
    }
    terms_mat[b, , ] <- tm
    eta_mat[b, ] <- as.numeric(stats::predict(fit, newdata = nd, type = "link",
                                              discrete = FALSE))
  }

  keep <- !is.na(eta_mat[, 1])
  tm_names <- dimnames(terms_mat)[[3]]

  # SD across resamples at each evaluation point, then averaged over points.
  sd_of <- function(M) mean(apply(M, 2, stats::sd, na.rm = TRUE), na.rm = TRUE)

  per_term <- vapply(seq_along(tm_names),
                     function(j) sd_of(terms_mat[keep, , j, drop = TRUE]),
                     numeric(1))
  names(per_term) <- tm_names

  # Blocks: the concurve groups. `level` is the pair the diagnostics name.
  is_level <- grepl("value_median|value_min|value_max|q05|q95", tm_names)
  is_pi    <- grepl("pi_minus|pi_plus", tm_names)
  is_iv    <- grepl("__", tm_names)

  block_sd <- function(sel) if (!any(sel)) NA_real_ else
    sd_of(apply(terms_mat[keep, , sel, drop = FALSE], c(1, 2), sum))

  cat(sprintf("\n=== %s / %s  (%d rows, %d bootstrap refits, %d eval points) ===\n",
              sg, md, nrow(d), sum(keep), nrow(nd)))
  cat("\n  bootstrap SD of the fitted contribution, averaged over eval points:\n\n")
  for (j in order(-per_term)) {
    cat(sprintf("    %-46s %.4f\n", tm_names[j], per_term[j]))
  }
  cat(sprintf("\n    %-46s %.4f\n", "SUM of level terms (the block)", block_sd(is_level)))
  if (any(is_pi)) cat(sprintf("    %-46s %.4f\n", "SUM of propensity terms", block_sd(is_pi)))
  if (any(is_iv)) cat(sprintf("    %-46s %.4f\n", "SUM of intervention terms", block_sd(is_iv)))
  cat(sprintf("    %-46s %.4f\n", "TOTAL linear predictor (what L is)", sd_of(eta_mat[keep, ])))

  lv <- per_term[is_level]
  cat(sprintf("\n  ratio: worst single level term / level block = %.2fx\n",
              max(lv) / block_sd(is_level)))
  cat(sprintf("  ratio: worst single term overall  / total eta  = %.2fx\n",
              max(per_term) / sd_of(eta_mat[keep, ])))
}

cat("\n  A large ratio means the DECOMPOSITION moves while the SUM does not:\n")
cat("  the model is well identified, its per-term attribution is not.\n")
