# tests/conditional_ablation.R -----------------------------------------------
# What do the two conditional-prior constructs COST and what do they BUY?
#
# Both are switchable in config, and both are BIJECTIONS of the raw pair given
# frozen parameters -- so neither removes information. What they change is which
# additive subspace the GAM is fitted in, and that is not free. This script
# fits every layer-1 spec twice, with the constructs off and on, on identical
# rows with an identical term count, and reports both sides of the trade:
#
#   A  deviance explained     what the reparameterisation costs or buys
#   B  attributability        bootstrap SD of each term's contribution against
#                             the SD of the total linear predictor. This is the
#                             quantity the interpretability claim rests on, and
#                             it is the reason the constructs exist.
#
# READ BOTH TABLES TOGETHER. A alone would suggest the constructs are neutral
# (mean dev_expl moves -0.4%); B alone would suggest they are a clear win. The
# honest summary is that they trade a little fit, very unevenly across signals,
# for a large gain in whether a fitted smooth means anything.
#
# A WARNING ABOUT TABLE A. Do NOT use it to pick which signals get a construct.
# Selecting covariate construction on the deviance change is selection on the
# response, and it is the same move CLAUDE.md forbids for diagnostic thresholds
# ("fixed in config before fitting, never chosen after seeing which models trip
# them"). An outcome-blind criterion is available and lives in
# tests/conditional_review.R section D -- residual decoupling -- which is what
# may legitimately drive a `magnitude_conditional_override` entry.
#
# AGGREGATES ONLY (hard rule 1): deviance explained, concurvity, and standard
# deviations across resamples. No row is ever emitted.
#
#   Rscript tests/conditional_ablation.R            both sections  ~12 min
#   Rscript tests/conditional_ablation.R 10         B with 10 resamples instead
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args   <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args)) as.integer(args[1]) else 25L
N_EVAL <- 400L

cfg_on  <- load_config("config/config.yml")
cfg_off <- cfg_on
cfg_off$magnitude_conditional <- FALSE
cfg_off$intensity_conditional <- FALSE

tabs  <- load_tables(cfg_on$paths$mimiciv, cfg_on, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg_on)
tr    <- folds$stay_id[folds$split == "train"]
pri_on  <- layer1_priors(tabs, folds, cfg_on,  verbose = FALSE)
pri_off <- layer1_priors(tabs, folds, cfg_off, verbose = FALSE)

CONFIGS <- list(off = list(cfg = cfg_off, pri = pri_off),
                on  = list(cfg = cfg_on,  pri = pri_on))

.frame <- function(sg, md, which) {
  k <- CONFIGS[[which]]
  signal_frame(sg, md, tabs, k$cfg, priors_for(k$pri, sg, "final"), stay_ids = tr)
}
.fit <- function(d) mgcv::bam(attr(d, "formula"), data = d, family = stats::binomial(),
                              method = "fREML", discrete = TRUE,
                              nthreads = cfg_on$bam$nthreads %||% 1L,
                              na.action = stats::na.fail)
.de <- function(b) (b$null.deviance - b$deviance) / b$null.deviance
.cc <- function(b) {
  m <- try(mgcv::concurvity(b, full = FALSE)$observed, silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  diag(m) <- NA
  max(m, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# A. deviance explained, off vs on
# ---------------------------------------------------------------------------
# `intv` models isolate `lambda` (they carry no magnitude term) and `meas`
# models isolate `delta` (they carry no intervention term), so the two
# constructs can be read separately off one table.
cat("\n=== A. conditional constructs OFF vs ON: deviance explained ===\n")
cat("    same rows, same term count. `meas` rows isolate delta; `intv` rows\n")
cat("    isolate lambda; `full` rows carry both.\n\n")
cat(sprintf("%-26s %9s %9s %9s   %8s %8s\n",
            "model", "de.off", "de.on", "change", "cc.off", "cc.on"))
acc <- list()
for (sg in cfg_on$signals) {
  for (md in models_of(sg, cfg_on)) {
    a <- .fit(.frame(sg, md, "off")); b <- .fit(.frame(sg, md, "on"))
    da <- .de(a); db <- .de(b)
    acc[[length(acc) + 1L]] <- data.frame(signal = sg, model = md,
                                          de_off = da, de_on = db,
                                          stringsAsFactors = FALSE)
    cat(sprintf("%-26s %9.5f %9.5f %+8.1f%%   %8.4f %8.4f\n",
                paste0(sg, "/", md), da, db, 100 * (db - da) / da, .cc(a), .cc(b)))
  }
}
A <- do.call(rbind, acc)
cat("\n")
for (md in c("meas", "intv", "full")) {
  z <- A[A$model == md, , drop = FALSE]
  cat(sprintf("  %-6s (%2d models)  mean de %.5f -> %.5f  (%+.1f%%)   %s\n",
              md, nrow(z), mean(z$de_off), mean(z$de_on),
              100 * (sum(z$de_on) - sum(z$de_off)) / sum(z$de_off),
              if (md == "meas") "<- delta only" else if (md == "intv") "<- lambda only" else "both"))
}
cat(sprintf("  %-6s (%2d models)  mean de %.5f -> %.5f  (%+.1f%%)\n", "ALL",
            nrow(A), mean(A$de_off), mean(A$de_on),
            100 * (sum(A$de_on) - sum(A$de_off)) / sum(A$de_off)))

# ---------------------------------------------------------------------------
# B. attributability: does a fitted smooth mean anything?
# ---------------------------------------------------------------------------
# Refit on bootstrap resamples and track, at a FIXED set of real evaluation
# rows, the spread of each term's contribution against the spread of the total
# linear predictor. A large ratio means the DECOMPOSITION moves while the SUM
# does not: the model is well identified, its per-term attribution is not.
#
# Evaluation points are real rows held constant across resamples; a synthetic
# grid would put the smooths off the data manifold where any two functions can
# disagree freely, and the answer would be rigged.
#
# The two cases are the ones docs/v2_findings_20260827.md SS5.1 already
# measured, so the before/after is directly comparable to the numbers there.
CASES <- list(list(sg = "creatinine", md = "meas"),
              list(sg = "mbp",        md = "full"))

boot_ratio <- function(sg, md, which) {
  d <- .frame(sg, md, which)
  f <- attr(d, "formula")
  set.seed(cfg_on$seed)
  nd <- d[sample.int(nrow(d), min(N_EVAL, nrow(d))), , drop = FALSE]

  terms_mat <- NULL
  eta_mat <- matrix(NA_real_, nrow = B, ncol = nrow(nd))
  for (b in seq_len(B)) {
    set.seed(cfg_on$seed + b)
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    fit <- try(mgcv::bam(f, data = d[idx, , drop = FALSE], family = stats::binomial(),
                         method = "fREML", discrete = TRUE,
                         nthreads = cfg_on$bam$nthreads %||% 1L,
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
  sd_of <- function(M) mean(apply(M, 2, stats::sd, na.rm = TRUE), na.rm = TRUE)
  nm <- dimnames(terms_mat)[[3]]
  per <- vapply(seq_along(nm), function(j) sd_of(terms_mat[keep, , j, drop = TRUE]),
                numeric(1))
  names(per) <- nm
  eta <- sd_of(eta_mat[keep, , drop = FALSE])
  list(worst = max(per), worst_name = nm[which.max(per)], eta = eta,
       ratio = max(per) / eta)
}

cat("\n=== B. attributability: worst single term / total linear predictor ===\n")
cat(sprintf("    %d bootstrap refits, %d fixed evaluation rows. Lower is better:\n",
            B, N_EVAL))
cat("    1.0 means a term is no less stable than the prediction it contributes to.\n\n")
cat(sprintf("%-22s %28s %8s   %28s %8s\n", "model",
            "worst term (off)", "ratio", "worst term (on)", "ratio"))
for (cs in CASES) {
  a <- boot_ratio(cs$sg, cs$md, "off")
  b <- boot_ratio(cs$sg, cs$md, "on")
  cat(sprintf("%-22s %28s %8.2fx   %28s %8.2fx\n",
              paste0(cs$sg, "/", cs$md),
              substr(a$worst_name, 1, 28), a$ratio,
              substr(b$worst_name, 1, 28), b$ratio))
}
cat("\n")
