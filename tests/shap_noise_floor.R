# tests/shap_noise_floor.R ---------------------------------------------------
# HOW UNSTABLE IS AN ATTRIBUTION WHEN NOTHING ABOUT THE SPECIFICATION CHANGES?
#
# The attribution ladder found that two defensible specifications name a
# different top-contributing signal for about a quarter of patients and a
# different top domain for about a fifth. Those numbers have no scale. Nobody
# knows whether 24% is alarming or ordinary, because there is nothing to compare
# them against.
#
# THIS SCRIPT BUILDS THE MISSING REFERENCE POINT. It runs the `xgb_feat` arm
# twice, changing NOTHING but the random seed, and measures the disagreement
# between the two runs with exactly the metrics the ladder used. That is pure
# estimation noise for a method with no specification choice at all. If SHAP
# disagrees with itself as much as our ladder disagrees across specifications,
# then attribution instability is a property of the problem rather than a defect
# of this design, and every number in the ladder has to be read against it.
#
# WHY `xgb_feat` AND NOT `xgb_raw`. `xgb_feat`'s 99 columns roll up into exactly
# 31 groups -- the 19 signals and the 12 interventions -- which is the same
# partition the L matrix uses, so the aggregation is exact rather than a
# judgement call. `xgb_raw` carries `n_obs` and a missingness channel that map
# onto signals only loosely, so its groups would not be commensurable.
#
# WHY SHAP IS THE RIGHT COMPARATOR. TreeSHAP is additive by construction, so
# summing contributions within a group is EXACT rather than an approximation --
# which is not true of most attribution methods. And for a binary:logistic
# booster the contributions are on the margin, which is the log-odds scale, the
# same units as an L. So the two are directly commensurable in units even where
# they are not commensurable in meaning.
#
# THE ONE SEMANTIC CAVEAT, STATED RATHER THAN BURIED. `xgb_feat` is ONE joint
# model over every covariate, so its per-signal SHAP is a CONDITIONAL
# contribution given all the others. Our `L` is one MARGINAL model per bundle,
# which is the whole reason layer 2 exists to discount the redundancy. The two
# are different objects and their VALUES should not be compared. Their
# STABILITY can be, which is all this script does.
#
# WHAT "A DIFFERENT SEED" CHANGES. `.xgb_fit1()` draws the early-stopping
# validation split under the seed and passes the seed to xgboost. So a re-seed
# re-runs the identical procedure with a different random draw, which is the
# right notion of estimation noise for this pipeline.
#
# DOES NOT TOUCH THE TARGETS STORE, on purpose: it reads the parquet and rebuilds
# folds and priors itself, so it can run beside a `tar_make()`. It is also
# INDEPENDENT OF `bam.gamma`, because `xgb_feat`'s design is built from the
# layer-1 COVARIATES -- pi_hat, delta, lambda -- none of which is fitted by
# `bam`. A GAM refit cannot change any number here.
#
# Aggregates only (hard rule 1). The SHAP matrices are row-level, are saved to
# the run directory the way `l_oof.rds` is, and are never printed.
#
#   Rscript tests/shap_noise_floor.R
#   Rscript tests/shap_noise_floor.R 2000        # the seed offset
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(xgboost); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args   <- commandArgs(trailingOnly = TRUE)
OFFSET <- if (length(args) >= 1L && nzchar(args[1])) as.integer(args[1]) else 1000L
DELTA_ABS <- c(0, 0.02, 0.05, 0.10, 0.25, 0.50)
KS <- c(1L, 3L)

cfg    <- load_config("config/config.yml")
tabs   <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds  <- assign_folds(tabs$cohort, cfg)
tr     <- folds$stay_id[folds$split == "train"]
y      <- as.integer(tabs$cohort$mortality[match(tr, tabs$cohort$stay_id)])
fold_k <- folds$fold[match(tr, folds$stay_id)]
grp_k  <- patient_group_of(tabs$cohort, cfg, tr)
priors <- layer1_priors(tabs, folds, cfg, verbose = FALSE)
domains <- load_domains("config/domains.csv")

run <- new_run("shapfloor", cfg, note = sprintf(
  "xgb_feat SHAP estimation noise floor, seed offset %d, fits no GAM", OFFSET))

# --- build the two out-of-fold SHAP matrices --------------------------------
#
# One design per fold, reused for both seeds, so the ONLY difference between the
# two runs is the random draw. Building it twice would also re-run the identical
# deterministic code and cost half the runtime again.
cat("\n=== out-of-fold SHAP, two seeds, identical designs ===\n\n")
t0 <- feats <- NULL; t0 <- start_timer()
S1 <- S2 <- NULL; ref <- NULL
for (f in sort(unique(fold_k))) {
  X  <- xgb_design_feat(tabs, cfg, priors, tr, role = "oof", fold = f)
  if (is.null(ref)) {
    ref <- colnames(X)
    S1 <- matrix(NA_real_, length(tr), length(ref), dimnames = list(as.character(tr), ref))
    S2 <- S1
  } else if (!identical(colnames(X), ref)) {
    stop("fold ", f, " produced a different column set", call. = FALSE)
  }
  ho <- fold_k == f
  for (which in 1:2) {
    sd_use <- cfg$seed + f + if (which == 2L) OFFSET else 0L
    b <- .xgb_fit1(X[!ho, , drop = FALSE], y[!ho], cfg, seed = sd_use,
                   group = grp_k[!ho])$booster
    # predcontrib returns one column per feature plus a trailing BIAS column.
    ctr <- stats::predict(b, xgboost::xgb.DMatrix(X[ho, , drop = FALSE], missing = NA),
                          predcontrib = TRUE)
    stopifnot(ncol(ctr) == length(ref) + 1L)
    if (which == 1L) S1[ho, ] <- ctr[, seq_along(ref), drop = FALSE]
    else             S2[ho, ] <- ctr[, seq_along(ref), drop = FALSE]
  }
  cat(sprintf("  fold %d: %d held-out rows, %d features\n", f, sum(ho), length(ref)))
}
stopifnot(!anyNA(S1), !anyNA(S2))

# --- roll up to the 31 groups ------------------------------------------------
grp <- sub("__.*$", "", ref)
gnm <- sort(unique(grp))
cat(sprintf("\n  %d features -> %d groups (%d signals, %d interventions)\n",
            length(ref), length(gnm), sum(gnm %in% unlist(cfg$signals)),
            sum(!gnm %in% unlist(cfg$signals))))
roll <- function(S) {
  G <- matrix(0, nrow(S), length(gnm), dimnames = list(rownames(S), gnm))
  for (g in gnm) {
    j <- which(grp == g)
    G[, g] <- if (length(j) == 1L) S[, j] else rowSums(S[, j, drop = FALSE])
  }
  G
}
G1 <- roll(S1); G2 <- roll(S2)
sig_cols <- intersect(gnm, unlist(cfg$signals))
save_object(run, list(shap_seed1 = G1, shap_seed2 = G2, stay_id = tr,
                      groups = gnm, signal_groups = sig_cols), "shap_oof_groups")

# --- aggregations, mirroring tests/attribution_ties.R ------------------------
dm  <- domains$domain[match(sig_cols, domains$signal)]
dnm <- sort(unique(dm))
to_dom <- function(M) {
  D <- matrix(0, nrow(M), length(dnm), dimnames = list(rownames(M), dnm))
  for (k in dnm) {
    j <- sig_cols[dm == k]
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE])
  }
  D
}
LV <- list(signal = function(M) M[, sig_cols, drop = FALSE],
           domain = function(M) to_dom(M),
           all31  = function(M) M)

prep <- function(M) {
  ab <- abs(M); n <- nrow(ab)
  ord <- t(apply(-ab, 1, order))
  list(ab = ab, n = n,
       kth = function(k) ab[cbind(seq_len(n), ord[, k])],
       inS = function(k) { S <- matrix(FALSE, n, ncol(ab))
         S[cbind(rep(seq_len(n), times = k), as.vector(ord[, seq_len(k)]))] <- TRUE; S })
}
agree_k <- function(pa, pb, k, delta) {
  Sa <- pa$inS(k); Sb <- pb$inS(k); ka <- pa$kth(k); kb <- pb$kth(k)
  !((rowSums(Sb & !Sa & (pa$ab < ka - delta)) > 0) |
    (rowSums(Sa & !Sb & (pb$ab < kb - delta)) > 0))
}

cat("\n=== A. SHAP against itself: agreement under a seed change alone ===\n\n")
cat("    This is the NOISE FLOOR. No specification changed. Any disagreement\n")
cat("    here is estimation noise in a method with no specification choice.\n")
TA <- list()
for (lv in names(LV)) {
  pa <- prep(LV[[lv]](G1)); pb <- prep(LV[[lv]](G2))
  cat(sprintf("\n  --- %s level (%d columns) ---\n\n", lv, ncol(pa$ab)))
  cat(sprintf("  %-6s", "k"))
  for (d in DELTA_ABS) cat(sprintf(" %9s", sprintf("d=%.2f", d)))
  cat("\n")
  for (k in KS) {
    r <- vapply(DELTA_ABS, function(d) mean(agree_k(pa, pb, k, d)), numeric(1))
    TA[[length(TA) + 1L]] <- data.frame(comparison = "shap_seed1_vs_seed2",
      level = lv, k = k, delta_kind = "absolute_nats", delta = DELTA_ABS,
      agree = round(r, 5), n_patients = pa$n, stringsAsFactors = FALSE)
    cat(sprintf("  top-%-2d", k)); for (v in r) cat(sprintf(" %9.4f", v)); cat("\n")
  }
}
save_table(run, do.call(rbind, TA), "shap_noise_agreement", subdir = "diagnostics")

cat("\n=== B. sign flip between the two seeds, gated ===\n\n")
TB <- list()
cat(sprintf("  %-8s %10s", "level", "cells"))
for (d in DELTA_ABS) cat(sprintf(" %9s", sprintf("d=%.2f", d)))
cat("\n")
for (lv in names(LV)) {
  a <- LV[[lv]](G1); b <- LV[[lv]](G2)
  fl <- sign(a) != sign(b); mx <- pmax(abs(a), abs(b))
  r <- vapply(DELTA_ABS, function(d) mean(fl & mx > d), numeric(1))
  TB[[length(TB) + 1L]] <- data.frame(comparison = "shap_seed1_vs_seed2", level = lv,
    n_cells = length(fl), delta = DELTA_ABS, flip = round(r, 5),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-8s %10d", lv, length(fl)))
  for (v in r) cat(sprintf(" %9.4f", v)); cat("\n")
}
save_table(run, do.call(rbind, TB), "shap_noise_signflip", subdir = "diagnostics")

cat("\n=== C. resolution and evidence budget of a SHAP explanation ===\n\n")
cat("    Directly comparable with the L arms' tables: the same statistics on\n")
cat("    the same patients, in the same log-odds units.\n\n")
cat(sprintf("  %-8s %-8s %10s %10s %10s %12s\n", "level", "seed", "budget",
            "hhi", "max share", "uniq@d=0.10"))
TC <- list()
for (lv in names(LV)) {
  for (s in 1:2) {
    M <- LV[[lv]](if (s == 1L) G1 else G2)
    ab <- abs(M); tot <- rowSums(ab); sh <- ab / pmax(tot, .Machine$double.eps)
    lead <- apply(ab, 1, max)
    TC[[length(TC) + 1L]] <- data.frame(level = lv, seed = s, n_columns = ncol(M),
      budget_median = round(stats::median(tot), 4),
      hhi_median = round(stats::median(rowSums(sh^2)), 5),
      max_share_median = round(stats::median(apply(sh, 1, max)), 5),
      frac_unique_leader_010 = round(mean(rowSums(ab >= lead - 0.10) == 1L), 5),
      stringsAsFactors = FALSE)
    z <- TC[[length(TC)]]
    cat(sprintf("  %-8s %-8d %10.4f %10.5f %10.5f %12.4f\n", lv, s,
                z$budget_median, z$hhi_median, z$max_share_median,
                z$frac_unique_leader_010))
  }
}
save_table(run, do.call(rbind, TC), "shap_budget", subdir = "diagnostics")

cat(sprintf("\n  %d boosters in %.1f minutes.\n", 2 * length(unique(fold_k)),
            t0()$elapsed_sec / 60))
cat("  READ THE LADDER AGAINST THIS. tests/attribution_ties.R reported, for\n")
cat("  cond vs cond_ti_all at domain level and top-1, 0.8121 at delta = 0 and\n")
cat("  0.8266 at delta = 0.10. If the seed-change floor above is comparable,\n")
cat("  attribution instability is a property of the problem and not of the\n")
cat("  specification choice.\n")

finalize_run(run, extra = list(seed_offset = OFFSET, n_features = length(ref),
                               n_groups = length(gnm)))
cat(sprintf("\nwritten: %s\n", run$path))
