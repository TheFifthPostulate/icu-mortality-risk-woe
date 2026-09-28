# tests/attribution_ties.R ---------------------------------------------------
# TIE-TOLERANT ATTRIBUTION AGREEMENT.
#
# tests/coupling_attribution.R reported exact top-1 and top-3 agreement between
# specifications, and tests/attribution_margin.R showed the disagreement is
# concentrated where the leader is not clear. This script closes that properly:
# it stops treating a reordering of two near-equal contributions as a
# disagreement at all.
#
# THE DEFINITION, AND WHY THIS ONE. Two specifications are said to agree at
# rank k with tie tolerance delta when NEITHER of them strongly prefers its own
# choice to the other's. Formally, with S_a and S_b the top-k sets and kth_a the
# k-th largest |L| under arm a:
#
#   a real disagreement exists if some g in S_b but not S_a has
#   |L_a[g]| < kth_a - delta   (arm a says b's pick is materially worse)
#   or symmetrically with a and b exchanged.
#
# Agreement otherwise. At delta = 0 this reduces to exact set equality, so the
# grid contains the original number and the reader can see the whole curve.
#
# It is SYMMETRIC on purpose. A one-sided version -- "is b's pick close under
# a?" -- would call it agreement whenever arm a happens to be flat, even if arm
# b is emphatic that a's pick is wrong. Requiring both directions means a
# disagreement is only forgiven when both models are genuinely undecided.
#
# TWO KINDS OF DELTA, because patients carry different amounts of evidence.
# ABSOLUTE delta is in nats and is the natural unit of the LLR. RELATIVE delta
# is a fraction of that patient's own total absolute evidence, which is the
# scale-free version and is fairer to patients with little evidence overall.
# Both are reported.
#
# THE RESOLUTION DIAGNOSTIC IS THE OTHER HALF, AND IT IS NOT A COMPARISON.
# `tie_set_size` counts, per patient, how many contributions sit within delta of
# their largest. That is a property of ONE arm and it bounds what ANY model can
# claim: if the typical patient has four signals within a tenth of a nat of the
# leader, then "the top signal for this patient" is not a supportable statement
# however stable it happens to be across specifications.
#
# FITS NOTHING. Reads the cached ladder from tests/coupling_attribution.R.
#
# Aggregates only (hard rule 1): agreement rates, counts and set sizes.
#
#   Rscript tests/attribution_ties.R out/runs/coupattr_20260905T184528
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args     <- commandArgs(trailingOnly = TRUE)
ladder_d <- if (length(args) >= 1L && nzchar(args[1])) args[1] else
              latest_run("coupattr", require_complete = FALSE)
lad <- readRDS(file.path(ladder_d, "tables", "l_oof_ladder.rds"))
A   <- lad$arms

cfg     <- tar_read(cfg)
tr      <- tar_read(train_ids)
domains <- load_domains("config/domains.csv")
sigs    <- unlist(cfg$signals)
stopifnot(identical(as.character(lad$stay_id), as.character(tr)))

run <- new_run("attrties", cfg, note = sprintf(
  "tie-tolerant attribution agreement, ladder from %s, fits nothing",
  basename(ladder_d)))

# Declared before anything is computed.
DELTA_ABS <- c(0, 0.02, 0.05, 0.10, 0.25, 0.50)     # nats
DELTA_REL <- c(0, 0.01, 0.02, 0.05, 0.10)           # share of the patient's own |L| total
KS        <- c(1L, 3L)

CMP <- list(
  c("cond",          "cond_ti_all"),
  c("cond",          "cond_ti_trend"),
  c("cond_ti_trend", "cond_ti_all"),
  c("meas",          "cond"),
  c("meas",          "cond_ti_all"),
  c("full",          "full_ti_all"),
  c("full",          "full_ti_trend"))
ARMS <- unique(unlist(CMP))

# --- aggregations -----------------------------------------------------------
dm  <- domains$domain[match(sigs, domains$signal)]
dnm <- sort(unique(dm))
to_dom <- function(M) {
  D <- matrix(0, nrow(M), length(dnm), dimnames = list(rownames(M), dnm))
  for (k in dnm) {
    j <- which(dm == k)
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE])
  }
  D
}
LEVELS <- list(signal = function(M) M, domain = to_dom)

# --- precompute, once per arm per level -------------------------------------
#
# `ord` is the column order by descending |L|, so the k-th largest value and the
# top-k membership both fall out by indexing rather than by a second pass over
# the rows. One `apply` per arm per level instead of one per comparison per k
# per delta, which is the difference between seconds and an hour.
prep <- function(M) {
  ab <- abs(M); n <- nrow(ab)
  ord <- t(apply(-ab, 1, order))
  kth <- function(k) ab[cbind(seq_len(n), ord[, k])]
  inS <- function(k) {
    S <- matrix(FALSE, n, ncol(ab))
    S[cbind(rep(seq_len(n), times = k), as.vector(ord[, seq_len(k)]))] <- TRUE
    S
  }
  list(ab = ab, ord = ord, kth = kth, inS = inS, total = rowSums(ab), n = n)
}

cat("\n=== precomputing orderings ===\n")
P <- list()
for (lv in names(LEVELS)) {
  P[[lv]] <- list()
  for (a in ARMS) P[[lv]][[a]] <- prep(LEVELS[[lv]](A[[a]]))
  cat(sprintf("  %s level: %d arms, %d columns\n", lv, length(ARMS),
              ncol(P[[lv]][[ARMS[1]]]$ab)))
}

#' Tie-tolerant agreement at rank k. `delta` is a scalar or a per-row vector.
agree_k <- function(pa, pb, k, delta) {
  Sa <- pa$inS(k); Sb <- pb$inS(k)
  ka <- pa$kth(k); kb <- pb$kth(k)
  # arm a is materially unhappy about something b put in its top-k
  viol_a <- rowSums(Sb & !Sa & (pa$ab < ka - delta)) > 0
  viol_b <- rowSums(Sa & !Sb & (pb$ab < kb - delta)) > 0
  !(viol_a | viol_b)
}

# --- A. absolute tie tolerance ----------------------------------------------
cat("\n=== A. agreement with an ABSOLUTE tie tolerance, in nats ===\n\n")
cat("    delta = 0 is exact set equality, the number reported previously.\n")
cat("    A disagreement is only counted when at least one arm says the other's\n")
cat("    pick is worse by more than delta.\n")
TA <- list()
for (lv in names(LEVELS)) {
  for (k in KS) {
    cat(sprintf("\n  --- %s level, top-%d ---\n\n", lv, k))
    cat(sprintf("  %-15s %-15s", "arm a", "arm b"))
    for (d in DELTA_ABS) cat(sprintf(" %9s", sprintf("d=%.2f", d)))
    cat("\n")
    for (p in CMP) {
      pa <- P[[lv]][[p[1]]]; pb <- P[[lv]][[p[2]]]
      r <- vapply(DELTA_ABS, function(d) mean(agree_k(pa, pb, k, d)), numeric(1))
      TA[[length(TA) + 1L]] <- data.frame(level = lv, k = k, arm_a = p[1], arm_b = p[2],
        delta_kind = "absolute_nats", delta = DELTA_ABS, agree = round(r, 5),
        n_patients = pa$n, stringsAsFactors = FALSE)
      cat(sprintf("  %-15s %-15s", p[1], p[2]))
      for (v in r) cat(sprintf(" %9.4f", v))
      cat("\n")
    }
  }
}
save_table(run, do.call(rbind, TA), "tie_agreement_absolute", subdir = "diagnostics")

# --- B. relative tie tolerance ----------------------------------------------
cat("\n=== B. agreement with a RELATIVE tie tolerance ===\n\n")
cat("    delta is a share of that patient's own total absolute evidence, so a\n")
cat("    patient with little evidence is not held to a nats threshold built for\n")
cat("    a patient with a lot.\n")
TB <- list()
for (lv in names(LEVELS)) {
  for (k in KS) {
    cat(sprintf("\n  --- %s level, top-%d ---\n\n", lv, k))
    cat(sprintf("  %-15s %-15s", "arm a", "arm b"))
    for (d in DELTA_REL) cat(sprintf(" %9s", sprintf("r=%.2f", d)))
    cat("\n")
    for (p in CMP) {
      pa <- P[[lv]][[p[1]]]; pb <- P[[lv]][[p[2]]]
      r <- vapply(DELTA_REL, function(d)
        mean(agree_k(pa, pb, k, d * pa$total)), numeric(1))
      TB[[length(TB) + 1L]] <- data.frame(level = lv, k = k, arm_a = p[1], arm_b = p[2],
        delta_kind = "relative_share", delta = DELTA_REL, agree = round(r, 5),
        n_patients = pa$n, stringsAsFactors = FALSE)
      cat(sprintf("  %-15s %-15s", p[1], p[2]))
      for (v in r) cat(sprintf(" %9.4f", v))
      cat("\n")
    }
  }
}
save_table(run, do.call(rbind, TB), "tie_agreement_relative", subdir = "diagnostics")

# --- C. resolution: how many contributions are in contention at all? --------
#
# NOT a comparison. A property of one arm, and the ceiling on what any per-patient
# explanation can claim.
cat("\n=== C. resolution: contributions within delta of the leader ===\n\n")
cat("    median and 90th percentile of the tie-set size, and the share of\n")
cat("    patients with an UNAMBIGUOUS leader (tie set of exactly one).\n\n")
TC <- list()
for (lv in names(LEVELS)) {
  ncol_lv <- ncol(P[[lv]][[ARMS[1]]]$ab)
  cat(sprintf("  --- %s level (%d columns) ---\n\n", lv, ncol_lv))
  cat(sprintf("  %-16s", "arm"))
  for (d in DELTA_ABS) cat(sprintf(" %14s", sprintf("d=%.2f med/uniq", d)))
  cat("\n")
  for (a in ARMS) {
    pa <- P[[lv]][[a]]; lead <- pa$kth(1L)
    cat(sprintf("  %-16s", a))
    for (d in DELTA_ABS) {
      sz <- rowSums(pa$ab >= lead - d)
      TC[[length(TC) + 1L]] <- data.frame(level = lv, arm = a, n_columns = ncol_lv,
        delta = d, tieset_median = stats::median(sz),
        tieset_mean = round(mean(sz), 3),
        tieset_p90 = unname(stats::quantile(sz, 0.90)),
        frac_unique_leader = round(mean(sz == 1L), 5), stringsAsFactors = FALSE)
      z <- TC[[length(TC)]]
      cat(sprintf(" %8.0f/%5.3f", z$tieset_median, z$frac_unique_leader))
    }
    cat("\n")
  }
  cat("\n")
}
save_table(run, do.call(rbind, TC), "tie_set_size", subdir = "diagnostics")

# --- D. the headline curve, one table --------------------------------------
cat("=== D. the headline: how fast does disagreement fall with tie tolerance? ===\n\n")
cat(sprintf("%-15s %-15s %-8s %-3s %9s %9s %9s %9s\n", "arm a", "arm b",
            "level", "k", "d=0", "d=0.05", "d=0.10", "d=0.25"))
TD <- do.call(rbind, TA)
for (p in CMP) {
  for (lv in c("domain", "signal")) {
    for (k in KS) {
      z <- TD[TD$arm_a == p[1] & TD$arm_b == p[2] & TD$level == lv & TD$k == k, ]
      g <- function(d) z$agree[z$delta == d]
      cat(sprintf("%-15s %-15s %-8s %-3d %9.4f %9.4f %9.4f %9.4f\n",
                  p[1], p[2], lv, k, g(0), g(0.05), g(0.10), g(0.25)))
    }
  }
}

finalize_run(run, extra = list(ladder = basename(ladder_d),
                               delta_abs = paste(DELTA_ABS, collapse = ","),
                               delta_rel = paste(DELTA_REL, collapse = ",")))
cat(sprintf("\nwritten: %s\n", run$path))
