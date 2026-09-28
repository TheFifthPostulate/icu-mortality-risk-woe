# tests/attr_metrics.R --------------------------------------------------------
# THE CHEAP CONSUMER. Reads cached attribution replicates, computes every metric
# at every contrast of the reproducibility hierarchy AS A DISTRIBUTION, writes
# the tables. FITS NOTHING, and it must stay that way: the whole point of the
# generator/consumer split in `docs/attribution_analysis_plan_20260906.md`
# section 10 is that a metric definition can be changed and re-run in minutes,
# without regenerating a single replicate.
#
# --- WHAT CHANGED ON 2026-09-06, AND WHY IT IS THE POINT OF THE FILE ---------
#
# EVERY LEVEL USED TO BE ONE NUMBER, COMPUTED FROM `mats[[1]]` AGAINST
# `mats[[2]]`. Forty replicates were generated per method and 778 of the 780
# pairs were thrown away. The comment justifying it said averaging over the
# pairs would be "a variance reduction on a number whose spread is itself the
# finding", which is right about the averaging and wrong about the remedy: the
# answer to "the spread is the finding" is to REPORT THE SPREAD, not to report
# one arbitrary draw from it.
#
# So every contrast is now a distribution, summarised by `attr_stat_summary()`
# with a fixed column set -- min, p05, p10, q1, median, mean, q3, p90, p95, max,
# IQR, sd, and two skew indicators. Three consequences, all of them things a
# single number could not say:
#
#   THE WORST CASE IS VISIBLE. `p95` is the robustness criterion and `max` is
#   the worst observed pair. They are reported together and only one of them is
#   ever a threshold: a maximum grows mechanically with B, so a criterion on it
#   would tighten every time more replicates were generated.
#
#   THE SKEW IS VISIBLE. `mean_minus_median` above zero means a heavy right
#   tail, which means the median understates how bad a bad resample is. That is
#   exactly the case a reproducibility claim has to survive.
#
#   TWO LEVELS CAN BE COMPARED WITHOUT DIVIDING. `attr_dominance()` reports
#   P(level 4 > level 3) and the median gap on the disagreement scale, which
#   has no denominator to explode and does not force two distributions with
#   opposite readings onto one axis.
#
# --- THE SIX CONTRASTS, AND WHERE EACH COMES FROM ---------------------------
#
#   L1   nothing differs      the nats budget, which sets every other scale
#   L2   seed only            SHAP, at the original sample and now also within
#                             each shared bag
#   L3   bag only             sampling noise, seed held. BOTH FAMILIES, ON THE
#                             SAME BAGS, which is what makes it the common
#                             currency
#   L3T  bag AND seed         operational reproducibility. Level-2 uncertainty
#                             propagated into level 3 BY DESIGN rather than by
#                             a variance-components model
#   L3P  posterior draw only  LLR estimation noise conditional on the sample. A
#                             DIFFERENT ESTIMAND from L3 and never pooled with it
#   L4   method only          specification, with every coordinate held. On the
#                             shared bags this is a DISTRIBUTION over bags, not
#                             the single number the two original fits gave
#
# --- THE CORRECTNESS GATE ---------------------------------------------------
#
# `--gate` recomputes the tables that `tests/coupling_attribution.R`,
# `tests/attribution_ties.R` and `tests/shap_noise_floor.R` already wrote, from
# THE SAME CACHED MATRICES, using only the consolidated library, and asserts
# they agree. If the consolidation changed a number, the consolidation is wrong
# and everything downstream of it is uninterpretable.
#
# THE GATE READS PRE-GAMMA RUNS ON PURPOSE, and no number it prints is a current
# result: those runs were built before `bam.gamma` became 1.5, which is exactly
# what makes them right here. The gate is a statement about CODE EQUIVALENCE and
# needs the inputs the old code actually ran on.
#
# THE THREE RUN PATHS MOVED INTO `config/attribution_eval.yml` ON 2026-09-06,
# with a CONTENT KEY beside each. They were script literals at lines 54 to 56 of
# this file, which made them undiffable and unauditable; and a pinned reference
# with no integrity check is a reference that can be regenerated underneath the
# gate without the gate noticing. The content key is a hash of every `.rds` in
# the run's `diagnostics/` and `tables/`, so a legacy run that moves stops the
# gate and names itself.
#
# Aggregates only (hard rule 1).
#
#   Rscript tests/attr_metrics.R --gate
#   Rscript tests/attr_metrics.R                       # latest attrgen run
#   Rscript tests/attr_metrics.R out/runs/attrgen_...
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
GATE <- "--gate" %in% args
pos  <- args[!grepl("^--", args)]
ecfg <- yaml::read_yaml("config/attribution_eval.yml")

# ============================================================================
# THE GATE
# ============================================================================

GATE_LADDER <- as.character(cfg_req(ecfg, "gate", "ladder"))
GATE_TIES   <- as.character(cfg_req(ecfg, "gate", "ties"))
GATE_FLOOR  <- as.character(cfg_req(ecfg, "gate", "floor"))

#' The integrity key of a pinned reference run.
#'
#' A DIFFERENT GUARD FROM THE LADDER'S DESIGN KEY, and the difference is the
#' whole reason both exist. A design key asks "was this built by the current
#' code?" and for these three runs the intended answer is NO. What must be
#' checked instead is that the REFERENCE ITSELF HAS NOT MOVED: if a legacy run
#' is regenerated, edited, or gains a table, the gate would go on comparing the
#' library against something other than what the config names, and would pass.
.run_content_key <- function(d) {
  fs <- c(sort(list.files(file.path(d, "diagnostics"), pattern = "[.]rds$",
                          full.names = TRUE)),
          sort(list.files(file.path(d, "tables"), pattern = "[.]rds$",
                          full.names = TRUE)))
  attr_key_hash(unname(tools::md5sum(fs)))
}

.guard_gate_run <- function(which, d) {
  if (!dir.exists(d)) {
    stop("config/attribution_eval.yml names gate.", which, " = ", d,
         "\nbut no such directory exists. The gate's references are PINNED, ",
         "not resolved with latest_run(): a gate that follows the newest run ",
         "is a comparison of the library with itself. Restore the run or ",
         "repin the key.", call. = FALSE)
  }
  want <- as.character(cfg_req(ecfg, "gate", "content_keys", which))
  got  <- .run_content_key(d)
  if (!identical(got, want)) {
    stop("gate.", which, " (", basename(d), ") has content key ", got,
         " but config/attribution_eval.yml declares ", want, ".\nThe pinned ",
         "reference has MOVED -- regenerated, edited, or a table added. The ",
         "gate would otherwise go on passing while comparing the library ",
         "against something the config does not name. Establish which branch ",
         "moved before repinning: a legacy run should not change.", call. = FALSE)
  }
  cat(sprintf("  %-7s %-34s content key %s  ok\n", which, basename(d), got))
  invisible(TRUE)
}

.legacy <- function(dir, name) {
  p <- file.path(dir, "diagnostics", paste0(name, ".rds"))
  if (!file.exists(p)) stop("gate: no legacy table ", name, " in ", dir, call. = FALSE)
  readRDS(p)
}

GATE_STATE <- new.env(parent = emptyenv())
GATE_STATE$rows <- list()

#' Compare a recomputed vector against the legacy one.
#'
#' `tol` is set from the LEGACY table's own rounding, not from a general sense
#' of what is close enough: `evidence_budget` rounded budgets to four decimals
#' and everything else to five, so demanding 1e-5 there would fail on the
#' rounding rather than on the statistic.
.cmp <- function(what, got, want, tol = 1e-5) {
  ok <- length(got) == length(want)
  d  <- if (ok) max(abs(as.numeric(got) - as.numeric(want)), na.rm = TRUE) else NA_real_
  pass <- ok && is.finite(d) && d <= tol
  GATE_STATE$rows[[length(GATE_STATE$rows) + 1L]] <- data.frame(
    check = what, n = length(want), max_abs_diff = d, tol = tol,
    pass = pass, stringsAsFactors = FALSE)
  cat(sprintf("  %-52s %6d values  max|d| = %.3e  %s\n", what, length(want),
              if (is.finite(d)) d else NA_real_, if (pass) "PASS" else "*** FAIL ***"))
  invisible(pass)
}

run_gate <- function() {
  cat("\n=== CORRECTNESS GATE: does R/14_attribution_eval.R reproduce the",
      "\n    scattered tables from the identical cached matrices? ===\n\n")
  cat("  the three references are PINNED IN CONFIG and PRE-GAMMA by design;\n")
  cat("  no number below is a current result.\n\n")
  .guard_gate_run("ladder", GATE_LADDER)
  .guard_gate_run("ties",   GATE_TIES)
  .guard_gate_run("floor",  GATE_FLOOR)
  cat("\n")

  lad <- readRDS(file.path(GATE_LADDER, "tables", "l_oof_ladder.rds"))
  A   <- lad$arms
  cfg <- tar_read(cfg)
  tr  <- tar_read(train_ids)
  domains <- load_domains("config/domains.csv")
  stopifnot(identical(as.character(lad$stay_id), as.character(tr)))

  # The measured mask. The one thing the gate cannot get from a run directory,
  # because `tests/coupling_attribution.R` never saved it. A parquet read, so
  # this branch is the only part of the file that is not instant.
  cat("  loading the measured mask from the targets store ...\n")
  meas_ok <- measured_matrix(tar_read(tabs), cfg, tr)

  DELTA_ABS <- c(0, 0.02, 0.05, 0.10, 0.25, 0.50)
  DELTA_REL <- c(0, 0.01, 0.02, 0.05, 0.10)
  TAUS      <- c(0, 0.05, 0.10, 0.25, 0.50)
  KS        <- c(1L, 3L)
  CMP <- list(c("cond", "cond_ti_all"), c("cond", "cond_ti_trend"),
              c("cond_ti_trend", "cond_ti_all"), c("meas", "cond"),
              c("meas", "cond_ti_all"), c("full", "full_ti_all"),
              c("full", "full_ti_trend"))
  ARMS <- unique(unlist(CMP))

  LV <- list(signal = function(M) M,
             domain = function(M) attr_to_domain(M, domains))
  P <- lapply(LV, function(f) lapply(stats::setNames(ARMS, ARMS),
                                     function(a) attr_prep(f(A[[a]]))))

  # --- 1. tie-tolerant agreement, both grids -------------------------------
  cat("\n  --- tie-tolerant agreement (attribution_ties.R) ---\n")
  for (kind in c("absolute", "relative")) {
    lg <- .legacy(GATE_TIES, paste0("tie_agreement_", kind))
    grid <- if (kind == "absolute") DELTA_ABS else DELTA_REL
    got <- want <- numeric(0)
    for (lv in names(LV)) for (k in KS) for (p in CMP) {
      pa <- P[[lv]][[p[1]]]; pb <- P[[lv]][[p[2]]]
      r <- vapply(grid, function(d) mean(attr_agree_k(
        pa, pb, k, if (kind == "relative") d * pa$total else d)), numeric(1))
      z <- lg[lg$level == lv & lg$k == k & lg$arm_a == p[1] & lg$arm_b == p[2], ]
      z <- z[match(grid, z$delta), ]
      got <- c(got, round(r, 5)); want <- c(want, z$agree)
    }
    .cmp(sprintf("attr_agree_k vs tie_agreement_%s", kind), got, want)
  }

  # --- 1b. the bitmask fast path against attr_agree_k(delta = 0) -----------
  #
  # ADDED 2026-09-06 WITH THE FAST PATH ITSELF. `attr_topk()` encodes each
  # patient's top-k set as an integer so that strict agreement over hundreds of
  # pairs is one vector comparison instead of two n x p logical masks. That is
  # a claim about tie-breaking and set semantics being identical to
  # `attr_agree_k(delta = 0)`, and a claim of that kind belongs in the gate
  # rather than in a comment: the distribution tables are built entirely on the
  # fast path, so a divergence here would move every headline number.
  got <- want <- numeric(0)
  for (lv in names(LV)) {
    TK <- lapply(stats::setNames(ARMS, ARMS), function(a) attr_topk(LV[[lv]](A[[a]]), KS))
    for (k in KS) for (p in CMP) {
      got  <- c(got,  attr_topk_agree(TK[[p[1]]], TK[[p[2]]], k))
      want <- c(want, mean(attr_agree_k(P[[lv]][[p[1]]], P[[lv]][[p[2]]], k, 0)))
    }
  }
  .cmp("attr_topk_agree vs attr_agree_k(delta = 0)", got, want, tol = 1e-12)

  # --- 1c. the boundary-tie counterexample, which the ladder cannot supply -
  #
  # ADDED 2026-09-09 (review finding A3). Check 1b passed on the anchor ladder
  # while the two functions DIFFERED in definition, because no patient there
  # happened to have an exact tie at the third-place boundary: a test on
  # historical inputs is a test of the inputs. This is the review's own
  # four-signal example -- A = (3, 2, 1, 1) against B = (1, 1, 3, 2) -- on
  # which the old bitmask said 0 and the tie-aware rule says 1, plus twenty
  # random tie-heavy matrices. `tests/attr_eval_unit.R` carries the same check
  # without needing the pinned runs.
  A_t <- rbind(c(3, 2, 1, 1), c(5, 4, 3, 2)); B_t <- rbind(c(1, 1, 3, 2), c(5, 4, 3, 2))
  got <- attr_topk_agree(attr_topk(A_t, 3L), attr_topk(B_t, 3L), 3L)
  want <- mean(attr_agree_k(attr_prep(A_t), attr_prep(B_t), 3L, 0))
  set.seed(20260909)
  for (trial in 1:20) {
    X_t <- matrix(sample(0:3, 300 * 11, TRUE), 300) * sample(c(-1, 1), 300 * 11, TRUE)
    Y_t <- X_t; Y_t[sample(length(X_t), 300)] <- sample(0:3, 300, TRUE)
    for (k in KS) {
      got  <- c(got,  attr_topk_agree(attr_topk(X_t, k), attr_topk(Y_t, k), k))
      want <- c(want, mean(attr_agree_k(attr_prep(X_t), attr_prep(Y_t), k, 0)))
    }
  }
  .cmp("attr_topk_agree vs attr_agree_k on boundary ties (synthetic)", got, want,
       tol = 1e-12)
  if (want[1] != 1) stop("gate: the tie counterexample no longer agrees under ",
                         "attr_agree_k; the definition itself has moved.", call. = FALSE)

  # --- 2. resolution --------------------------------------------------------
  lg <- .legacy(GATE_TIES, "tie_set_size")
  got <- want <- numeric(0)
  for (lv in names(LV)) for (a in ARMS) {
    t <- attr_tieset(P[[lv]][[a]], DELTA_ABS)
    z <- lg[lg$level == lv & lg$arm == a, ]; z <- z[match(DELTA_ABS, z$delta), ]
    got <- c(got, t$tieset_median, t$tieset_mean, t$tieset_p90, t$frac_unique_leader)
    want <- c(want, z$tieset_median, z$tieset_mean, z$tieset_p90, z$frac_unique_leader)
  }
  .cmp("attr_tieset vs tie_set_size", got, want)

  # --- 3. gated sign flip ---------------------------------------------------
  cat("\n  --- ladder tables (coupling_attribution.R) ---\n")
  lg <- .legacy(GATE_LADDER, "sign_flip_gated")
  got <- want <- numeric(0)
  for (p in CMP) {
    f <- attr_sign_flip(A[[p[1]]], A[[p[2]]], keep = meas_ok, taus = TAUS)
    z <- lg[lg$a == p[1] & lg$b == p[2], ]; z <- z[match(TAUS, z$tau), ]
    got <- c(got, f$flip, f$n_cells[1]); want <- c(want, z$flip, z$n_cells[1])
  }
  .cmp("attr_sign_flip vs sign_flip_gated", got, want)

  # --- 4. within-patient ranking -------------------------------------------
  lg <- .legacy(GATE_LADDER, "attribution_ranking")
  got <- want <- numeric(0)
  for (p in CMP) {
    r <- attr_rank_agreement(A[[p[1]]], A[[p[2]]])
    z <- lg[lg$a == p[1] & lg$b == p[2], ]
    got <- c(got, r$top1_agree, r$top3_agree, r$rho_median, r$rho_p10)
    want <- c(want, z$top1_agree, z$top3_agree, z$rho_median, z$rho_p10)
  }
  .cmp("attr_rank_agreement vs attribution_ranking", got, want)

  # --- 5. nats budget -------------------------------------------------------
  lg <- .legacy(GATE_LADDER, "evidence_budget")
  got <- want <- numeric(0)
  for (a in names(A)) {
    b <- attr_budget(A[[a]]); z <- lg[lg$arm == a, ]
    got <- c(got, b$budget_median, b$budget_p90, b$hhi_median, b$max_share_median)
    want <- c(want, z$budget_median, z$budget_p90, z$hhi_median, z$max_share_median)
  }
  # 1e-4: `evidence_budget` rounded the two budget columns to four decimals.
  .cmp("attr_budget vs evidence_budget", got, want, tol = 1e-4)

  # --- 6. domain level ------------------------------------------------------
  lg <- .legacy(GATE_LADDER, "domain_attribution")
  dmeas <- attr_to_domain(meas_ok * 1, domains) > 0
  got <- want <- numeric(0)
  for (p in CMP) {
    Da <- attr_to_domain(A[[p[1]]], domains); Db <- attr_to_domain(A[[p[2]]], domains)
    f <- attr_sign_flip(Da, Db, keep = dmeas, taus = c(0, 0.10))
    r <- attr_rank_agreement(Da, Db)
    z <- lg[lg$a == p[1] & lg$b == p[2], ]
    got <- c(got, f$flip, r$top1_agree, r$rho_median)
    want <- c(want, z$flip_t0, z$flip_t10, z$top1_agree, z$rho_median)
  }
  .cmp("attr_to_domain + flip/rank vs domain_attribution", got, want)

  # --- 7. the SHAP floor ----------------------------------------------------
  cat("\n  --- SHAP floor tables (shap_noise_floor.R) ---\n")
  sh <- readRDS(file.path(GATE_FLOOR, "tables", "shap_oof_groups.rds"))
  G1 <- sh$shap_seed1; G2 <- sh$shap_seed2
  sc <- sh$signal_groups
  SLV <- list(signal = function(M) M[, sc, drop = FALSE],
              domain = function(M) attr_to_domain(M[, sc, drop = FALSE], domains),
              all31  = function(M) M)

  lg <- .legacy(GATE_FLOOR, "shap_noise_agreement")
  got <- want <- numeric(0)
  for (lv in names(SLV)) {
    pa <- attr_prep(SLV[[lv]](G1)); pb <- attr_prep(SLV[[lv]](G2))
    for (k in KS) {
      r <- vapply(DELTA_ABS, function(d) mean(attr_agree_k(pa, pb, k, d)), numeric(1))
      z <- lg[lg$level == lv & lg$k == k, ]; z <- z[match(DELTA_ABS, z$delta), ]
      got <- c(got, round(r, 5)); want <- c(want, z$agree)
    }
  }
  .cmp("attr_agree_k vs shap_noise_agreement", got, want)

  lg <- .legacy(GATE_FLOOR, "shap_noise_signflip")
  got <- want <- numeric(0)
  for (lv in names(SLV)) {
    f <- attr_sign_flip(SLV[[lv]](G1), SLV[[lv]](G2), keep = NULL, taus = DELTA_ABS)
    z <- lg[lg$level == lv, ]; z <- z[match(DELTA_ABS, z$delta), ]
    got <- c(got, f$flip, f$n_cells[1]); want <- c(want, z$flip, z$n_cells[1])
  }
  .cmp("attr_sign_flip vs shap_noise_signflip", got, want)

  lg <- .legacy(GATE_FLOOR, "shap_budget")
  got <- want <- numeric(0)
  for (lv in names(SLV)) for (s in 1:2) {
    M <- SLV[[lv]](if (s == 1L) G1 else G2)
    b <- attr_budget(M); u <- attr_tieset(attr_prep(M), 0.10)
    z <- lg[lg$level == lv & lg$seed == s, ]
    got <- c(got, b$budget_median, b$hhi_median, b$max_share_median, u$frac_unique_leader)
    want <- c(want, z$budget_median, z$hhi_median, z$max_share_median,
              z$frac_unique_leader_010)
  }
  .cmp("attr_budget + attr_tieset vs shap_budget", got, want, tol = 1e-4)

  # --- the verdict ----------------------------------------------------------
  res <- do.call(rbind, GATE_STATE$rows)
  cat("\n=== GATE VERDICT ===\n\n")
  cat(sprintf("  %d checks, %d passed, %d failed\n\n", nrow(res), sum(res$pass),
              sum(!res$pass)))
  if (!all(res$pass)) {
    print(res[!res$pass, ], row.names = FALSE)
    stop("the consolidated library does not reproduce the scattered tables. ",
         "STOP -- plan section 10.5 step 5. Everything downstream of this is ",
         "uninterpretable until it is resolved.", call. = FALSE)
  }
  cat("  PASS. R/14_attribution_eval.R reproduces every legacy table bitwise\n")
  cat("  (to the precision each was saved at) from the identical inputs.\n")
  res
}

if (GATE) {
  res <- run_gate()
  run <- new_run("attrgate", tar_read(cfg), note =
    "correctness gate: R/14_attribution_eval.R against the scattered tables")
  save_table(run, res, "gate_checks", subdir = "diagnostics")
  finalize_run(run, extra = list(ladder = basename(GATE_LADDER),
                                 ties = basename(GATE_TIES),
                                 floor = basename(GATE_FLOOR),
                                 all_pass = all(res$pass)))
  cat(sprintf("\nwritten: %s\n", run$path))
  quit(save = "no", status = 0L)
}

# ============================================================================
# THE CONSUMER
# ============================================================================
#
# REWRITTEN 2026-09-09 AGAINST `docs/attribution_eval_review_20260909.md`. The
# response, finding by finding, is `docs/attribution_eval_review_response_20260909.md`.
# What changed here, keyed to the review:
#
#   A1/A2  the store must carry a FINGERPRINT (design key plus booster
#          settings, design config subset, priors, generator seeds and draw
#          layout) that matches the live design. A store without one is
#          refused with the command that stamps it.
#   A3     the fast top-k path is now the tie-aware definition (library).
#   A4     SHAP's discrimination axis reads the STORED OUT-OF-FOLD MARGIN,
#          never the truncated 19-group row sum; the row-sum AUROC is reported
#          for every method under its own name, `auroc_attr_sum`.
#   A6     the noise tolerance is calibrated PER AGGREGATION LEVEL, and a
#          missing calibration is an `unscored` row, never delta = 0.
#   A7     the magnitude tables evaluate every declared tolerance family --
#          each noise quantile, the common absolute grid, the relative grid --
#          and every declared sign tau; L2 has a calibrated block beside L3
#          and L4; the resolution quantile is a config key that is validated.
#   A8     leader tables carry both populations (all scored, leader differs);
#          share metrics are ORIENTED before dominance and paired deltas.
#   A9     the manifest and every matrix are validated; planned-against-present
#          coverage is reported and a partial store is refused; comparative
#          tables use ONE COMMON BAG POPULATION, declared per row.
#   A10    the default store must be marked complete.
#   A12    leader cells are pooled incrementally and the bag-outer level-4
#          pass streams its per-bag cells to scratch inside the run directory.

ALLOW_INCOMPLETE <- "--allow-incomplete" %in% args
ALLOW_PARTIAL    <- "--allow-partial" %in% args
gen_d <- if (length(pos)) pos[1] else latest_run("attrgen", require_complete = TRUE)
if (is.null(gen_d) || !dir.exists(gen_d)) {
  stop("no COMPLETE attrgen run to consume. Run `Rscript tests/attr_replicates.R` ",
       "first, or pass a run directory (add --allow-incomplete for one whose ",
       "manifest is not marked complete). To check the library itself instead, ",
       "run with --gate.", call. = FALSE)
}
# A GENERATOR RUN THAT NEVER FINISHED IS NOT A STORE (review finding A10).
# `latest_run(require_complete = FALSE)` was the default until 2026-09-09 and
# on this machine it resolved to a 39-file test store written AFTER the
# 642-file one, so a bare `Rscript tests/attr_metrics.R` would have consumed
# the wrong directory without a word.
.gm <- tryCatch(read_manifest(gen_d), error = function(e) NULL)
if (is.null(.gm) || !identical(.gm$status, "complete")) {
  if (!ALLOW_INCOMPLETE) {
    stop("the generator run ", basename(gen_d), " is not marked complete ",
         "(status: ", .gm$status %||% "no manifest", "). Its replicate set may ",
         "be anything. Pass --allow-incomplete to consume it anyway; every ",
         "table is then labelled store_status = incomplete.", call. = FALSE)
  }
  cat("\n  *** consuming an INCOMPLETE generator run (--allow-incomplete)\n")
}
STORE_STATUS <- if (identical(.gm$status, "complete")) "complete" else "incomplete"

cfg     <- tar_read(cfg)
folds   <- tar_read(folds)
tr      <- tar_read(train_ids)
priors  <- tar_read(priors)
domains <- load_domains("config/domains.csv")
sigs    <- as.character(unlist(cfg$signals))
ids_ch  <- as.character(tr)

man  <- utils::read.csv(file.path(gen_d, "manifest_replicates.csv"),
                        stringsAsFactors = FALSE)
meas <- qs2::qs_read(file.path(gen_d, "measured.qs2"))
sub  <- cfg_req(ecfg, "storage", "subdir")

# --- THE STALENESS GUARD ON THE GENERATOR RUN -------------------------------
#
# A GUARD ON THE PRODUCER PROTECTS THE PRODUCER'S OWN REUSE AND NOTHING ELSE,
# which is the lesson `tests/attribution_vs_shap.R` records: the generator
# refuses a mismatched LADDER, but nothing stopped this script from consuming a
# replicate store built under a different design. On 2026-09-06 a table was
# produced whose SHAP half was current and whose L half predated a `gamma`
# change; nothing errored and the only way to find out was to read two manifests
# by hand. The generator stamps `design.qs2`; this recomputes the key from the
# live config and the live fold assignment and STOPS on any differing field.
#
# AND THE FINGERPRINT BESIDE IT (review findings A1, A2). The key covers the
# GAM design only. The fingerprint covers what the key does not -- the booster
# settings, the design subset of config.yml, the priors, the generator's
# seeds, draw count and draw layout -- and the same refusal applies.
cat("\n=== staleness guard ===\n\n")
.dp <- file.path(gen_d, "design.qs2")
if (!file.exists(.dp)) {
  stop("the generator run ", basename(gen_d), " carries no design.qs2, so it ",
       "was written before 2026-09-06 and its design cannot be checked. ",
       "Re-run tests/attr_replicates.R against it -- it is content-addressed ",
       "and resumable, so nothing is refitted.", call. = FALSE)
}
.dg <- qs2::qs_read(.dp)
if (is.null(.dg$k_ti)) {
  stop("the generator run ", basename(gen_d), " carries a design but no `k_ti`. ",
       "A `%||% 5L` here would be a third place the interaction basis dimension ",
       "is written down, silently substituting 5 for whatever the store used. ",
       "Re-run the generator.", call. = FALSE)
}
.fold_vec <- folds$fold[match(tr, folds$stay_id)]
.want <- attr_design_key(cfg, .fold_vec, .dg$k_ti)
.diff <- attr_design_diff(.dg$design_key, .want)
if (nrow(.diff)) {
  cat("\n!!! REPLICATE STORE DESIGN MISMATCH !!!\n\n")
  print(.diff[, c("field", "hash_a", "hash_b", "n_elements_a", "n_elements_b")],
        row.names = FALSE)
  stop("the replicate store in ", basename(gen_d), " was built under a ",
       "different design from the one config/config.yml and the targets store ",
       "now describe (", nrow(.diff), " field(s) differ, listed above). Every ",
       "number below would be from a pipeline that no longer exists. Re-run ",
       "tests/attr_replicates.R.", call. = FALSE)
}
stopifnot(identical(as.character(.dg$stay_id), ids_ch))
cat(sprintf("  replicate store design key %s matches the live design.\n",
            attr_key_hash(.want)))
if (is.null(.dg$fingerprint)) {
  stop("the replicate store in ", basename(gen_d), " carries no FINGERPRINT ",
       "(written before 2026-09-09). Its booster settings, priors, draw seed ",
       "and draw count cannot be checked against the live design. Stamp it, ",
       "with evidence, by resuming the generator on its ladder route:\n",
       "  Rscript tests/attr_replicates.R --levels spec --resume ", gen_d, "\n",
       "That verifies every LLR ladder replicate bitwise against the targets ",
       "cache, refits and verifies the SHAP ladder (about two minutes), and ",
       "writes the stamp beside that evidence.", call. = FALSE)
}
.fp <- attr_design_fingerprint(cfg, .fold_vec, .dg$k_ti, priors, ecfg,
                               design_key = .want)
.fd <- attr_fingerprint_diff(.dg$fingerprint, .fp)
if (nrow(.fd)) {
  cat("\n!!! REPLICATE STORE FINGERPRINT MISMATCH !!!\n\n")
  print(.fd[, c("field", "hash_a", "hash_b", "routes_affected")], row.names = FALSE)
  # ROUTE-SELECTIVE (2026-09-09): refuse only if a replicate the differing
  # field(s) invalidate is actually on disk. The generator quarantines and
  # regenerates those on resume and re-stamps the fingerprint.
  .st <- attr_stale_routes(.fd)
  .stale <- if (.st$all) rep(TRUE, nrow(man)) else attr_stale_mask(man, .st$rules)
  if (any(.stale)) {
    print(table(method = man$method[.stale], route = man$route[.stale]))
    stop("the replicate store's fingerprint differs from the live design in ",
         nrow(.fd), " field(s) (above), and ", sum(.stale), " of ", nrow(man),
         " replicate(s) on disk are on the routes those fields invalidate. The ",
         "design KEY is unchanged, so nothing else would have noticed. Resume ",
         "the generator on the store -- it quarantines and regenerates exactly ",
         "those routes and re-stamps the fingerprint:\n",
         "  Rscript tests/attr_replicates.R --resume ", gen_d, call. = FALSE)
  }
  cat("  the differing field(s) invalidate no replicate on disk; proceeding.\n")
} else {
  cat(sprintf("  replicate store fingerprint %s matches the live design.\n",
              attr_key_hash(.fp)))
}

# THE MANIFEST AND THE MASK (review finding A9): coordinates present, unique
# and integral; every row the same shape; the mask in the store's row order.
attr_validate_manifest(man, "tests/attr_metrics.R")
if (!identical(rownames(meas), ids_ch) || !identical(colnames(meas), sigs)) {
  stop("the store's measured mask is not in stay_id row order or signal ",
       "column order.", call. = FALSE)
}
cat("  manifest carries unique, integral replicate coordinates.\n")

# --- COVERAGE: what the live plan wants against what the store holds -------
#
# ADDED 2026-09-09 (review finding A9). The plan is recomputed from the LIVE
# config, so a count that grew since generation shows as MISSING rather than
# as a quietly smaller distribution. A tombstoned bag is accounted for, not
# missing; anything else absent makes the store PARTIAL, and a partial store
# is refused unless `--allow-partial` says the reader knows.
.plan <- attr_replicate_plan(ecfg, n_folds = as.integer(cfg_req(cfg, "n_folds")),
                             n_pass_fits = 0L)
.tomb_p <- file.path(gen_d, "llr_bootstrap_excluded.csv")
TOMB_BAGS <- integer(0)
if (file.exists(.tomb_p)) {
  .tomb <- utils::read.csv(.tomb_p, stringsAsFactors = FALSE)
  .keys_ok <- c(attr_key_hash(.dg$design_key), .dg$migrated_from)
  TOMB_BAGS <- sort(unique(as.integer(.tomb$boot_id[.tomb$design %in% .keys_ok])))
}
COV <- attr_replicate_coverage(.plan, man, tombstoned = TOMB_BAGS)
cat("\n=== coverage: planned / present / tombstoned / missing ===\n\n")
print(COV, row.names = FALSE)
if (any(COV$missing > 0L)) {
  if (!ALLOW_PARTIAL) {
    stop("the store is PARTIAL: ", sum(COV$missing), " planned replicate(s) are ",
         "neither present nor tombstoned (rows with missing > 0 above). Every ",
         "distribution would be computed on an undeclared subset. Finish the ",
         "generator, or pass --allow-partial to proceed with store_status = ",
         "partial on every table.", call. = FALSE)
  }
  STORE_STATUS <- "partial"
  cat("\n  *** PARTIAL store consumed under --allow-partial\n")
}
if (length(TOMB_BAGS)) {
  cat(sprintf("\n  tombstoned bags (excluded whole, with cause): %s\n",
              paste(TOMB_BAGS, collapse = ", ")))
}

cat(sprintf("\n=== consuming %s (store_status = %s) ===\n\n", basename(gen_d),
            STORE_STATUS))
cat(sprintf("  %d replicates across %d methods and %d route(s)\n",
            nrow(man), length(unique(man$method)), length(unique(man$route))))
print(table(man$method, man$route))

#' Read one replicate and refuse anything but the object every metric assumes.
load_rep <- function(key) {
  M <- qs2::qs_read(file.path(gen_d, sub, paste0(key, ".qs2")))
  attr_check_replicate(M, ids_ch, sigs, key)
  M
}

run <- new_run("attrmetrics", cfg, note = sprintf(
  "attribution reproducibility hierarchy as DISTRIBUTIONS, replicates from %s, fits nothing",
  basename(gen_d)))
save_table(run, COV, "replicate_coverage", subdir = "diagnostics")
save_table(run, ATTR_ESTIMAND_NOTES, "design_notes", subdir = "diagnostics")

DELTA_ABS <- as.numeric(cfg_req(ecfg, "delta", "common_absolute"))
DELTA_REL <- as.numeric(cfg_req(ecfg, "delta", "relative"))
TAUS      <- as.numeric(cfg_req(ecfg, "delta", "sign_taus"))
KS        <- as.integer(cfg_req(ecfg, "aggregation", "top_k"))
LVLS      <- as.character(cfg_req(ecfg, "aggregation", "levels"))
NOISE_RT  <- as.character(cfg_req(ecfg, "delta", "noise_route"))
NOISE_Q   <- as.numeric(cfg_req(ecfg, "delta", "noise_quantiles"))
MAXP      <- as.integer(cfg_req(ecfg, "delta", "max_pairs"))
CFRC      <- as.numeric(cfg_req(ecfg, "delta", "cell_frac"))
DSED      <- as.integer(cfg_req(ecfg, "delta", "sample_seed"))
D_MAXP    <- as.integer(cfg_req(ecfg, "distribution", "max_pairs"))
D_MAXR    <- as.integer(cfg_req(ecfg, "distribution", "max_reps_delta"))
D_SEED    <- as.integer(cfg_req(ecfg, "distribution", "pair_seed"))
COL_SH    <- as.numeric(cfg_req(ecfg, "delta", "collapse_shares"))
COL_MAXP  <- as.integer(cfg_req(ecfg, "distribution", "max_pairs_collapse"))
L3_STRATUM <- as.character(cfg_req(ecfg, "selection", "stability_stratum"))
STAB_C     <- as.character(cfg_req(ecfg, "selection", "stability_contrast"))
RES_Q      <- as.numeric(cfg_req(ecfg, "selection", "resolution_noise_quantile"))

# THE SUPPORTED OUTPUT GRID, VALIDATED BEFORE A MATRIX IS READ (review finding
# A7). Reading a key is not evidence that changing it changes a result; every
# key read above is used below, and the two that other keys depend on are
# checked here so a misconfiguration is a configuration error and not a
# silent `unscored`.
if (!RES_Q %in% NOISE_Q) {
  abort_values(paste0("config: selection.resolution_noise_quantile must be one ",
                      "of delta.noise_quantiles (", paste(NOISE_Q, collapse = ", "),
                      ")"), RES_Q)
}
if (!0 %in% DELTA_ABS) {
  abort_values("config: delta.common_absolute must include 0, the strict number",
               paste(DELTA_ABS, collapse = ","))
}
if (!NOISE_RT %in% names(ATTR_ROUTE_FAMILY)) {
  abort_values("config: delta.noise_route is not a generation route", NOISE_RT)
}

to_level <- function(M, lv) if (lv == "signal") M else attr_to_domain(M, domains)
keep_at  <- function(lv) if (lv == "signal") meas else attr_to_domain(meas * 1, domains) > 0

# `keep` applies to the LLR arms only. A SHAP cell for an unmeasured signal is
# not an assigned zero -- the tree routes a missing feature through its default
# direction and gives it a real contribution -- so masking SHAP would remove
# genuine attributions and flatter our own agreement scores.
keep_for <- function(method, lv) if (attr_method_family(method) == "llr") keep_at(lv) else NULL
keep_pair <- function(a, b, lv) {
  if (attr_method_family(a) == "llr" && attr_method_family(b) == "llr") keep_at(lv) else NULL
}

# THE PRIMARY ROUTE, READ FROM CONFIG AND NOT WRITTEN DOWN HERE. Both families
# bootstrap on the SAME shared bags, which is what makes level 3 a like-for-like
# comparison; one route serves both. The LLR posterior draws stay in their own
# contrast (L3P) and are never the primary quantity.
prim_route <- function(m) NOISE_RT

# ============================================================================
# LEVEL 1: the nats budget, which sets every other scale
# ============================================================================
cat("\n=== L1: the evidence budget ===\n\n")
cat("    Not a defect being measured. The SCALE against which every tolerance\n")
cat("    at every other contrast has to be read: an absolute tolerance of 0.25\n")
cat("    nats forgives twice as much of a SHAP disagreement as of an L one.\n\n")
spec <- man[man$route == "ladder", ]
REP1 <- list()
for (i in seq_len(nrow(spec))) REP1[[spec$method[i]]] <- load_rep(spec$key[i])
if (!length(REP1)) stop("no `ladder` replicates in the store: level 1 and every ",
                        "anchor comparison need the arms as fitted.", call. = FALSE)

# THE METHOD SET IS THE CONFIGURED ONE, NOT WHATEVER HAPPENS TO BE PRESENT
# (review finding A9). `intersect(configured, present)` let a method vanish
# from every comparison without a word.
SCORE <- as.character(cfg_req(ecfg, "methods"))
attr_check_methods(SCORE)
.absent <- setdiff(SCORE, names(REP1))
if (length(.absent)) {
  if (!ALLOW_PARTIAL) {
    abort_values(paste0("configured method(s) have no `ladder` replicate in the ",
                        "store, so they would silently drop out of every ",
                        "comparison. Generate them, or pass --allow-partial"),
                 .absent)
  }
  cat(sprintf("  *** method(s) absent from the store and DROPPED: %s\n",
              paste(.absent, collapse = ", ")))
  SCORE <- setdiff(SCORE, .absent)
  STORE_STATUS <- "partial"
}

cat(sprintf("  %-20s %-8s %10s %10s %10s %10s\n", "method", "agg", "budget",
            "budget_p90", "hhi", "max_share"))
B1 <- list()
for (m in SCORE) for (lv in LVLS) {
  b <- attr_budget(to_level(REP1[[m]], lv))
  B1[[length(B1) + 1L]] <- cbind(method = m, agg = lv, contrast = "L1", b)
  cat(sprintf("  %-20s %-8s %10.4f %10.4f %10.5f %10.5f\n", m, lv,
              b$budget_median, b$budget_p90, b$hhi_median, b$max_share_median))
}
save_table(run, do.call(rbind, B1), "l1_budget", subdir = "diagnostics")

cat("\n--- resolution: how many contributions are in contention at all ---\n\n")
cat("    A property of ONE arm and the ceiling on what any per-patient\n")
cat("    explanation can claim. This is what catches degenerate stability.\n\n")
R1 <- list()
for (m in SCORE) for (lv in LVLS) {
  R1[[length(R1) + 1L]] <- cbind(method = m, agg = lv,
                                 attr_tieset(attr_prep(to_level(REP1[[m]], lv)), DELTA_ABS))
}
R1 <- do.call(rbind, R1)
save_table(run, R1, "l1_resolution", subdir = "diagnostics")

# ============================================================================
# THE COMMON BAG POPULATION
# ============================================================================
#
# ADDED 2026-09-09 (review finding A9). An LLR bag is excluded whole when any
# spec-fold cannot be fitted; a SHAP bag never is. So SHAP's level 3 rested
# on 40 bags and every LLR arm's on 38, and any comparison between them --
# the selection table, the dominance table -- compared stability under
# different success criteria. Every COMPARATIVE quantity below is computed
# on the bags every scored method succeeded on, and the row says so in
# `bag_population = common`. A method's own full population is reported
# beside it as `own` wherever the two differ, so nothing is thrown away.
boot_man <- man[man$route == "bootstrap", ]
BAGS_BY_M <- lapply(stats::setNames(SCORE, SCORE), function(m)
  sort(unique(boot_man$boot_id[boot_man$method == m])))
COMMON_BAGS <- if (length(BAGS_BY_M)) Reduce(intersect, BAGS_BY_M) else integer(0)
in_common <- function(b) b == 0L | b %in% COMMON_BAGS
.planned_bags <- vapply(SCORE, function(m)
  sum(.plan$method == m & .plan$route == "bootstrap"), integer(1))
BAGPOP <- data.frame(
  method = SCORE, n_bags_planned = .planned_bags,
  n_bags_succeeded = lengths(BAGS_BY_M)[SCORE],
  n_bags_tombstoned = vapply(SCORE, function(m)
    if (attr_method_family(m) == "llr") length(TOMB_BAGS) else 0L, integer(1)),
  n_bags_common = length(COMMON_BAGS),
  common_bags = paste(COMMON_BAGS, collapse = ","),
  store_status = STORE_STATUS, stringsAsFactors = FALSE)
save_table(run, BAGPOP, "bag_population", subdir = "diagnostics")
cat(sprintf("\n  common bag population: %d bags shared by every scored method\n",
            length(COMMON_BAGS)))
for (m in SCORE) cat(sprintf("    %-20s %2d succeeded of %2d planned\n", m,
                             length(BAGS_BY_M[[m]]), .planned_bags[[m]]))

# ============================================================================
# THE NOISE-CALIBRATED DELTA, PER METHOD, ROUTE AND AGGREGATION LEVEL
# ============================================================================
#
# Computed first, because every noise-calibrated table below is evaluated at
# it. Grouped by (method, route) and NEVER pooled across routes: a posterior
# draw and a bootstrap refit estimate different quantities.
#
# PER LEVEL AS OF 2026-09-09 (review finding A6). The displacement was
# measured on signal matrices only and the resulting scalar was applied to
# domain SUMS under the same `noise_calibrated` label. Summing signed
# contributions can cancel or amplify noise, so the domain-level displacement
# is not the signal-level one -- two replicates (1, 1) and (2, 0) move every
# signal cell by 1 and their one-domain sum by 0. The tolerance is now
# measured on the matrices it will be applied to. A method with no
# calibration on a level gets NA here and an `unscored` row below, never a
# zero.
cat("\n=== the noise-calibrated tolerance, per method, route and level ===\n\n")
DELTA_NOISE <- stats::setNames(vector("list", length(LVLS)), LVLS)
NOISE <- list()
samp <- man[man$route %in% c("seed", "bootstrap", "bootstrap_seeded", "posterior") &
            man$method %in% SCORE, ]
samp <- samp[in_common(samp$boot_id), ]
for (g in unique(paste(samp$method, samp$route, sep = "|"))) {
  z <- samp[paste(samp$method, samp$route, sep = "|") == g, ]
  m <- z$method[1]; rt <- z$route[1]
  if (nrow(z) < 2L) next
  mats <- lapply(z$key, load_rep)
  is_prim <- identical(rt, prim_route(m))
  for (lv in LVLS) {
    ml   <- lapply(mats, to_level, lv = lv)
    disp <- attr_displacement(ml, keep = keep_for(m, lv),
                              max_pairs = MAXP, cell_frac = CFRC, seed = DSED)
    dp <- attr_delta_policy(disp, NOISE_Q)
    NOISE[[length(NOISE) + 1L]] <- cbind(method = m, route = rt, agg = lv,
                                         b = nrow(z), calibrates = is_prim, dp)
    if (is_prim) DELTA_NOISE[[lv]][[m]] <- stats::setNames(dp$delta, sprintf("q%g", dp$q))
    cat(sprintf("  %-20s %-18s %-7s B=%2d  q50 = %.4f  q90 = %.4f nats%s\n", m, rt, lv,
                nrow(z), dp$delta[dp$q == 0.5][1], dp$delta[dp$q == 0.9][1],
                if (is_prim) "   <- calibrates delta" else ""))
    rm(ml)
  }
  rm(mats); gc(verbose = FALSE)
}
if (length(NOISE)) save_table(run, do.call(rbind, NOISE), "noise_calibrated_delta",
                              subdir = "diagnostics")

#' The calibrated tolerance for one method at one level and quantile, or NA.
noise_delta <- function(m, lv, q) {
  v <- DELTA_NOISE[[lv]][[m]]
  if (is.null(v)) return(NA_real_)
  d <- v[sprintf("q%g", q)]
  if (is.null(d) || is.na(d)) NA_real_ else unname(d)
}

# ============================================================================
# WITHIN-METHOD CONTRASTS, AS DISTRIBUTIONS
# ============================================================================
#
# THE FAST PATH AND WHY IT MATTERS. Top-k agreement at delta = 0 is two
# bitwise tests per patient once each replicate's top-k set and boundary-tie
# set are encoded (`attr_topk()`), so EVERY available pair can be used and the
# cost is O(B) encodings rather than O(B^2) matrix passes. The gate asserts
# the encoding reproduces `attr_agree_k(delta = 0)` exactly, boundary ties
# included.
CP <- attr_contrast_pairs(man[man$method %in% SCORE, ])
CP$common <- in_common(CP$boot_i) & in_common(CP$boot_j)
if (!nrow(CP)) {
  cat("\n  no within-method pairs in this store. Generate stage `sample`.\n")
}
cat(sprintf("\n=== within-method contrasts: %d pairs across %d (method, family, contrast, stratum) cells ===\n\n",
            nrow(CP), length(unique(paste(CP$method, CP$route_family, CP$code,
                                          CP$held)))))
for (cd in sort(unique(CP$code))) {
  cat(sprintf("  %-4s %s\n", cd, ATTR_CONTRAST_LABELS[[cd]]))
}

#' The pairs of one cell, in the requested population, bounded and seeded.
cell_pairs <- function(z, pop, bound = D_MAXP) {
  if (pop == "common") z <- z[z$common, , drop = FALSE]
  if (nrow(z) > bound) z <- z[with_seed(D_SEED, sort(sample.int(nrow(z), bound))), , drop = FALSE]
  z
}
#' Which populations a cell is reported in: `common` always, `own` only where
#' it differs from `common`.
cell_pops <- function(z) if (all(z$common)) "common" else c("common", "own")

DIST <- list(); RAW <- list()
for (lv in LVLS) {
  cat(sprintf("\n--- aggregation level `%s` ---\n\n", lv))
  # One pass over the store per aggregation level. Encodings are two integer
  # vectors per k, so the whole store's encodings are a few hundred megabytes
  # and the matrices themselves are never all held.
  TK <- list()
  for (i in seq_len(nrow(man))) {
    if (!man$method[i] %in% SCORE) next
    TK[[man$key[i]]] <- attr_topk(to_level(load_rep(man$key[i]), lv), KS)
  }
  # A CELL IS (method, route family, contrast, STRATUM). Once `bootstrap` and
  # `bootstrap_seeded` pair with each other, one method has TWO L3
  # distributions -- one at each held seed -- and two L2 distributions, one on
  # the original sample and one on the resamples. Keying without `held` would
  # pool them.
  cells <- unique(CP[, c("method", "route_family", "code", "held")])
  for (ci in seq_len(nrow(cells))) {
    m <- cells$method[ci]; rt <- cells$route_family[ci]
    cd <- cells$code[ci]; hd <- cells$held[ci]
    z_all <- CP[CP$method == m & CP$route_family == rt & CP$code == cd &
                CP$held == hd, , drop = FALSE]
    for (pop in cell_pops(z_all)) {
      z <- cell_pairs(z_all, pop)
      if (!nrow(z)) next
      for (k in KS) {
        v <- vapply(seq_len(nrow(z)), function(q)
          1 - attr_topk_agree(TK[[man$key[z$i[q]]]], TK[[man$key[z$j[q]]]], k),
          numeric(1))
        DIST[[length(DIST) + 1L]] <- cbind(
          data.frame(contrast = cd, hierarchy_level = z$hierarchy_level[1],
                     method_a = m, method_b = m, route = rt, held = hd,
                     agg = lv, k = k, delta_kind = "absolute_nats", delta = 0,
                     bag_population = pop, stringsAsFactors = FALSE),
          attr_stat_summary(v))
        # The per-pair values are kept for the dominance and paired-delta tables.
        RAW[[length(RAW) + 1L]] <- data.frame(
          contrast = cd, method_a = m, method_b = m, route = rt, held = hd,
          agg = lv, metric = sprintf("top%d", k), k = k, bag_population = pop,
          boot_i = z$boot_i, boot_j = z$boot_j, value = v,
          stringsAsFactors = FALSE)
      }
    }
  }
  rm(TK); gc(verbose = FALSE)
}
DIST_W <- if (length(DIST)) do.call(rbind, DIST) else NULL
RAW_W  <- if (length(RAW))  do.call(rbind, RAW)  else NULL
if (!is.null(DIST_W)) {
  DIST_W$anchor_boot0 <- NA_real_
  cat("\n  strict top-1 disagreement (1 - agreement, delta = 0), signal level, common bags\n\n")
  cat(sprintf("  %-4s %-20s %-9s %-16s %5s %7s %7s %7s %7s %7s %7s\n",
              "code", "method", "family", "held", "pairs", "min", "median",
              "mean", "p90", "p95", "max"))
  z <- DIST_W[DIST_W$agg == "signal" & DIST_W$k == 1L & DIST_W$bag_population == "common", ]
  z <- z[order(z$contrast, z$median), ]
  for (i in seq_len(nrow(z))) cat(sprintf(
    "  %-4s %-20s %-9s %-16s %5d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
    z$contrast[i], z$method_a[i], z$route[i], z$held[i], z$n[i], z$min[i],
    z$median[i], z$mean[i], z$p90[i], z$p95[i], z$max[i]))
  cat("\n    `mean` above `median` is a right tail: the median understates how\n")
  cat("    bad a bad resample is. `p95` is the criterion; `max` is a diagnostic\n")
  cat("    that grows mechanically with the number of pairs.\n")
}

# ============================================================================
# LEVEL 4: specification, as a DISTRIBUTION over the shared bags
# ============================================================================
#
# Arm A and arm B are fitted on the IDENTICAL bag b, so D(A_b, B_b) holds the
# sample fixed and the variation caused by which patients entered bag b is
# shared and differences out. Across the COMMON bags that is a distribution;
# the anchor (`boot_id = 0`, the two arms as fitted on all of train) is
# reported beside it and is one draw from the distribution, not its summary.
cat("\n=== L4: specification, as a distribution over the common bags ===\n\n")
CMP4 <- utils::combn(SCORE, 2L, simplify = FALSE)

L4 <- list(); L4RAW <- list()
for (lv in LVLS) {
  use <- rbind(spec[spec$method %in% SCORE, names(boot_man)],
               boot_man[boot_man$method %in% SCORE & in_common(boot_man$boot_id), ])
  TK <- list()
  for (i in seq_len(nrow(use))) {
    TK[[use$key[i]]] <- attr_topk(to_level(load_rep(use$key[i]), lv), KS)
  }
  for (p in CMP4) {
    ka <- use[use$method == p[1], ]; kb <- use[use$method == p[2], ]
    # PAIRED ON `boot_id`. `attr_pair_contrast()` refuses any pair that differs
    # in a coordinate as well as in the method.
    shared <- intersect(ka$boot_id, kb$boot_id)
    if (!length(shared)) next
    for (k in KS) {
      v <- vapply(shared, function(b) {
        A <- ka$key[ka$boot_id == b][1]; B <- kb$key[kb$boot_id == b][1]
        co <- attr_pair_contrast(
          list(method = p[1], boot_id = b, seed_id = 0L, draw_id = 0L),
          list(method = p[2], boot_id = b, seed_id = 0L, draw_id = 0L))
        stopifnot(identical(co$code, "L4"))
        1 - attr_topk_agree(TK[[A]], TK[[B]], k)
      }, numeric(1))
      anchor <- if (0L %in% shared) v[match(0L, shared)] else NA_real_
      resamp <- v[shared != 0L]
      L4[[length(L4) + 1L]] <- cbind(
        data.frame(contrast = "L4", hierarchy_level = 4L, method_a = p[1],
                   method_b = p[2], route = "refit", held = "boot,seed,draw",
                   agg = lv, k = k, delta_kind = "absolute_nats", delta = 0,
                   bag_population = "common",
                   anchor_boot0 = round(anchor, 6), stringsAsFactors = FALSE),
        attr_stat_summary(resamp))
      if (length(resamp)) {
        bb <- shared[shared != 0L]
        L4RAW[[length(L4RAW) + 1L]] <- data.frame(
          contrast = "L4", method_a = p[1], method_b = p[2], route = "refit",
          held = "boot,seed,draw", agg = lv, metric = sprintf("top%d", k), k = k,
          bag_population = "common",
          boot_i = bb, boot_j = bb, value = resamp, stringsAsFactors = FALSE)
      }
    }
  }
  rm(TK); gc(verbose = FALSE)
}
L4_W    <- if (length(L4)) do.call(rbind, L4) else NULL
L4RAW_W <- if (length(L4RAW)) do.call(rbind, L4RAW) else NULL

ALLDIST <- if (is.null(DIST_W)) L4_W else
  if (is.null(L4_W)) DIST_W else rbind(DIST_W, L4_W)
if (!is.null(ALLDIST)) {
  ALLDIST$store_status <- STORE_STATUS
  save_table(run, ALLDIST, "disagreement_distributions", subdir = "diagnostics")
}

cat(sprintf("  %-20s %-20s %-7s %2s %8s %8s %8s %8s %8s %8s\n", "method a",
            "method b", "agg", "k", "anchor", "median", "mean", "p90", "p95", "max"))
if (!is.null(L4_W)) {
  z <- L4_W[L4_W$agg == "signal" & L4_W$k == 1L, ]
  z <- z[order(z$median), ]
  for (i in seq_len(nrow(z))) cat(sprintf(
    "  %-20s %-20s %-7s %2d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
    z$method_a[i], z$method_b[i], z$agg[i], z$k[i], z$anchor_boot0[i],
    z$median[i], z$mean[i], z$p90[i], z$p95[i], z$max[i]))
}

# ============================================================================
# COSINE DISSIMILARITY AND LEADER COLLAPSE, PER PAIR, OVER THE FULL PAIR SET
# ============================================================================
#
# Cosine needs no ordering -- three vector passes per pair, one once norms are
# cached -- so it runs over every pair the rank metrics use; the leader
# collapse is `max.col`, `rowSums` and one comparison and is bounded
# separately at `max_pairs_collapse`. The value stored is `1 - cosine`, a
# DISSIMILARITY, so larger means more disagreement like every other metric.
#
# LEADER CELLS ARE POOLED INCREMENTALLY AND SUMMARISED PER CELL (review
# finding A12). The first version kept every pair's full 41,250-row cells for
# every cell until the end of the run. Now each pair's scored entries are
# appended to the cell's pool and the pool is summarised and dropped as soon
# as its cell is complete. The level-4 loop is bag-outer for the reason
# section 32 of the plan gives (each (method, bag) matrix is read once), so
# its per-bag cells are streamed to scratch INSIDE the run directory and
# pooled per pair at the end; the scratch is deleted before the run closes.
cat("\n=== cosine dissimilarity and leader collapse, per pair ===\n\n")
COSRAW <- list(); LCD <- list()
LC_NAMES <- names(attr_leader_pair_scalars(attr_leader_cells(
  matrix(c(1, 2, 3, 4), 2), matrix(c(4, 3, 2, 1), 2)), shares = COL_SH))
TMP <- file.path(run$path, "tmp_leader_cells")
dir.create(TMP, showWarnings = FALSE)
emit_lcd <- function(hdr, pool) {
  LCD[[length(LCD) + 1L]] <<- cbind(site = "mimic", hdr,
                                    attr_leader_distribution(pool, COL_SH))
}
for (lv in LVLS) {
  # --- within method: every contrast the store supports, common bags --------
  cells <- unique(CP[, c("method", "route_family", "code", "held")])
  for (ci in seq_len(nrow(cells))) {
    m <- cells$method[ci]; rt <- cells$route_family[ci]
    cd <- cells$code[ci]; hd <- cells$held[ci]
    z <- cell_pairs(CP[CP$method == m & CP$route_family == rt & CP$code == cd &
                       CP$held == hd, , drop = FALSE], "common")
    if (!nrow(z)) next
    keys <- unique(c(man$key[z$i], man$key[z$j]))
    M <- lapply(stats::setNames(keys, keys), function(k) to_level(load_rep(k), lv))
    N <- lapply(M, attr_row_norm)
    kp <- keep_for(m, lv)
    v <- vapply(seq_len(nrow(z)), function(q) {
      ka <- man$key[z$i[q]]; kb <- man$key[z$j[q]]
      stats::median(attr_cosine_dissim(M[[ka]], M[[kb]], N[[ka]], N[[kb]]),
                    na.rm = TRUE)
    }, numeric(1))
    COSRAW[[length(COSRAW) + 1L]] <- data.frame(
      contrast = cd, method_a = m, method_b = m, route = rt, held = hd,
      agg = lv, metric = "cosine_dissim", k = NA_integer_, bag_population = "common",
      boot_i = z$boot_i, boot_j = z$boot_j, value = v, stringsAsFactors = FALSE)
    zc <- if (nrow(z) <= COL_MAXP) z else
      z[with_seed(D_SEED, sort(sample.int(nrow(z), COL_MAXP))), , drop = FALSE]
    pool <- .pool_leader_new()
    cvm <- vapply(seq_len(nrow(zc)), function(q) {
      ka <- man$key[zc$i[q]]; kb <- man$key[zc$j[q]]
      cl <- attr_leader_cells(M[[ka]], M[[kb]], keep = kp)
      pool <<- .pool_leader_add(pool, cl)
      attr_leader_pair_scalars(cl, shares = COL_SH)
    }, numeric(length(LC_NAMES)))
    emit_lcd(data.frame(contrast = cd, held = hd, method_a = m, method_b = m,
                        agg = lv, route = rt, n_pairs = nrow(zc),
                        bag_population = "common", masked = !is.null(kp),
                        stringsAsFactors = FALSE), pool)
    rm(pool)
    for (mt in LC_NAMES) {
      COSRAW[[length(COSRAW) + 1L]] <- data.frame(
        contrast = cd, method_a = m, method_b = m, route = rt, held = hd,
        agg = lv, metric = mt, k = NA_integer_, bag_population = "common",
        boot_i = zc$boot_i, boot_j = zc$boot_j,
        value = cvm[mt, ], stringsAsFactors = FALSE)
    }
    rm(M, N); gc(verbose = FALSE)
  }
  # --- level 4, paired on the bag, bag-outer, streamed --------------------
  #
  # BOUNDED BY `distribution.max_pairs_collapse` LIKE THE WITHIN-METHOD PATH
  # (re-review, A7): the bag set is one observation per bag, so the bound is
  # on bags, drawn under `pair_seed` when it binds. At 38 common bags against
  # a bound of 150 it does not bind today.
  acc <- list(); cc <- list()
  bags4 <- if (length(COMMON_BAGS) <= COL_MAXP) COMMON_BAGS else
    sort(with_seed(D_SEED, sample(COMMON_BAGS, COL_MAXP)))
  for (b in bags4) {
    zb <- boot_man[boot_man$boot_id == b & boot_man$method %in% SCORE, , drop = FALSE]
    M <- lapply(stats::setNames(zb$method, zb$method),
                function(mm) to_level(load_rep(zb$key[zb$method == mm][1]), lv))
    N <- lapply(M, attr_row_norm)
    for (pp in CMP4) {
      if (is.null(M[[pp[1]]]) || is.null(M[[pp[2]]])) next
      tag <- paste(pp[1], pp[2], sep = "|")
      acc[[tag]] <- c(acc[[tag]], stats::median(
        attr_cosine_dissim(M[[pp[1]]], M[[pp[2]]], N[[pp[1]]], N[[pp[2]]]),
        na.rm = TRUE))
      names(acc[[tag]])[length(acc[[tag]])] <- as.character(b)
      cl4 <- attr_leader_cells(M[[pp[1]]], M[[pp[2]]], keep = keep_pair(pp[1], pp[2], lv))
      cc[[tag]] <- rbind(cc[[tag]], c(bag = b, attr_leader_pair_scalars(cl4, shares = COL_SH)))
      qs2::qs_save(.pool_leader_extract(cl4),
                   file.path(TMP, sprintf("%s__%s__%02d.qs2", lv, gsub("|", "__", tag, fixed = TRUE), b)))
      rm(cl4)
    }
    rm(M, N); gc(verbose = FALSE)
  }
  for (tag in names(acc)) {
    pp <- strsplit(tag, "|", fixed = TRUE)[[1]]
    fs <- file.path(TMP, sprintf("%s__%s__%02d.qs2", lv, gsub("|", "__", tag, fixed = TRUE),
                                 as.integer(names(acc[[tag]]))))
    pool <- .pool_leader_new()
    for (f in fs) pool <- .pool_leader_add(pool, qs2::qs_read(f))
    unlink(fs)
    emit_lcd(data.frame(contrast = "L4", held = "shared_bags", method_a = pp[1],
                        method_b = pp[2], agg = lv, route = "refit",
                        n_pairs = length(fs), bag_population = "common",
                        masked = !is.null(keep_pair(pp[1], pp[2], lv)),
                        stringsAsFactors = FALSE), pool)
    rm(pool)
    bb <- as.integer(names(acc[[tag]]))
    COSRAW[[length(COSRAW) + 1L]] <- data.frame(
      contrast = "L4", method_a = pp[1], method_b = pp[2], route = "refit",
      held = "boot,seed,draw", agg = lv, metric = "cosine_dissim",
      k = NA_integer_, bag_population = "common",
      boot_i = bb, boot_j = bb, value = unname(acc[[tag]]), stringsAsFactors = FALSE)
    if (!is.null(cc[[tag]])) {
      mm <- cc[[tag]]; bg <- as.integer(mm[, "bag"])
      for (mt in LC_NAMES) {
        COSRAW[[length(COSRAW) + 1L]] <- data.frame(
          contrast = "L4", method_a = pp[1], method_b = pp[2], route = "refit",
          held = "boot,seed,draw", agg = lv, metric = mt,
          k = NA_integer_, bag_population = "common", boot_i = bg, boot_j = bg,
          value = unname(mm[, mt]), stringsAsFactors = FALSE)
      }
    }
  }
}
unlink(TMP, recursive = TRUE)
COS_W <- if (length(COSRAW)) do.call(rbind, COSRAW) else NULL
if (!is.null(COS_W)) {
  # `metric` IS PART OF THE SPLIT KEY: cosine dissimilarity and the leader
  # scalars share this table and are on different scales.
  CD <- do.call(rbind, lapply(split(COS_W, paste(COS_W$metric, COS_W$contrast,
                                                 COS_W$method_a, COS_W$method_b,
                                                 COS_W$agg, COS_W$held)),
    function(g) {
      stopifnot(length(unique(g$metric)) == 1L)
      cbind(g[1, c("contrast", "method_a", "method_b", "route", "held",
                   "agg", "metric", "bag_population")],
            data.frame(orientation = attr_metric_orientation(g$metric[1]),
                       stringsAsFactors = FALSE),
            attr_stat_summary(g$value))
    }))
  CD$store_status <- STORE_STATUS
  save_table(run, CD, "pairwise_metric_distributions", subdir = "diagnostics")
  z <- CD[CD$agg == "signal" & CD$contrast %in% c("L3", "L4") & CD$metric == "cosine_dissim", ]
  z <- z[order(z$contrast, z$median), ]
  cat(sprintf("  %-4s %-20s %-20s %5s %9s %9s %9s\n", "code", "method a",
              "method b", "pairs", "median", "p95", "max"))
  for (i in seq_len(nrow(z))) cat(sprintf("  %-4s %-20s %-20s %5d %9.5f %9.5f %9.5f\n",
    z$contrast[i], z$method_a[i], z$method_b[i], z$n[i], z$median[i],
    z$p95[i], z$max[i]))
  cat("\n    Values are 1 - cosine, so LARGER IS MORE DISAGREEMENT, matching the\n")
  cat("    top-k columns. Subtract from 1 to read as a similarity.\n")
}

# ============================================================================
# THE TOLERANCE GRID, COSINE AND SIGN-FLIP DISTRIBUTIONS
# ============================================================================
#
# Bounded by replicates rather than by pairs, because these need the full
# magnitude matrices. EVERY DECLARED TOLERANCE FAMILY IS EVALUATED (review
# finding A7): the noise-calibrated tolerance at each declared quantile, the
# common absolute grid, and the relative grid on each arm's own total; and
# every declared sign tau. `attr_agree_k_grid()` makes the grid one pass per
# (pair, k). A tolerance that could not be calibrated for a method and level
# produces a row with `status = unscored_no_calibration` and no numbers.
#
# THREE BLOCKS: L3 within method on the common bags (the primary stability
# quantity), L2 within method on the original sample (SHAP only, by
# construction), and L4 paired on the common bags.
cat(sprintf("\n=== tolerance grid, cosine and sign flip, over %d replicates per cell ===\n\n",
            D_MAXR))
delta_grid <- function(methods, lv) {
  rows <- list()
  for (q in NOISE_Q) {
    d <- vapply(methods, function(m) noise_delta(m, lv, q), numeric(1))
    rows[[length(rows) + 1L]] <- data.frame(
      delta_kind = sprintf("noise_calibrated_q%g", q),
      delta = if (anyNA(d)) NA_real_ else max(d), relative = FALSE,
      stringsAsFactors = FALSE)
  }
  for (d in DELTA_ABS) rows[[length(rows) + 1L]] <- data.frame(
    delta_kind = "common_absolute", delta = d, relative = FALSE, stringsAsFactors = FALSE)
  for (d in DELTA_REL) rows[[length(rows) + 1L]] <- data.frame(
    delta_kind = "relative_own_total", delta = d, relative = TRUE, stringsAsFactors = FALSE)
  do.call(rbind, rows)
}
#' Agreement over the grid for one (pair, k): NA where the grid row is unscored.
grid_agree <- function(pa, pb, k, G) {
  out <- rep(NA_real_, nrow(G))
  a_ok <- !G$relative & !is.na(G$delta); r_ok <- G$relative
  if (any(a_ok)) out[a_ok] <- attr_agree_k_grid(pa, pb, k, G$delta[a_ok])
  if (any(r_ok)) out[r_ok] <- attr_agree_k_grid(pa, pb, k, G$delta[r_ok], relative = TRUE)
  out
}
#' Emit the grid rows, the cosine row and the sign-flip rows for one cell
#' from per-pair matrices `AG[[k]]` (pairs x grid), `cs` and `FL` (pairs x taus).
emit_aux <- function(hdr, G, AG, cs, FL) {
  rows <- list()
  for (k in KS) for (g in seq_len(nrow(G))) {
    v <- if (is.na(G$delta[g])) numeric(0) else 1 - AG[[as.character(k)]][, g]
    rows[[length(rows) + 1L]] <- cbind(hdr, data.frame(
      k = k, delta_kind = G$delta_kind[g], delta = round(G$delta[g], 6),
      metric = "disagree",
      status = if (is.na(G$delta[g])) "unscored_no_calibration" else "scored",
      stringsAsFactors = FALSE), attr_stat_summary(v))
  }
  rows[[length(rows) + 1L]] <- cbind(hdr, data.frame(
    k = NA_integer_, delta_kind = "none", delta = NA_real_, metric = "cosine_median",
    status = "scored", stringsAsFactors = FALSE), attr_stat_summary(cs))
  for (t in seq_along(TAUS)) rows[[length(rows) + 1L]] <- cbind(hdr, data.frame(
    k = NA_integer_, delta_kind = "sign_tau", delta = TAUS[t], metric = "sign_flip",
    status = "scored", stringsAsFactors = FALSE), attr_stat_summary(FL[, t]))
  do.call(rbind, rows)
}
#' One within-method cell: replicates `z` (manifest rows, in order).
aux_within <- function(m, z, contrast, route_label, held, lv) {
  z <- z[seq_len(min(nrow(z), D_MAXR)), , drop = FALSE]
  M <- lapply(z$key, function(k) to_level(load_rep(k), lv))
  pr <- lapply(M, attr_prep)
  pix <- attr_pair_index(length(M), D_MAXP, D_SEED)
  G <- delta_grid(m, lv); kp <- keep_for(m, lv)
  AG <- stats::setNames(lapply(KS, function(k) t(vapply(seq_len(ncol(pix)), function(q)
    grid_agree(pr[[pix[1, q]]], pr[[pix[2, q]]], k, G), numeric(nrow(G))))),
    as.character(KS))
  cs <- vapply(seq_len(ncol(pix)), function(q)
    stats::median(attr_cosine(M[[pix[1, q]]], M[[pix[2, q]]]), na.rm = TRUE), numeric(1))
  FL <- t(vapply(seq_len(ncol(pix)), function(q)
    attr_sign_flip(M[[pix[1, q]]], M[[pix[2, q]]], keep = kp, taus = TAUS)$flip,
    numeric(length(TAUS))))
  rm(M, pr); gc(verbose = FALSE)
  emit_aux(data.frame(contrast = contrast, method_a = m, method_b = m,
                      route = route_label, held = held, agg = lv,
                      bag_population = "common", n_replicates = nrow(z),
                      stringsAsFactors = FALSE), G, AG, cs, FL)
}
AUXD <- list()
for (lv in LVLS) {
  # --- L3, within method, on the common bags, at the primary stratum -------
  for (m in SCORE) {
    z <- boot_man[boot_man$method == m & boot_man$boot_id > 0L &
                  in_common(boot_man$boot_id) & boot_man$seed_id == 0L, , drop = FALSE]
    z <- z[order(z$boot_id), , drop = FALSE]
    if (nrow(z) < 2L) next
    AUXD[[length(AUXD) + 1L]] <- aux_within(m, z, "L3", "bootstrap", L3_STRATUM, lv)
  }
  # --- L2, within method, on the original sample (SHAP only by construction) -
  for (m in SCORE) {
    if (!attr_has_seed_noise(m)) next
    z <- man[man$method == m & man$route %in% c("ladder", "seed"), , drop = FALSE]
    z <- z[order(z$seed_id), , drop = FALSE]
    if (nrow(z) < 2L) next
    AUXD[[length(AUXD) + 1L]] <- aux_within(m, z, "L2", "seed", "original_sample", lv)
  }
  # --- L4, paired on the bag, bag-outer ------------------------------------
  bags <- COMMON_BAGS[seq_len(min(length(COMMON_BAGS), D_MAXR))]
  acc <- list()
  for (b in bags) {
    z <- boot_man[boot_man$boot_id == b & boot_man$method %in% SCORE, , drop = FALSE]
    MM <- PP <- list()
    for (i in seq_len(nrow(z))) {
      MM[[z$method[i]]] <- to_level(load_rep(z$key[i]), lv)
      PP[[z$method[i]]] <- attr_prep(MM[[z$method[i]]])
    }
    for (p in CMP4) {
      if (is.null(MM[[p[1]]]) || is.null(MM[[p[2]]])) next
      tag <- paste(p[1], p[2], sep = "|")
      if (is.null(acc[[tag]])) {
        acc[[tag]] <- list(G = delta_grid(p, lv), cs = numeric(0), FL = NULL,
                           AG = stats::setNames(rep(list(NULL), length(KS)), as.character(KS)))
      }
      kp <- keep_pair(p[1], p[2], lv)
      for (k in KS) {
        acc[[tag]]$AG[[as.character(k)]] <- rbind(acc[[tag]]$AG[[as.character(k)]],
          grid_agree(PP[[p[1]]], PP[[p[2]]], k, acc[[tag]]$G))
      }
      acc[[tag]]$cs <- c(acc[[tag]]$cs,
        stats::median(attr_cosine(MM[[p[1]]], MM[[p[2]]]), na.rm = TRUE))
      acc[[tag]]$FL <- rbind(acc[[tag]]$FL,
        attr_sign_flip(MM[[p[1]]], MM[[p[2]]], keep = kp, taus = TAUS)$flip)
    }
    rm(MM, PP); gc(verbose = FALSE)
  }
  for (tag in names(acc)) {
    p <- strsplit(tag, "|", fixed = TRUE)[[1]]; a <- acc[[tag]]
    AUXD[[length(AUXD) + 1L]] <- emit_aux(
      data.frame(contrast = "L4", method_a = p[1], method_b = p[2],
                 route = "bootstrap", held = "boot,seed,draw", agg = lv,
                 bag_population = "common", n_replicates = length(a$cs),
                 stringsAsFactors = FALSE), a$G, a$AG, a$cs, a$FL)
  }
}
if (length(AUXD)) {
  AUXD <- do.call(rbind, AUXD); AUXD$store_status <- STORE_STATUS
  save_table(run, AUXD, "magnitude_metric_distributions", subdir = "diagnostics")
  z <- AUXD[AUXD$agg == "signal" & AUXD$metric == "disagree" & AUXD$k == 1L &
            AUXD$delta_kind == sprintf("noise_calibrated_q%g", RES_Q), ]
  z <- z[order(z$contrast, z$median), ]
  cat(sprintf("  top-1 disagreement at the noise-calibrated q%g tolerance, signal level\n\n", RES_Q))
  cat(sprintf("  %-4s %-20s %-20s %8s %-9s %7s %7s\n", "code", "method a", "method b",
              "delta", "status", "median", "p95"))
  for (i in seq_len(nrow(z))) cat(sprintf("  %-4s %-20s %-20s %8.4f %-9s %7.4f %7.4f\n",
    z$contrast[i], z$method_a[i], z$method_b[i], z$delta[i],
    if (z$status[i] == "scored") "scored" else "UNSCORED", z$median[i], z$p95[i]))
}

# ============================================================================
# LEADER COLLAPSE: the pooled cells above, plus the anchor fits
# ============================================================================
#
# Every row is a distribution over patients (or patient-pair observations, per
# `unit`) with the fixed `attr_stat_summary()` columns, in TWO POPULATIONS:
# every scored patient, and the patients whose leader differs between the two
# arms. `tests/attr_external.R` writes the same schema with `site = "eicu"`.
cat("\n=== leader collapse: where does one arm's leader land in the other? ===\n\n")
for (lv in LVLS) for (p in CMP4) {
  if (is.null(REP1[[p[1]]]) || is.null(REP1[[p[2]]])) next
  kp <- keep_pair(p[1], p[2], lv)
  LCD[[length(LCD) + 1L]] <- cbind(
    site = "mimic",
    data.frame(contrast = "L4", held = "boot=0", method_a = p[1], method_b = p[2],
               agg = lv, route = "refit", n_pairs = 1L, bag_population = "anchor",
               masked = !is.null(kp), stringsAsFactors = FALSE),
    attr_leader_patient_summary(to_level(REP1[[p[1]]], lv), to_level(REP1[[p[2]]], lv),
                                keep = kp, shares = COL_SH))
}
if (length(LCD)) {
  LCD <- do.call(rbind, LCD); LCD$store_status <- STORE_STATUS
  save_table(run, LCD, "leader_collapse_distributions", subdir = "diagnostics")
  z <- LCD[LCD$agg == "signal" & LCD$metric == "rank_displacement" &
           LCD$direction == "a_leader_in_b" & LCD$population == "leader_differs" &
           (LCD$held == "boot=0" | (LCD$contrast == "L3" & LCD$held == L3_STRATUM)), ]
  z <- z[order(z$contrast != "L3", z$p95, z$median), ]
  cat("  population: patients whose leader DIFFERS between the two arms\n\n")
  cat(sprintf("  %-4s %-20s %-20s %7s %6s %6s %6s %6s %9s\n", "code", "method a",
              "method b", "differs", "median", "q3", "p90", "p95", "bottom3rd"))
  for (i in seq_len(nrow(z))) cat(sprintf(
    "  %-4s %-20s %-20s %7.4f %6.0f %6.0f %6.0f %6.0f %9.5f\n", z$contrast[i],
    z$method_a[i], z$method_b[i], z$frac_leader_differs[i], z$median[i], z$q3[i],
    z$p90[i], z$p95[i], z$frac_bottom_third[i]))
  cat("\n    Rank of method a's leading signal inside method b's ranking, among\n")
  cat("    patients whose leaders differ (`differs` is that share). The\n")
  cat("    `all_scored` rows in the table are the unconditional version.\n")
}

# ============================================================================
# SPECIFICATION AGAINST SAMPLING, WITHOUT DIVIDING
# ============================================================================
cat("\n=== L4 against L3: is the specification effect bigger than the sampling",
    "\n    noise it has to clear? ===\n\n")
cat("    `p_dominates` is P(L4 pair > L3 pair) + 0.5 P(equal), over all cross\n")
cat("    pairs, on the COMMON bag population and with every metric ORIENTED so\n")
cat("    that larger is more disagreement (share metrics enter as 1 - share).\n\n")
ALLRAW  <- do.call(rbind, Filter(Negate(is.null), list(RAW_W, COS_W)))
ALL4RAW <- do.call(rbind, Filter(Negate(is.null), list(L4RAW_W, COS_W)))
ALL4RAW <- if (is.null(ALL4RAW)) NULL else
  ALL4RAW[ALL4RAW$contrast == "L4" & ALL4RAW$bag_population == "common", , drop = FALSE]
l3raw <- if (is.null(ALLRAW)) data.frame() else
  ALLRAW[ALLRAW$contrast == "L3" & ALLRAW$route == "refit" &
         ALLRAW$held == L3_STRATUM & ALLRAW$bag_population == "common", , drop = FALSE]
DOM <- list()
if (!is.null(ALL4RAW) && nrow(l3raw) > 0L) {
  METRICS <- unique(ALL4RAW$metric)
  for (lv in LVLS) for (mt in METRICS) {
    z4 <- ALL4RAW[ALL4RAW$agg == lv & ALL4RAW$metric == mt, , drop = FALSE]
    if (!nrow(z4)) next
    k <- z4$k[1]
    u4 <- unique(z4[, c("method_a", "method_b")])
    for (i in seq_len(nrow(u4))) {
      u <- u4[i, ]
      x <- z4$value[z4$method_a == u$method_a & z4$method_b == u$method_b]
      for (side in c("a", "b")) {
        mm <- if (side == "a") u$method_a else u$method_b
        y <- l3raw$value[l3raw$method_a == mm & l3raw$agg == lv & l3raw$metric == mt]
        if (!length(y)) next
        ox <- attr_orient(x, mt); oy <- attr_orient(y, mt)
        DOM[[length(DOM) + 1L]] <- cbind(
          data.frame(agg = lv, metric = mt, metric_oriented = ox$metric,
                     orientation = ox$orientation, k = k,
                     method_a = u$method_a, method_b = u$method_b,
                     l3_reference = mm, bag_population = "common",
                     stringsAsFactors = FALSE),
          attr_dominance(ox$value, oy$value),
          attr_stability_ratio(ox$value, oy$value))
      }
    }
  }
}
if (length(DOM)) {
  DOM <- do.call(rbind, DOM); DOM$store_status <- STORE_STATUS
  save_table(run, DOM, "l4_against_l3", subdir = "diagnostics")
  for (mt in unique(DOM$metric_oriented)) {
    z <- DOM[DOM$agg == "signal" & DOM$metric_oriented == mt, ]
    z <- z[order(-z$p_dominates), ]
    cat(sprintf("\n  --- metric `%s`, signal level ---\n\n", mt))
    cat(sprintf("  %-20s %-20s %-20s %11s %11s %10s\n", "method a", "method b",
                "L3 reference", "p_dominates", "median gap", "p95 gap"))
    for (i in seq_len(nrow(z))) cat(sprintf(
      "  %-20s %-20s %-20s %11.3f %11.5f %10.5f\n", z$method_a[i], z$method_b[i],
      z$l3_reference[i], z$p_dominates[i], z$median_gap[i], z$p95_gap[i]))
  }
  cat("\n    `ratio_median` is the DEMOTED level-4-over-level-3 statistic, kept\n")
  cat("    so numbers reported before 2026-09-06 stay comparable. It decides\n")
  cat("    nothing; see the `selection:` block of config/attribution_eval.yml.\n")
} else {
  DOM <- NULL
  cat("  no L3 or no L4 replicates yet; generate stage `sample`.\n")
}

# ============================================================================
# THE PAIRED CROSS-METHOD LEVEL-3 COMPARISON
# ============================================================================
cat("\n=== L3, SHAP against each LLR arm, PAIRED ON THE SHARED BAG PAIRS ===\n\n")
cat("    Under IDENTICAL training-sample perturbations, whose patient-level\n")
cat("    explanation moves less? Both arms see the same bags, so the variation\n")
cat("    caused by which patients entered a bag is differenced out. Share\n")
cat("    metrics enter as 1 - share so that positive always means SHAP moved more.\n\n")
PAIRED <- list()
shap_m <- intersect(SCORE, ATTR_SHAP_ARMS)
llr_m  <- setdiff(SCORE, ATTR_SHAP_ARMS)
if (length(shap_m) && nrow(l3raw) > 0L) {
  for (sm in shap_m) for (lm in llr_m) for (lv in LVLS)
    for (mt in unique(l3raw$metric)) {
      a <- l3raw[l3raw$method_a == sm & l3raw$agg == lv & l3raw$metric == mt, ]
      b <- l3raw[l3raw$method_a == lm & l3raw$agg == lv & l3raw$metric == mt, ]
      if (!nrow(a) || !nrow(b)) next
      oa <- attr_orient(a$value, mt); ob <- attr_orient(b$value, mt)
      a$value <- oa$value; b$value <- ob$value
      PAIRED[[length(PAIRED) + 1L]] <- cbind(
        data.frame(shap = sm, llr = lm, agg = lv, metric = mt,
                   metric_oriented = oa$metric, orientation = oa$orientation,
                   k = a$k[1], bag_population = "common", stringsAsFactors = FALSE),
        attr_paired_delta(a, b))
    }
}
if (length(PAIRED)) {
  PAIRED <- do.call(rbind, PAIRED); PAIRED$store_status <- STORE_STATUS
  save_table(run, PAIRED, "l3_shap_minus_llr_paired", subdir = "diagnostics")
  for (mt in unique(PAIRED$metric_oriented)) {
    z <- PAIRED[PAIRED$agg == "domain" & PAIRED$metric_oriented == mt, ]
    if (!nrow(z)) next
    cat(sprintf("\n  --- metric `%s`, domain level ---\n\n", mt))
    cat(sprintf("  %-16s %-20s %6s %10s %10s %10s %9s\n", "shap", "llr",
                "pairs", "median", "p05", "p95", "frac_pos"))
    for (i in seq_len(nrow(z))) {
      cat(sprintf("  %-16s %-20s %6d %10.5f %10.5f %10.5f %9.4f\n",
                  z$shap[i], z$llr[i], z$n_shared_bag_pairs[i],
                  z$median[i], z$p05[i], z$p95[i], z$frac_positive[i]))
    }
  }
  cat("\n    Positive means SHAP disagrees with itself MORE than the LLR arm\n")
  cat("    does, on the identical pair of resamples.\n")
} else {
  PAIRED <- NULL
  cat("  needs both families at contrast L3 on shared bags.\n")
}

# ============================================================================
# THE SELECTION RULE
# ============================================================================
cat("\n=== the selection rule ===\n\n")
cat("    THE STABILITY AXIS IS INTRINSIC: a method's OWN level-3 disagreement\n")
cat("    distribution on the common bags, at its median and its p95. Read the\n")
cat("    PROVENANCE block of config/attribution_eval.yml before quoting a\n")
cat("    verdict: the rule remains EVALUATIVE ONLY.\n\n")
y_tr <- as.integer(tar_read(y_train))
SEL <- data.frame(method = SCORE, stringsAsFactors = FALSE)

# THE DISCRIMINATION AXIS READS THE MODEL'S SCORE, NOT THE ATTRIBUTION SUM
# (review finding A4). For an LLR arm the row sum of its L matrix IS the
# layer-1 score by construction (`llr_sum`, `llr_cond`, ...). For SHAP the
# stored matrix is 19 signal groups with the BIAS and the intervention groups
# dropped, so its row sum is neither the booster's margin nor its ranking; the
# generator stores the out-of-fold margin beside the ladder replicate and it
# is read here. The row-sum AUROC is still reported for EVERY method as
# `auroc_attr_sum`, named for what it is.
SEL$auroc_attr_sum <- vapply(SCORE, function(m) .auroc(rowSums(REP1[[m]]), y_tr), numeric(1))
SEL$auprc_attr_sum <- vapply(SCORE, function(m) .auprc(rowSums(REP1[[m]]), y_tr), numeric(1))
SEL$discrimination_source <- ifelse(attr_method_family(SCORE) == "llr",
                                    "row_sum_is_layer1_score", "stored_oof_margin")
SEL$auroc <- SEL$auroc_attr_sum; SEL$auprc <- SEL$auprc_attr_sum
SEL$margin_vs_targets_max_abs <- NA_real_
for (m in intersect(SCORE, ATTR_SHAP_ARMS)) {
  mp <- file.path(gen_d, sub, paste0(spec$key[spec$method == m][1], "__margin.qs2"))
  if (!file.exists(mp)) {
    if (!ALLOW_PARTIAL) {
      stop("the store carries no out-of-fold margin for `", m, "`, so its ",
           "discrimination cannot be scored from the model's prediction (the ",
           "19-group row sum is not it -- review finding A4). Stamp it:\n",
           "  Rscript tests/attr_replicates.R --levels spec --resume ", gen_d, "\n",
           "or pass --allow-partial to leave the axis unscored.", call. = FALSE)
    }
    SEL$auroc[SEL$method == m] <- NA_real_; SEL$auprc[SEL$method == m] <- NA_real_
    SEL$discrimination_source[SEL$method == m] <- "unscored_no_margin"
    next
  }
  mg <- qs2::qs_read(mp)
  # THE MARGIN IS CHECKED AS A REPLICATE IS (re-review): one finite value per
  # training stay, in the store's row order, from the same replicate key.
  if (!is.numeric(mg$eta) || length(mg$eta) != length(ids_ch) ||
      !identical(names(mg$eta), ids_ch) || !all(is.finite(mg$eta))) {
    stop("margin file for `", m, "` is not one finite logit per training stay ",
         "in stay_id order.", call. = FALSE)
  }
  if (!identical(mg$replicate_key, spec$key[spec$method == m][1])) {
    stop("margin file for `", m, "` names a different replicate key from the ",
         "ladder replicate in the manifest.", call. = FALSE)
  }
  SEL$auroc[SEL$method == m] <- .auroc(mg$eta, y_tr)
  SEL$auprc[SEL$method == m] <- .auprc(mg$eta, y_tr)
  # Cross-check against the pipeline's own OOF booster, which is fitted from
  # the same design, folds and seeds and should therefore be the same
  # booster: an informational gap, not a stop, because the STORE's margin is
  # the one that belongs to the store's replicate.
  tg <- tryCatch(tar_read(oof_xgb_feat), error = function(e) NULL)
  if (!is.null(tg) && length(tg$score) == length(mg$eta)) {
    SEL$margin_vs_targets_max_abs[SEL$method == m] <-
      max(abs(mg$eta - (tg$score + logit(mean(y_tr)))))
  }
}
SEL$auprc_lift <- round(SEL$auprc / mean(y_tr), 4)

# RESOLUTION IS EVALUATED AT EACH METHOD'S OWN NOISE-CALIBRATED TOLERANCE at
# the configured quantile, signal level.
res_at_noise <- function(m) {
  d <- noise_delta(m, "signal", RES_Q)
  if (is.na(d)) return(c(delta = NA_real_, res = NA_real_))
  t <- attr_tieset(attr_prep(REP1[[m]]), d)
  c(delta = d, res = t$frac_unique_leader[1])
}
RES <- vapply(SEL$method, res_at_noise, numeric(2))
SEL$resolution_delta   <- RES["delta", ]
SEL$frac_unique_leader <- RES["res", ]
save_table(run, data.frame(
  method = SEL$method, agg = "signal",
  delta_source = sprintf("noise_calibrated_q%g", RES_Q), delta = round(SEL$resolution_delta, 6),
  frac_unique_leader = SEL$frac_unique_leader, stringsAsFactors = FALSE),
  "resolution_at_noise_delta", subdir = "diagnostics")

#' A method's OWN level-3 distribution on its primary route, common bags,
#' signal level, top-1.
own_l3 <- function(m, what) {
  if (is.null(DIST_W)) return(NA_real_)
  z <- DIST_W[DIST_W$contrast == STAB_C & DIST_W$method_a == m &
              DIST_W$route == "refit" & DIST_W$held == L3_STRATUM &
              DIST_W$agg == "signal" & DIST_W$k == 1L &
              DIST_W$bag_population == "common", , drop = FALSE]
  if (nrow(z) > 1L) {
    abort_values(paste0("own_l3: more than one ", STAB_C, " cell for `", m,
                        "` at stratum ", L3_STRATUM), unique(z$held))
  }
  if (!nrow(z)) NA_real_ else z[[what]][1]
}
SEL$l3_median <- vapply(SEL$method, own_l3, numeric(1), what = "median")
SEL$l3_p95    <- vapply(SEL$method, own_l3, numeric(1), what = "p95")
SEL$l3_max    <- vapply(SEL$method, own_l3, numeric(1), what = "max")
SEL$l3_n      <- vapply(SEL$method, own_l3, numeric(1), what = "n")
SEL$bag_population <- "common"
SEL$store_status   <- STORE_STATUS

if (all(is.na(SEL$l3_median))) {
  cat(sprintf("  no `%s` replicates: the stability axis cannot be scored and the\n",
              STAB_C))
  cat("  verdict is WITHHELD -- two axes out of three is not the rule.\n\n")
  print(SEL[order(-SEL$auroc), c("method", "auroc", "auprc", "auprc_lift",
                                 "discrimination_source", "frac_unique_leader")],
        row.names = FALSE)
} else {
  SC <- attr_selection_score(SEL, cfg_req(ecfg, "selection"))
  SC$l3_max <- SEL$l3_max[match(SC$method, SEL$method)]
  SC$l3_n_pairs <- SEL$l3_n[match(SC$method, SEL$method)]
  SC$auprc      <- round(SEL$auprc[match(SC$method, SEL$method)], 5)
  SC$auprc_lift <- SEL$auprc_lift[match(SC$method, SEL$method)]
  SC$auroc_attr_sum <- round(SEL$auroc_attr_sum[match(SC$method, SEL$method)], 5)
  SC$discrimination_source <- SEL$discrimination_source[match(SC$method, SEL$method)]
  SC$bag_population <- "common"; SC$store_status <- STORE_STATUS
  print(SC, row.names = FALSE)
  save_table(run, SC, "selection_verdict", subdir = "diagnostics")
}
save_table(run, SEL, "selection_inputs", subdir = "diagnostics")

finalize_run(run, extra = list(
  generator = basename(gen_d),
  store_status = STORE_STATUS,
  store_fingerprint = attr_key_hash(.fp),
  n_replicates = nrow(man),
  n_planned = sum(COV$planned), n_missing = sum(COV$missing),
  n_tombstoned = sum(COV$tombstoned),
  common_bags = paste(COMMON_BAGS, collapse = ","),
  n_common_bags = length(COMMON_BAGS),
  n_within_pairs = nrow(CP),
  contrasts = paste(sort(unique(CP$code)), collapse = ","),
  stability_contrast = STAB_C,
  resolution_noise_quantile = RES_Q,
  selection_status = cfg_req(ecfg, "selection", "status"),
  level4_role = cfg_req(ecfg, "selection", "level4_role"),
  topk_definition = "tie_tolerant_at_delta_0 (attr_agree_k semantics, boundary ties forgiven)",
  resample_kind = ATTR_ESTIMAND_NOTES$value[ATTR_ESTIMAND_NOTES$field == "resample_kind"],
  methods = paste(SCORE, collapse = ",")))
cat(sprintf("\nwritten: %s\n", run$path))
