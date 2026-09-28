# tests/attribution_vs_shap.R ------------------------------------------------
# THE LADDER AGAINST ITS NOISE FLOOR, ON A COMMON SCALE.
#
# tests/shap_noise_floor.R measured how much a SHAP attribution disagrees with
# ITSELF when nothing changes but the random seed. tests/attribution_ties.R
# measured how much our L attribution changes across specifications. Putting the
# two numbers side by side is the point of both, and doing it naively is wrong.
#
# THE SCALE CONFOUND, AND WHY IT MATTERS. A SHAP explanation and an L
# explanation do not carry the same total evidence. MEASURED: the median patient
# has 3.24 nats of |SHAP| spread over 19 signals against 6.84 nats of |L_cond|.
# So an absolute tolerance of 0.25 nats forgives twice as much of a SHAP
# disagreement as it does of an L disagreement, and the absolute-delta curves
# are not comparable. That difference is not an artifact either -- it is the
# marginal-versus-conditional distinction. Our L's are one marginal model per
# bundle, so shared evidence is counted in every column that carries it, which
# is precisely the redundancy layer 2's Sigma-inverse exists to discount. SHAP
# over one joint model splits that shared evidence between the columns.
#
# So the comparison is made on the RELATIVE tolerance, where delta is a share of
# each patient's own total absolute evidence UNDER THAT METHOD. The absolute
# curves are still reported, labelled as not scale-matched.
#
# FITS NOTHING. Reads the SHAP matrices and the L ladder that the two upstream
# scripts saved.
#
# THE STALENESS GUARD, AND THE RUN THAT FORCED IT. On 2026-09-06 this script
# was run against `coupattr_20260905T184528`, a ladder built at 18:45 the
# previous evening -- BEFORE the 00:14 rebuild that set `bam.gamma` to 1.5 and
# refitted all 258 GAMs. Every SHAP row in the resulting table was current,
# because `xgb_feat`'s design is built from `pi_hat`, `delta` and `lambda` and
# none of those is fitted by `bam`. Every L row was from the previous pipeline.
# Nothing errored and nothing looked wrong: the two halves of one table were
# simply from two different models, and the only way to find out was to read a
# manifest and compare two timestamps by hand.
#
# `tests/coupling_attribution.R` already refused a mismatched cache through a
# fingerprint. The DEFECT WAS THAT THE CONSUMER HAD NO SUCH CHECK -- a guard on
# the producer protects the producer's own reuse and nothing else. The ladder
# now carries `design_key` (see `attr_design_key()`), this script recomputes it
# from the live config and the live fold assignment, and a mismatch is a STOP
# that names the differing field rather than a warning that scrolls past.
#
# Aggregates only (hard rule 1).
#
#   Rscript tests/attribution_vs_shap.R <shapfloor_run> <coupattr_run>
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
shap_d  <- if (length(args) >= 1L && nzchar(args[1])) args[1] else
             latest_run("shapfloor", require_complete = FALSE)
lad_d   <- if (length(args) >= 2L && nzchar(args[2])) args[2] else
             latest_run("coupattr", require_complete = FALSE)

sh  <- readRDS(file.path(shap_d, "tables", "shap_oof_groups.rds"))
lad <- readRDS(file.path(lad_d,  "tables", "l_oof_ladder.rds"))
stopifnot(identical(as.character(sh$stay_id), as.character(lad$stay_id)))

# The design the LADDER must have been built under. Read from the targets store
# rather than from `config/config.yml` alone, because the fold assignment is
# part of the design and only the store has it without a data load.
cfg     <- tar_read(cfg)
folds   <- tar_read(folds)
tr      <- tar_read(train_ids)
domains <- load_domains("config/domains.csv")
sigs    <- unlist(cfg$signals)

.guard_ladder <- function(lad, lad_d) {
  if (is.null(lad$design_key)) {
    stop("the ladder in ", basename(lad_d), " carries no `design_key`, so it ",
         "was written by a version of tests/coupling_attribution.R from before ",
         "2026-09-06 and predates the gamma rebuild. Re-run\n",
         "  Rscript tests/coupling_attribution.R\n",
         "with no --reuse and point this script at the new run.", call. = FALSE)
  }
  # `lad$k_ti` with NO FALLBACK. `design_key` and `k_ti` were stamped onto the
  # ladder by the same change, so a ladder that has one has the other, and a
  # `%||% 5L` here would be a third place the interaction basis dimension is
  # written down -- silently substituting 5 for whatever a future ladder used
  # and reporting the design as matching.
  if (is.null(lad$k_ti)) {
    stop("the ladder in ", basename(lad_d), " carries a design_key but no ",
         "`k_ti`; it was written by an inconsistent version of ",
         "tests/coupling_attribution.R. Re-run it.", call. = FALSE)
  }
  want <- attr_design_key(cfg, folds$fold[match(tr, folds$stay_id)], lad$k_ti)
  d <- attr_design_diff(lad$design_key, want)
  if (nrow(d)) {
    cat("\n!!! LADDER DESIGN MISMATCH !!!\n\n")
    print(d[, c("field", "hash_a", "hash_b", "n_elements_a", "n_elements_b")],
          row.names = FALSE)
    stop("the ladder in ", basename(lad_d), " was built under a different ",
         "design from the one config/config.yml and the targets store now ",
         "describe (", nrow(d), " field(s) differ, listed above). Mixing it ",
         "with a current SHAP floor produces a table whose two halves come ",
         "from different pipelines. Re-run tests/coupling_attribution.R with ",
         "no --reuse.", call. = FALSE)
  }
  cat(sprintf("  ladder design key %s matches the live design.\n",
              attr_key_hash(want)))
  invisible(TRUE)
}
cat("\n=== staleness guard ===\n\n")
.guard_ladder(lad, lad_d)

# The SHAP floor needs no design key: `xgb_feat`'s design is built from the
# layer-1 COVARIATES -- pi_hat, delta, lambda -- none of which is fitted by
# `bam`, so a GAM refit cannot move a number in it. What CAN move it is a
# change to those covariate constructions or to the fold assignment, and both
# are recorded in the floor run's own config snapshot. Checked here rather than
# assumed.
.guard_floor <- function(shap_d) {
  m <- read_manifest(shap_d)
  # Compared as HASHES OF THE IN-MEMORY CONFIG, not of the YAML snapshot: the
  # manifest's `config_hash` was computed by `write_manifest()` on the object
  # `tests/shap_noise_floor.R` held, and `load_config()` here rebuilds that same
  # object. Hashing the round-tripped YAML instead would differ on type coercion
  # alone and the guard would cry wolf on every run.
  live <- .hash(load_config("config/config.yml"))
  got  <- as.character(m$config_hash %||% "absent")
  if (!identical(live, got)) {
    cat(sprintf("  NOTE: the SHAP floor in %s was built under a different config\n",
                basename(shap_d)))
    cat("        snapshot than the live one. That is expected after a `bam`-only\n")
    cat("        change, which cannot reach a booster; it is NOT expected after a\n")
    cat("        change to pi_hat, delta, lambda, the folds or the xgb settings.\n")
    cat(sprintf("        floor config %s vs live %s -- check before quoting.\n", got, live))
  } else {
    cat(sprintf("  SHAP floor config snapshot %s matches the live config.\n", live))
  }
  invisible(TRUE)
}
.guard_floor(shap_d)

run <- new_run("attrshap", cfg, note = sprintf(
  "L ladder against the SHAP noise floor, scale-matched. shap=%s ladder=%s",
  basename(shap_d), basename(lad_d)))

DELTA_ABS <- c(0, 0.02, 0.05, 0.10, 0.25, 0.50)
DELTA_REL <- c(0, 0.01, 0.02, 0.05, 0.10)
KS <- c(1L, 3L)

# --- put every arm on the same 19 signal columns and 11 domains -------------
sig_cols <- intersect(sh$signal_groups, sigs)
dm  <- domains$domain[match(sig_cols, domains$signal)]
dnm <- sort(unique(dm))
to_dom <- function(M) {
  D <- matrix(0, nrow(M), length(dnm), dimnames = list(rownames(M), dnm))
  for (k in dnm) { j <- sig_cols[dm == k]
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE]) }
  D
}
sub <- function(M) M[, sig_cols, drop = FALSE]

PAIRS <- list(
  list(lab = "shap_seed1 vs shap_seed2", kind = "noise floor",
       a = sub(sh$shap_seed1),      b = sub(sh$shap_seed2)),
  list(lab = "cond vs cond_ti_all",     kind = "specification",
       a = sub(lad$arms$cond),      b = sub(lad$arms$cond_ti_all)),
  list(lab = "cond vs cond_ti_trend",   kind = "specification",
       a = sub(lad$arms$cond),      b = sub(lad$arms$cond_ti_trend)),
  list(lab = "meas vs cond",            kind = "specification",
       a = sub(lad$arms$meas),      b = sub(lad$arms$cond)))

prep <- function(M) {
  ab <- abs(M); n <- nrow(ab); ord <- t(apply(-ab, 1, order))
  list(ab = ab, n = n, total = rowSums(ab),
       kth = function(k) ab[cbind(seq_len(n), ord[, k])],
       inS = function(k) { S <- matrix(FALSE, n, ncol(ab))
         S[cbind(rep(seq_len(n), times = k), as.vector(ord[, seq_len(k)]))] <- TRUE; S })
}
agree_k <- function(pa, pb, k, delta) {
  Sa <- pa$inS(k); Sb <- pb$inS(k); ka <- pa$kth(k); kb <- pb$kth(k)
  !((rowSums(Sb & !Sa & (pa$ab < ka - delta)) > 0) |
    (rowSums(Sa & !Sb & (pb$ab < kb - delta)) > 0))
}

cat("\n=== evidence budgets, which is why the scales differ ===\n\n")
cat(sprintf("%-26s %-8s %12s %12s\n", "arm pair (first arm)", "level",
            "median budget", "median hhi"))
BB <- list()
for (P in PAIRS) {
  for (lv in c("signal", "domain")) {
    M <- if (lv == "signal") P$a else to_dom(P$a)
    ab <- abs(M); tot <- rowSums(ab); shr <- ab / pmax(tot, .Machine$double.eps)
    BB[[length(BB) + 1L]] <- data.frame(pair = P$lab, level = lv,
      budget_median = round(stats::median(tot), 4),
      hhi_median = round(stats::median(rowSums(shr^2)), 5), stringsAsFactors = FALSE)
    cat(sprintf("%-26s %-8s %12.4f %12.5f\n", P$lab, lv,
                BB[[length(BB)]]$budget_median, BB[[length(BB)]]$hhi_median))
  }
}
save_table(run, do.call(rbind, BB), "budget_by_method", subdir = "diagnostics")

emit <- function(kind, grid, rel) {
  cat(sprintf("\n=== %s tie tolerance %s ===\n", kind,
              if (rel) "(SCALE-MATCHED: share of each method's own budget)"
              else     "(NOT scale-matched; read the relative table instead)"))
  R <- list()
  for (lv in c("signal", "domain")) {
    for (k in KS) {
      cat(sprintf("\n  --- %s level, top-%d ---\n\n", lv, k))
      cat(sprintf("  %-26s %-14s", "comparison", "kind"))
      for (d in grid) cat(sprintf(" %9s", sprintf("%s%.2f", if (rel) "r=" else "d=", d)))
      cat("\n")
      for (P in PAIRS) {
        A <- if (lv == "signal") P$a else to_dom(P$a)
        B <- if (lv == "signal") P$b else to_dom(P$b)
        pa <- prep(A); pb <- prep(B)
        r <- vapply(grid, function(d)
          mean(agree_k(pa, pb, k, if (rel) d * pa$total else d)), numeric(1))
        R[[length(R) + 1L]] <- data.frame(comparison = P$lab, kind = P$kind,
          level = lv, k = k, delta_kind = if (rel) "relative_share" else "absolute_nats",
          delta = grid, agree = round(r, 5), stringsAsFactors = FALSE)
        cat(sprintf("  %-26s %-14s", P$lab, P$kind))
        for (v in r) cat(sprintf(" %9.4f", v)); cat("\n")
      }
    }
  }
  do.call(rbind, R)
}
RA <- emit("ABSOLUTE", DELTA_ABS, FALSE)
RR <- emit("RELATIVE", DELTA_REL, TRUE)
save_table(run, rbind(RA, RR), "ladder_vs_noise_floor", subdir = "diagnostics")

cat("\n=== how to read it ===\n\n")
cat("  A specification disagreement that is NO LARGER than the seed-change\n")
cat("  floor is not evidence that the specification matters -- it is what any\n")
cat("  attribution of this kind does. A specification disagreement that clearly\n")
cat("  EXCEEDS the floor is a real effect of the modelling choice.\n")

finalize_run(run, extra = list(
  shap_run = basename(shap_d), ladder_run = basename(lad_d),
  ladder_design_key = attr_key_hash(lad$design_key),
  guard = "ladder design key checked against the live config and folds"))
cat(sprintf("\nwritten: %s\n", run$path))
