# tests/attr_hospital_disagreement.R --------------------------------------------
# THE eICU REPRODUCIBILITY ARM, PER HOSPITAL. FITS NOTHING. Added 2026-10-06 at
# the author's request (paper/plans/plan_eicu_reproducibility_arm.md).
#
# The pooled metrics run (`tests/attr_metrics.R --site eicu`) gives, for every
# pair of replicates, one disagreement value over all 95,507 eICU stays. This
# script computes the same quantities WITHIN EACH eICU HOSPITAL, so that the
# reproducibility of an explanation can be read the way the paper reads
# discrimination: a median hospital, a 5th-95th percentile spread, and a range.
#
# WHAT IS COMPUTED, AND ON WHICH PAIRS. Everything mirrors the pooled run, so a
# per-hospital number and a pooled number are the same quantity on a different
# population of stays:
#   pairs       each within-method cell's pairs are the pooled run's own:
#               `attr_contrast_pairs()` on the store's manifest, the common
#               bags, bounded at `distribution.max_pairs` under `pair_seed`
#               (`cell_pairs()` of tests/attr_metrics.R). Method contrasts (L4)
#               are the common bags, one pair of methods per bag, plus the
#               anchor (the frozen bundle) reported separately.
#   top-k       strict top-1 and top-3 disagreement (the tie-tolerant
#               definition at delta = 0, `attr_topk()` encoding) on every pair.
#               At the noise threshold, on the pooled run's thresholded design:
#               the first `distribution.max_reps_delta` replicates (16) of the
#               L3 bootstrap stratum and of the common bags for L4, all their
#               pairs, each method's eICU-calibrated noise threshold (median,
#               bootstrap route, signal level), and for a pair of methods the
#               larger of the two thresholds.
#   leader      rank displacement and leader collapse (`attr_leader_cells()`),
#               both directions, on the pooled run's leader subsample
#               (`distribution.max_pairs_collapse` pairs per cell under
#               `pair_seed`; every common bag for L4). The measured mask applies
#               when both members are evidence arms, as in the pooled run.
#   level       signal level only (the level the paper reports).
#
# PER HOSPITAL, THEN ACROSS HOSPITALS. For each pair, a metric is summarized
# over the stays of each hospital (a share, or a median / 95th percentile of a
# displacement); for each hospital, over the pairs of the cell (median); and
# across hospitals, by the median, mean, 5th / 25th / 75th / 95th percentile,
# minimum and maximum. The same metric over all eligible stays pooled is
# reported beside it. For the sampling floors, SHAP's floor is compared with
# each evidence arm's WITHIN each hospital on the same bag pairs.
#
# THE HOSPITALS. The hospital table, the join rules (missing or conflicting
# assignments stop the run; duplicates with one assignment are de-duplicated;
# unmatched stays follow `hospital.unmatched_policy`) and the eligibility floors
# (`hospital.min_stays`, `hospital.min_events`, and as many survivors as the
# event floor) are those of `run/external.R` and `group_metrics()`, so the
# eligible hospitals are the ones the per-hospital AUROC uses. The outcome
# decides ELIGIBILITY ONLY; no metric here reads it.
#
# HARD RULE 1. Hospital identifiers are row-level under the PhysioNet DUA: the
# per-hospital table is written as an .rds object and never printed or written
# as CSV. Everything printed is a count or a summary across hospitals.
#
#   Rscript tests/attr_hospital_disagreement.R
#   Rscript tests/attr_hospital_disagreement.R --store out/runs/attrextgen_... --metrics out/runs/attrmetricsext_...
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml); library(qs2); library(arrow) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/attr_external_common.R")

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
STORE   <- .opt("--store", latest_run("attrextgen", require_complete = TRUE))
METRICS <- .opt("--metrics", latest_run("attrmetricsext", require_complete = TRUE))
MEM_LIMIT_GB <- 12
invisible(gc(reset = TRUE))
mem_max_gb <- function() { g <- gc(verbose = FALSE); sum(g[, ncol(g)]) / 1024 }
mem_check <- function(where) {
  mx <- mem_max_gb()
  cat(sprintf("  [memory] %-34s R max used %.2f GB\n", where, mx))
  if (mx > MEM_LIMIT_GB) stop("R's maximum memory passed ", MEM_LIMIT_GB, " GB at ", where, call. = FALSE)
}
t_all <- start_timer()

# --- the store, the metrics run, and the chain between them ---------------------
.need <- function(cond, msg) if (!isTRUE(cond)) stop("hospital disagreement: ", msg, call. = FALSE)
.need(!is.null(STORE) && dir.exists(STORE), "no complete eICU store (attrextgen_*)")
.need(!is.null(METRICS) && dir.exists(METRICS), "no complete eICU metrics run (attrmetricsext_*)")
sm <- read_manifest(STORE); mm <- read_manifest(METRICS)
.need(identical(sm$status, "complete"), paste(basename(STORE), "is not marked complete"))
.need(identical(mm$status, "complete") && identical(mm$site, "eicu"),
      paste(basename(METRICS), "is not a complete eICU metrics run"))
.need(identical(mm$generator, basename(STORE)),
      sprintf("the metrics run consumed %s, not %s", mm$generator, basename(STORE)))
dg <- qs2::qs_read(file.path(STORE, "design.qs2"))
.need(identical(dg$site, "eicu") && identical(dg$generator, EICU_GEN_VERSION), "not an eICU whole-bag store")

ecfg <- yaml::read_yaml("config/attribution_eval.yml")
xcfg <- yaml::read_yaml("config/external.yml")
bundle <- load_bundle(dg$bundle, verbose = FALSE)
cfg    <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
sigs   <- as.character(unlist(cfg$signals))
METHODS <- as.character(dg$methods)
ids_ch <- as.character(dg$stay_id)
n_e    <- length(ids_ch)
meas   <- qs2::qs_read(file.path(STORE, "measured.qs2"))
.need(identical(rownames(meas), ids_ch) && identical(colnames(meas), sigs), "measured mask out of order")
man <- utils::read.csv(file.path(STORE, "manifest_replicates.csv"), stringsAsFactors = FALSE)
attr_validate_manifest(man, "tests/attr_hospital_disagreement.R")
# `attr_contrast_pairs()` returns row indices into the manifest it is given, and
# they are read against `man` below; the store holds its methods only, so the
# two are the same table. Asserted rather than assumed.
.need(all(man$method %in% METHODS), "the store's manifest holds a method outside its design")
SUB <- cfg_req(ecfg, "storage", "subdir")
load_rep <- function(key) {
  M <- qs2::qs_read(file.path(STORE, SUB, paste0(key, ".qs2")))
  attr_check_replicate(M, ids_ch, sigs, key)
  M
}

KS        <- as.integer(cfg_req(ecfg, "aggregation", "top_k"))
D_MAXP    <- as.integer(cfg_req(ecfg, "distribution", "max_pairs"))
D_MAXR    <- as.integer(cfg_req(ecfg, "distribution", "max_reps_delta"))
D_SEED    <- as.integer(cfg_req(ecfg, "distribution", "pair_seed"))
COL_MAXP  <- as.integer(cfg_req(ecfg, "distribution", "max_pairs_collapse"))
COL_SH    <- as.numeric(cfg_req(ecfg, "delta", "collapse_shares"))
NOISE_RT  <- as.character(cfg_req(ecfg, "delta", "noise_route"))
L3_STRATUM <- as.character(cfg_req(ecfg, "selection", "stability_stratum"))
RES_Q     <- as.numeric(cfg_req(ecfg, "selection", "resolution_noise_quantile"))

# The eICU-calibrated noise thresholds, read from the metrics run rather than
# recomputed, so the per-hospital and pooled numbers share one threshold.
ND <- readRDS(file.path(METRICS, "diagnostics", "noise_calibrated_delta.rds"))
ND <- ND[ND$route == NOISE_RT & ND$agg == "signal" & ND$calibrates %in% TRUE &
           abs(ND$q - RES_Q) < 1e-12, , drop = FALSE]
DELTA <- stats::setNames(ND$delta, ND$method)
.need(all(METHODS %in% names(DELTA)), "a method has no calibrated noise threshold in the metrics run")

# --- the hospitals: the join and the floors of run/external.R ---------------------
hp   <- cfg_req(xcfg, "hospital")
gcol <- hp$group_col %||% "hospitalid"
hz   <- as.data.frame(arrow::read_parquet(cfg_req(xcfg, "hospital_table")))
.need(all(c("stay_id", gcol) %in% names(hz)), "hospital table lacks stay_id or the group column")
.need(!anyNA(hz[[gcol]]), "hospital table has missing hospital assignments")
ids_h <- as.character(hz$stay_id)
dup <- duplicated(ids_h) | duplicated(ids_h, fromLast = TRUE)
if (any(dup)) {
  conflict <- tapply(as.character(hz[[gcol]][dup]), ids_h[dup], function(v) length(unique(v)) > 1L)
  .need(!any(conflict), "hospital table assigns some stays to more than one hospital")
  hz <- hz[!duplicated(ids_h), , drop = FALSE]; ids_h <- as.character(hz$stay_id)
}
grp <- as.character(hz[[gcol]])[match(ids_ch, ids_h)]
n_unmatched <- sum(is.na(grp))
if (n_unmatched && identical(hp$unmatched_policy, "error")) {
  stop("hospital disagreement: ", n_unmatched, " stay(s) have no hospital and the policy is `error`.",
       call. = FALSE)
}
# Eligibility needs the outcome (the event floor); it is read for that alone.
tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
.need(identical(as.character(tabs_e$cohort$stay_id), ids_ch), "eICU cohort rows differ from the store's")
y <- as.integer(tabs_e$cohort$mortality)
rm(tabs_e); invisible(gc(verbose = FALSE))
HN <- tapply(!is.na(grp), grp, sum); HK <- tapply(y, grp, sum)
min_n <- as.integer(hp$min_stays); min_k <- as.integer(hp$min_events)
elig_h <- names(HN)[HN >= min_n & HK >= min_k & (HN - HK) >= min_k]
in_e   <- !is.na(grp) & grp %in% elig_h
H      <- factor(grp[in_e], levels = sort(elig_h))
nH     <- nlevels(H)
ELIG <- data.frame(n_stays_scored = n_e, n_unmatched = n_unmatched,
                   n_hospitals = length(HN), n_hospitals_eligible = nH,
                   n_stays_eligible = sum(in_e), min_stays = min_n, min_events = min_k,
                   unmatched_policy = hp$unmatched_policy, stringsAsFactors = FALSE)
cat(sprintf("\n=== per-hospital disagreement: %s (metrics %s) ===\n\n", basename(STORE), basename(METRICS)))
cat(sprintf("  hospitals: %d in the table, %d eligible (>= %d stays, >= %d deaths and survivors); %d eligible stays; %d unmatched\n",
            length(HN), nH, min_n, min_k, sum(in_e), n_unmatched))
rm(y)

# --- per-pair, per-hospital summaries ---------------------------------------------
#' Share of each hospital's eligible stays where `v` is TRUE (or the mean of a
#' numeric `v`), among the stays where `use` is TRUE; NA for a hospital with none.
h_share <- function(v, use = NULL) {
  v <- as.numeric(v[in_e]); u <- if (is.null(use)) rep(TRUE, length(v)) else use[in_e]
  num <- tapply(v[u], H[u], sum); den <- tapply(rep(1, sum(u)), H[u], sum)
  out <- as.numeric(num / den); out[is.na(den) | den == 0] <- NA_real_; out
}
h_quant <- function(v, use, pr) {
  v <- as.numeric(v[in_e]); u <- use[in_e]
  as.numeric(tapply(v[u], H[u], function(z) if (length(z)) stats::quantile(z, pr, names = FALSE) else NA_real_))
}
pooled_share <- function(v, use = NULL) {
  u <- in_e & (if (is.null(use)) TRUE else use); if (!any(u)) NA_real_ else mean(as.numeric(v[u]))
}
pooled_quant <- function(v, use, pr) {
  u <- in_e & use; if (!any(u)) NA_real_ else stats::quantile(as.numeric(v[u]), pr, names = FALSE)
}
#' Per-patient strict top-k disagreement from two `attr_topk()` encodings: the
#' per-patient form of `attr_topk_agree()` (identical definition).
topk_disagree <- function(ta, tb, k) {
  kk <- as.character(as.integer(k))
  Sa <- ta[[kk]]$top; Sb <- tb[[kk]]$top; Ta <- ta[[kk]]$tie; Tb <- tb[[kk]]$tie
  (bitwAnd(bitwAnd(Sb, bitwNot(Sa)), bitwNot(Ta)) != 0L) |
    (bitwAnd(bitwAnd(Sa, bitwNot(Sb)), bitwNot(Tb)) != 0L)
}
#' Leader metrics of one pair, per hospital and pooled, as a named list of
#' (per-hospital vector, pooled scalar).
leader_metrics <- function(A, B, keep) {
  cl <- attr_leader_cells(A, B, keep = keep)
  ok <- cl$ok; ch <- ok & !cl$same_leader
  out <- list(
    frac_leader_differs = list(h_share(!cl$same_leader, ok), pooled_share(!cl$same_leader, ok)),
    disp_ab_median  = list(h_quant(cl$disp_a, ok, 0.5),  pooled_quant(cl$disp_a, ok, 0.5)),
    disp_ab_p95     = list(h_quant(cl$disp_a, ok, 0.95), pooled_quant(cl$disp_a, ok, 0.95)),
    disp_ba_median  = list(h_quant(cl$disp_b, ok, 0.5),  pooled_quant(cl$disp_b, ok, 0.5)),
    disp_ba_p95     = list(h_quant(cl$disp_b, ok, 0.95), pooled_quant(cl$disp_b, ok, 0.95)),
    disp_ab_changed_median = list(h_quant(cl$disp_a, ch, 0.5), pooled_quant(cl$disp_a, ch, 0.5)),
    disp_ab_changed_p95    = list(h_quant(cl$disp_a, ch, 0.95), pooled_quant(cl$disp_a, ch, 0.95)),
    disp_ba_changed_median = list(h_quant(cl$disp_b, ch, 0.5), pooled_quant(cl$disp_b, ch, 0.5)),
    disp_ba_changed_p95    = list(h_quant(cl$disp_b, ch, 0.95), pooled_quant(cl$disp_b, ch, 0.95)))
  for (s in COL_SH) {
    out[[sprintf("collapse_ab_%g", s)]] <- list(h_share(cl$sh_a_in_b < s, ok), pooled_share(cl$sh_a_in_b < s, ok))
    out[[sprintf("collapse_ba_%g", s)]] <- list(h_share(cl$sh_b_in_a < s, ok), pooled_share(cl$sh_b_in_a < s, ok))
  }
  out
}
keep_pair <- function(a, b) if (attr_method_family(a) == "llr" && attr_method_family(b) == "llr") meas else NULL

# Accumulators: one row per (cell, metric, pair), holding the per-hospital
# vector and the pooled value; and the SHAP-minus-LLR floor differences.
ACC <- new.env(parent = emptyenv()); ACC$rows <- list()
acc_add <- function(cell, metric, pair_id, hv, pooled) {
  ACC$rows[[length(ACC$rows) + 1L]] <- list(cell = cell, metric = metric, pair = pair_id,
                                            h = hv, pooled = pooled)
}
cell_id <- function(contrast, a, b, route, held) paste(contrast, a, b, route, held, sep = "|")

# --- WITHIN-METHOD CELLS: the pooled run's pairs -----------------------------------
in_common_of <- function(man) {
  boot <- man[man$route == "bootstrap", ]
  bags <- lapply(stats::setNames(METHODS, METHODS), function(m) sort(unique(boot$boot_id[boot$method == m])))
  Reduce(intersect, bags)
}
COMMON_BAGS <- in_common_of(man)
in_common <- function(b) b == 0L | b %in% COMMON_BAGS
CP <- attr_contrast_pairs(man[man$method %in% METHODS, ])
CP$common <- in_common(CP$boot_i) & in_common(CP$boot_j)
cells <- unique(CP[, c("method", "route_family", "code", "held")])
cat(sprintf("  common bags: %d; within-method cells: %d\n", length(COMMON_BAGS), nrow(cells)))

for (ci in seq_len(nrow(cells))) {
  m <- cells$method[ci]; rt <- cells$route_family[ci]; cd <- cells$code[ci]; hd <- cells$held[ci]
  z <- CP[CP$method == m & CP$route_family == rt & CP$code == cd & CP$held == hd & CP$common, , drop = FALSE]
  if (!nrow(z)) next
  if (nrow(z) > D_MAXP) z <- z[with_seed(D_SEED, sort(sample.int(nrow(z), D_MAXP))), , drop = FALSE]
  zc <- if (nrow(z) <= COL_MAXP) z else z[with_seed(D_SEED, sort(sample.int(nrow(z), COL_MAXP))), , drop = FALSE]
  cid <- cell_id(cd, m, m, rt, hd)
  t0 <- start_timer()
  keys <- unique(c(man$key[z$i], man$key[z$j]))
  M  <- lapply(stats::setNames(keys, keys), load_rep)
  TK <- lapply(M, attr_topk, ks = KS)
  kp <- keep_pair(m, m)
  for (q in seq_len(nrow(z))) {
    ka <- man$key[z$i[q]]; kb <- man$key[z$j[q]]
    pid <- paste(z$boot_i[q], z$boot_j[q], man$seed_id[z$i[q]], man$seed_id[z$j[q]],
                 man$draw_id[z$i[q]], man$draw_id[z$j[q]], sep = ":")
    for (k in KS) {
      d <- topk_disagree(TK[[ka]], TK[[kb]], k)
      acc_add(cid, sprintf("top%d_strict", k), pid, h_share(d), pooled_share(d))
    }
  }
  for (q in seq_len(nrow(zc))) {
    ka <- man$key[zc$i[q]]; kb <- man$key[zc$j[q]]
    pid <- paste(zc$boot_i[q], zc$boot_j[q], man$seed_id[zc$i[q]], man$seed_id[zc$j[q]],
                 man$draw_id[zc$i[q]], man$draw_id[zc$j[q]], sep = ":")
    lm <- leader_metrics(M[[ka]], M[[kb]], kp)
    for (mt in names(lm)) acc_add(cid, mt, pid, lm[[mt]][[1]], lm[[mt]][[2]])
  }
  rm(M, TK); invisible(gc(verbose = FALSE))
  cat(sprintf("  %-4s %-14s %-9s %-16s %4d top-k pairs, %3d leader pairs (%.1f min)\n",
              cd, m, rt, hd, nrow(z), nrow(zc), t0()$elapsed_sec / 60))
}
mem_check("within-method cells")

# --- THE NOISE-THRESHOLD TOP-K: the pooled run's thresholded design --------------------
boot_man <- man[man$route == "bootstrap", ]
for (m in METHODS) {
  z <- boot_man[boot_man$method == m & boot_man$boot_id > 0L & in_common(boot_man$boot_id) &
                  boot_man$seed_id == 0L, , drop = FALSE]
  z <- z[order(z$boot_id), , drop = FALSE]
  z <- z[seq_len(min(nrow(z), D_MAXR)), , drop = FALSE]
  PP <- lapply(z$key, function(k) attr_prep(load_rep(k)))
  pix <- attr_pair_index(length(PP), D_MAXP, D_SEED)
  cid <- cell_id("L3", m, m, "refit", L3_STRATUM)
  for (q in seq_len(ncol(pix))) {
    a <- pix[1, q]; b <- pix[2, q]
    pid <- paste(z$boot_id[a], z$boot_id[b], 0, 0, 0, 0, sep = ":")
    for (k in KS) {
      d <- !attr_agree_k(PP[[a]], PP[[b]], k, DELTA[[m]])
      acc_add(cid, sprintf("top%d_noise", k), pid, h_share(d), pooled_share(d))
    }
  }
  rm(PP); invisible(gc(verbose = FALSE))
  cat(sprintf("  L3   %-14s noise threshold %.4f nats: %d pairs of the first %d bags\n",
              m, DELTA[[m]], ncol(pix), nrow(z)))
}
mem_check("noise-threshold floors")

# --- METHOD CONTRASTS (L4): every common bag, and the anchor -----------------------------
CMP <- utils::combn(METHODS, 2L, simplify = FALSE)
bags_n <- COMMON_BAGS[seq_len(min(length(COMMON_BAGS), D_MAXR))]
for (b in c(0L, COMMON_BAGS)) {
  zb <- if (b == 0L) man[man$route == "ladder" & man$method %in% METHODS, , drop = FALSE] else
    boot_man[boot_man$boot_id == b & boot_man$seed_id == 0L & boot_man$method %in% METHODS, , drop = FALSE]
  M  <- lapply(stats::setNames(zb$method, zb$method), function(mm) load_rep(zb$key[zb$method == mm][1]))
  TK <- lapply(M, attr_topk, ks = KS)
  PP <- if (b %in% bags_n) lapply(M, attr_prep) else NULL
  for (p in CMP) {
    cid <- cell_id(if (b == 0L) "L4_anchor" else "L4", p[1], p[2], "refit", "boot,seed,draw")
    pid <- as.character(b)
    for (k in KS) {
      d <- topk_disagree(TK[[p[1]]], TK[[p[2]]], k)
      acc_add(cid, sprintf("top%d_strict", k), pid, h_share(d), pooled_share(d))
      if (!is.null(PP)) {
        d2 <- !attr_agree_k(PP[[p[1]]], PP[[p[2]]], k, max(DELTA[[p[1]]], DELTA[[p[2]]]))
        acc_add(cid, sprintf("top%d_noise", k), pid, h_share(d2), pooled_share(d2))
      }
    }
    lm <- leader_metrics(M[[p[1]]], M[[p[2]]], keep_pair(p[1], p[2]))
    for (mt in names(lm)) acc_add(cid, mt, pid, lm[[mt]][[1]], lm[[mt]][[2]])
  }
  rm(M, TK, PP); invisible(gc(verbose = FALSE))
}
cat(sprintf("  L4   %d common bags and the anchor, %d pairs of methods\n", length(COMMON_BAGS), length(CMP)))
mem_check("method contrasts")

# --- ASSEMBLY ---------------------------------------------------------------------
R <- ACC$rows
meta <- data.frame(cell = vapply(R, `[[`, "", "cell"), metric = vapply(R, `[[`, "", "metric"),
                   pair = vapply(R, `[[`, "", "pair"), pooled = vapply(R, `[[`, 0, "pooled"),
                   stringsAsFactors = FALSE)
HM <- do.call(rbind, lapply(R, `[[`, "h"))          # rows: (cell, metric, pair); cols: hospitals
colnames(HM) <- levels(H)
rm(R, ACC); invisible(gc(verbose = FALSE))

# Per hospital: the median over the pairs of a cell. Then across hospitals.
grp_cm <- paste(meta$cell, meta$metric, sep = "##")
PH <- do.call(rbind, lapply(split(seq_len(nrow(meta)), grp_cm), function(ix) {
  apply(HM[ix, , drop = FALSE], 2, stats::median, na.rm = TRUE)
}))
PH[!is.finite(PH)] <- NA_real_
cm <- do.call(rbind, strsplit(rownames(PH), "##", fixed = TRUE))
cc <- do.call(rbind, strsplit(cm[, 1], "|", fixed = TRUE))
across <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(n_hospitals = 0, median = NA, mean = NA, q05 = NA, q25 = NA, q75 = NA,
                           q95 = NA, min = NA, max = NA))
  c(n_hospitals = length(v), median = stats::median(v), mean = mean(v),
    q05 = stats::quantile(v, 0.05, names = FALSE), q25 = stats::quantile(v, 0.25, names = FALSE),
    q75 = stats::quantile(v, 0.75, names = FALSE), q95 = stats::quantile(v, 0.95, names = FALSE),
    min = min(v), max = max(v))
}
pooled_med <- tapply(meta$pooled, grp_cm, stats::median, na.rm = TRUE)
n_pairs    <- tapply(meta$pair, grp_cm, length)
SUMM <- data.frame(contrast = cc[, 1], method_a = cc[, 2], method_b = cc[, 3], route = cc[, 4],
                   held = cc[, 5], metric = cm[, 2],
                   n_pairs = as.integer(n_pairs[rownames(PH)]),
                   pooled_eligible_median = as.numeric(pooled_med[rownames(PH)]),
                   t(apply(PH, 1, across)), stringsAsFactors = FALSE, row.names = NULL)
SUMM <- SUMM[order(SUMM$metric, SUMM$contrast, SUMM$method_a, SUMM$method_b, SUMM$held), ]

# The paired floor comparison, within each hospital: SHAP's sampling floor minus
# each evidence arm's, on the same bag pairs (L3 bootstrap stratum).
shap <- intersect(METHODS, ATTR_SHAP_ARMS); llr <- setdiff(METHODS, ATTR_SHAP_ARMS)
PAIRED <- list(); PAIRED_H <- list()
#' The unordered bag pair of a pair id, so the two families' pairs match
#' whatever order each listed its two bags in.
bagpair <- function(p) {
  b <- do.call(rbind, lapply(strsplit(p, ":", fixed = TRUE), function(v) as.integer(v[1:2])))
  paste(pmin(b[, 1], b[, 2]), pmax(b[, 1], b[, 2]), sep = ":")
}
for (sm_ in shap) for (lm_ in llr) for (mt in c("top1_strict", "top3_strict", "top1_noise", "top3_noise")) {
  ia <- which(meta$cell == cell_id("L3", sm_, sm_, "refit", L3_STRATUM) & meta$metric == mt)
  ib <- which(meta$cell == cell_id("L3", lm_, lm_, "refit", L3_STRATUM) & meta$metric == mt)
  pa <- bagpair(meta$pair[ia]); pb <- bagpair(meta$pair[ib])
  sh <- intersect(pa, pb)
  if (!length(sh)) next
  D <- HM[ia[match(sh, pa)], , drop = FALSE] - HM[ib[match(sh, pb)], , drop = FALSE]
  med_h  <- apply(D, 2, stats::median, na.rm = TRUE)
  fpos_h <- apply(D, 2, function(v) mean(v > 0, na.rm = TRUE))
  PAIRED_H[[paste(sm_, lm_, mt)]] <- data.frame(hospital = levels(H), shap = sm_, llr = lm_, metric = mt,
                                                median_diff = med_h, frac_pairs_shap_higher = fpos_h,
                                                stringsAsFactors = FALSE)
  dp <- meta$pooled[ia[match(sh, pa)]] - meta$pooled[ib[match(sh, pb)]]
  PAIRED[[length(PAIRED) + 1L]] <- data.frame(
    shap = sm_, llr = lm_, metric = mt, n_bag_pairs = length(sh),
    n_hospitals = sum(is.finite(med_h)),
    hospitals_shap_higher = sum(med_h > 0, na.rm = TRUE),
    hospitals_equal = sum(med_h == 0, na.rm = TRUE),
    hospitals_llr_higher = sum(med_h < 0, na.rm = TRUE),
    t(across(med_h))[, -1, drop = FALSE],
    pooled_eligible_median_diff = stats::median(dp, na.rm = TRUE),
    pooled_eligible_frac_shap_higher = mean(dp > 0, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
PAIRED <- do.call(rbind, PAIRED)

# --- OUTPUT -------------------------------------------------------------------------
run <- new_run("attrhosp", cfg, note = sprintf(
  "eICU attribution reproducibility per hospital: %s, %s; fits nothing", basename(STORE), basename(METRICS)))
save_table(run, ELIG, "hospital_eligibility", subdir = "diagnostics")
save_table(run, SUMM, "hospital_disagreement_summary", subdir = "tables")
save_table(run, PAIRED, "hospital_paired_floor", subdir = "tables")
# ROW-LEVEL UNDER THE DUA: hospital identifiers. Objects only, never printed or CSV.
PERH <- data.frame(cell = rep(sub("##.*$", "", rownames(PH)), times = ncol(PH)),
                   metric = rep(sub("^.*##", "", rownames(PH)), times = ncol(PH)),
                   hospital = rep(colnames(PH), each = nrow(PH)),
                   value = as.vector(PH), stringsAsFactors = FALSE)
save_object(run, PERH, "hospital_per_hospital_medians")
save_object(run, do.call(rbind, PAIRED_H), "hospital_paired_floor_per_hospital")

show <- SUMM[SUMM$metric %in% c("top1_strict", "top1_noise") &
               SUMM$contrast %in% c("L3", "L4", "L4_anchor", "L2", "L3T", "L3P"),
             c("contrast", "method_a", "method_b", "held", "metric", "pooled_eligible_median",
               "n_hospitals", "median", "q05", "q95", "min", "max")]
cat("\n=== top-1 disagreement per hospital (median over pairs within a hospital), across hospitals ===\n\n")
print(show, row.names = FALSE, digits = 3)
cat("\n=== sampling floor, SHAP minus each evidence arm, within each hospital ===\n\n")
print(PAIRED[, c("shap", "llr", "metric", "n_bag_pairs", "n_hospitals", "hospitals_shap_higher",
                 "hospitals_equal", "hospitals_llr_higher", "median", "q05", "q95")],
      row.names = FALSE, digits = 3)
.mx <- mem_max_gb()
finalize_run(run, extra = list(store = basename(STORE), metrics = basename(METRICS),
                               n_hospitals_eligible = nH, n_stays_eligible = sum(in_e),
                               common_bags = paste(COMMON_BAGS, collapse = ","),
                               methods = paste(METHODS, collapse = ","), level = "signal",
                               r_max_used_gb = round(.mx, 2)))
cat(sprintf("\nwritten: %s  (%.1f min, R max used %.2f GB)\n", run$path, t_all()$elapsed_sec / 60, .mx))
