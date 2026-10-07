# tests/attr_collapse_unmeasured.R ----------------------------------------------
# WHY DOES A SHAP LEADER COLLAPSE IN THE PAIRED WEIGHT OF EVIDENCE? FITS NOTHING.
# Added 2026-10-06 at the author's request.
#
# The per-hospital run (`tests/attr_hospital_disagreement.R`) found that the
# share of SHAP leaders that carry under 2% of the paired weight of evidence
# (leader collapse, paired WoE vs SHAP, the 2% rule of `attr_leader_cells()`)
# is 38% pooled at eICU and spreads from about 8% to 82% across hospitals,
# against 15% at MIMIC-IV. The hypothesis under test: SHAP assigns a
# contribution to an UNMEASURED channel (the booster routes a missing feature
# through its default direction), so SHAP's leader can be a channel whose
# paired weight of evidence is zero by assignment -- a collapse by
# construction, more frequent where a site or a hospital charts less.
#
# THE DECOMPOSITION. For every scored patient of a pair (the paired WoE as
# member a, SHAP as member b, unmasked, exactly as the pooled run forms this
# contrast), SHAP's leader is classified as a measured or an unmeasured channel
# of that stay. Reported, for each pair and then as the median over the common
# bags (and for the anchor):
#   collapse                 share of SHAP leaders that collapse (< 2%)
#   leader_unmeasured        share of SHAP leaders that are unmeasured channels
#   collapse_given_unmeas    collapse rate among unmeasured SHAP leaders
#   collapse_given_meas      collapse rate among measured SHAP leaders
#   collapses_unmeasured     share of all collapses whose leader is unmeasured
# at MIMIC-IV (the internal store) and at eICU (the whole-bag store), pooled;
# at eICU also per eligible hospital (the join and floors of run/external.R),
# with the across-hospital Spearman correlation between the collapse rate and
# the unmeasured-leader share; and, at the anchor, which channels the
# unmeasured SHAP leaders are.
#
# HARD RULE 1: hospital identifiers stay inside an .rds object; everything
# printed is a count, a share or a summary across hospitals.
#
#   Rscript tests/attr_collapse_unmeasured.R
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
EICU_STORE <- .opt("--eicu", "out/runs/attrextgen_20261006T082355")
MIMIC_STORE <- .opt("--mimic", "out/runs/attrgen_20260909T115047")
HOSP_RUN   <- .opt("--hosp", "out/runs/attrhosp_20261006T133213")
A_M <- "llr_full"; B_M <- "shap_xgb_feat"
ecfg <- yaml::read_yaml("config/attribution_eval.yml")
xcfg <- yaml::read_yaml("config/external.yml")
SUB  <- cfg_req(ecfg, "storage", "subdir")
SHARE <- 0.02

#' The per-patient decomposition of one pair (paired WoE `A`, SHAP `B`), with the
#' leader and shares of `attr_leader_cells()` (unmasked, as the pooled run forms
#' an LLR-against-SHAP pair).
decompose <- function(A, B, meas) {
  cl <- attr_leader_cells(A, B, keep = NULL)
  lb <- max.col(abs(B), ties.method = "first")        # SHAP's leader, as in attr_leader_cells
  unm <- !meas[cbind(seq_len(nrow(B)), lb)]
  list(ok = cl$ok, collapse = cl$sh_b_in_a < SHARE, unmeasured = unm, leader = lb)
}
summ <- function(d, use = d$ok) {
  c(n = sum(use),
    collapse = mean(d$collapse[use]),
    leader_unmeasured = mean(d$unmeasured[use]),
    collapse_given_unmeas = if (any(use & d$unmeasured)) mean(d$collapse[use & d$unmeasured]) else NA,
    collapse_given_meas = if (any(use & !d$unmeasured)) mean(d$collapse[use & !d$unmeasured]) else NA,
    collapses_unmeasured = if (any(use & d$collapse)) mean(d$unmeasured[use & d$collapse]) else NA)
}

#' One store: the anchor and every common bag.
run_store <- function(store) {
  man  <- utils::read.csv(file.path(store, "manifest_replicates.csv"), stringsAsFactors = FALSE)
  meas <- qs2::qs_read(file.path(store, "measured.qs2"))
  rd   <- function(key) qs2::qs_read(file.path(store, SUB, paste0(key, ".qs2")))
  ba <- sort(unique(man$boot_id[man$method == A_M & man$route == "bootstrap"]))
  bb <- sort(unique(man$boot_id[man$method == B_M & man$route == "bootstrap" & man$seed_id == 0L]))
  bags <- intersect(ba, bb)
  key_of <- function(m, b) if (b == 0L) man$key[man$method == m & man$route == "ladder"][1] else
    man$key[man$method == m & man$route == "bootstrap" & man$boot_id == b & man$seed_id == 0L][1]
  out <- list()
  for (b in c(0L, bags)) {
    A <- rd(key_of(A_M, b)); B <- rd(key_of(B_M, b))
    stopifnot(identical(dimnames(A), dimnames(meas)), identical(dimnames(B), dimnames(meas)))
    out[[as.character(b)]] <- decompose(A, B, meas)
  }
  list(dec = out, meas = meas, bags = bags)
}

# --- MIMIC-IV, pooled -------------------------------------------------------------
cat("\n=== SHAP-leader collapse in the paired WoE, decomposed by measured / unmeasured SHAP leader ===\n")
MI <- run_store(MIMIC_STORE)
SITE_ROWS <- list()
site_rows <- function(site, R, use_fun = function(d) d$ok) {
  S <- t(vapply(R$dec, function(d) summ(d, use_fun(d)), numeric(6)))
  anchor <- S["0", ]; bags <- S[rownames(S) != "0", , drop = FALSE]
  rbind(data.frame(site = site, population = "anchor", n_bags = 0L, t(anchor), check.names = FALSE),
        data.frame(site = site, population = "median over common bags", n_bags = nrow(bags),
                   t(apply(bags, 2, stats::median)), check.names = FALSE))
}
SITE_ROWS[["mimic"]] <- site_rows("MIMIC-IV train (out of fold)", MI)

# --- eICU, pooled and per hospital ---------------------------------------------------
EI <- run_store(EICU_STORE)
dg <- qs2::qs_read(file.path(EICU_STORE, "design.qs2"))
bundle <- load_bundle(dg$bundle, verbose = FALSE)
cfg <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
ids_ch <- rownames(EI$meas)
hp <- cfg_req(xcfg, "hospital"); gcol <- hp$group_col %||% "hospitalid"
hz <- as.data.frame(arrow::read_parquet(cfg_req(xcfg, "hospital_table")))
stopifnot(!anyNA(hz[[gcol]]))
ids_h <- as.character(hz$stay_id)
if (anyDuplicated(ids_h)) {
  dup <- duplicated(ids_h) | duplicated(ids_h, fromLast = TRUE)
  stopifnot(!any(tapply(as.character(hz[[gcol]][dup]), ids_h[dup], function(v) length(unique(v)) > 1L)))
  hz <- hz[!duplicated(ids_h), , drop = FALSE]; ids_h <- as.character(hz$stay_id)
}
grp <- as.character(hz[[gcol]])[match(ids_ch, ids_h)]
tabs_e <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
stopifnot(identical(as.character(tabs_e$cohort$stay_id), ids_ch))
y <- as.integer(tabs_e$cohort$mortality); rm(tabs_e)
HN <- tapply(!is.na(grp), grp, sum); HK <- tapply(y, grp, sum); rm(y)
min_n <- as.integer(hp$min_stays); min_k <- as.integer(hp$min_events)
elig <- names(HN)[HN >= min_n & HK >= min_k & (HN - HK) >= min_k]
in_e <- !is.na(grp) & grp %in% elig
# Tied to the per-hospital run: the same eligible set, by count.
.he <- readRDS(file.path(HOSP_RUN, "diagnostics", "hospital_eligibility.rds"))
stopifnot(length(elig) == .he$n_hospitals_eligible, sum(in_e) == .he$n_stays_eligible)
cat(sprintf("  eICU: %d eligible hospitals, %d eligible stays (as in %s)\n",
            length(elig), sum(in_e), basename(HOSP_RUN)))

SITE_ROWS[["eicu"]] <- site_rows("eICU (all stays)", EI)
SITE_ROWS[["eicu_elig"]] <- site_rows("eICU (eligible hospitals' stays)", EI, function(d) d$ok & in_e)
SITE <- do.call(rbind, SITE_ROWS); rownames(SITE) <- NULL

# Per hospital: each metric per bag, then the median over the common bags.
H <- factor(grp[in_e], levels = sort(elig))
per_bag_h <- lapply(EI$dec[names(EI$dec) != "0"], function(d) {
  u <- d$ok[in_e]
  f <- function(v, w) as.numeric(tapply(v[w], H[w], mean))
  cbind(collapse = f(d$collapse[in_e], u),
        leader_unmeasured = f(d$unmeasured[in_e], u),
        collapse_given_meas = f(d$collapse[in_e], u & !d$unmeasured[in_e]),
        collapse_given_unmeas = f(d$collapse[in_e], u & d$unmeasured[in_e]),
        n_unmeasured_channels = as.numeric(tapply(rowSums(!EI$meas[in_e, , drop = FALSE]), H, mean)))
})
PH <- apply(simplify2array(per_bag_h), c(1, 2), stats::median, na.rm = TRUE)
rownames(PH) <- levels(H)
across <- function(v) { v <- v[is.finite(v)]
  c(n_hospitals = length(v), median = stats::median(v), q05 = stats::quantile(v, 0.05, names = FALSE),
    q25 = stats::quantile(v, 0.25, names = FALSE), q75 = stats::quantile(v, 0.75, names = FALSE),
    q95 = stats::quantile(v, 0.95, names = FALSE), min = min(v), max = max(v)) }
HS <- data.frame(metric = colnames(PH), t(apply(PH, 2, across)), check.names = FALSE, row.names = NULL)
COR <- data.frame(
  x = c("leader_unmeasured", "n_unmeasured_channels", "leader_unmeasured"),
  y = c("collapse", "collapse", "collapse_given_meas"),
  spearman = c(stats::cor(PH[, "leader_unmeasured"], PH[, "collapse"], method = "spearman", use = "complete.obs"),
               stats::cor(PH[, "n_unmeasured_channels"], PH[, "collapse"], method = "spearman", use = "complete.obs"),
               stats::cor(PH[, "leader_unmeasured"], PH[, "collapse_given_meas"], method = "spearman", use = "complete.obs")),
  n_hospitals = nrow(PH), stringsAsFactors = FALSE)

# Which channels are SHAP's unmeasured leaders, at the anchor (counts and shares).
chan <- function(R, use = NULL, site) {
  d <- R$dec[["0"]]; u <- d$ok & d$unmeasured & (if (is.null(use)) TRUE else use)
  tb <- table(factor(colnames(R$meas)[d$leader[u]], levels = colnames(R$meas)))
  tb <- tb[tb > 0]
  data.frame(site = site, channel = names(tb), n_unmeasured_leaders = as.integer(tb),
             share_of_unmeasured_leaders = round(as.numeric(tb) / sum(tb), 4),
             collapse_rate = vapply(names(tb), function(ch)
               mean(d$collapse[u & colnames(R$meas)[d$leader] == ch]), numeric(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}
CH <- rbind(chan(MI, NULL, "MIMIC-IV"), chan(EI, in_e, "eICU (eligible)"))
CH <- CH[order(CH$site, -CH$n_unmeasured_leaders), ]

# --- OUTPUT ---------------------------------------------------------------------------
run <- new_run("attrcollapse", cfg, note = "SHAP-leader collapse in the paired WoE, by measured/unmeasured leader; fits nothing")
save_table(run, SITE, "collapse_decomposition_pooled", subdir = "tables")
save_table(run, HS, "collapse_decomposition_across_hospitals", subdir = "tables")
save_table(run, COR, "collapse_hospital_correlations", subdir = "tables")
save_table(run, CH, "unmeasured_shap_leaders_by_channel", subdir = "tables")
save_object(run, data.frame(hospital = rownames(PH), PH, row.names = NULL), "collapse_per_hospital")  # IDs: object only
cat("\n--- pooled (paired WoE as a, SHAP as b; collapse = SHAP leader under 2% of the paired WoE) ---\n\n")
print(SITE, row.names = FALSE, digits = 3)
cat("\n--- eICU, per eligible hospital (median over the common bags), across hospitals ---\n\n")
print(HS, row.names = FALSE, digits = 3)
cat("\n--- across hospitals, Spearman ---\n\n"); print(COR, row.names = FALSE, digits = 3)
cat("\n--- unmeasured SHAP leaders by channel, at the anchor ---\n\n"); print(CH, row.names = FALSE, digits = 3)
finalize_run(run, extra = list(eicu_store = basename(EICU_STORE), mimic_store = basename(MIMIC_STORE),
                               hosp_run = basename(HOSP_RUN), n_common_bags_eicu = length(EI$bags),
                               n_common_bags_mimic = length(MI$bags), share = SHARE))
cat(sprintf("\nwritten: %s\n", run$path))
