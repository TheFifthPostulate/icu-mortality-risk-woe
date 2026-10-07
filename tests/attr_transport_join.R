# tests/attr_transport_join.R ---------------------------------------------------
# THE eICU REPRODUCIBILITY ARM: TRANSPORT JOIN. FITS NOTHING. Added 2026-10-06
# (paper/plans/plan_eicu_reproducibility_arm.md).
#
# Sets every distribution of an eICU metrics run (`tests/attr_metrics.R --site
# eicu`, prefix `attrmetricsext`) beside the same distribution of a MIMIC-IV
# metrics run (prefix `attrmetrics`), row for row on the tables' key columns,
# for the eICU store's three methods. The two runs describe the same experiment
# only if the eICU store mirrors the internal store the MIMIC-IV run consumed,
# so that chain is checked before anything is joined.
#
# A NOISE-CALIBRATED TOLERANCE IS A VALUE, NOT A KEY. Each site calibrates its
# own tolerance from its own bag replicates, so the `delta` of a
# `noise_calibrated_*` row differs between the sites and is carried as a value
# (`delta_mimic`, `delta_eicu`); the fixed grids (absolute, relative, sign taus)
# keep their delta in the key.
#
#   Rscript tests/attr_transport_join.R
#   Rscript tests/attr_transport_join.R --mimic out/runs/attrmetrics_... --eicu out/runs/attrmetricsext_...
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/attr_external_common.R")

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(nm, default = NULL) {
  i <- which(args == nm)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
MIMIC <- .opt("--mimic", "out/runs/attrmetrics_20260909T165755")
EICU  <- .opt("--eicu", latest_run("attrmetricsext", require_complete = TRUE))
if (is.null(EICU) || !dir.exists(EICU)) stop("no complete eICU metrics run (attrmetricsext_*).", call. = FALSE)
if (!dir.exists(MIMIC)) stop("no such MIMIC-IV metrics run: ", MIMIC, call. = FALSE)

# --- THE CHAIN ------------------------------------------------------------------
mm <- read_manifest(MIMIC); em <- read_manifest(EICU)
.need <- function(cond, msg) if (!isTRUE(cond)) stop("transport join: ", msg, call. = FALSE)
.need(identical(mm$status, "complete"), paste(basename(MIMIC), "is not marked complete"))
.need(identical(em$status, "complete"), paste(basename(EICU), "is not marked complete"))
.need(is.null(mm$site) || identical(mm$site, "mimic"), paste(basename(MIMIC), "is not a MIMIC-IV run"))
.need(identical(em$site, "eicu"), paste(basename(EICU), "is not an eICU run"))
.need(identical(mm$store_status, "complete") && identical(em$store_status, "complete"),
      "both runs must have consumed a complete store")
.eg <- qs2::qs_read(file.path(dirname(EICU), em$generator, "design.qs2"))
.need(identical(.eg$site, "eicu") && identical(.eg$generator, EICU_GEN_VERSION),
      paste(em$generator, "is not an eICU whole-bag store"))
.need(identical(.eg$internal_store, mm$generator),
      sprintf("the eICU store mirrors %s but the MIMIC-IV run consumed %s", .eg$internal_store, mm$generator))
METHODS <- as.character(.eg$methods)
.mm_methods <- strsplit(mm$methods %||% "", ",", fixed = TRUE)[[1]]
.em_methods <- strsplit(em$methods %||% "", ",", fixed = TRUE)[[1]]
.need(setequal(.em_methods, METHODS), "the eICU run's methods are not its store's")
.need(all(METHODS %in% .mm_methods), "the MIMIC-IV run does not cover the eICU methods")
cat(sprintf("\n=== transport join: %s (MIMIC-IV) and %s (eICU) ===\n", basename(MIMIC), basename(EICU)))
cat(sprintf("  chain: eICU store %s mirrors %s, which %s consumed; methods %s\n\n",
            em$generator, .eg$internal_store, basename(MIMIC), paste(METHODS, collapse = ", ")))

# --- THE TABLES AND THEIR KEYS ------------------------------------------------------
KEYS <- list(
  disagreement_distributions     = c("contrast", "hierarchy_level", "method_a", "method_b", "route",
                                     "held", "agg", "k", "delta_kind", "delta_key", "bag_population"),
  pairwise_metric_distributions  = c("contrast", "method_a", "method_b", "route", "held", "agg",
                                     "metric", "bag_population", "orientation"),
  magnitude_metric_distributions = c("contrast", "method_a", "method_b", "route", "held", "agg",
                                     "bag_population", "k", "delta_kind", "delta_key", "metric"),
  leader_collapse_distributions  = c("contrast", "held", "method_a", "method_b", "agg", "route",
                                     "bag_population", "masked", "metric", "direction", "population",
                                     "unit", "tie_policy"),
  l3_shap_minus_llr_paired       = c("shap", "llr", "agg", "metric", "metric_oriented", "orientation",
                                     "k", "bag_population"),
  l4_against_l3                  = c("agg", "metric", "metric_oriented", "orientation", "k", "method_a",
                                     "method_b", "l3_reference", "bag_population"),
  l1_budget                      = c("method", "agg", "contrast"),
  l1_resolution                  = c("method", "agg", "delta"),
  resolution_at_noise_delta      = c("method", "agg", "delta_source"),
  noise_calibrated_delta         = c("method", "route", "agg", "calibrates", "q"),
  bag_population                 = c("method"))
DROP <- c("site", "store_status")
METHOD_COLS <- c("method", "method_a", "method_b", "shap", "llr", "l3_reference")

read_tab <- function(run, name) {
  p <- file.path(run, "diagnostics", paste0(name, ".rds"))
  if (!file.exists(p)) stop("transport join: ", basename(run), " has no table ", name, call. = FALSE)
  x <- readRDS(p)
  names(x) <- make.unique(names(x))          # leader_collapse_distributions repeats `n_pairs`
  if ("delta_kind" %in% names(x)) {
    x$delta_key <- ifelse(grepl("^noise_calibrated", x$delta_kind), x$delta_kind,
                          sprintf("%s=%s", x$delta_kind, format(x$delta, digits = 10, trim = TRUE)))
  }
  # Restricted to the eICU store's methods: every method column present is
  # empty, missing, or one of them.
  keep <- rep(TRUE, nrow(x))
  for (mc in intersect(METHOD_COLS, names(x))) {
    v <- as.character(x[[mc]])
    keep <- keep & (is.na(v) | !nzchar(v) | v %in% METHODS)
  }
  x[keep, setdiff(names(x), DROP), drop = FALSE]
}

run <- new_run("attrtransport", list(), note = sprintf(
  "eICU against MIMIC-IV attribution distributions: %s and %s; fits nothing",
  basename(MIMIC), basename(EICU)))
COVJ <- list(); J <- list()
for (nm in names(KEYS)) {
  a <- read_tab(MIMIC, nm); b <- read_tab(EICU, nm)
  k <- KEYS[[nm]]
  miss <- setdiff(k, intersect(names(a), names(b)))
  if (length(miss)) stop("transport join: table ", nm, " lacks key column(s) ",
                         paste(miss, collapse = ", "), call. = FALSE)
  for (s in list(list(a, "MIMIC-IV"), list(b, "eICU"))) {
    if (anyDuplicated(s[[1]][, k, drop = FALSE])) {
      stop("transport join: table ", nm, " has duplicated keys at ", s[[2]],
           "; the key set is incomplete.", call. = FALSE)
    }
  }
  vals <- intersect(setdiff(names(a), k), setdiff(names(b), k))
  j <- merge(a[, c(k, vals), drop = FALSE], b[, c(k, vals), drop = FALSE], by = k, all = TRUE,
             suffixes = c("_mimic", "_eicu"))
  both <- !is.na(j[[paste0(vals[1], "_mimic")]]) & !is.na(j[[paste0(vals[1], "_eicu")]])
  COVJ[[nm]] <- data.frame(table = nm, rows_mimic = nrow(a), rows_eicu = nrow(b),
                           rows_joined = sum(both), rows_mimic_only = sum(!is.na(j[[paste0(vals[1], "_mimic")]]) & !both),
                           rows_eicu_only = sum(!is.na(j[[paste0(vals[1], "_eicu")]]) & !both),
                           stringsAsFactors = FALSE)
  J[[nm]] <- j
  save_table(run, j, nm, subdir = "tables")
}
COVJ <- do.call(rbind, COVJ)
save_table(run, COVJ, "join_coverage", subdir = "diagnostics")
cat("=== join coverage ===\n\n"); print(COVJ, row.names = FALSE)

# --- THE HEADLINE: top-1 disagreement at signal level, common bags ------------------
d <- J$disagreement_distributions
h <- d[d$agg == "signal" & d$k == 1L & d$bag_population == "common", , drop = FALSE]
h <- h[, c("contrast", "method_a", "method_b", "route", "held",
           "n_mimic", "median_mimic", "p95_mimic", "anchor_boot0_mimic",
           "n_eicu", "median_eicu", "p95_eicu", "anchor_boot0_eicu")]
h <- h[order(h$contrast, h$method_a, h$method_b, h$held), ]
save_table(run, h, "headline_top1_signal_common", subdir = "tables")
cat("\n=== top-1 disagreement, signal level, common bags: MIMIC-IV against eICU ===\n\n")
print(h, row.names = FALSE, digits = 4)
p <- J$l3_shap_minus_llr_paired
p <- p[p$agg == "signal" & p$k %in% 1L, c("shap", "llr", "metric_oriented", "n_shared_bag_pairs_mimic",
                                          "frac_positive_mimic", "median_mimic",
                                          "n_shared_bag_pairs_eicu", "frac_positive_eicu", "median_eicu")]
cat("\n=== floor against floor, SHAP minus each LLR arm, signal level, top-1 ===\n\n")
print(p, row.names = FALSE, digits = 4)

.g <- gc(verbose = FALSE)
finalize_run(run, extra = list(mimic_run = basename(MIMIC), eicu_run = basename(EICU),
                               eicu_store = em$generator, internal_store = .eg$internal_store,
                               methods = paste(METHODS, collapse = ","),
                               r_max_used_gb = round(sum(.g[, ncol(.g)]) / 1024, 2)))
cat(sprintf("\nwritten: %s\n", run$path))
