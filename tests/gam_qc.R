# tests/gam_qc.R ---------------------------------------------------------------
# QC OF THE 64 FINAL LAYER-1 GAMs IN THE BUNDLE OF RECORD.
#
# Read docs/gam_qc_plan_20260908.md first. This script is five of that plan's
# six QC families in one pass over every fitted spec; the sixth, bootstrap
# stability, is tests/gam_qc_bootstrap.R and reads what this writes.
#
#   FIT             convergence read off the fREML optimiser (gradient norm,
#                   Hessian definiteness, smoothing parameters at a bound),
#                   EDF per term and its saturation ratio, mgcv::k.check on ALL
#                   rows with a pinned seed, and a k-doubling refit for every
#                   saturated s() term.
#   RESIDUALS       what a 0/1 outcome allows: randomised quantile residuals
#                   against N(0,1), binned residuals over fitted probability,
#                   leftover structure against every smooth covariate,
#                   calibration in the large, and a separation count.
#   DEPENDENCE      not recomputed. The concurvity numbers in diag_final are
#                   JOINED against the calibrated null tests/concurvity_null.R
#                   already saved, so every spec carries its excess over null.
#   GENERALISATION  the held-out deviance explained and log-loss of each spec
#                   from the OUT-OF-FOLD L it produced, beside its in-sample
#                   deviance explained (the optimism), its out-of-fold AUROC and
#                   AUPRC as a single column, the agreement between the final
#                   fit's in-sample L and the out-of-fold L on the same stays,
#                   and the fold-to-fold range of every diagnostic.
#   SURFACE         the EMPIRICAL SUPPORT MASK: for every smooth, a grid over
#                   the covariate's training distribution with the count of
#                   stays and of deaths in each bin, and a supported flag when
#                   both clear the floors declared in config/gam_qc.yml. For
#                   every ti() term, the same over the product grid.
#   INTERPRETABILITY the evidence-response curve of every smooth with its
#                   Vc-corrected band, summarised INSIDE the supported region
#                   only: range in nats, the share of the curve's variation that
#                   lies outside support, extrema count, the fraction of the
#                   supported grid resolved from zero, and the monotone
#                   direction against the design's declared expectation.
#
# THE OBJECT UNDER TEST IS THE BUNDLE. Each model's fitting frame is REBUILT
# through the same signal_frame() the apply path uses, from the bundle's frozen
# priors and the training ids, and re-attached as `$model`: a stripped bam
# object keeps its residuals, its response and its `Vc` but not the frame, and
# k.check needs the frame. Two guards hold the rebuilt frame to the fitted one
# -- the row count and the response must match the object bitwise, and when the
# targets store is present every frame is compared to the unstripped copy it
# holds. A guard that fails names the branch: the training ids, the frame
# builder or the bundle moved, and the script stops rather than QC-ing a frame
# the model was not fitted on.
#
# WHAT IS FITTED HERE, AND WHY THAT IS NOT A HARD RULE 8 BREACH. The k-doubling
# refits are the one fit in the script. They reach no score, no L and no bundle;
# each exists only to be compared, on a grid, with the frozen fit it doubles,
# and is discarded. Same status as the diagnostic smooth in R/04c.
#
# TWO MGCV FACTS THIS SCRIPT DEPENDS ON (docs/attribution_analysis_plan section
# 11). `unconditional = TRUE` is a no-op on a `bam`, so every band here comes
# from assigning `Vp <- Vc` before predict(). And a `bam(discrete = TRUE)` fit
# is in a model-scoped basis, so every comparison between two fits goes through
# predict(type = "terms") on a shared grid and never through coefficients.
#
# AGGREGATES ONLY (hard rule 1). The console gets one line per spec of fitted
# statistics and counts. The long tables hold grids, curve values, bin counts
# and event counts -- summaries over 40,000 stays, never a stay. The figures
# draw curves and bin-count histograms, never a rug of individual values.
#
#   Rscript tests/gam_qc.R                       # all 64 specs; MEASURED ~40 min
#   Rscript tests/gam_qc.R --plan-only           # what would run, no data read
#   Rscript tests/gam_qc.R --signals mbp,lactate # a subset
#   Rscript tests/gam_qc.R --models meas,full    # a subset of models
#   Rscript tests/gam_qc.R --no-refit --no-kcheck --no-figs   # the fast pass
#   Rscript tests/gam_qc.R --external out/runs/external_...   # join eICU arms
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(qs2)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
source("tests/gam_qc_common.R")   # helpers shared with tests/gam_qc_bootstrap.R, defined once

# --- arguments ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
ONLY_SG   <- .opt("--signals"); if (!is.null(ONLY_SG)) ONLY_SG <- strsplit(ONLY_SG, ",")[[1]]
ONLY_MD   <- .opt("--models");  if (!is.null(ONLY_MD)) ONLY_MD <- strsplit(ONLY_MD, ",")[[1]]
EXT_D     <- .opt("--external", NA_character_)
NO_REFIT  <- "--no-refit"  %in% args
NO_KCHECK <- "--no-kcheck" %in% args
NO_FIGS   <- "--no-figs"   %in% args
PLAN_ONLY <- "--plan-only" %in% args

# --- config: the QC plan, then the bundle of record ---------------------------
qcfg <- yaml::read_yaml("config/gam_qc.yml")
BUNDLE_P <- as.character(cfg_req(qcfg, "bundle"))
RUN_D    <- dirname(BUNDLE_P)
STORE    <- as.character(cfg_req(qcfg, "store"))
NULL_D   <- as.character(cfg_req(qcfg, "concurvity_null"))
G   <- cfg_req(qcfg, "grid");      S <- cfg_req(qcfg, "support")
KC  <- cfg_req(qcfg, "kcheck");    RS <- cfg_req(qcfg, "residuals")
SPB <- cfg_req(qcfg, "sp_bounds"); RK <- cfg_req(qcfg, "refit_k")
Z   <- as.numeric(cfg_req(qcfg, "band", "z"))
RHO_MONO <- as.numeric(cfg_req(qcfg, "interp", "monotone_rho_min"))
FLAT_NATS <- as.numeric(cfg_req(qcfg, "interp", "flat_range_nats"))
REPORT <- cfg_req(qcfg, "report")
for (k in c("frac_var_outside_high", "resolved_low")) cfg_req(REPORT, k)
for (k in c("n_points", "q_lo", "q_hi", "max_distinct", "n_points_2d")) cfg_req(G, k)
for (k in c("min_rows", "min_events", "min_rows_2d", "min_events_2d",
            "window_half_steps", "window_half_steps_2d")) cfg_req(S, k)
W1 <- as.integer(S$window_half_steps); W2 <- as.integer(S$window_half_steps_2d)
for (k in c("subsample", "n_rep", "seed")) cfg_req(KC, k)
for (k in c("n_bins", "n_bins_covariate", "z_out", "seed", "sep_eps")) cfg_req(RS, k)
for (k in c("huge", "tiny", "shrunk_edf")) cfg_req(SPB, k)
for (k in c("enabled", "multiplier")) cfg_req(RK, k)
DO_REFIT <- isTRUE(cfg_flag(RK, "enabled")) && !NO_REFIT

cfg_local <- load_config("config/config.yml")
bundle    <- load_bundle(BUNDLE_P, cfg = cfg_local, strict = TRUE)
cfg       <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
DIAG_THR  <- cfg_req(cfg, "diagnostics")     # the FROZEN thresholds, never restated
# Each read once here with no fallback: a `$` on an absent key would yield NULL
# and turn every comparison below into a zero-length logical that data.frame()
# rejects forty minutes in, with a message that names none of this.
for (k in c("k_index_min", "k_index_p_max", "edf_ratio_max", "concurvity_max", "dev_expl_min"))
  cfg_req(DIAG_THR, k, what = "the bundle's frozen diagnostics block must carry every threshold")

# The run directory that produced the bundle is the only place its out-of-fold
# and final diagnostics live; an absent table names a moved directory.
.run_tab <- function(sub, nm) {
  p <- file.path(RUN_D, sub, paste0(nm, ".rds"))
  if (!file.exists(p)) stop("gam_qc: the bundle's run directory has no ", sub, "/",
                            nm, ".rds (", RUN_D, "). The bundle and its run were ",
                            "separated; QC needs both.", call. = FALSE)
  readRDS(p)
}
DF <- .run_tab("diagnostics", "diag_final")
DO <- .run_tab("diagnostics", "diag_oof")

# --- the specs ---------------------------------------------------------------
keys <- names(bundle$models)
SPECS <- data.frame(key = keys,
                    signal = sub("/.*$", "", keys), model = sub("^.*/", "", keys),
                    stringsAsFactors = FALSE)
# The bundle carries FITTED specs only; an aliased or assigned cell has no
# model. Asserted, because a bundle that shipped an alias as a model would be
# QC-ed twice under two names.
bad <- SPECS[!is.na(mapply(spec_source, SPECS$signal, SPECS$model,
                          MoreArgs = list(cfg = cfg))), , drop = FALSE]
if (nrow(bad)) abort_values("gam_qc: the bundle carries a model for a cell spec_source() says is not fitted", bad$key)
if (!is.null(ONLY_SG)) SPECS <- SPECS[SPECS$signal %in% ONLY_SG, , drop = FALSE]
if (!is.null(ONLY_MD)) SPECS <- SPECS[SPECS$model  %in% ONLY_MD, , drop = FALSE]
if (!nrow(SPECS)) stop("gam_qc: no spec selected", call. = FALSE)

n_1d <- sum(vapply(SPECS$key, function(k) sum(vapply(bundle$models[[k]]$smooth,
              function(s) s$dim == 1L, logical(1))), integer(1)))
n_2d <- sum(vapply(SPECS$key, function(k) sum(vapply(bundle$models[[k]]$smooth,
              function(s) s$dim == 2L, logical(1))), integer(1)))
# Counted from the OBJECTS and over 1-D terms only, which is exactly what the
# refit block triggers on. `diag_final`'s `edf_ratio_max` also covers ti terms
# and would over-count the plan.
n_ref <- if (DO_REFIT) sum(vapply(SPECS$key, function(k) {
  b <- bundle$models[[k]]
  any(vapply(b$smooth, function(s) s$dim == 1L &&
        sum(b$edf[s$first.para:s$last.para]) / length(s$first.para:s$last.para) >
          DIAG_THR$edf_ratio_max, logical(1)))
}, logical(1))) else 0L
cat("\n=== plan ===\n")
cat(sprintf("  bundle                    : %s\n", BUNDLE_P))
cat(sprintf("  specs                     : %d of %d in the bundle\n", nrow(SPECS), length(keys)))
cat(sprintf("  1-D smooths / ti surfaces : %d / %d\n", n_1d, n_2d))
cat(sprintf("  k.check                   : %s\n", if (NO_KCHECK) "skipped" else
            sprintf("all rows (subsample %s), n_rep %d, seed %d",
                    if (as.integer(KC$subsample) == 0L) "off" else KC$subsample,
                    as.integer(KC$n_rep), as.integer(KC$seed))))
cat(sprintf("  k-doubling refits         : %s\n", if (!DO_REFIT) "skipped" else
            sprintf("%d spec(s) with edf/k' > %.2f (the bundle's edf_ratio_max)",
                    n_ref, DIAG_THR$edf_ratio_max)))
cat(sprintf("  support floors            : 1-D %d rows & %d deaths over +/-%d grid steps; 2-D %d & %d over +/-%d cells\n",
            S$min_rows, S$min_events, W1, S$min_rows_2d, S$min_events_2d, W2))
cat(sprintf("  grid                      : %d points over [q%.3f, q%.3f]; atoms if <= %d distinct\n",
            G$n_points, G$q_lo, G$q_hi, G$max_distinct))
if (PLAN_ONLY) quit(save = "no", status = 0L)

# --- data -------------------------------------------------------------------
tabs  <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
tr    <- folds$stay_id[folds$split == "train"]
y_of  <- function(ids) as.integer(tabs$cohort$mortality[match(as.character(ids),
                                                             as.character(tabs$cohort$stay_id))])

run <- new_run("gamqc", cfg_local, note = sprintf(
  "GAM QC of %d final spec(s) from %s", nrow(SPECS), basename(RUN_D)))
log_msg(run, "bundle: ", BUNDLE_P)
save_table(run, verify_bundle(bundle, cfg = cfg_local, strict = FALSE),
           "bundle_checks", subdir = "diagnostics")

# ==============================================================================
# 0. PROVENANCE: is this bundle the fit the store and diag_final describe?
# ==============================================================================
RAW <- NULL; store_coef <- NA_integer_; store_n <- NA_integer_
if (dir.exists(STORE) && requireNamespace("targets", quietly = TRUE)) {
  fm <- tryCatch(targets::tar_read_raw("final_models", store = STORE), error = function(e) NULL)
  if (!is.null(fm)) {
    store_n <- length(fm)
    store_coef <- sum(vapply(keys, function(k) !is.null(fm[[k]]) &&
                        identical(stats::coef(bundle$models[[k]]), stats::coef(fm[[k]])),
                        logical(1)))
    rm(fm)
  }
  RAW <- tryCatch(targets::tar_read_raw("layer1_final_raw", store = STORE)$models,
                  error = function(e) NULL)
}
cat("\n=== 0. provenance ===\n\n")
cat(sprintf("  models in bundle          : %d;  diag_final rows: %d;  diag_oof rows: %d\n",
            length(keys), nrow(DF), nrow(DO)))
if (is.na(store_coef)) {
  cat("  targets store             : absent or unreadable. Frame identity is\n")
  cat("                              held by row count and response only.\n")
} else {
  cat(sprintf("  store final_models        : %d models; coefficients identical to the bundle on %d of %d\n",
              store_n, store_coef, length(keys)))
  if (store_coef != length(keys)) {
    cat("  *** THE STORE AND THE BUNDLE DESCRIBE DIFFERENT FITS. The store was rebuilt\n")
    cat("      after this bundle was exported. Frame identity against the store is\n")
    cat("      therefore not evidence about this bundle and is NOT used below.\n")
    RAW <- NULL
  }
}

# ==============================================================================
# helpers
# ==============================================================================

# The grid, the neutral row and the window sum live in tests/gam_qc_common.R,
# shared with the bootstrap script. These are the config-bound aliases.
grid_of <- function(x) qc_grid_of(x, G$n_points, G$q_lo, G$q_hi, G$max_distinct)
grid_2d <- function(x) qc_grid_of(x, G$n_points_2d, G$q_lo, G$q_hi, G$max_distinct)

#' Count stays and deaths around each grid point. `n` and `events` are the
#' bin's own counts; `supported` is judged on the WINDOWED counts (see
#' config/gam_qc.yml support.window_half_steps). An atom grid has no window.
support_of <- function(x, y, gr, min_rows, min_events, h) {
  bin <- findInterval(x, gr$edges, rightmost.closed = TRUE)
  ok  <- bin >= 1L & bin <= length(gr$grid)
  n   <- tabulate(bin[ok], nbins = length(gr$grid))
  e   <- tabulate(bin[ok & y == 1L], nbins = length(gr$grid))
  hh  <- if (gr$kind == "atoms") 0L else h
  nw  <- qc_window_sum(n, hh); ew <- qc_window_sum(e, hh)
  data.frame(grid_x = gr$grid, n = n, events = e, n_window = nw, events_window = ew,
             supported = nw >= min_rows & ew >= min_events, stringsAsFactors = FALSE)
}

#' Which way the design expects a smooth to run. +1: risk rises with the
#' covariate; -1: risk falls; NA: no declared expectation. A curve running the
#' other way is REPORTED beside it and never corrected (CLAUDE.md: fit
#' unconstrained, report monotonicity violations).
#'
#' THE DECLARATION IS THE EXCURSION SIDE, NOT THE VARIABLE NAME (finding F36).
#' The first version expected every `value_min` / `q05` term to lower risk as
#' it rises, which is right for a low-side signal and wrong for a high-side
#' one: a rising MINIMUM lactate means the whole window sat high. And an
#' unpaired signal declares no side, so its two level terms carry no
#' expectation at all. Under the name rule the full pass reported 30 smooths
#' against; under the side rule 19, and the 11 that dropped were all of the
#' first kind. The propensity coordinate on the excursion side is expected to
#' raise risk; the other tail is undeclared. Intervention intensities are
#' expected to raise risk as a PROPENSITY reading (CLAUDE.md: sick patients
#' receive vasopressors), which is a general expectation rather than a per-
#' intervention declaration -- `expectation_source` says which kind each is.
expected_dir <- function(v, signal) {
  side <- excursion_side_of(signal, cfg)
  if (is.null(side) || is.na(side) || !nzchar(side)) side <- NA_character_
  # A LIST, not c(): c(1L, "x") is a character vector and "1" == 1 is TRUE in
  # R by coercion, which is the kind of comparison that works until it does not.
  if (grepl("__", v, fixed = TRUE)) return(list(dir = 1L, src = "propensity"))
  if (v == "trend") return(list(dir = NA_integer_, src = "none"))
  if (v == "pi_minus") return(list(dir = if (identical(side, "low"))  1L else NA_integer_,
                                   src = if (identical(side, "low"))  "excursion_side" else "none"))
  if (v == "pi_plus")  return(list(dir = if (identical(side, "high")) 1L else NA_integer_,
                                   src = if (identical(side, "high")) "excursion_side" else "none"))
  base <- if (grepl("_delta$", v)) delta_base_of(v) else v
  if (base %in% c("q05", "value_min", "q95", "value_max", "value_median")) {
    return(list(dir = if (identical(side, "low")) -1L else if (identical(side, "high")) 1L else NA_integer_,
                src = if (is.na(side)) "none" else "excursion_side"))
  }
  list(dir = NA_integer_, src = "none")
}

#' Summaries of one curve, computed INSIDE the supported region.
curve_summary <- function(f, se, x, sup) {
  out <- list(n_grid = length(f), n_supported = sum(sup),
              range_all = max(f) - min(f), range_in = NA_real_,
              frac_var_outside = NA_real_, frac_curv_outside = NA_real_,
              n_extrema_in = NA_integer_, slope_rho_in = NA_real_,
              resolved_frac_in = NA_real_, se_med_in = NA_real_, se_med_out = NA_real_,
              max_abs_in = NA_real_, max_abs_out = NA_real_)
  if (sum(sup) < 3L) return(out)
  i <- which(sup); fi <- f[i]; xi <- x[i]
  out$range_in <- max(fi) - min(fi)
  d1 <- abs(diff(f)); step_out <- !(sup[-1] & sup[-length(sup)])
  out$frac_var_outside <- if (sum(d1) > 0) sum(d1[step_out]) / sum(d1) else NA_real_
  if (length(f) >= 3L) {
    d2 <- abs(diff(f, differences = 2L))
    n <- length(sup); mid_out <- !(sup[1:(n - 2)] & sup[2:(n - 1)] & sup[3:n])
    out$frac_curv_outside <- if (sum(d2) > 0) sum(d2[mid_out]) / sum(d2) else NA_real_
  }
  out$n_extrema_in <- qc_n_extrema(fi)
  out$slope_rho_in <- if (stats::sd(fi) > 0 && stats::sd(xi) > 0)
    stats::cor(xi, fi, method = "spearman") else NA_real_
  out$resolved_frac_in <- mean(abs(fi) > Z * se[i])
  out$se_med_in  <- stats::median(se[i])
  out$se_med_out <- if (any(!sup)) stats::median(se[!sup]) else NA_real_
  out$max_abs_in  <- max(abs(fi))
  out$max_abs_out <- if (any(!sup)) max(abs(f[!sup])) else NA_real_
  out
}

#' Randomised quantile residuals for a 0/1 outcome (Dunn and Smyth).
rqr_of <- function(y, p, seed) {
  with_seed(seed, {
    u  <- stats::runif(length(y))
    lo <- ifelse(y == 1L, 1 - p, 0); hi <- ifelse(y == 1L, 1, 1 - p)
    stats::qnorm(pmin(pmax(lo + u * (hi - lo), 1e-12), 1 - 1e-12))
  })
}

#' Equal-count bins of fitted probability, TIE-AWARE: observed against expected.
#'
#' Rows with the identical fitted probability always share a bin (finding
#' F37). The first version split them by `order()`, which is a stable sort and
#' therefore splits a tie block by ROW ORDER -- and the frame's row order is
#' correlated with outcome, so the sub-bins of one calibrated block showed
#' observed rates from 0.04 to 0.21 around a common fitted 0.099, and
#' `creatinine/intv` reported a Hosmer-Lemeshow of 1156 that fell to 23 once
#' the block was kept whole. Bins that a tie block swallows are collapsed, so
#' a model with a large point mass has fewer than `nb` bins and says so.
binned_of <- function(y, p, nb) {
  r <- rank(round(p, 10), ties.method = "min")
  g <- ceiling(r * nb / length(p)); g <- match(g, sort(unique(g)))
  do.call(rbind, lapply(sort(unique(g)), function(b) {
    i <- which(g == b); n <- length(i); ex <- mean(p[i]); ob <- mean(y[i])
    se <- sqrt(ex * (1 - ex) / n)
    data.frame(bin = b, n = n, p_mean = round(ex, 5), obs_rate = round(ob, 5),
               gap = round(ob - ex, 5), z = round((ob - ex) / se, 3), stringsAsFactors = FALSE)
  }))
}

#' Double k on the named s() terms of a formula. String surgery on the emitted
#' term, so the refit differs from the frozen fit in exactly one number.
formula_with_k <- function(f, vars, mult) {
  lab <- attr(stats::terms(f), "term.labels")
  for (v in vars) {
    j <- grep(sprintf("^s\\(%s,", v), lab)
    if (length(j) != 1L) next
    a <- .call_args(lab[j]); k_new <- as.integer(a$k) * mult
    lab[j] <- sub("k = [0-9]+", sprintf("k = %d", k_new), lab[j])
  }
  stats::as.formula(paste(all.vars(f)[1], "~", paste(lab, collapse = " + ")), env = globalenv())
}

#' Draw every 1-D curve of one spec: band, unsupported shading, bin counts.
plot_curves <- function(run, key, cv, tm) {
  terms <- unique(cv$term); n <- length(terms)
  if (!n) return(invisible(NULL))
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
  save_fig(run, paste0("curves_", gsub("/", "_", key)), width = 4 * nc, height = 3.2 * nr, dpi = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2.5, 0.8), mgp = c(2.1, 0.7, 0))
  for (t in terms) {
    z <- cv[cv$term == t, ]; z <- z[order(z$grid_x), ]
    ti <- tm[tm$term == t, ][1, ]
    lo <- z$fit - Z * z$se; hi <- z$fit + Z * z$se
    yl <- range(c(lo, hi), na.rm = TRUE); yl[1] <- yl[1] - 0.22 * diff(yl)
    graphics::plot(z$grid_x, z$fit, type = "n", ylim = yl, xlab = z$variable[1],
                   ylab = "partial effect (nats)",
                   main = sprintf("%s  edf %.1f  k-idx %.2f", t, ti$edf, ti$k_index),
                   cex.main = 0.85)
    # unsupported regions, shaded
    r <- rle(!z$supported); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    for (i in which(r$values)) {
      x0 <- if (starts[i] > 1L) (z$grid_x[starts[i] - 1L] + z$grid_x[starts[i]]) / 2 else graphics::par("usr")[1]
      x1 <- if (ends[i] < nrow(z)) (z$grid_x[ends[i]] + z$grid_x[ends[i] + 1L]) / 2 else graphics::par("usr")[2]
      graphics::rect(x0, yl[1], x1, yl[2] + diff(yl), col = "#00000014", border = NA)
    }
    # bin counts along the bottom (counts, never rows)
    h <- z$n / max(z$n, 1) * 0.18 * diff(yl)
    w <- if (nrow(z) > 1L) min(diff(z$grid_x)) * 0.9 else 1
    graphics::rect(z$grid_x - w / 2, yl[1], z$grid_x + w / 2, yl[1] + h, col = "grey75", border = NA)
    graphics::polygon(c(z$grid_x, rev(z$grid_x)), c(lo, rev(hi)), col = "#3B6FB620", border = NA)
    if (z$kind[1] == "atoms") {
      graphics::points(z$grid_x, z$fit, pch = 19, cex = 0.7, col = "#1F3F7A")
      graphics::segments(z$grid_x, lo, z$grid_x, hi, col = "#1F3F7A")
    } else graphics::lines(z$grid_x, z$fit, lwd = 2, col = "#1F3F7A")
    graphics::abline(h = 0, lty = 3, col = "grey40")
  }
  invisible(NULL)
}

#' Draw every ti surface of one spec, masked outside support.
plot_surfaces <- function(run, key, sf) {
  terms <- unique(sf$term); n <- length(terms)
  if (!n) return(invisible(NULL))
  nc <- min(4L, n); nr <- ceiling(n / nc)
  save_fig(run, paste0("surface_", gsub("/", "_", key)), width = 4 * nc, height = 3.6 * nr, dpi = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2.5, 0.8), mgp = c(2.1, 0.7, 0))
  for (t in terms) {
    z <- sf[sf$term == t, ]
    gx <- sort(unique(z$x1)); gy <- sort(unique(z$x2))
    M <- matrix(NA_real_, length(gx), length(gy))
    M[cbind(match(z$x1, gx), match(z$x2, gy))] <- ifelse(z$supported, z$fit, NA_real_)
    if (all(is.na(M))) { graphics::plot.new(); graphics::title(paste(t, "(no supported cell)")); next }
    graphics::image(gx, gy, M, col = grDevices::hcl.colors(24, "RdBu", rev = TRUE),
                    xlab = z$var1[1], ylab = z$var2[1], main = t, cex.main = 0.8)
    if (length(gx) > 2L && length(gy) > 2L)
      try(graphics::contour(gx, gy, M, add = TRUE, col = "grey20", lwd = 0.6), silent = TRUE)
  }
  invisible(NULL)
}

# ==============================================================================
# THE PASS OVER EVERY SPEC
# ==============================================================================
cat("\n=== per spec ===\n")
cat("  conv: converged | |g|: max |gradient| of the fREML objective | H+: Hessian positive-definite\n")
cat("  edf: total | sat: max edf/k' (term) | k-idx: min k-index (p) on all rows\n")
cat("  ks: randomised-quantile-residual KS vs N(0,1) | bins: binned-residual bins out of 20\n")
cat("  dir-: smooths running against the declared expectation | out: max share of curve variation outside support\n\n")

FIT <- list(); TERMS <- list(); RESID <- list(); BINS <- list(); STRUCT <- list()
CUR1 <- list(); CSUM <- list(); SUR2 <- list(); SSUM <- list()
REFIT <- list(); FRAMES <- list()
FRAME_CACHE <- new.env(parent = emptyenv())   # one frame per (signal, column set); see qc_frame_for()
t_all <- start_timer()

for (.i in seq_len(nrow(SPECS))) {
  key <- SPECS$key[.i]; sg <- SPECS$signal[.i]; md <- SPECS$model[.i]
  b   <- bundle$models[[key]]
  pri <- priors_for(bundle$priors, sg, "final")
  d   <- qc_frame_for(FRAME_CACHE, sg, md, tabs, cfg, pri, stay_ids = tr)
  f   <- attr(d, "formula")
  y   <- as.integer(d$mortality)

  # --- the frame identity guards -----------------------------------------
  if (nrow(d) != length(b$y) || !identical(as.integer(b$y), y)) {
    stop(sprintf("gam_qc [%s]: the rebuilt frame (%d rows) does not match the fitted response (%d rows%s). ",
                 key, nrow(d), length(b$y),
                 if (nrow(d) == length(b$y)) ", response differs" else ""),
         "One of: the training ids, the measured-subset rule, or the bundle moved. ",
         "Name the branch before touching anything.", call. = FALSE)
  }
  frame_same <- NA
  if (!is.null(RAW) && !is.null(RAW[[key]]$model)) {
    m0 <- RAW[[key]]$model; cols <- intersect(names(m0), names(d))
    frame_same <- isTRUE(all.equal(as.data.frame(m0[, cols, drop = FALSE]),
                                   as.data.frame(d[, cols, drop = FALSE]),
                                   check.attributes = FALSE, tolerance = 0))
    if (!frame_same) {
      stop(sprintf("gam_qc [%s]: the rebuilt frame differs from the frame the store fitted on, ",
                   "while the coefficients are identical. The frame builder or a prior moved.", key),
           call. = FALSE)
    }
    RAW[[key]] <- NULL   # 86 MB of store copies; each is needed exactly once
  }
  vars <- setdiff(names(d), c("stay_id", "mortality"))
  b$model <- as.data.frame(d[, c("mortality", vars), drop = FALSE])
  bV <- b; if (!is.null(bV$Vc)) bV$Vp <- bV$Vc
  p_hat <- as.numeric(b$fitted.values); eta <- as.numeric(b$linear.predictors)

  # --- A. fit ---------------------------------------------------------------
  grad <- b$outer.info$grad; hess <- b$outer.info$hess
  grad_norm <- if (length(grad)) max(abs(grad)) else NA_real_
  # The fREML Hessian over log smoothing parameters. A term shrunk to nothing
  # sits at a huge sp where the objective is FLAT, so its eigenvalue is ~0 and
  # can come out marginally negative from round-off; that is not a saddle. PD
  # is therefore judged against the largest eigenvalue, and the flat directions
  # are counted separately.
  hev <- if (length(hess)) eigen(hess, symmetric = TRUE, only.values = TRUE)$values else numeric(0)
  hess_min <- if (length(hev)) min(hev) else NA_real_
  hess_pd  <- if (length(hev)) min(hev) > -1e-6 * max(abs(hev)) else NA
  n_flat   <- if (length(hev)) sum(abs(hev) < 1e-6 * max(abs(hev))) else NA_integer_
  sp <- b$sp
  trows <- lapply(seq_along(b$smooth), function(j) {
    s <- b$smooth[[j]]; idx <- s$first.para:s$last.para
    e <- sum(b$edf[idx]); e1 <- sum(b$edf1[idx])
    # mgcv's smooth objects are NOT one schema. `bs.dim` is a field of a
    # univariate `mgcv.smooth` and does not exist on a `tensor.smooth`, whose
    # basis sizes live on its margins (finding F24 in the plan: the first
    # version read `s$bs.dim` on every smooth and data.frame() refused the
    # zero-length column on the first ti term).
    bsd <- if (!is.null(s$bs.dim)) as.integer(s$bs.dim) else
      as.integer(prod(vapply(s$margin, function(m) m$bs.dim, numeric(1))))
    data.frame(key = key, signal = sg, model = md, term = s$label, dim = s$dim,
               variable = paste(s$term, collapse = " x "),
               k_prime = length(idx), bs_dim = bsd, edf = round(e, 3), edf1 = round(e1, 3),
               edf_ratio = round(e / length(idx), 4),
               shrunk_out = e < as.numeric(SPB$shrunk_edf),
               saturated = e / length(idx) > DIAG_THR$edf_ratio_max,
               stringsAsFactors = FALSE)
  })
  TM <- do.call(rbind, trows)
  # THE SCHEMA IS FIXED BEFORE ANY BRANCH. The k.check block below fills these
  # on success and leaves them on failure or under --no-kcheck; without the
  # initialisation one failed k.check would give that spec's term table fewer
  # columns and rbind() would refuse the whole table at the end of the run.
  TM$k_index <- NA_real_; TM$k_p <- NA_real_; TM$k_prime_mgcv <- NA_integer_
  # smoothing parameters: one per penalty; a ti carries two. Reported per model.
  fit_row <- data.frame(
    key = key, signal = sg, model = md, n_rows = nrow(d), n_events = sum(y),
    n_coef = length(stats::coef(b)), n_smooth = length(b$smooth), rank = b$rank,
    converged = isTRUE(b$converged), boundary = isTRUE(b$boundary),
    # `bam` records its iteration count at the top level; `outer.info` carries
    # only the gradient and Hessian for fREML.
    n_iter = if (!is.null(b$iter)) as.integer(b$iter) else NA_integer_,
    grad_norm = signif(grad_norm, 3), hess_min_eig = signif(hess_min, 3),
    hess_pd = hess_pd, n_hess_flat = n_flat,
    n_sp = length(sp), sp_min = signif(min(sp), 3), sp_max = signif(max(sp), 3),
    n_sp_huge = sum(sp > as.numeric(SPB$huge)), n_sp_tiny = sum(sp < as.numeric(SPB$tiny)),
    edf_total = round(sum(b$edf), 3), edf1_total = round(sum(b$edf1), 3),
    n_terms_shrunk_out = sum(TM$shrunk_out), n_terms_saturated = sum(TM$saturated),
    edf_ratio_max = max(TM$edf_ratio), edf_worst = TM$term[which.max(TM$edf_ratio)],
    dev_expl = round((b$null.deviance - b$deviance) / b$null.deviance, 5),
    aic = round(b$aic, 2), frame_identical_to_store = frame_same,
    stringsAsFactors = FALSE)

  # --- B. k.check on all rows, seeded -------------------------------------
  kc_min <- NA_real_; kc_p <- NA_real_; kc_worst <- NA_character_
  fit_row$kcheck_note <- if (NO_KCHECK) "skipped (--no-kcheck)" else ""
  if (!NO_KCHECK) {
    sub <- as.integer(KC$subsample); if (sub == 0L) sub <- nrow(d) + 1L
    kc <- with_seed(as.integer(KC$seed),
                    tryCatch(mgcv::k.check(b, subsample = sub, n.rep = as.integer(KC$n_rep)),
                             error = function(e) e))
    if (inherits(kc, "error")) {
      fit_row$kcheck_note <- conditionMessage(kc)
    } else {
      m <- match(TM$term, rownames(kc))
      TM$k_index <- round(unname(kc[m, "k-index"]), 4); TM$k_p <- round(unname(kc[m, "p-value"]), 4)
      TM$k_prime_mgcv <- as.integer(unname(kc[m, "k'"]))
      j <- which.min(TM$k_index)
      kc_min <- TM$k_index[j]; kc_p <- TM$k_p[j]; kc_worst <- TM$term[j]
    }
  }
  # `diag_final` divides EDF by k.check's k'; this table divides by the
  # coefficient count. They agree for every s() term (both 9 at k = 10) and the
  # count of terms where they do not is reported rather than assumed zero.
  TM$k_prime_agrees <- is.na(TM$k_prime_mgcv) | TM$k_prime_mgcv == TM$k_prime
  fit_row$k_index_min <- kc_min; fit_row$k_index_p <- kc_p; fit_row$k_worst <- kc_worst
  fit_row$k_flag <- isTRUE(!is.na(kc_min) && !is.na(kc_p) &&
                           kc_min < DIAG_THR$k_index_min && kc_p < DIAG_THR$k_index_p_max)
  fit_row$n_kprime_disagree <- sum(!TM$k_prime_agrees)
  r0 <- DF[DF$signal == sg & DF$model == md, , drop = FALSE]
  fit_row$k_index_recorded <- if (nrow(r0)) r0$k_index_min[1] else NA_real_
  fit_row$k_worst_recorded <- if (nrow(r0)) r0$k_worst[1] else NA_character_
  fit_row$dev_expl_recorded <- if (nrow(r0)) r0$dev_expl[1] else NA_real_

  # --- C. residual sanity -------------------------------------------------
  zq <- rqr_of(y, p_hat, as.integer(RS$seed))
  ks <- suppressWarnings(stats::ks.test(zq, "pnorm"))
  m2 <- mean(zq^2); m3 <- mean(zq^3); m4 <- mean(zq^4)
  bn <- binned_of(y, p_hat, as.integer(RS$n_bins))
  hl <- {
    O <- bn$obs_rate * bn$n; E <- bn$p_mean * bn$n
    sum((O - E)^2 / (E * (1 - bn$p_mean)))
  }
  # THE DOMINANT TIE BLOCK, measured once as itself. Every stay sitting on an
  # intervention's point mass at zero gets the SAME fitted probability, and in
  # an `intv` model that is most of the cohort. Equal-count bins split that one
  # block across many bins, so a single miscalibration at the atom shows up as
  # many bins out (first seen on mbp/intv: 11 of 20 bins) and reads as
  # widespread. The block's own observed-minus-expected is the honest number.
  tb  <- table(round(p_hat, 10))
  blk <- round(p_hat, 10) == as.numeric(names(tb)[which.max(tb)])
  # `blk_frac`, not `atom_frac`: the 1-D loop below has an `atom_frac` of its
  # own (a covariate's point mass) and the residual row must not be able to
  # pick that one up if the two blocks are ever reordered.
  blk_frac <- mean(blk); atom_z <- NA_real_; atom_obs <- NA_real_; atom_exp <- NA_real_
  if (sum(blk) >= as.integer(S$min_rows)) {
    atom_obs <- mean(y[blk]); atom_exp <- mean(p_hat[blk])
    atom_z <- (atom_obs - atom_exp) / sqrt(atom_exp * (1 - atom_exp) / sum(blk))
  }
  # leftover structure: the mean RQR in equal-count bins of each smooth covariate,
  # standardised by sqrt(n_bin). k.check asks the same question with neighbour
  # differences; this asks it on the scale a reader can picture.
  st <- do.call(rbind, lapply(unique(unlist(lapply(b$smooth, `[[`, "term"))), function(v) {
    x <- d[[v]]; nbv <- as.integer(RS$n_bins_covariate)
    o <- order(x); grp <- ceiling(seq_along(o) * nbv / length(o))
    mm <- vapply(seq_len(nbv), function(g) mean(zq[o[grp == g]]), numeric(1))
    nn <- vapply(seq_len(nbv), function(g) sum(grp == g), numeric(1))
    data.frame(key = key, signal = sg, model = md, variable = v,
               # the mean RQR per covariate bin, in SD units (the effect size),
               # and the same scaled by sqrt(n) (the test statistic). At 40,000
               # rows the second is large for tiny effects; read the first.
               max_abs_binned_mean = round(max(abs(mm)), 4),
               max_abs_binned_z = round(max(abs(mm) * sqrt(nn)), 3),
               spearman_z_x = round(stats::cor(x, zq, method = "spearman"), 4),
               stringsAsFactors = FALSE)
  }))
  # In-sample calibration slope of y on the fitted log-odds. Penalisation
  # shrinks the linear predictor, so a slope above 1 in-sample is the
  # signature of shrinkage (under-confidence), not of a fitting failure.
  cal <- llr_calibration(eta - logit(pri$p_bar), y, pri$p_bar)
  eps <- as.numeric(RS$sep_eps)
  res_row <- data.frame(
    key = key, signal = sg, model = md,
    rqr_ks = round(unname(ks$statistic), 4), rqr_ks_p = signif(ks$p.value, 3),
    rqr_var = round(m2, 4), rqr_skew = round(m3 / m2^1.5, 4), rqr_kurt = round(m4 / m2^2 - 3, 4),
    rqr_frac_abs_gt3 = round(mean(abs(zq) > 3), 5),
    cal_in_large = signif(mean(p_hat) - mean(y), 3),
    n_bins = nrow(bn), bins_out = sum(abs(bn$z) > as.numeric(RS$z_out)), bins_max_abs_z = max(abs(bn$z)),
    # the worst gap in PROBABILITY units over bins of at least min_rows: the
    # effect size the z and the Hosmer-Lemeshow statistic scale with n
    bins_max_abs_gap = round(max(abs(bn$gap[bn$n >= as.integer(S$min_rows)])), 5),
    hosmer_lemeshow = round(hl, 2),
    atom_block_frac = round(blk_frac, 4), atom_block_obs = round(atom_obs, 5),
    atom_block_exp = round(atom_exp, 5), atom_block_z = round(atom_z, 3),
    struct_max_abs_mean = max(st$max_abs_binned_mean),
    struct_max_abs_z = max(st$max_abs_binned_z), struct_worst = st$variable[which.max(st$max_abs_binned_z)],
    cal_slope_in = cal$slope, cal_slope_in_lo = cal$slope_lo, cal_slope_in_hi = cal$slope_hi,
    frac_separated = round(mean(p_hat < eps | p_hat > 1 - eps), 5),
    eta_abs_max = round(max(abs(eta)), 3),
    pearson_dispersion = round(sum(stats::residuals(b, type = "pearson")^2) / b$df.residual, 4),
    stringsAsFactors = FALSE)
  bn$key <- key; bn$signal <- sg; bn$model <- md

  # --- D. support and curves, 1-D ------------------------------------------
  nr0 <- qc_neutral_row(d, vars)
  cv <- list(); cs <- list()
  for (j in seq_along(b$smooth)) {
    s <- b$smooth[[j]]; if (s$dim != 1L) next
    v <- s$term[1]; x <- d[[v]]
    gr <- grid_of(x)
    su <- support_of(x, y, gr, as.integer(S$min_rows), as.integer(S$min_events), W1)
    nd <- nr0[rep(1L, length(gr$grid)), , drop = FALSE]; nd[[v]] <- gr$grid
    pr <- stats::predict(bV, newdata = nd, type = "terms", se.fit = TRUE, discrete = FALSE)
    jj <- match(s$label, colnames(pr$fit))
    if (is.na(jj)) stop("gam_qc [", key, "]: predict(type = 'terms') has no column ", s$label, call. = FALSE)
    fit <- as.numeric(pr$fit[, jj]); se <- as.numeric(pr$se.fit[, jj])
    tab <- table(x); atom_val <- as.numeric(names(tab)[which.max(tab)]); atom_frac <- max(tab) / length(x)
    cv[[length(cv) + 1L]] <- data.frame(
      key = key, signal = sg, model = md, term = s$label, variable = v, kind = gr$kind,
      grid_x = gr$grid, fit = round(fit, 6), se = round(se, 6),
      n = su$n, events = su$events, n_window = su$n_window, events_window = su$events_window,
      supported = su$supported, stringsAsFactors = FALSE)
    cm <- curve_summary(fit, se, gr$grid, su$supported)
    ex <- expected_dir(v, sg); ed <- ex$dir
    # A term shrunk to nothing is a flat line with round-off for a slope, so
    # its direction is not a finding either way; the same for a curve whose
    # supported range is below a thousandth of a nat.
    shrunk <- isTRUE(TM$shrunk_out[TM$term == s$label][1])
    flat   <- is.na(cm$range_in) || cm$range_in < FLAT_NATS
    if (flat) cm$n_extrema_in <- NA_integer_    # round-off wobble is not an extremum
    cs[[length(cs) + 1L]] <- cbind(
      data.frame(key = key, signal = sg, model = md, term = s$label, variable = v, kind = gr$kind,
                 n_distinct = length(tab), atom_value = atom_val, atom_frac = round(atom_frac, 4),
                 frac_rows_in_grid = round(sum(su$n) / length(x), 4),
                 frac_rows_supported = round(sum(su$n[su$supported]) / length(x), 4),
                 shrunk_out = shrunk, stringsAsFactors = FALSE),
      as.data.frame(lapply(cm, function(z) if (is.numeric(z)) round(z, 5) else z)),
      data.frame(expected_dir = ed, expectation_source = ex$src,
                 # SIX STATES, because "runs against the declared direction"
                 # is one of several things a curve can do and the others are
                 # not failures: a U-shape has no direction, a shrunk-out term
                 # has no curve, and a term with no declared expectation
                 # cannot contradict one.
                 direction = if (shrunk) "shrunk_out" else if (flat) "flat"
                             else if (is.na(cm$slope_rho_in)) "flat"
                             else if (abs(cm$slope_rho_in) < RHO_MONO) "non_monotone"
                             else if (is.na(ed)) "undeclared"
                             else if (sign(cm$slope_rho_in) == ed) "with" else "against",
                 dir_ok = if (is.na(ed) || is.na(cm$slope_rho_in) || shrunk || flat ||
                              abs(cm$slope_rho_in) < RHO_MONO) NA
                          else sign(cm$slope_rho_in) == ed,
                 stringsAsFactors = FALSE))
  }
  CV <- if (length(cv)) do.call(rbind, cv) else NULL
  CS <- if (length(cs)) do.call(rbind, cs) else NULL

  # --- E. support and surfaces, 2-D (the ti terms) --------------------------
  sf <- list(); ss <- list()
  for (j in seq_along(b$smooth)) {
    s <- b$smooth[[j]]; if (s$dim != 2L) next
    v1 <- s$term[1]; v2 <- s$term[2]; x1 <- d[[v1]]; x2 <- d[[v2]]
    gr1 <- grid_2d(x1); gr2 <- grid_2d(x2)
    b1 <- findInterval(x1, gr1$edges, rightmost.closed = TRUE); b2 <- findInterval(x2, gr2$edges, rightmost.closed = TRUE)
    ok <- b1 >= 1L & b1 <= length(gr1$grid) & b2 >= 1L & b2 <= length(gr2$grid)
    cell <- (b2 - 1L) * length(gr1$grid) + b1
    nC <- length(gr1$grid) * length(gr2$grid)
    n  <- tabulate(cell[ok], nbins = nC); e <- tabulate(cell[ok & y == 1L], nbins = nC)
    # windowed over +/- W2 cells in each margin (no window along an atom margin)
    h1 <- if (gr1$kind == "atoms") 0L else W2; h2 <- if (gr2$kind == "atoms") 0L else W2
    win2 <- function(v) {
      n1 <- length(gr1$grid); n2 <- length(gr2$grid)
      M <- matrix(v, n1, n2)
      # apply() drops to a vector when a margin has one point; the explicit
      # matrix() puts the dimensions back before the second pass.
      M <- matrix(apply(M, 2, qc_window_sum, h = h1), n1, n2)
      M <- t(matrix(apply(M, 1, qc_window_sum, h = h2), n2, n1))
      as.numeric(M)
    }
    nw <- win2(n); ew <- win2(e)
    supp <- nw >= as.integer(S$min_rows_2d) & ew >= as.integer(S$min_events_2d)
    ndg <- expand.grid(x1 = gr1$grid, x2 = gr2$grid)      # x1 varies fastest, matching `cell`
    nd <- nr0[rep(1L, nrow(ndg)), , drop = FALSE]; nd[[v1]] <- ndg$x1; nd[[v2]] <- ndg$x2
    pr <- stats::predict(bV, newdata = nd, type = "terms", se.fit = TRUE, discrete = FALSE)
    jj <- match(s$label, colnames(pr$fit))
    if (is.na(jj)) stop("gam_qc [", key, "]: predict(type = 'terms') has no column ", s$label, call. = FALSE)
    fit <- as.numeric(pr$fit[, jj]); se <- as.numeric(pr$se.fit[, jj])
    sf[[length(sf) + 1L]] <- data.frame(
      key = key, signal = sg, model = md, term = s$label, var1 = v1, var2 = v2,
      x1 = ndg$x1, x2 = ndg$x2, fit = round(fit, 6), se = round(se, 6),
      n = n, events = e, supported = supp, stringsAsFactors = FALSE)
    fi <- fit[supp]
    ss[[length(ss) + 1L]] <- data.frame(
      key = key, signal = sg, model = md, term = s$label, var1 = v1, var2 = v2,
      n_cells = nC, n_supported = sum(supp), frac_cells_supported = round(mean(supp), 4),
      frac_rows_supported = round(sum(n[supp]) / length(x1), 4),
      range_in = if (sum(supp) >= 2L) round(max(fi) - min(fi), 5) else NA_real_,
      range_all = round(max(fit) - min(fit), 5),
      sd_in = if (sum(supp) >= 2L) round(stats::sd(fi), 5) else NA_real_,
      sd_out = if (sum(!supp) >= 2L) round(stats::sd(fit[!supp]), 5) else NA_real_,
      max_abs_in = if (sum(supp)) round(max(abs(fi)), 5) else NA_real_,
      max_abs_out = if (sum(!supp)) round(max(abs(fit[!supp])), 5) else NA_real_,
      resolved_frac_in = if (sum(supp)) round(mean(abs(fi) > Z * se[supp]), 4) else NA_real_,
      se_med_in = if (sum(supp)) round(stats::median(se[supp]), 5) else NA_real_,
      stringsAsFactors = FALSE)
  }
  SF <- if (length(sf)) do.call(rbind, sf) else NULL
  SS <- if (length(ss)) do.call(rbind, ss) else NULL

  # --- F. the k-doubling refit for saturated s() terms ----------------------
  if (DO_REFIT && any(TM$saturated & TM$dim == 1L)) {
    sat <- TM$variable[TM$saturated & TM$dim == 1L]
    mult <- as.integer(RK$multiplier)
    for (v in sat) {
      nd_ <- length(unique(d[[v]]))
      k_old <- TM$bs_dim[TM$variable == v & TM$dim == 1L][1]
      k_new <- k_old * mult
      if (k_new >= nd_) {
        REFIT[[length(REFIT) + 1L]] <- data.frame(
          key = key, signal = sg, model = md, variable = v, k_old = k_old, k_new = k_new,
          n_distinct = nd_, status = "skipped: k_new >= distinct values (saturation is structural)",
          dev_expl_old = fit_row$dev_expl, dev_expl_new = NA_real_, edf_old = TM$edf[TM$variable == v][1],
          edf_new = NA_real_, edf_ratio_new = NA_real_, k_index_new = NA_real_,
          shape_sd_in = NA_real_, shape_max_abs_in = NA_real_, curve_cor_in = NA_real_,
          stringsAsFactors = FALSE)
        next
      }
      f2 <- formula_with_k(f, v, mult)
      b2 <- try(.bam_fit(f2, d, cfg), silent = TRUE)
      if (inherits(b2, "try-error")) {
        REFIT[[length(REFIT) + 1L]] <- data.frame(
          key = key, signal = sg, model = md, variable = v, k_old = k_old, k_new = k_new,
          n_distinct = nd_, status = paste0("refit failed: ", substr(conditionMessage(attr(b2, "condition")), 1, 120)),
          dev_expl_old = fit_row$dev_expl, dev_expl_new = NA_real_, edf_old = TM$edf[TM$variable == v][1],
          edf_new = NA_real_, edf_ratio_new = NA_real_, k_index_new = NA_real_,
          shape_sd_in = NA_real_, shape_max_abs_in = NA_real_, curve_cor_in = NA_real_,
          stringsAsFactors = FALSE)
        next
      }
      s2 <- b2$smooth[[which(vapply(b2$smooth, function(s) identical(s$term[1], v) && s$dim == 1L, logical(1)))[1]]]
      idx2 <- s2$first.para:s2$last.para; edf2 <- sum(b2$edf[idx2])
      kc2 <- NA_real_
      if (!NO_KCHECK) {
        k2 <- with_seed(as.integer(KC$seed), tryCatch(mgcv::k.check(b2, subsample = nrow(d) + 1L,
                                                                    n.rep = as.integer(KC$n_rep)),
                                                       error = function(e) NULL))
        if (!is.null(k2) && s2$label %in% rownames(k2)) kc2 <- round(unname(k2[s2$label, "k-index"]), 4)
      }
      z0 <- CV[CV$variable == v, ]
      b2V <- b2; if (!is.null(b2V$Vc)) b2V$Vp <- b2V$Vc
      nd <- nr0[rep(1L, nrow(z0)), , drop = FALSE]; nd[[v]] <- z0$grid_x
      p2 <- stats::predict(b2V, newdata = nd, type = "terms", discrete = FALSE)
      f_new <- as.numeric(p2[, match(s2$label, colnames(p2))])
      i <- which(z0$supported)
      dif <- (f_new - z0$fit)[i]; dif <- dif - mean(dif)
      REFIT[[length(REFIT) + 1L]] <- data.frame(
        key = key, signal = sg, model = md, variable = v, k_old = k_old, k_new = k_new,
        n_distinct = nd_, status = "ok",
        dev_expl_old = fit_row$dev_expl,
        dev_expl_new = round((b2$null.deviance - b2$deviance) / b2$null.deviance, 5),
        edf_old = TM$edf[TM$variable == v][1], edf_new = round(edf2, 3),
        edf_ratio_new = round(edf2 / length(idx2), 4), k_index_new = kc2,
        shape_sd_in = if (length(i) >= 3L) round(stats::sd(dif), 5) else NA_real_,
        shape_max_abs_in = if (length(i)) round(max(abs(dif)), 5) else NA_real_,
        curve_cor_in = if (length(i) >= 3L && stats::sd(z0$fit[i]) > 0 && stats::sd(f_new[i]) > 0)
          round(stats::cor(z0$fit[i], f_new[i]), 5) else NA_real_,
        stringsAsFactors = FALSE)
      rm(b2, b2V)
    }
  }

  # --- collect -----------------------------------------------------------------
  FIT[[.i]] <- fit_row; TERMS[[.i]] <- TM; RESID[[.i]] <- res_row; BINS[[.i]] <- bn; STRUCT[[.i]] <- st
  if (!is.null(CV)) { CUR1[[.i]] <- CV; CSUM[[.i]] <- CS }
  if (!is.null(SF)) { SUR2[[.i]] <- SF; SSUM[[.i]] <- SS }
  FRAMES[[key]] <- list(n = nrow(d), stay_id = d$stay_id, eta = eta)   # for section G, in memory only

  if (!NO_FIGS) {
    if (!is.null(CV)) plot_curves(run, key, CV, TM)
    if (!is.null(SF)) plot_surfaces(run, key, SF)
  }

  n_dir <- if (!is.null(CS)) sum(!CS$dir_ok, na.rm = TRUE) else 0L
  out_v <- if (!is.null(CS)) max(CS$frac_var_outside, na.rm = TRUE) else NA_real_
  cat(sprintf("  %-28s conv=%s |g|=%.1e H+=%s edf=%5.1f sat=%.2f(%s) k-idx=%s ks=%.3f bins=%d dir-=%d out=%s\n",
              key, if (fit_row$converged) "y" else "N", fit_row$grad_norm,
              if (isTRUE(fit_row$hess_pd)) "y" else "N", fit_row$edf_total, fit_row$edf_ratio_max,
              gsub("^s\\(|^ti\\(|\\)$", "", fit_row$edf_worst),
              if (is.na(kc_min)) "-" else sprintf("%.3f(p=%.2f)", kc_min, kc_p),
              res_row$rqr_ks, res_row$bins_out, n_dir,
              if (is.na(out_v)) "-" else sprintf("%.2f", out_v)))
  rm(b, bV, d)
}
cat(sprintf("\n  %d spec(s) in %.1f min\n", nrow(SPECS), t_all()$elapsed_sec / 60))

FIT <- do.call(rbind, FIT); TERMS <- do.call(rbind, TERMS); RESID <- do.call(rbind, RESID)
BINS <- do.call(rbind, BINS); STRUCT <- do.call(rbind, STRUCT)
CUR1 <- do.call(rbind, Filter(Negate(is.null), CUR1)); CSUM <- do.call(rbind, Filter(Negate(is.null), CSUM))
SUR2 <- do.call(rbind, Filter(Negate(is.null), SUR2)); SSUM <- do.call(rbind, Filter(Negate(is.null), SSUM))
REFIT <- if (length(REFIT)) do.call(rbind, REFIT) else NULL
rownames(FIT) <- rownames(TERMS) <- rownames(RESID) <- NULL

# ==============================================================================
# G. GENERALISATION: what the out-of-fold L says about each spec
# ==============================================================================
cat("\n=== G. generalisation: held-out deviance from the out-of-fold L ===\n\n")
LO <- .run_tab("tables", "l_oof")
GEN <- list()
for (.i in seq_len(nrow(SPECS))) {
  key <- SPECS$key[.i]; sg <- SPECS$signal[.i]; md <- SPECS$model[.i]
  z <- LO[LO$signal == sg & LO$model == md & LO$role == "oof", , drop = FALSE]
  if (!nrow(z)) next
  pb <- DO[DO$signal == sg & DO$model == md, c("fold", "p_bar")]
  pbf <- pb$p_bar[match(z$fold, pb$fold)]
  if (anyNA(pbf)) stop("gam_qc [", key, "]: no out-of-fold p_bar for some fold", call. = FALSE)
  yz <- y_of(z$stay_id)
  p  <- inv_logit(z$l + logit(pbf))
  ll <- function(p) -(yz * log(p) + (1 - yz) * log(1 - p))
  d_mod <- 2 * sum(ll(p)); d_null <- 2 * sum(ll(pbf))
  fr <- FRAMES[[key]]
  m  <- match(as.character(z$stay_id), as.character(fr$stay_id))
  l_in <- fr$eta[m] - logit(priors_for(bundle$priors, sg, "final")$p_bar)
  r0 <- FIT[FIT$key == key, ]
  GEN[[length(GEN) + 1L]] <- data.frame(
    key = key, signal = sg, model = md, n_oof = nrow(z), n_events_oof = sum(yz),
    dev_expl_in = r0$dev_expl,
    dev_expl_oof = round(1 - d_mod / d_null, 5),
    optimism = round(r0$dev_expl - (1 - d_mod / d_null), 5),
    logloss_oof = round(mean(ll(p)), 5), logloss_null = round(mean(ll(pbf)), 5),
    auroc_oof = round(.auroc(z$l, yz), 5), auprc_oof = round(.auprc(z$l, yz), 5),
    auprc_lift_oof = round(.auprc(z$l, yz) / mean(yz), 3),
    l_oof_sd = round(stats::sd(z$l), 5), l_in_sd = round(stats::sd(l_in, na.rm = TRUE), 5),
    l_cor_in_vs_oof = round(stats::cor(l_in, z$l, use = "complete.obs"), 5),
    l_mad_in_vs_oof = round(mean(abs(l_in - z$l), na.rm = TRUE), 5),
    stringsAsFactors = FALSE)
}
GEN <- do.call(rbind, GEN)
if (nrow(GEN) < nrow(SPECS)) {
  cat(sprintf("  *** %d spec(s) have no out-of-fold L rows in %s and get no held-out number: %s\n",
              nrow(SPECS) - nrow(GEN), basename(RUN_D),
              paste(setdiff(SPECS$key, GEN$key), collapse = ", ")))
}
# fold-to-fold range of the out-of-fold diagnostics, per spec
DOs <- DO[paste(DO$signal, DO$model, sep = "/") %in% SPECS$key, , drop = FALSE]
FS <- do.call(rbind, lapply(split(DOs, paste(DOs$signal, DOs$model, sep = "/")), function(z) {
  # NA, not -Inf, for a column with no value: the two single-smooth `intv`
  # specs carry no concurvity pair, so their range is undefined.
  rg <- function(col) { v <- z[[col]][!is.na(z[[col]])]; if (length(v)) round(max(v) - min(v), 5) else NA_real_ }
  data.frame(key = paste(z$signal[1], z$model[1], sep = "/"), signal = z$signal[1], model = z$model[1],
             n_folds = nrow(z), n_unconverged = sum(!z$converged),
             dev_expl_range = rg("dev_expl"), edf_total_range = rg("edf_total"),
             k_index_range = rg("k_index_min"), concurvity_range = rg("concurvity_max"),
             l_sd_range = rg("l_sd"), stringsAsFactors = FALSE)
}))
FS <- FS[FS$key %in% SPECS$key, ]; rownames(FS) <- NULL
cat(sprintf("%-28s %8s %8s %8s %7s %7s %7s\n", "spec", "dev_in", "dev_oof", "optim", "auroc", "lift", "cor_L"))
for (i in order(-GEN$optimism)) cat(sprintf("%-28s %8.4f %8.4f %+8.4f %7.4f %7.3f %7.4f\n",
  GEN$key[i], GEN$dev_expl_in[i], GEN$dev_expl_oof[i], GEN$optimism[i], GEN$auroc_oof[i],
  GEN$auprc_lift_oof[i], GEN$l_cor_in_vs_oof[i]))
cat(sprintf("\n  optimism (in-sample minus held-out deviance explained): median %+.4f, max %+.4f (%s)\n",
            stats::median(GEN$optimism), max(GEN$optimism), GEN$key[which.max(GEN$optimism)]))
cat(sprintf("  correlation of the final fit's in-sample L with the out-of-fold L: median %.4f, min %.4f (%s)\n",
            stats::median(GEN$l_cor_in_vs_oof), min(GEN$l_cor_in_vs_oof), GEN$key[which.min(GEN$l_cor_in_vs_oof)]))

# arm-level reference points, joined rather than recomputed
ARMS <- tryCatch(.run_tab("tables", "train_ref_arms"), error = function(e) NULL)
if (is.na(EXT_D)) EXT_D <- latest_run("external") %||% NA_character_
TRANS <- NULL
if (!is.na(EXT_D) && file.exists(file.path(EXT_D, "tables", "transport.rds"))) {
  TRANS <- readRDS(file.path(EXT_D, "tables", "transport.rds"))
  cat(sprintf("\n  arm-level transport table joined from %s\n", basename(EXT_D)))
} else {
  cat(sprintf("\n  no arm-level transport table joined (%s)\n",
              if (is.na(EXT_D)) "no completed external run found" else
                paste0(EXT_D, " has no tables/transport.rds")))
}

# ==============================================================================
# H. DEPENDENCE: join the recorded concurvity against its calibrated null
# ==============================================================================
cat("\n=== H. dependence: concurvity against the calibrated null ===\n\n")
CN <- NULL
if (dir.exists(NULL_D) && file.exists(file.path(NULL_D, "diagnostics", "concurvity_null.rds"))) {
  CN <- readRDS(file.path(NULL_D, "diagnostics", "concurvity_null.rds"))
}
DEP <- DF[paste(DF$signal, DF$model, sep = "/") %in% SPECS$key,
          c("signal", "model", "n_smooth", "concurvity_max", "concurvity_obs", "concurvity_est", "concurvity_term")]
DEP$key <- paste(DEP$signal, DEP$model, sep = "/")
DEP$over_threshold <- !is.na(DEP$concurvity_max) & DEP$concurvity_max > DIAG_THR$concurvity_max
if (!is.null(CN)) {
  for (m in c("worst", "observed")) {
    z <- CN[CN$measure == m, ]
    k <- match(DEP$key, paste(z$signal, z$model, sep = "/"))
    DEP[[paste0(m, "_null_used")]]  <- z$null_used[k]
    DEP[[paste0(m, "_excess")]]     <- z$excess_used[k]
    DEP[[paste0(m, "_null_source")]] <- z$null_source[k]
    DEP[[paste0(m, "_pair")]]       <- z$worst_pair[k]
  }
  cat(sprintf("  null run: %s\n", basename(NULL_D)))
  cat(sprintf("  specs with a null row     : %d of %d (a single-smooth model has no pair)\n",
              sum(!is.na(DEP$worst_excess)), nrow(DEP)))
  cat(sprintf("  `worst` over the flag %.2f : %d of %d;  excess over its EXACT permutation null > 0: %d\n",
              DIAG_THR$concurvity_max, sum(DEP$over_threshold), nrow(DEP), sum(DEP$worst_excess > 0, na.rm = TRUE)))
  cat(sprintf("  `observed` excess over the synthetic null > 0: %d (approximate null, labelled as such)\n",
              sum(DEP$observed_excess > 0, na.rm = TRUE)))
  cat(sprintf("  worst-pair excess, `worst` measure: median %+.4f, max %+.4f (%s)\n",
              stats::median(DEP$worst_excess, na.rm = TRUE), max(DEP$worst_excess, na.rm = TRUE),
              DEP$key[which.max(DEP$worst_excess)]))
} else {
  cat("  NO CALIBRATED NULL FOUND at ", NULL_D, ". Raw concurvity is reported and is not\n")
  cat("  interpretable on its own at this design (see tests/concurvity_null.R).\n")
}
rownames(DEP) <- NULL

# ==============================================================================
# I. THE SPEC SUMMARY, AND THE FLAGS
# ==============================================================================
term_agg <- function(df, key, col, fun) if (is.null(df)) NA else {
  z <- df[[col]][df$key == key]; z <- z[!is.na(z)]; if (length(z)) fun(z) else NA }
SUM <- FIT[, c("key", "signal", "model", "n_rows", "n_events", "converged", "hess_pd", "grad_norm",
               "n_sp_huge", "n_sp_tiny", "edf_total", "n_terms_shrunk_out", "n_terms_saturated",
               "edf_ratio_max", "edf_worst", "k_index_min", "k_index_p", "k_worst", "k_flag",
               "k_index_recorded", "dev_expl")]
SUM <- merge(SUM, RESID[, c("key", "rqr_ks", "rqr_skew", "rqr_kurt", "cal_in_large", "bins_out",
                            "hosmer_lemeshow", "atom_block_frac", "atom_block_z", "cal_slope_in",
                            "struct_max_abs_mean", "struct_max_abs_z", "struct_worst", "frac_separated")],
             by = "key", all.x = TRUE)
SUM <- merge(SUM, GEN[, c("key", "dev_expl_oof", "optimism", "auroc_oof", "auprc_lift_oof", "l_cor_in_vs_oof")],
             by = "key", all.x = TRUE)
SUM <- merge(SUM, DEP[, intersect(c("key", "concurvity_max", "concurvity_obs", "worst_excess", "observed_excess"),
                                  names(DEP))], by = "key", all.x = TRUE)
SUM$n_terms_1d <- vapply(SUM$key, function(k) sum(CSUM$key == k), integer(1))
SUM$n_terms_2d <- vapply(SUM$key, function(k) if (is.null(SSUM)) 0L else sum(SSUM$key == k), integer(1))
SUM$out_var_max     <- vapply(SUM$key, function(k) term_agg(CSUM, k, "frac_var_outside", max), numeric(1))
SUM$out_curv_max    <- vapply(SUM$key, function(k) term_agg(CSUM, k, "frac_curv_outside", max), numeric(1))
SUM$resolved_min    <- vapply(SUM$key, function(k) term_agg(CSUM, k, "resolved_frac_in", min), numeric(1))
SUM$n_dir_contra    <- vapply(SUM$key, function(k) sum(!CSUM$dir_ok[CSUM$key == k], na.rm = TRUE), integer(1))
SUM$n_dir_declared  <- vapply(SUM$key, function(k) sum(!is.na(CSUM$dir_ok[CSUM$key == k])), integer(1))
SUM$surf_cells_supported_min <- vapply(SUM$key, function(k) term_agg(SSUM, k, "frac_cells_supported", min), numeric(1))
# Always present, NA when no refit ran, so spec_summary has one column set
# whatever the flags were.
SUM$refit_k_shape_max <- vapply(SUM$key, function(k) term_agg(REFIT, k, "shape_sd_in", max), numeric(1))
SUM$refit_k_dev_gain_max <- vapply(SUM$key, function(k) {
  if (is.null(REFIT)) return(NA_real_)
  z <- REFIT[REFIT$key == k & REFIT$status == "ok", ]; if (!nrow(z)) NA_real_ else max(z$dev_expl_new - z$dev_expl_old) }, numeric(1))
SUM <- SUM[match(SPECS$key, SUM$key), ]; rownames(SUM) <- NULL

# --- write -------------------------------------------------------------------
save_table(run, data.frame(
  item = c("bundle", "bundle_run", "bundle_design_hash", "n_models", "n_specs_qc",
           "store_coef_identical", "store_frames_checked", "mgcv_version", "concurvity_null", "external_run"),
  value = c(BUNDLE_P, basename(RUN_D), .hash(bundle$cfg), length(keys), nrow(SPECS),
            if (is.na(store_coef)) "store absent" else sprintf("%d/%d", store_coef, length(keys)),
            if (all(is.na(FIT$frame_identical_to_store))) "not checked (store absent or stale)"
            else sprintf("%d/%d identical", sum(FIT$frame_identical_to_store, na.rm = TRUE), nrow(FIT)),
            as.character(utils::packageVersion("mgcv")), basename(NULL_D),
            if (is.na(EXT_D)) "" else basename(EXT_D)),
  stringsAsFactors = FALSE), "provenance", subdir = "diagnostics")
save_table(run, FIT,   "fit_qc",        subdir = "diagnostics")
save_table(run, TERMS, "fit_terms",     subdir = "diagnostics")
save_table(run, RESID, "residual_qc",   subdir = "diagnostics")
save_table(run, BINS,  "binned_residuals", subdir = "diagnostics")
save_table(run, STRUCT, "residual_structure", subdir = "diagnostics")
if (!is.null(REFIT)) save_table(run, REFIT, "refit_k", subdir = "diagnostics")
save_table(run, DEP,   "concurvity_join", subdir = "diagnostics")
save_table(run, GEN,   "generalisation")
save_table(run, FS,    "fold_stability")
if (!is.null(ARMS))  save_table(run, ARMS,  "arms_oof")
if (!is.null(TRANS)) save_table(run, TRANS, "arms_transport")
save_table(run, CUR1, "curves_1d", csv = FALSE)          # long; .rds only
save_table(run, CUR1[, c("key", "term", "variable", "kind", "grid_x", "n", "events", "supported")],
           "support_1d", csv = FALSE)
save_table(run, CSUM, "curve_summary")
if (!is.null(SUR2)) {
  save_table(run, SUR2, "surfaces_2d", csv = FALSE)
  save_table(run, SSUM, "surface_summary")
}
save_table(run, SUM, "spec_summary")

# ==============================================================================
# THE READING
# ==============================================================================
cat("\n=== summary over ", nrow(SUM), " spec(s) ===\n\n", sep = "")
cat("FIT\n")
cat(sprintf("  converged: %d;  Hessian PD (to 1e-6 of its largest eigenvalue): %d;  flat directions: %d over %d spec(s);  max |gradient|: %.2e\n",
            sum(SUM$converged), sum(SUM$hess_pd, na.rm = TRUE), sum(FIT$n_hess_flat, na.rm = TRUE),
            sum(FIT$n_hess_flat > 0, na.rm = TRUE), max(SUM$grad_norm, na.rm = TRUE)))
if (any(!is.na(FIT$n_iter))) cat(sprintf("  iterations: %s\n", paste(range(FIT$n_iter, na.rm = TRUE), collapse = "-")))
cat(sprintf("  in-sample calibration slope (1 = the fitted log-odds are taken at face value): median %.3f, range %.3f-%.3f\n",
            stats::median(SUM$cal_slope_in), min(SUM$cal_slope_in), max(SUM$cal_slope_in)))
cat(sprintf("  smoothing parameters at a bound: %d huge (> %.0e, term shrunk to nothing), %d tiny (< %.0e)\n",
            sum(SUM$n_sp_huge), as.numeric(SPB$huge), sum(SUM$n_sp_tiny), as.numeric(SPB$tiny)))
cat(sprintf("  terms shrunk out (edf < %.1f): %d of %d;  saturated (edf/k' > %.2f): %d of %d over %d spec(s)\n",
            as.numeric(SPB$shrunk_edf), sum(TERMS$shrunk_out), nrow(TERMS), DIAG_THR$edf_ratio_max,
            sum(TERMS$saturated), nrow(TERMS), sum(SUM$n_terms_saturated > 0)))
if (!NO_KCHECK) {
  cat(sprintf("  k-index on ALL rows: min %.3f;  flagged (index < %.2f and p < %.2f): %d spec(s)\n",
              min(SUM$k_index_min, na.rm = TRUE), DIAG_THR$k_index_min, DIAG_THR$k_index_p_max, sum(SUM$k_flag)))
  dk <- SUM$k_index_min - SUM$k_index_recorded
  cat(sprintf("  against diag_final's recorded k-index (one unseeded 5,000-row subsample): median |diff| %.4f, max %.4f;\n",
              stats::median(abs(dk), na.rm = TRUE), max(abs(dk), na.rm = TRUE)))
  cat(sprintf("    worst term agrees on %d of %d spec(s)\n", sum(SUM$k_worst == FIT$k_worst_recorded, na.rm = TRUE), nrow(SUM)))
  cat(sprintf("  k' from k.check against the coefficient count: disagree on %d of %d terms (%d ti)\n",
              sum(!TERMS$k_prime_agrees), nrow(TERMS), sum(!TERMS$k_prime_agrees & TERMS$dim == 2L)))
}
if (!is.null(REFIT)) {
  ok <- REFIT[REFIT$status == "ok", ]
  cat(sprintf("  k-doubling refits: %d done, %d skipped (k bounded by distinct values), %d failed\n",
              nrow(ok), sum(grepl("^skipped", REFIT$status)), sum(grepl("^refit failed", REFIT$status))))
  if (nrow(ok)) {
    cat(sprintf("    deviance-explained gain: median %+.5f, max %+.5f;  curve shape change inside support (sd, nats): median %.4f, max %.4f (%s/%s)\n",
                stats::median(ok$dev_expl_new - ok$dev_expl_old), max(ok$dev_expl_new - ok$dev_expl_old),
                stats::median(ok$shape_sd_in, na.rm = TRUE), max(ok$shape_sd_in, na.rm = TRUE),
                ok$key[which.max(ok$shape_sd_in)], ok$variable[which.max(ok$shape_sd_in)]))
    cat(sprintf("    edf/k' after doubling: median %.3f, max %.3f -- a ratio that halves says the basis was binding, one that does not says it was not\n",
                stats::median(ok$edf_ratio_new), max(ok$edf_ratio_new)))
  }
}
cat("RESIDUALS\n")
cat(sprintf("  RQR KS vs N(0,1): median %.4f, max %.4f (%s);  |skew| max %.3f;  excess kurtosis max %.3f\n",
            stats::median(SUM$rqr_ks), max(SUM$rqr_ks), SUM$key[which.max(SUM$rqr_ks)],
            max(abs(SUM$rqr_skew)), max(SUM$rqr_kurt)))
cat(sprintf("  binned residuals: spec(s) with any bin beyond %.0f SE: %d;  max Hosmer-Lemeshow over %d bins: %.1f (%s)\n",
            as.numeric(RS$z_out), sum(SUM$bins_out > 0), as.integer(RS$n_bins),
            max(SUM$hosmer_lemeshow), SUM$key[which.max(SUM$hosmer_lemeshow)]))
if (any(!is.na(SUM$atom_block_z))) {
  cat(sprintf("  dominant tie block (the point mass): share of rows median %.2f, max %.2f;  its |obs - exp| z: median %.2f, max %.2f (%s)\n",
              stats::median(SUM$atom_block_frac), max(SUM$atom_block_frac),
              stats::median(abs(SUM$atom_block_z), na.rm = TRUE), max(abs(SUM$atom_block_z), na.rm = TRUE),
              SUM$key[which.max(abs(SUM$atom_block_z))]))
}
cat(sprintf("  leftover structure vs a smooth covariate, max |binned mean RQR| in SD units: median %.3f, max %.3f (%s on %s);  as z: max %.1f\n",
            stats::median(SUM$struct_max_abs_mean), max(SUM$struct_max_abs_mean),
            SUM$key[which.max(SUM$struct_max_abs_mean)], SUM$struct_worst[which.max(SUM$struct_max_abs_mean)],
            max(SUM$struct_max_abs_z)))
cat(sprintf("  calibration in the large |mean p - mean y|: max %.2e;  separated fitted values: max fraction %.5f\n",
            max(abs(SUM$cal_in_large)), max(SUM$frac_separated)))
cat("GENERALISATION\n")
cat(sprintf("  held-out deviance explained: median %.4f, min %.4f (%s);  optimism median %+.4f, max %+.4f\n",
            stats::median(SUM$dev_expl_oof, na.rm = TRUE), min(SUM$dev_expl_oof, na.rm = TRUE),
            SUM$key[which.min(SUM$dev_expl_oof)], stats::median(SUM$optimism, na.rm = TRUE),
            max(SUM$optimism, na.rm = TRUE)))
cat(sprintf("  spec(s) whose held-out deviance explained is below the bundle's dev_expl_min %.3f: %d\n",
            DIAG_THR$dev_expl_min, sum(SUM$dev_expl_oof < DIAG_THR$dev_expl_min, na.rm = TRUE)))
cat("SUPPORT AND INTERPRETABILITY\n")
cat(sprintf("  1-D smooths: %d;  supported share of the grid: median %.2f, min %.2f;  rows inside supported bins: median %.3f\n",
            nrow(CSUM), stats::median(CSUM$n_supported / CSUM$n_grid), min(CSUM$n_supported / CSUM$n_grid),
            stats::median(CSUM$frac_rows_supported)))
cat(sprintf("  share of curve variation OUTSIDE support: median %.3f, max %.3f (%s %s);  smooths above %.2f: %d\n",
            stats::median(CSUM$frac_var_outside, na.rm = TRUE), max(CSUM$frac_var_outside, na.rm = TRUE),
            CSUM$key[which.max(CSUM$frac_var_outside)], CSUM$term[which.max(CSUM$frac_var_outside)],
            as.numeric(REPORT$frac_var_outside_high),
            sum(CSUM$frac_var_outside > as.numeric(REPORT$frac_var_outside_high), na.rm = TRUE)))
cat(sprintf("  resolved from zero inside support (|f| > %.2f SE): median %.2f of the supported grid; smooths under %.2f: %d\n",
            Z, stats::median(CSUM$resolved_frac_in, na.rm = TRUE), as.numeric(REPORT$resolved_low),
            sum(CSUM$resolved_frac_in < as.numeric(REPORT$resolved_low), na.rm = TRUE)))
dt <- table(factor(CSUM$direction, levels = c("with", "against", "non_monotone", "undeclared", "flat", "shrunk_out")))
cat(sprintf("  direction inside support: with the declaration %d, against it %d, non-monotone (|rho| < %.2f) %d, undeclared %d, flat %d, shrunk out %d\n",
            dt[["with"]], dt[["against"]], RHO_MONO, dt[["non_monotone"]], dt[["undeclared"]], dt[["flat"]], dt[["shrunk_out"]]))
if (dt[["against"]] > 0) {
  z <- CSUM[CSUM$direction == "against", ]
  cat("    against (reported, never corrected):\n")
  cat(sprintf("    %s\n", paste(sprintf("%s %s (rho %+.2f, range %.2f nats)", z$key, z$variable, z$slope_rho_in, z$range_in),
                                collapse = "\n    ")))
}
if (dt[["non_monotone"]] > 0) {
  z <- CSUM[CSUM$direction == "non_monotone", ]
  cat(sprintf("    non-monotone: %s\n", paste(sprintf("%s %s", z$key, z$variable), collapse = "; ")))
}
if (!is.null(SSUM)) {
  cat(sprintf("  ti surfaces: %d;  supported share of cells: median %.2f, min %.2f;  rows inside supported cells: median %.3f\n",
              nrow(SSUM), stats::median(SSUM$frac_cells_supported), min(SSUM$frac_cells_supported),
              stats::median(SSUM$frac_rows_supported)))
  cat(sprintf("  surface range inside support: median %.3f nats, max %.3f (%s %s);  resolved share median %.2f\n",
              stats::median(SSUM$range_in, na.rm = TRUE), max(SSUM$range_in, na.rm = TRUE),
              SSUM$key[which.max(SSUM$range_in)], SSUM$term[which.max(SSUM$range_in)],
              stats::median(SSUM$resolved_frac_in, na.rm = TRUE)))
  cat(sprintf("  surface max |f| outside support / inside: median ratio %.2f -- the size of the extrapolated part relative to the supported part\n",
              stats::median(SSUM$max_abs_out / SSUM$max_abs_in, na.rm = TRUE)))
}

finalize_run(run, extra = list(n_specs = nrow(SPECS), n_terms_1d = nrow(CSUM),
                               n_terms_2d = if (is.null(SSUM)) 0L else nrow(SSUM),
                               kcheck = !NO_KCHECK, refit_k = DO_REFIT,
                               concurvity_null = basename(NULL_D)))
cat(sprintf("\nwritten: %s\n", run$path))
cat("  NEXT: Rscript tests/gam_qc_bootstrap.R --qc ", run$path, "\n\n", sep = "")
