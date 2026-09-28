# run/internal.R -------------------------------------------------------------
# THE ONLY RUNNER THAT FITS ANYTHING (hard rule 8).
#
# Two steps and no analysis. `targets::tar_make()` computes and caches the whole
# training graph; `export_run()` snapshots the named values into a dated,
# immutable run directory. The split is hard rule 7: the moment a run directory
# becomes a target, or a target's value carries a timestamp, the cache is gone
# and 258 GAMs recompute to redraw a figure.
#
# This runner does NOT touch MIMIC test. Test is spent once, by
# `run/test_look.R`, which loads the bundle this run writes and fits nothing.
# Keeping the two apart means nobody spends the test look by typing the command
# that rebuilds a plot.
#
# NOTHING ROW-LEVEL IS PRINTED (hard rule 1). The console gets counts,
# eigenvalues and check results. `l_oof` and `oof_scores` are exported as
# objects, so they land in the run directory as .rds and are never rendered.
#
#   Rscript run/internal.R
#   Rscript run/internal.R --dry          # what would rebuild, and nothing else
#   Rscript run/internal.R --note "..."   # free text into the manifest
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(targets); library(mgcv); library(arrow); library(yaml)
  library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
dry  <- "--dry" %in% args
note <- {
  i <- match("--note", args)
  if (!is.na(i) && length(args) > i) args[i + 1L] else NULL
}

rc  <- yaml::read_yaml("config/internal.yml")
cfg <- load_config(rc$config)

# THE RUNNER'S CONFIG MUST BE THE GRAPH'S CONFIG (plumbing review F1,
# 2026-09-08). `_targets.R` tracks its own config, pairing and domain files as
# `format = "file"` targets and never reads `internal.yml`, so until now
# `internal.yml$config` could name any valid design and this runner would print
# that design's budget, export that design into the manifest and describe a
# bundle fitted under a different one. Nothing compared the two. Two checks,
# neither of which restates the graph's path: before `tar_make()`, the file the
# runner loaded must be the file the graph tracks (read from the store when one
# exists); after it, the object the runner holds must be identical to the
# `cfg` target, and the graph's object is what reaches the export.
r1 <- function(nm) targets::tar_read_raw(nm, store = rc$store)
# A `format = "file"` target has no object file, only a meta record, so the
# existence test is the read itself; a fresh store has neither and the object
# check after `tar_make()` covers that case.
graph_cfg_file <- if (dir.exists(rc$store)) {
  tryCatch(r1("config_file"), error = function(e) NULL)
} else NULL
if (!is.null(graph_cfg_file)) {
  if (normalizePath(rc$config, mustWork = TRUE) !=
      normalizePath(graph_cfg_file, mustWork = TRUE)) {
    stop("run/internal.R: config/internal.yml names `", rc$config, "` but the ",
         "graph fits from `", graph_cfg_file, "`. The runner would describe one ",
         "design and export another. Point internal.yml at the graph's file.",
         call. = FALSE)
  }
}

# --- what is about to happen ------------------------------------------------
cat("\n=== internal run: fit MIMIC-IV train, build the bundle ===\n\n")
print(layer1_budget(cfg))
cat("\n")

out <- targets::tar_outdated(store = rc$store, reporter = "silent")
if (!length(out)) {
  cat("  every target is up to date; tar_make() will recompute nothing.\n")
} else {
  cat(sprintf("  %d target(s) will rebuild: %s\n", length(out),
              paste(utils::head(out, 15), collapse = ", ")))
  if (length(out) > 15) cat(sprintf("  ... and %d more\n", length(out) - 15))
}
if (dry) {
  cat("\n  --dry: stopping before tar_make().\n\n")
  quit(save = "no")
}

# --- compute ----------------------------------------------------------------
# Timing is measured HERE, at the call site, and printed. It never enters a
# target's value (hard rule 7).
t0 <- Sys.time()
targets::tar_make(store = rc$store)
cat(sprintf("\n  tar_make(): %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# --- the config that was actually fitted ------------------------------------
# Second half of the F1 guard: the object, not the path. `identical()` rather
# than a design-key comparison, because the manifest snapshots the whole thing.
cfg_graph <- r1("cfg")
if (!identical(cfg, cfg_graph)) {
  stop("run/internal.R: the config this runner loaded is not the `cfg` target ",
       "the graph fitted under. Same file path, different object -- a pairing ",
       "or domain file the graph tracks and the runner did not, or an edit ",
       "between the two loads. Nothing is exported.", call. = FALSE)
}
cfg <- cfg_graph

# --- read back the things worth reading -------------------------------------

cat("\n=== the branch point (L_full, zero-filled, 19 signals) ===\n\n")
print(r1("eigen_full_zero")[1:6, ], row.names = FALSE)
bp <- r1("branch_pt")
cat("\n")
report_branch_point(bp)

cat("\n=== layer-1 diagnostics: what tripped a threshold ===\n\n")
tri <- r1("triage_oof")
if (!nrow(tri)) cat("  nothing tripped a configured threshold.\n") else
  print(utils::head(tri, 20), row.names = FALSE)

cat("\n=== out-of-fold discrimination, all arms (MIMIC train) ===\n\n")
tr <- r1("train_ref")
print(tr$arms[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
                  "auprc", "auprc_lift", "spearman", "cal_slope")],
      row.names = FALSE)

cat("\n=== per-signal discrimination, L_full (top 8 by AUROC) ===\n")
cat("  `auprc_lift` is AUPRC over the event rate; 1.0 is a column that has\n")
cat("  learned nothing. Read it beside AUROC: at a 12% event rate a column can\n")
cat("  order the cohort acceptably and add nothing where the deaths are.\n\n")
print(utils::head(r1("signal_auroc_full"), 8), row.names = FALSE)

cat("\n=== bundle ===\n\n")
b <- r1("bundle")
print(b)
cat("\n")
print(r1("bundle_checks"), row.names = FALSE)

# --- the reporting pass -----------------------------------------------------
# THE RISK-ORDERING FIGURES AND THE PER-ARM MONOTONICITY TABLES.
#
# WHY THIS IS A CALLBACK AND NOT A TARGET. `plot_risk_curve()` writes a PNG, and
# a PNG needs a directory. Hard rule 9 forbids anything in `R/` from building a
# path and hard rule 7 forbids a target's value from carrying a timestamp, so
# the graph can compute `monotonicity()` -- and does, inside `train_ref` -- but
# it has nowhere to put a figure. The consequence until 2026-09-05 was that
# run/test_look.R and run/external.R each produced a full set of risk-ordering
# curves and per-arm monotonicity tables while the internal run produced none,
# and nothing errored to say so: the graph reported `spearman` and `rate_ratio`
# as two columns of `train_ref` and the curves behind them existed nowhere.
#
# `export_run()` owns the run directory and hands it to this function, so the
# library layer still constructs no path.
#
# SELF-BINNED ONLY, AND DELIBERATELY. `cutpoints` is frozen FROM these very
# out-of-fold scores, so binning them on it would reproduce the self-binned
# curve to within tied values. The frozen-versus-self contrast is a transport
# question and it belongs at the apply sites, where the two genuinely differ.
report_oof <- function(run) {
  sc <- r1("oof_scores")
  yt <- r1("y_train")
  pb <- r1("p_bar_cohort_train")
  nbins <- cfg$metrics$n_bins %||% 20L
  nboot <- cfg$metrics$n_boot %||% 200L

  out <- score_arms(run, sc, yt, pb, breaks = NULL, suffix = "_oof",
                    n_bins = nbins, n_boot = nboot, seed = cfg$seed,
                    group = r1("group_of_train"))
  save_table(run, out$summary, "score_summary_oof")

  # The graph computed these numbers once already, inside `train_ref`, from the
  # same scores with the same seed and the same bootstrap count, so they must
  # agree exactly. Checking rather than assuming is nearly free, and it is the
  # only thing standing between a reporting pass that reads `oof_scores` and one
  # that reads something subtly else -- a stale target, or an arm added in one
  # place and not the other.
  ref <- tr$arms
  cmp <- merge(data.frame(label = sub("_oof$", "", out$summary$label),
                          auroc_here = out$summary$auroc,
                          auprc_here = out$summary$auprc,
                          stringsAsFactors = FALSE),
               data.frame(label = ref$label, auroc_ref = ref$auroc,
                          auprc_ref = ref$auprc, stringsAsFactors = FALSE),
               by = "label", all = TRUE)
  bad <- cmp[is.na(cmp$auroc_here) | is.na(cmp$auroc_ref) |
               abs(cmp$auroc_here - cmp$auroc_ref) > 1e-9 |
               abs(cmp$auprc_here - cmp$auprc_ref) > 1e-9, , drop = FALSE]
  if (nrow(bad)) {
    print(bad, row.names = FALSE)
    stop("report_oof: the reporting pass disagrees with `train_ref` on ",
         nrow(bad), " arm(s). Same scores, same seed, same n_boot -- they ",
         "cannot differ, so one of the two is not reading what it claims to.",
         call. = FALSE)
  }
  log_msg(run, sprintf("reporting pass: %d arm(s); figures written; agrees with train_ref",
                       nrow(cmp)))

  mono <- do.call(rbind, lapply(names(out$reports), function(nm)
    cbind(arm = nm, out$reports[[nm]]$mono, stringsAsFactors = FALSE)))
  save_table(run, mono, "monotonicity_oof")

  cat("\n=== monotonicity of the out-of-fold risk ordering (MIMIC train) ===\n")
  cat("  `n_sig_inversions` is the one to read: adjacent bins that go the wrong\n")
  cat("  way AND whose Wilson intervals do not overlap. A handful of\n")
  cat("  overlapping inversions in 20 bins is sampling noise.\n\n")
  print(mono[, c("arm", "n_bins", "spearman", "n_inversions", "n_sig_inversions",
                 "worst_inversion", "rate_bottom", "rate_top", "rate_ratio")],
        row.names = FALSE)
  invisible(NULL)
}

# --- snapshot ---------------------------------------------------------------
# export_run() is a PLAIN FUNCTION and never a target. It is what reads the
# clock; the graph never does. `after = report_oof` runs inside the run's
# lifetime, so a failure there leaves the manifest at status "running", which is
# the correct statement about a run whose reporting step did not finish.
ex <- rc$export
path <- export_run(
  "internal", cfg,
  tables      = as.character(unlist(ex$tables      %||% character(0))),
  diagnostics = as.character(unlist(ex$diagnostics %||% character(0))),
  objects     = as.character(unlist(ex$objects     %||% character(0))),
  bundle      = ex$bundle,
  store       = rc$store,
  after       = report_oof,
  note        = note %||% sprintf(
    "MIMIC-IV train: %d stays, %d fits, PC1 = %s. Bundle written for test and eICU.",
    tr$n, layer1_budget(cfg)$fits_total,
    round(bp$pc1[bp$model == "full" & bp$fill == "zero" &
                   bp$gcs == "all three"][1], 4)))

cat(sprintf("\n  run directory: %s\n", path))
cat(sprintf("  bundle:        %s\n", file.path(path, "bundle.qs2")))
cat("\n  NEXT: put that bundle path into config/external.yml (and into\n")
cat("  config/internal.yml under test_look.bundle). Both are explicit paths on\n")
cat("  purpose -- a published result must not change because a newer run\n")
cat("  appeared.\n\n")
