# run/attribution.R -----------------------------------------------------------
# THE ATTRIBUTABILITY ARM, MIMIC-ONLY. Steps 3 and 4 of
# `docs/v2_attributability_plan_20260902.md` section 12.
#
# WHAT THIS PRODUCES, AND WHY IT IS A RESULT ON ITS OWN. For each of the five
# arms, how much does its per-signal attribution move when the estimator is
# refitted on perturbed training data? That is the NOISE FLOOR, and without it
# an agreement score is a number with no referent: "the two sites agree at
# 0.71" means nothing until you know that one site's own resampling produces
# 0.68. The floors are MIMIC-only, so none of this waits on eICU, and the
# ordering they induce over methods -- the plan's section 10 predicts
# `llr_sum`, then SHAP on `xgb_l`, then `xgb_feat`, then `xgb_raw` -- is
# publishable by itself.
#
# HARD RULE 8, AND THE ONE PLACE THIS RUNNER FITS. `refit_replicate()` refits
# alpha, delta, lambda, the 43 GAMs and the three boosters on a resampled
# subset of MIMIC TRAIN. That is a fit, and it is legitimate here for a reason
# that must travel with the result: a replicate is never scored, never reported
# as a performance number, and never enters a bundle. It exists solely as the
# second argument to an agreement metric. `refit_replicate()` returns a bundle
# built by `build_bundle()` from the SAME cfg the real bundle froze, so
# `attribution_set()`'s design-stamp guard passes for the right reason rather
# than by being bypassed, and `layer2`, `sigma`, `cutpoints` and `train_ref`
# are deliberately left empty so a replicate cannot be scored even by accident.
#
# THE EVALUATION SET NEVER MOVES. It is drawn once, from the split named in
# config/attribution.yml, and every replicate of every perturbation attributes
# on those same stays. Resampling it too would confound sampling of the
# evaluation cohort with sampling of the estimator, and the estimator is the
# only thing under study.
#
# AGGREGATES ONLY (hard rule 1). The attribution matrices are n-by-19 and
# row-level. They stay in memory and go nowhere: what leaves this runner is
# correlations, quantiles, shares and counts.
#
#   Rscript run/attribution.R
#   Rscript run/attribution.R --dry            # what would run, and its cost
#   Rscript run/attribution.R --floors bootstrap
#   Rscript run/attribution.R --note "..."
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
dry  <- "--dry" %in% args
note <- { i <- match("--note", args)
          if (!is.na(i) && length(args) > i) args[i + 1L] else NULL }
only <- { i <- match("--floors", args)
          if (!is.na(i) && length(args) > i)
            trimws(strsplit(args[i + 1L], ",")[[1]]) else NULL }

ac <- yaml::read_yaml("config/attribution.yml")

# EVERY SETTING READ WITH NO FALLBACK. A threshold written down twice can be
# edited apart, and the divergence stays invisible until the config key is
# deleted or misspelled (audit findings F9 to F11). An absent key here is a
# loud stop naming the path, not a quiet substitution.
eval_split  <- cfg_req(ac, "eval", "split")
eval_n      <- cfg_req(ac, "eval", "n_stays")  # an integer, or the string "all"
eval_seed   <- cfg_req(ac, "eval", "seed")
methods     <- as.character(unlist(cfg_req(ac, "methods")))
iv_handling <- cfg_req(ac, "intervention_handling")
mag_floor   <- cfg_req(ac, "mag_floor")
min_signals <- as.integer(cfg_req(ac, "min_signals"))
domain_bins <- as.integer(cfg_req(ac, "domain_bins"))
floor_seed  <- cfg_req(ac, "floor_seed")
l_source    <- cfg_req(ac, "replicate_l_source")
bundle_path <- cfg_req(ac, "bundle")

floors <- Filter(function(f) isTRUE(f$run) &&
                   (is.null(only) || f$perturbation %in% only),
                 cfg_req(ac, "floors"))

# --- load --------------------------------------------------------------------
cfg_local <- load_config("config/config.yml")
bundle    <- load_bundle(bundle_path, cfg = cfg_local, strict = TRUE)
cfg       <- bundle_cfg(bundle, cfg_local$paths$mimiciv)

run <- new_run("attribution", cfg_local, note = note %||% sprintf(
  "Attributability floors on MIMIC %s. bundle %s. Replicates fit; nothing is scored.",
  eval_split, basename(dirname(bundle_path))))
log_msg(run, "bundle: ", bundle_path)
save_table(run, verify_bundle(bundle, cfg = cfg_local, strict = FALSE),
           "bundle_checks", subdir = "diagnostics")
save_table(run, bundle_summary(bundle), "bundle_contents", subdir = "diagnostics")

tabs  <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)

# --- the evaluation set, drawn once and then frozen --------------------------
pool <- folds$stay_id[folds$split == eval_split]
if (!length(pool)) {
  abort_values("run/attribution.R: no stays in the requested split", eval_split)
}
ids <- if (identical(as.character(eval_n), "all")) sort(pool) else {
  k <- suppressWarnings(as.integer(eval_n))
  if (is.na(k) || k < 1L) {
    abort_values("config/attribution.yml `eval.n_stays` must be a positive integer or \"all\"",
                 as.character(eval_n))
  }
  sort(with_seed(eval_seed, sample(pool, min(k, length(pool)))))
}

log_msg(run, sprintf("evaluation set: %d of %d %s stays (seed %s), FROZEN",
                     length(ids), length(pool), eval_split, eval_seed))
if (identical(eval_split, "train")) {
  log_msg(run, "NOTE: a train evaluation set is partly IN-SAMPLE for every ",
          "replicate, which inflates agreement. See config/attribution.yml.")
}

cat("\n=== attributability arm: floors on MIMIC ", eval_split, " ===\n\n", sep = "")
cat(sprintf("  evaluation set        %d stays, frozen\n", length(ids)))
cat(sprintf("  methods               %s\n", paste(methods, collapse = ", ")))
cat(sprintf("  intervention handling %s\n", iv_handling))
cat(sprintf("  magnitude floor       %.3f log-odds\n", mag_floor))
cat(sprintf("  replicate L source    %s%s\n", l_source,
            if (identical(l_source, "in_sample"))
              "   *** floor sits too high; declare it ***" else ""))
for (f in floors) {
  cat(sprintf("  floor                 %-14s B = %d  (%d replicate pairs, %d refits)\n",
              f$perturbation, f$B, f$B, 2L * f$B))
}
if (!length(floors)) cat("  floor                 none enabled\n")
if (dry) {
  cat("\n  --dry: stopping before any fit.\n\n")
  quit(save = "no")
}

# --- the reference attribution set, on the bundle itself ---------------------
# ONE pass of layer 1 serves every method (`attribution_set()` evaluates the 43
# smooths once and derives all five matrices from that), which is what makes a
# floor replicate affordable at all.
cat("\n=== reference attributions, from the bundle ===\n\n")
tm <- start_timer()
ref <- attribution_set(bundle, tabs, cfg, ids, methods = methods,
                       intervention_handling = iv_handling, verbose = TRUE)
el <- tm()
log_msg(run, sprintf("attribution_set: %.1f min elapsed, %.1f min CPU",
                     el$elapsed_sec / 60, el$cpu_sec / 60))

cells <- attribution_cells(ref$measured)
cov <- data.frame(
  method = names(ref$A),
  n_stays = vapply(ref$A, nrow, integer(1)),
  n_signals = vapply(ref$A, ncol, integer(1)),
  dropped_frac = round(unlist(ref$dropped_frac[names(ref$A)]), 5),
  frac_cells_used = round(cells$frac_used, 5),
  stringsAsFactors = FALSE, row.names = NULL)
print(cov, row.names = FALSE)
save_table(run, cov, "attribution_coverage", subdir = "diagnostics")
cat("\n  `frac_cells_used` excludes assigned zeros: an unmeasured signal gets\n")
cat("  L = 0 by design, and two fits agree perfectly on it without any\n")
cat("  estimation having happened. `dropped_frac` is the share of total |SHAP|\n")
cat(sprintf("  sent to intervention groups and discarded under `%s`.\n", iv_handling))

# --- descriptive: how far apart are the methods on ONE fit? ------------------
# NOT the arm's result. The arm compares ONE method across TWO fits; this
# compares methods on one fit, and it is here so that the floors have context
# and so a metric stuck at 1 could not pass unnoticed.
cat("\n=== method against method, on the same fit (orientation only) ===\n\n")
# The two `llr_cond` rows were added 2026-09-05 with the arm. `llr_sum` against
# `llr_cond` is the pure propensity block seen as an attribution difference, and
# `llr_cond` against `xgb_l` is the aggregator contrast asked on the arm the
# interpretability claim is actually about.
pairs <- list(c("llr_sum", "xgb_l"), c("llr_sum", "xgb_feat"),
              c("llr_sum", "xgb_raw"), c("xgb_l", "xgb_raw"),
              c("llr_sum", "llr_meas"),
              c("llr_sum", "llr_cond"), c("llr_cond", "xgb_l"))
ma <- do.call(rbind, lapply(pairs, function(p) {
  if (!all(p %in% names(ref$A))) return(NULL)
  a <- attribution_agreement(ref$A[[p[1]]], ref$A[[p[2]]], cells$keep,
                             mag_floor = mag_floor)
  cbind(a = p[1], b = p[2], a$overall)
}))
print(ma, row.names = FALSE)
save_table(run, ma, "method_agreement")

# --- metric five: domain composition ----------------------------------------
cat("\n=== domain composition within risk quartiles (llr_sum) ===\n\n")
dom <- load_domains(cfg_req(cfg_local, "paths", "domains"))
dc <- attribution_domains(ref$A$llr_sum, rowSums(ref$A$llr_sum), dom,
                          n_bins = domain_bins)
print(utils::head(dc[dc$bin == max(dc$bin), ], 6), row.names = FALSE)
save_table(run, dc, "domain_composition")
cat("\n  Shown: the top quartile only. The full table is in the run directory.\n")

# --- the floors --------------------------------------------------------------
fs <- list(); fo <- list(); fp <- list(); fr <- list(); ff <- list(); ft <- list()
for (f in floors) {
  cat(sprintf("\n=== floor: %s, B = %d ===\n\n", f$perturbation, f$B))
  fl <- attribution_floor(
    bundle, tabs, cfg, eval_ids = ids, folds = folds,
    perturbation = f$perturbation, B = as.integer(f$B), methods = methods,
    mag_floor = mag_floor, seed = floor_seed,
    intervention_handling = iv_handling, min_signals = min_signals,
    l_source = l_source, verbose = TRUE)

  su <- floor_summary(fl)
  if (nrow(su)) {
    su <- cbind(perturbation = f$perturbation, su)
    print(su, row.names = FALSE)
  } else {
    cat("  NO REPLICATE SUCCEEDED. The failure table below is the result.\n")
  }
  log_msg(run, sprintf("%s: %d attempted, %d failed, %.1f min elapsed",
                       f$perturbation, fl$n_attempted, fl$n_failed, fl$elapsed_min))

  # TAGGED ONLY IF NON-EMPTY. `cbind(perturbation = "x", NULL)` returns a
  # one-row frame out of nothing, which would put a fictitious row into the
  # export of a floor where every replicate failed (finding F20).
  tag <- function(d) if (is.null(d) || !nrow(d)) NULL else
    cbind(perturbation = f$perturbation, d)
  if (nrow(su)) fs[[length(fs) + 1L]] <- su
  fo[[length(fo) + 1L]] <- tag(fl$overall)
  fp[[length(fp) + 1L]] <- tag(fl$per_signal)
  fr[[length(fr) + 1L]] <- tag(fl$ranking)
  ff[[length(ff) + 1L]] <- tag(fl$failures)
  ft[[length(ft) + 1L]] <- tag(fl$timing)

  if (fl$n_failed > 0L) {
    cat(sprintf("\n  %d of %d replicates FAILED and are counted, not dropped.\n",
                fl$n_failed, fl$n_attempted))
    cat("  A resample changes the size of every signal's measured subset, and a\n")
    cat("  thin signal can fall below its declared basis dimension. A floor\n")
    cat("  computed only over the replicates that worked is conditioned on the\n")
    cat("  estimator having worked, which is what makes it too narrow.\n")
  }
}

if (length(floors)) {
  # The failure and timing tables are written even when nothing succeeded,
  # because in that case they ARE the result (finding F20).
  bind <- function(l) { l <- Filter(Negate(is.null), l)
                        if (length(l)) do.call(rbind, l) else NULL }
  if (length(fs)) save_table(run, do.call(rbind, fs), "floor_summary")
  if (!is.null(bind(fo))) save_table(run, bind(fo), "floor_overall")
  if (!is.null(bind(fr))) save_table(run, bind(fr), "floor_ranking")
  if (!is.null(bind(fp)))
    save_table(run, bind(fp), "floor_per_signal", subdir = "diagnostics")
  save_table(run, bind(ft) %||%
    data.frame(perturbation = character(0), replicate = character(0),
               status = character(0), elapsed_min = numeric(0),
               cpu_min = numeric(0)),
    "floor_timing", subdir = "diagnostics")
  save_table(run, bind(ff) %||%
    data.frame(perturbation = character(0), replicate = character(0),
               message = character(0)),
    "floor_failures", subdir = "diagnostics")

  cat("\n=== how to read this ===\n\n")
  cat("  The floor is the distribution of the agreement statistic between two\n")
  cat("  fits of the SAME method on perturbed training data from ONE site. A\n")
  cat("  cross-site number is read against it in one sentence: the MIMIC-eICU\n")
  cat("  sign agreement is x, the within-site floor is y with interquartile\n")
  cat("  range a to b, so the site effect is or is not larger than what one\n")
  cat("  site's own sampling produces.\n\n")
  cat("  The DISJOINT HALF is the correct null for a transport claim, because\n")
  cat("  MIMIC and eICU share no patients either. The BOOTSTRAP floor sits too\n")
  cat("  high: two resamples share about 63% of their distinct patients, so\n")
  cat("  reading a cross-site number against it biases toward declaring a site\n")
  cat("  effect that is really sample-to-sample variation.\n")
} else {
  cat("\n  No floor was enabled. Set `run: true` in config/attribution.yml.\n")
}

write_manifest(run, "ok", note = note, extra = list(
  bundle_path = bundle_path,
  eval_split = eval_split, eval_n = length(ids), eval_seed = eval_seed,
  methods = as.list(methods),
  intervention_handling = iv_handling, replicate_l_source = l_source,
  mag_floor = mag_floor, min_signals = min_signals,
  floors_run = as.list(vapply(floors, function(f) f$perturbation, character(1))),
  floor_B = as.list(vapply(floors, function(f) as.integer(f$B), integer(1))),
  floor_seed = floor_seed,
  frac_cells_used = round(cells$frac_used, 5)))

cat(sprintf("\n  run directory: %s\n\n", run$path))
