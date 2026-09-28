# tests/invariance_suite.R ---------------------------------------------------
# THE INVARIANCE SUITE of docs/v2_audit_plan_20260903.md section 3, items 1, 3
# and 5. Audits A to F find today's drift; these assertions prevent tomorrow's.
#
# Item 2 is Audit A and lives in `tests/audit_a_invalidation.R`. Item 4 is
# already implemented in `tests/attribution_smoke.R`. Item 6 is deferred to
# section 6 step 5, where the plan puts it.
#
# ITEM 1 — DETERMINISM. Two consecutive `tar_make()` runs on an unchanged tree
#   must produce bitwise-identical values for every target. The FULL test costs
#   a complete rebuild and cannot be run as a side effect of an audit, so it is
#   split:
#     1a  the necessary condition, free: `tar_outdated()` is empty on an
#         unchanged tree. If this fails, something in the graph is not stable
#         under a no-op and nothing further need be measured.
#     1b  the sufficient test, as a two-phase harness: `--snapshot` records
#         every target's data hash from `tar_meta()`; after a rebuild,
#         `--compare` reports every target whose hash moved. The script never
#         calls `tar_destroy()` itself. Deliberately: destroying a store that
#         took hours to build must be a decision somebody makes, not something
#         a test does on the way past.
#
# ITEM 3 — APPLY-PATH CONSISTENCY. The claim under test is that the fitting
#   site and an apply site traverse the same code path. Applying the bundle to
#   TRAINING rows must reproduce what the graph's own objects produce on those
#   same rows. A divergence here would silently corrupt both the eICU result
#   and every attribution, and it would look like a finding rather than a bug.
#
#   WHAT IS ACTUALLY COMPARED, because the naive comparison is empty. The apply
#   path (`apply_one`) and the fitting path (`fit_one`'s predict branch) call
#   the same `signal_frame(stage = "predict")` and the same `predict.gam`, so
#   comparing those two expressions proves nothing. What CAN differ is the
#   OBJECTS they are handed:
#
#     cfg      the graph uses `load_config()`; an apply site uses
#              `bundle_cfg(bundle, paths)`, the FROZEN design. Every level
#              rule, tail rule, basis dimension and conditional-prior switch
#              the frame builder reads comes from there.
#     priors   the graph uses the full container including fold rows; the
#              bundle carries `priors_final()`, a restriction.
#     models   the graph holds live gam objects; the bundle holds
#              `strip_gam()`ed ones.
#     p_bar    the graph's per-signal training prior against the bundle's.
#
#   So the test drives BOTH paths over the same training stays -- one from the
#   graph's objects, one from the bundle's -- and asserts the L values are
#   identical, spec by spec. That is the comparison with content in it.
#
#   `--refit N` additionally re-fits N specs from scratch and checks the fresh
#   fit predicts identically to the stripped one in the bundle, which is the
#   one thing the object comparison cannot reach.
#
# ITEM 5 — DESIGN-HASH AGREEMENT. `bundle_cfg()` must be the only route to a
#   scoreable cfg, and the frozen design must be exactly `BUNDLE_DESIGN_KEYS`
#   of the live config -- no more (a leaked `paths` would make the hash
#   site-dependent) and no less (a missing key leaves a hole an apply site
#   fills from somewhere else). Checked as a hash, as a key-set comparison in
#   both directions, as a static count of where the stamp is written, and
#   empirically by handing `apply_bundle()` a config it must refuse.
#
# AGGREGATES ONLY (hard rule 1). L values are compared in bulk; the maxima and
# counts are printed and no stay identifier or row ever is.
#
#   Rscript tests/invariance_suite.R                  items 1a, 3, 5   ~10 min
#   Rscript tests/invariance_suite.R --items=5        one item
#   Rscript tests/invariance_suite.R --items=3 --refit=2
#   Rscript tests/invariance_suite.R --snapshot       item 1b, phase 1
#   Rscript tests/invariance_suite.R --compare        item 1b, phase 2
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
  library(targets)
})
root <- "."
for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)
source(file.path(root, "tests", "audit_common.R"))

args   <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^", flag, "="), "", hit[1])
}
items   <- strsplit(argval("--items", "1,3,5"), ",")[[1]]
n_refit <- as.integer(argval("--refit", "0"))
SNAP    <- file.path(root, "tests", "invariance_determinism_snapshot.csv")

res <- list()
note <- function(id, claim, ok, detail = "") {
  res[[length(res) + 1L]] <<- data.frame(
    item = id, claim = claim,
    verdict = if (isTRUE(ok)) "PASS" else if (is.na(ok)) "SKIP" else "FAIL",
    detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%-4s] %-6s %s%s\n", id,
              if (isTRUE(ok)) "PASS" else if (is.na(ok)) "SKIP" else "FAIL",
              claim, if (nzchar(detail)) paste0("  --  ", detail) else ""))
}

# --- item 1b phases, which short-circuit everything else ---------------------

if ("--snapshot" %in% args || "--compare" %in% args) {
  m <- targets::tar_meta(fields = c("name", "data", "type"))
  m <- m[m$type %in% c("stem", "branch", "pattern"), c("name", "data")]
  m <- m[order(m$name), ]
  if ("--snapshot" %in% args) {
    utils::write.csv(m, SNAP, row.names = FALSE)
    cat(sprintf("\nitem 1b: snapshot of %d target hashes written to %s\n", nrow(m), SNAP))
    cat("Now rebuild (`targets::tar_destroy(); targets::tar_make()`), then re-run\n")
    cat("this script with --compare. Any target whose hash moved is a hidden\n")
    cat("clock, RNG or ordering dependence.\n")
  } else {
    if (!file.exists(SNAP)) stop("no snapshot at ", SNAP, "; run --snapshot first",
                                 call. = FALSE)
    old <- utils::read.csv(SNAP, colClasses = "character")
    j <- merge(old, m, by = "name", suffixes = c("_before", "_after"), all = TRUE)
    # A target present in one snapshot and not the other counts as moved: it
    # means the graph itself changed between the two runs, which invalidates
    # the comparison rather than passing it.
    j$moved <- is.na(j$data_before) | is.na(j$data_after) |
               (j$data_before != j$data_after)
    cat(sprintf("\nitem 1b: %d target(s) compared, %d moved\n", nrow(j), sum(j$moved)))
    if (any(j$moved)) print(j[j$moved, c("name", "data_before", "data_after")],
                            row.names = FALSE)
  }
  quit(save = "no")
}

cat("\n=== INVARIANCE SUITE ===\n")

# --- item 1a -----------------------------------------------------------------

if ("1" %in% items) {
  cat("\nitem 1 — determinism\n")
  o <- targets::tar_outdated(callr_function = NULL)
  note("1a", "tar_outdated() is empty on an unchanged tree", length(o) == 0L,
       if (length(o)) paste0(length(o), " outdated: ",
                             paste(utils::head(o, 10), collapse = ", ")) else "")
  note("1b", "full two-run bitwise comparison", NA,
       sprintf("harness only; run --snapshot, rebuild, then --compare. Snapshot %s.",
               if (file.exists(SNAP)) "exists" else "not yet taken"))
}

# --- shared load, for items 3 and 5 -----------------------------------------

need_bundle <- any(c("3", "5") %in% items)
bundle <- NULL
if (need_bundle) {
  rc <- yaml::read_yaml(file.path(root, "config", "internal.yml"))
  # `--bundle=` overrides `test_look.bundle`. The config path is explicit on
  # purpose -- a published result must not change because a newer run appeared
  # -- so a freshly built bundle must be testable WITHOUT editing it, or the
  # test would force the very edit it is supposed to inform.
  bundle_path <- argval("--bundle", rc$test_look$bundle)
  if (is.null(bundle_path) || !file.exists(bundle_path)) {
    note("3/5", "bundle available", FALSE,
         paste0("no bundle at ", bundle_path %||% "<unset>",
                "; set test_look.bundle in config/internal.yml"))
    items <- setdiff(items, c("3", "5"))
  } else {
    cfg_local <- load_config(rc$config)
    bundle    <- load_bundle(bundle_path, cfg = cfg_local, strict = FALSE, verbose = FALSE)
    cat("\nbundle: ", bundle_path, "\n", sep = "")
  }
}

# --- item 5 ------------------------------------------------------------------

if ("5" %in% items && !is.null(bundle)) {
  cat("\nitem 5 — design-hash agreement\n")

  # 5a. The frozen design IS the live config restricted to the design keys.
  live_design <- cfg_local[intersect(BUNDLE_DESIGN_KEYS, names(cfg_local))]
  note("5a", "hash(bundle$cfg) == hash(local config restricted to design keys)",
       identical(.hash(bundle$cfg), .hash(live_design)),
       sprintf("bundle %s, local %s", .hash(bundle$cfg), .hash(live_design)))

  # 5b/5c. Key set, both directions. A leaked key would make the hash
  # site-dependent; a missing key leaves a hole an apply site fills elsewhere.
  note("5b", "every BUNDLE_DESIGN_KEY is present in the frozen design",
       !length(setdiff(BUNDLE_DESIGN_KEYS, names(bundle$cfg))),
       paste("missing:", paste(setdiff(BUNDLE_DESIGN_KEYS, names(bundle$cfg)),
                               collapse = ", ")))
  extra <- setdiff(names(bundle$cfg), BUNDLE_DESIGN_KEYS)
  note("5c", "the frozen design carries NOTHING beyond the design keys",
       !length(extra), paste("extra:", paste(extra, collapse = ", ")))

  # 5d. Per-key equality, so a failure of 5a names the key rather than a hash.
  diffs <- unlist(Filter(Negate(is.null), lapply(BUNDLE_DESIGN_KEYS, function(k)
    if (isTRUE(all.equal(bundle$cfg[[k]], cfg_local[[k]]))) NULL else k)))
  note("5d", "every design key is value-identical to the live config",
       !length(diffs), paste("differs in:", paste(diffs, collapse = ", ")))

  # 5e. `bundle_cfg()` is the ONLY writer of the stamp. Static, over comment-
  # stripped source, so a stamp assigned in a second place is a failure even if
  # it happens to compute the same value today.
  #
  # SCOPED TO THE NON-TEST LAYERS, and that scoping is the point rather than a
  # convenience: check 5h below forges a stamp deliberately, so a search that
  # included `tests/` would report this file and the check would fail on its own
  # evidence. A test forging a stamp is how the guard gets exercised; a LIBRARY
  # or RUNNER forging one is the failure.
  src <- scan_sources(root)
  wr_all <- src[grepl("[$]\\s*\\.bundle_design\\s*<-", src$code, perl = TRUE), ]
  writers <- wr_all[wr_all$layer != "tests", , drop = FALSE]
  in_tests <- wr_all[wr_all$layer == "tests", , drop = FALSE]
  note("5e", "`.bundle_design` is assigned in exactly one place outside tests/",
       nrow(writers) == 1L && grepl("10_bundle[.]R$", writers$path[1]),
       sprintf("%s%s",
               paste(paste0(sub("^[.]/", "", writers$path), ":", writers$line), collapse = ", "),
               if (nrow(in_tests))
                 paste0("  (plus ", nrow(in_tests), " deliberate forgery/forgeries in tests/)")
               else ""))

  # 5f. Every `cfg$<key>` R/ reads is either frozen or a documented exclusion.
  refs <- cfg_refs(src, "lib")
  unfrozen <- setdiff(refs, c(BUNDLE_DESIGN_KEYS, "paths", "pairing"))
  note("5f", "every cfg key R/ reads is frozen (paths and pairing excepted)",
       !length(unfrozen), paste("unfrozen:", paste(unfrozen, collapse = ", ")))

  # 5g. `apply_bundle()` must REFUSE a config that did not come from
  # `bundle_cfg()`. Checked by doing it, because the guard is the load-bearing
  # part of hard rule 8 and a guard nobody has tripped is a guard nobody knows
  # works.
  refused <- tryCatch({
    apply_bundle(bundle, tabs = NULL, cfg = cfg_local, stay_ids = character(0),
                 arms = "llr_sum", verbose = FALSE)
    FALSE
  }, error = function(e) grepl("did not come from bundle_cfg", conditionMessage(e)))
  note("5g", "apply_bundle() refuses a cfg that did not come from bundle_cfg()",
       isTRUE(refused))

  # 5h. And it must refuse a stamp forged from a DIFFERENT bundle's design.
  forged <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
  forged$.bundle_design <- .hash(list(not = "this bundle"))
  refused2 <- tryCatch({
    apply_bundle(bundle, tabs = NULL, cfg = forged, stay_ids = character(0),
                 arms = "llr_sum", verbose = FALSE)
    FALSE
  }, error = function(e) grepl("different bundle", conditionMessage(e)))
  note("5h", "apply_bundle() refuses a stamp from a different bundle",
       isTRUE(refused2))
}

# --- item 3 ------------------------------------------------------------------

if ("3" %in% items && !is.null(bundle)) {
  cat("\nitem 3 — apply-path consistency on TRAINING rows\n")

  cfg_bundle <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
  tabs <- load_tables(cfg_bundle$paths, cfg_bundle, site = "mimic", verbose = FALSE)

  # The graph's own objects. `tar_read()` values are row-level and go nowhere
  # near a console (hard rule 1).
  priors_graph <- targets::tar_read(priors, store = file.path(root, "_targets"))
  models_graph <- targets::tar_read(final_models, store = file.path(root, "_targets"))
  folds_graph  <- targets::tar_read(folds, store = file.path(root, "_targets"))
  train_ids    <- folds_graph$stay_id[folds_graph$split == "train"]

  # The fold partition an apply site would compute for itself, from the FROZEN
  # seed and split, must be the partition the graph used. If this fails, every
  # comparison below is on different rows and nothing else is interpretable.
  folds_apply <- assign_folds(tabs$cohort, cfg_bundle)
  same_split <- identical(folds_apply$split, folds_graph$split) &&
                identical(folds_apply$stay_id, folds_graph$stay_id)
  note("3a", "the frozen seed/split reproduces the graph's train/test partition",
       same_split,
       sprintf("%d train stays", length(train_ids)))

  # p_bar, per signal, bundle against graph.
  pg <- priors_final(priors_graph)$signal
  pb <- bundle$priors$signal
  m  <- merge(pg[, c("signal", "p_bar")], pb[, c("signal", "p_bar")],
              by = "signal", suffixes = c("_graph", "_bundle"))
  dp <- max(abs(m$p_bar_graph - m$p_bar_bundle))
  note("3b", "p_bar_train is identical in the bundle and in the graph",
       nrow(m) == nrow(pg) && dp == 0, sprintf("max |diff| = %.3e over %d signals", dp, nrow(m)))

  # The L's themselves. Two full apply passes over the training stays: one
  # driven by the graph's cfg/priors/models, one by the bundle's.
  cat("  running the graph-driven apply pass...\n")
  a_graph <- apply_layer1(models_graph, tabs, cfg_local,
                          priors_final(priors_graph), train_ids, verbose = FALSE)
  cat("  running the bundle-driven apply pass...\n")
  a_bund  <- apply_layer1(bundle$models, tabs, cfg_bundle,
                          bundle$priors, train_ids, verbose = FALSE)

  kg <- paste(a_graph$l$signal, a_graph$l$model, a_graph$l$stay_id)
  kb <- paste(a_bund$l$signal,  a_bund$l$model,  a_bund$l$stay_id)
  note("3c", "the two passes score the identical (signal, model, stay) set",
       setequal(kg, kb),
       sprintf("graph %d rows, bundle %d rows, %d only in one",
               length(kg), length(kb), length(union(setdiff(kg, kb), setdiff(kb, kg)))))

  if (setequal(kg, kb)) {
    lb <- a_bund$l$l[match(kg, kb)]
    dl <- max(abs(a_graph$l$l - lb))
    note("3d", "every L is bitwise identical between the two passes", dl == 0,
         sprintf("max |diff| = %.3e over %d L values", dl, length(kg)))

    # Per-spec, so a failure names the model rather than the pipeline.
    per <- do.call(rbind, lapply(split(seq_along(kg),
                                       paste(a_graph$l$signal, a_graph$l$model)),
      function(ix) data.frame(spec = paste(a_graph$l$signal[ix[1]],
                                           a_graph$l$model[ix[1]]),
                              n = length(ix),
                              max_abs_diff = max(abs(a_graph$l$l[ix] - lb[ix])),
                              stringsAsFactors = FALSE)))
    bad <- per[per$max_abs_diff > 0, , drop = FALSE]
    if (nrow(bad)) print(bad[order(-bad$max_abs_diff), ], row.names = FALSE)
  }

  # Optional: a fresh fit against the stripped bundle model. This is the one
  # thing an object comparison cannot reach -- whether `strip_gam()` and the
  # store round trip left prediction unchanged.
  if (n_refit > 0L) {
    jobs <- layer1_jobs(cfg_bundle)
    jobs <- jobs[jobs$fit & jobs$role == "final", , drop = FALSE]
    # Spread across the job table rather than taking the first N: the specs
    # differ in whether they are paired, whether they carry an ordinal `delta`,
    # and whether `full` is an alias of `meas`, and taking the first N would
    # sample one corner of that.
    ix <- unique(round(seq(1, nrow(jobs), length.out = min(n_refit, nrow(jobs)))))
    pick <- jobs[ix, , drop = FALSE]
    rows <- list()
    for (i in seq_len(nrow(pick))) {
      sg <- pick$signal[i]; md <- pick$model[i]
      key <- paste(sg, md, sep = "/")
      cat(sprintf("  refitting %s ...\n", key))
      pri <- priors_for(priors_final(priors_graph), sg, "final", NA_integer_)
      fr  <- fit_one(sg, md, tabs, cfg_local, pri,
                     fit_ids = train_ids, predict_ids = train_ids,
                     keep_model = FALSE, role = "final", fold = NA_integer_)
      ap  <- apply_one(bundle$models[[key]], sg, md, tabs, cfg_bundle,
                       priors_for(bundle$priors, sg, "final", NA_integer_), train_ids)
      j <- match(fr$l$stay_id, ap$stay_id)
      rows[[i]] <- data.frame(spec = key, n = nrow(fr$l),
                              max_abs_diff = if (anyNA(j)) NA_real_
                                             else max(abs(fr$l$l - ap$l[j])),
                              stringsAsFactors = FALSE)
    }
    rt <- do.call(rbind, rows)
    note("3e", sprintf("a fresh fit predicts as the bundle's stripped model (%d spec(s))",
                       nrow(rt)),
         all(!is.na(rt$max_abs_diff)) && max(rt$max_abs_diff) < 1e-9,
         sprintf("max |diff| = %.3e", suppressWarnings(max(rt$max_abs_diff, na.rm = TRUE))))
    print(rt, row.names = FALSE)
  }
}

# --- report ------------------------------------------------------------------

out <- do.call(rbind, res)
md <- c("# Invariance suite — items 1, 3, 5", "",
        "Generated by `tests/invariance_suite.R`. Aggregates only (hard rule 1).", "",
        md_table(out), "")
audit_write(md, "invariance_suite.md", root)

cat("\n=== SUMMARY ===\n")
print(out[, c("item", "verdict", "claim")], row.names = FALSE)
if (any(out$verdict == "FAIL")) {
  cat("\nFAILURES:\n")
  print(out[out$verdict == "FAIL", ], row.names = FALSE)
}
