# tests/attr_store_migrate.R --------------------------------------------------
# MIGRATE THE ATTRIBUTION REPLICATE STORE ACROSS THE 2026-09-07 SPEC CHANGE,
# WITH THE EVIDENCE ATTACHED. Fits nothing. Runs in seconds.
#
# --- WHAT HAPPENED, AND WHY THIS IS NOT A REGENERATION -----------------------
#
# `full_ti_trend` and `full_ti_all` joined `LAYER1_MODELS`, so `layer1_jobs()`
# now enumerates 64 fitted specs where it enumerated 43. `attr_design_key()`
# hashes every one of those formulas into its `formulas` field, so the key
# changed, so `tests/attr_metrics.R`'s staleness guard refuses a store built
# before the change -- all 678 files of it, including the 40 shared-bag LLR
# bootstrap replicates that cost about five hours.
#
# THE STORED L VALUES ARE UNCHANGED, and the reason is structural rather than
# hopeful. Adding a spec does not alter any OTHER spec's formula, its fitting
# rows, its folds or its `bam` settings; `p_bar_train` is per (signal, fold) and
# not per model; layer 2 and Sigma see signals, not alternative specs of one
# signal; and the interaction arms in the store were fitted from
# `update(build_formula(sg, "full", cfg), . ~ . + <ti>)`, which
# `build_formula(sg, "full_ti_*", cfg)` now reproduces term for term. So the
# `formulas` field reports a difference that is TRUE about the design and FALSE
# about the data.
#
# That makes the field over-broad in one direction and it stays that way on
# purpose: narrowing it to "the formulas of the specs this store holds" would
# stop it catching a spec being REMOVED, which is a change that does invalidate
# a store. The remedy for over-breadth is not a looser key, it is a migration
# that PROVES the claim the key cannot express.
#
# --- WHAT THIS SCRIPT PROVES -------------------------------------------------
#
#   1. Every field of the design key EXCEPT `formulas` is unchanged. If `bam`,
#      `signals`, `paired`, `ti_all`, `ti_trend`, `k_ti` or `folds` differs, the
#      store really is stale and the migration stops.
#   2. `formulas` differs ONLY by ADDITION. Every spec the old key described is
#      present in the new key with a byte-identical formula. A spec whose
#      formula CHANGED, or one that disappeared, stops the migration.
#   3. The eight ladder replicates on disk -- the anchor `(0, 0, 0)` matrices
#      the whole level-4 comparison is read against -- are BITWISE identical to
#      what the pipeline now produces from `l_mats_zero`. Not close, identical:
#      `max(abs(stored - live))` must be exactly 0.
#
# (3) is the load-bearing one. The 634 expensive replicates were produced by the
# same code from the same formulas on resampled rows, so a bitwise match on the
# unresampled anchor is the strongest available evidence that they are equally
# valid short of regenerating them, which would take five hours to reproduce
# numbers this asserts are already right.
#
# The pattern is `tests/refactor_identity.R`'s: assert that a refactor changed
# no number, and record the assertion. It is the standard this project holds
# elsewhere and it is what makes a key rewrite a migration rather than a
# `--force`.
#
# --- WHAT IT WRITES ----------------------------------------------------------
#
# On success it rewrites `<store>/design.qs2` with the new key and writes
# `<store>/migration_20260907.csv` recording every comparison and its result,
# so the store carries its own evidence. `design_pre20260907.qs2` keeps the old
# key beside it -- the same convention `manifest_replicates_pre20260906.csv`
# already uses -- because a migration that destroys what it migrated from
# cannot be audited afterwards.
#
#   Rscript tests/attr_store_migrate.R                 # dry run, changes nothing
#   Rscript tests/attr_store_migrate.R --write         # rewrite the key
#   Rscript tests/attr_store_migrate.R --store <dir>   # a store other than the
#                                                      # one config names
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(targets); library(qs2); library(yaml); library(digest)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

.args  <- commandArgs(trailingOnly = TRUE)
.opt   <- function(flag, default = NA_character_) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[i + 1L]
}
WRITE <- "--write" %in% .args

ecfg  <- yaml::read_yaml("config/attribution_eval.yml")
# RESOLVED EXACTLY AS `tests/attr_metrics.R` RESOLVES IT: a positional argument
# first, `latest_run("attrgen")` otherwise. Two scripts that disagree about
# which store is "the" store would migrate one and read another, and nothing
# would say so.
.pos  <- setdiff(.args, c("--write", "--store"))
.pos  <- .pos[!.pos %in% .opt("--store", character(0))]
STORE <- .opt("--store", if (length(.pos)) .pos[1] else NA_character_)
if (is.na(STORE) || !nzchar(STORE)) {
  STORE <- latest_run("attrgen", require_complete = FALSE)
}
if (is.null(STORE) || is.na(STORE) || !dir.exists(STORE)) {
  stop("no attrgen store found. Pass the run directory as the first argument ",
       "or with --store.", call. = FALSE)
}

cat("\n=== attribution store migration, 2026-09-07 spec change ===\n\n")
cat(sprintf("  store : %s\n", STORE))
cat(sprintf("  mode  : %s\n\n", if (WRITE) "WRITE" else "dry run (pass --write to commit)"))

# --- the two keys ------------------------------------------------------------
cfg    <- tar_read(cfg)
folds  <- tar_read(folds)
tr     <- tar_read(train_ids)
fold_k <- folds$fold[match(tr, folds$stay_id)]

dg <- qs2::qs_read(file.path(STORE, "design.qs2"))
if (is.null(dg$design_key)) {
  stop("the store carries no `design_key` and predates 2026-09-06. It cannot ",
       "be migrated, only regenerated.", call. = FALSE)
}
if (is.null(dg$k_ti)) {
  stop("the store carries a design but no `k_ti`; it was written by an ",
       "inconsistent version of tests/attr_replicates.R.", call. = FALSE)
}
if (!identical(as.character(dg$stay_id), as.character(tr))) {
  stop("the store's `stay_id` vector is not the current training set. That is ",
       "a cohort change, not a spec change, and nothing here can migrate it.",
       call. = FALSE)
}

OLD <- dg$design_key
# KEYED TO THE STORE'S OWN `k_ti`, not to today's config value. They are equal
# at present -- the move of `k_ti` into config/config.yml deliberately kept the
# value at 5 -- and asserting that separately below is a stronger statement than
# assuming it here.
NEW <- attr_design_key(cfg, fold_k, dg$k_ti)

d <- attr_design_diff(OLD, NEW)
cat("=== 1. which design fields differ ===\n\n")
if (!nrow(d)) {
  cat("  none. The store is already current and there is nothing to migrate.\n\n")
} else {
  print(d[, c("field", "hash_a", "hash_b")], row.names = FALSE)
}

other <- d$field[d$field != "formulas"]
if (length(other)) {
  cat("\n")
  abort_values(paste0("field(s) besides `formulas` differ, so the store is ",
                      "genuinely stale and the migration must not proceed"), other)
}
cat(sprintf("\n  `k_ti` in the store: %d; in config/config.yml: %d\n",
            as.integer(dg$k_ti), as.integer(cfg_req(cfg, "k_ti"))))
if (!identical(as.integer(dg$k_ti), as.integer(cfg_req(cfg, "k_ti")))) {
  stop("the store's `k_ti` differs from the config's. The cross terms in the ",
       "store are not the cross terms the pipeline now fits.", call. = FALSE)
}

# --- 2. `formulas` differs only by addition ----------------------------------
cat("\n=== 2. is `formulas` a pure addition? ===\n\n")
fo <- OLD$formulas; fn <- NEW$formulas
gone    <- setdiff(names(fo), names(fn))
added   <- setdiff(names(fn), names(fo))
shared  <- intersect(names(fo), names(fn))
changed <- shared[!vapply(shared, function(k) identical(fo[[k]], fn[[k]]), logical(1))]
cat(sprintf("  specs before        : %d\n", length(fo)))
cat(sprintf("  specs now           : %d\n", length(fn)))
cat(sprintf("  added               : %d\n", length(added)))
cat(sprintf("  removed             : %d\n", length(gone)))
cat(sprintf("  shared but CHANGED  : %d\n", length(changed)))
if (length(gone))    abort_values("spec(s) REMOVED from the design; a store keyed to them is stale", gone)
if (length(changed)) abort_values("spec(s) whose formula CHANGED; the stored L for these is stale", changed)
cat(sprintf("\n  added specs: %s\n", paste(sort(added), collapse = ", ")))

# --- 3. the bitwise check ----------------------------------------------------
cat("\n=== 3. do the stored ladder replicates equal the pipeline's L? ===\n\n")
cat("  The anchor (boot 0, seed 0, draw 0) of every LLR arm, against\n")
cat("  `l_mats_zero` as the graph now computes it. Exact equality required.\n\n")

Lz   <- tar_read(l_mats_zero)
sigs <- as.character(unlist(cfg$signals))
need <- ATTR_BASE_ARMS
miss <- setdiff(need, names(Lz))
if (length(miss)) {
  abort_values(paste0("`l_mats_zero` is missing L matri(ces) -- the interaction ",
                      "models must be in LAYER1_MODELS and the graph must be up ",
                      "to date. Run targets::tar_make()"), miss)
}
LIVE <- attr_derive_arms(Lz[need])

man <- utils::read.csv(file.path(STORE, "manifest_replicates.csv"),
                       stringsAsFactors = FALSE)
SUB <- as.character(cfg_req(ecfg, "storage", "subdir"))

rows <- list()
for (m in paste0("llr_", ATTR_LLR_ARMS)) {
  z <- man[man$method == m & man$route == "ladder", , drop = FALSE]
  if (!nrow(z)) {
    rows[[length(rows) + 1L]] <- data.frame(
      method = m, present = FALSE, max_abs_diff = NA_real_, identical = NA,
      stringsAsFactors = FALSE)
    cat(sprintf("  %-24s  no ladder replicate in the manifest\n", m))
    next
  }
  if (nrow(z) > 1L) {
    stop("more than one `ladder` replicate for `", m, "`; the manifest is not ",
         "content-addressed as it claims.", call. = FALSE)
  }
  M <- qs2::qs_read(file.path(STORE, SUB, paste0(z$key[1], ".qs2")))
  L <- LIVE[[attr_arm_of(m)]][, sigs, drop = FALSE]
  if (!identical(dim(M), dim(L))) {
    stop("shape mismatch for `", m, "`: stored ", paste(dim(M), collapse = "x"),
         ", live ", paste(dim(L), collapse = "x"), call. = FALSE)
  }
  if (!identical(colnames(M), colnames(L))) {
    stop("column order differs for `", m, "`; the comparison would be between ",
         "different signals.", call. = FALSE)
  }
  dm <- max(abs(M - L))
  rows[[length(rows) + 1L]] <- data.frame(
    method = m, present = TRUE, max_abs_diff = dm, identical = (dm == 0),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-24s  max |stored - live| = %.3e   %s\n", m, dm,
              if (dm == 0) "IDENTICAL" else "*** DIFFERS ***"))
}
CHK <- do.call(rbind, rows)

bad <- CHK$method[isTRUE(any(CHK$present)) & !is.na(CHK$identical) & !CHK$identical]
if (length(bad)) {
  cat("\n")
  abort_values(paste0("the pipeline's L differs from the stored ladder for ",
                      "arm(s) below. The 634 expensive replicates were built ",
                      "from the same code on resampled rows, so a difference ",
                      "here means the store no longer describes the design. ",
                      "Regenerate rather than migrate"), bad)
}
if (!any(CHK$present, na.rm = TRUE)) {
  stop("no ladder replicate was found for any arm, so there is nothing to ",
       "verify against and the migration has no evidence. Refusing.",
       call. = FALSE)
}

n_ok <- sum(CHK$identical, na.rm = TRUE)
cat(sprintf("\n  %d of %d LLR arms verified bitwise identical.\n",
            n_ok, nrow(CHK)))
cat(sprintf("  %d replicate file(s) in the store are covered by this migration.\n",
            nrow(man)))

# --- commit ------------------------------------------------------------------
cat("\n=== verdict ===\n\n")
if (!WRITE) {
  cat("  DRY RUN. Every check passed. Re-run with --write to rewrite\n")
  cat(sprintf("  %s and record the evidence.\n\n", file.path(STORE, "design.qs2")))
  quit(save = "no", status = 0L)
}

file.copy(file.path(STORE, "design.qs2"),
          file.path(STORE, "design_pre20260907.qs2"), overwrite = FALSE)
# `naming_design_key` IS THE KEY THE FILES STAY NAMED UNDER, and it was
# missing from the first migration (found 2026-09-09): `attr_replicate_key()`
# hashes the whole design key into every file name, so rewriting `design_key`
# alone left the generator unable to find a single existing replicate on its
# next resume. The files are not renamed here; the generator hashes names
# under this key and checks the live design against `design_key`. A store
# that already carries a naming key keeps it -- a second migration must not
# move it to the first migration's key.
qs2::qs_save(list(design_key = NEW, stay_id = tr, k_ti = as.integer(dg$k_ti),
                  naming_design_key = dg$naming_design_key %||% OLD,
                  fingerprint = dg$fingerprint,
                  fingerprint_evidence = dg$fingerprint_evidence,
                  migrated_from = attr_key_hash(OLD),
                  migrated_on = "2026-09-07",
                  migration = "full_ti_trend/full_ti_all entered LAYER1_MODELS"),
             file.path(STORE, "design.qs2"))
CHK$old_key <- attr_key_hash(OLD)
CHK$new_key <- attr_key_hash(NEW)
CHK$n_specs_added <- length(added)
CHK$n_replicates_covered <- nrow(man)
utils::write.csv(CHK, file.path(STORE, "migration_20260907.csv"), row.names = FALSE)

cat(sprintf("  design.qs2 rewritten: %s -> %s\n",
            attr_key_hash(OLD), attr_key_hash(NEW)))
cat(sprintf("  old key kept at      %s\n",
            file.path(STORE, "design_pre20260907.qs2")))
cat(sprintf("  evidence written to  %s\n\n",
            file.path(STORE, "migration_20260907.csv")))
