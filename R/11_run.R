# R/11_run.R -----------------------------------------------------------------
# Run directories, manifests, and save conventions.
#
# This is the ONLY file in R/ permitted to construct a path or read the clock
# (CLAUDE.md hard rule 9). Every other library file receives a `run` object and
# writes through the helpers here. That restriction is what makes the MIMIC-test
# and eICU code paths provably identical: neither can differ by a path.
#
# A run directory is immutable once finalised. `targets` caches computation in
# _targets/ and overwrites by design; that is a cache. This is the output.
# ----------------------------------------------------------------------------

RUN_SUBDIRS <- c("tables", "figs", "diagnostics")

# --- creation ---------------------------------------------------------------

#' Create a new, unique run directory.
#'
#' @param prefix   short run kind, e.g. "internal", "external", "survival".
#' @param config   the resolved config list; snapshotted into the manifest.
#' @param root     parent directory for all runs.
#' @param note     optional free-text label recorded in the manifest.
#' @return a `run` object (list) to be threaded through the pipeline.
new_run <- function(prefix, config, root = "out/runs", note = NULL) {
  stopifnot(is.character(prefix), length(prefix) == 1L, nzchar(prefix))
  if (!grepl("^[a-z0-9_]+$", prefix)) {
    stop("`prefix` must be lowercase alphanumeric/underscore, got: ", prefix)
  }

  started <- Sys.time()
  stamp   <- format(started, "%Y%m%dT%H%M%S")

  # Two runs launched in the same second must not collide. Suffix rather than
  # overwrite: losing a previous run's outputs is the failure this guards.
  base <- file.path(root, paste0(prefix, "_", stamp))
  path <- base
  i <- 1L
  while (dir.exists(path)) {
    path <- paste0(base, "-", i)
    i <- i + 1L
    if (i > 100L) stop("could not allocate a unique run directory under ", root)
  }

  dir.create(path, recursive = TRUE)
  for (d in RUN_SUBDIRS) dir.create(file.path(path, d))

  run <- list(
    prefix    = prefix,
    id        = basename(path),
    path      = path,
    started   = started,
    config    = config,
    # KEPT ON THE RUN OBJECT (plumbing review F10, 2026-09-08). The note was
    # written into the "running" manifest and then dropped: `finalize_run()`
    # rewrote the manifest with its own `note = NULL`, so every completed
    # manifest in out/runs/ carries `note: ''` whatever the runner said.
    note      = note,
    log_file  = file.path(path, "log.txt")
  )
  class(run) <- "llr_run"

  # Written immediately, with status "running", so an interrupted run leaves a
  # directory that says so rather than one that merely looks incomplete.
  write_manifest(run, status = "running", note = note)
  log_msg(run, "run started: ", run$id)
  run
}

print.llr_run <- function(x, ...) {
  cat("<llr_run>", x$id, "\n  path: ", x$path, "\n", sep = "")
  invisible(x)
}

# --- paths ------------------------------------------------------------------

#' Path inside a run directory. The only path constructor in the library.
run_path <- function(run, ...) {
  stopifnot(inherits(run, "llr_run"))
  file.path(run$path, ...)
}

# --- logging ----------------------------------------------------------------

#' Append a timestamped line to the run log and echo it to the console.
log_msg <- function(run, ...) {
  stopifnot(inherits(run, "llr_run"))
  line <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(line, "\n", sep = "", file = run$log_file, append = TRUE)
  message(line)
  invisible(line)
}

# --- saving -----------------------------------------------------------------

#' Save a data frame to tables/, as .rds and (by default) a .csv alongside.
#'
#' The .rds is canonical — it round-trips types and factor levels. The .csv
#' exists purely so results can be eyeballed without starting R, and is never
#' read back by the pipeline.
save_table <- function(run, x, name, csv = TRUE, subdir = "tables") {
  stopifnot(inherits(run, "llr_run"), is.data.frame(x))
  if (!subdir %in% RUN_SUBDIRS) abort_values("unknown run subdirectory", subdir)
  name <- .check_name(name)
  saveRDS(x, run_path(run, subdir, paste0(name, ".rds")))
  if (csv) {
    utils::write.csv(x, run_path(run, subdir, paste0(name, ".csv")),
                     row.names = FALSE, na = "")
  }
  log_msg(run, "table  ", subdir, "/", name, "  [", nrow(x), " x ", ncol(x), "]")
  invisible(x)
}

#' Save an arbitrary R object to a named subdirectory as .rds.
save_object <- function(run, x, name, subdir = "tables") {
  stopifnot(inherits(run, "llr_run"))
  name <- .check_name(name)
  saveRDS(x, run_path(run, subdir, paste0(name, ".rds")))
  log_msg(run, "object ", subdir, "/", name)
  invisible(x)
}

#' Save the fitted-model bundle. qs2, because it carries ~190 gam objects.
#'
#' Only `run/internal.R` ever calls this — external and survival load a bundle
#' and fit nothing (CLAUDE.md hard rule 8). The path returned here is what goes
#' into config/external.yml.
save_bundle <- function(run, bundle, name = "bundle") {
  stopifnot(inherits(run, "llr_run"))
  if (run$prefix != "internal") {
    stop("only an internal run may write a bundle; this run is '", run$prefix,
         "'. External and survival runs fit nothing (hard rule 8).")
  }
  p <- run_path(run, paste0(.check_name(name), ".qs2"))
  qs2::qs_save(bundle, p)
  log_msg(run, "bundle ", basename(p), "  (",
          round(file.size(p) / 1024^2, 1), " MB)")
  invisible(p)
}

#' Open a PNG device writing into figs/. Caller closes with dev.off().
save_fig <- function(run, name, width = 7, height = 5, dpi = 150) {
  stopifnot(inherits(run, "llr_run"))
  p <- run_path(run, "figs", paste0(.check_name(name), ".png"))
  grDevices::png(p, width = width * dpi, height = height * dpi, res = dpi)
  invisible(p)
}

# --- manifest ---------------------------------------------------------------

#' Write manifest.yml: everything needed to say what this run was.
#'
#' The config is snapshotted in full rather than referenced, because config/
#' changes and the run directory must stay self-describing. The hashes let you
#' tell at a glance whether two runs used identical inputs.
write_manifest <- function(run, status, note = NULL, extra = NULL) {
  stopifnot(inherits(run, "llr_run"))

  m <- list(
    run_id       = run$id,
    prefix       = run$prefix,
    status       = status,
    started_at   = format(run$started, "%Y-%m-%d %H:%M:%S %Z"),
    note         = note %||% run$note %||% "",
    r_version    = as.character(getRversion()),
    platform     = R.version$platform,
    hostname     = unname(Sys.info()[["nodename"]]),
    # `xgboost` produces three of the ten arms and was not in this list;
    # `digest` produces every hash in this manifest (plumbing review F10).
    packages     = .pkg_versions(c("mgcv", "arrow", "qs2", "yaml", "data.table", "targets",
                                   "xgboost", "digest")),
    config_hash  = .hash(run$config),
    vcs          = .git_sha(),
    source_hashes = .source_hashes(),
    # The graph and the runners, hashed beside the library. `source_hashes`
    # covers R/ only, and a Git SHA plus a dirty flag cannot reconstruct an
    # uncommitted edit to `_targets.R` or to the script that assembled this run
    # (plumbing review F10).
    orchestration_hashes = .orchestration_hashes(),
    config       = run$config
  )
  if (!is.null(extra)) m <- utils::modifyList(m, extra)

  yaml::write_yaml(m, run_path(run, "manifest.yml"))
  invisible(m)
}

#' Close a run: stamp elapsed time and flip status to "complete".
#'
#' A run directory without status "complete" was interrupted and its outputs
#' should not be trusted or cited.
finalize_run <- function(run, note = NULL, extra = NULL) {
  stopifnot(inherits(run, "llr_run"))
  elapsed <- as.numeric(difftime(Sys.time(), run$started, units = "mins"))
  write_manifest(run, status = "complete", note = note, extra = c(
    list(finished_at  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
         elapsed_mins = round(elapsed, 2)),
    extra
  ))
  log_msg(run, "run complete in ", round(elapsed, 2), " min")
  invisible(run)
}

#' Stamp a run as FAILED at a named stage, keeping every output already written.
#'
#' EXTERNAL RUNNER REVIEW E1 (2026-09-09). `run/external.R` caught every error
#' from the severity arm, logged "skipped", and finalised the run as
#' `complete`, so a completed manifest no longer established that its enabled
#' arms had succeeded, and the failure reason survived only in the log. The
#' manifest now has a third status beside "running" and "complete": "failed",
#' with the stage and the condition message, written BEFORE the error is
#' re-raised. Primary outputs saved earlier in the run stay on disk; the
#' manifest says which stage did not finish rather than implying they all did.
#'
#' Status semantics, for anyone reading a run directory:
#'   running    interrupted or still in progress; do not cite
#'   failed     a named stage raised; earlier outputs exist, the run is partial
#'   complete   every requested stage finished
fail_run <- function(run, stage, reason, extra = NULL) {
  stopifnot(inherits(run, "llr_run"), is.character(stage), length(stage) == 1L)
  elapsed <- as.numeric(difftime(Sys.time(), run$started, units = "mins"))
  reason <- substr(paste(as.character(reason), collapse = " "), 1L, 2000L)
  write_manifest(run, status = "failed", extra = c(
    list(failed_at      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
         elapsed_mins   = round(elapsed, 2),
         failed_stage   = stage,
         failure_reason = reason),
    extra
  ))
  log_msg(run, "run FAILED at stage `", stage, "`: ", reason)
  invisible(run)
}

#' Evaluate one stage of a runner; on error, stamp the manifest and re-raise.
#'
#' `expr` is a promise and is forced inside the handler, so any condition it
#' raises is caught here. Nothing is swallowed: the error propagates after the
#' manifest has been written, so the console shows the same message and the
#' process exits non-zero. `extra` is the manifest content known at the time
#' of failure (counts, arm statuses), so a failed manifest still describes the
#' stages that did complete.
run_stage <- function(run, stage, expr, extra = NULL) {
  stopifnot(inherits(run, "llr_run"))
  tryCatch(expr, error = function(e) {
    fail_run(run, stage, conditionMessage(e), extra = extra)
    stop(e)
  })
}

# --- export -----------------------------------------------------------------

#' Snapshot targets values into a dated run directory.
#'
#' This is a PLAIN FUNCTION, never a target (hard rule 7). `tar_make()` computes
#' and caches; `export_run()` snapshots. Calling it from inside the graph would
#' put `Sys.time()` into a target's value, changing its hash on every run and
#' invalidating everything downstream.
#'
#' Target names are read with `tar_read()` and routed by argument:
#'   tables      -> tables/<name>.rds + .csv
#'   diagnostics -> diagnostics/<name>.rds + .csv
#'   objects     -> tables/<name>.rds only (non-rectangular values)
#'   bundle      -> <name>.qs2, internal runs only
#'
#' `after` is a callback receiving the open `run` object once every named target
#' has been snapshotted and before the run is finalised. It exists for ONE
#' reason: a figure cannot be a target. Hard rule 9 forbids anything in `R/`
#' from building a path, and hard rule 7 forbids a target's value from carrying
#' a timestamp, so the graph can compute a monotonicity table but has nowhere to
#' put a PNG. Before 2026-09-05 the consequence was that the internal run
#' produced no risk-ordering figures at all while both apply sites produced a
#' full set, and the asymmetry was invisible because nothing errored -- the
#' plots simply did not exist. The callback closes that hole without giving the
#' graph a filesystem: `export_run()` owns the directory, the callback receives
#' it, and `R/` still constructs no path.
#'
#' It runs INSIDE the run's lifetime, so an error in it leaves the manifest at
#' status "running", which is correct: a run whose reporting step failed is
#' incomplete and must not be cited. Nothing is swallowed.
#'
#' @param label short run kind: "internal", "external", "survival"
#' @param cfg   resolved config, snapshotted into the manifest
#' @param after optional `function(run)` called after the exports and before
#'   finalize_run(). Use it for figures and for tables that need a `run`.
#' @examples
#' \dontrun{
#'   targets::tar_make()
#'   export_run("internal", cfg,
#'              tables      = c("eigen", "l_matrix_summary"),
#'              diagnostics = "gam_diagnostics",
#'              bundle      = "bundle")
#' }
export_run <- function(label, cfg,
                       tables      = character(0),
                       diagnostics = character(0),
                       objects     = character(0),
                       bundle      = NULL,
                       root = "out/runs", note = NULL, store = "_targets",
                       after = NULL) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("export_run() needs the `targets` package", call. = FALSE)
  }
  if (!dir.exists(store)) {
    stop("no targets store at '", store, "'. Run targets::tar_make() first.", call. = FALSE)
  }

  run <- new_run(label, cfg, root = root, note = note)
  read1 <- function(nm) targets::tar_read_raw(nm, store = store)

  for (nm in tables)      save_table(run, as.data.frame(read1(nm)), nm)
  for (nm in diagnostics) {
    v <- read1(nm)
    if (is.data.frame(v)) {
      save_table(run, .flatten_list_cols(as.data.frame(v)), nm, subdir = "diagnostics")
    } else {
      # A NAMED LIST OF TABLES gets one CSV per element as well as the .rds.
      # `severity_diag` is the case that forced this: it is a list of five
      # aggregate tables -- APACHE point coverage, SOFA organ coverage, native
      # agreement -- and writing it only as an .rds means the numbers exist but
      # nobody reads them. It was not exported at all until 2026-09-03 (audit
      # finding F4). The .rds is still written, so nothing is lost to the
      # flattening.
      save_object(run, v, nm, subdir = "diagnostics")
      if (is.list(v) && !is.null(names(v)) && length(v) &&
          all(vapply(v, function(z) is.null(z) || is.data.frame(z), logical(1)))) {
        for (el in names(v)) {
          if (is.null(v[[el]]) || !nrow(v[[el]])) next
          save_table(run, .flatten_list_cols(as.data.frame(v[[el]])),
                     paste0(nm, "__", el), subdir = "diagnostics")
        }
      }
    }
  }
  for (nm in objects)     save_object(run, read1(nm), nm)
  if (!is.null(bundle))   save_bundle(run, read1(bundle), "bundle")

  # Provenance. The config snapshot and the R/ source hashes already go into
  # manifest.yml via write_manifest(); tar_meta() adds the per-target record of
  # what was actually computed, and when, which is what makes the snapshot
  # reproducible rather than merely dated.
  meta <- .flatten_list_cols(as.data.frame(targets::tar_meta(store = store)))
  save_table(run, meta, "tar_meta")

  # Figures and anything else that needs the open run directory. See the note
  # on `after` above for why this is a callback rather than a target.
  if (!is.null(after)) {
    if (!is.function(after)) stop("export_run: `after` must be a function", call. = FALSE)
    log_msg(run, "after: reporting step")
    after(run)
  }

  finalize_run(run, extra = list(
    exported = list(tables = as.list(tables), diagnostics = as.list(diagnostics),
                    objects = as.list(objects), bundle = bundle %||% ""),
    targets_store = store,
    n_targets     = nrow(meta)
  ))
  run$path
}

#' tar_meta() carries list columns (`path`, `children`, ...) which write.csv
#' cannot render. Collapse them so the CSV is readable without changing the .rds.
.flatten_list_cols <- function(df) {
  for (j in seq_along(df)) {
    if (is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(z) paste(as.character(z), collapse = ";"),
                        character(1))
    }
  }
  df
}

# --- reading a previous run -------------------------------------------------

#' Read the manifest of a finished run. Used by external/survival to record
#' which internal run produced the bundle they are applying.
read_manifest <- function(run_dir) {
  p <- file.path(run_dir, "manifest.yml")
  if (!file.exists(p)) stop("no manifest.yml in ", run_dir)
  yaml::read_yaml(p)
}

#' Most recent completed run of a given prefix, or NULL. Convenience for
#' interactive work only — runners take an explicit path from their config, so
#' that a result can never silently change because a new run appeared.
latest_run <- function(prefix, root = "out/runs", require_complete = TRUE) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(paste0("^", prefix, "_"), basename(dirs))]
  if (!length(dirs)) return(NULL)
  if (require_complete) {
    ok <- vapply(dirs, function(d) {
      m <- tryCatch(read_manifest(d), error = function(e) NULL)
      !is.null(m) && identical(m$status, "complete")
    }, logical(1))
    dirs <- dirs[ok]
    if (!length(dirs)) return(NULL)
  }
  sort(dirs, decreasing = TRUE)[1L]
}

# --- internal helpers -------------------------------------------------------

# `%||%` was defined here as well until 2026-09-03, byte-identical to the one
# in R/00_utils.R. Removed: `tar_source("R")` sources in filename order, so the
# copy in force was this one, and every file-based attribution of `%||%` named
# the run layer -- which made 27 targets look as though they reached
# `R/11_run.R` and buried the one that genuinely did (audit C4b). Verified
# identical as an AST, as a deparse, and behaviourally before removal.

.check_name <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!grepl("^[a-z0-9_]+$", name)) {
    stop("output name must be lowercase alphanumeric/underscore, got: ", name)
  }
  name
}

.hash <- function(x) substr(digest::digest(x, algo = "md5"), 1L, 12L)

#' Commit SHA, or an explicit statement that there isn't one.
#'
#' This project is not currently a git repository. `system("git rev-parse HEAD",
#' intern = TRUE)` there returns character(0) with a warning, and writing that
#' out produces an EMPTY git_sha.txt — a snapshot that looks provenanced but is
#' not, which is worse than one that admits it has no SHA. So the absence is
#' recorded as a value, and .source_hashes() carries the real provenance until
#' the repo is initialised.
.git_sha <- function() {
  sha <- suppressWarnings(tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)))
  status <- attr(sha, "status")
  if (!length(sha) || (!is.null(status) && status != 0L)) {
    return(list(system = "none", sha = NA_character_,
                note = "not a git repository; see source_hashes"))
  }
  dirty <- suppressWarnings(tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
    error = function(e) NA))
  list(system = "git", sha = sha[1L], dirty = dirty)
}

.pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) "absent")
  })
  stats::setNames(out, pkgs)
}

#' Hash every R source file, for the RUN MANIFEST.
#'
#' The path listing happens here, in the one file allowed to build paths (hard
#' rule 9); the hashing itself is `source_hashes()` in R/00_utils.R. That split
#' was made on 2026-09-03: `build_bundle()` used to call this function, which
#' gave the `bundle` target an undeclared filesystem input and reached into the
#' run layer from the library layer (audit finding F3). `_targets.R` now tracks
#' `R/` as a `format = "file"` target and passes the hashes into
#' `build_bundle()` as an argument, so the bundle's provenance is a declared
#' dependency rather than a hidden read.
#'
#' This wrapper stays because the run manifest legitimately wants the same
#' thing at export time, and export is not a target.
.source_hashes <- function(dir = "R") {
  source_hashes(list.files(dir, pattern = "\\.R$", full.names = TRUE))
}

#' Hash the orchestration: `_targets.R`, every runner under run/, and every
#' declaration under config/. Same hasher, same one-file-builds-paths rule.
#' Files that do not exist are simply not listed, so a run directory made by a
#' script outside run/ still gets a manifest.
#'
#' THE CONFIG FILES JOINED ON 2026-09-09 (external runner review E3). The
#' manifest's `config` block snapshots whatever the runner passed to
#' `new_run()`, and for the external runner that used to be the LOCAL MIMIC
#' config rather than `config/external.yml`, so two external runs on different
#' eICU inputs or hospital thresholds could carry the same `config_hash` and no
#' hash anywhere covered the file that actually steered the run. Hashing the
#' whole config/ directory here closes that independently of what a runner
#' chooses to snapshot. Names are basenames, as for run/ and R/; the two
#' directories share no filename.
.orchestration_hashes <- function() {
  paths <- c(if (file.exists("_targets.R")) "_targets.R",
             list.files("run", pattern = "\\.R$", full.names = TRUE),
             list.files("config", pattern = "\\.(yml|yaml|csv)$", full.names = TRUE))
  if (!length(paths)) return(list())
  source_hashes(paths)
}

#' Fingerprint a set of input files by size and full MD5, for a run manifest.
#'
#' CONTENT IDENTITY WITHOUT CONTENT (external runner review E3, hard rule 1).
#' A path in a manifest is not an identity: the file behind it can be
#' re-extracted or replaced, and a bundle path likewise. The MD5 of the bytes
#' is, and it reveals nothing about any row. Files that do not exist are
#' reported as such rather than dropped, so a manifest cannot quietly describe
#' fewer inputs than the config declared.
#'
#' @param paths named character vector or list of file paths
#' @return a list, one element per name: path, exists, bytes, md5
.input_fingerprint <- function(paths) {
  paths <- unlist(paths, use.names = TRUE)
  if (!length(paths)) return(list())
  nm <- names(paths) %||% basename(paths)
  nm[!nzchar(nm)] <- basename(paths)[!nzchar(nm)]
  out <- lapply(seq_along(paths), function(i) {
    p <- as.character(paths[[i]])
    ex <- file.exists(p)
    list(path = p, exists = ex,
         bytes = if (ex) as.numeric(file.size(p)) else NA_real_,
         md5   = if (ex) unname(as.character(tools::md5sum(p))) else NA_character_)
  })
  stats::setNames(out, nm)
}
