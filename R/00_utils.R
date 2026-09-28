# R/00_utils.R ---------------------------------------------------------------
# Small shared helpers. No paths, no data access, and no clock: `start_timer()`
# below measures a DURATION off `proc.time()` and never a datetime (hard rule 9).
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

#' Evaluate `expr` under a fixed RNG seed, restoring the caller's stream after.
#'
#' Fold assignment must be reproducible without the pipeline silently reseeding
#' the global stream for everything downstream.
with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

#' Percentage, formatted for logs.
pct <- function(x, digits = 2) paste0(round(100 * x, digits), "%")

#' Stop with a message listing offending values, truncated so a bad load cannot
#' dump 50k rows into the console (hard rule 1).
#'
#' USE ONLY FOR CONTROLLED VOCABULARIES — column names, signal and intervention
#' names, factor levels, config keys. Never for anything row-level: not
#' `stay_id`, not `subject_id`, not a measurement, not a free-text column. Those
#' are under a PhysioNet DUA and must never reach a console, a log, or a model
#' context. For row-level problems report a COUNT and let the user inspect
#' locally — see check_folds() in R/03_folds.R for the pattern.
abort_values <- function(msg, values, max_show = 10L) {
  v <- unique(as.character(values))
  shown <- paste(utils::head(v, max_show), collapse = ", ")
  more  <- if (length(v) > max_show) paste0(" ... and ", length(v) - max_show, " more") else ""
  stop(msg, ": ", shown, more, call. = FALSE)
}

#' Are two character vectors the same set?
same_set <- function(a, b) setequal(as.character(a), as.character(b))

#' Hash a set of source files, so a run can be traced to the code that made it.
#'
#' TAKES PATHS, BUILDS NONE. That is the whole reason it lives here rather than
#' in `R/11_run.R`: `build_bundle()` used to call `R/11_run.R`'s
#' `.source_hashes()`, which does `list.files("R")` itself, and that made two
#' things true at once (audit finding F3). It breached hard rule 9, because a
#' library-layer function was reaching into the one file allowed to build
#' paths. And more seriously it gave the `bundle` target an UNDECLARED
#' FILESYSTEM INPUT: the contents of `R/` decided part of the target's value
#' and `targets` could not see it, so the value was not a function of its
#' declared dependencies.
#'
#' With the paths passed in, `_targets.R` can track them as a `format = "file"`
#' target and the input becomes declared. The provenance still travels inside
#' the bundle, which is where it belongs.
#'
#' Stands in for a commit hash; this project is not a git repository.
source_hashes <- function(paths) {
  fs <- sort(as.character(paths))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(list())
  stats::setNames(
    lapply(fs, function(f) unname(substr(tools::md5sum(f), 1L, 12L))),
    basename(fs)
  )
}

#' A config value that MUST be declared, with no fallback.
#'
#' WHY THIS EXISTS RATHER THAN `cfg$a$b %||% <literal>`. The `%||%` idiom is
#' good for optional things and wrong for settings that decide what model gets
#' fitted, because it writes the value down TWICE -- once in `config.yml` and
#' once as the fallback literal -- and nothing makes the two agree. They can be
#' edited apart, and the divergence is then invisible until the config key is
#' deleted or misspelled, at which point the pipeline quietly fits a different
#' model instead of failing. Audit E found four such disagreements on
#' 2026-09-03, of which `bam.nthreads` was measured to change fitted L values
#' at the 1e-13 level (audit findings F9, F10).
#'
#' `cfg_req()` removes the second declaration entirely. There is nothing to
#' disagree with, and an absent key is a loud stop naming the path.
#'
#' @param cfg  a config list
#' @param ...  the path, one key per argument: `cfg_req(cfg, "bam", "method")`
#' @param what optional context for the error message
cfg_req <- function(cfg, ..., what = NULL) {
  path <- as.character(c(...))
  v <- cfg
  for (i in seq_along(path)) {
    if (!is.list(v) || is.null(v[[path[i]]])) {
      stop("config is missing required key `", paste(path, collapse = "."), "`",
           if (i > 1L) paste0(" (stopped at `", paste(path[seq_len(i)], collapse = "."), "`)") else "",
           if (!is.null(what)) paste0(": ", what) else "",
           ".\nThis setting has no default on purpose -- a fallback here would ",
           "be a second declaration of the same thing, and the two can be ",
           "edited apart without anything noticing.", call. = FALSE)
    }
    v <- v[[path[i]]]
  }
  v
}

#' A required config flag, coerced to a single TRUE/FALSE.
#'
#' Same argument as `cfg_req()`, and the reason it is separate is that
#' `isTRUE(cfg$x %||% FALSE)` is the worst version of the pattern: an absent or
#' misspelled key silently turns a construct OFF, and a construct that is off
#' while the methods section says it is on is precisely the failure this
#' project has already had once with `delta`.
cfg_flag <- function(cfg, ..., what = NULL) {
  v <- cfg_req(cfg, ..., what = what)
  if (!is.logical(v) || length(v) != 1L || is.na(v)) {
    stop("config key `", paste(as.character(c(...)), collapse = "."),
         "` must be a single TRUE or FALSE, got: ",
         paste(utils::head(as.character(v), 5), collapse = ", "), call. = FALSE)
  }
  v
}

#' Time a long-running arm without reading the clock.
#'
#' WHY THIS IS NOT A HARD RULE 9 BREACH, and the distinction is the whole
#' point. `Sys.time()` reads the CALENDAR: it returns a datetime, which is
#' exactly the thing rule 9 keeps out of `R/` so that a run's identity enters
#' through the run object and through nothing else. `proc.time()` reads two
#' counters this R process has kept since it started -- CPU consumed and
#' seconds elapsed -- and cannot produce a datetime at all. The difference of
#' two of them is a DURATION. Nothing here can stamp when a run happened, so
#' MIMIC-test and eICU still provably execute the same code path, which is what
#' rule 9 is for.
#'
#' HARD RULE 7 IS A SEPARATE QUESTION AND STILL BINDS. A duration does differ
#' between two runs of identical code, so it must never enter a target's value.
#' `tests/audit_c_targets.R` keeps `proc.time` in its CLOCK list precisely so
#' that a target whose closure reaches this is reported. Print a duration, log
#' it, write it into a run directory -- never cache it.
#'
#' `user.child` and `sys.child` are NA on Windows, so CPU is counted for this
#' process only. Nothing is lost by that here: `bam(nthreads = )` and `xgboost`
#' both use threads inside this process rather than child processes.
#'
#' @return a function of no arguments, returning elapsed and CPU seconds since
#'   the timer was created.
start_timer <- function() {
  t0 <- proc.time()
  function() {
    d <- proc.time() - t0
    list(elapsed_sec = unname(d[["elapsed"]]),
         cpu_sec     = unname(d[["user.self"]] + d[["sys.self"]]))
  }
}
