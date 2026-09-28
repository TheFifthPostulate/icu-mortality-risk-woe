# R/03_folds.R ---------------------------------------------------------------
# Split and fold assignment. Defined once; everything downstream keys to it.
#
# Two nested levels:
#   1. split  - train / test, 80/20. Test is touched ONCE, for the APACHE II /
#               SOFA
#               discrimination comparison.
#   2. folds  - 5 CV folds INSIDE train, for the out-of-fold L's.
#
# Both are grouped by patient and stratified on mortality. Grouping is applied
# FIRST: a patient's stratum label is the max of the stratum column over their
# stays, and patients are assigned whole, so no repeat stay can straddle a
# boundary.
#
# EACH LEVEL DECLARES ITS OWN GROUPING AND STRATIFICATION, and as of 2026-09-03
# it actually reads them. Four config keys govern this file:
#
#   split.group_by        grouping for the train/test split
#   split.stratify_on     stratum column for the train/test split
#   fold_group_by         grouping for the 5 CV folds inside train
#   fold_stratify_on      stratum column for those folds
#
# They are separate keys because the split and the fold assignment are two
# separate resampling operations, and each is entitled to its own answer. Until
# 2026-09-03 three of the four were declared and read by nothing: `mortality`
# was hard-coded as the stratum at both levels, and `split.group_by` was used
# for both levels' grouping. The declarations agreed with what the code did, so
# no number was ever wrong -- but a control that reads as being in force and is
# not is the failure class this pipeline is built against (audit finding F5,
# docs/v2_audit_findings_20260903.md).
#
# The fold-level keys DEFAULT TO THE SPLIT-LEVEL ONES rather than to a literal,
# so the two levels agree unless somebody deliberately separates them, and
# there is no second place where "mortality" is written down.
#
# Measured 2026-08-22: in MIMIC-IV every patient appears exactly once
# (n_subject == n_stay == n_hadm == 51,563), so grouping is a no-op here and
# stratification is exact. The machinery still runs; eICU may differ.
#
# No paths, no clock (hard rule 9). Seed comes from config.
# ----------------------------------------------------------------------------

#' Resolve a declared resampling column against the cohort.
#'
#' Applies the ONE documented alias and nothing else: spec §5.1 was revised to
#' table-defined names, so the patient key is `subject_id`, and a config still
#' saying `patient_id` is accepted. Anything else that does not resolve is an
#' error rather than a silent substitution -- the whole point of wiring these
#' keys is that a typo must be loud.
.resolve_resample_col <- function(name, cohort, key) {
  nm <- as.character(name)[1]
  if (!nzchar(nm) || is.na(nm)) {
    stop("config `", key, "` is empty; name a cohort column", call. = FALSE)
  }
  if (nm %in% names(cohort)) return(nm)
  if (identical(nm, "patient_id") && "subject_id" %in% names(cohort)) return("subject_id")
  abort_values(paste0("config `", key, "` names a column the cohort does not have"), nm)
}

#' A stratum label as an integer, with the reasons it might not be one stated.
#'
#' `.allocate()` shuffles within `sort(unique(y))` and balances stays inside
#' each level, so the stratum must be a small set of ordered labels with no
#' NA. A character stratum would work mechanically and is refused anyway: the
#' patient-level label is `max()` over the patient's stays, and the "worst" of
#' a set of character labels is a lexicographic accident rather than a
#' clinical statement.
.strat_label <- function(v, key) {
  if (anyNA(v)) stop("config `", key, "` column has NA; a stratum must be complete",
                     call. = FALSE)
  if (is.logical(v)) return(as.integer(v))
  if (!is.numeric(v)) {
    stop("config `", key, "` names a non-numeric column. A stratum label is ",
         "reduced to the patient level with max(), which is only meaningful ",
         "for an ordered numeric or logical column.", call. = FALSE)
  }
  as.integer(v)
}

#' The group-level table `.allocate()` consumes: one row per group, in
#' `aggregate()`'s sorted order.
#'
#' ROW ORDER IS LOAD-BEARING AND IS WHY THIS USES `aggregate()` RATHER THAN
#' ANYTHING TIDIER. `.allocate()` shuffles `which(y == lab)` under a fixed
#' seed, so the assignment a group receives depends on where it sits in this
#' table. Changing how the table is built reshuffles every fold and moves every
#' number downstream, whether or not anything else changed.
.group_table <- function(grp, strat) {
  t <- stats::aggregate(list(y = strat), by = list(grp = grp), FUN = max)
  t$n_stays <- as.integer(table(grp)[t$grp])
  t
}

#' Assign train/test and CV folds for a cohort.
#'
#' @param cohort data frame with stay_id, subject_id, mortality, plus whatever
#'   columns the four resampling keys name
#' @param cfg    from load_config()
#' @return data frame: stay_id, subject_id, mortality, split, fold
#'         `fold` is NA for test rows -- they belong to no fold by construction.
#'         `mortality` is carried unconditionally rather than being the
#'         stratum column, because downstream needs the OUTCOME (`y_train`,
#'         the balance report) whatever the split was stratified on.
assign_folds <- function(cohort, cfg) {
  need <- c("stay_id", "subject_id", "mortality")
  miss <- setdiff(need, names(cohort))
  if (length(miss)) abort_values("cohort missing columns needed for folds", miss)
  if (anyNA(cohort$mortality)) stop("cohort: `mortality` has NA", call. = FALSE)
  if (any(duplicated(cohort$stay_id))) stop("cohort: `stay_id` is not unique", call. = FALSE)

  d <- data.frame(
    stay_id    = cohort$stay_id,
    subject_id = cohort$subject_id,
    mortality  = as.integer(cohort$mortality),
    stringsAsFactors = FALSE
  )

  ck <- resample_cols(cohort, cfg)

  test_frac <- cfg$split$test_fraction %||% 0.2
  n_folds   <- cfg$n_folds %||% 5L

  # --- level 1: train / test ------------------------------------------------
  # Patient-level stratum label: max over the patient's stays. A patient who
  # died in any admission counts as an event for allocation purposes.
  g1 <- as.character(cohort[[ck$split_group]])
  y1 <- .strat_label(cohort[[ck$split_stratify]], "split.stratify_on")
  p1 <- .group_table(g1, y1)
  p1$split <- with_seed(cfg$seed, .allocate(p1, k = NULL, frac = test_frac))
  d$split <- p1$split[match(g1, p1$grp)]

  # --- level 2: folds inside train -----------------------------------------
  # A separate seed stream, derived from the same root, so that changing
  # test_fraction cannot silently reshuffle the folds and vice versa.
  #
  # The fold table is built from the TRAIN ROWS rather than by filtering the
  # split table, so that `fold_group_by` and `fold_stratify_on` can genuinely
  # differ from the split-level ones. When they do not differ -- the only
  # configuration this project has ever run -- the two constructions give the
  # identical table: `aggregate()` sorts by group either way, a patient wholly
  # inside train contributes the same stratum label and the same stay count,
  # and the row order that `.allocate()` depends on is therefore preserved.
  # That equivalence was verified against the pre-wiring code before this
  # change was made, fold for fold over all 51,563 stays.
  is_tr <- d$split == "train"
  g2 <- as.character(cohort[[ck$fold_group]])
  y2 <- .strat_label(cohort[[ck$fold_stratify]], "fold_stratify_on")
  p2 <- .group_table(g2[is_tr], y2[is_tr])
  p2$fold <- with_seed(cfg$seed + 1L, .allocate(p2, k = n_folds, frac = NULL))

  # Assigned to TRAIN ROWS ONLY. Under a fold grouping that differs from the
  # split grouping a test stay can share a fold group with a train stay, and a
  # test row must carry no fold whatever the grouping says.
  d$fold <- NA_integer_
  d$fold[is_tr] <- as.integer(p2$fold[match(g2[is_tr], p2$grp)])

  d[order(d$stay_id), ]
}

#' The four resampling columns, resolved against a cohort.
#'
#' One function so that `assign_folds()` and `check_folds()` cannot disagree
#' about which column was used for what -- which is the same argument as
#' `priors_for()` returning one container rather than four arguments.
#'
#' The fold-level keys default to the split-level ones, not to a literal. There
#' is exactly one place in this project where "mortality" and "patient_id"
#' appear as resampling defaults, and it is here.
resample_cols <- function(cohort, cfg) {
  sg <- cfg$split$group_by    %||% "patient_id"
  sy <- cfg$split$stratify_on %||% "mortality"
  list(
    split_group    = .resolve_resample_col(sg, cohort, "split.group_by"),
    split_stratify = .resolve_resample_col(sy, cohort, "split.stratify_on"),
    fold_group     = .resolve_resample_col(cfg$fold_group_by    %||% sg, cohort,
                                           "fold_group_by"),
    fold_stratify  = .resolve_resample_col(cfg$fold_stratify_on %||% sy, cohort,
                                           "fold_stratify_on"))
}

#' The patient of every stay, as the dependence unit every resampling reads.
#'
#' ONE RESOLVER (2026-09-09). The cluster bootstrap (R/09), the boosters'
#' inner early-stopping split (R/09b), the attribution refits (R/13b) and the
#' replicate scripts all need "the patient of each stay", and each had begun
#' to write `tabs$cohort[[resample_cols(...)$...]][match(...)]` for itself --
#' six copies, two of which read `split_group` and four `fold_group`. The
#' unit is `fold_group_by`: the grouping the CV folds and the attribution bags
#' already declare, and which defaults to `split.group_by` unless somebody
#' deliberately separates the two levels. Returned as character and aligned
#' to `stay_ids`; a stay with no grouping value is an error, never NA.
#'
#' AGGREGATES ONLY on print: this vector is row-level and is never printed.
patient_group_of <- function(cohort, cfg, stay_ids) {
  col <- resample_cols(cohort, cfg)$fold_group
  v <- cohort[[col]][match(as.character(stay_ids), as.character(cohort$stay_id))]
  if (anyNA(v)) {
    stop(sprintf("patient_group_of: %d stay(s) have no `%s` value", sum(is.na(v)), col),
         call. = FALSE)
  }
  as.character(v)
}

#' Hold the four resampling declarations to the cohort, at BOTH sites.
#'
#' Same idiom as `check_signal_tails()` and the rest of the design-check family:
#' a static declaration re-checked against the data, so that a column present at
#' MIMIC and absent or degenerate at eICU is a loud failure rather than a
#' silently different partition. Registered in `design_checks`.
#'
#' What it reports per key: the resolved column, its distinct-value count, and
#' whether it is usable. A stratum with ONE level is legal but is not
#' stratification, and a grouping with as many groups as rows is legal but is
#' not grouping; both are reported so the difference is visible rather than
#' inferred.
#'
#' AGGREGATES ONLY (hard rule 1): counts, never identifiers.
check_resampling_cols <- function(tabs, cfg, strict = TRUE) {
  cohort <- tabs$cohort
  ck <- resample_cols(cohort, cfg)
  rows <- lapply(names(ck), function(k) {
    col <- ck[[k]]
    v <- cohort[[col]]
    is_strat <- grepl("stratify$", k)
    data.frame(
      key = switch(k, split_group = "split.group_by",
                      split_stratify = "split.stratify_on",
                      fold_group = "fold_group_by",
                      fold_stratify = "fold_stratify_on"),
      role = if (is_strat) "stratum" else "grouping",
      column = col,
      n_distinct = length(unique(v)),
      n_rows = length(v),
      n_na = sum(is.na(v)),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$ok <- out$n_na == 0L &
            ifelse(out$role == "stratum", out$n_distinct >= 2L, out$n_distinct >= 1L)
  out$note <- ifelse(out$n_na > 0L, "has NA",
              ifelse(out$role == "stratum" & out$n_distinct < 2L,
                     "single-valued: not stratification",
              ifelse(out$role == "grouping" & out$n_distinct == out$n_rows,
                     "one group per row: grouping is a no-op here", "")))
  if (any(!out$ok)) {
    bad <- paste0(out$key[!out$ok], " (", out$column[!out$ok], ": ",
                  out$note[!out$ok], ")")
    msg <- paste0("resampling declarations unusable at this site: ",
                  paste(bad, collapse = "; "),
                  ". Fix config/config.yml or the extraction, not R/03_folds.R.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  out
}

#' Allocate patients to groups, balancing stays within each stratum.
#'
#' Shuffle patients inside a stratum, then walk them in order assigning each to
#' whichever bucket currently holds the fewest STAYS. Balancing on stays rather
#' than patients is what keeps the split at the requested fraction when patients
#' contribute unequal numbers of stays.
#'
#' @param k    number of equal buckets (folds), or NULL for a fractional split
#' @param frac size of the "test" bucket, or NULL when k is given
.allocate <- function(pat, k = NULL, frac = NULL) {
  stopifnot(xor(is.null(k), is.null(frac)))
  out <- character(nrow(pat))

  for (lab in sort(unique(pat$y))) {
    idx <- which(pat$y == lab)
    idx <- idx[sample.int(length(idx))]          # shuffle within stratum
    sizes <- pat$n_stays[idx]

    if (!is.null(frac)) {
      target <- frac * sum(sizes)
      cum <- cumsum(sizes)
      # Take patients until the next one would overshoot the target by more
      # than it undershoots -- i.e. round to the nearest achievable boundary.
      n_take <- sum(cum <= target)
      if (n_take < length(cum) &&
          abs(cum[n_take + 1L] - target) < abs(target - c(0, cum)[n_take + 1L])) {
        n_take <- n_take + 1L
      }
      out[idx] <- "train"
      if (n_take > 0L) out[idx[seq_len(n_take)]] <- "test"
    } else {
      load <- numeric(k)
      assign_to <- integer(length(idx))
      for (i in seq_along(idx)) {
        j <- which.min(load)
        assign_to[i] <- j
        load[j] <- load[j] + sizes[i]
      }
      out[idx] <- as.character(assign_to)
    }
  }
  out
}

#' One-row-per-bucket summary. Cheap, and it is how you confirm at a glance
#' that stratification and grouping actually held.
fold_summary <- function(folds) {
  f <- folds
  f$bucket <- ifelse(f$split == "test", "test", paste0("train/fold", f$fold))
  s <- do.call(rbind, lapply(split(f, f$bucket), function(z) data.frame(
    bucket     = z$bucket[1],
    n_stays    = nrow(z),
    n_patients = length(unique(z$subject_id)),
    deaths     = sum(z$mortality),
    mortality  = round(mean(z$mortality), 4),
    stringsAsFactors = FALSE)))
  s <- s[order(s$bucket), ]
  rownames(s) <- NULL
  s
}

#' Assertions the split must satisfy. Called by the runner after assignment.
check_folds <- function(folds, cfg, tol = 0.02, cohort = NULL) {
  n_folds <- cfg$n_folds %||% 5L

  if (anyNA(folds$split)) stop("folds: some stays have no split", call. = FALSE)
  te <- folds$split == "test"
  if (any(!is.na(folds$fold[te]))) stop("folds: a test stay was given a fold", call. = FALSE)
  if (anyNA(folds$fold[!te])) stop("folds: a train stay has no fold", call. = FALSE)

  got <- sort(unique(folds$fold[!te]))
  if (!identical(as.integer(got), seq_len(n_folds))) {
    abort_values("folds: fold ids are not 1..n_folds", got)
  }

  # WHAT IS CHECKED DEPENDS ON WHAT WAS DECLARED, which is why `cohort` is
  # worth passing. Without it this function can only check the columns the
  # folds table happens to carry -- `subject_id` and `mortality` -- and if the
  # config declared different resampling columns it would be asserting
  # something other than what was done. With it, the check reads the same four
  # keys `assign_folds()` read, through the same `resample_cols()`.
  ck <- if (!is.null(cohort)) resample_cols(cohort, cfg) else NULL
  key <- function(col_role, fallback) {
    if (is.null(ck)) return(folds[[fallback]])
    cohort[[ck[[col_role]]]][match(folds$stay_id, cohort$stay_id)]
  }

  # No group may straddle a boundary -- the whole point of grouping. Both
  # levels are checked, because they can be grouped differently: a split group
  # must not straddle train/test, and a fold group must not straddle folds.
  # Reported as COUNTS, never as identifiers (hard rule 1).
  g_split <- as.character(key("split_group", "subject_id"))
  n_bad <- sum(tapply(folds$split, g_split, function(z) length(unique(z))) > 1L)
  if (n_bad) {
    stop(sprintf(paste0("folds: %d split-group(s) straddle train/test. ",
                        "Inspect locally; ids are not printed."), n_bad), call. = FALSE)
  }
  g_fold <- as.character(key("fold_group", "subject_id"))[!te]
  n_bad <- sum(tapply(folds$fold[!te], g_fold, function(z) length(unique(z))) > 1L)
  if (n_bad) {
    stop(sprintf(paste0("folds: %d fold-group(s) straddle folds. ",
                        "Inspect locally; ids are not printed."), n_bad), call. = FALSE)
  }

  # Stratification held, within tolerance -- on the column that was actually
  # stratified on. The two levels can name different columns, so each is
  # checked against its own.
  y_split <- .strat_label(key("split_stratify", "mortality"), "split.stratify_on")
  if (abs(mean(y_split[te]) - mean(y_split)) > tol) {
    stop(sprintf("folds: split stratum rate %.4f in test against %.4f overall, off by more than %.2f",
                 mean(y_split[te]), mean(y_split), tol), call. = FALSE)
  }
  y_fold <- .strat_label(key("fold_stratify", "mortality"), "fold_stratify_on")[!te]
  by_fold <- tapply(y_fold, folds$fold[!te], mean)
  off <- names(by_fold)[abs(by_fold - mean(y_fold)) > tol]
  if (length(off)) {
    abort_values(sprintf("folds: fold stratum rate deviates >%.0f%% from %.4f in fold(s)",
                         100 * tol, mean(y_fold)), off)
  }

  frac <- mean(te)
  target <- cfg$split$test_fraction %||% 0.2
  if (abs(frac - target) > tol) {
    stop(sprintf("folds: test fraction %.4f differs from target %.4f by more than %.2f",
                 frac, target, tol), call. = FALSE)
  }
  invisible(TRUE)
}
