# R/07_layer2.R --------------------------------------------------------------
# The L matrix, its correlation structure, and the eigenspectrum.
#
# THIS IS THE BRANCH POINT, not a step. The question it answers is whether the
# 19 L's carry more than one direction of information. If one component
# dominates, layer 2 is weakly motivated and the paper's claims change — so
# nothing downstream of this file should be built before its output is read.
#
# Three matrices, all pivots of layer-1 OUTPUT (hard rule 2 — this is the only
# genuinely wide object in the project):
#
#   L_full   the primary joint measurement-intervention LLR
#   L_intv   the pure intervention-propensity LLR
#   L_cond   L_full - L_intv, the measurement-deviation LLR conditional on
#            class and intervention context. Derived, never fitted.
#
# and each is built twice, because the answer can differ:
#
#   fill = "zero"  unmeasured stays get L = 0, the value the design assigns
#                  them (spec §5.5). This is the matrix layer 2 would use.
#   fill = "na"    unmeasured stays are missing, and cor() runs pairwise.
#                  Co-ordered labs correlate through co-missingness ALONE. If
#                  that rather than physiology is driving PC1, the branch-point
#                  conclusion changes — which is why this is not optional.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). Nothing here reads the filesystem or
# Sys.time(); the run object carries both.
#
# AGGREGATES ONLY on print (hard rule 1). The L matrix itself is row-level and
# must never be printed, summarised per stay, or written anywhere but a run
# directory. Every function below that prints emits eigenvalues, counts or
# correlations — never a row.
# ----------------------------------------------------------------------------

# --- assembly ---------------------------------------------------------------

#' Which stays have a measured value for each signal.
#'
#' The distinction the zero/NA fill turns on. `n_obs > 0` is the same predicate
#' `signal_frame()` applies, read from the same column, so the two cannot drift.
#'
#' @return logical matrix, stays x signals
measured_matrix <- function(tabs, cfg, stay_ids) {
  sf <- tabs$signal_features
  M <- matrix(FALSE, nrow = length(stay_ids), ncol = length(cfg$signals),
              dimnames = list(as.character(stay_ids), unlist(cfg$signals)))
  for (sg in cfg$signals) {
    z <- sf[sf$signal == sg, c("stay_id", "n_obs"), drop = FALSE]
    i <- match(as.character(stay_ids), as.character(z$stay_id))
    M[, sg] <- !is.na(i) & z$n_obs[i] > 0
  }
  M
}

#' Pivot the long out-of-fold L's into one stays x signals matrix.
#'
#' Handles the two kinds of row layer 1 never fits, and READS THE SAME
#' `spec_source()` `layer1_jobs()` reads rather than restating its rules:
#'
#'   alias   the 7 unpaired signals have no `full` fit, because `full` and
#'           `meas` are the same formula there; `full` reads `meas`, and so do
#'           both interaction models. `full_ti_trend` additionally reads `full`
#'           on the three paired signals with no `trend` covariate.
#'   assign  the 7 unpaired signals have no `intv` fit either: the formula
#'           would be `mortality ~ 1`, so L_intv = 0 exactly on every MEASURED
#'           stay. Assigned, never fitted.
#'
#' An alias can now CHAIN -- an unpaired signal's `full_ti_all` points at
#' `meas`, which is fitted, but nothing in `spec_source()` guarantees a
#' one-step hop -- so the source is resolved in a loop with a depth bound
#' rather than by a single lookup. A cycle would otherwise hang the pivot after
#' 215 fits.
#'
#' The assignment rule is why `fill` cannot be applied blindly. Under
#' fill = "na" an unpaired signal's `intv` column is 0 where measured and NA
#' where not — a constant column, which is correct and which cor() will report
#' as NA. That is a fact about the design, not a bug: those 7 signals have no
#' intervention, so they contribute no intervention evidence.
#'
#' @param l_long  the `l` element of run_layer1(), long
#' @param model   "full", "meas" or "intv"
#' @param fill    "zero" or "na", for stays with no measured value
#' @return numeric matrix, stays x signals
l_matrix <- function(l_long, model = LAYER1_MODELS, tabs, cfg, stay_ids,
                     fill = c("zero", "na")) {
  model <- match.arg(model)
  fill  <- match.arg(fill)

  meas_ok <- measured_matrix(tabs, cfg, stay_ids)
  ids <- as.character(stay_ids)
  M <- matrix(if (fill == "zero") 0 else NA_real_,
              nrow = length(ids), ncol = length(cfg$signals),
              dimnames = list(ids, unlist(cfg$signals)))

  for (sg in cfg$signals) {
    # The alias walk lives in R/05 (`resolve_spec_source()`), shared with the
    # nested cross-fit, so the pivot and the fit enumeration read one rule.
    src <- resolve_spec_source(sg, model, cfg)

    if (identical(src, "zero")) {
      # assigned: 0 on measured stays, fill elsewhere
      M[meas_ok[, sg], sg] <- 0
      next
    }

    z <- l_long[l_long$signal == sg & l_long$model == src, , drop = FALSE]
    if (!nrow(z)) {
      stop(sprintf("l_matrix: no L rows for signal '%s' model '%s'", sg, src),
           call. = FALSE)
    }
    if (any(duplicated(z$stay_id))) {
      # Out-of-fold means each stay is predicted exactly once. More than one row
      # is a fold-scoping bug, and it would average away silently.
      stop(sprintf("l_matrix: %d duplicate stay(s) for '%s'/'%s' — a stay was ",
                   "predicted in more than one fold. Ids not printed (hard rule 1).",
                   sum(duplicated(z$stay_id)), sg, src), call. = FALSE)
    }
    i <- match(ids, as.character(z$stay_id))
    M[!is.na(i), sg] <- z$l[i[!is.na(i)]]

    # A stay that is measured but has no L is a scoping failure, not a fill case.
    gap <- meas_ok[, sg] & is.na(i)
    if (any(gap)) {
      stop(sprintf("l_matrix: %d measured stay(s) have no L for '%s'/'%s'",
                   sum(gap), sg, src), call. = FALSE)
    }
  }
  M
}

#' Every matrix at one fill, plus the derived conditionals.
#'
#' `L_cond = L_full - L_intv` is computed here rather than fitted. On the 7
#' unpaired signals L_intv is 0, so L_cond equals L_full there — expected, and
#' the same statement as `meas == full` for those signals.
#'
#' THE SAME SUBTRACTION IS APPLIED TO BOTH INTERACTION MODELS, and it is the
#' same subtraction rather than an analogous one. `full_ti_*` adds cross terms
#' to `full`'s two blocks and adds nothing to `intv`'s, so the term sets still
#' partition: `L_cond_ti_* = L_full_ti_* - L_intv` is the measurement-deviation
#' LLR conditional on intervention context, WITH treatment response expressible.
#' `p_bar` cancels for exactly the reason it cancels in the additive case --
#' both models are fitted on the same rows with the same prior.
#'
#' Six matrices where there were four, and `attr_derive_arms()` in R/14 builds
#' the identical eight arms from five bases. The two are the same arithmetic in
#' two places, which is tolerable only because the replicate store predates the
#' spec change; `tests/attr_replicates.R` now reads these matrices directly.
#'
#' `meas` IS REQUIRED BY DEFAULT, as of 2026-09-08 (plumbing review F2). Until
#' then it was wrapped in `try(..., silent = TRUE)` and dropped on ANY error,
#' on the argument that it was information-only. That stopped being true on
#' 2026-09-01 when `llr_meas` became a scored arm and the APACHE II
#' counterpart, and the `try()` was hiding the wrong thing anyway: an ABSENT
#' `meas` (a branch-point script that never fitted it) and a CORRUPT `meas`
#' (a measured stay with no L, or a stay predicted in two folds -- exactly what
#' `l_matrix()` stops for) both came out as "no meas arm". Reproduced with
#' synthetic predictions: delete one `mbp/meas` row and the primary graph
#' produced seven matrices and silently no eighth.
#'
#' THE MECHANISM: a broad error catch on a required input converts a
#' correctness failure into an absent output, and `oof_scores` then dropped the
#' NULL arm, so the run succeeded with nine arms and nothing said so.
#'
#' `require_meas = FALSE` is the explicit legacy mode for a caller that
#' deliberately fitted no `meas` (tests/branch_point.R can). Even then an error
#' from `l_matrix()` PROPAGATES: absence is a mode, corruption never is.
l_matrices <- function(l_long, tabs, cfg, stay_ids, fill = c("zero", "na"),
                       require_meas = TRUE) {
  fill <- match.arg(fill)
  get1 <- function(md) l_matrix(l_long, md, tabs, cfg, stay_ids, fill = fill)

  out <- list(full = get1("full"), intv = get1("intv"))
  out$cond <- out$full - out$intv
  for (md in LAYER1_TI_MODELS) {
    out[[md]] <- get1(md)
    out[[sub("^full", "cond", md)]] <- out[[md]] - out$intv
  }
  paired <- cfg$signals[vapply(cfg$signals, function(sg)
    length(interventions_of(sg, cfg)) > 0L, logical(1))]
  has_meas <- any(l_long$model == "meas" & l_long$signal %in% paired)
  if (!has_meas && isTRUE(require_meas)) {
    stop("l_matrices: no `meas` predictions for any paired signal. `llr_meas` ",
         "is a scored arm and the APACHE II counterpart, so the primary graph ",
         "requires it; pass require_meas = FALSE only from a script that ",
         "deliberately fitted no `meas` model.", call. = FALSE)
  }
  if (has_meas) out$meas <- get1("meas")
  attr(out, "fill") <- fill
  out
}

# --- correlation and spectrum ------------------------------------------------

#' Sigma: the correlation matrix of the L's, estimated on out-of-fold values.
#'
#' Out-of-fold L's are systematically noisier than final-fit L's, so the two are
#' not identically distributed. Estimating here and applying to L_test is the
#' correct choice and every alternative is worse; it is a limitation to state in
#' the methods, not a bug to fix (CLAUDE.md).
#'
#' Constant columns are dropped rather than left as NA, and reported: a constant
#' column is an unpaired signal's `intv`, which carries no information by
#' construction. Keeping it would make every downstream eigen() fail.
#'
#' @return correlation matrix, with attr "dropped" naming any constant columns
l_correlation <- function(M) {
  keep <- vapply(seq_len(ncol(M)), function(j) {
    v <- M[, j]; v <- v[!is.na(v)]
    length(v) > 1L && stats::sd(v) > 0
  }, logical(1))
  dropped <- colnames(M)[!keep]
  use <- if (anyNA(M)) "pairwise.complete.obs" else "everything"
  C <- stats::cor(M[, keep, drop = FALSE], use = use)

  # Under pairwise-complete a cell is NA when two signals share no measured
  # stay. It should not happen at 41,000 stays, but eigen() dies on an NA and
  # this runs AFTER 215 fits — so drop the worst offender until the matrix is
  # clean rather than lose the run at the last step. Anything dropped here is
  # reported, never silent.
  while (anyNA(C) && ncol(C) > 2L) {
    j <- which.max(colSums(is.na(C)))
    dropped <- c(dropped, colnames(C)[j])
    C <- C[-j, -j, drop = FALSE]
  }
  attr(C, "dropped") <- dropped
  attr(C, "use") <- use
  C
}

#' The eigenspectrum of a correlation matrix.
#'
#' NEGATIVE EIGENVALUES ARE REPORTED, NOT CLIPPED (statistical review S9,
#' 2026-09-09). A pairwise-complete correlation matrix -- the NA-filled
#' sensitivity reading -- estimates every cell on a different stay population
#' and need not be positive semidefinite, and the deficit is not confined to
#' roundoff: a six-row synthetic missingness pattern gives eigenvalues 2, 2, -1.
#' Until now every negative eigenvalue was replaced by zero and the rest were
#' renormalised, so the clipped spectrum summed to 4 where the trace was 3 and
#' the table looked like an ordinary PCA. `prop` is now the eigenvalue over the
#' TRACE, which is the number of signals for a correlation matrix and is
#' preserved whether or not the matrix is definite, so the shares still sum to
#' one and a negative share is visible as such. The attributes carry the
#' minimum eigenvalue and the negative spectral mass as a fraction of the
#' trace; `spectrum_summary()` puts both in its row. Read an indefinite
#' matrix's `pc1` as a redundancy summary of an inconsistent estimate, never as
#' an explained-variance fraction. The zero-filled primary matrix is complete
#' and therefore definite; this concerns the pairwise sensitivity arm.
#'
#' @return data frame: component, eigenvalue, prop, cumprop, with attributes
#'   `min_eigenvalue`, `neg_mass`, `psd`
eigenspectrum <- function(C) {
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  tr <- sum(ev)
  out <- data.frame(component = seq_along(ev),
                    eigenvalue = round(ev, 5),
                    prop = round(ev / tr, 5),
                    cumprop = round(cumsum(ev) / tr, 5),
                    stringsAsFactors = FALSE)
  attr(out, "min_eigenvalue") <- min(ev)
  attr(out, "neg_mass") <- sum(-ev[ev < 0]) / tr
  attr(out, "psd") <- min(ev) >= -1e-8 * max(abs(ev))
  out
}

#' The branch point, as one row of numbers.
#'
#' `pc1` is the headline: the share of variance in the first component. But a
#' single number is a poor summary of "does one component dominate", so three
#' others are reported beside it, and they disagree in informative ways:
#'
#'   n80 / n90   how many components to reach 80% / 90%. Blunt but readable.
#'   pr          participation ratio, (sum L)^2 / sum(L^2). A CONTINUOUS
#'               effective dimensionality — it does not depend on a threshold,
#'               and it is the one to quote if only one is quoted. pr near 1
#'               means one direction; pr near S means isotropic.
#'   kaiser      components with eigenvalue > 1. On a correlation matrix that is
#'               "explains more than one signal's worth of variance". Included
#'               because it is conventional, not because it is better.
#'
#' NOT a decision rule. There is no threshold at which layer 2 becomes
#' motivated; the numbers inform a judgement that has to be made and stated.
spectrum_summary <- function(C, label = NA_character_) {
  e <- eigenspectrum(C)
  ev <- e$eigenvalue
  data.frame(
    label      = label,
    n_signals  = ncol(C),
    n_dropped  = length(attr(C, "dropped") %||% character(0)),
    use        = attr(C, "use") %||% NA_character_,
    pc1        = e$prop[1],
    pc12       = round(sum(e$prop[1:min(2, length(ev))]), 5),
    n80        = which(e$cumprop >= 0.80)[1],
    n90        = which(e$cumprop >= 0.90)[1],
    pr         = round(sum(ev)^2 / sum(ev^2), 3),
    kaiser     = sum(ev > 1),
    mean_abs_r = round(mean(abs(C[upper.tri(C)]), na.rm = TRUE), 4),
    max_abs_r  = round(max(abs(C[upper.tri(C)]), na.rm = TRUE), 4),
    # The definiteness reading (review S9). `psd` FALSE means the row above
    # summarises an indefinite estimate: `pc1`, `pr` and `n80` are then
    # descriptions of that spectrum, not variance fractions of a covariance.
    min_eig    = round(attr(e, "min_eigenvalue"), 5),
    neg_mass   = round(attr(e, "neg_mass"), 5),
    psd        = isTRUE(attr(e, "psd")),
    stringsAsFactors = FALSE)
}

# --- the reading -------------------------------------------------------------

#' Collapse the GCS triple to one representative.
#'
#' The three GCS models share nine sedation terms, so three columns of the L
#' matrix are built partly from the same covariates. That could inflate PC1 for
#' a DESIGN reason rather than a physiological one, which would be the wrong
#' basis for a branch-point decision.
#'
#' `L_cond` addresses the same worry more directly — it removes the shared
#' intervention term by construction — so the two readings are complementary,
#' not alternatives: if collapsing GCS moves the spectrum but L_cond does not,
#' the shared terms are not what is driving it.
gcs_collapse <- function(M, keep = "gcs_motor") {
  gcs <- c("gcs_motor", "gcs_eyes", "gcs_verbal")
  drop <- setdiff(intersect(gcs, colnames(M)), keep)
  if (!length(drop)) return(M)
  M[, setdiff(colnames(M), drop), drop = FALSE]
}

#' The whole branch-point reading, as one table.
#'
#' Every combination of {full, cond, intv} x {zero, na} x {19 signals, GCS
#' collapsed}, which is the 2x3 of CLAUDE.md's frozen decision widened by the
#' sensitivity reading in next_steps.md. All of it is cor() and eigen() once the
#' fits exist — seconds, no refitting.
#'
#' @param mats OPTIONAL named list `list(zero = , na = )` of already-built
#'   matrix sets. Supply them and this function pivots nothing.
#'
#'   WHY THE ARGUMENT EXISTS. Without it this function calls `l_matrices()` once
#'   per fill, which recomputes exactly what the `l_mats_zero` and `l_mats_na`
#'   targets already hold -- twelve pivots of a 41,250 x 19 matrix per run since
#'   the interaction models took the set from four matrices to six. The wasted
#'   half-minute is the small part. The real cost is that "the L matrices" was
#'   computed in two places, and nothing asserted the two agreed: the branch
#'   point could have described a different matrix from the one every other
#'   target reads, with no symptom, the moment the two call sites drifted in an
#'   argument. Passing them in makes the identity structural instead.
#'
#'   The fallback path is kept because `tests/branch_point.R` calls this with a
#'   long table and no target store.
#' @return data frame, one row per reading, most-dominant PC1 first
branch_point <- function(l_long, tabs, cfg, stay_ids, models = c("full", "cond", "intv"),
                         mats = NULL) {
  if (!is.null(mats)) {
    miss <- setdiff(c("zero", "na"), names(mats))
    if (length(miss)) {
      abort_values("branch_point: `mats` must carry both fills", miss)
    }
  }
  rows <- list()
  for (fl in c("zero", "na")) {
    mm <- if (is.null(mats)) l_matrices(l_long, tabs, cfg, stay_ids, fill = fl)
          else mats[[fl]]
    # ASSERTED, NOT ASSUMED. A supplied set must carry the fill it claims, or
    # the table would label an NA-filled reading "zero" and every PC1 in the
    # co-missingness half of the branch point would be the wrong number under
    # the right name.
    if (!is.null(mats) && !identical(attr(mm, "fill"), fl)) {
      stop("branch_point: the supplied `mats$", fl, "` is stamped fill = ",
           attr(mm, "fill") %||% "<none>", call. = FALSE)
    }
    mats_fl <- mm
    for (md in models) {
      M <- mats_fl[[md]]
      if (is.null(M)) next
      for (cl in c(FALSE, TRUE)) {
        Mi <- if (cl) gcs_collapse(M) else M
        C <- l_correlation(Mi)
        if (ncol(C) < 2L) next
        s <- spectrum_summary(C, label = sprintf("%s/%s/%s", md, fl,
                                                 if (cl) "gcs1" else "gcs3"))
        s$model <- md; s$fill <- fl; s$gcs <- if (cl) "collapsed" else "all three"
        rows[[length(rows) + 1L]] <- s
      }
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(-out$pc1), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Print the branch point in the form the decision is actually made from.
#'
#' Says what the numbers are and what they do not settle. Deliberately refuses
#' to emit a verdict: "does one component dominate" is a judgement about whether
#' 19 marginal models are measuring 19 things or one thing, and no threshold
#' decides it.
report_branch_point <- function(bp, top_n = 3L) {
  cat("\n=== BRANCH POINT: does one component dominate? ===\n\n")
  print(bp[, c("label", "n_signals", "pc1", "pc12", "n80", "n90", "pr",
               "kaiser", "mean_abs_r")], row.names = FALSE)

  f <- bp[bp$model == "full" & bp$fill == "zero" & bp$gcs == "all three", ]
  if (nrow(f) == 1L) {
    cat(sprintf("\n  primary reading (L_full, zero-filled, 19 signals):\n"))
    cat(sprintf("    PC1 = %.1f%% of variance, participation ratio = %.2f of %d\n",
                100 * f$pc1, f$pr, f$n_signals))
    cat(sprintf("    %d component(s) reach 80%%, %d reach 90%%\n", f$n80, f$n90))
  }
  cat("\n  Read the SPREAD across rows, not any single row. A PC1 that moves\n")
  cat("  between zero-fill and pairwise is co-missingness, not physiology.\n")
  cat("  A PC1 that moves between L_full and L_cond is the shared intervention\n")
  cat("  term, which is exactly what Sigma-inverse is meant to discount.\n")
  cat("  No threshold decides this. State the judgement and the numbers.\n\n")
  invisible(bp)
}

# --- redundancy, measured rather than asserted -------------------------------

#' How near-collinear are the columns that share covariates?
#'
#' The GCS triple share nine sedation terms; creatinine and urine output share
#' `diuretic`. CLAUDE.md's position is that this is inherent to a set of
#' MARGINAL models and that Sigma-inverse is what discounts it. With L_intv
#' computed explicitly that stops being an assertion: the three GCS `intv`
#' columns should be near-collinear, and how near is a number.
#'
#' @return data frame: pair, r_full, r_cond, r_intv
shared_covariate_pairs <- function(l_long, tabs, cfg, stay_ids) {
  mats <- l_matrices(l_long, tabs, cfg, stay_ids, fill = "zero")
  pairs <- list(
    c("gcs_motor", "gcs_eyes"), c("gcs_motor", "gcs_verbal"),
    c("gcs_eyes", "gcs_verbal"),
    c("creatinine", "urine_output_rate"),
    c("mbp", "heart_rate"),
    c("spo2", "resp_rate")
  )
  rows <- lapply(pairs, function(p) {
    r <- function(md) {
      M <- mats[[md]]
      if (is.null(M) || !all(p %in% colnames(M))) return(NA_real_)
      a <- M[, p[1]]; b <- M[, p[2]]
      if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
      round(stats::cor(a, b), 4)
    }
    data.frame(pair = paste(p, collapse = " ~ "),
               r_full = r("full"), r_cond = r("cond"), r_intv = r("intv"),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # The interesting column is the drop from full to cond: it is the part of the
  # correlation that was the shared intervention block rather than physiology.
  out$drop_full_to_cond <- round(out$r_full - out$r_cond, 4)
  out
}

# --- summaries safe to print -------------------------------------------------

#' Per-signal summary of an L matrix. Aggregates only (hard rule 1).
l_summary <- function(M) {
  data.frame(
    signal   = colnames(M),
    n        = apply(M, 2, function(v) sum(!is.na(v))),
    n_zero   = apply(M, 2, function(v) sum(v == 0, na.rm = TRUE)),
    mean     = round(apply(M, 2, mean, na.rm = TRUE), 4),
    sd       = round(apply(M, 2, stats::sd, na.rm = TRUE), 4),
    q05      = round(apply(M, 2, stats::quantile, 0.05, na.rm = TRUE), 4),
    q95      = round(apply(M, 2, stats::quantile, 0.95, na.rm = TRUE), 4),
    stringsAsFactors = FALSE, row.names = NULL)
}
