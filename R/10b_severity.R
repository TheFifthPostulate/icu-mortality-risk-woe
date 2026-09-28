# R/10b_severity.R -----------------------------------------------------------
# The severity comparators as a TRANSPORTABLE object, and the one apply path
# both sites use for them.
#
# WHY THIS FILE EXISTS. R/09c and R/09d already compute APACHE II and SOFA, and
# neither of them fits anything: the point tables are frozen in code, so a score
# recomputed at eICU is the same construction as the one recomputed at MIMIC.
# That is the whole reason those tables live in R rather than in SQL.
#
# ONE THING IN THE ARM IS FITTED, AND IT IS THE REASON THIS FILE IS NEEDED.
# An integer severity score is not on the log-odds scale, so every cell is
# mapped onto it by a logistic recalibration before it can be compared against
# an L. At the training site `recalibrate_oof()` does that fold-wise. AN APPLY
# SITE HAS NO FOLDS, and re-deriving the mapping from the apply site's own
# outcomes would fit at eICU -- which is hard rule 8, and it would flatter the
# baseline rather than us. The transportable object is therefore the
# MIMIC-fitted intercept and slope, frozen into the bundle here and applied
# unchanged at both apply sites.
#
# THE RESTRICTION IS PART OF THE FROZEN OBJECT TOO. A stay scored on four of
# twelve APACHE variables is a fragment, not a severity score, so
# `apache_complete_mask()` and `sofa_complete_mask()` gate the comparison. The
# floors are declared in `config/apache` and `config/sofa` before any run and
# are frozen into the bundle's design, so an apply site cannot quietly use a
# different one. The mask is applied to EVERY cell including our own arms:
# comparing a restricted baseline against an unrestricted proposed method is a
# different and much weaker experiment.
#
# WHAT IS DELIBERATELY NOT TRANSPORTED.
#
#   the native scores  Each site ships its own (eICU has APACHE IVa / APS;
#                      MIMIC-IV ships none), so there is no shared quantity to
#                      freeze. They are reported RAW, as a rank-only sanity
#                      check, and labelled site-local.
#   llr_plus_demo      The age / chronic-health cell in tests/metrics_severity.R
#                      needs a fold-wise glm over `admission_class`, whose factor
#                      levels are not guaranteed to match across sites. It stays
#                      a MIMIC-internal descriptive cell.
#
# AGGREGATES ONLY (hard rule 1). Everything returned here that is row-level --
# the scores, the mask, the per-stay point frames -- goes to a run directory.
# The diagnostic tables are counts and means and are safe to print.
# ----------------------------------------------------------------------------

#' The severity cells that are frozen and travel.
#'
#' `apache2_aps` is the acute physiology score alone and is the honest
#' counterpart to `llr_meas`, which carries no intervention term either.
#' `apache2_total` adds the age and chronic-health points. `sofa` is the
#' counterpart to `llr_sum`, because SOFA bundles vasopressor dose into its
#' cardiovascular component and is therefore a measurement-intervention bundle
#' rather than pure physiology (docs/severity_baselines.md).
SEVERITY_CELLS <- c("apache2_aps", "apache2_total", "sofa")

#' Pull the severity settings out of a config, with the same defaults
#' `tests/metrics_severity.R` uses.
#'
#' Read from the BUNDLE's frozen design at an apply site, never from
#' `config/config.yml` -- `bundle_cfg()` is what makes that true, and `apache`
#' and `sofa` are in `BUNDLE_DESIGN_KEYS` so that it can be.
severity_settings <- function(cfg) {
  list(
    min_vars     = as.integer(cfg$apache$min_vars_present %||% 10L),
    arf_doubling = isTRUE(cfg$apache$arf_doubling),
    gcs_source   = as.character(cfg$apache$gcs_source %||% "native"),
    min_organs   = as.integer(cfg$sofa$min_organs_present %||% 4L),
    resp_support = as.character(cfg$sofa$resp_support %||% "invasive"),
    organs       = as.character(unlist(cfg$sofa$organs %||% SOFA_ORGANS))
  )
}

#' Recompute both severity scores on a fixed stay set. FITS NOTHING.
#'
#' Identical at both sites by construction: the point tables are in code, the
#' settings come from the frozen design, and the only site-specific input is the
#' parquet. This is called once inside the training graph to fit the
#' recalibration, and once at each apply site to score with it.
#'
#' @param restrict apply the coverage floors. FALSE returns an all-TRUE mask so
#'   the unrestricted numbers can be reported beside the restricted ones -- the
#'   two bracket the truth in opposite directions and neither is primary alone.
#' @return list with the aligned severity rows, the two score frames, the mask,
#'   and the raw scores per cell. Row-level; never print it.
severity_raw <- function(tabs, cfg, stay_ids, settings = NULL, restrict = TRUE) {
  st <- settings %||% severity_settings(cfg)
  if (is.null(tabs$severity)) {
    stop("severity_raw: no `severity` table. Add `severity:` to this site's ",
         "`paths` block and re-run the extraction.", call. = FALSE)
  }
  ids <- as.character(stay_ids)
  ci  <- match(ids, as.character(tabs$cohort$stay_id))
  ai  <- match(ids, as.character(tabs$severity$stay_id))
  if (anyNA(ci)) stop("severity_raw: a scored stay has no cohort row", call. = FALSE)
  if (anyNA(ai)) stop("severity_raw: a scored stay has no severity row", call. = FALSE)

  ap  <- tabs$severity[ai, , drop = FALSE]
  age <- tabs$cohort$age[ci]

  ap2 <- apache2_score(ap, age = age, arf_doubling = st$arf_doubling,
                       gcs_source = st$gcs_source)
  sf  <- sofa_score(ap, resp_support = st$resp_support, organs = st$organs)

  # `sf` rather than `ap` for the SOFA mask: the extraction's
  # `sofa_n_organs_present` counts six, while `sofa_score()` recounts over the
  # declared organ set, and respiration is excluded from that set. Using the raw
  # column would apply a 4-of-6 floor to a five-organ score.
  # THE SCORE'S OWN VARIABLE COUNT CONTROLS APACHE INCLUSION, as of 2026-09-09
  # (statistical review, severity coverage). `apache2_score()` recounts
  # `n_vars_present` under the GCS and oxygenation rules actually applied, and
  # the extraction's `ap2_n_vars_present` can differ from it when those rules
  # change. The floor is applied to the count the scored points were built
  # from -- the same policy `sofa_complete_mask()` already follows with `sf` --
  # and the number of stays on which the two counts disagree is reported in
  # `severity_diagnostics()` rather than assumed to be zero.
  n_ext <- ap[["ap2_n_vars_present"]]
  n_sc  <- ap2$n_vars_present
  n_vars_disagree <- if (is.null(n_ext)) NA_integer_ else
    sum(!is.na(n_ext) & !is.na(n_sc) & n_ext != n_sc)
  keep <- if (isTRUE(restrict)) {
    apache_complete_mask(data.frame(ap2_n_vars_present = n_sc),
                         min_vars = st$min_vars) &
      sofa_complete_mask(sf, min_organs = st$min_organs)
  } else rep(TRUE, length(ids))

  # The length assertion is not defensive clutter. A mask built from a missing
  # column is zero-length, recycles to nothing, reports "0 dropped", and makes
  # every downstream cell declare itself entirely NA far away from the cause.
  # See the companion note in `sofa_complete_mask()`, which is where that
  # actually happened.
  if (length(keep) != length(ids)) {
    stop(sprintf("severity_raw: coverage mask is length %d against %d stays.",
                 length(keep), length(ids)), call. = FALSE)
  }
  if (!any(keep)) stop("severity_raw: the coverage floors retained no stays", call. = FALSE)

  nat <- function(v) if (!is.null(ap[[v]])) as.numeric(ap[[v]]) else rep(NA_real_, length(ids))

  list(stay_id  = ids,
       settings = st,
       ap       = ap,
       ap2      = ap2,
       sf       = sf,
       keep     = keep,
       n_drop   = sum(!keep),
       n_vars_disagree = n_vars_disagree,
       raw      = list(apache2_aps   = as.numeric(ap2$aps),
                       apache2_total = as.numeric(ap2$total),
                       sofa          = as.numeric(sf$sofa)),
       native   = list(aps_native  = nat("aps_native"),
                       sofa_native = nat("sofa_native")))
}

# --- the fitted half, which happens once, at MIMIC ---------------------------

#' Fit the transportable recalibration: one intercept and slope per cell.
#'
#' `glm(y ~ raw)` on the FULL training set, restricted rows only. The full-train
#' analogue of what `recalibrate_oof()` does fold-wise, and the object an apply
#' site cannot re-derive without fitting on its own outcomes.
#'
#' `p_bar` is carried beside the coefficients and is the event rate of the
#' RESTRICTED TRAINING rows, not of the cohort and not of the apply site. It is
#' the constant every recalibrated score is centred on, so that a severity cell
#' and an L arrive at an apply site centred on the same kind of quantity.
#'
#' @return data.frame, one row per cell, plus the attributes the apply path
#'   needs. A non-positive slope is recorded rather than repaired: it would mean
#'   the score runs the wrong way at the training site, which is a finding.
fit_severity_recal <- function(sev, y, cells = SEVERITY_CELLS) {
  if (length(y) != length(sev$keep)) {
    stop("fit_severity_recal: `y` is not aligned to the severity rows", call. = FALSE)
  }
  k  <- sev$keep
  yk <- as.integer(y[k])
  p_bar <- mean(yk)
  rows <- lapply(cells, function(nm) {
    x <- sev$raw[[nm]][k]
    ok <- !is.na(x) & !is.na(yk)
    if (sum(ok) < 100L) {
      return(data.frame(cell = nm, intercept = NA_real_, slope = NA_real_,
                        n = sum(ok), converged = FALSE, stringsAsFactors = FALSE))
    }
    fit <- stats::glm(yk[ok] ~ x[ok], family = stats::binomial())
    cf  <- stats::coef(fit)
    data.frame(cell = nm, intercept = unname(cf[1]), slope = unname(cf[2]),
               n = sum(ok), converged = isTRUE(fit$converged),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "p_bar") <- p_bar
  out
}

#' Apply a frozen recalibration to one cell's raw score.
#'
#' Deliberately tiny and deliberately not vectorised over cells: every call site
#' names the cell it is mapping, so a mis-join between coefficient rows and
#' score vectors is a missing-row error rather than a silent recycle.
recal_apply <- function(recal, cell, x, p_bar) {
  i <- match(cell, recal$cell)
  if (is.na(i)) {
    abort_values("recal_apply: the bundle carries no recalibration for cell", cell)
  }
  if (!is.finite(recal$intercept[i]) || !is.finite(recal$slope[i])) {
    stop("recal_apply: cell `", cell, "` has no fitted recalibration. An absent ",
         "mapping must fail rather than pass the raw integer score through on ",
         "a log-odds scale.", call. = FALSE)
  }
  recal$intercept[i] + recal$slope[i] * x - logit(p_bar)
}

#' Assemble the severity slot the bundle carries.
#'
#' A plain constructor: it stores, it does not estimate. Everything in it was
#' produced by a target upstream, which is what lets `verify_severity()` be a
#' pure check.
make_severity_bundle <- function(settings, recal, p_bar, n_train, frac_kept,
                                 cutpoints = list(), train_ref = NULL,
                                 symmetric = NULL) {
  list(settings  = settings,
       symmetric = symmetric,
       recal     = recal,
       p_bar     = p_bar,
       n_train   = n_train,
       frac_kept = frac_kept,
       cutpoints = cutpoints,
       train_ref = train_ref)
}

#' Check rows for the severity slot, in the shape `verify_bundle()` returns.
#'
#' Absent is legitimate and is reported as a skip rather than a failure: a
#' bundle built before the severity extraction existed is a valid bundle, and
#' the apply runners already refuse to score a severity arm that is not there.
verify_severity <- function(bundle, cfg = NULL) {
  row <- function(check, ok, detail) data.frame(check = check, ok = isTRUE(ok),
                                                detail = detail, stringsAsFactors = FALSE)
  sv <- bundle$severity
  if (is.null(sv)) {
    return(row("severity", TRUE, "absent (severity arms not frozen into this bundle)"))
  }
  rc <- sv$recal
  miss <- setdiff(SEVERITY_CELLS, rc$cell)
  fin  <- is.finite(rc$intercept) & is.finite(rc$slope)
  pos  <- fin & rc$slope > 0
  conv <- vapply(rc$converged %||% rep(NA, nrow(rc)), isTRUE, logical(1))
  out <- list(
    row("severity_recal", !length(miss) && all(fin),
        sprintf("%d/%d cell(s)%s", sum(fin), length(SEVERITY_CELLS),
                if (length(miss)) paste0("; missing: ", paste(miss, collapse = ", ")) else "")),
    # Finite coefficients are not acceptance (statistical review S10). The
    # full-train recalibration records its convergence flag; a bundle whose
    # cell did not converge fails verification rather than travelling.
    row("severity_recal_converged", all(conv),
        if (all(conv)) "all cells converged"
        else paste("not converged:", paste(rc$cell[!conv], collapse = ", "))),
    # A negative slope means the severity score runs the wrong way at the
    # training site. That is a finding about the extraction, not something to
    # repair here, so it is reported as a failed check and the run continues.
    row("severity_slope_sign", all(pos),
        if (all(pos)) "all slopes positive"
        else paste("non-positive slope:", paste(rc$cell[!pos], collapse = ", "))),
    row("severity_p_bar", is.finite(sv$p_bar) && sv$p_bar > 0 && sv$p_bar < 1,
        sprintf("%.4f on %d restricted train stays (%.3f of train)",
                sv$p_bar %||% NA_real_, sv$n_train %||% NA_integer_,
                sv$frac_kept %||% NA_real_))
  )
  if (!is.null(cfg)) {
    want <- severity_settings(cfg)
    same <- identical(want, sv$settings)
    out[[length(out) + 1L]] <- row(
      "severity_settings_match_local", same,
      if (same) "identical"
      else paste("frozen:", paste(names(sv$settings), unlist(lapply(sv$settings, paste,
                 collapse = "/")), sep = "=", collapse = " ")))
  }
  do.call(rbind, out)
}

# --- the apply path ----------------------------------------------------------

#' Score the severity arms at an apply site with a frozen bundle. FITS NOTHING.
#'
#' The severity counterpart to `apply_bundle()`, and the same argument for its
#' existence: MIMIC-test and eICU go through this one function, so "the same
#' code path" is a property of the code rather than a claim about two runner
#' scripts that happen to look alike.
#'
#' @return list(scores, keep, raw, native, sev). `scores` and `keep` are
#'   row-level and aligned to `stay_ids`.
apply_severity <- function(bundle, tabs, cfg, stay_ids, verbose = TRUE) {
  sv <- bundle$severity
  if (is.null(sv)) {
    stop("apply_severity: this bundle carries no severity slot. Rebuild the ",
         "internal graph with the severity targets before scoring the APACHE ",
         "II / SOFA arms; an absent arm must fail rather than silently shorten ",
         "the comparison table.", call. = FALSE)
  }
  stamp <- cfg$.bundle_design
  if (is.null(stamp) || !identical(stamp, .hash(bundle$cfg))) {
    stop("apply_severity: `cfg` did not come from bundle_cfg(bundle, paths). ",
         "The coverage floors and the GCS source are part of the frozen design ",
         "(hard rule 8).", call. = FALSE)
  }
  if (verbose) message("apply_severity: ", length(stay_ids), " stays")

  sev <- severity_raw(tabs, cfg, stay_ids, settings = sv$settings, restrict = TRUE)
  sc  <- lapply(SEVERITY_CELLS, function(nm)
    stats::setNames(recal_apply(sv$recal, nm, sev$raw[[nm]], sv$p_bar), sev$stay_id))
  names(sc) <- SEVERITY_CELLS

  list(scores = sc, keep = sev$keep, raw = sev$raw, native = sev$native, sev = sev)
}

# --- the tables that have to sit beside the numbers --------------------------

#' Coverage, point contributions and the native sanity checks, in one call.
#'
#' All aggregates. The coverage table is the one that has to be read at both
#' sites side by side: if the baseline is thinner at one site than the other,
#' part of any cross-site difference is the baseline's coverage rather than the
#' method's transportability.
severity_diagnostics <- function(sev) {
  k <- sev$keep
  ap2 <- sev$ap2; sf <- sev$sf
  pt_cols   <- grep("^pt_", names(ap2), value = TRUE)
  org_cols  <- grep("^sofa_(resp|coag|liver|cardio|cns|renal)", names(sf), value = TRUE)

  points <- data.frame(
    variable     = sub("^pt_", "", pt_cols),
    mean_points  = round(vapply(pt_cols, function(v) mean(ap2[[v]][k], na.rm = TRUE), numeric(1)), 4),
    frac_nonzero = round(vapply(pt_cols, function(v) mean(ap2[[v]][k] > 0, na.rm = TRUE), numeric(1)), 4),
    stringsAsFactors = FALSE, row.names = NULL)

  organs <- data.frame(
    organ        = sub("^sofa_", "", org_cols),
    mean_score   = round(vapply(org_cols, function(v) mean(sf[[v]][k], na.rm = TRUE), numeric(1)), 4),
    frac_nonzero = round(vapply(org_cols, function(v) mean(sf[[v]][k] > 0, na.rm = TRUE), numeric(1)), 4),
    stringsAsFactors = FALSE, row.names = NULL)

  natv <- sev$native$aps_native
  agree <- if (any(!is.na(natv))) {
    apache_native_agreement(sev$raw$apache2_total[k], natv[k])
  } else NULL
  sofa_agree <- sofa_native_agreement(sf[k, , drop = FALSE], sev$ap[k, , drop = FALSE],
                                      scored = sev$settings$organs)

  list(coverage = apache_coverage(sev$ap),
       points = points, organs = organs,
       apache_vs_native = agree, sofa_vs_native = sofa_agree,
       kept = data.frame(n = length(k), n_kept = sum(k), n_drop = sum(!k),
                         frac_kept = round(mean(k), 5),
                         min_vars = sev$settings$min_vars,
                         # Stays where the extraction's APACHE variable count
                         # and the score's own recount disagree. The recount
                         # controls inclusion; this says how often it mattered.
                         n_vars_count_disagree = sev$n_vars_disagree %||% NA_integer_,
                         min_organs = sev$settings$min_organs,
                         gcs_source = sev$settings$gcs_source,
                         resp_support = sev$settings$resp_support,
                         organs = paste(sev$settings$organs, collapse = "+"),
                         stringsAsFactors = FALSE))
}

# --- the arm, as one function both apply sites call --------------------------

#' The pre-registered severity contrasts.
#'
#' APACHE II is pure physiology, so its counterpart is `llr_meas`. SOFA bundles
#' vasopressor dose into its cardiovascular component, so its counterpart is
#' `llr_sum`. THE PAIRING IS NOT INTERCHANGEABLE: swapping them would compare a
#' physiology score against a measurement-intervention bundle and attribute the
#' difference to the method rather than to the feature set
#' (docs/severity_baselines.md, and CLAUDE.md on why SOFA is lineage rather than
#' a third baseline).
#'
#' `apache2_total` against `llr_meas` is the third row and is the harder of the
#' two APACHE contrasts, because the total adds age and chronic-health points
#' that no L carries. It is reported so the margin cannot be read as depending
#' on the softer comparison.
SEVERITY_CONTRASTS <- list(
  c("llr_meas", "apache2_aps"),
  c("llr_meas", "apache2_total"),
  c("llr_sum",  "sofa")
)

#' Score the APACHE II / SOFA comparison at an apply site, end to end.
#'
#' ONE FUNCTION, BOTH SITES. `run/test_look.R` and `run/external.R` call this
#' with different rows and nothing else differs, which is the same argument
#' `apply_bundle()` makes one level down: "the same code path" becomes a
#' property of the code rather than a claim about two scripts that look alike.
#'
#' THE COVERAGE RESTRICTION IS APPLIED TO EVERY CELL, our own arms included. A
#' restricted baseline against an unrestricted proposed method is a different
#' and much weaker experiment, and the shift in `llr_sum` between this table and
#' the unrestricted one is the size of the selection -- a number to report, not
#' to assume away.
#'
#' THE NATIVE SCORES ARE RANK-ONLY. Each site ships its own and there is nothing
#' shared to recalibrate, so they are scored raw and their AUROC and AUPRC are
#' the only cells that mean anything. They are a sanity check on the
#' recomputation, never a result.
#'
#' @param scores  the bundle arms from `apply_bundle()`, aligned to `stay_ids`
#' @param l_full  the zero-filled `full` L matrix, for the domain table. NULL
#'   skips it.
#' @return list(summary, contrasts, transport, keep, scores, diagnostics)
#' @param group   patient id per row of `stay_ids`; the bootstrap unit
#'   (statistical review S4). NULL resamples stays and the tables say so.
severity_arm <- function(run, bundle, tabs, cfg, stay_ids, y, scores, group = NULL,
                         l_full = NULL, domains = NULL,
                         n_bins = 20L, n_boot = 200L, seed = 1L,
                         verbose = TRUE) {
  sev <- apply_severity(bundle, tabs, cfg, stay_ids, verbose = verbose)
  k   <- sev$keep
  if (length(y) != length(k)) {
    stop("severity_arm: `y` is not aligned to `stay_ids`", call. = FALSE)
  }
  yk    <- as.integer(y[k])
  p_bar <- mean(yk)
  if (!is.null(group) && length(group) != length(k)) {
    stop("severity_arm: `group` is not aligned to `stay_ids`", call. = FALSE)
  }
  gk <- if (is.null(group)) NULL else group[k]

  # Every cell on the identical rows. The order is fixed so a table read at two
  # sites has its rows in the same places.
  cells <- c(lapply(scores, function(v) as.numeric(v)[k]),
             lapply(sev$scores, function(v) as.numeric(v)[k]))
  for (nm in names(sev$native)) {
    v <- sev$native[[nm]]
    if (any(!is.na(v[k]))) cells[[nm]] <- v[k]
  }

  log_msg(run, sprintf(paste0("severity arm: %d of %d stays kept (min_vars %d, ",
                              "min_organs %d); %d cell(s); restricted event rate %.4f"),
                       sum(k), length(k), sev$sev$settings$min_vars,
                       sev$sev$settings$min_organs, length(cells), p_bar))

  sc <- score_arms(run, cells, yk, p_bar, breaks = NULL, suffix = "_sev",
                   n_bins = n_bins, n_boot = n_boot, seed = seed, group = gk)
  summ <- sc$summary
  save_table(run, summ, "severity_score_summary")

  # --- the SYMMETRIC restriction ------------------------------------------
  # The cell above restricts to stays where the BASELINE is scoreable. This one
  # additionally requires OUR inputs to be present, which is the same argument
  # applied in the direction it was not. Both are reported; neither is primary
  # alone. See `severity_symmetric` in config for why, and for the eICU Glasgow
  # coverage asymmetry that motivated it.
  sym <- bundle$severity$symmetric %||% list(enabled = FALSE)
  sym_out <- NULL
  if (isTRUE(sym$enabled)) {
    M    <- measured_matrix(tabs, cfg, stay_ids)
    need <- intersect(as.character(unlist(sym$require_signals %||% character(0))),
                      colnames(M))
    miss <- setdiff(as.character(unlist(sym$require_signals %||% character(0))),
                    colnames(M))
    if (length(miss)) {
      # A required signal that is not a column is a design/extraction mismatch,
      # not a stay that fails the rule. Silently treating it as satisfied would
      # turn the restriction off without saying so.
      abort_values("severity_arm: `require_signals` names a signal that is not modelled", miss)
    }
    ok_named <- if (length(need)) rowSums(M[, need, drop = FALSE]) == length(need)
                else rep(TRUE, nrow(M))
    floor_n  <- as.integer(sym$min_signals_measured %||% 0L)
    ok_count <- if (floor_n > 0L) rowSums(M) >= floor_n else rep(TRUE, nrow(M))
    k2 <- k & ok_named & ok_count

    if (!any(k2)) {
      log_msg(run, "symmetric restriction retained no stays; cell skipped")
    } else {
      yk2 <- as.integer(y[k2]); p2 <- mean(yk2)
      gk2 <- if (is.null(group)) NULL else group[k2]
      cells2 <- c(lapply(scores, function(v) as.numeric(v)[k2]),
                  lapply(sev$scores, function(v) as.numeric(v)[k2]))
      log_msg(run, sprintf(paste0("symmetric cell: %d of %d stays (%.1f%% of the ",
                                  "baseline-restricted set); require %s; min_signals %d"),
                           sum(k2), sum(k), 100 * sum(k2) / max(sum(k), 1L),
                           if (length(need)) paste(need, collapse = "+") else "none",
                           floor_n))
      sc2 <- score_arms(run, cells2, yk2, p2, breaks = NULL, suffix = "_symsev",
                        n_bins = n_bins, n_boot = n_boot, seed = seed, group = gk2)
      save_table(run, sc2$summary, "severity_score_summary_symmetric")
      have2 <- names(cells2)
      pr2 <- Filter(function(q) all(q %in% have2), SEVERITY_CONTRASTS)
      ct2 <- if (length(pr2)) arm_contrasts(cells2, yk2, pr2, n_boot = n_boot, seed = seed,
                                            group = gk2) else NULL
      if (!is.null(ct2)) save_table(run, ct2, "severity_contrasts_symmetric")
      sym_out <- list(summary = sc2$summary, contrasts = ct2, keep = k2,
                      n = sum(k2), frac_of_restricted = sum(k2) / max(sum(k), 1L),
                      p_bar = p2, require_signals = need,
                      min_signals_measured = floor_n)
    }
  }

  # Paired, because every cell is scored on identical rows and marginal
  # intervals settle nothing when they are.
  have <- names(cells)
  pairs <- Filter(function(p) all(p %in% have), SEVERITY_CONTRASTS)
  ct <- if (length(pairs)) arm_contrasts(cells, yk, pairs, n_boot = n_boot, seed = seed,
                                         group = gk) else NULL
  if (!is.null(ct)) save_table(run, ct, "severity_contrasts")

  # Against what the training site recorded for itself on ITS restricted rows.
  tr <- bundle$severity$train_ref
  transport <- if (!is.null(tr)) {
    m <- merge(
      data.frame(cell = sub("_sev$", "", summ$label), auroc_site = summ$auroc,
                 auprc_site = summ$auprc, stringsAsFactors = FALSE),
      data.frame(cell = tr$arms$label, auroc_mimic_train = tr$arms$auroc,
                 auprc_mimic_train = tr$arms$auprc, stringsAsFactors = FALSE),
      by = "cell", all.x = TRUE)
    m$d_auroc  <- round(m$auroc_site - m$auroc_mimic_train, 5)
    # AUPRC MOVED WITH IT AS OF 2026-09-07. Both columns were already in the
    # table and only the AUROC difference was taken, which is the failure mode
    # `score_metrics()`'s own header warns about: at a roughly 10% event rate a
    # score can retain its ranking and lose its precision where the deaths are,
    # and an AUROC ratio alone cannot see that. The two ratios disagreeing is
    # the finding rather than a defect.
    #
    # RENAMED 2026-09-09 (external runner review E8). `retained` and
    # `retained_auprc` were raw quotients that credit chance: a cell at AUROC
    # 0.5 against a training 0.8 read as 62.5% "retained". They are now
    # `auroc_ratio` / `auprc_ratio`, and `auroc_retained_above_chance` is the
    # retention measured from 0.5 under the declared near-chance policy. The
    # absolute AUROCs and `d_auroc` remain the primary transport columns.
    m$d_auprc  <- round(m$auprc_site - m$auprc_mimic_train, 5)
    m$auroc_ratio <- round(m$auroc_site / m$auroc_mimic_train, 4)
    m$auroc_retained_above_chance <- retention_above_chance(m$auroc_site, m$auroc_mimic_train)
    m$auprc_ratio <- round(m$auprc_site / m$auprc_mimic_train, 4)
    m
  } else NULL
  if (!is.null(transport)) save_table(run, transport, "severity_transport")

  dg <- severity_diagnostics(sev$sev)
  save_table(run, dg$kept,     "severity_coverage_kept", subdir = "diagnostics")
  save_table(run, dg$coverage, "apache_coverage",        subdir = "diagnostics")
  save_table(run, dg$points,   "apache2_points_by_variable", subdir = "diagnostics")
  save_table(run, dg$organs,   "sofa_score_by_organ",    subdir = "diagnostics")
  if (!is.null(dg$apache_vs_native)) {
    save_table(run, dg$apache_vs_native, "apache2_vs_native", subdir = "diagnostics")
  }
  if (!is.null(dg$sofa_vs_native)) {
    save_table(run, dg$sofa_vs_native, "sofa_vs_native", subdir = "diagnostics")
  }

  # The domain-level table, where each `D_k` meets its SOFA organ. Equal
  # weights, because the layer-2 weight vector does not exist yet; `D_k` reduces
  # to the unweighted sum at w = 1 and the comparison sharpens rather than
  # changes shape when the weights land.
  dom <- NULL
  if (!is.null(l_full) && !is.null(domains)) {
    Mk  <- l_full[k, , drop = FALSE]
    dom <- sofa_domain_table(Mk, sev$sev$sf[k, , drop = FALSE], yk, domains,
                             organs = sev$sev$settings$organs)
    save_table(run, dom, "domain_vs_sofa_organ")
    save_table(run, sofa_uncovered_domains(Mk, yk, domains),
               "domains_sofa_cannot_express")
  }

  list(summary = summ, contrasts = ct, transport = transport,
       keep = k, scores = sev$scores, raw = sev$raw,
       diagnostics = dg, domains = dom, p_bar = p_bar, symmetric = sym_out)
}

#' Print the severity arm. Aggregates only.
report_severity_arm <- function(sa, site_label = "this site") {
  cat("\n=== THE SEVERITY ARM: APACHE II and SOFA, ", site_label, " ===\n", sep = "")
  kp <- sa$diagnostics$kept
  cat(sprintf(paste0("\n  %d of %d stays kept (%.1f%%) by the pre-declared coverage\n",
                     "  floors: APACHE >= %d of 12 variables AND SOFA >= %d organs.\n",
                     "  The restriction is applied to EVERY cell, ours included.\n\n"),
              kp$n_kept, kp$n, 100 * kp$frac_kept, kp$min_vars, kp$min_organs))
  print(sa$summary[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
                       "auprc", "auprc_lift")], row.names = FALSE)

  if (!is.null(sa$contrasts)) {
    cat("\n=== the pre-registered contrasts, paired on identical rows ===\n")
    cat("  APACHE II is pure physiology so it meets `llr_meas`; SOFA bundles\n")
    cat("  vasopressor dose so it meets `llr_sum`. The pairing is not\n")
    cat("  interchangeable. Interval and `auroc_p` are the paired bootstrap\n")
    cat("  resampled by `boot_unit`; `delong_p` is the iid reference, a\n")
    cat("  different assumption, printed separately below.\n\n")
    print(sa$contrasts[, c("a", "b", "auroc_a", "auroc_b", "d_auroc",
                           "auroc_lo", "auroc_hi", "auroc_p", "boot_unit",
                           "d_auprc", "auprc_p")], row.names = FALSE)
    cat("\n  DeLong, observation-level iid (sensitivity, not cluster-adjusted):\n\n")
    print(sa$contrasts[, c("a", "b", "d_auroc", "delong_se", "delong_z", "delong_p")],
          row.names = FALSE)
  }
  if (!is.null(sa$transport)) {
    cat("\n=== against MIMIC train, out of fold, on ITS restricted rows ===\n")
    cat("  `auroc_ratio` is the raw quotient and credits chance;\n")
    cat("  `auroc_retained_above_chance` measures from 0.5 and is NA when the\n")
    cat("  training cell is within 0.05 of chance. Read `d_auroc` first.\n\n")
    print(sa$transport[, c("cell", "auroc_mimic_train", "auroc_site", "d_auroc",
                           "auroc_ratio", "auroc_retained_above_chance",
                           "auprc_mimic_train", "auprc_site", "auprc_ratio")],
          row.names = FALSE)
  }
  if (!is.null(sa$symmetric)) {
    sm <- sa$symmetric
    cat("\n=== THE SYMMETRIC CELL: both sides' inputs required ===\n")
    cat(sprintf(paste0("\n  %d stays (%.1f%% of the baseline-restricted set).",
                       " Additionally requires\n  %s to be measured%s.\n"),
                sm$n, 100 * sm$frac_of_restricted,
                if (length(sm$require_signals)) paste(sm$require_signals, collapse = ", ") else "nothing",
                if (sm$min_signals_measured > 0)
                  sprintf(", and at least %d of 19 signals", sm$min_signals_measured) else ""))
    cat("  Read it BESIDE the cell above, never instead of it: this one selects\n")
    cat("  the hospitals that chart our inputs and is hospital-confounded, while\n")
    cat("  that one lets the baseline use inputs we do not have. The two bracket\n")
    cat("  the truth in opposite directions.\n\n")
    print(sm$summary[, c("label", "n", "n_events", "auroc", "auroc_lo",
                         "auroc_hi", "auprc")], row.names = FALSE)
    if (!is.null(sm$contrasts)) {
      cat("\n  the same pre-registered contrasts, on those stays:\n\n")
      print(sm$contrasts[, c("a", "b", "auroc_a", "auroc_b", "d_auroc",
                             "auroc_lo", "auroc_hi", "auroc_p", "boot_unit",
                             "delong_p")], row.names = FALSE)
    }
  }

  cat("\n=== APACHE II input coverage here ===\n")
  cat("  Read beside the other site's. If the baseline is thinner at one site,\n")
  cat("  part of any cross-site difference is the BASELINE's coverage rather\n")
  cat("  than the method's transportability.\n\n")
  print(sa$diagnostics$coverage[, c("variable", "frac_present")], row.names = FALSE)
  if (!is.null(sa$domains)) {
    cat("\n=== domain-level: each D_k against its SOFA organ ===\n")
    cat("  AUPRC beside AUROC because a SOFA organ is a 0-4 integer: it splits\n")
    cat("  the cohort into five blocks and cannot order within one, which AUROC\n")
    cat("  averages over and AUPRC does not. `signals` is in the saved table.\n\n")
    print(sa$domains[, c("sofa_organ", "domain", "n_signals", "auroc_D",
                         "auroc_organ", "delta", "auprc_D", "auprc_organ",
                         "delta_auprc", "auprc_lift_D", "auprc_lift_organ",
                         "spearman")], row.names = FALSE)
  }
  invisible(sa)
}
