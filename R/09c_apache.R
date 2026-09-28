# R/09c_apache.R -------------------------------------------------------------
# The APACHE comparison arm: a baseline for the FEATURE SET, not for the model.
#
# WHAT THIS ARM ANSWERS, and why it is not the XGBoost arm again. R/09b holds
# the learner's information fixed and varies the representation. This file holds
# the representation problem fixed and varies WHO DESIGNED THE FEATURES. APACHE
# II sees the same measurements over the same 24-hour window and reduces each to
# a worst value mapped through a published point table. We shrink a proportion,
# fit a smooth, and take a log-likelihood ratio. Two answers to one question,
# and the difference between them is attributable to the construction rather
# than to the data, the window, or the learner.
#
# THE POINT TABLES LIVE HERE AND NOWHERE ELSE. The extraction (sql/*/v2_06_*)
# stops at raw worst-in-24h values precisely so that one scoring function serves
# both sites. A point table duplicated in two SQL files is a point table that
# will diverge, and hard rule 5 says no code branches on site. `apache2_score()`
# takes a data frame and does not know or care which site produced it.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). NO FITTING except the cross-fitted
# recalibration at the end, which is explicitly a scoring-scale conversion and
# is documented as such.
#
# AGGREGATES ONLY (hard rule 1). Scores are row-level and never printed.
#
# SOURCE. Knaus WA, Draper EA, Wagner DP, Zimmerman JE. APACHE II: a severity of
# disease classification system. Crit Care Med 1985;13(10):818-829. The bands
# below are transcribed from that paper's Table 1 and are not tuned, adjusted,
# or re-fitted. If a band here is wrong the arm is wrong; check_apache_tables()
# asserts the structural properties that a transcription error would break.
# ----------------------------------------------------------------------------

# --- the band scorer --------------------------------------------------------

#' Map a value to points through an ascending band table.
#'
#' `breaks` are the LOWER edges of bands 2..K, so `points` has one more element
#' than `breaks`. A value below every break gets `points[1]`; a value at or
#' above the last break gets `points[K]`. This is `findInterval()`'s exact
#' semantics, which is why it is used rather than `cut()` — `cut()` would need
#' an infinite outer edge and a decision about closed ends on every call, and
#' getting that wrong on one variable is invisible.
#'
#' NA in, NA out. The caller decides what a missing variable is worth; see
#' `apache2_score()`, where the decision is stated rather than defaulted.
.ap2_band <- function(x, breaks, points) {
  stopifnot(length(points) == length(breaks) + 1L, !is.unsorted(breaks))
  out <- points[findInterval(x, breaks) + 1L]
  out[is.na(x)] <- NA_real_
  as.numeric(out)
}

#' The twelve APACHE II physiologic bands, plus age.
#'
#' Each entry is (breaks, points) for `.ap2_band()`. Eleven of the twelve are
#' two-sided and U-shaped in the value, which is why the extraction emits a min
#' and a max for each and `apache2_score()` scores both ends and keeps the
#' worse. GCS and the oxygenation branch are the two exceptions and are handled
#' separately below.
AP2_BANDS <- list(
  # temperature, degrees C, rectal in the original. We use the charted value.
  temp       = list(breaks = c(30, 32, 34, 36, 38.5, 39, 41),
                    points = c(4, 3, 2, 1, 0, 1, 3, 4)),
  # mean arterial pressure, mmHg
  mbp        = list(breaks = c(50, 70, 110, 130, 160),
                    points = c(4, 2, 0, 2, 3, 4)),
  # heart rate, beats/min
  hr         = list(breaks = c(40, 55, 70, 110, 140, 180),
                    points = c(4, 3, 2, 0, 2, 3, 4)),
  # respiratory rate, breaths/min
  rr         = list(breaks = c(6, 10, 12, 25, 35, 50),
                    points = c(4, 2, 1, 0, 1, 3, 4)),
  # arterial pH
  ph         = list(breaks = c(7.15, 7.25, 7.33, 7.5, 7.6, 7.7),
                    points = c(4, 3, 2, 0, 1, 3, 4)),
  # serum bicarbonate, mmol/L. APACHE II's OWN stated substitute for pH when no
  # arterial gas is available — not an invention of ours. Used only where pH is
  # missing, and the substitution rate is reported.
  hco3       = list(breaks = c(15, 18, 22, 32, 41, 52),
                    points = c(4, 3, 2, 0, 1, 3, 4)),
  # serum sodium, mmol/L
  sodium     = list(breaks = c(111, 120, 130, 150, 155, 160, 180),
                    points = c(4, 3, 2, 0, 1, 2, 3, 4)),
  # serum potassium, mmol/L
  potassium  = list(breaks = c(2.5, 3.0, 3.5, 5.5, 6.0, 7.0),
                    points = c(4, 2, 1, 0, 1, 3, 4)),
  # serum creatinine, mg/dL. Doubled under acute renal failure; see
  # `arf_doubling` in apache2_score().
  creatinine = list(breaks = c(0.6, 1.5, 2.0, 3.5),
                    points = c(2, 0, 2, 3, 4)),
  # haematocrit, percent
  hematocrit = list(breaks = c(20, 30, 46, 50, 60),
                    points = c(4, 2, 0, 1, 2, 4)),
  # white cell count, x10^3/mm^3
  wbc        = list(breaks = c(1, 3, 15, 20, 40),
                    points = c(4, 2, 0, 1, 2, 4)),
  # A-aDO2, mmHg. Scored only when FiO2 >= 0.5.
  aado2      = list(breaks = c(200, 350, 500),
                    points = c(0, 2, 3, 4)),
  # PaO2, mmHg. Scored only when FiO2 < 0.5.
  pao2       = list(breaks = c(55, 61, 71),
                    points = c(4, 3, 1, 0)),
  # age, years. NOT part of the APS — added on top for the total only.
  age        = list(breaks = c(45, 55, 65, 75),
                    points = c(0, 2, 3, 5, 6))
)

#' Structural assertions on the band tables.
#'
#' Not a test of whether the numbers match Knaus 1985 — nothing in R can check
#' that. It checks the properties a transcription slip breaks: point vectors one
#' longer than their breaks, ascending breaks, points inside [0, 6], and a
#' zero-point band present in every two-sided variable. A table with no zero
#' band would score every patient as abnormal and the arm would look absurdly
#' strong for a reason that has nothing to do with APACHE.
#'
#' Called by `apache2_score()` on every invocation. It costs microseconds.
check_apache_tables <- function() {
  for (nm in names(AP2_BANDS)) {
    b <- AP2_BANDS[[nm]]
    if (length(b$points) != length(b$breaks) + 1L) {
      stop("AP2_BANDS[", nm, "]: points must be one longer than breaks", call. = FALSE)
    }
    if (is.unsorted(b$breaks)) {
      stop("AP2_BANDS[", nm, "]: breaks must be ascending", call. = FALSE)
    }
    if (any(b$points < 0) || any(b$points > 6)) {
      stop("AP2_BANDS[", nm, "]: points outside [0, 6]", call. = FALSE)
    }
    if (!nm %in% c("aado2", "age") && !any(b$points == 0)) {
      stop("AP2_BANDS[", nm, "]: no zero-point band; every value would score",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# --- the score --------------------------------------------------------------

#' Alveolar-arterial oxygen gradient, at sea level on room-air barometrics.
#'
#'   A-aDO2 = FiO2 x (P_atm - P_H2O) - PaCO2 / R - PaO2
#'          = FiO2 x 713 - PaCO2 / 0.8 - PaO2
#'
#' `fio2` arrives as a PERCENT at both sites (MIMIC's `bg.fio2` and eICU's
#' `apacheapsvar.fio2` are both 21-100) and is converted to a fraction here,
#' once, for both. Doing the conversion in either SQL file would be the exact
#' kind of site-local arithmetic hard rule 5 exists to prevent.
.ap2_aado2 <- function(fio2_pct, paco2, pao2) {
  f <- fio2_pct / 100
  pmax(f * 713 - paco2 / 0.8 - pao2, 0)
}

#' APACHE II, from the extracted worst-in-24h values.
#'
#' THE WORST-END RULE. Eleven variables are two-sided, so the extraction hands
#' us a min and a max and this function scores both and keeps whichever earns
#' more points. That is APACHE II's own definition of "worst value in the first
#' 24 hours": worst means most points, not most extreme.
#'
#' MISSING VARIABLES SCORE ZERO, and this is the single most important caveat in
#' the arm. It is the standard retrospective convention and it is what every
#' published APACHE-in-EHR study does, but it means a stay with no arterial gas
#' scores 0 on oxygenation exactly as a stay with a normal gas does. The
#' direction of the bias is knowable and it runs AGAINST the baseline: thin
#' coverage makes APACHE II look weaker than it is, which flatters our arm. Two
#' things follow, and both are done rather than promised. `n_vars_present`
#' travels with every row so the arm can be restricted to well-covered stays as
#' a sensitivity analysis, and `tests/metrics_severity.R` runs that restriction by
#' default rather than on request.
#'
#' @param ap   a data frame with the `ap2_*` columns of CONTRACT$severity
#' @param age  optional numeric, aligned to `ap`; needed for `total` only
#' @param arf_doubling apply APACHE II's acute-renal-failure doubling of the
#'   creatinine points. FALSE by default; see the config comment for why.
#' @param gcs_source which GCS column feeds the score. MEASURED 2026-08-30,
#'   audit A4: at MIMIC the native and reconstructed totals differ on 41.4% of
#'   stays with a mean gap of 3.80 points, while eICU -- which masks nothing --
#'   differs by 0.12. So this is the verbal convention, not hourly binning, and
#'   3.8 points is 3.8 APACHE II points on two fifths of the cohort. It is the
#'   LARGEST SINGLE LEVER in the arm, ahead of `arf_doubling` and the coverage
#'   floors, and the direction is knowable: a higher GCS means fewer points,
#'   a weaker baseline, and therefore a result that flatters our method.
#'     "native"     each site's canonical total. PRIMARY, because it is the
#'                  published convention and a baseline should be computed the
#'                  way the score normally is. Also the weakest-baseline
#'                  choice, which is why the other two must be reported.
#'     "ours"       our hourly components with verbal imputed 1 -- the
#'                  pessimistic bound, asserting an unassessable patient is
#'                  unresponsive.
#'     "ours_vnorm" the same with verbal imputed 5 -- the optimistic bound, and
#'                  the standard prospective convention.
#'   The two `ours` variants BRACKET the baseline and the paper reports the
#'   bracket. Declared in `config/apache`, never switched after seeing a
#'   result.
#' @return data frame, one row per input row: one `pt_*` column per variable,
#'   plus `aps` (the 12-variable acute physiology score), `age_points`,
#'   `chronic_points`, `total`, and `n_vars_present`.
apache2_score <- function(ap, age = NULL, arf_doubling = FALSE,
                          gcs_source = c("native", "ours", "ours_vnorm")) {
  check_apache_tables()
  gcs_source <- match.arg(gcs_source)
  gcs_col <- switch(gcs_source,
                    native     = "ap2_gcs_min_native",
                    ours       = "ap2_gcs_min_ours",
                    ours_vnorm = "ap2_gcs_min_ours_vnorm")

  need <- c("ap2_temp_min", "ap2_temp_max", "ap2_mbp_min", "ap2_mbp_max",
            "ap2_hr_min", "ap2_hr_max", "ap2_rr_min", "ap2_rr_max",
            "ap2_pao2", "ap2_paco2", "ap2_fio2", "ap2_ph_min", "ap2_ph_max",
            "ap2_hco3_min", "ap2_hco3_max", "ap2_sodium_min", "ap2_sodium_max",
            "ap2_potassium_min", "ap2_potassium_max",
            "ap2_creatinine_min", "ap2_creatinine_max",
            "ap2_hematocrit_min", "ap2_hematocrit_max",
            "ap2_wbc_min", "ap2_wbc_max", gcs_col,
            if (isTRUE(arf_doubling)) c("ap2_urine_ml_24h", "ap2_dialysis"))
  miss <- setdiff(need, names(ap))
  if (length(miss)) abort_values("apache2_score: missing input columns", miss)

  # Worst of the two ends, per two-sided variable.
  worst <- function(var, lo, hi) {
    b <- AP2_BANDS[[var]]
    a <- .ap2_band(ap[[lo]], b$breaks, b$points)
    z <- .ap2_band(ap[[hi]], b$breaks, b$points)
    pmax(a, z, na.rm = TRUE)   # NA only where BOTH ends are NA
  }

  p <- data.frame(
    pt_temp       = worst("temp",       "ap2_temp_min",       "ap2_temp_max"),
    pt_mbp        = worst("mbp",        "ap2_mbp_min",        "ap2_mbp_max"),
    pt_hr         = worst("hr",         "ap2_hr_min",         "ap2_hr_max"),
    pt_rr         = worst("rr",         "ap2_rr_min",         "ap2_rr_max"),
    pt_sodium     = worst("sodium",     "ap2_sodium_min",     "ap2_sodium_max"),
    pt_potassium  = worst("potassium",  "ap2_potassium_min",  "ap2_potassium_max"),
    pt_hematocrit = worst("hematocrit", "ap2_hematocrit_min", "ap2_hematocrit_max"),
    pt_wbc        = worst("wbc",        "ap2_wbc_min",        "ap2_wbc_max"),
    stringsAsFactors = FALSE)

  # Creatinine, with the optional acute-renal-failure doubling.
  pt_creat <- worst("creatinine", "ap2_creatinine_min", "ap2_creatinine_max")
  if (isTRUE(arf_doubling)) {
    # APACHE II doubles the creatinine points for ACUTE renal failure. The
    # operational definition used here is pre-specified: creatinine at or above
    # 1.5 mg/dL, urine output under 410 mL in the window, and no renal
    # replacement therapy (which would indicate chronic dialysis). Every one of
    # those three is imperfect in retrospective data — the urine total is
    # summed over covered hours only and therefore under-states, and `rrt` in
    # the first 24 hours does not distinguish chronic from acute — which is why
    # the flag defaults to FALSE and why turning it on is a config decision
    # made before a run rather than after seeing one.
    arf <- !is.na(ap$ap2_creatinine_max) & ap$ap2_creatinine_max >= 1.5 &
           !is.na(ap$ap2_urine_ml_24h)  & ap$ap2_urine_ml_24h < 410 &
           (is.na(ap$ap2_dialysis) | ap$ap2_dialysis == 0L)
    pt_creat[arf] <- pt_creat[arf] * 2
  }
  p$pt_creatinine <- pt_creat

  # Acid-base: pH where available, bicarbonate where it is not. APACHE II's own
  # substitution rule, applied in APACHE II's own direction.
  pt_ph   <- worst("ph",   "ap2_ph_min",   "ap2_ph_max")
  pt_hco3 <- worst("hco3", "ap2_hco3_min", "ap2_hco3_max")
  p$pt_acidbase <- ifelse(is.na(pt_ph), pt_hco3, pt_ph)

  # Oxygenation: the branch, not a value. FiO2 >= 0.5 scores the gradient,
  # below that scores PaO2 directly. A stay with a PaO2 but no FiO2 is scored on
  # the PaO2 branch, which is the conservative reading — the gradient branch
  # awards more points at the same PaO2.
  f  <- ap$ap2_fio2
  hi <- !is.na(f) & f >= 50
  aa <- .ap2_aado2(f, ap$ap2_paco2, ap$ap2_pao2)
  p$pt_oxygenation <- ifelse(
    hi,
    .ap2_band(aa, AP2_BANDS$aado2$breaks, AP2_BANDS$aado2$points),
    .ap2_band(ap$ap2_pao2, AP2_BANDS$pao2$breaks, AP2_BANDS$pao2$points))

  # GCS: points are 15 minus the score, with no bands at all. The one variable
  # where APACHE II does not discretise, and the one our design splits into
  # three separately modelled components.
  p$pt_gcs <- 15 - as.numeric(ap[[gcs_col]])

  pt_cols <- grep("^pt_", names(p), value = TRUE)
  P <- as.matrix(p[, pt_cols, drop = FALSE])

  # How many of the twelve contributed a real value, BEFORE the zero-fill. This
  # is the column that lets a reader tell "normal" from "not measured", and it
  # is why the zero-fill below is safe to state rather than dangerous to hide.
  p$n_vars_present <- rowSums(!is.na(P))
  P[is.na(P)] <- 0
  p$aps <- rowSums(P)
  p[pt_cols] <- P

  # Age and chronic health. NOT part of the APS: APACHE II's total adds them on
  # top, which is exactly the decomposition docs/v2_state_20260828.md section
  # 5.3 wants, because the LLR design deliberately excludes both.
  p$age_points <- if (is.null(age)) NA_real_ else
    .ap2_band(as.numeric(age), AP2_BANDS$age$breaks, AP2_BANDS$age$points)

  chronic <- rep(0, nrow(ap))
  if (all(c("chronic_immunocompromised", "chronic_severe_organ", "admission_class")
          %in% names(ap))) {
    has <- (!is.na(ap$chronic_immunocompromised) & ap$chronic_immunocompromised == 1L) |
           (!is.na(ap$chronic_severe_organ)      & ap$chronic_severe_organ == 1L)
    # 5 points for a non-operative or emergency post-operative admission,
    # 2 for an elective post-operative one.
    chronic <- ifelse(has,
                      ifelse(ap$admission_class == "elective_postop", 2, 5),
                      0)
  }
  p$chronic_points <- chronic
  p$total <- p$aps + ifelse(is.na(p$age_points), 0, p$age_points) + p$chronic_points
  p
}

# --- putting a points score on the log-odds scale ---------------------------

#' Cross-fitted logistic recalibration of an unscaled score.
#'
#' WHY THIS IS NEEDED AND WHAT IT DOES NOT DO. `score_report()` in R/09 requires
#' a LOG-ODDS score for its calibration column to mean anything: the slope-1
#' reference is the Bayes statement `logit p = logit(p_bar) + score`. An APACHE
#' II point total is not on that scale and never was — it is a count of points.
#' Feeding it in raw would produce a "calibration slope" of about 0.1 that says
#' nothing about APACHE and everything about the units.
#'
#' So the points are mapped to log-odds by a one-covariate logistic regression,
#' CROSS-FITTED on the same five folds the L's use. Fitting it in-sample would
#' hand the baseline an advantage the proposed method does not get, which is the
#' opposite of the asymmetry a comparison arm should carry.
#'
#' WHAT IT DOES TO DISCRIMINATION, stated precisely because the obvious claim is
#' slightly wrong. The transformation is monotone WITHIN a fold, so it cannot
#' reorder any two stays in the same fold. It is not one global transformation:
#' each fold gets its own intercept and slope, so two stays in DIFFERENT folds
#' can swap. The effect is small — on a 4,000-row synthetic check with five
#' folds it moved AUROC by 0.0015 — but it is not zero, and asserting it is zero
#' would be the kind of plausible-sounding claim this project is built to avoid.
#' `tests/metrics_severity.R` therefore MEASURES the raw-versus-recalibrated AUROC
#' gap and prints it, rather than assuming it away. If that gap is ever material
#' relative to the effect being reported, the recalibration is doing work and
#' the raw score is the number to quote.
#'
#' @param score  the raw score, any scale
#' @param y      0/1 outcome
#' @param fold_of fold index per row; the same vector the L's were built with
#' @param p_bar  the prior the returned score is centred on
#' @return numeric, on the log-odds scale, centred so that `logit(p_bar) + out`
#'   is the model's predicted log-odds
recalibrate_oof <- function(score, y, fold_of, p_bar) {
  stopifnot(length(score) == length(y), length(score) == length(fold_of))
  out <- rep(NA_real_, length(score))
  for (f in sort(unique(fold_of))) {
    te <- fold_of == f
    tr <- !te & !is.na(score) & !is.na(y)
    if (!any(tr) || !any(te)) next
    fit <- stats::glm(y[tr] ~ score[tr], family = stats::binomial())
    cf  <- stats::coef(fit)
    # ACCEPTANCE, NOT JUST A RECORD (statistical review S10, 2026-09-09). A
    # recalibration that did not converge, or that returned a non-finite
    # coefficient, must not map a score onto the log-odds scale and enter the
    # training reference as if it had.
    if (!isTRUE(fit$converged) || !all(is.finite(cf))) {
      stop(sprintf(paste0("recalibrate_oof: the fold %s recalibration %s. A ",
                          "severity cell is scored only through an accepted fit."),
                   f, if (!isTRUE(fit$converged)) "did not converge"
                      else "returned a non-finite coefficient"), call. = FALSE)
    }
    out[te] <- cf[1] + cf[2] * score[te]
  }
  out - logit(p_bar)
}

# --- the arm ----------------------------------------------------------------

#' Restrict to stays whose APACHE II is actually an APACHE II.
#'
#' A stay scored on four of twelve variables is not a severity score, it is a
#' fragment, and averaging fragments into the baseline degrades the baseline in
#' the direction that flatters the proposed method. `min_vars` is declared in
#' config before a run, never chosen after seeing which cut helps.
#'
#' Returns a logical mask rather than a subset, so the caller keeps the row
#' alignment that every other arm depends on.
#' The explicit presence check is not defensive clutter. Reading a missing
#' column returns NULL, `is.na(NULL)` is `logical(0)`, and a zero-length mask
#' propagates silently: `sum(!keep)` reports 0 dropped, `y[keep]` is empty, and
#' `all(is.na(numeric(0)))` is TRUE, so every downstream cell reports itself as
#' entirely NA and the run dies far away on an unrelated `rbind`. See the
#' companion note in `sofa_complete_mask()`, which is where that actually
#' happened.
apache_complete_mask <- function(ap, min_vars = 10L) {
  if (is.null(ap[["ap2_n_vars_present"]])) {
    stop("apache_complete_mask: `ap2_n_vars_present` is absent. A missing ",
         "column here returns a ZERO-LENGTH mask that silently empties every ",
         "downstream cell.", call. = FALSE)
  }
  !is.na(ap$ap2_n_vars_present) & ap$ap2_n_vars_present >= min_vars
}

#' Coverage of the APACHE inputs, one row per variable. Aggregates only.
#'
#' The table that has to sit beside any AUROC gap. If the baseline is thinner at
#' one site than the other, part of any cross-site difference is the baseline's
#' coverage and not the method's transportability.
apache_coverage <- function(ap) {
  vars <- c(temp = "ap2_temp_min", mbp = "ap2_mbp_min", hr = "ap2_hr_min",
            rr = "ap2_rr_min", oxygenation = "ap2_pao2", ph = "ap2_ph_min",
            hco3 = "ap2_hco3_min", sodium = "ap2_sodium_min",
            potassium = "ap2_potassium_min", creatinine = "ap2_creatinine_min",
            hematocrit = "ap2_hematocrit_min", wbc = "ap2_wbc_min",
            gcs_native = "ap2_gcs_min_native", gcs_ours = "ap2_gcs_min_ours")
  data.frame(
    variable = names(vars),
    n        = nrow(ap),
    n_present = vapply(vars, function(v)
      if (v %in% names(ap)) sum(!is.na(ap[[v]])) else NA_integer_, integer(1)),
    frac_present = round(vapply(vars, function(v)
      if (v %in% names(ap)) mean(!is.na(ap[[v]])) else NA_real_, numeric(1)), 4),
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Does the recomputed APACHE II agree with the site's native APS?
#'
#' A SANITY CHECK, not a result. APS III and APACHE IVa are different scores
#' from APACHE II, so the two must correlate strongly without agreeing. A
#' Spearman correlation below about 0.7 means the point table, a unit
#' assumption, or a sentinel guard is wrong, and the arm must not be reported
#' until it is found. Audit A5 in both SQL files points here.
apache_native_agreement <- function(ap2_total, aps_native) {
  ok <- !is.na(ap2_total) & !is.na(aps_native)
  data.frame(
    n         = sum(ok),
    spearman  = round(stats::cor(ap2_total[ok], aps_native[ok], method = "spearman"), 4),
    pearson   = round(stats::cor(ap2_total[ok], aps_native[ok], method = "pearson"), 4),
    stringsAsFactors = FALSE)
}
