# R/02_validate.R ------------------------------------------------------------
# The nine checks of docs/canonical_variable_spec.md §8, run on every table load
# at both sites.
#
# All checks run and all failures are collected before anything stops. Failing
# on the first would mean nine load-fix-reload cycles to find nine problems.
#
# Check 9 (schema equality across sites) is the one that catches the failures
# that matter: a missing column throws loudly, but a silently different factor
# level produces a plausible-looking transportability result that is actually a
# plumbing artifact, and nothing downstream flags it.
# ----------------------------------------------------------------------------

# --- the column contract ----------------------------------------------------
# Exact column sets the loader must emit, per spec §5.1-5.4 as revised
# 2026-08-22. `site` is included because the loader stamps it.

CONTRACT <- list(
  cohort = c(
    "site", "subject_id", "hadm_id", "stay_id", "intime", "outtime",
    "landmark_time", "admittime", "dischtime", "deathtime", "mortality",
    "discharge_location", "gender", "age", "race", "insurance",
    "first_careunit", "last_careunit", "anchor_year_group", "icu_los_days",
    "weight_kg", "qc_death_no_timestamp", "qc_discharge_hospice",
    "qc_weight_missing", "event_after_24h", "duration_hours_from_24h",
    "time_to_death_from_icu_hours"
  ),
  signal_features = c(
    "site", "stay_id", "signal", "signal_class", "ref_low", "ref_high",
    "n_obs", "dx_n_masked", "k_low", "k_mid", "k_high", "q05", "q95",
    # `value_min` / `value_max` lost the `dx_` prefix in SQL rev 3 (2026-08-26):
    # they are modelled for the signals whose measurement density cannot support
    # a quantile, and `dx_` is refused by the formula builder's assertion.
    "value_min", "value_max", "value_median", "dx_n_hours_early",
    "dx_n_hours_late", "dx_first_hour_low", "dx_first_hour_high",
    "dx_first_hour_obs", "dx_last_hour_obs", "dx_trend_span", "trend",
    "dx_trend_defined", "dx_value_missing",
    # Replicate structure for the `delta` within-stay variance component,
    # added 2026-09-02. Sufficient statistics, not covariates: `^se_` is
    # refused by the formula builder. `se_within_ss` and `se_within_n` pool
    # exactly across stays to give s_e = SUM(ss) / SUM(n - 1), which replaces
    # inferring s_e from heteroscedasticity in n_obs -- an inference measured
    # to be unidentified on 19 of 38 delta rows.
    "se_within_n", "se_within_mean", "se_within_ss", "se_hour_range_mean",
    "se_n_raw"
  ),
  intervention_features = c(
    "site", "stay_id", "intervention", "shape", "ever_active", "exposure_frac",
    "n_hours", "total_amount", "peak_intensity", "dx_nee_peak", "n_agents",
    "max_concurrent_agents", "has_norepinephrine", "first_hour",
    "present_at_admission"
  ),
  # `o_flag` is still IN THE PARQUET and so belongs in this contract: CONTRACT
  # describes the file schema, not the model. It was dropped from the design on
  # 2026-08-25 and is read but never coerced, never modelled, and never enters a
  # frame. Removing it here would make the validator report the extraction as
  # wrong when it is the design that changed.
  ordering = c("site", "stay_id", "signal", "intervention", "excursion_side", "o_flag"),

  # The severity-score comparators (sql/*/v2_06_severity_*.sql, scored by
  # R/09c_apache.R and R/09d_sofa.R). ONE table carrying TWO scores, because
  # they share seven inputs and splitting them would duplicate every one.
  # OPTIONAL, like `ordering`: it is a comparison arm and no model frame reads
  # it, so a site that has not run the extraction yet loads and fits normally.
  # It is in the contract anyway, because the failure this catches is the one
  # that matters — a column silently absent at one site would produce a
  # baseline that is weaker at eICU for a plumbing reason and would read as a
  # transportability result.
  #
  # NOTHING HERE MAY EVER ENTER A FORMULA. These are the comparator's inputs,
  # not ours. The `ap2_` prefix is not in R/05's refusal guard because the
  # table never reaches a model frame in the first place; if that ever changes,
  # add it to the guard rather than relying on this comment.
  severity = c(
    "site", "stay_id",
    # the site's own published acute physiology score. Different scores at the
    # two sites — APS III at MIMIC, APACHE IVa at eICU — which is why
    # `aps_native_version` travels with the value and why the two are never
    # differenced or pooled.
    "aps_native", "aps_native_version", "aps_native_prob", "severity_total_native",
    # the recomputed APACHE II inputs. Worst-in-first-24h, min and max, scored
    # by one function at both sites.
    "ap2_temp_min", "ap2_temp_max", "ap2_mbp_min", "ap2_mbp_max",
    "ap2_hr_min", "ap2_hr_max", "ap2_rr_min", "ap2_rr_max",
    "ap2_pao2", "ap2_paco2", "ap2_fio2", "ap2_pf_min",
    "ap2_ph_min", "ap2_ph_max",
    "ap2_hco3_min", "ap2_hco3_max", "ap2_sodium_min", "ap2_sodium_max",
    "ap2_potassium_min", "ap2_potassium_max",
    "ap2_creatinine_min", "ap2_creatinine_max",
    "ap2_hematocrit_min", "ap2_hematocrit_max",
    "ap2_wbc_min", "ap2_wbc_max",
    "ap2_gcs_min_native", "ap2_gcs_min_ours", "ap2_gcs_min_ours_vnorm",
    # acute renal failure inputs, for the optional creatinine doubling
    "ap2_urine_ml_24h", "ap2_dialysis",
    # secondary and site-approximate; see the SQL headers
    "chronic_immunocompromised", "chronic_severe_organ", "admission_class",
    # provenance and coverage
    "ap2_gap_source", "ap2_n_vars_present",
    # SOFA's native score and its six organ subscores. Structurally NULL at
    # eICU, which ships no SOFA of any kind; it is the RECOMPUTED SOFA that
    # transports, validated against these at MIMIC first. The subscores are the
    # point - the total alone would make SOFA a second APACHE arm, where the six
    # organs are what compare against the six domains of config/domains.csv
    # that carry a `sofa_organ`.
    "sofa_native", "sofa_native_respiration", "sofa_native_coagulation",
    "sofa_native_liver", "sofa_native_cardiovascular", "sofa_native_cns",
    "sofa_native_renal",
    # the SOFA inputs the APACHE block does not already carry. Cardiovascular
    # reuses ap2_mbp_min, CNS reuses ap2_gcs_min_native, and renal reuses
    # ap2_creatinine_max and ap2_urine_ml_24h. Sharing those is what guarantees
    # the two arms are computed from identical numbers, so any difference
    # between them is construction and never input.
    #
    # RESPIRATION IS THE EXCEPTION and has its own column, `ap2_pf_min`, above.
    # It cannot share APACHE's oxygenation triple: APACHE II scores PaO2 and
    # SOFA scores the P/F RATIO, and the gas that minimises one is not the gas
    # that minimises the other. Sharing them is what broke the SOFA respiration
    # component on 2026-08-31 (33.8% agreement, Spearman 0.086).
    "sofa_platelet_min", "sofa_bilirubin_max",
    "sofa_vent_invasive", "sofa_vent_noninvasive", "sofa_vasopressor",
    "sofa_nee_peak", "sofa_inotrope", "sofa_n_organs_present"
  )
)

# --- driver -----------------------------------------------------------------

#' Run all nine checks. Returns a report data frame; stops if anything FAILed.
#'
#' @param tabs   list from load_tables()
#' @param cfg    from load_config()
#' @param strict if FALSE, report but do not stop. For triage only — the
#'               pipeline always runs strict.
validate_tables <- function(tabs, cfg, strict = TRUE) {
  r <- .reporter()

  .check_columns(tabs, r)
  .check_value_sets(tabs, cfg, r)
  .check_class_shape(tabs, cfg, r)
  .check_pair_grid(tabs, cfg, r)
  .check_required(tabs, cfg, r)
  .check_counts(tabs, cfg, r)
  .check_ranges(tabs, cfg, r)
  .check_forbidden_guard(cfg, r)
  .check_schema_equality(tabs, r)

  rep <- r$result()
  n_fail <- sum(rep$status == "FAIL")
  n_skip <- sum(rep$status == "SKIP")
  message(sprintf("validator: %d checks, %d passed, %d failed, %d skipped",
                  nrow(rep), sum(rep$status == "PASS"), n_fail, n_skip))
  if (n_fail > 0L) {
    msg <- paste0("  [", rep$check[rep$status == "FAIL"], "] ",
                  rep$detail[rep$status == "FAIL"], collapse = "\n")
    if (strict) stop("validation failed:\n", msg, call. = FALSE)
    warning("validation failed:\n", msg, call. = FALSE)
  }
  invisible(rep)
}

.reporter <- function() {
  rows <- list()
  add <- function(check, table, ok, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, table = table,
      status = if (is.na(ok)) "SKIP" else if (isTRUE(ok)) "PASS" else "FAIL",
      detail = detail, stringsAsFactors = FALSE)
    invisible(NULL)
  }
  list(add = add, result = function() do.call(rbind, rows))
}

.fmt <- function(v, n = 8L) {
  v <- unique(as.character(v))
  paste0(paste(utils::head(v, n), collapse = ", "),
         if (length(v) > n) paste0(" ... +", length(v) - n) else "")
}

# --- 1. columns -------------------------------------------------------------

.check_columns <- function(tabs, r) {
  for (nm in names(CONTRACT)) {
    if (!nm %in% names(tabs)) {
      # `ordering` is audit-only and `apache` is a comparison arm; both are
      # legitimately absent. The other three are not.
      r$add("1 columns", nm,
            if (nm %in% c("ordering", "severity")) NA else FALSE, "table absent")
      next
    }
    got <- names(tabs[[nm]]); want <- CONTRACT[[nm]]
    extra <- setdiff(got, want); missing <- setdiff(want, got)
    ok <- !length(extra) && !length(missing)
    r$add("1 columns", nm, ok, if (ok) "" else paste0(
      if (length(missing)) paste0("missing: ", .fmt(missing)) else "",
      if (length(missing) && length(extra)) "; " else "",
      if (length(extra)) paste0("unexpected: ", .fmt(extra)) else ""))
  }
}

# --- 2. value sets ----------------------------------------------------------

.check_value_sets <- function(tabs, cfg, r) {
  # signal ⊆ the 19
  for (nm in intersect(c("signal_features", "ordering"), names(tabs))) {
    v <- setdiff(unique(as.character(tabs[[nm]]$signal)), cfg$signals)
    r$add("2 value sets", paste0(nm, "$signal"), !length(v),
          if (length(v)) paste0("outside the 19: ", .fmt(v)) else "")
  }
  # intervention ⊆ the 18 extracted
  for (nm in intersect(c("intervention_features", "ordering"), names(tabs))) {
    v <- setdiff(unique(as.character(tabs[[nm]]$intervention)), cfg$interventions_extracted)
    r$add("2 value sets", paste0(nm, "$intervention"), !length(v),
          if (length(v)) paste0("outside the extracted set: ", .fmt(v)) else "")
  }
  # every modelled intervention must actually be present in the extraction
  present <- unique(as.character(tabs$intervention_features$intervention))
  v <- setdiff(cfg$interventions_modelled, present)
  r$add("2 value sets", "modelled ⊆ present", !length(v),
        if (length(v)) paste0("modelled but absent from data: ", .fmt(v)) else "")
}

# --- 3. signal_class and shape ----------------------------------------------

.check_class_shape <- function(tabs, cfg, r) {
  s <- unique(as.data.frame(tabs$signal_features)[, c("signal", "signal_class")])
  s$signal <- as.character(s$signal); s$signal_class <- as.character(s$signal_class)
  dup <- s$signal[duplicated(s$signal)]
  r$add("3 class/shape", "signal_class single-valued", !length(dup),
        if (length(dup)) paste0("varies within: ", .fmt(dup)) else "")

  spec <- unique(cfg$pairing[, c("signal", "signal_class")])
  m <- merge(s, spec, by = "signal", suffixes = c("_data", "_spec"))
  bad <- m$signal[m$signal_class_data != m$signal_class_spec]
  r$add("3 class/shape", "signal_class vs pairing.csv", !length(bad),
        if (length(bad)) paste0("disagree: ", .fmt(bad)) else "")

  iv <- unique(as.data.frame(tabs$intervention_features)[, c("intervention", "shape")])
  iv$intervention <- as.character(iv$intervention); iv$shape <- as.character(iv$shape)
  dup <- iv$intervention[duplicated(iv$intervention)]
  r$add("3 class/shape", "shape single-valued", !length(dup),
        if (length(dup)) paste0("varies within: ", .fmt(dup)) else "")

  want <- unlist(cfg$intervention_shape)
  bad <- iv$intervention[iv$shape != want[iv$intervention]]
  bad <- bad[!is.na(bad)]
  r$add("3 class/shape", "shape vs config map", !length(bad),
        if (length(bad)) paste0("disagree with intervention_shape: ", .fmt(bad)) else "")
}

# --- 4. the pair grid, cross-checked -----------------------------------------
#
# `o_flag` was dropped from the design on 2026-08-25, and with it the level
# check that used to head this section. What survives is the part that was
# always the more valuable half: `ordering` encodes the signal x intervention
# grid independently of pairing.csv, and two independent definitions of one
# thing will disagree eventually (v2_handoff.md §5 — this is the class of check
# that found the phenylephrine bug).
#
# SKIPs rather than fails when the table is absent. `ordering` is audit-only
# now, so a site that does not ship it is a configuration choice, not an error.

.check_pair_grid <- function(tabs, cfg, r) {
  if (is.null(tabs$ordering)) {
    r$add("4 pair grid", "ordering table present", NA,
          "not loaded — audit-only since o_flag was dropped")
    return(invisible(NULL))
  }

  # ordering$excursion_side duplicates pairing.csv; pairing.csv is authoritative
  o <- unique(as.data.frame(tabs$ordering)[, c("signal", "intervention", "excursion_side")])
  o[] <- lapply(o, as.character)
  p <- cfg$pairing[cfg$pairing$intervention != "", c("signal", "intervention", "excursion_side")]
  m <- merge(o, p, by = c("signal", "intervention"), suffixes = c("_data", "_spec"), all.x = TRUE)
  bad <- m[!is.na(m$excursion_side_spec) & m$excursion_side_data != m$excursion_side_spec, ]
  r$add("4 pair grid", "excursion_side vs pairing.csv", !nrow(bad),
        if (nrow(bad)) paste0("disagree: ", .fmt(paste0(bad$signal, "/", bad$intervention))) else "")

  # THE TWO DIRECTIONS ARE NOT THE SAME FAILURE, and collapsing them was wrong.
  #
  # This check used to assert an exact set match. That is right while the two
  # tables are meant to be two encodings of one grid, and it stopped being right
  # on 2026-09-01, when `inotrope`, `rrt` and `rrt` were paired with `mbp`,
  # `creatinine` and `urine_output_rate` (canonical_variable_spec.md §4.1).
  # `ordering` is AUDIT-ONLY since `o_flag` was dropped on 2026-08-25 and is a
  # snapshot of the grid the EXTRACTION was written against; the design is now
  # deliberately ahead of it, and re-extracting it would serve no model.
  #
  #   ordering has a pair pairing.csv lacks   -> FATAL. The extraction knows
  #       about a pair the design does not, which means the design table is
  #       stale or a pair was silently dropped. This is the direction that hides
  #       a real mistake, and it stays fatal.
  #
  #   pairing.csv has a pair ordering lacks   -> REPORTED, not fatal. The design
  #       moved ahead of an audit-only snapshot. It is safe ONLY because the
  #       modelled columns come from `intervention_features`, which carries every
  #       intervention in `interventions_extracted` on every stay -- check 5
  #       asserts that separately, and `.frame_intervention()` errors on a stay
  #       absent from it. Nothing reads `ordering` to build a model frame.
  #
  # eICU never had this table, so check 4 already skips there. Reporting rather
  # than failing here also makes the two sites behave the same way, which is the
  # spirit of hard rule 5.
  got  <- paste0(o$signal, "/", o$intervention)
  want <- paste0(p$signal, "/", p$intervention)
  extra_in_data <- setdiff(got, want)
  ahead_in_spec <- setdiff(want, got)

  r$add("4 pair grid", "ordering ⊆ pairing.csv", !length(extra_in_data),
        if (length(extra_in_data))
          paste0("the extraction has pair(s) the design does not: ",
                 .fmt(extra_in_data)) else "")

  r$add("4 pair grid", "pairing.csv ahead of ordering",
        if (length(ahead_in_spec)) NA else TRUE,
        if (length(ahead_in_spec))
          paste0("design ahead of the audit-only ordering snapshot (safe; ",
                 "intervention_features carries these): ", .fmt(ahead_in_spec))
        else "")
}

# --- 5. required-by-class / by-shape / when-measured ------------------------

.check_required <- function(tabs, cfg, r) {
  s <- tabs$signal_features
  measured <- s$n_obs > 0L

  # Level terms: non-NULL exactly where measured, NULL exactly where not.
  # Both directions are contract violations (spec §5.2, §5.5).
  for (v in c("q05", "q95", "value_median")) {
    n_bad_meas   <- sum(measured & is.na(s[[v]]))
    n_bad_unmeas <- sum(!measured & !is.na(s[[v]]))
    ok <- n_bad_meas == 0L && n_bad_unmeas == 0L
    r$add("5 required", paste0("signal_features$", v), ok, if (ok) "" else
      sprintf("NA on %d measured row(s); non-NA on %d unmeasured row(s)",
              n_bad_meas, n_bad_unmeas))
  }

  # dx_value_missing is the measured mask and must agree with n_obs exactly.
  n <- sum((s$n_obs == 0L) != (s$dx_value_missing == 1L))
  r$add("5 required", "dx_value_missing == (n_obs == 0)", n == 0L,
        if (n) sprintf("%d row(s) disagree", n) else "")

  # Always-required columns.
  for (v in c("n_obs", "k_low", "k_mid", "k_high", "ref_low", "ref_high", "trend",
              "dx_trend_defined")) {
    n <- sum(is.na(s[[v]]))
    r$add("5 required", paste0("signal_features$", v), n == 0L,
          if (n) sprintf("%d NA", n) else "")
  }

  # Required-by-shape.
  iv <- tabs$intervention_features
  is_state <- iv$shape == "state"; is_event <- iv$shape == "event"
  chk <- list(
    list("exposure_frac", is_state, "state"),
    list("n_hours",       is_event, "event"),
    list("total_amount",  is_event, "event")
  )
  for (c3 in chk) {
    n <- sum(c3[[2]] & is.na(iv[[c3[[1]]]]))
    r$add("5 required", paste0("intervention_features$", c3[[1]]), n == 0L,
          if (n) sprintf("%d NA on shape=%s", n, c3[[3]]) else "")
  }

  # n_agents / max_concurrent_agents: vasopressor and inotrope only (spec §5.3
  # as revised). Required there, and required-absent everywhere else.
  wa <- unlist(cfg$interventions_with_agent_counts)
  in_wa <- as.character(iv$intervention) %in% wa
  for (v in c("n_agents", "max_concurrent_agents")) {
    n_bad_in  <- sum(in_wa & is.na(iv[[v]]))
    n_bad_out <- sum(!in_wa & !is.na(iv[[v]]))
    ok <- n_bad_in == 0L && n_bad_out == 0L
    r$add("5 required", paste0("intervention_features$", v), ok, if (ok) "" else
      sprintf("%d NA where required; %d non-NA where inapplicable", n_bad_in, n_bad_out))
  }

  # first_hour: present exactly when ever_active.
  n <- sum((iv$ever_active == 1L) & is.na(iv$first_hour))
  r$add("5 required", "first_hour when ever_active", n == 0L,
        if (n) sprintf("%d active row(s) with NA first_hour", n) else "")
}

# --- 6. count identities ----------------------------------------------------

.check_counts <- function(tabs, cfg, r) {
  s <- tabs$signal_features
  n <- sum(s$n_obs != s$k_low + s$k_mid + s$k_high, na.rm = TRUE)
  r$add("6 counts", "n_obs == k_low+k_mid+k_high", n == 0L,
        if (n) sprintf("%d row(s) violate", n) else "")
  n <- sum(s$n_obs > 24L, na.rm = TRUE)
  r$add("6 counts", "n_obs <= 24", n == 0L,
        if (n) sprintf("%d row(s) exceed 24", n) else "")

  # The signal grid must be dense: every stay x every signal.
  n_stay <- length(unique(tabs$cohort$stay_id))
  n_sig  <- length(unique(as.character(s$signal)))
  ok <- nrow(s) == n_stay * n_sig
  r$add("6 counts", "signal grid dense", ok,
        if (!ok) sprintf("%d rows != %d stays x %d signals", nrow(s), n_stay, n_sig) else "")

  # KEY INTEGRITY, NOT JUST A ROW COUNT (plumbing review F11, 2026-09-08). The
  # count above cannot tell a duplicated (stay, signal) pair from a missing one
  # -- one of each leaves the count exact -- and every feature constructor
  # downstream reaches a row by `match()`, which takes the FIRST of a duplicate
  # and says nothing. So: every key exactly once, every stay in the cohort, and
  # the signal set equal to the DECLARED one rather than to whatever the table
  # happens to contain. The same three for the intervention grid, whose density
  # is reported rather than enforced because the spec does not require it.
  co_ids <- tabs$cohort$stay_id
  key <- paste(s$stay_id, as.character(s$signal))
  n <- sum(duplicated(key))
  r$add("6 counts", "signal grid (stay, signal) unique", n == 0L,
        if (n) sprintf("%d duplicated key(s)", n) else "")
  n <- sum(!s$stay_id %in% co_ids)
  r$add("6 counts", "signal grid stay_id in cohort", n == 0L,
        if (n) sprintf("%d row(s) with a stay_id outside the cohort", n) else "")
  miss <- setdiff(cfg$signals, unique(as.character(s$signal)))
  r$add("6 counts", "signal grid covers declared signals", !length(miss),
        if (length(miss)) paste0("absent: ", .fmt(miss)) else "")
  per <- tabulate(match(s$stay_id, co_ids), nbins = length(co_ids))
  n <- sum(per != length(cfg$signals))
  r$add("6 counts", "signal grid complete per stay", n == 0L,
        if (n) sprintf("%d cohort stay(s) without exactly %d signal rows", n, length(cfg$signals)) else "")

  iv <- tabs$intervention_features
  key <- paste(iv$stay_id, as.character(iv$intervention))
  n <- sum(duplicated(key))
  r$add("6 counts", "intervention grid (stay, intervention) unique", n == 0L,
        if (n) sprintf("%d duplicated key(s)", n) else "")
  n <- sum(!iv$stay_id %in% co_ids)
  r$add("6 counts", "intervention grid stay_id in cohort", n == 0L,
        if (n) sprintf("%d row(s) with a stay_id outside the cohort", n) else "")
  per <- tabulate(match(iv$stay_id, co_ids), nbins = length(co_ids))
  r$add("6 counts", "intervention rows per stay (reported, not enforced)", TRUE,
        sprintf("distinct counts: %s", paste(sort(unique(per)), collapse = ", ")))

  .check_severity_counts(tabs, r)
}

#' The severity comparators' counts. SKIPs entirely when the table is absent.
#'
#' Three things, and the third is the one worth having. One row per cohort stay
#' catches the failure mode the eICU port is most exposed to, because
#' `apachePatientResult` carries one row per APACHE version and a careless join
#' silently doubles the cohort. The `-1` sentinel scan catches the other:
#' eICU codes "not measured" as -1 throughout `apacheApsVar`, and a -1
#' potassium scores four points on APACHE II's low band, which would inflate
#' the baseline for every unmeasured stay in a way no downstream number would
#' reveal. Coverage is reported as a count and never as a pass/fail, because
#' what counts as thin is a config decision (`apache.min_vars_present`) made
#' before a run.
.check_severity_counts <- function(tabs, r) {
  if (is.null(tabs$severity)) {
    r$add("6 counts", "severity", NA, "not loaded — comparison arms, optional")
    return(invisible(NULL))
  }
  ap <- tabs$severity

  n_stay <- length(unique(tabs$cohort$stay_id))
  ok <- nrow(ap) == n_stay && !anyDuplicated(ap$stay_id)
  r$add("6 counts", "severity one row per stay", ok,
        if (!ok) sprintf("%d rows, %d duplicated stay_id, against %d cohort stays",
                         nrow(ap), sum(duplicated(ap$stay_id)), n_stay) else "")

  # No sentinel may survive. Every one of these is a physiologic quantity that
  # cannot be negative; a negative value means a -1 guard was missed in SQL.
  num <- grep("^(ap2_|sofa_|aps_native$|aps_native_prob$|severity_total_native$)",
              names(ap), value = TRUE)
  num <- num[vapply(ap[num], is.numeric, logical(1))]
  neg <- num[vapply(num, function(v) any(ap[[v]] < 0, na.rm = TRUE), logical(1))]
  r$add("6 counts", "severity no negative sentinels", !length(neg),
        if (length(neg)) paste0("negative values in: ", .fmt(neg)) else "")

  n_thin <- sum(ap$ap2_n_vars_present < 10L, na.rm = TRUE)
  r$add("6 counts", "apache coverage (reported, not enforced)", TRUE,
        sprintf("mean %.2f of 12 variables present; %d stay(s) below 10; %d with no native score",
                mean(ap$ap2_n_vars_present, na.rm = TRUE), n_thin,
                sum(is.na(ap$aps_native))))

  # SOFA organ subscores are bounded 0-4 and the total 0-24. A value outside
  # those means the derived concept's column mapping is wrong - most likely a
  # release renamed the subscores and the query silently picked up a different
  # column. This is SQL audit A7, re-asserted where it cannot be skipped.
  org <- grep("^sofa_native_", names(ap), value = TRUE)
  bad <- org[vapply(org, function(v) any(ap[[v]] > 4L, na.rm = TRUE), logical(1))]
  if (any(ap$sofa_native > 24L, na.rm = TRUE)) bad <- c(bad, "sofa_native")
  r$add("6 counts", "sofa subscores in range", !length(bad),
        if (length(bad)) paste0("out of range: ", .fmt(bad)) else "")

  r$add("6 counts", "sofa coverage (reported, not enforced)", TRUE,
        sprintf("mean %.2f of 6 organs present; %d stay(s) below 5; %d with no native SOFA",
                mean(ap$sofa_n_organs_present, na.rm = TRUE),
                sum(ap$sofa_n_organs_present < 5L, na.rm = TRUE),
                sum(is.na(ap$sofa_native))))
}

# --- 7. ranges --------------------------------------------------------------

.check_ranges <- function(tabs, cfg, r) {
  iv <- tabs$intervention_features
  n <- sum(!is.na(iv$exposure_frac) & (iv$exposure_frac < 0 | iv$exposure_frac > 1))
  r$add("7 ranges", "exposure_frac in [0,1]", n == 0L, if (n) sprintf("%d outside", n) else "")
  n <- sum(!is.na(iv$first_hour) & (iv$first_hour < 0L | iv$first_hour > 23L))
  r$add("7 ranges", "first_hour in [0,23]", n == 0L, if (n) sprintf("%d outside", n) else "")

  co <- tabs$cohort
  cap <- cfg$guards$age_max %||% 90
  n <- sum(co$age > cap, na.rm = TRUE)
  r$add("7 ranges", paste0("age <= ", cap), n == 0L, if (n) sprintf("%d exceed", n) else "")

  lo <- cfg$guards$weight_kg_min %||% 30; hi <- cfg$guards$weight_kg_max %||% 300
  n <- sum(!is.na(co$weight_kg) & (co$weight_kg < lo | co$weight_kg > hi))
  r$add("7 ranges", sprintf("weight_kg in [%g,%g] or NA", lo, hi), n == 0L,
        if (n) sprintf("%d outside after guard", n) else "")

  n <- sum(!co$mortality %in% c(0L, 1L))
  r$add("7 ranges", "mortality in {0,1}", n == 0L, if (n) sprintf("%d outside", n) else "")

  n <- sum(duplicated(co$stay_id))
  r$add("7 ranges", "stay_id unique in cohort", n == 0L, if (n) sprintf("%d duplicated", n) else "")
}

# --- 8. no dx_ / qc_ in any formula -----------------------------------------

.check_forbidden_guard <- function(cfg, r) {
  # The guard itself lives in R/05_formula.R and fires at formula-build time.
  # What is verifiable here is that it is loaded and actually rejects a known
  # forbidden term — an assertion about the assertion (hard rule 4).
  if (!exists("assert_no_forbidden", mode = "function")) {
    r$add("8 forbidden", "guard available", FALSE,
          "assert_no_forbidden() not found; source R/05_formula.R")
    return(invisible(NULL))
  }
  probes <- c("dx_trend_defined", "qc_discharge_hospice", "dx_value_min")
  caught <- vapply(probes, function(p) {
    inherits(try(assert_no_forbidden(p), silent = TRUE), "try-error")
  }, logical(1))
  r$add("8 forbidden", "guard rejects dx_/qc_", all(caught),
        if (!all(caught)) paste0("guard let through: ", .fmt(probes[!caught])) else "")

  allowed <- inherits(try(assert_no_forbidden(c("q05", "value_median")), silent = TRUE), "try-error")
  r$add("8 forbidden", "guard admits legal terms", !allowed,
        if (allowed) "guard rejected a legal term" else "")
}

# --- 9. schema equality across sites ----------------------------------------

.check_schema_equality <- function(tabs, r) {
  sites <- unique(unlist(lapply(tabs, function(x) as.character(unique(x$site)))))
  if (length(sites) < 2L) {
    r$add("9 schema equality", paste(sites, collapse = "+"), NA,
          "only one site loaded; run validate_schema_equality() when eICU exists")
    return(invisible(NULL))
  }
  r$add("9 schema equality", paste(sites, collapse = "+"), NA,
        "multi-site frame: use validate_schema_equality(tabs_a, tabs_b)")
}

#' The post-loader schema of a table list, compact enough to travel in a bundle.
#'
#' EXTERNAL RUNNER REVIEW E7 (2026-09-09). `validate_schema_equality()` took
#' two fully loaded table lists, so the external runner materialised the whole
#' MIMIC cohort beside eICU to compare column names, classes and factor levels,
#' and it compared against whatever MIMIC extraction was on disk at apply time
#' rather than against the tables the models were fitted on. This signature is
#' exactly the information that comparison used -- names, first class, factor
#' levels with `site` excluded -- taken AFTER the loader's coercions, so it
#' describes the frame the formulas saw. The internal graph freezes it into the
#' bundle and the apply runners compare against that.
#'
#' AGGREGATES ONLY (hard rule 1): column names, classes and level sets. No
#' value of any row enters the signature.
schema_signature <- function(tabs) {
  stopifnot(is.list(tabs), !is.null(names(tabs)))
  lapply(tabs, function(x) {
    facs <- setdiff(names(x)[vapply(x, is.factor, logical(1))], "site")
    list(names   = names(x),
         classes = vapply(x, function(z) class(z)[1], character(1)),
         levels  = lapply(stats::setNames(facs, facs), function(v) levels(x[[v]])))
  })
}

#' Compare two schema signatures: names, classes and factor levels per table.
#'
#' The comparison `validate_schema_equality()` has always made, lifted onto
#' signatures so one side can be a frozen bundle slot. Tables present on only
#' one side are reported as a SKIP row rather than silently ignored: eICU has
#' no `ordering` table by design, and a reader should see that stated.
#'
#' @param label_a,label_b names for the two sides, used in messages only
compare_schema_signatures <- function(sig_a, sig_b, strict = TRUE,
                                      label_a = "a", label_b = "b") {
  r <- .reporter()
  only_a <- setdiff(names(sig_a), names(sig_b))
  only_b <- setdiff(names(sig_b), names(sig_a))
  if (length(only_a) || length(only_b)) {
    r$add("9 tables", "all", NA, paste0(
      if (length(only_a)) paste0("only in ", label_a, ": ", .fmt(only_a)) else "",
      if (length(only_a) && length(only_b)) "; " else "",
      if (length(only_b)) paste0("only in ", label_b, ": ", .fmt(only_b)) else ""))
  }
  common_tabs <- intersect(names(sig_a), names(sig_b))
  if (!length(common_tabs)) {
    r$add("9 tables", "all", FALSE, "no table in common")
  }
  for (nm in common_tabs) {
    a <- sig_a[[nm]]; b <- sig_b[[nm]]

    d <- c(setdiff(a$names, b$names), setdiff(b$names, a$names))
    r$add("9 names", nm, !length(d), if (length(d)) .fmt(d) else "")

    common <- intersect(a$names, b$names)
    ta <- a$classes[common]; tb <- b$classes[common]
    bad <- common[ta != tb]
    r$add("9 types", nm, !length(bad), if (length(bad)) paste0(
      .fmt(paste0(bad, " (", ta[bad], " vs ", tb[bad], ")"))) else "")

    # `site` is never in a signature's level set; it is the one column that
    # must differ between sites.
    facs <- intersect(intersect(names(a$levels), names(b$levels)), common)
    bad <- facs[!vapply(facs, function(v) identical(a$levels[[v]], b$levels[[v]]), logical(1))]
    r$add("9 levels", nm, !length(bad), if (length(bad)) .fmt(bad) else "")
  }
  rep <- r$result()
  n_fail <- sum(rep$status == "FAIL")
  message(sprintf("schema equality (%s vs %s): %d checks, %d failed",
                  label_a, label_b, nrow(rep), n_fail))
  if (n_fail && strict) {
    stop("schema equality failed:\n",
         paste0("  [", rep$check[rep$status == "FAIL"], " ", rep$table[rep$status == "FAIL"],
                "] ", rep$detail[rep$status == "FAIL"], collapse = "\n"), call. = FALSE)
  }
  invisible(rep)
}

#' Check 9 proper: two sites' tables must match in column names, types, and
#' factor levels. Call with the MIMIC and eICU table lists. Implemented as the
#' signature comparison, so the live and the frozen path cannot disagree.
validate_schema_equality <- function(tabs_a, tabs_b, strict = TRUE) {
  la <- unique(unlist(lapply(tabs_a, function(x) as.character(unique(x$site)))))
  lb <- unique(unlist(lapply(tabs_b, function(x) as.character(unique(x$site)))))
  compare_schema_signatures(schema_signature(tabs_a), schema_signature(tabs_b),
                            strict = strict,
                            label_a = if (length(la) == 1L) la else "a",
                            label_b = if (length(lb) == 1L) lb else "b")
}
