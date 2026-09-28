-- =====================================================================
-- v2 STEP 4: eICU GCS components -> hourly lattice
--
-- PORT OF sql/mimiciv/v2_04_hourly_gcs_mimiciv.sql. Emits the SAME SCHEMA
-- as v2_hourly_eicu, so step 5 unions the two and applies masking and
-- thresholds uniformly.
--
-- Reads `v2_nursecharting_slice_eicu` from step 3a, NOT `nursecharting`.
-- One scan of that 8M-row table serves both files; re-running this one
-- costs a few hundred MB. Same argument as 3a existing at all.
--
-- WHY COMPONENTS, NOT TOTAL GCS: unchanged. Verbal is unscorable in
-- intubated patients, so total GCS would be NULL for exactly the sickest
-- subgroup. eICU charts 'GCS Total' at 3,248,821 rows (audit 06a) and
-- the design correctly does not use it - the slice in 3a does not even
-- extract it.
--
-- THE FACTS / DECISIONS LINE, unchanged:
--   Stays here (recorded facts): component scale bounds.
--   Moves to step 5 (decisions): impairment thresholds, sedation and
--   paralytic masking.
--
-- =====================================================================
-- THE ONE REAL LOSS: THERE IS NO `gcs_unable` ANALOGUE IN eICU
-- =====================================================================
-- MIMIC's `mimiciv_derived.gcs.gcs_unable` is a clinician's DOCUMENTED
-- statement that verbal could not be assessed, and the MIMIC file drops
-- verbal at exactly those timepoints. That is data, not a modelling
-- choice, which is why it lives in step 4 rather than step 5.
--
-- eICU has nothing equivalent (reconciliation section L). The closest
-- candidates and why neither works:
--   * `apacheapsvar.meds` - complete for all 171,177 rows, but it is
--     stay-level and day-1-worst, so it cannot drive PER-TIMEPOINT
--     masking, which is the whole point of the v1 -> v2 fix.
--   * `nursecharting` Sedation Scale/Score - 1,450,629 time-stamped rows,
--     but a sedation score is a measurement, not a statement that GCS was
--     unassessable, and using it would fold a decision into the facts
--     layer.
--
-- WHAT IS ACTUALLY DONE. Nothing is substituted. Verbal is extracted
-- unfiltered, and the non-numeric rows are counted in
-- `v2_gcs_unable_eicu` below so that the size of the gap is a reported
-- number rather than a silent one.
--
-- *** CORRECTED 2026-09-04, AND THE ORIGINAL CLAIM IS LEFT BELOW SO THE
-- ERROR IS VISIBLE. *** This paragraph used to say that eICU charts
-- non-numeric verbal values, naming '1T', 'Unable to score' and 'Unable
-- to assess', and called `v2_gcs_unable_eicu` a PARTIAL ANALOGUE of
-- MIMIC's flag on that basis. Audit round 6, query U2, grouped the
-- non-numeric rows by their actual value: there is exactly ONE value per
-- component in the cohort window and it is BLANK. Verbal has 15,731 such
-- rows over 2,508 stays, motor 12,801 over 1,955, eyes 12,814 over
-- 1,958. No '1T' appears, and neither does any 'Unable' string.
--
-- So `gcs_unable_frac` counts BLANK CELLS, which is missingness, and it
-- is not an analogue of anything. U2b sharpens it further: of 14,988
-- stay-hours with a non-numeric verbal, 78.1% also have a non-numeric
-- motor AND eyes, so the construct behind a blank is that the Glasgow
-- assessment was not recorded this hour, and not that verbal could not
-- be assessed in this intubated patient.
--
-- NOTHING NUMERICAL DEPENDS ON THIS. `v2_gcs_unable_eicu` is built here
-- and read by no R file, no runner and no target, so no result changes.
-- What changes is what may be written in the methods: eICU has no
-- partial unassessability construct, not a small one.
--
-- PRE-SPECIFY THE CONSEQUENCE: eICU verbal will include timepoints a
-- MIMIC clinician would have marked unassessable, and those are
-- concentrated in ventilated patients. If verbal is scored 1 rather than
-- left blank, eICU's verbal distribution gains a spike at 1 that MIMIC's
-- flag removes. Check 4 below measures it directly against
-- `invasive_vent`. If the spike is large, the honest options are to
-- report gcs_verbal's transport result with that caveat attached, or to
-- drop gcs_verbal from the cross-site comparison and keep motor and eyes
-- - which stay assessable at both sites. DECIDE BEFORE FITTING.
--
-- TREND COVERAGE, the other site difference worth knowing in advance.
-- Audit 20a: eICU GCS coverage is mean 6.4 hours per stay and `trend` is
-- definable for 49.7% of stays. Audit 02_02a measured MIMIC: eyes 7.6
-- hours / 85.3% definable, motor 7.58 / 85.0%, verbal 6.25 / 61.1%.
-- The COVERAGE is comparable; the DEFINABLE FRACTION is not - 0.497 at
-- eICU against 0.85 at MIMIC for motor and eyes. Since `trend = 0` when
-- undefined, the point mass at exactly zero is roughly 15% at MIMIC and
-- 50% at eICU, and `s(trend)` would see a materially different input
-- distribution. This is a real transport threat and it needs a decision
-- (reconciliation section L): accept and report the definable fraction
-- per site, or restrict GCS trend. It is NOT fixed in this file.
-- =====================================================================

CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_hourly_gcs_eicu`
PARTITION BY RANGE_BUCKET(hour_bin, GENERATE_ARRAY(0, 24, 1))
CLUSTER BY stay_id, signal AS

WITH slice AS (
  SELECT stay_id, item, offset_min, value
  FROM `eicu-ext.eicu_ext_data.v2_nursecharting_slice_eicu`
  WHERE item IN ('gcs_motor', 'gcs_eyes', 'gcs_verbal')
    AND value IS NOT NULL
),

-- Long format. Component bounds are structural (motor 1-6, eyes 1-4,
-- verbal 1-5) - these are scale definitions, not plausibility judgements,
-- so they belong here rather than in signal_spec. VERBATIM from MIMIC.
--
-- The verbal branch is the ONE place the two files differ: MIMIC adds
-- `AND COALESCE(gcs_unable, 0) = 0`. eICU has no such column, and no
-- substitute is invented. See the header.
component_values AS (
  SELECT stay_id, offset_min, 'gcs_motor' AS signal, value
  FROM slice
  WHERE item = 'gcs_motor' AND value BETWEEN 1 AND 6

  UNION ALL
  SELECT stay_id, offset_min, 'gcs_eyes', value
  FROM slice
  WHERE item = 'gcs_eyes' AND value BETWEEN 1 AND 4

  UNION ALL
  SELECT stay_id, offset_min, 'gcs_verbal', value
  FROM slice
  WHERE item = 'gcs_verbal' AND value BETWEEN 1 AND 5
)

SELECT
  cv.stay_id,
  cv.signal,
  'dense'                                          AS signal_class,
  CAST(FLOOR(cv.offset_min / 60) AS INT64)         AS hour_bin,
  COUNT(*)                                         AS n_raw_in_hour,
  MIN(cv.value)                                    AS hr_min,
  APPROX_QUANTILES(cv.value, 2)[OFFSET(1)]         AS hr_med,
  MAX(cv.value)                                    AS hr_max
FROM component_values cv
GROUP BY cv.stay_id, cv.signal, hour_bin
HAVING hour_bin BETWEEN 0 AND 23;


-- =====================================================================
-- Stay-level airway diagnostic. Does not fit the hourly schema and is not
-- a modelling input - kept for the cross-check against invasive_vent from
-- step 2b, and as the closest available comparator to MIMIC's
-- gcs_unable_frac.
--
-- READ THE COLUMN NAMES AS DEFINED HERE, NOT AS AT MIMIC. `n_unable` is
-- NOT a documented-unassessable count: it is the count of verbal rows
-- whose charted value is non-numeric, which is a DIFFERENT and almost
-- certainly SMALLER construct. Comparing this fraction to MIMIC's
-- gcs_unable_frac is informative about the size of the gap; it is not a
-- like-for-like substitution and must not be presented as one.
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_gcs_unable_eicu` AS
SELECT
  c.stay_id,
  COUNTIF(s.item = 'gcs_verbal')                                    AS n_assessments,
  COUNTIF(s.item = 'gcs_verbal' AND s.value IS NULL)                AS n_unable,
  SAFE_DIVIDE(COUNTIF(s.item = 'gcs_verbal' AND s.value IS NULL),
              NULLIF(COUNTIF(s.item = 'gcs_verbal'), 0))            AS gcs_unable_frac,
  MAX(IF(s.item = 'gcs_verbal' AND s.value IS NULL, 1, 0))          AS gcs_unable_ever
FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu` c
LEFT JOIN `eicu-ext.eicu_ext_data.v2_nursecharting_slice_eicu` s
  ON c.stay_id = s.stay_id
GROUP BY c.stay_id;


-- =====================================================================
-- AUDITS
-- =====================================================================
-- 1. COVERAGE PER COMPONENT, against MIMIC (audit 02_02a): eyes 7.60
--    hours / 85.3% trend-definable, motor 7.58 / 85.0%, verbal 6.25 /
--    61.1%. Audit 20a predicts ~6.4 hours and ~49.7% here.
--    THE DEFINABLE FRACTION IS THE NUMBER THAT MATTERS - it sets the
--    height of the point mass at trend = 0, and that mass is a mixture of
--    "flat" and "not estimable" whose composition differs by site.
-- SELECT signal, AVG(n_hours) AS mean_hours,
--        APPROX_QUANTILES(n_hours, 4) AS quartiles,
--        AVG(CAST(n_early >= 1 AND n_late >= 1
--                 AND (max_h - min_h) >= 6 AS INT64)) AS frac_trend_definable
-- FROM (SELECT stay_id, signal, COUNT(*) AS n_hours,
--              COUNTIF(hour_bin < 12) AS n_early, COUNTIF(hour_bin >= 12) AS n_late,
--              MIN(hour_bin) AS min_h, MAX(hour_bin) AS max_h
--       FROM `...v2_hourly_gcs_eicu` GROUP BY stay_id, signal)
-- GROUP BY signal;
--
--    NOTE, as at MIMIC: if mean_hours comes back below ~4, reclassify
--    these as 'sparse' in step 5 - min/max only, no pi. 'dense' is a
--    hypothesis this audit tests, and it is a hypothesis about the SITE,
--    so it must be tested separately here.
--
-- 2. VERBAL vs MOTOR COVERAGE. At MIMIC verbal is materially lower than
--    motor because gcs_unable fires. At eICU, with no such filter, the
--    two should be NEARLY EQUAL. If verbal is much lower here anyway,
--    something else is dropping it and that needs explaining before the
--    signal is used.
-- SELECT signal, COUNT(DISTINCT stay_id) AS n_stays, COUNT(*) AS n_stay_hours
-- FROM `...v2_hourly_gcs_eicu` GROUP BY 1;
--
-- 3. VALUE DISTRIBUTION - confirms the ordinal scales are intact and
--    shows where the step-5 impairment thresholds land.
-- SELECT signal, hr_med, COUNT(*) AS n
-- FROM `...v2_hourly_gcs_eicu` GROUP BY 1, 2 ORDER BY 1, 2;
--
-- 4. *** THE UNASSESSABLE-CONTAMINATION CHECK. Run this one. ***
--    Verbal = 1 among ventilated stays, eICU vs MIMIC. A spike here that
--    MIMIC does not have is intubated patients scored 1 instead of
--    flagged unassessable, and it inflates gcs_verbal's apparent
--    impairment in exactly the sickest subgroup - which is where the
--    transport comparison is most load-bearing.
-- WITH vent AS (SELECT stay_id FROM `...v2_intervention_features_eicu`
--               WHERE intervention = 'invasive_vent' AND ever_active = 1)
-- SELECT v.stay_id IS NOT NULL AS vented, g.hr_med, COUNT(*) AS n
-- FROM `...v2_hourly_gcs_eicu` g LEFT JOIN vent v USING (stay_id)
-- WHERE g.signal = 'gcs_verbal' GROUP BY 1, 2 ORDER BY 1, 2;
--
-- 5. NON-NUMERIC VERBAL RATE - the partial gcs_unable analogue. Compare
--    to MIMIC's gcs_unable_frac, and read the caveat in the side-table
--    header before quoting the comparison.
-- SELECT AVG(gcs_unable_frac), AVG(gcs_unable_ever)
-- FROM `...v2_gcs_unable_eicu` WHERE n_assessments > 0;
--
--    And the vocabulary itself, which is worth seeing once:
-- SELECT value_raw, COUNT(*) FROM `...v2_nursecharting_slice_eicu`
-- WHERE item = 'gcs_verbal' AND value IS NULL GROUP BY 1 ORDER BY 2 DESC;
