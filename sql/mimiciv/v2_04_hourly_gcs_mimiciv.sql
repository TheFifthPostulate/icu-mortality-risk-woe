-- =====================================================================
-- v2 STEP 4: MIMIC-IV GCS components -> hourly lattice
--
-- REWRITTEN to emit the SAME SCHEMA as v2_hourly_mimiciv, so step 5 can
-- simply UNION the two and apply masking + thresholds uniformly.
--
-- WHY THE EARLIER STAY-LEVEL VERSION HAD TO GO
--   It aggregated straight to n_obs / K counts / quantiles. But GCS is the
--   signal that MOST needs step-5 masking: under deep sedation or
--   neuromuscular blockade a low GCS is not a low score, it is an ABSENT
--   measurement. Masking has to happen before aggregation, because dropping
--   hours changes n_obs, K, q05 and q95. Aggregating here would have locked
--   in unmasked values.
--
-- THE FACTS / DECISIONS LINE
--   Stays here (recorded facts):
--     - gcs_unable filtering. The clinician DOCUMENTED that verbal could not
--       be assessed. That is data, not a modelling choice.
--   Moves to step 5 (decisions):
--     - impairment thresholds (motor < 5, eyes < 3, verbal < 3)
--     - sedation / paralytic masking
--   Same principle as reference ranges for the other signals.
--
-- PER-TIMEPOINT VERBAL MASKING (unchanged from the earlier version, and the
--   main fix over v1): v1 masked verbal for the ENTIRE stay if gcs_unable
--   was ever 1, discarding 16h of valid assessments from a patient
--   extubated at hour 8. Here the flag is applied row by row.
--
-- WHY COMPONENTS, NOT TOTAL GCS
--   Verbal is unscorable in intubated patients, so total GCS would be NULL
--   for exactly the sickest subgroup. Motor and eyes stay assessable.
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_gcs_mimiciv`
PARTITION BY RANGE_BUCKET(hour_bin, GENERATE_ARRAY(0, 24, 1))
CLUSTER BY stay_id, signal AS

WITH cohort AS (
  SELECT stay_id, intime, landmark_time
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

gcs_raw AS (
  SELECT
    c.stay_id,
    c.intime,
    g.charttime,
    SAFE_CAST(g.gcs_motor  AS FLOAT64) AS gcs_motor,
    SAFE_CAST(g.gcs_verbal AS FLOAT64) AS gcs_verbal,
    SAFE_CAST(g.gcs_eyes   AS FLOAT64) AS gcs_eyes,
    SAFE_CAST(g.gcs_unable AS INT64)   AS gcs_unable
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
    ON c.stay_id = g.stay_id
   AND g.charttime >= c.intime
   AND g.charttime <  c.landmark_time
),

-- Long format. Component bounds are structural (motor 1-6, eyes 1-4,
-- verbal 1-5) - these are scale definitions, not plausibility judgements,
-- so they belong here rather than in signal_spec.
component_values AS (
  SELECT stay_id, intime, charttime, 'gcs_motor' AS signal, gcs_motor AS value
  FROM gcs_raw
  WHERE gcs_motor BETWEEN 1 AND 6

  UNION ALL
  SELECT stay_id, intime, charttime, 'gcs_eyes', gcs_eyes
  FROM gcs_raw
  WHERE gcs_eyes BETWEEN 1 AND 4

  UNION ALL
  -- Verbal dropped AT THE TIMEPOINT the patient was not assessable.
  SELECT stay_id, intime, charttime, 'gcs_verbal', gcs_verbal
  FROM gcs_raw
  WHERE gcs_verbal BETWEEN 1 AND 5
    AND COALESCE(gcs_unable, 0) = 0
)

SELECT
  cv.stay_id,
  cv.signal,
  'dense'                                          AS signal_class,
  TIMESTAMP_DIFF(cv.charttime, cv.intime, HOUR)    AS hour_bin,
  COUNT(*)                                         AS n_raw_in_hour,
  MIN(cv.value)                                    AS hr_min,
  APPROX_QUANTILES(cv.value, 2)[OFFSET(1)]         AS hr_med,
  MAX(cv.value)                                    AS hr_max
FROM component_values cv
GROUP BY cv.stay_id, cv.signal, hour_bin
HAVING hour_bin BETWEEN 0 AND 23;


-- =====================================================================
-- Stay-level airway diagnostic. Does not fit the hourly schema and is not
-- a modelling input - kept for the cross-check against the ventilation
-- table in step 2. Large disagreement means one of the two airway
-- definitions is wrong.
-- =====================================================================
CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_gcs_unable_mimiciv` AS
SELECT
  c.stay_id,
  COUNT(g.charttime)                                        AS n_assessments,
  COUNTIF(g.gcs_unable = 1)                                 AS n_unable,
  SAFE_DIVIDE(COUNTIF(g.gcs_unable = 1), COUNT(g.charttime)) AS gcs_unable_frac,
  MAX(COALESCE(g.gcs_unable, 0))                            AS gcs_unable_ever
FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv` c
LEFT JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
  ON c.stay_id = g.stay_id
 AND g.charttime >= c.intime
 AND g.charttime <  c.landmark_time
GROUP BY c.stay_id;


-- =====================================================================
-- AUDITS
-- =====================================================================
-- 1. Coverage per component. Verbal n_obs should be materially lower than
--    motor/eyes, and the gap should track gcs_unable_frac. If verbal
--    matches motor exactly, the masking is not firing.
-- SELECT signal, AVG(n_hours) AS mean_hours,
--        APPROX_QUANTILES(n_hours, 4) AS quartiles,
--        COUNTIF(n_hours = 0) AS n_never
-- FROM (SELECT stay_id, signal, COUNT(*) AS n_hours
--       FROM `...v2_hourly_gcs_mimiciv` GROUP BY stay_id, signal)
-- GROUP BY signal;
--
--    NOTE: if mean_hours comes back low (say < 4), reclassify these as
--    'sparse' in step 5 - min/max only, no pi. GCS is charted with neuro
--    checks, typically q1-4h, so 'dense' is a hypothesis this audit tests.
--
-- 2. Value distribution - confirms the ordinal scales are intact and shows
--    where impairment thresholds will actually land.
-- SELECT signal, hr_med, COUNT(*) AS n
-- FROM `...v2_hourly_gcs_mimiciv` GROUP BY signal, hr_med ORDER BY signal, hr_med;
--
--    Use this to set the step-5 thresholds empirically. If "below component
--    maximum" would flag ~90% of hours, that threshold does not discriminate
--    and motor < 5 / eyes < 3 / verbal < 3 are the better cutpoints.
--
-- 3. Airway cross-check against step 2 ventilation:
-- SELECT u.gcs_unable_ever, COUNT(DISTINCT v.stay_id) AS n_vented
-- FROM `...v2_gcs_unable_mimiciv` u
-- LEFT JOIN `...v2_interventions_mimiciv` v
--   ON u.stay_id = v.stay_id AND v.intervention = 'vent_invasive'
-- GROUP BY u.gcs_unable_ever;
