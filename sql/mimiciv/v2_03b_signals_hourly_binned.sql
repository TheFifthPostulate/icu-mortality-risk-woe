-- =====================================================================
-- v2 STEP 3b: hourly signal lattice
--
-- Reads the MATERIALISED raw table from step 3a, not physionet-data, so
-- re-running with different plausibility bounds costs a few hundred MB
-- instead of ~12 GB. Bounds and signal_class are declared once, in
-- signal_spec - that table IS the pre-specification. Commit before fitting.
--
-- ONE EXCEPTION: urine output still reads physionet directly. The
-- urine_output derived table is small, and the rate needs weight from the
-- cohort anyway, so materialising it would buy nothing.
--
-- WHY BOUNDS GO ON RAW VALUES, NOT ON QUANTILES
--   hr_min is the MIN of raw readings in an hour, so a single artifact sets
--   the hourly minimum. q05 over ~24 hourly values then lands at roughly the
--   1st order statistic - i.e. the artifact. Quantiles give almost NO
--   protection once the artifact has been absorbed at the binning step.
--
--   Bounds are PHYSIOLOGICALLY IMPOSSIBLE limits, not clinically unusual
--   ones: WBC 407 (leukemic hyperleukocytosis) and platelet 2360
--   (myeloproliferative) are real and retained. Fixed a priori, identical
--   for train / test / eICU - no test-set contamination.
--
-- STILL NO REFERENCE RANGES, STILL NO RRT/SEDATION MASKING - step 5.
-- SPARSE BY DESIGN - grid completion happens at stay level in step 5.
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_mimiciv`
PARTITION BY RANGE_BUCKET(hour_bin, GENERATE_ARRAY(0, 24, 1))
CLUSTER BY stay_id, signal AS

WITH cohort AS (
  SELECT stay_id, intime, landmark_time, weight_kg
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

-- =====================================================================
-- DESIGN TABLE: plausibility bounds + variable template class
-- =====================================================================
signal_spec AS (
  SELECT 'heart_rate'      AS signal, 'dense'  AS signal_class,  20.0 AS lo,  300.0 AS hi UNION ALL
  SELECT 'mbp',                       'dense',                   20.0,        200.0      UNION ALL
  SELECT 'resp_rate',                 'dense',                    1.0,         70.0      UNION ALL
  SELECT 'spo2',                      'dense',                   20.0,        100.0      UNION ALL
  SELECT 'temperature',               'dense',                   25.0,         43.0      UNION ALL
  SELECT 'glucose',                   'dense',                   10.0,       2000.0      UNION ALL
  SELECT 'sodium',                    'sparse',                  90.0,        200.0      UNION ALL
  SELECT 'bicarbonate',               'sparse',                   2.0,         60.0      UNION ALL
  SELECT 'creatinine',                'sparse',                   0.1,         30.0      UNION ALL
  SELECT 'bun',                       'sparse',                   1.0,        250.0      UNION ALL
  SELECT 'wbc',                       'sparse',                   0.1,        500.0      UNION ALL
  SELECT 'hemoglobin',                'sparse',                   2.0,         25.0      UNION ALL
  SELECT 'platelet',                  'sparse',                   5.0,       3000.0      UNION ALL
  SELECT 'lactate',                   'sparse',                   0.1,         30.0      UNION ALL
  SELECT 'bilirubin_total',           'sparse',                   0.1,         60.0      UNION ALL
  -- Applies to the RATE (mL/kg/hr) after weight normalisation, enforced
  -- inside hourly_uo rather than on raw charted volume.
  SELECT 'urine_output_rate',         'rate',                     0.0,         20.0
),

-- ---------------------------------------------------------------------
-- Unit normalisation BEFORE bounds.
-- Defensive Fahrenheit conversion: a no-op if derived.vitalsign is
-- uniformly Celsius (audit B), and the correct fix if it is not. Doing
-- this as CONVERSION rather than exclusion means mixed-unit rows are
-- recovered instead of silently deleted.
-- ---------------------------------------------------------------------
normalised AS (
  SELECT
    stay_id,
    signal,
    charttime,
    CASE
      WHEN signal = 'temperature' AND value BETWEEN 90 AND 110
        THEN (value - 32.0) * 5.0 / 9.0
      ELSE value
    END AS value,
    source
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_raw_values_mimiciv`
),

-- Bounds applied on raw values, before any aggregation.
raw_vals AS (
  SELECT n.stay_id, n.signal, n.charttime, n.value
  FROM normalised n
  JOIN signal_spec s ON n.signal = s.signal
  WHERE n.value BETWEEN s.lo AND s.hi
),

hourly_point AS (
  SELECT
    r.stay_id,
    r.signal,
    TIMESTAMP_DIFF(r.charttime, c.intime, HOUR) AS hour_bin,
    COUNT(*)                                    AS n_raw_in_hour,
    MIN(r.value)                                AS hr_min,
    APPROX_QUANTILES(r.value, 2)[OFFSET(1)]     AS hr_med,
    MAX(r.value)                                AS hr_max
  FROM raw_vals r
  JOIN cohort c ON r.stay_id = c.stay_id
  GROUP BY r.stay_id, r.signal, hour_bin
),

-- ---------------------------------------------------------------------
-- Urine output: a RATE, not a point measurement.
--   * weight_kg constrained to 30-300. Absurd observed rates (max 1600
--     mL/kg/hr) traced to a near-zero denominator.
--   * hourly net volume floored at 0. mimic-code's urine_output subtracts
--     GU irrigant volume, so irrigated patients produce genuinely negative
--     net output (observed min -93.9). Net output cannot be negative.
--   * rates above 20 mL/kg/hr dropped - cumulative-total charting events
--     misread as hourly increments.
-- Stays with absent or implausible weight yield no UO rows -> n_obs = 0
-- downstream. That is correct: unavailable, not zero.
-- ---------------------------------------------------------------------
hourly_uo_raw AS (
  SELECT
    c.stay_id,
    'urine_output_rate'                             AS signal,
    TIMESTAMP_DIFF(u.charttime, c.intime, HOUR)     AS hour_bin,
    COUNT(*)                                        AS n_raw_in_hour,
    SAFE_DIVIDE(GREATEST(SUM(u.urineoutput), 0.0),
                c.weight_kg)                        AS uo_rate
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_derived.urine_output` u
    ON c.stay_id = u.stay_id
   AND u.charttime >= c.intime AND u.charttime < c.landmark_time
  WHERE u.urineoutput IS NOT NULL
    AND c.weight_kg BETWEEN 30 AND 300
  GROUP BY c.stay_id, hour_bin, c.weight_kg
),

hourly_uo AS (
  SELECT
    stay_id, signal, hour_bin, n_raw_in_hour,
    uo_rate AS hr_min,
    uo_rate AS hr_med,
    uo_rate AS hr_max
  FROM hourly_uo_raw
  WHERE uo_rate IS NOT NULL
    AND uo_rate BETWEEN 0.0 AND 20.0
),

unioned AS (
  SELECT * FROM hourly_point
  UNION ALL
  SELECT * FROM hourly_uo
)

SELECT
  u.stay_id,
  u.signal,
  s.signal_class,
  u.hour_bin,
  u.n_raw_in_hour,
  u.hr_min,
  u.hr_med,
  u.hr_max
FROM unioned u
JOIN signal_spec s ON u.signal = s.signal
WHERE u.hour_bin BETWEEN 0 AND 23;


-- =====================================================================
-- AUDITS
-- =====================================================================
-- 1. Hours covered per signal - the n that shrinkage acts on.
-- SELECT signal, signal_class,
--        APPROX_QUANTILES(n_hours, 4) AS quartiles, AVG(n_hours) AS mean_hours
-- FROM (SELECT stay_id, signal, ANY_VALUE(signal_class) AS signal_class,
--              COUNT(*) AS n_hours
--       FROM `...v2_hourly_mimiciv` GROUP BY stay_id, signal)
-- GROUP BY signal, signal_class ORDER BY signal_class, signal;
--
-- 2. Range sanity - every min/max must now sit inside signal_spec bounds.
--    UO specifically: min >= 0, max <= 20.
-- SELECT signal, MIN(hr_min), MAX(hr_max)
-- FROM `...v2_hourly_mimiciv` GROUP BY signal ORDER BY signal;
--
-- 3. Glucose density after pooling - should be well above the serum-only
--    2.0 hours. If not, the itemid list missed a variant.
-- SELECT APPROX_QUANTILES(n_hours, 4), AVG(n_hours) FROM (
--   SELECT stay_id, COUNT(*) AS n_hours FROM `...v2_hourly_mimiciv`
--   WHERE signal = 'glucose' GROUP BY stay_id);
--
-- 4. Rejection audit - now trivial against the 3a table.
-- SELECT r.signal, COUNT(*) AS n_raw,
--        COUNTIF(r.value < s.lo OR r.value > s.hi) AS n_rejected,
--        SAFE_DIVIDE(COUNTIF(r.value < s.lo OR r.value > s.hi), COUNT(*))
--          AS frac_rejected
-- FROM `...v2_raw_values_mimiciv` r
-- JOIN (<paste signal_spec>) s USING (signal)
-- GROUP BY r.signal ORDER BY frac_rejected DESC;
--
-- 5. Temperature conversion check - how many rows were F and got converted?
-- SELECT COUNTIF(value BETWEEN 90 AND 110) AS n_fahrenheit, COUNT(*) AS n
-- FROM `...v2_raw_values_mimiciv` WHERE signal = 'temperature';
