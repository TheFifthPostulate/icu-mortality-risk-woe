-- =====================================================================
-- v2 STEP 3b: hourly signal lattice, eICU
--
-- PORT OF sql/mimiciv/v2_03b_signals_hourly_binned.sql. `signal_spec`,
-- the normalisation rule, the bounds-before-aggregation ordering and the
-- urine-output guards are COPIED VERBATIM. This is the file where that
-- matters most: bounds are pre-specified, identical for train / test /
-- eICU, and any divergence here is a site effect manufactured in SQL.
--
-- If you change a bound, change it in BOTH files in the same commit.
--
-- ONE EXCEPTION, same as at MIMIC: urine output does not come from the 3a
-- table. It needs weight from the cohort and its source table is small,
-- so materialising it would buy nothing.
--
-- ---------------------------------------------------------------------
-- WHY BOUNDS GO ON RAW VALUES, NOT ON QUANTILES (unchanged, and it bites
-- harder here). hr_min is the MIN of raw readings in an hour, so a single
-- artifact sets the hourly minimum. At eICU's 11.5 readings/hour there
-- are ~10x more chances per hour for an artifact to land than at MIMIC's
-- 1.2, and audit 19b measured the consequence directly: the top ventile
-- of the hourly min-to-median gap for heart rate is 152 bpm. The bounds
-- in `signal_spec` are the only thing standing between that and
-- `value_min`.
--
-- ---------------------------------------------------------------------
-- URINE OUTPUT: cellpath is a RULE, not a list
-- ---------------------------------------------------------------------
-- Audit 10 shows the same failure mode as the MIMIC cumulative-total
-- problem, in a different costume: the `%urine%` sweep returns counts and
-- flags masquerading as volumes -
--   Urine Count (3,106 stays, p95 = 57)
--   Urine Occurrence (601 stays, p95 = 5)
--   urine incontinence / Urine Incontinence (32 stays, p95 = 1)
--   Mixed Urine/Stool Volume (536 stays) - not urine alone
-- The dominant path is '...|I&O|Output (ml)|Urine' (145,186 stays), but
-- there are 25+ spelling variants of Foley and nephrostomy paths behind
-- it. A 25-item literal whitelist would be brittle to a re-extraction; an
-- INCLUDE rule plus an EXCLUDE rule is auditable and generalises. Audit 3
-- prints what the rule actually matched - run it, do not assume.
-- =====================================================================

CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_hourly_eicu`
PARTITION BY RANGE_BUCKET(hour_bin, GENERATE_ARRAY(0, 24, 1))
CLUSTER BY stay_id, signal AS

WITH cohort AS (
  SELECT stay_id, intime, landmark_time, weight_kg
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

-- =====================================================================
-- DESIGN TABLE: plausibility bounds + variable template class
-- VERBATIM COPY of the MIMIC file. Do not edit one without the other.
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
  SELECT 'urine_output_rate',         'rate',                     0.0,         20.0
),

-- ---------------------------------------------------------------------
-- Unit normalisation BEFORE bounds. VERBATIM from MIMIC.
-- At eICU this is not merely defensive: audit 17a found ~1.9% of
-- temperature rows in Fahrenheit form (64,539 rows above 45 C against
-- 3.34M in the Celsius band) under the loose label match. Step 3a takes
-- only the 'Temperature (C)' valname, which should make this a no-op -
-- audit 5 says whether it did. CONVERSION rather than exclusion means
-- mixed-unit rows are recovered instead of silently deleted.
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
  FROM `eicu-ext.eicu_ext_data.v2_raw_values_eicu`
),

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
-- Urine output: a RATE, not a point measurement. Guards VERBATIM from
-- MIMIC:
--   * weight_kg constrained to 30-300 (already applied in step 1)
--   * hourly net volume floored at 0
--   * rates above 20 mL/kg/hr dropped - cumulative-total charting events
--     misread as hourly increments
-- Stays with absent or implausible weight yield no UO rows -> n_obs = 0
-- downstream. That is correct: unavailable, not zero. eICU has 8.3% such
-- stays against MIMIC's much smaller share, so the unmeasured fraction
-- for this one signal differs by site by construction. Report it.
-- ---------------------------------------------------------------------
uo_raw AS (
  SELECT
    io.patientunitstayid AS stay_id,
    io.intakeoutputoffset AS offset_min,
    io.cellvaluenumeric   AS volume_ml
  FROM `physionet-data.eicu_crd.intakeoutput` io
  JOIN cohort c ON io.patientunitstayid = c.stay_id
  WHERE io.intakeoutputoffset BETWEEN 0 AND 1439
    AND io.cellvaluenumeric IS NOT NULL
    -- INCLUDE rule
    AND LOWER(io.cellpath) LIKE '%|output (ml)|%'
    AND (LOWER(io.cellpath) LIKE '%urine%'
      OR LOWER(io.cellpath) LIKE '%foley%'
      OR LOWER(io.cellpath) LIKE '%neph%')
    -- EXCLUDE rule: counts, occurrences, incontinence flags, mixed output
    AND LOWER(io.cellpath) NOT LIKE '%count%'
    AND LOWER(io.cellpath) NOT LIKE '%occurrence%'
    AND LOWER(io.cellpath) NOT LIKE '%incontinen%'
    AND LOWER(io.cellpath) NOT LIKE '%stool%'
),

hourly_uo_raw AS (
  SELECT
    c.stay_id,
    'urine_output_rate'                             AS signal,
    CAST(FLOOR(u.offset_min / 60) AS INT64)         AS hour_bin,
    COUNT(*)                                        AS n_raw_in_hour,
    SAFE_DIVIDE(GREATEST(SUM(u.volume_ml), 0.0),
                c.weight_kg)                        AS uo_rate
  FROM cohort c
  JOIN uo_raw u ON c.stay_id = u.stay_id
  WHERE c.weight_kg BETWEEN 30 AND 300
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
-- 1. HOURS COVERED PER SIGNAL - the n that shrinkage acts on, and the
--    single most important cross-site table in this extraction. Compare
--    every row to the MIMIC figures in audit 02_02b: heart_rate 22.71,
--    mbp 22.00, spo2 22.38, resp_rate 22.54, temperature 7.67, glucose
--    6.63, and 1.49-2.85 for the sparse labs.
--    A dense signal at ~24 here against ~22 at MIMIC is fine. A sparse
--    lab materially below MIMIC's ~2.2 is not, and it changes what
--    pi_hat's confidence weighting means at the two sites.
-- SELECT signal, signal_class,
--        APPROX_QUANTILES(n_hours, 4) AS quartiles, AVG(n_hours) AS mean_hours,
--        COUNT(*) AS n_stays_with_any
-- FROM (SELECT stay_id, signal, ANY_VALUE(signal_class) AS signal_class,
--              COUNT(*) AS n_hours
--       FROM `...v2_hourly_eicu` GROUP BY stay_id, signal)
-- GROUP BY signal, signal_class ORDER BY signal_class, signal;
--
-- 2. RANGE SANITY - every min/max must sit inside signal_spec bounds.
--    UO specifically: min >= 0, max <= 20.
-- SELECT signal, MIN(hr_min), MAX(hr_max)
-- FROM `...v2_hourly_eicu` GROUP BY signal ORDER BY signal;
--
-- 3. WHAT THE URINE CELLPATH RULE ACTUALLY MATCHED. Run this and read
--    every row. It is the check that the include/exclude rule did what
--    the header says, and it is cheap.
-- SELECT cellpath, COUNT(*) AS n, COUNT(DISTINCT patientunitstayid) AS n_stays,
--        APPROX_QUANTILES(cellvaluenumeric, 20)[OFFSET(19)] AS p95
-- FROM `physionet-data.eicu_crd.intakeoutput`
-- WHERE LOWER(cellpath) LIKE '%|output (ml)|%'
--   AND (LOWER(cellpath) LIKE '%urine%' OR LOWER(cellpath) LIKE '%foley%'
--     OR LOWER(cellpath) LIKE '%neph%')
--   AND LOWER(cellpath) NOT LIKE '%count%' AND LOWER(cellpath) NOT LIKE '%occurrence%'
--   AND LOWER(cellpath) NOT LIKE '%incontinen%' AND LOWER(cellpath) NOT LIKE '%stool%'
-- GROUP BY 1 ORDER BY n_stays DESC;
--
-- 4. REJECTION AUDIT against the 3a table, per signal. Compare the
--    rejected fraction to MIMIC's. A materially higher rejection rate for
--    a dense vital at eICU is the 5-minute interface's artifact tail, and
--    it is expected; a higher rate for a LAB is a vocabulary problem.
-- SELECT r.signal, COUNT(*) AS n_raw,
--        COUNTIF(r.value < s.lo OR r.value > s.hi) AS n_rejected,
--        SAFE_DIVIDE(COUNTIF(r.value < s.lo OR r.value > s.hi), COUNT(*)) AS frac
-- FROM `...v2_raw_values_eicu` r
-- JOIN (<paste signal_spec>) s USING (signal)
-- GROUP BY r.signal ORDER BY frac DESC;
--
-- 5. TEMPERATURE CONVERSION - how many rows were F and got converted?
--    Should be near zero given step 3a's exact (C) valname.
-- SELECT COUNTIF(value BETWEEN 90 AND 110) AS n_fahrenheit, COUNT(*) AS n
-- FROM `...v2_raw_values_eicu` WHERE signal = 'temperature';
--
-- 6. URINE OUTPUT COVERAGE vs MIMIC. Audit 16 put eICU urine coverage at
--    0.698 stay-weighted with 40 hospitals below 0.5, before the weight
--    requirement is applied. The number after both filters is what the
--    model actually sees.
-- SELECT COUNT(DISTINCT stay_id) FROM `...v2_hourly_eicu`
-- WHERE signal = 'urine_output_rate';
