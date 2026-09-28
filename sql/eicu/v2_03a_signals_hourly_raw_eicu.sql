-- =====================================================================
-- v2 STEP 3a: raw long-format signal values, eICU (materialised once)
--
-- PORT OF sql/mimiciv/v2_03a_signals_hourly_raw.sql. Same reason for
-- existing: BigQuery bills BYTES SCANNED, and re-running step 3b with
-- different plausibility bounds must not re-pay the scan.
--
-- The expensive reads at eICU are `vitalperiodic` (45.4M rows) and
-- `nursecharting` (temperature 6.3M rows, GCS 8.2M rows). MIMIC's single
-- expensive block was `chartevents` for point-of-care glucose; here it is
-- two tables, and `nursecharting` is read by BOTH this file and step 04.
-- So this file materialises a NURSECHARTING SLICE first, and step 04
-- reads that slice rather than re-scanning. One scan, two consumers.
-- That is the same argument that produced 3a in the first place.
--
-- NO BOUNDS, NO BINNING, NO MASKING. Facts only.
--
-- ---------------------------------------------------------------------
-- SITE SUBSTITUTIONS, per signal
-- ---------------------------------------------------------------------
--   heart_rate   vitalperiodic.heartrate
--   resp_rate    vitalperiodic.respiration
--   spo2         vitalperiodic.sao2
--   mbp          vitalaperiodic.noninvasivemean   <- SEE DECISION B BELOW
--   temperature  nursecharting 'Temperature (C)'  <- SEE DECISION C BELOW
--                 UNION vitalperiodic.temperature
--   labs         `lab`, EXACT labname whitelist   <- SEE DECISION I BELOW
--   glucose      `lab` 'glucose' + 'bedside glucose'
--   urine_output_rate  NOT HERE - it needs weight from the cohort, so it
--                 stays in 3b exactly as at MIMIC.
--
-- ---------------------------------------------------------------------
-- DECISION B: NON-INVASIVE PREFERRED, ARTERIAL AS PER-HOUR FALLBACK.
-- REVISED 2026-08-28 - the first version of this file said non-invasive
-- ONLY, on an argument that measurement showed was overstated. Both the
-- old argument and the correction are kept below, because the correction
-- is the reason the rule is affordable.
--
-- THE RULE NOW: for each stay-hour, take the cuff readings if the hour has
-- any; take the arterial readings only for hours with no cuff reading at
-- all. `source` distinguishes them ('vitalaperiodic' vs
-- 'vitalperiodic_arterial') so the mixture is measurable. MIMIC's 03a
-- applies the identical rule as COALESCE(mbp_ni, mbp).
--
-- NOBODY LOSES MBP. The 2,230 arterial-only stays keep the signal instead
-- of falling to n_obs = 0 -> L = 0, and those are the sickest quarter of
-- the cohort - exactly the patients the score exists to identify.
--
-- WHY NOT SIMPLY POOL, and why that objection was overstated.
-- ---------------------------------------------------------------------
-- Stay-level availability in the LOS-filtered eICU cohort (audit 18):
--   non-invasive only  80,009 stays,  7.3% mortality
--   both               26,605 stays, 12.7% mortality
--   invasive only       2,230 stays, 10.6% mortality
--   neither             1,684 stays,  7.4% mortality
--
-- Arterial-line presence is a SEVERITY MARKER (12.7% vs 7.3%), and
-- invasive MAP is charted at 11.4 readings/hour against non-invasive
-- 2.5 (audit 02_04). Free pooling therefore appears to hand the sicker
-- quarter of the cohort ~5x the measurement density of everyone else -
-- a density confound CORRELATED WITH THE OUTCOME, which `n_obs` would
-- carry straight into pi_hat's shrinkage denominator.
--
-- *** THE CORRECTION: THE 5x NEVER REACHES THE MODEL. ***
-- `n_obs` counts COVERED HOURS, not raw readings. Step 3b collapses each
-- hour to `hr_med` before step 5 counts anything, so n_obs is capped at
-- 24 whatever the interface rate. What survives binning is the COVERAGE
-- differential - an arterial line fills every hour, a cuff misses some -
-- and that is worth about three hours out of twenty-four, not a factor
-- of five.
--
-- MEASURED, and this is the number that settled it: eICU non-invasive-
-- only `mbp` gets 20.50 mean covered hours (audit 01 of this run) against
-- MIMIC's pooled 21.73 (mimiciv audit 03_01). The two sites' `mbp`
-- coverage already agrees to within 1.2 hours, which is closer than
-- `heart_rate` (23.40 vs 22.22). The density argument was real but small,
-- and it does not justify deleting a signal from the sickest 2% of the
-- cohort.
--
-- THE PER-HOUR RESTRICTION IS WHAT KEEPS IT SMALL. Arterial readings are
-- admitted only for hours with NO cuff reading, so an a-line never adds
-- density to an hour the cuff already covers. Free pooling would have
-- re-introduced the 5x inside the hour; this does not.
--
-- REMAINING COST, to state rather than to fix: arterial and oscillometric
-- MAP are not the same instrument, so `mbp` is an instrument mixture
-- whose composition correlates with severity. It is now the SAME mixture
-- rule at both sites, which is what the transport comparison needs -
-- MIMIC's `mbp` was already such a mixture and nobody had said so.
--
-- *** COMPANION CHANGE AT MIMIC: MADE 2026-08-28. ***
-- `mimiciv_derived.vitalsign.mbp` pools arterial and cuff with arterial
-- winning. `v2_03a_signals_hourly_raw.sql` now reads
-- COALESCE(mbp_ni, mbp), which is the same preference rule as above.
-- Audit F in that file measures MIMIC's arterial-only share and carries a
-- pre-specified rule for switching both sites to strict non-invasive if
-- that share turns out to be small. RUN IT BEFORE THE RE-EXTRACTION.
--
-- ---------------------------------------------------------------------
-- DECISION C: TEMPERATURE IS RECOVERED FROM nursecharting
-- ---------------------------------------------------------------------
-- `vitalperiodic.temperature` covers only 0.112 of stays stay-weighted
-- (audit 16), 164 of 166 hospitals below 0.5 - which would have cost the
-- unpaired control group one of its seven members, and the unpaired
-- group carries the paper's most important internal contrast.
-- Audit 02_03a/03b found it: `nursecharting` has Temperature at
-- 6,267,541 rows over 187,371 stays, recovering 105,222 of the cohort
-- (~96%). Both sources are unioned and `source` distinguishes them, the
-- same way MIMIC pools chemistry / bg / poc glucose.
--
-- UNITS: nursecharting carries BOTH 'Temperature (C)' and
-- 'Temperature (F)' as separate valnames for the same reading. Only the
-- (C) valname is taken. The defensive Fahrenheit conversion still runs in
-- 3b, unchanged from MIMIC, and is a no-op if this holds.
--
-- ---------------------------------------------------------------------
-- DECISION I: EXACT labname WHITELIST, never LIKE
-- ---------------------------------------------------------------------
-- Audit 04 found exactly the contamination the MIMIC glucose audit
-- predicted: '%creatinine%' pulls urinary creatinine (p50 82.75 mg/dL),
-- '%glucose%' pulls glucose-CSF, '%wbc%' pulls seven body-fluid variants
-- with p95 up to 270,300, '%bilirubin%' pulls direct bilirubin,
-- '%urea%' pulls 24h urine urea nitrogen. A single one of those would
-- destroy a signal's reference-range counts silently.
--
-- UNITS ALIGN WITH MIMIC and need no rescaling. Despite the name,
-- 'platelets x 1000' reads K/mcL with p50 = 197 and 'WBC x 1000' reads
-- K/mcL with p50 = 10.07 - already MIMIC's units. The name is the trap,
-- not the data.
--
-- 'bedside glucose' is the point-of-care analogue: 3.18M rows and 5.5
-- covered hours per stay against serum glucose's 1.6 (audit 20b). The
-- MIMIC decision to add fingerstick glucose transfers cleanly and
-- matters just as much here.
-- =====================================================================

DECLARE epoch TIMESTAMP DEFAULT TIMESTAMP '2014-01-01 00:00:00';


-- =====================================================================
-- 3a-0. NURSECHARTING SLICE. Cohort x window x the four labels this
-- project needs. Read once here; step 04 reads THIS, not nursecharting.
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_nursecharting_slice_eicu`
CLUSTER BY item, stay_id AS
SELECT
  n.patientunitstayid AS stay_id,
  CASE
    WHEN n.nursingchartcelltypevallabel = 'Temperature'
     AND n.nursingchartcelltypevalname  = 'Temperature (C)'  THEN 'temperature'
    WHEN n.nursingchartcelltypevallabel = 'Glasgow coma score'
     AND n.nursingchartcelltypevalname  = 'Motor'            THEN 'gcs_motor'
    WHEN n.nursingchartcelltypevallabel = 'Glasgow coma score'
     AND n.nursingchartcelltypevalname  = 'Eyes'             THEN 'gcs_eyes'
    WHEN n.nursingchartcelltypevallabel = 'Glasgow coma score'
     AND n.nursingchartcelltypevalname  = 'Verbal'           THEN 'gcs_verbal'
  END                                                        AS item,
  n.nursingchartoffset                                       AS offset_min,
  -- RAW string kept alongside the parsed number. eICU charts non-numeric
  -- verbal values ('1T', 'Unable to score') that SAFE_CAST drops; those
  -- rows are the closest thing eICU has to MIMIC's `gcs_unable` flag and
  -- step 04 uses the raw string to count them. Do not drop this column.
  n.nursingchartvalue                                        AS value_raw,
  SAFE_CAST(n.nursingchartvalue AS FLOAT64)                  AS value
FROM `physionet-data.eicu_crd.nursecharting` n
JOIN `eicu-ext.eicu_ext_data.v2_cohort_eicu` c
  ON n.patientunitstayid = c.stay_id
WHERE n.nursingchartoffset BETWEEN 0 AND 1439
  AND ((n.nursingchartcelltypevallabel = 'Temperature'
        AND n.nursingchartcelltypevalname = 'Temperature (C)')
    OR (n.nursingchartcelltypevallabel = 'Glasgow coma score'
        AND n.nursingchartcelltypevalname IN ('Motor', 'Eyes', 'Verbal')));


-- =====================================================================
-- 3a-1. RAW VALUES. Same schema as v2_raw_values_mimiciv:
--       stay_id, signal, charttime, value, source
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_raw_values_eicu`
CLUSTER BY signal, stay_id AS

WITH cohort AS (
  SELECT stay_id
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

-- ---------------------------------------------------------------------
-- VITALS from the 5-minute interface (vitalperiodic)
-- ---------------------------------------------------------------------
vp AS (
  SELECT v.patientunitstayid AS stay_id, v.observationoffset AS offset_min,
         v.heartrate, v.respiration, v.sao2, v.temperature
  FROM `physionet-data.eicu_crd.vitalperiodic` v
  JOIN cohort c ON v.patientunitstayid = c.stay_id
  WHERE v.observationoffset BETWEEN 0 AND 1439
),

-- ---------------------------------------------------------------------
-- LABS. Exact whitelist, one scan.
-- ---------------------------------------------------------------------
lb AS (
  SELECT l.patientunitstayid AS stay_id, l.labresultoffset AS offset_min,
         l.labname, l.labresult
  FROM `physionet-data.eicu_crd.lab` l
  JOIN cohort c ON l.patientunitstayid = c.stay_id
  WHERE l.labresultoffset BETWEEN 0 AND 1439
    AND l.labresult IS NOT NULL
    AND l.labname IN ('sodium', 'bicarbonate', 'creatinine', 'BUN',
                      'glucose', 'bedside glucose', 'Hgb',
                      'platelets x 1000', 'WBC x 1000',
                      'total bilirubin', 'lactate')
)

SELECT stay_id, 'heart_rate' AS signal,
       TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE) AS charttime,
       CAST(heartrate AS FLOAT64) AS value, 'vitalperiodic' AS source
FROM vp WHERE heartrate IS NOT NULL

UNION ALL
SELECT stay_id, 'resp_rate', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(respiration AS FLOAT64), 'vitalperiodic'
FROM vp WHERE respiration IS NOT NULL

UNION ALL
SELECT stay_id, 'spo2', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(sao2 AS FLOAT64), 'vitalperiodic'
FROM vp WHERE sao2 IS NOT NULL

UNION ALL
-- DECISION B, REVISED 2026-08-28: NON-INVASIVE PREFERRED, NOT NON-INVASIVE
-- ONLY. The cuff channel is taken wherever it exists; the arterial channel
-- fills only the hours it does not. `source` records which, so the mixture
-- is measurable rather than assumed. MIMIC's 03a now applies the identical
-- rule via COALESCE(mbp_ni, mbp). See the header.
SELECT a.patientunitstayid, 'mbp',
       TIMESTAMP_ADD(epoch, INTERVAL a.observationoffset MINUTE),
       CAST(a.noninvasivemean AS FLOAT64), 'vitalaperiodic'
FROM `physionet-data.eicu_crd.vitalaperiodic` a
JOIN cohort c ON a.patientunitstayid = c.stay_id
WHERE a.observationoffset BETWEEN 0 AND 1439
  AND a.noninvasivemean IS NOT NULL

UNION ALL
-- Arterial fallback. Restricted to the STAY-HOURS with no cuff reading, so
-- an arterial line never adds density to an hour the cuff already covers -
-- which is what keeps the severity-correlated coverage differential to the
-- ~3 hours measured in the header rather than reintroducing the 5x.
-- DELETE THIS BLOCK for strict non-invasive; nothing else changes.
SELECT v.patientunitstayid, 'mbp',
       TIMESTAMP_ADD(epoch, INTERVAL v.observationoffset MINUTE),
       CAST(v.systemicmean AS FLOAT64), 'vitalperiodic_arterial'
FROM `physionet-data.eicu_crd.vitalperiodic` v
JOIN cohort c ON v.patientunitstayid = c.stay_id
LEFT JOIN (
  SELECT DISTINCT patientunitstayid,
         CAST(FLOOR(observationoffset / 60) AS INT64) AS hour_bin
  FROM `physionet-data.eicu_crd.vitalaperiodic`
  WHERE observationoffset BETWEEN 0 AND 1439 AND noninvasivemean IS NOT NULL
) ni
  ON v.patientunitstayid = ni.patientunitstayid
 AND CAST(FLOOR(v.observationoffset / 60) AS INT64) = ni.hour_bin
WHERE v.observationoffset BETWEEN 0 AND 1439
  AND v.systemicmean IS NOT NULL
  AND ni.patientunitstayid IS NULL

UNION ALL
-- DECISION C. Temperature, primary source.
SELECT s.stay_id, 'temperature',
       TIMESTAMP_ADD(epoch, INTERVAL s.offset_min MINUTE),
       s.value, 'nursecharting'
FROM `eicu-ext.eicu_ext_data.v2_nursecharting_slice_eicu` s
WHERE s.item = 'temperature' AND s.value IS NOT NULL

UNION ALL
-- Temperature, secondary source. Adds the ~12k stays charted on the
-- periodic interface but not by nursing.
SELECT stay_id, 'temperature', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(temperature AS FLOAT64), 'vitalperiodic'
FROM vp WHERE temperature IS NOT NULL

-- ---------------------------------------------------------------------
-- LABS -> the eleven sparse signals plus glucose
-- ---------------------------------------------------------------------
UNION ALL
SELECT stay_id, 'sodium', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'sodium'

UNION ALL
SELECT stay_id, 'bicarbonate', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'bicarbonate'

UNION ALL
SELECT stay_id, 'creatinine', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'creatinine'

UNION ALL
SELECT stay_id, 'bun', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'BUN'

UNION ALL
SELECT stay_id, 'wbc', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'WBC x 1000'

UNION ALL
SELECT stay_id, 'hemoglobin', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'Hgb'

UNION ALL
SELECT stay_id, 'platelet', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'platelets x 1000'

UNION ALL
SELECT stay_id, 'lactate', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'lactate'

UNION ALL
SELECT stay_id, 'bilirubin_total', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'lab'
FROM lb WHERE labname = 'total bilirubin'

UNION ALL
-- Serum glucose. `source` mirrors MIMIC's 'chemistry'.
SELECT stay_id, 'glucose', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'chemistry'
FROM lb WHERE labname = 'glucose'

UNION ALL
-- Point-of-care glucose. `source` mirrors MIMIC's 'poc'.
SELECT stay_id, 'glucose', TIMESTAMP_ADD(epoch, INTERVAL offset_min MINUTE),
       CAST(labresult AS FLOAT64), 'poc'
FROM lb WHERE labname = 'bedside glucose';


-- =====================================================================
-- AUDITS
-- =====================================================================
-- A. LABNAME VERIFICATION - run on `lab` alone first, it is the cheap
--    scan. Confirm the eleven strings above are the only ones wanted and
--    that no new variant has appeared.
-- SELECT labname, unit, COUNT(*) FROM `physionet-data.eicu_crd.lab`
-- WHERE LOWER(labname) LIKE '%glucose%' OR LOWER(labname) LIKE '%creatinine%'
--    OR LOWER(labname) LIKE '%wbc%'     OR LOWER(labname) LIKE '%bilirubin%'
-- GROUP BY 1,2 ORDER BY 3 DESC;
--
-- B. TEMPERATURE UNITS - is 'Temperature (C)' really Celsius?
--    Audit 02_03b saw ~990k rows above 45 under the LOOSE label match;
--    under the exact (C) valname this should be near zero. If it is not,
--    the 3b conversion catches it, which is why it stays.
-- SELECT COUNTIF(value BETWEEN 25 AND 43) AS n_celsius,
--        COUNTIF(value BETWEEN 90 AND 110) AS n_fahrenheit, COUNT(*) AS n
-- FROM `...v2_raw_values_eicu` WHERE signal = 'temperature';
--
-- C. GLUCOSE DENSITY BY SOURCE. 'bedside glucose' should dominate, as
--    fingerstick does at MIMIC.
-- SELECT source, COUNT(*) AS n_rows, COUNT(DISTINCT stay_id) AS n_stays,
--        COUNT(*) / COUNT(DISTINCT stay_id) AS rows_per_stay
-- FROM `...v2_raw_values_eicu` WHERE signal = 'glucose' GROUP BY source;
--
-- D. THE DENSITY COMPARISON THAT MATTERS. Re-run discovery Q4 against
--    THIS table and compare `mean_min_gap_per_iqr` to the MIMIC column of
--    the same name (audit 02_01). MIMIC is 0.000-0.020 across every
--    signal; eICU was 0.13-0.45 for the vitalperiodic signals and 0.093
--    for non-invasive MAP.
--    WHY THIS IS NO LONGER DESIGN-CHANGING: config/level_terms_by_class
--    gives dense and rate signals `quantile` (q05/q95 of hourly MEDIANS)
--    and only sparse signals `extreme` (value_min/value_max of raw
--    readings). The min/max path is therefore taken ONLY by the labs and
--    the three GCS components, which chart at ~1 reading/hour at BOTH
--    sites. The site-dependence reconciliation section A identified lands
--    on a statistic no dense signal uses. Re-measure it anyway and put
--    the table in the supplement - the argument is only as good as the
--    numbers behind it.
--
-- E. TEMPERATURE COVERAGE, the check that decides whether the unpaired
--    control group keeps its seventh member.
-- SELECT COUNT(DISTINCT stay_id) FROM `...v2_raw_values_eicu`
-- WHERE signal = 'temperature';   -- expect ~106,000 of 109,200
--
-- F. TABLE SIZE. eICU vitalperiodic is 45M rows before filtering; the
--    materialised output should still be low single-digit GB.
-- SELECT table_id, ROUND(size_bytes/POW(1024,3), 3) AS gb, row_count
-- FROM `eicu-ext.eicu_ext_data.__TABLES__`
-- WHERE table_id IN ('v2_raw_values_eicu', 'v2_nursecharting_slice_eicu');
--
-- G. REJECTION AUDIT, per signal, against the 3b bounds.
-- SELECT signal, COUNT(*) AS n_raw, MIN(value), MAX(value)
-- FROM `...v2_raw_values_eicu` GROUP BY signal ORDER BY signal;
