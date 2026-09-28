-- =====================================================================
-- v2 STEP 3a: raw long-format signal values (materialised once)
--
-- WHY THIS EXISTS
--   BigQuery bills BYTES SCANNED IN REFERENCED COLUMNS, not rows returned.
--   The POC glucose block reads mimiciv_3_1_icu.chartevents (several hundred
--   million rows); the LIKE filter removes rows from the RESULT but not from
--   the SCAN. That single block took step 3 from ~1 GB to ~12 GB, and every
--   bounds tweak re-paid it.
--
--   This materialises the expensive read ONCE. Step 3b then reads this table
--   instead of physionet-data, so re-running with different plausibility
--   bounds costs a few hundred MB rather than 12 GB.
--
-- ON OUTPUT SIZE - this does NOT write 12 GB.
--   12 GB is the amount SCANNED. What lands here is roughly:
--     ~5 dense vitals  x ~25 rows/stay
--     ~10 sparse labs  x ~2-3 rows/stay
--     POC glucose      x ~5-8 rows/stay
--   For a ~50k cohort that is on the order of 8-10M rows at ~36 bytes each
--   uncompressed -> a few hundred MB stored, less after columnar
--   compression. At BigQuery active-storage rates that is cents per month.
--   Set an expiration below if you want it to clean itself up.
--
-- TWO SCAN OPTIMISATIONS APPLIED
--   1. LITERAL itemids instead of a JOIN to d_items. A join forces a full
--      scan; a literal IN list lets BigQuery prune if chartevents is
--      clustered on itemid. VERIFY THE ITEMIDS FIRST (audit A below) -
--      a wrong or missing itemid truncates coverage with no error.
--   2. bg.glucose added. mimiciv_derived.bg is ALREADY being scanned for
--      lactate, so blood-gas glucose costs nothing extra and is drawn far
--      more often than serum chemistry in ICU patients. Check audit C: if
--      bg + chemistry alone gets glucose to adequate density, you may not
--      need the chartevents read at all.
--
-- NO BOUNDS, NO BINNING, NO MASKING. Facts only.
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_raw_values_mimiciv`
CLUSTER BY signal, stay_id
-- OPTIONS (expiration_timestamp = TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 90 DAY))
AS

WITH cohort AS (
  SELECT stay_id, hadm_id, intime, landmark_time
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
)

-- ---------------------------------------------------------------------
-- VITALS  (mimiciv_derived.vitalsign, stay-keyed)
-- ---------------------------------------------------------------------
SELECT c.stay_id, 'heart_rate' AS signal, v.charttime,
       CAST(v.heart_rate AS FLOAT64) AS value, 'vitalsign' AS source
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
  ON c.stay_id = v.stay_id
 AND v.charttime >= c.intime AND v.charttime < c.landmark_time
WHERE v.heart_rate IS NOT NULL

UNION ALL
-- ---------------------------------------------------------------------
-- MBP SOURCE RULE, changed 2026-08-28. NON-INVASIVE PREFERRED PER READING.
--
-- WAS: `v.mbp`, which mimic-code already coalesces over the arterial and
-- the oscillometric channels with arterial winning. That made `mbp` a
-- source mixture whose composition correlates with severity, and it did
-- not match the eICU extraction, where the two channels live in separate
-- tables (`vitalperiodic.systemicmean` vs `vitalaperiodic.noninvasivemean`)
-- and must be chosen between explicitly.
--
-- IS: COALESCE(mbp_ni, mbp) — the cuff reading when the hour has one, the
-- arterial reading only when it does not. eICU's 03a does the identical
-- COALESCE over its two tables, so the SOURCE PREFERENCE is now the same
-- rule at both sites, which is the whole point.
--
-- WHY NOT STRICT NON-INVASIVE (i.e. `v.mbp_ni` alone). It was the first
-- proposal, and it is one token away — delete the COALESCE. It was
-- rejected because it deletes the signal for arterial-only stays
-- (n_obs = 0 -> L = 0), and arterial-only is not rare in the sickest
-- subgroup: when an a-line is in, cuff readings are often not charted at
-- all. At eICU that group is 2,230 stays (2.0%, audit 18). At MIMIC the
-- number is UNMEASURED — run audit F below before treating this as
-- settled. If MIMIC's arterial-only share is also ~2%, strict
-- non-invasive is defensible and marginally cleaner; if it is 10%+, it
-- zeroes out `mbp` for exactly the patients the score exists to identify.
--
-- WHY POOLING WAS NOT SIMPLY KEPT. The objection in
-- sql/eicu/eicu_extraction_decisions.md §1 was that arterial lines chart
-- at ~11.4 readings/hour against the cuff's ~2.5 (eICU audit 02_04), and
-- arterial-line presence is a severity marker (12.7% vs 7.3% mortality),
-- so pooling hands the sicker quarter of the cohort ~5x the measurement
-- density of everyone else — a density confound correlated with the
-- outcome, which `n_obs` carries into pi_hat's shrinkage denominator.
--
-- THAT ARGUMENT IS TRUE BUT WAS OVERSTATED, and the correction is worth
-- recording because it is what makes this rule affordable. `n_obs` counts
-- COVERED HOURS, not raw readings: step 3b collapses each hour to
-- `hr_med` before step 5 counts anything, so n_obs is capped at 24
-- whatever the interface rate. The 5x never reaches the model. What does
-- reach it is the COVERAGE differential — an arterial line fills every
-- hour, a cuff misses some. MEASURED: eICU non-invasive-only mbp gets
-- 20.50 mean covered hours; MIMIC's pooled mbp gets 21.73. The whole
-- effect is worth about three hours out of twenty-four, not a factor of
-- five, and COALESCE keeps it to that.
-- ---------------------------------------------------------------------
SELECT c.stay_id, 'mbp', v.charttime,
       CAST(COALESCE(v.mbp_ni, v.mbp) AS FLOAT64), 'vitalsign'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
  ON c.stay_id = v.stay_id
 AND v.charttime >= c.intime AND v.charttime < c.landmark_time
WHERE COALESCE(v.mbp_ni, v.mbp) IS NOT NULL

UNION ALL
SELECT c.stay_id, 'resp_rate', v.charttime, CAST(v.resp_rate AS FLOAT64), 'vitalsign'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
  ON c.stay_id = v.stay_id
 AND v.charttime >= c.intime AND v.charttime < c.landmark_time
WHERE v.resp_rate IS NOT NULL

UNION ALL
SELECT c.stay_id, 'spo2', v.charttime, CAST(v.spo2 AS FLOAT64), 'vitalsign'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
  ON c.stay_id = v.stay_id
 AND v.charttime >= c.intime AND v.charttime < c.landmark_time
WHERE v.spo2 IS NOT NULL

UNION ALL
-- Expected CELSIUS. No unit filter here - audit B tells you whether that
-- assumption holds. If a large share falls outside 25-43, the fix is
-- CONVERSION in step 3b, not exclusion.
SELECT c.stay_id, 'temperature', v.charttime, CAST(v.temperature AS FLOAT64), 'vitalsign'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` v
  ON c.stay_id = v.stay_id
 AND v.charttime >= c.intime AND v.charttime < c.landmark_time
WHERE v.temperature IS NOT NULL

-- ---------------------------------------------------------------------
-- CHEMISTRY PANEL  (one draw -> sodium, bicarbonate, creatinine, BUN, glucose)
-- hadm_id join, time-filtered to the ICU window, so pre-ICU draws are out.
-- ---------------------------------------------------------------------
UNION ALL
SELECT c.stay_id, 'sodium', ch.charttime, CAST(ch.sodium AS FLOAT64), 'chemistry'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
  ON c.hadm_id = ch.hadm_id
 AND ch.charttime >= c.intime AND ch.charttime < c.landmark_time
WHERE ch.sodium IS NOT NULL

UNION ALL
SELECT c.stay_id, 'bicarbonate', ch.charttime, CAST(ch.bicarbonate AS FLOAT64), 'chemistry'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
  ON c.hadm_id = ch.hadm_id
 AND ch.charttime >= c.intime AND ch.charttime < c.landmark_time
WHERE ch.bicarbonate IS NOT NULL

UNION ALL
SELECT c.stay_id, 'creatinine', ch.charttime, CAST(ch.creatinine AS FLOAT64), 'chemistry'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
  ON c.hadm_id = ch.hadm_id
 AND ch.charttime >= c.intime AND ch.charttime < c.landmark_time
WHERE ch.creatinine IS NOT NULL

UNION ALL
SELECT c.stay_id, 'bun', ch.charttime, CAST(ch.bun AS FLOAT64), 'chemistry'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
  ON c.hadm_id = ch.hadm_id
 AND ch.charttime >= c.intime AND ch.charttime < c.landmark_time
WHERE ch.bun IS NOT NULL

UNION ALL
SELECT c.stay_id, 'glucose', ch.charttime, CAST(ch.glucose AS FLOAT64), 'chemistry'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.chemistry` ch
  ON c.hadm_id = ch.hadm_id
 AND ch.charttime >= c.intime AND ch.charttime < c.landmark_time
WHERE ch.glucose IS NOT NULL

-- ---------------------------------------------------------------------
-- CBC PANEL  (one draw -> WBC, hemoglobin, platelet)
-- ---------------------------------------------------------------------
UNION ALL
SELECT c.stay_id, 'wbc', cbc.charttime, CAST(cbc.wbc AS FLOAT64), 'cbc'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
  ON c.hadm_id = cbc.hadm_id
 AND cbc.charttime >= c.intime AND cbc.charttime < c.landmark_time
WHERE cbc.wbc IS NOT NULL

UNION ALL
SELECT c.stay_id, 'hemoglobin', cbc.charttime, CAST(cbc.hemoglobin AS FLOAT64), 'cbc'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
  ON c.hadm_id = cbc.hadm_id
 AND cbc.charttime >= c.intime AND cbc.charttime < c.landmark_time
WHERE cbc.hemoglobin IS NOT NULL

UNION ALL
SELECT c.stay_id, 'platelet', cbc.charttime, CAST(cbc.platelet AS FLOAT64), 'cbc'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.complete_blood_count` cbc
  ON c.hadm_id = cbc.hadm_id
 AND cbc.charttime >= c.intime AND cbc.charttime < c.landmark_time
WHERE cbc.platelet IS NOT NULL

-- ---------------------------------------------------------------------
-- BLOOD GAS  -> lactate and glucose. bg is scanned once for both.
-- Lactate is INFORMATIVELY MISSING (drawn when shock is suspected);
-- its presence is itself a signal. Note in limitations, do not impute.
-- ---------------------------------------------------------------------
UNION ALL
SELECT c.stay_id, 'lactate', bg.charttime, CAST(bg.lactate AS FLOAT64), 'bg'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.bg` bg
  ON c.hadm_id = bg.hadm_id
 AND bg.charttime >= c.intime AND bg.charttime < c.landmark_time
WHERE bg.lactate IS NOT NULL

UNION ALL
SELECT c.stay_id, 'glucose', bg.charttime, CAST(bg.glucose AS FLOAT64), 'bg'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.bg` bg
  ON c.hadm_id = bg.hadm_id
 AND bg.charttime >= c.intime AND bg.charttime < c.landmark_time
WHERE bg.glucose IS NOT NULL

-- ---------------------------------------------------------------------
-- LIVER: total bilirubin
-- ---------------------------------------------------------------------
UNION ALL
SELECT c.stay_id, 'bilirubin_total', e.charttime,
       CAST(e.bilirubin_total AS FLOAT64), 'enzyme'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_derived.enzyme` e
  ON c.hadm_id = e.hadm_id
 AND e.charttime >= c.intime AND e.charttime < c.landmark_time
WHERE e.bilirubin_total IS NOT NULL

-- ---------------------------------------------------------------------
-- POINT-OF-CARE GLUCOSE  <-- THE EXPENSIVE BLOCK (~11 GB of the ~12 GB)
--
-- Literal itemids, NOT a d_items join: the join forces a full scan, the
-- literal list allows cluster pruning if chartevents is clustered on itemid.
--
-- ***RUN AUDIT A BEFORE TRUSTING THESE ITEMIDS.***
--   220621 = Glucose (serum)
--   225664 = Glucose finger stick
--   226537 = Glucose (whole blood)
-- If audit C shows chemistry + bg already give adequate glucose density,
-- DELETE this block entirely and save the 11 GB.
-- ---------------------------------------------------------------------
UNION ALL
SELECT c.stay_id, 'glucose', ce.charttime, CAST(ce.valuenum AS FLOAT64), 'poc'
FROM cohort c
JOIN `physionet-data.mimiciv_3_1_icu.chartevents` ce
  ON c.stay_id = ce.stay_id
 AND ce.charttime >= c.intime AND ce.charttime < c.landmark_time
WHERE ce.valuenum IS NOT NULL
  AND ce.itemid IN (220621, 225664, 226537, 228388);


-- =====================================================================
-- AUDITS
-- =====================================================================
-- A. ITEMID VERIFICATION - run this FIRST, on d_items only (tiny scan).
--    Confirm the three literals above are right and nothing is missing.
-- SELECT itemid, label, abbreviation, unitname
-- FROM `physionet-data.mimiciv_3_1_icu.d_items`
-- WHERE LOWER(label) LIKE '%glucose%' ORDER BY itemid;
--
-- B. TEMPERATURE UNITS - is derived.vitalsign.temperature really Celsius?
-- SELECT COUNTIF(value BETWEEN 25 AND 43) AS n_celsius,
--        COUNTIF(value BETWEEN 90 AND 110) AS n_fahrenheit,
--        COUNT(*) AS n_total
-- FROM `...v2_raw_values_mimiciv` WHERE signal = 'temperature';
--
-- C. GLUCOSE DENSITY BY SOURCE - does the chartevents read earn its 11 GB?
-- SELECT source, COUNT(*) AS n_rows,
--        COUNT(DISTINCT stay_id) AS n_stays,
--        COUNT(*) / COUNT(DISTINCT stay_id) AS rows_per_stay
-- FROM `...v2_raw_values_mimiciv` WHERE signal = 'glucose' GROUP BY source;
--
-- D. TABLE SIZE - confirm this is a few hundred MB, not 12 GB.
-- SELECT table_name, ROUND(size_bytes/POW(1024,3), 3) AS gb, row_count
-- FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.__TABLES__`
-- WHERE table_id = 'v2_raw_values_mimiciv';
--
-- E. REJECTION AUDIT (now trivial - the unbounded values are a real table).
-- SELECT signal, COUNT(*) AS n_raw, MIN(value), MAX(value)
-- FROM `...v2_raw_values_mimiciv` GROUP BY signal ORDER BY signal;
--
-- F. *** MBP SOURCE AUDIT - RUN THIS BEFORE THE RE-EXTRACTION. ***
--    It decides whether the COALESCE above should instead be strict
--    `mbp_ni`, and it is the MIMIC analogue of eICU audit 18. Three
--    numbers matter:
--      1. how many cohort stays have arterial mbp but NO cuff mbp
--         (strict non-invasive would give these n_obs = 0 -> L = 0)
--      2. their mortality, against the cuff-covered group - if it is
--         elevated the way eICU's is (12.7% vs 7.3%), the arterial-only
--         group is the sick tail and deleting its mbp is not a neutral
--         act
--      3. how much COALESCE actually changes coverage: mean covered
--         hours under `mbp_ni` alone vs under COALESCE(mbp_ni, mbp)
-- WITH v AS (
--   SELECT c.stay_id, c.mortality,
--          COUNTIF(vs.mbp_ni IS NOT NULL)                             AS n_ni,
--          COUNTIF(vs.mbp    IS NOT NULL)                             AS n_any,
--          COUNTIF(vs.mbp IS NOT NULL AND vs.mbp_ni IS NULL)          AS n_art_only
--   FROM `...v2_cohort_mimiciv` c
--   LEFT JOIN `physionet-data.mimiciv_3_1_derived.vitalsign` vs
--     ON c.stay_id = vs.stay_id
--    AND vs.charttime >= c.intime AND vs.charttime < c.landmark_time
--   GROUP BY 1, 2)
-- SELECT
--   COUNT(*)                                              AS n_stays,
--   COUNTIF(n_ni = 0 AND n_any > 0)                       AS n_arterial_only,
--   SAFE_DIVIDE(COUNTIF(n_ni = 0 AND n_any > 0), COUNT(*)) AS frac_arterial_only,
--   AVG(IF(n_ni = 0 AND n_any > 0, mortality, NULL))      AS mort_arterial_only,
--   AVG(IF(n_ni > 0, mortality, NULL))                    AS mort_cuff_covered,
--   AVG(IF(n_ni > 0 AND n_art_only > 0, mortality, NULL)) AS mort_both
-- FROM v;
--
--    DECISION RULE, pre-specified so it is not chosen after seeing the
--    branch point: if frac_arterial_only <= 0.03, either rule is
--    defensible and strict `mbp_ni` is marginally cleaner because it
--    removes the source mixture entirely. Above that, keep the COALESCE
--    and report `mbp` as an instrument mixture at both sites.
--
-- G. COLUMN EXISTENCE. `mbp_ni` is a mimic-code vitalsign column, but
--    confirm it before running the block above - a typo would silently
--    become "arterial only" via the COALESCE.
-- SELECT column_name FROM
--   `physionet-data.mimiciv_3_1_derived.INFORMATION_SCHEMA.COLUMNS`
-- WHERE table_name = 'vitalsign' ORDER BY ordinal_position;
