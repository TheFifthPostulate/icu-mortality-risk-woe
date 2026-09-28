-- =====================================================================
-- v2 STEP 1: eICU-CRD cohort
--
-- PORT OF sql/mimiciv/v2_01_cohort_mimiciv.sql. Every design decision is
-- carried over unchanged; only the source columns differ. Where eICU
-- cannot express a MIMIC concept the substitution is named in the
-- SITE SUBSTITUTIONS block below and recorded in
-- sql/eicu/eicu_extraction_decisions.md. Nothing is substituted silently.
--
-- OUTPUT CONTRACT. This table must emit EXACTLY the 26 columns of
-- R/02_validate.R CONTRACT$cohort (minus `site`, which the loader stamps),
-- in the same types. Check 9 compares types across sites, so a STRING
-- where MIMIC has INT64 fails the run rather than producing a
-- transportability artifact. Do not add columns here - see the
-- `v2_cohort_hospital_eicu` side table at the bottom for hospitalid.
--
-- ---------------------------------------------------------------------
-- SITE SUBSTITUTIONS
-- ---------------------------------------------------------------------
--   subject_id          DENSE_RANK() over `uniquepid`. eICU's patient key
--                       is a STRING ('002-33870'); MIMIC's is INT64 and
--                       check 9 compares types. The rank is deterministic
--                       for a fixed cohort and is used ONLY for fold
--                       grouping (config: split.group_by = patient_id),
--                       which needs an equivalence class, not an identity.
--   hadm_id             patienthealthsystemstayid (the hospital stay).
--   stay_id             patientunitstayid.
--   TIMESTAMPS          eICU HAS NO ABSOLUTE TIME. Every offset is minutes
--                       from unit admission. All six timestamp columns are
--                       synthesised from ONE shared epoch, so `intime` is
--                       identical for every stay. This is safe because
--                       nothing downstream compares absolute times ACROSS
--                       stays - every consumer takes a difference within a
--                       stay (hour_bin, landmark, duration). It would NOT
--                       be safe for any calendar-time analysis, and there
--                       is none in the design.
--   age                 '> 89' -> 90, matching the MIMIC-side LEAST(age,90).
--                       Blank age (95 stays) is non-numeric and fails the
--                       >= 18 filter, so it drops with the paediatric rows.
--   gender              'Female'/'Male' -> 'F'/'M'. The loader errors on
--                       anything else, so stays with a blank / 'Unknown' /
--                       'Other' gender are EXCLUDED here. Count them in
--                       check C3 and report the exclusion.
--   race                `ethnicity`.
--   insurance           DOES NOT EXIST in eICU. Emitted as a typed NULL so
--                       the column set matches. Never modelled at either
--                       site.
--   first/last_careunit `unittype`, twice. eICU has one unit type per unit
--                       stay, so the transfer distinction MIMIC carries
--                       does not exist. Leave-one-unit-out validation is
--                       therefore coarser at eICU; state it.
--   anchor_year_group   CAST(hospitaldischargeyear AS STRING). Takes
--                       exactly two values (2014, 2015), so temporal
--                       validation at eICU is a two-level contrast, not
--                       MIMIC's four.
--   weight_kg           `patient.admissionweight`, same 30-300 guard.
--                       `infusiondrug.patientweight` is NOT a fallback:
--                       null in 4.36M of 4.80M rows (audit 17c).
--   qc_death_no_timestamp  STRUCTURALLY ALWAYS 0. eICU takes `dischtime`
--                       and `deathtime` from the SAME field
--                       (hospitaldischargeoffset), so a death with no
--                       timestamp cannot survive the non-null filter that
--                       MIMIC also applies. The column and its expression
--                       are kept verbatim; only the reachable value set
--                       differs.
--   qc_discharge_hospice   STRUCTURALLY ALWAYS 0. eICU's
--                       `hospitaldischargelocation` has NO hospice
--                       category anywhere (audit 13a/13b: zero across all
--                       158 hospitals). The expression is kept identical
--                       so no code branches on site; the consequence is
--                       that the pre-registered hospice-as-death
--                       sensitivity analysis is MIMIC-INTERNAL ONLY, and
--                       that the outcome label is differentially
--                       misclassified between development and validation.
--                       MIMIC has 1,207 identifiable hospice survivors
--                       (23% of its deaths); eICU has an unknown number
--                       inside 'Other External' (9,757) and 'Other'
--                       (7,888). Part of the 10.3% vs 8.78% mortality gap
--                       is therefore LABEL DEFINITION, not case mix. Say
--                       this before a reviewer says it.
--
-- ---------------------------------------------------------------------
-- FIRST STAY PER PATIENT - the one rule that needed a new decision
-- ---------------------------------------------------------------------
--   MIMIC: ROW_NUMBER() PARTITION BY subject_id ORDER BY intime.
--   eICU has no intime, and `unitvisitnumber` RESETS inside each
--   patienthealthsystemstayid, so 18,881 patients with >= 2 hospital
--   stays contribute multiple rows at unitvisitnumber = 1 - about 17,000
--   rows past the intended one-per-patient rule (audit 12c).
--   The only ordering information in the whole dataset is
--   `hospitaldischargeyear`, which has two values and cannot break
--   within-year ties (audit 12d).
--
--   RULE, pre-specified and deterministic:
--     PARTITION BY uniquepid
--     ORDER BY hospitaldischargeyear, patienthealthsystemstayid,
--              unitvisitnumber, patientunitstayid
--   Earliest discharge year, then lowest hospital-stay id, then the first
--   unit visit inside it. The tie-break is ARBITRARY and is documented as
--   arbitrary. Repeat-hospitalisation patients are NOT dropped: dropping
--   them would select against chronic illness, which interacts with
--   mortality in exactly the direction that would flatter the validation.
-- =====================================================================

DECLARE epoch TIMESTAMP DEFAULT TIMESTAMP '2014-01-01 00:00:00';

CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_cohort_eicu` AS

WITH base AS (
  SELECT
    p.patientunitstayid,
    p.patienthealthsystemstayid,
    p.uniquepid,
    p.gender,
    p.ethnicity,
    p.unittype,
    p.unitvisitnumber,
    p.hospitaldischargeyear,
    p.hospitaladmitoffset,
    p.hospitaldischargeoffset,
    p.hospitaldischargelocation,
    p.hospitaldischargestatus,
    p.unitdischargeoffset,
    p.admissionweight,

    -- '> 89' is the only non-numeric age string (7,081 stays); blanks (95)
    -- become NULL and fail the adult filter.
    CASE WHEN TRIM(p.age) = '> 89' THEN 90
         ELSE SAFE_CAST(p.age AS INT64) END      AS age_years,

    ROW_NUMBER() OVER (
      PARTITION BY p.uniquepid
      ORDER BY p.hospitaldischargeyear,
               p.patienthealthsystemstayid,
               p.unitvisitnumber,
               p.patientunitstayid
    ) AS rn

  FROM `physionet-data.eicu_crd.patient` p
),

first_stay AS (
  SELECT * FROM base WHERE rn = 1
),

filtered AS (
  SELECT *
  FROM first_stay
  WHERE age_years >= 18
    -- outcome must be known; blank status is 1,751 stays (audit 13a)
    AND hospitaldischargestatus IN ('Alive', 'Expired')
    AND hospitaldischargeoffset IS NOT NULL
    AND unitdischargeoffset     IS NOT NULL
    -- 24h of actual ICU observation must exist. This is the dominant cut
    -- (30%, audit 12a) and belongs in the CONSORT diagram with mortality
    -- at each rung, not buried.
    AND unitdischargeoffset     >= 1440
    AND hospitaldischargeoffset >= 1440
    -- landmark design: exclude deaths inside the observation window.
    -- Removes only 18 stays once the LOS filter has run; keep it for
    -- correctness, do not present it as a meaningful exclusion.
    AND NOT (hospitaldischargestatus = 'Expired' AND hospitaldischargeoffset < 1440)
    -- the loader factors gender to {F, M} and errors on anything else
    AND gender IN ('Female', 'Male')
)

SELECT
  DENSE_RANK() OVER (ORDER BY f.uniquepid)                        AS subject_id,
  f.patienthealthsystemstayid                                     AS hadm_id,
  f.patientunitstayid                                             AS stay_id,

  epoch                                                           AS intime,
  TIMESTAMP_ADD(epoch, INTERVAL f.unitdischargeoffset MINUTE)     AS outtime,
  TIMESTAMP_ADD(epoch, INTERVAL 24 HOUR)                          AS landmark_time,
  TIMESTAMP_ADD(epoch, INTERVAL f.hospitaladmitoffset MINUTE)     AS admittime,
  TIMESTAMP_ADD(epoch, INTERVAL f.hospitaldischargeoffset MINUTE) AS dischtime,
  CASE WHEN f.hospitaldischargestatus = 'Expired'
       THEN TIMESTAMP_ADD(epoch, INTERVAL f.hospitaldischargeoffset MINUTE)
       ELSE NULL END                                              AS deathtime,

  CASE WHEN f.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS mortality,
  f.hospitaldischargelocation                                     AS discharge_location,

  -- Stratifiers / covariates. NONE of these enter the primary score.
  CASE WHEN f.gender = 'Female' THEN 'F' ELSE 'M' END             AS gender,
  LEAST(f.age_years, 90)                                          AS age,
  f.ethnicity                                                     AS race,
  CAST(NULL AS STRING)                                            AS insurance,
  f.unittype                                                      AS first_careunit,
  f.unittype                                                      AS last_careunit,
  CAST(f.hospitaldischargeyear AS STRING)                         AS anchor_year_group,
  f.unitdischargeoffset / 1440.0                                  AS icu_los_days,
  -- Guard applied in SQL, as at MIMIC. R/01_load.R re-asserts it.
  CASE WHEN f.admissionweight BETWEEN 30 AND 300
       THEN CAST(f.admissionweight AS FLOAT64) ELSE NULL END      AS weight_kg,

  -- ---------------------------------------------------------------
  -- Data-quality flags. Expressions IDENTICAL to the MIMIC file; two of
  -- the three are structurally unreachable at eICU (see header).
  -- ---------------------------------------------------------------
  CASE WHEN f.hospitaldischargestatus = 'Expired'
        AND f.hospitaldischargeoffset IS NULL THEN 1 ELSE 0 END
    AS qc_death_no_timestamp,

  CASE WHEN LOWER(f.hospitaldischargelocation) LIKE '%hospice%' THEN 1 ELSE 0 END
    AS qc_discharge_hospice,

  CASE WHEN f.admissionweight IS NULL
         OR NOT (f.admissionweight BETWEEN 30 AND 300) THEN 1 ELSE 0 END
    AS qc_weight_missing,

  -- ---------------------------------------------------------------
  -- Survival endpoints, clock starting AT THE LANDMARK
  -- ---------------------------------------------------------------
  CASE WHEN f.hospitaldischargestatus = 'Expired'
        AND f.hospitaldischargeoffset >= 1440 THEN 1 ELSE 0 END
    AS event_after_24h,

  -- Both branches are the same expression because eICU's death time and
  -- discharge time are the same field. Written out in full anyway, so the
  -- structure matches the MIMIC file line for line and a future eICU
  -- release that separates them needs no restructuring.
  -- FLOOR, not CAST-round: TIMESTAMP_DIFF(..., HOUR) at MIMIC TRUNCATES, and
  -- CAST(x/60 AS INT64) in BigQuery rounds half away from zero. Both offsets
  -- here are non-negative after the landmark filter, so FLOOR reproduces
  -- MIMIC's arithmetic exactly.
  CAST(
    CASE WHEN f.hospitaldischargestatus = 'Expired'
         THEN FLOOR((f.hospitaldischargeoffset - 1440) / 60)
         ELSE FLOOR((f.hospitaldischargeoffset - 1440) / 60)
    END AS INT64)                                                 AS duration_hours_from_24h,

  CASE WHEN f.hospitaldischargestatus = 'Expired'
       THEN CAST(FLOOR(f.hospitaldischargeoffset / 60) AS INT64)
       ELSE NULL END                                              AS time_to_death_from_icu_hours

FROM filtered f;


-- =====================================================================
-- SIDE TABLE: hospital identity, for the leave-one-hospital-out
-- heterogeneity analysis (reconciliation section N).
--
-- DELIBERATELY NOT IN THE COHORT TABLE. CONTRACT$cohort is an exact
-- column set and validator check 1 fails on any extra column, so adding
-- `hospitalid` to the cohort parquet is a SPEC CHANGE to discuss, not a
-- convenience. Until that decision is made this table carries it, is not
-- listed in config/external.yml, and is never loaded by R/01_load.R.
--
-- Hospital inclusion criterion is UNDECIDED and must be pre-specified
-- before fitting. The two defensible options (reconciliation N):
--   (a) exclude hospitals with GCS or urine coverage < 0.5 (costs ~40,000
--       stays, leaves ~70,000 - still larger than MIMIC);
--   (b) keep all hospitals and report leave-one-hospital-out
--       heterogeneity with per-hospital coverage as the explanatory
--       covariate.
-- (b) is the stronger paper and turns reconciliation C and D from
-- limitations into the mechanism that explains the heterogeneity. This
-- extraction implements NEITHER - it keeps every hospital and hands the
-- decision downstream, which is what (b) needs and what (a) can still be
-- applied on top of.
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_cohort_hospital_eicu` AS
SELECT
  c.stay_id,
  p.hospitalid,
  p.unittype,
  p.unitadmitsource,
  p.apacheadmissiondx
FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu` c
JOIN `physionet-data.eicu_crd.patient` p
  ON c.stay_id = p.patientunitstayid;


-- =====================================================================
-- POST-EXTRACTION CHECKS - run these, do not skip
-- =====================================================================
-- C1. Headline, against the audit-12 ladder (109,200 stays, 8.78%).
--     A material shortfall means the gender filter or the first-stay rule
--     cost more than expected - quantify before accepting the number.
-- SELECT
--   COUNT(*)                    AS n_stays,
--   COUNT(DISTINCT subject_id)  AS n_patients,
--   SUM(mortality)              AS n_deaths,
--   AVG(mortality)              AS mortality_rate,
--   SUM(qc_death_no_timestamp)  AS n_death_no_timestamp,   -- must be 0
--   SUM(qc_discharge_hospice)   AS n_hospice,              -- must be 0
--   SUM(qc_weight_missing)      AS n_weight_missing,
--   MIN(age), MAX(age)
-- FROM `...v2_cohort_eicu`;
--
-- C2. n_stays MUST equal n_patients. If it does not, the first-stay rule
--     is not doing what the header says and the fold grouping is unsound.
--     Audit 12c predicts ~17,000 extra rows if the rule is wrong, so this
--     is a loud failure, not a subtle one.
--
-- C3. Gender exclusion count - the one stay-level exclusion this port
--     adds that MIMIC does not have. Report it in the CONSORT diagram.
-- SELECT gender, COUNT(*) FROM `physionet-data.eicu_crd.patient`
-- WHERE unitvisitnumber = 1 AND unitdischargeoffset >= 1440 GROUP BY 1;
--
-- C4. Attrition ladder, to sit beside audit 12a. Rebuild it against THIS
--     file's filter chain, not the discovery query's, because the two
--     differ in the gender and first-stay rules.
--
-- C5. Weight coverage. Audit 03a: 16,718/200,859 null (8.3%), 255
--     implausible. Stays with no usable weight yield NO urine-output rows
--     downstream (n_obs = 0), which is correct - unavailable, not zero -
--     but the fraction is a site difference worth reporting because
--     urine_output_rate is a modelled signal.
-- SELECT AVG(qc_weight_missing) FROM `...v2_cohort_eicu`;
--
-- C6. Timestamp sanity. Every stay must have intime = the epoch,
--     landmark_time = intime + 24h, and outtime >= landmark_time.
-- SELECT COUNTIF(outtime < landmark_time)   AS bad_outtime,
--        COUNTIF(dischtime < landmark_time) AS bad_dischtime,
--        COUNT(DISTINCT intime)             AS n_distinct_intime  -- must be 1
-- FROM `...v2_cohort_eicu`;
