-- =====================================================================
-- v2 STEP 1: MIMIC-IV cohort
--
-- Design decisions (frozen before analysis):
--   * First ICU stay PER PATIENT (PARTITION BY subject_id) -> one row per
--     person, full independence. NOT per hospitalization.
--   * Landmark at intime + 24h.
--   * Requires ICU outtime >= landmark. This is the key change from v1:
--     v1 only required HOSPITAL discharge after the landmark, which admitted
--     patients with <24h of ICU observation into a 24h feature window.
--   * Deaths before the landmark excluded (landmark design; documented as a
--     selection effect - the earliest deaths are removed).
--   * Age from mimiciv_derived.age (age at THIS admission), not anchor_age
--     (age at the anchor year).
--
-- Added vs v1: first_careunit, anchor_year_group, weight, race.
--   - first_careunit    -> leave-one-unit-out internal-external validation
--   - anchor_year_group -> temporal validation
--   - weight            -> required for urine output (mL/kg/hr)
--   - race              -> subgroup analysis ONLY, never enters the score
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv` AS

WITH base AS (
  SELECT
    icu.subject_id,
    icu.hadm_id,
    icu.stay_id,
    icu.intime,
    icu.outtime,
    icu.first_careunit,
    icu.last_careunit,
    icu.los                              AS icu_los_days,

    adm.admittime,
    adm.dischtime,
    adm.deathtime,
    adm.hospital_expire_flag             AS mortality,
    adm.discharge_location,
    adm.insurance,
    adm.race,

    pat.gender,
    pat.anchor_year_group,
    age.age,

    ROW_NUMBER() OVER (
      PARTITION BY icu.subject_id
      ORDER BY icu.intime
    ) AS rn

  FROM `physionet-data.mimiciv_3_1_icu.icustays` icu
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm
    ON icu.hadm_id = adm.hadm_id
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` pat
    ON icu.subject_id = pat.subject_id
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.age` age
    ON icu.hadm_id = age.hadm_id
),

first_stay AS (
  SELECT
    *,
    TIMESTAMP_ADD(intime, INTERVAL 24 HOUR) AS landmark_time
  FROM base
  WHERE rn = 1
),

filtered AS (
  SELECT *
  FROM first_stay
  WHERE age >= 18
    AND dischtime IS NOT NULL
    AND outtime   IS NOT NULL
    -- 24h of actual ICU observation must exist
    AND outtime   >= landmark_time
    AND dischtime >= landmark_time
    -- landmark design: exclude deaths inside the observation window
    AND (deathtime IS NULL OR deathtime >= landmark_time)
),

wt AS (
  SELECT
    stay_id,
    COALESCE(weight_admit, weight) AS weight_kg
  FROM `physionet-data.mimiciv_3_1_derived.first_day_weight`
)

SELECT
  f.subject_id,
  f.hadm_id,
  f.stay_id,

  f.intime,
  f.outtime,
  f.landmark_time,
  f.admittime,
  f.dischtime,
  f.deathtime,

  f.mortality,
  f.discharge_location,

  -- Stratifiers / covariates. NONE of these enter the primary score.
  f.gender,
  LEAST(f.age, 90) age,
  f.race,
  f.insurance,
  f.first_careunit,
  f.last_careunit,
  f.anchor_year_group,
  f.icu_los_days,
  CASE WHEN w.weight_kg >= 30 AND w.weight_kg <= 300 THEN w.weight_kg ELSE NULL END AS weight_kg,

  -- ---------------------------------------------------------------
  -- Data-quality flags (inspect these before trusting the cohort)
  -- ---------------------------------------------------------------
  -- Died in hospital per flag but no deathtime recorded: such rows CANNOT be
  -- excluded by the landmark rule and may include pre-landmark deaths.
  CASE WHEN f.mortality = 1 AND f.deathtime IS NULL THEN 1 ELSE 0 END
    AS qc_death_no_timestamp,

  -- Discharged to hospice: coded alive, effectively a death. Hospice referral
  -- rates vary by site, so this biases outcome rates site-dependently.
  CASE WHEN LOWER(f.discharge_location) LIKE '%hospice%' THEN 1 ELSE 0 END
    AS qc_discharge_hospice,

  CASE WHEN w.weight_kg IS NULL THEN 1 ELSE 0 END AS qc_weight_missing,

  -- ---------------------------------------------------------------
  -- Survival endpoints, clock starting AT THE LANDMARK
  -- ---------------------------------------------------------------
  CASE
    WHEN f.deathtime IS NOT NULL AND f.deathtime >= f.landmark_time THEN 1
    ELSE 0
  END AS event_after_24h,

  CASE
    WHEN f.deathtime IS NOT NULL AND f.deathtime >= f.landmark_time
      THEN TIMESTAMP_DIFF(f.deathtime, f.landmark_time, HOUR)
    ELSE TIMESTAMP_DIFF(f.dischtime, f.landmark_time, HOUR)
  END AS duration_hours_from_24h,

  CASE
    WHEN f.deathtime IS NOT NULL
      THEN TIMESTAMP_DIFF(f.deathtime, f.intime, HOUR)
    ELSE NULL
  END AS time_to_death_from_icu_hours

FROM filtered f
LEFT JOIN wt w
  ON f.stay_id = w.stay_id;


-- =====================================================================
-- POST-EXTRACTION CHECKS - run these, do not skip
-- =====================================================================
-- SELECT
--   COUNT(*)                             AS n_stays,
--   SUM(mortality)                       AS n_deaths,
--   AVG(mortality)                       AS mortality_rate,
--   SUM(qc_death_no_timestamp)           AS n_death_no_timestamp,
--   SUM(qc_discharge_hospice)            AS n_hospice,
--   SUM(qc_weight_missing)               AS n_weight_missing,
--   MIN(age), MAX(age)
-- FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`;
--
-- -- Cohort attrition: quantify what the ICU LOS filter removed.
-- -- Report this in the paper; it is the single largest selection effect.
