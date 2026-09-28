-- =====================================================================
-- v2 STEP 4: MIMIC-IV interventions (hourly lattice, long format)
--
-- OUTPUT SHAPE
--   One row per (stay_id, intervention, hour_bin) where the intervention was
--   active, plus an intensity value where dose is meaningful. Deliberately
--   NOT aggregated to stay level: step 5 needs hour-resolved state to
--   (a) censor invalidated measurements and (b) compute the ordering flag O.
--
-- TWO INTERVENTION SHAPES - they are not interchangeable
--   STATE     (vasopressor, vent, RRT, sedation, insulin infusion)
--             -> exposure = fraction of hours active; intensity = peak dose
--   EVENT     (transfusion units, diuretic doses)
--             -> exposure = count over the fixed 24h window
--
-- NO DOSE EQUIVALENCE FOR SEDATION. Unlike vasopressors, propofol /
-- midazolam / dexmedetomidine / fentanyl have no accepted potency
-- equivalence. Exposure only. Benzodiazepine and non-benzodiazepine are
-- kept separate because benzo sedation carries worse outcomes in the
-- delirium literature.
--
-- SCHEMA CAVEAT: verify column names before running. Run
--   SELECT table_name, column_name
--   FROM `physionet-data.mimiciv_3_1_derived.INFORMATION_SCHEMA.COLUMNS`
--   WHERE table_name IN ('norepinephrine_equivalent_dose','ventilation',
--                        'rrt','ventilator_setting');
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_interventions_mimiciv` AS

WITH cohort AS (
  SELECT *
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

-- Hour spine: every stay gets all 24 slots, so absence is explicit.
hours AS (
  SELECT
    c.stay_id,
    h AS hour_bin,
    TIMESTAMP_ADD(c.intime, INTERVAL h     HOUR) AS hr_start,
    TIMESTAMP_ADD(c.intime, INTERVAL h + 1 HOUR) AS hr_end
  FROM cohort c, UNNEST(GENERATE_ARRAY(0, 23)) AS h
),

-- Label lookup. Matching on d_items labels rather than hard-coded itemids:
-- itemids drift between MIMIC releases, labels are self-documenting, and a
-- silently missed variant truncates coverage without any error.
items AS (
  SELECT itemid, label, LOWER(label) AS l
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
),

-- ---------------------------------------------------------------------
-- VASOPRESSORS - norepinephrine-equivalent dose (STATE + intensity)
-- Uses the derived table so NEE conversion is not hand-rolled.
-- ---------------------------------------------------------------------
vaso AS (
  SELECT
    p.stay_id,
    'vasopressor' AS intervention,
    p.hour_bin,
    MAX(n.norepinephrine_equivalent_dose) AS intensity
  FROM (
    -- PRESENCE: authoritative, from inputevents label matching.
    SELECT DISTINCT stay_id, hour_bin
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_agents_mimiciv`
    WHERE intervention = 'vasopressor'
  ) p
  LEFT JOIN cohort c
    ON p.stay_id = c.stay_id
  -- INTENSITY: LEFT JOIN, so a stay with no computable NEE keeps its
  -- presence row and simply carries a NULL dose. Under the old INNER-join
  -- semantics that stay vanished entirely.
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.norepinephrine_equivalent_dose` n
    ON p.stay_id    = n.stay_id
   AND n.starttime  <  TIMESTAMP_ADD(c.intime, INTERVAL p.hour_bin + 1 HOUR)
   AND n.endtime    >  TIMESTAMP_ADD(c.intime, INTERVAL p.hour_bin     HOUR)
   AND n.norepinephrine_equivalent_dose > 0
  GROUP BY 1, 2, 3
),

-- ---------------------------------------------------------------------
-- VENTILATION (STATE, ordinal). Invasive is kept distinct from NIV/HFNC:
-- collapsing them loses the clinically decisive difference.
-- ---------------------------------------------------------------------
vent AS (
  SELECT
    h.stay_id,
    CASE
      WHEN MAX(CASE WHEN v.ventilation_status IN ('InvasiveVent','Tracheostomy')
                    THEN 1 ELSE 0 END) = 1 THEN 'invasive_vent'
      WHEN MAX(CASE WHEN v.ventilation_status IN ('NonInvasiveVent','HFNC')
                    THEN 1 ELSE 0 END) = 1 THEN 'noninvasive_vent'
      ELSE 'supplemental_o2'
    END AS intervention,
    h.hour_bin,
    CAST(NULL AS FLOAT64) AS intensity
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_derived.ventilation` v
    ON h.stay_id = v.stay_id
   AND v.starttime <  h.hr_end
   AND v.endtime   >  h.hr_start
  WHERE v.ventilation_status IS NOT NULL
    AND v.ventilation_status <> 'None'
  GROUP BY h.stay_id, h.hour_bin
),

-- Ventilator settings: FiO2 and PEEP as separate intensity channels.
-- FiO2 stays on the INTERVENTION side only. Do not also build P/F ratio -
-- that would put FiO2 on both sides of the oxygenation pair.
vent_settings AS (
  SELECT
    h.stay_id,
    'fio2' AS intervention,
    h.hour_bin,
    MAX(vs.fio2) AS intensity
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_derived.ventilator_setting` vs
    ON h.stay_id = vs.stay_id
   AND vs.charttime >= h.hr_start
   AND vs.charttime <  h.hr_end
  WHERE vs.fio2 IS NOT NULL
  GROUP BY 1, 2, 3

  UNION ALL

  SELECT
    h.stay_id,
    'peep',
    h.hour_bin,
    MAX(vs.peep)
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_derived.ventilator_setting` vs
    ON h.stay_id = vs.stay_id
   AND vs.charttime >= h.hr_start
   AND vs.charttime <  h.hr_end
  WHERE vs.peep IS NOT NULL
  GROUP BY 1, 2, 3
),

-- ---------------------------------------------------------------------
-- RRT (STATE). Critical for step 5: this is the mask that INVALIDATES
-- creatinine and urine output, not merely attenuates them.
-- ---------------------------------------------------------------------
rrt AS (
  SELECT
    h.stay_id,
    'rrt' AS intervention,
    h.hour_bin,
    CAST(NULL AS FLOAT64) AS intensity
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_derived.rrt` r
    ON h.stay_id = r.stay_id
   AND r.charttime >= h.hr_start
   AND r.charttime <  h.hr_end
  WHERE COALESCE(r.dialysis_active, 0) = 1
  GROUP BY 1, 2, 3
),

-- ---------------------------------------------------------------------
-- INFUSION-BASED interventions from inputevents (STATE).
-- Sedation split benzo / non-benzo. No dose equivalence attempted.
-- ---------------------------------------------------------------------
infusions AS (
  SELECT
    h.stay_id,
    CASE
      WHEN i.l LIKE '%midazolam%' OR i.l LIKE '%lorazepam%'
        THEN 'sedation_benzo'
      WHEN i.l LIKE '%dexmedetomidine%' OR i.l LIKE '%precedex%'
        THEN 'sedation_dexmed'
      WHEN i.l LIKE '%propofol%'
        THEN 'sedation_propofol'
      WHEN i.l LIKE '%fentanyl%' OR i.l LIKE '%morphine%'
           OR i.l LIKE '%hydromorphone%'
        THEN 'opioid'
      WHEN i.l LIKE '%insulin%'
        THEN 'insulin'
      WHEN i.l LIKE '%cisatracurium%' OR i.l LIKE '%rocuronium%'
           OR i.l LIKE '%vecuronium%'
        THEN 'paralytic'
      WHEN i.l LIKE '%dobutamine%' OR i.l LIKE '%milrinone%'
         THEN 'inotrope'
    END AS intervention,
    h.hour_bin,
    MAX(ie.rate) AS intensity
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie
    ON h.stay_id = ie.stay_id
   AND ie.starttime <  h.hr_end
   AND ie.endtime   >  h.hr_start
  JOIN items i
    ON ie.itemid = i.itemid
  WHERE ie.amount > 0
  -- Route restriction, insulin ONLY. Subcutaneous basal and sliding-scale
    -- insulin is routine glycemic maintenance, not an ICU intervention, and
    -- it corrupts the glucose/insulin O flag. Every other drug here keeps
    -- both routes: benzo and paralytic boluses drive GCS masking and must
    -- not be dropped.
    AND NOT (i.l LIKE '%insulin%' AND ie.rate IS NULL)
    AND (i.l LIKE '%midazolam%' OR i.l LIKE '%lorazepam%'
      OR i.l LIKE '%propofol%'  OR i.l LIKE '%dexmedetomidine%'
      OR i.l LIKE '%precedex%'
      OR i.l LIKE '%fentanyl%'  OR i.l LIKE '%morphine%'
      OR i.l LIKE '%hydromorphone%'
      OR i.l LIKE '%insulin%'
      OR i.l LIKE '%cisatracurium%' OR i.l LIKE '%rocuronium%'
      OR i.l LIKE '%vecuronium%'
      OR i.l LIKE '%dobutamine%' OR i.l LIKE '%milrinone%')
  GROUP BY 1, 2, 3
),

-- ---------------------------------------------------------------------
-- EVENT-BASED interventions (discrete). intensity = amount given in
-- that hour; exposure downstream is a COUNT over the window, not a fraction.
-- ---------------------------------------------------------------------
events AS (
  SELECT
    h.stay_id,
    CASE
      WHEN i.l LIKE '%furosemide%' OR i.l LIKE '%bumetanide%'
        THEN 'diuretic'
      WHEN i.l LIKE '%platelet%'
        THEN 'transfusion_platelet'
      WHEN i.l LIKE '%packed red blood%' OR i.l LIKE '%prbc%'
        THEN 'transfusion_prbc'
      WHEN i.l LIKE '%fresh frozen plasma%' OR i.l LIKE '%ffp%'
        THEN 'transfusion_ffp'
    END AS intervention,
    h.hour_bin,
    SUM(ie.amount) AS intensity
  FROM hours h
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie
    ON h.stay_id = ie.stay_id
   AND ie.starttime >= h.hr_start
   AND ie.starttime <  h.hr_end
  JOIN items i
    ON ie.itemid = i.itemid
  WHERE ie.amount > 0
    AND (i.l LIKE '%furosemide%' OR i.l LIKE '%bumetanide%'
      OR i.l LIKE '%platelet%'
      OR i.l LIKE '%packed red blood%' OR i.l LIKE '%prbc%'
      OR i.l LIKE '%fresh frozen plasma%' OR i.l LIKE '%ffp%')
  GROUP BY 1, 2, 3
),

combined AS (
  SELECT stay_id, intervention, hour_bin, intensity, 'state' AS shape FROM vaso
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'state' FROM vent
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'state' FROM vent_settings
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'state' FROM rrt
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'state' FROM infusions
  WHERE intervention IS NOT NULL
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'event' FROM events
  WHERE intervention IS NOT NULL
)

SELECT
  c.stay_id,
  c.intervention,
  c.hour_bin,
  c.shape,
  c.intensity,
  1 AS active
FROM combined c;


-- =====================================================================
-- 5b. INTERVENTION FEATURES  (unchanged from rev 1)
--
-- Two shapes, not interchangeable:
--   shape = 'state' -> use exposure_frac and peak_intensity
--   shape = 'event' -> use n_hours and total_amount; exposure_frac is
--                      MEANINGLESS for these (diuretic 0.064 means ~1.5
--                      hours with a dose, not 6% exposure)
-- =====================================================================
CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_intervention_features_mimiciv` AS

WITH cohort AS (
  SELECT stay_id
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

iv AS (
  SELECT *
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_interventions_mimiciv`
),

intervention_list AS (
  SELECT DISTINCT intervention, shape FROM iv
),

grid AS (
  SELECT c.stay_id, l.intervention, l.shape
  FROM cohort c CROSS JOIN intervention_list l
),

-- per-stay agent identity
agent_ids AS (
  SELECT
    stay_id,
    intervention,
    COUNT(DISTINCT agent)                    AS n_agents,
    MAX(IF(agent = 'norepinephrine', 1, 0))  AS has_norepinephrine
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_agents_mimiciv`
  GROUP BY 1, 2
),

-- per-hour concurrency, then reduced to a stay-level max.
-- Two stages because COUNT(DISTINCT ...) OVER (...) is rejected by BigQuery.
agent_concurrency AS (
  SELECT stay_id, intervention, MAX(n_concurrent) AS max_concurrent_agents
  FROM (
    SELECT stay_id, intervention, hour_bin, COUNT(DISTINCT agent) AS n_concurrent
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_agents_mimiciv`
    GROUP BY 1, 2, 3
  )
  GROUP BY 1, 2
),

rolled AS (
  SELECT
    g.stay_id,
    g.intervention,
    g.shape,
    -- COUNT DISTINCT, not COUNT: guards against any duplicate
    -- (stay, intervention, hour_bin) arriving from the UNION ALL in 02
    COUNT(DISTINCT iv.hour_bin)                              AS n_hours,
    SUM(iv.intensity)                                        AS total_amount,
    MAX(iv.intensity)                                        AS raw_peak,
    MIN(iv.hour_bin)                                         AS first_hour,
    CASE WHEN COUNT(iv.hour_bin) = 0 THEN 0 ELSE 1 END       AS ever_active,
    CASE WHEN MIN(iv.hour_bin) = 0   THEN 1 ELSE 0 END       AS present_at_admission
  FROM grid g
  LEFT JOIN iv
    ON g.stay_id = iv.stay_id AND g.intervention = iv.intervention
  GROUP BY g.stay_id, g.intervention, g.shape
)

SELECT
  'mimic'                                                    AS site,
  r.stay_id,
  r.intervention,
  r.shape,
  r.ever_active,

  -- ACTIVITY MEASURES: 0 when never active, never NULL. A GAM needs the
  -- unexposed patients evaluated at zero, not dropped as missing.
  IF(r.shape = 'state', r.n_hours / 24.0, NULL)              AS exposure_frac,
  IF(r.shape = 'event', r.n_hours,        NULL)              AS n_hours,
  IF(r.shape = 'event', COALESCE(r.total_amount, 0), NULL)   AS total_amount,

  -- DOSE MEASURES: NULL when never active. "The dose he was not given" is
  -- undefined, not zero. See the R-SIDE DECISION note below — this is the
  -- one column whose NULL handling must be settled deliberately rather
  -- than defaulted.
  IF(r.intervention IN ('vasopressor','inotrope'), NULL, r.raw_peak)
                                                             AS peak_intensity,
  IF(r.intervention = 'vasopressor', r.raw_peak, NULL)       AS dx_nee_peak,

  -- AGENT MEASURES: 0 for the multi-agent interventions when unexposed,
  -- NULL for every other intervention where the concept does not apply.
  IF(r.intervention IN ('vasopressor','inotrope'),
     COALESCE(ai.n_agents, 0),             NULL)             AS n_agents,
  IF(r.intervention IN ('vasopressor','inotrope'),
     COALESCE(ac.max_concurrent_agents, 0), NULL)            AS max_concurrent_agents,
  -- vasopressor only. A priori justification: norepinephrine is the
  -- guideline first-line agent for shock at both sites, so its presence
  -- marks that the patient was TREATED AS being in shock.
  IF(r.intervention = 'vasopressor',
     COALESCE(ai.has_norepinephrine, 0),   NULL)             AS has_norepinephrine,

  r.first_hour,
  r.present_at_admission
FROM rolled r
LEFT JOIN agent_ids         ai ON r.stay_id = ai.stay_id AND r.intervention = ai.intervention
LEFT JOIN agent_concurrency ac ON r.stay_id = ac.stay_id AND r.intervention = ac.intervention;


-- =====================================================================
-- CONVENIENCE ROLL-UP (step 5 will redo this with the ordering flag)
-- =====================================================================
-- CREATE OR REPLACE TABLE `...v2_intervention_summary_mimiciv` AS
-- SELECT
--   stay_id,
--   intervention,
--   ANY_VALUE(shape)                        AS shape,
--   COUNT(DISTINCT hour_bin) / 24.0         AS exposure_frac,   -- states
--   COUNT(DISTINCT hour_bin)                AS n_hours,
--   SUM(intensity)                          AS total_amount,    -- events
--   MAX(intensity)                          AS peak_intensity,
--   MIN(hour_bin)                           AS first_hour,
--   CASE WHEN MIN(hour_bin) = 0 THEN 1 ELSE 0 END AS present_at_admission
-- FROM `...v2_interventions_mimiciv`
-- GROUP BY stay_id, intervention;
--
-- CHECKS
-- 1. Prevalence per intervention. Vasopressor ~25-35%, invasive vent ~40%,
--    RRT ~5-10% are the rough expectations; large deviations mean the label
--    matching missed variants.
-- 2. Label audit - confirm nothing was silently dropped:
--    SELECT DISTINCT label FROM `physionet-data.mimiciv_3_1_icu.d_items`
--    WHERE LOWER(label) LIKE '%propofol%' OR LOWER(label) LIKE '%fentanyl%';
-- 3. Cross-check vent_invasive against gcs_unable_ever from step 3.
