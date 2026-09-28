-- =====================================================================
-- v2 STEP 2b: eICU interventions (hourly lattice, long format)
--             + the stay-level intervention feature table
--
-- PORT OF sql/mimiciv/v2_02b_interventions_mimiciv.sql. Output shape,
-- column names, types and the state/event split are identical. What
-- changes is where each intervention comes from and how an INTERVAL is
-- reconstructed, because eICU has no table with the interval semantics
-- MIMIC's `inputevents` and `mimiciv_derived.ventilation` provide.
--
-- ---------------------------------------------------------------------
-- WHAT ACTUALLY REACHES A MODEL, and why that governs this file
-- ---------------------------------------------------------------------
-- R/05_formula.R FORBIDDEN_VARS excludes `peak_intensity`, `dx_nee_peak`,
-- `first_hour`, `max_concurrent_agents` and `ever_active`. The intervention
-- columns that DO enter a layer-1 formula are exactly:
--
--   state  ->  {iv}__exposure_frac, {iv}__present_at_admission
--              (+ {iv}__n_agents for vasopressor / inotrope, which
--               `intensity_conditional` may replace with {iv}__lambda)
--   event  ->  {iv}__n_hours, {iv}__total_amount, {iv}__present_at_admission
--
-- THIS RETIRES THE HEADLINE eICU RISK. Reconciliation section F called
-- norepinephrine-equivalent dose "the worst problem, and it lands on the
-- paper's core claim", on the premise that peak_intensity is a fitted GAM
-- term. It is not, and has not been since `dx_nee_peak` took the `dx_`
-- prefix: roughly half of eICU norepinephrine exposure being charted in
-- uninterpretable ml/hr costs a DIAGNOSTIC column, not a model term. NEE
-- is still reconstructed below, three-tier and provenance-tagged, because
-- it is worth reporting and because the exposure-only ablation section F
-- proposes is then already the primary analysis rather than a fallback.
--
-- What DOES now carry the transport risk is `exposure_frac`, i.e. how many
-- of the 24 hours each intervention is judged active. That is entirely a
-- function of the interval reconstruction rules below. They are the part
-- of this file to argue about.
--
-- ---------------------------------------------------------------------
-- THREE INTERVAL RECONSTRUCTION RULES, all pre-specified here
-- ---------------------------------------------------------------------
--   R1  INFUSION CARRY-FORWARD (infusiondrug). A charted rate holds until
--       the next charted row for the same drug, or `infusion_carry_forward_min`
--       minutes, whichever is sooner. Identical to step 2a; see that file's
--       header for the argument and check A6 for the sensitivity analysis.
--
--   R2  BOLUS OCCUPANCY (medication). A bolus occupies the single hour
--       containing `drugstartoffset`. This mirrors MIMIC, where an
--       `inputevents` bolus has starttime and endtime minutes apart and
--       therefore occupies one hour. `drugstopoffset` is NOT used as an
--       interval end: audit 14b found 1,082,358 rows with
--       drugstopoffset <= drugstartoffset, so it is an order-expiry field,
--       not an administration end.
--
--   R3  TREATMENT ONSET-TO-WINDOW-END (treatment). A `treatment` row is a
--       care-plan entry, not an interval, and is charted per review rather
--       than per hour. An intervention recorded there is treated as active
--       from its first offset in the window to the end of the window.
--       This is the WEAKEST of the three and it BIASES exposure_frac
--       UPWARD relative to MIMIC. It is used only where nothing better
--       exists, and for invasive_vent it is the fallback, not the primary
--       route: `respiratorycare` carries real ventstart/ventend offsets
--       for 31,273 of the 36,927 stays in the union (audit 06a), and the
--       treatment route supplies only the remainder. Check B6 measures
--       what the fallback does to the exposure distribution.
--
-- ---------------------------------------------------------------------
-- VENTILATION ASCERTAINMENT - a known, quantified shortfall
-- ---------------------------------------------------------------------
-- MIMIC observes 41.4% invasive ventilation. The union rule below gives
-- 33.7% at eICU (audit 06a), with vented mortality 16.8% against 4.7%
-- unvented, so the ascertainment is at least behaving correctly. The 7.7
-- point gap is partly case mix (eICU includes many non-tertiary sites)
-- and partly under-ascertainment; audit 06b shows it is NOT uniform -
-- per-hospital prevalence runs from 0.000 to 0.844. That heterogeneity is
-- an input to the hospital inclusion decision (reconciliation section N)
-- and is not repaired here.
--
-- ---------------------------------------------------------------------
-- STATIC INTERVENTION LIST - a deliberate divergence from the MIMIC file
-- ---------------------------------------------------------------------
-- MIMIC builds `intervention_list` as SELECT DISTINCT FROM the hourly
-- table, so an intervention with zero rows would silently vanish from the
-- feature table. At MIMIC all 18 are present and it never fires. At eICU
-- an intervention could plausibly be empty, and validator check 2
-- ("modelled subset of present") would then fail with a message about the
-- model rather than about the extraction. The list below is STATIC and
-- matches config/config.yml interventions_extracted exactly, so every
-- stay gets all 18 rows and an absent intervention shows up as
-- ever_active = 0 everywhere - which is a readable fact, not a missing
-- row. The MIMIC file should adopt the same list for symmetry; that is a
-- one-CTE change and it alters no MIMIC number.
-- =====================================================================

DECLARE infusion_carry_forward_min INT64 DEFAULT 60;

CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_interventions_eicu` AS

WITH cohort AS (
  SELECT stay_id, weight_kg
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

-- Hour spine: every stay gets all 24 slots, so absence is explicit.
hours AS (
  SELECT
    c.stay_id,
    h            AS hour_bin,
    h * 60       AS hr_start_min,
    (h + 1) * 60 AS hr_end_min
  FROM cohort c, UNNEST(GENERATE_ARRAY(0, 23)) AS h
),

-- =====================================================================
-- INFUSION SOURCE (R1). One pass over `infusiondrug`, classified to an
-- intervention, with the canonical-unit flag that gates intensity.
--
-- CANONICAL UNITS. eICU charts the same drug under several units
-- (audit 08: Propofol appears as (ml/hr), (mcg/kg/min), (mg/kg/min) and
-- bare). MAX() across mixed units is meaningless, so `peak_intensity` is
-- computed ONLY from rows whose drugname carries the intervention's
-- canonical unit - chosen to match MIMIC's dominant `rateuom` for that
-- drug class. PRESENCE uses every unit variant, because exposure is
-- unit-free. A stay charted only in ml/hr therefore keeps its exposure
-- and carries a NULL peak_intensity, exactly as MIMIC's LEFT JOIN
-- semantics give a vasopressor stay with no computable NEE.
-- =====================================================================
inf_classified AS (
  SELECT
    i.patientunitstayid AS stay_id,
    i.infusionoffset    AS offset_min,
    LOWER(i.drugname)   AS dn,
    SAFE_CAST(i.drugrate AS FLOAT64) AS rate,
    CASE
      WHEN LOWER(i.drugname) LIKE '%midazolam%' OR LOWER(i.drugname) LIKE '%versed%'
        OR LOWER(i.drugname) LIKE '%lorazepam%' OR LOWER(i.drugname) LIKE '%ativan%'
        THEN 'sedation_benzo'
      WHEN LOWER(i.drugname) LIKE '%dexmedetomidine%'
        OR LOWER(i.drugname) LIKE '%dexmetetomidine%'
        OR LOWER(i.drugname) LIKE '%precedex%'
        THEN 'sedation_dexmed'
      WHEN LOWER(i.drugname) LIKE '%propofol%' OR LOWER(i.drugname) LIKE '%diprivan%'
        THEN 'sedation_propofol'
      WHEN LOWER(i.drugname) LIKE '%fentanyl%' OR LOWER(i.drugname) LIKE '%morphine%'
        OR LOWER(i.drugname) LIKE '%hydromorphone%' OR LOWER(i.drugname) LIKE '%dilaudid%'
        THEN 'opioid'
      WHEN LOWER(i.drugname) LIKE '%insulin%'
        THEN 'insulin'
      WHEN LOWER(i.drugname) LIKE '%cisatracurium%' OR LOWER(i.drugname) LIKE '%nimbex%'
        OR LOWER(i.drugname) LIKE '%rocuronium%'    OR LOWER(i.drugname) LIKE '%vecuronium%'
        OR LOWER(i.drugname) LIKE '%norcuron%'
        THEN 'paralytic'
      WHEN LOWER(i.drugname) LIKE '%dobutamine%' OR LOWER(i.drugname) LIKE '%milrinone%'
        THEN 'inotrope'
      WHEN LOWER(i.drugname) LIKE '%furosemide%' OR LOWER(i.drugname) LIKE '%lasix%'
        OR LOWER(i.drugname) LIKE '%bumetanide%'  OR LOWER(i.drugname) LIKE '%bumex%'
        THEN 'diuretic'
    END AS intervention
  FROM `physionet-data.eicu_crd.infusiondrug` i
  JOIN cohort c ON i.patientunitstayid = c.stay_id
  WHERE i.infusionoffset BETWEEN 0 AND 1439
    AND SAFE_CAST(i.drugrate AS FLOAT64) > 0
),

-- Epidural / PCA / combination bags are analgesia by a different route and
-- a different pharmacology; MIMIC's inputevents opioid set is systemic.
-- Removing them here keeps `opioid` the same construct at both sites.
inf_filtered AS (
  SELECT *
  FROM inf_classified
  WHERE intervention IS NOT NULL
    AND NOT (intervention = 'opioid'
             AND (dn LIKE '%epidural%' OR dn LIKE '%bupiv%' OR dn LIKE '%pca%'))
),

inf_points AS (
  SELECT
    stay_id, intervention, offset_min,
    MAX(rate) AS rate,
    -- Canonical unit per intervention. See header.
    MAX(CASE
      WHEN intervention = 'sedation_propofol' AND dn LIKE '%(mcg/kg/min)%' THEN rate
      WHEN intervention = 'sedation_benzo'    AND dn LIKE '%(mg/hr)%'      THEN rate
      WHEN intervention = 'sedation_dexmed'   AND dn LIKE '%(mcg/kg/hr)%'  THEN rate
      WHEN intervention = 'opioid'            AND dn LIKE '%(mcg/hr)%'     THEN rate
      WHEN intervention = 'insulin'           AND dn LIKE '%(units/hr)%'   THEN rate
      WHEN intervention = 'paralytic'         AND dn LIKE '%(mcg/kg/min)%' THEN rate
      WHEN intervention = 'inotrope'          AND dn LIKE '%(mcg/kg/min)%' THEN rate
      WHEN intervention = 'diuretic'          AND dn LIKE '%(mg/hr)%'      THEN rate
      ELSE NULL END)                          AS rate_canonical
  FROM inf_filtered
  GROUP BY 1, 2, 3
),

inf_intervals AS (
  SELECT
    stay_id, intervention, rate, rate_canonical,
    offset_min AS start_min,
    LEAST(
      COALESCE(
        LEAD(offset_min) OVER (PARTITION BY stay_id, intervention ORDER BY offset_min),
        offset_min + infusion_carry_forward_min),
      offset_min + infusion_carry_forward_min) AS end_min
  FROM inf_points
),

inf_hourly AS (
  SELECT
    h.stay_id,
    v.intervention,
    h.hour_bin,
    MAX(v.rate_canonical) AS intensity
  FROM hours h
  JOIN inf_intervals v
    ON h.stay_id  = v.stay_id
   AND v.start_min <  h.hr_end_min
   AND v.end_min   >  h.hr_start_min
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- BOLUS SOURCE (R2). `medication`, for the drug classes MIMIC keeps both
-- routes for. MIMIC's comment is explicit: "benzo and paralytic boluses
-- drive GCS masking and must not be dropped."
--
-- NO ROUTE FILTER for the masking drugs. Route filtering exists to
-- exclude the routes MIMIC's `inputevents` structurally lacks - oral
-- diuretic and subcutaneous insulin - and there is no oral rocuronium.
-- Filtering here on routeadmin would instead drop the junk-route rows
-- audit 14a shows ('ZPYXVEND', 'DEVICE', 'MISC', '.ROUTE' account for
-- 2,214 rocuronium rows), and losing those weakens GCS masking for no
-- principled reason.
--
-- INSULIN IS ABSENT from this CTE, deliberately: MIMIC restricts insulin
-- to infusions only (`NOT (insulin AND rate IS NULL)`), because
-- subcutaneous basal and sliding-scale insulin is routine glycaemic
-- maintenance, not an ICU intervention. eICU's SubQ insulin is large
-- (INSULIN-LISPRO SubQ 9,105 stays) and would swamp the construct.
-- =====================================================================
med_bolus AS (
  SELECT
    m.patientunitstayid AS stay_id,
    CAST(FLOOR(m.drugstartoffset / 60) AS INT64) AS hour_bin,
    CASE
      WHEN LOWER(m.drugname) LIKE '%midazolam%' OR LOWER(m.drugname) LIKE '%versed%'
        OR LOWER(m.drugname) LIKE '%lorazepam%' OR LOWER(m.drugname) LIKE '%ativan%'
        THEN 'sedation_benzo'
      WHEN LOWER(m.drugname) LIKE '%dexmedetomidine%' OR LOWER(m.drugname) LIKE '%precedex%'
        THEN 'sedation_dexmed'
      WHEN LOWER(m.drugname) LIKE '%propofol%' OR LOWER(m.drugname) LIKE '%diprivan%'
        THEN 'sedation_propofol'
      WHEN LOWER(m.drugname) LIKE '%fentanyl%' OR LOWER(m.drugname) LIKE '%morphine%'
        OR LOWER(m.drugname) LIKE '%hydromorphone%' OR LOWER(m.drugname) LIKE '%dilaudid%'
        THEN 'opioid'
      WHEN LOWER(m.drugname) LIKE '%cisatracurium%' OR LOWER(m.drugname) LIKE '%nimbex%'
        OR LOWER(m.drugname) LIKE '%rocuronium%'    OR LOWER(m.drugname) LIKE '%vecuronium%'
        OR LOWER(m.drugname) LIKE '%norcuron%'      OR LOWER(m.drugname) LIKE '%succinylcholine%'
        THEN 'paralytic'
    END AS intervention
  FROM `physionet-data.eicu_crd.medication` m
  JOIN cohort c ON m.patientunitstayid = c.stay_id
  WHERE m.drugstartoffset BETWEEN 0 AND 1439
    AND COALESCE(m.drugordercancelled, 'No') != 'Yes'   -- 205,273 cancelled orders
    AND NOT (LOWER(m.drugname) LIKE '%epidural%' OR LOWER(m.drugname) LIKE '%bupiv%')
),

bolus_hourly AS (
  SELECT stay_id, intervention, hour_bin, CAST(NULL AS FLOAT64) AS intensity
  FROM med_bolus
  WHERE intervention IS NOT NULL
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- VASOPRESSOR (STATE + diagnostic intensity)
--
-- PRESENCE from the agents table, exactly as at MIMIC: it is the
-- authoritative label-matched source and it is already hour-resolved.
--
-- INTENSITY is norepinephrine-equivalent dose in mcg/kg/min, three-tier
-- and DIAGNOSTIC ONLY (`dx_nee_peak`, refused by the formula builder):
--   tier 1  drugname carries (mcg/min) or (mcg/kg/min) -> direct
--   tier 2  drugname carries (ml/hr) AND a concentration is recoverable
--           from this stay's `medication` rows -> rate * mcg_per_mL / 60
--   tier 3  neither -> NULL. The presence row survives; only the dose is
--           absent, which is what a LEFT JOIN gives at MIMIC too.
-- Audit 05b: tier 2 recovers 2,629 of the 4,838 ml/hr-only stays (54.3%).
-- Audit 05c shows the recovered concentrations are BIMODAL at 32 and
-- 16 mcg/mL, so a single pre-specified default for tier 3 is NOT
-- defensible and none is applied - that was the pre-registered criterion
-- in discovery Q5c and it came back negative.
--
-- NEE weights follow mimiciv_derived.norepinephrine_equivalent_dose:
--   norepinephrine 1, epinephrine 1, phenylephrine /10, dopamine /100,
--   vasopressin (units/min) * 2.5.
-- =====================================================================
nee_conc AS (
  SELECT
    patientunitstayid AS stay_id,
    -- One concentration per stay: the modal recovered value, broken by
    -- taking the most frequently ordered. Deterministic and reported.
    ARRAY_AGG(mcg_per_ml ORDER BY n_orders DESC, mcg_per_ml LIMIT 1)[OFFSET(0)] AS mcg_per_ml
  FROM (
    SELECT
      patientunitstayid,
      SAFE_DIVIDE(
        SAFE_CAST(REGEXP_EXTRACT(LOWER(drugname), r'([0-9]+\.?[0-9]*)\s*mg') AS FLOAT64) * 1000,
        SAFE_CAST(REGEXP_EXTRACT(LOWER(drugname), r'([0-9]+\.?[0-9]*)\s*m[l]') AS FLOAT64)
      ) AS mcg_per_ml,
      COUNT(*) AS n_orders
    FROM `physionet-data.eicu_crd.medication`
    WHERE (LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%levophed%')
      AND COALESCE(drugordercancelled, 'No') != 'Yes'
    GROUP BY 1, 2
  )
  WHERE mcg_per_ml IS NOT NULL AND mcg_per_ml > 0
  GROUP BY 1
),

nee_points AS (
  SELECT
    i.patientunitstayid AS stay_id,
    i.infusionoffset    AS offset_min,
    LOWER(i.drugname)   AS dn,
    SAFE_CAST(i.drugrate AS FLOAT64) AS rate,
    c.weight_kg,
    nc.mcg_per_ml
  FROM `physionet-data.eicu_crd.infusiondrug` i
  JOIN cohort c ON i.patientunitstayid = c.stay_id
  LEFT JOIN nee_conc nc ON i.patientunitstayid = nc.stay_id
  WHERE i.infusionoffset BETWEEN 0 AND 1439
    AND SAFE_CAST(i.drugrate AS FLOAT64) > 0
    AND (LOWER(i.drugname) LIKE '%norepinephrine%' OR LOWER(i.drugname) LIKE '%levophed%'
      OR LOWER(i.drugname) LIKE '%epinephrine%'    OR LOWER(i.drugname) LIKE '%adrenalin%'
      OR LOWER(i.drugname) LIKE '%phenylephrine%'  OR LOWER(i.drugname) LIKE '%neosynephrine%'
      OR LOWER(i.drugname) LIKE '%vasopressin%'
      OR LOWER(i.drugname) LIKE '%dopamine%')
),

nee_dose AS (
  SELECT
    stay_id,
    CAST(FLOOR(offset_min / 60) AS INT64) AS hour_bin,
    SUM(
      -- to mcg/kg/min first, then weight the molecule
      (CASE
         WHEN dn LIKE '%(mcg/kg/min)%' THEN rate
         WHEN dn LIKE '%(mcg/min)%'    THEN SAFE_DIVIDE(rate, weight_kg)
         WHEN dn LIKE '%(units/min)%'  THEN rate            -- vasopressin, see below
         WHEN dn LIKE '%(ml/hr)%'      THEN SAFE_DIVIDE(rate * mcg_per_ml / 60.0, weight_kg)
         ELSE NULL
       END)
      *
      (CASE
         WHEN dn LIKE '%norepinephrine%' OR dn LIKE '%levophed%'      THEN 1.0
         WHEN dn LIKE '%epinephrine%'    OR dn LIKE '%adrenalin%'     THEN 1.0
         WHEN dn LIKE '%phenylephrine%'  OR dn LIKE '%neosynephrine%' THEN 0.1
         WHEN dn LIKE '%dopamine%'                                    THEN 0.01
         WHEN dn LIKE '%vasopressin%'                                 THEN 2.5
         ELSE NULL
       END)
    ) AS nee
  FROM nee_points
  GROUP BY 1, 2
),

vaso AS (
  SELECT
    p.stay_id,
    'vasopressor' AS intervention,
    p.hour_bin,
    MAX(n.nee) AS intensity
  FROM (
    -- PRESENCE: authoritative, from the agents table.
    SELECT DISTINCT stay_id, hour_bin
    FROM `eicu-ext.eicu_ext_data.v2_agents_eicu`
    WHERE intervention = 'vasopressor'
  ) p
  -- INTENSITY: LEFT JOIN, so a stay with no computable NEE keeps its
  -- presence row and simply carries a NULL dose.
  LEFT JOIN nee_dose n
    ON p.stay_id = n.stay_id AND p.hour_bin = n.hour_bin
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- VENTILATION (STATE, ordinal, R3 for the treatment routes)
-- Invasive kept distinct from NIV; collapsing them loses the clinically
-- decisive difference. Precedence per hour: invasive > NIV > supplemental,
-- matching the MIMIC CASE.
-- =====================================================================
-- REWRITTEN 2026-09-04 after audit round 6. This is the FOURTH version of
-- the eICU ventilation rule and the first one whose duration has been
-- validated against a known truth. The three earlier versions and why each
-- failed are recorded below, because the sequence is the argument.
--
-- v1 (original). Mapped every non-increasing `ventendoffset` to 1440.
--   MEASURED: that branch is taken on 505,467 of 582,356 rows (86.8%), and
--   it produced `invasive_vent__exposure_frac` = 0.978 mean with every
--   decile from the tenth up at exactly 1.0, against MIMIC's 0.661 with a
--   genuine spread. The rule was not recovering unknown end times; it was
--   overwriting known ones.
--
-- v2. Dropped rows with unusable end times and let those stays fall to the
--   `treatment` route. Exposure improved only to 0.838 and PREVALENCE
--   COLLAPSED from 32.2% to 8.1%, with mortality among the apparently
--   non-ventilated rising from 4.7% to 8.6% — thousands of genuinely
--   ventilated, sicker patients reclassified. Never trade ascertainment
--   for duration.
--
-- v3 (previous). Separated the two questions: PRESENCE from
--   `ventstartoffset`, which is reliable, and DURATION from a valid
--   `ventendoffset` interval where one exists. Presence was fixed —
--   prevalence 34.0% with the cleanest mortality split of the three,
--   17.7% against 5.0%. Duration was not: 84.1% of the 19,397 stays with
--   at least one valid interval STILL came out at an exposure of exactly
--   1.0, because the intervals tile the window rather than describing
--   episodes. Audit round 4 concluded from this that `respiratorycare`
--   does not encode extubation at all.
--
-- v4 (this version). Round 4's conclusion was drawn from one column and
-- was too strong, and round 6 found two things it had missed.
--
--   FIRST, A GOLD STANDARD. `priorventendoffset` describes a COMPLETED
--   prior episode. 6,880 cohort stays carry a valid one inside the window,
--   its values spread across the whole day (deciles 76 to 1,220 minutes),
--   and only 2.6% sit at the window end. That is an event time, not a
--   fallback showing through, so for those stays the true end is known.
--   Coverage of 21.2% is far too low to BE the duration source, but it is
--   enough to score every candidate against a truth. (V4, V4b: only 1.8%
--   of them are re-intubated, so it is a genuine stopping time.)
--
--   SECOND, THE RIGHT INSTRUMENT. Audit round 4 rejected ventilator
--   settings after counting distinct PEEP-charted hours (median 6) and
--   comparing that with MIMIC's median DURATION of 17 hours. A count of
--   charted hours is a coverage statistic: respiratory therapy charts on
--   rounds, and MEASURED, barely half the hours inside a stay's charting
--   span carry a value. The matching quantity is the SPAN, and at the
--   median these stays have 6 charted hours across a 19-hour span.
--
--   SCORED AGAINST THE GOLD STANDARD (V7, V7b), in hours of error:
--     onset-to-window-end, i.e. v1/v3's fallback   +13.61   r = —
--     last respiratory-care review                  +4.33   r = 0.358
--     strict ventilator-charting span               +3.57   r = 0.850
--   The variable this project shipped until today ran 13.6 hours long on
--   the only stays where the truth is knowable.
--
-- THE RULE NOW. Presence is decided exactly as in v3 and is UNTOUCHED, so
-- prevalence and the mortality split cannot move. Only the END changes,
-- and it is taken from the first of these that exists:
--
--   1. the last minute of strict ventilator charting   73.9% of stays
--   2. the last in-window respiratory care review       7.7%
--   3. onset to window end, as before                  18.4%
--
-- MEASURED CONSEQUENCE (V9, simulated from source): mean exposure 0.744
-- against 0.979, and the share at exactly 1.0 falls from 90.5% to 18.6%,
-- of which 98% is rung 3. Deciles run 0.168, 0.422, 0.683, 0.842, 0.903,
-- 0.942, 0.973, 0.996. MIMIC's comparator is 0.661 with 25.0% at 1.0, and
-- was not used to choose anything here — see the label list below.
--
-- WHAT IS NOT USED, AND WHY.
--
--   `ventendoffset` no longer contributes at all. It is 86.8%
--   order-artefact, and where valid it tiles the window. V10 settles the
--   last question about it: of the 5,979 stays that fall to rung 3, ZERO
--   have a valid in-window `ventendoffset` interval, so adding it as a
--   fourth rung would move no stay. `vent_respcare_valid` is deleted
--   rather than left unused.
--
--   `priorventendoffset` is NOT a rung, although it is the most accurate
--   value where it exists. A stay has one precisely when its ventilation
--   ENDED inside the window, so the gold stays are the short-ventilation
--   patients — mean exposure 0.4328 against 0.7440 pooled. Using it would
--   put unbiased values at the low end of the covariate's range and
--   +3.6-hour ones at the high end, which is a range-dependent distortion
--   rather than a shift. The covariate is read by a smooth fitted at
--   MIMIC, so what has to survive transport is its SHAPE: a uniform bias
--   is a translation that a monotone smooth absorbs while preserving
--   patient ordering, and a distortion that stretches one end of the
--   range is not. One estimator for every stay, one meaning for the
--   covariate. The gold-standard variant is reported as a sensitivity and
--   differs by 0.032 on the mean (audit round 6, V9, hierarchy H1).
--
-- =====================================================================
-- Every stay with a plausible ventilation START. THIS DECIDES PRESENCE
-- and is unchanged from v3: no stay is lost because its end times are
-- unusable, which is the v2 failure.
vent_respcare_any AS (
  SELECT
    r.patientunitstayid                 AS stay_id,
    MIN(GREATEST(r.ventstartoffset, 0)) AS start_min
  FROM `physionet-data.eicu_crd.respiratorycare` r
  JOIN cohort c ON r.patientunitstayid = c.stay_id
  WHERE r.ventstartoffset IS NOT NULL
    AND r.ventstartoffset BETWEEN -1440 AND 1440
  GROUP BY 1
),

-- RUNG 1. The last minute at which a ventilator was being charted.
--
-- THE LABEL LIST IS THE WHOLE DESIGN AND IT WAS CHOSEN ON MORTALITY.
-- `respiratorycharting` has three categories: `respFlowPtVentData`
-- (measured from the ventilator), `respFlowSettings` (ventilator and
-- oxygen settings) and `respFlowCareData` (a 60-label hospital-specific
-- tail). Kept here is everything that is charted ONLY when a ventilator
-- is delivering and measuring breaths.
--
-- DELIBERATELY EXCLUDED, and each exclusion costs coverage on purpose:
--   `FiO2` (37,564 stays) and `LPM O2` (19,708) — charted for anyone on
--     supplemental oxygen.
--   `SaO2` (16,694) — pulse oximetry. `Humidifier Temp` (5,583) — also
--     high-flow nasal cannula.
--   `CPAP` and `PEEP/CPAP` — non-invasive by name, and this project
--     models `noninvasive_vent` separately.
--   `RR (patient)` (26,749) and `Total RR` (18,170) — the tempting ones,
--     because together they lift coverage from 73.9% to 80.5%. MEASURED
--     (V6): 18,135 stays carry that charting without being flagged
--     ventilated, only 1,939 of them are flagged non-invasive, and their
--     mortality is 0.0529 against 0.0481 for stays with neither flag nor
--     charting. They are not ventilated patients we are missing; somebody
--     charted a respiratory rate for them. The strict set's 1,785
--     unflagged stays sit at 0.1048, between the never-ventilated rate
--     and the ventilated 0.1765, which is what a mixture of real missed
--     ventilation and non-invasive support should look like.
--   `PEEP` — adds 35 stays over that and changes nothing.
--
-- *** A TRAP, NAMED, BECAUSE IT WAS ONE NUMBER AWAY FROM SPRINGING. ***
-- The rejected wider set puts 0.2502 of ventilated stays at exactly 1.0
-- and MIMIC's share is 0.250. Had the label sets been compared on
-- resemblance to MIMIC instead of on the mortality of the stays they
-- recruit, the contaminated set would have won on a coincidence. The
-- acceptance criteria were fixed before any of these queries ran.
vent_setting_span AS (
  SELECT
    r.patientunitstayid    AS stay_id,
    MAX(r.respchartoffset) AS last_setting_min
  FROM `physionet-data.eicu_crd.respiratorycharting` r
  JOIN cohort c ON r.patientunitstayid = c.stay_id
  WHERE r.respchartoffset BETWEEN 0 AND 1439
    AND r.respcharttypecat IN ('respFlowSettings', 'respFlowPtVentData')
    AND r.respchartvaluelabel IN (
      -- measured from the ventilator
      'Exhaled MV', 'Exhaled TV (patient)', 'Exhaled TV (machine)',
      'Peak Insp. Pressure', 'Plateau Pressure', 'Mean Airway Pressure',
      'Compliance', 'A1: High Exhaled Vt',
      -- set on the ventilator
      'Vent Rate', 'Tidal Volume (set)', 'TV/kg IBW',
      'Pressure Support', 'Pressure Control', 'Peak Flow',
      'Flow Sensitivity', 'Pressure to Trigger PS')
  GROUP BY 1
),

-- RUNG 2. The last respiratory care review inside the window.
--
-- THE WINDOW RESTRICTION IS LOAD-BEARING. eICU stays run for days and
-- respiratory care keeps reviewing throughout: MEASURED, 73.4% of these
-- stays' `respiratorycare` rows lie past minute 1439 and 64.8% of the
-- stays have at least one. Without `respcarestatusoffset BETWEEN 0 AND
-- 1439` a review written on day three sets the end of day one. The audit
-- made exactly that mistake once (round 6, V3 against V3b) and it moved
-- the share at an exposure of 1.0 from 0.17% to 64%.
vent_last_review AS (
  SELECT
    r.patientunitstayid          AS stay_id,
    MAX(r.respcarestatusoffset)  AS last_review_min
  FROM `physionet-data.eicu_crd.respiratorycare` r
  JOIN cohort c ON r.patientunitstayid = c.stay_id
  WHERE r.ventstartoffset IS NOT NULL
    AND r.ventstartoffset BETWEEN -1440 AND 1440
    AND r.respcarestatusoffset BETWEEN 0 AND 1439
  GROUP BY 1
),

-- ONE INTERVAL PER STAY: start from presence, end from the hierarchy.
--
-- THE `start_min + 60` FLOOR IS A PRESENCE GUARD AND NOT A FUDGE. A
-- stay's last ventilator charting can precede its recorded ventilation
-- start, which would give a backwards interval, zero ventilated hours,
-- and `ever_active = 0` — silently converting a duration fix into an
-- ascertainment change, which is the v2 failure wearing a different
-- coat. Presence is decided by the start, so the onset hour is
-- ventilated by definition and the interval is never shorter than it.
vent_respcare AS (
  SELECT
    a.stay_id,
    a.start_min,
    LEAST(GREATEST(COALESCE(s.last_setting_min, lr.last_review_min, 1440),
                   a.start_min + 60), 1440) AS end_min
  FROM vent_respcare_any a
  LEFT JOIN vent_setting_span s  USING (stay_id)
  LEFT JOIN vent_last_review  lr USING (stay_id)
),

-- The second source (R3): stays `respiratorycare` does not mention at
-- all. They get the same duration hierarchy, minus rung 2, which cannot
-- exist for them. Applying one rule wherever the evidence exists is what
-- keeps the covariate meaning one thing.
vent_treatment AS (
  SELECT
    t.patientunitstayid                 AS stay_id,
    MIN(GREATEST(t.treatmentoffset, 0)) AS start_min
  FROM `physionet-data.eicu_crd.treatment` t
  JOIN cohort c ON t.patientunitstayid = c.stay_id
  WHERE t.treatmentoffset BETWEEN 0 AND 1439
    AND LOWER(t.treatmentstring) LIKE '%mechanical ventilation%'
    AND LOWER(t.treatmentstring) NOT LIKE '%non-invasive%'
  GROUP BY 1
),

vent_treatment_iv AS (
  SELECT
    t.stay_id,
    t.start_min,
    LEAST(GREATEST(COALESCE(s.last_setting_min, 1440),
                   t.start_min + 60), 1440) AS end_min
  FROM vent_treatment t
  LEFT JOIN vent_setting_span s USING (stay_id)
),

-- WHAT TO CHECK AFTER RE-RUNNING, and the first two are vetoes.
--
--   PREVALENCE must stay at 0.340 and non-ventilated mortality at 0.050.
--     Presence logic is untouched, so any movement means the duration
--     change leaked into ascertainment and the run is wrong.
--   `ever_active` must remain 1 for every stay it was 1 for. The
--     `start_min + 60` floor guarantees it; check rather than assume.
--   MEAN `invasive_vent__exposure_frac` should land near 0.74 and the
--     share at exactly 1.0 near 0.19, against 0.979 and 0.905 today.
--
-- EXPECT SMALL DIFFERENCES FROM V9's SIMULATION RATHER THAN NONE. V9
-- computed ends straight from the source tables; this builds `n_hours`
-- as COUNT(DISTINCT hour_bin) over the hour grid, so binning rounds.
-- V9 also used minute 0 as the start for treatment-route stays where
-- this uses `MIN(treatmentoffset)`, and the presence floor above gives a
-- minimum of one hour where V9 allowed zero. A few points of difference
-- are expected; a large one means something else changed.
vent_invasive AS (
  SELECT stay_id, start_min, end_min FROM vent_respcare
  UNION ALL
  SELECT t.stay_id, t.start_min, t.end_min
  FROM vent_treatment_iv t
  LEFT JOIN (SELECT DISTINCT stay_id FROM vent_respcare) r USING (stay_id)
  WHERE r.stay_id IS NULL
),

vent_niv AS (
  SELECT
    t.patientunitstayid AS stay_id,
    MIN(GREATEST(t.treatmentoffset, 0)) AS start_min,
    1440 AS end_min
  FROM `physionet-data.eicu_crd.treatment` t
  JOIN cohort c ON t.patientunitstayid = c.stay_id
  WHERE t.treatmentoffset BETWEEN 0 AND 1439
    AND (LOWER(t.treatmentstring) LIKE '%non-invasive ventilation%'
      OR LOWER(t.treatmentstring) LIKE '%cpap/peep therapy%')
  GROUP BY 1
),

vent_o2 AS (
  SELECT
    t.patientunitstayid AS stay_id,
    MIN(GREATEST(t.treatmentoffset, 0)) AS start_min,
    1440 AS end_min
  FROM `physionet-data.eicu_crd.treatment` t
  JOIN cohort c ON t.patientunitstayid = c.stay_id
  WHERE t.treatmentoffset BETWEEN 0 AND 1439
    AND LOWER(t.treatmentstring) LIKE '%oxygen therapy%'
  GROUP BY 1
),

vent_flags AS (
  SELECT
    h.stay_id,
    h.hour_bin,
    MAX(IF(iv.stay_id IS NOT NULL, 1, 0)) AS f_inv,
    MAX(IF(nv.stay_id IS NOT NULL, 1, 0)) AS f_niv,
    MAX(IF(o2.stay_id IS NOT NULL, 1, 0)) AS f_o2
  FROM hours h
  LEFT JOIN vent_invasive iv
    ON h.stay_id = iv.stay_id AND iv.start_min < h.hr_end_min AND iv.end_min > h.hr_start_min
  LEFT JOIN vent_niv nv
    ON h.stay_id = nv.stay_id AND nv.start_min < h.hr_end_min AND nv.end_min > h.hr_start_min
  LEFT JOIN vent_o2 o2
    ON h.stay_id = o2.stay_id AND o2.start_min < h.hr_end_min AND o2.end_min > h.hr_start_min
  GROUP BY 1, 2
),

vent AS (
  SELECT
    stay_id,
    CASE WHEN f_inv = 1 THEN 'invasive_vent'
         WHEN f_niv = 1 THEN 'noninvasive_vent'
         ELSE 'supplemental_o2' END AS intervention,
    hour_bin,
    CAST(NULL AS FLOAT64)           AS intensity
  FROM vent_flags
  WHERE f_inv = 1 OR f_niv = 1 OR f_o2 = 1
),

-- =====================================================================
-- VENTILATOR SETTINGS: FiO2 and PEEP as separate intensity channels.
-- FiO2 stays on the INTERVENTION side only. Do not also build P/F ratio -
-- that would put FiO2 on both sides of the oxygenation pair.
--
-- UNIT RULE (audit 17b): eICU FiO2 is overwhelmingly percent-form
-- (3,066,451 rows vs 7,889 fraction) and MIXED WITHIN A SINGLE LABEL, so
-- a per-value rule is required: value <= 1.0 -> x100. The ambiguity at
-- exactly 1.0 is harmless, since both readings mean 100%. Emitted on the
-- PERCENT scale to match mimiciv_derived.ventilator_setting.fio2.
-- (canonical_variable_spec.md section 4 says "fraction (0-1)"; the spec
-- doc is wrong and the MIMIC derived table is right, same class of error
-- as `map` vs `mbp`. This only ever reaches `peak_intensity`, which no
-- formula may carry, so nothing is fitted on either scale.)
-- =====================================================================
vent_settings AS (
  SELECT
    r.patientunitstayid AS stay_id,
    'fio2'              AS intervention,
    CAST(FLOOR(r.respchartoffset / 60) AS INT64) AS hour_bin,
    MAX(CASE WHEN SAFE_CAST(r.respchartvalue AS FLOAT64) <= 1.0
             THEN SAFE_CAST(r.respchartvalue AS FLOAT64) * 100.0
             ELSE SAFE_CAST(r.respchartvalue AS FLOAT64) END) AS intensity
  FROM `physionet-data.eicu_crd.respiratorycharting` r
  JOIN cohort c ON r.patientunitstayid = c.stay_id
  WHERE r.respchartoffset BETWEEN 0 AND 1439
    AND r.respchartvaluelabel IN ('FiO2', 'FIO2 (%)',
                                  'Set Fraction of Inspired Oxygen (FIO2)')
    AND SAFE_CAST(r.respchartvalue AS FLOAT64) BETWEEN 0.21 AND 100
  GROUP BY 1, 2, 3

  UNION ALL

  SELECT
    r.patientunitstayid,
    'peep',
    CAST(FLOOR(r.respchartoffset / 60) AS INT64),
    MAX(SAFE_CAST(r.respchartvalue AS FLOAT64))
  FROM `physionet-data.eicu_crd.respiratorycharting` r
  JOIN cohort c ON r.patientunitstayid = c.stay_id
  WHERE r.respchartoffset BETWEEN 0 AND 1439
    AND r.respchartvaluelabel IN ('PEEP', 'PEEP/CPAP')
    AND SAFE_CAST(r.respchartvalue AS FLOAT64) BETWEEN 0 AND 50
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- RRT (STATE, R3). Critical for step 5: this is the mask that INVALIDATES
-- creatinine and urine output, and the mask mode is 'from_onset', so the
-- onset-to-window-end rule is EXACTLY right here rather than an
-- approximation - only `first_hour` is consumed.
--
-- Whitelist is the dialysis subtree. ACCESS PROCEDURES ARE EXCLUDED:
-- 'insertion of venous catheter for hemodialysis' (2,292 + 1,087 + 775 +
-- 537 + 176 + 148 rows) and 'arteriovenous shunt' / 'dialysis access
-- surgery' record that a line was placed, not that blood was cleared. A
-- creatinine measured after line insertion but before the first session
-- is still valid, and masking it would delete real evidence.
-- =====================================================================
rrt_onset AS (
  SELECT
    t.patientunitstayid AS stay_id,
    MIN(GREATEST(t.treatmentoffset, 0)) AS start_min
  FROM `physionet-data.eicu_crd.treatment` t
  JOIN cohort c ON t.patientunitstayid = c.stay_id
  WHERE t.treatmentoffset BETWEEN 0 AND 1439
    AND LOWER(t.treatmentstring) LIKE '%dialysis%'
    AND LOWER(t.treatmentstring) NOT LIKE '%insertion of%'
    AND LOWER(t.treatmentstring) NOT LIKE '%arteriovenous shunt%'
    AND LOWER(t.treatmentstring) NOT LIKE '%access surgery%'
    AND LOWER(t.treatmentstring) NOT LIKE '%insertion of catheter%'
  GROUP BY 1
),

rrt AS (
  SELECT
    h.stay_id,
    'rrt'                 AS intervention,
    h.hour_bin,
    CAST(NULL AS FLOAT64) AS intensity
  FROM hours h
  JOIN rrt_onset o
    ON h.stay_id = o.stay_id
   AND o.start_min < h.hr_end_min
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- STATE infusions: sedation (split benzo / dexmed / propofol), opioid,
-- insulin, paralytic, inotrope. No dose equivalence attempted for
-- sedation - exposure only, unchanged from the MIMIC design.
-- =====================================================================
infusions AS (
  SELECT stay_id, intervention, hour_bin, MAX(intensity) AS intensity
  FROM (
    SELECT stay_id, intervention, hour_bin, intensity FROM inf_hourly
    WHERE intervention IN ('sedation_benzo', 'sedation_dexmed', 'sedation_propofol',
                           'opioid', 'insulin', 'paralytic', 'inotrope')
    UNION ALL
    SELECT stay_id, intervention, hour_bin, intensity FROM bolus_hourly
  )
  GROUP BY 1, 2, 3
),

-- =====================================================================
-- EVENT-BASED interventions (discrete). intensity = amount given in that
-- hour; exposure downstream is a COUNT over the window, not a fraction.
--
-- DIURETIC comes mainly from `medication`, which is the correct source
-- for an event-shaped intervention: audit 14a shows furosemide is
-- overwhelmingly a bolus drug (LASIX IV 6,205 stays, FUROSEMIDE IV 4,194,
-- IV push 3,708) while `infusiondrug` sees 701 + 542. Both routes are
-- included and both contribute milligrams delivered in the hour:
--   medication -> parsed leading number from `dosage`
--   infusiondrug (mg/hr) -> the rate, i.e. mg over that hour
-- ROUTE IS FILTERED HERE and only here, to IV forms, because oral
-- furosemide (FUROSEMIDE 40 MG PO TABS, 3,834 stays) is maintenance
-- therapy and MIMIC's `inputevents` structurally cannot contain it.
--
-- NO POTENCY EQUIVALENCE between furosemide and bumetanide, matching
-- MIMIC exactly, which sums `ie.amount` across both. 1 mg bumetanide is
-- roughly 40 mg furosemide, so `total_amount` is a mixture; it is the
-- SAME mixture at both sites, which is what the comparison needs. State
-- it as a limitation of the term, not of the transport.
--
-- TRANSFUSIONS come from `intakeoutput`, not `treatment`: audit 15b vs
-- 15a is 17,575 vs 1,733 stays for pRBC and 3,668 vs 602 for platelets,
-- and only intakeoutput carries the VOLUME that `total_amount` needs.
-- Cellpaths are an exact whitelist. Audit 15a's sweep also caught
-- antiplatelet agents (aspirin, clopidogrel), which are not transfusions
-- and are matched by nothing below.
-- =====================================================================
diuretic_events AS (
  SELECT
    m.patientunitstayid AS stay_id,
    'diuretic'          AS intervention,
    CAST(FLOOR(m.drugstartoffset / 60) AS INT64) AS hour_bin,
    SUM(SAFE_CAST(REGEXP_EXTRACT(m.dosage, r'([0-9]+\.?[0-9]*)') AS FLOAT64)) AS intensity
  FROM `physionet-data.eicu_crd.medication` m
  JOIN cohort c ON m.patientunitstayid = c.stay_id
  WHERE m.drugstartoffset BETWEEN 0 AND 1439
    AND COALESCE(m.drugordercancelled, 'No') != 'Yes'
    AND (LOWER(m.drugname) LIKE '%furosemide%' OR LOWER(m.drugname) LIKE '%lasix%'
      OR LOWER(m.drugname) LIKE '%bumetanide%' OR LOWER(m.drugname) LIKE '%bumex%')
    AND (UPPER(m.routeadmin) LIKE 'IV%' OR UPPER(m.routeadmin) LIKE 'INTRAVEN%')
  GROUP BY 1, 2, 3

  UNION ALL

  SELECT
    stay_id, 'diuretic', hour_bin, MAX(intensity)
  FROM inf_hourly
  WHERE intervention = 'diuretic'
  GROUP BY 1, 2, 3
),

transfusion_events AS (
  SELECT
    io.patientunitstayid AS stay_id,
    CASE
      WHEN io.cellpath IN (
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|pRBCs',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|Volume-Transfuse red blood cells',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|Volume (ml)-Transfuse - Leukoreduced Packed RBCs',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|PRBC',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Generic Intake (ml)|PRBC')
        THEN 'transfusion_prbc'
      WHEN io.cellpath IN (
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|Platelets',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|Volume-Transfuse platelet pheresis')
        THEN 'transfusion_platelet'
      WHEN io.cellpath IN (
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|FFP',
        'flowsheet|Flowsheet Cell Labels|I&O|Intake (ml)|Blood Products (ml)|Volume-Transfuse plasma')
        THEN 'transfusion_ffp'
    END AS intervention,
    CAST(FLOOR(io.intakeoutputoffset / 60) AS INT64) AS hour_bin,
    SUM(io.cellvaluenumeric) AS intensity
  FROM `physionet-data.eicu_crd.intakeoutput` io
  JOIN cohort c ON io.patientunitstayid = c.stay_id
  WHERE io.intakeoutputoffset BETWEEN 0 AND 1439
    AND io.cellvaluenumeric > 0
  GROUP BY 1, 2, 3
),

events AS (
  SELECT stay_id, intervention, hour_bin, SUM(intensity) AS intensity
  FROM (
    SELECT * FROM diuretic_events
    UNION ALL
    SELECT stay_id, intervention, hour_bin, intensity FROM transfusion_events
    WHERE intervention IS NOT NULL
  )
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
  UNION ALL
  SELECT stay_id, intervention, hour_bin, intensity, 'event' FROM events
)

SELECT
  c.stay_id,
  c.intervention,
  c.hour_bin,
  c.shape,
  c.intensity,
  1 AS active
FROM combined c
WHERE c.hour_bin BETWEEN 0 AND 23;


-- =====================================================================
-- 2b-ii. INTERVENTION FEATURES
--
-- Byte-identical in schema to the MIMIC file. The ONLY structural change
-- is the static `intervention_list` - see the header.
--
-- Two shapes, not interchangeable:
--   shape = 'state' -> exposure_frac and peak_intensity
--   shape = 'event' -> n_hours and total_amount; exposure_frac is
--                      MEANINGLESS for these
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_intervention_features_eicu` AS

WITH cohort AS (
  SELECT stay_id
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

iv AS (
  SELECT *
  FROM `eicu-ext.eicu_ext_data.v2_interventions_eicu`
),

-- STATIC. Must match config/config.yml interventions_extracted exactly.
-- An intervention with no rows at eICU then appears with ever_active = 0
-- on every stay, instead of disappearing and failing validator check 2
-- with a message about the model rather than the extraction.
intervention_list AS (
  SELECT 'diuretic'             AS intervention, 'event' AS shape UNION ALL
  SELECT 'fio2',                  'state' UNION ALL
  SELECT 'inotrope',              'state' UNION ALL
  SELECT 'insulin',               'state' UNION ALL
  SELECT 'invasive_vent',         'state' UNION ALL
  SELECT 'noninvasive_vent',      'state' UNION ALL
  SELECT 'opioid',                'state' UNION ALL
  SELECT 'paralytic',             'state' UNION ALL
  SELECT 'peep',                  'state' UNION ALL
  SELECT 'rrt',                   'state' UNION ALL
  SELECT 'sedation_benzo',        'state' UNION ALL
  SELECT 'sedation_dexmed',       'state' UNION ALL
  SELECT 'sedation_propofol',     'state' UNION ALL
  SELECT 'supplemental_o2',       'state' UNION ALL
  SELECT 'transfusion_ffp',       'event' UNION ALL
  SELECT 'transfusion_platelet',  'event' UNION ALL
  SELECT 'transfusion_prbc',      'event' UNION ALL
  SELECT 'vasopressor',           'state'
),

grid AS (
  SELECT c.stay_id, l.intervention, l.shape
  FROM cohort c CROSS JOIN intervention_list l
),

agent_ids AS (
  SELECT
    stay_id,
    intervention,
    COUNT(DISTINCT agent)                    AS n_agents,
    MAX(IF(agent = 'norepinephrine', 1, 0))  AS has_norepinephrine
  FROM `eicu-ext.eicu_ext_data.v2_agents_eicu`
  GROUP BY 1, 2
),

agent_concurrency AS (
  SELECT stay_id, intervention, MAX(n_concurrent) AS max_concurrent_agents
  FROM (
    SELECT stay_id, intervention, hour_bin, COUNT(DISTINCT agent) AS n_concurrent
    FROM `eicu-ext.eicu_ext_data.v2_agents_eicu`
    GROUP BY 1, 2, 3
  )
  GROUP BY 1, 2
),

rolled AS (
  SELECT
    g.stay_id,
    g.intervention,
    g.shape,
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
  'eicu'                                                     AS site,
  r.stay_id,
  r.intervention,
  r.shape,
  r.ever_active,

  IF(r.shape = 'state', r.n_hours / 24.0, NULL)              AS exposure_frac,
  IF(r.shape = 'event', r.n_hours,        NULL)              AS n_hours,
  IF(r.shape = 'event', COALESCE(r.total_amount, 0), NULL)   AS total_amount,

  IF(r.intervention IN ('vasopressor','inotrope'), NULL, r.raw_peak)
                                                             AS peak_intensity,
  IF(r.intervention = 'vasopressor', r.raw_peak, NULL)       AS dx_nee_peak,

  IF(r.intervention IN ('vasopressor','inotrope'),
     COALESCE(ai.n_agents, 0),             NULL)             AS n_agents,
  IF(r.intervention IN ('vasopressor','inotrope'),
     COALESCE(ac.max_concurrent_agents, 0), NULL)            AS max_concurrent_agents,
  IF(r.intervention = 'vasopressor',
     COALESCE(ai.has_norepinephrine, 0),   NULL)             AS has_norepinephrine,

  r.first_hour,
  r.present_at_admission
FROM rolled r
LEFT JOIN agent_ids         ai ON r.stay_id = ai.stay_id AND r.intervention = ai.intervention
LEFT JOIN agent_concurrency ac ON r.stay_id = ac.stay_id AND r.intervention = ac.intervention;


-- =====================================================================
-- CHECKS - the cross-site prevalence gate is the important one
-- =====================================================================
-- B1. PREVALENCE AGAINST MIMIC. This is the gate that decides whether the
--     label matching missed variants. MIMIC v2 observed: vasopressor
--     ~25-35%, invasive_vent 41.4%, rrt ~5-10%.
--     Audit 06a predicts invasive_vent 33.7% here; a figure far from that
--     means the union rule was mis-implemented.
-- SELECT intervention, AVG(ever_active) AS prevalence,
--        AVG(IF(ever_active = 1, exposure_frac, NULL)) AS mean_exposure_if_active,
--        AVG(IF(ever_active = 1, n_hours, NULL))       AS mean_n_hours_if_active
-- FROM `...v2_intervention_features_eicu` GROUP BY 1 ORDER BY 2 DESC;
--
-- B2. All 18 must appear, each on all 109,200 stays.
-- SELECT COUNT(DISTINCT intervention) AS n_iv, COUNT(*) / COUNT(DISTINCT stay_id) AS rows_per_stay
-- FROM `...v2_intervention_features_eicu`;
--
-- B3. SHAPE MAP must agree with config/config.yml intervention_shape.
--     The validator checks it on load; check it here first, it is cheaper.
-- SELECT intervention, ANY_VALUE(shape) FROM `...v2_intervention_features_eicu` GROUP BY 1;
--
-- B4. NEE TIER COVERAGE. Report this; it is the honest version of
--     reconciliation section F now that the column is diagnostic only.
-- SELECT COUNTIF(ever_active = 1)                          AS n_vaso_stays,
--        COUNTIF(ever_active = 1 AND dx_nee_peak IS NULL)   AS n_no_dose,
--        SAFE_DIVIDE(COUNTIF(ever_active = 1 AND dx_nee_peak IS NULL),
--                    COUNTIF(ever_active = 1))              AS frac_no_dose
-- FROM `...v2_intervention_features_eicu` WHERE intervention = 'vasopressor';
--
-- B5. AGENT POOL. Must be 5 / 2. See step 2a check A2 - this is the same
--     assertion viewed from the feature table.
-- SELECT intervention, MAX(n_agents) FROM `...v2_intervention_features_eicu`
-- WHERE intervention IN ('vasopressor','inotrope') GROUP BY 1;
--
-- B6. THE R3 BIAS, measured. Split invasive_vent exposure_frac by which
--     route supplied the stay. If the treatment-only stays pile up at
--     exposure_frac = 1.0 while the respiratorycare stays spread out, the
--     onset-to-window-end rule is manufacturing exposure and the fallback
--     route should be reported separately - or dropped, at the cost of
--     prevalence.
--     READ IT DIFFERENTLY AFTER THE 2026-09-04 REWRITE. The split that
--     matters is no longer respiratorycare against treatment but WHICH
--     RUNG supplied the end. Audit round 6 V9 predicts 73.9% of stays on
--     the charting span at mean exposure 0.709, 7.7% on the last review
--     at 0.480, and 18.4% still on onset-to-window-end at 0.994 with
--     99.2% at exactly 1.0. That last group IS the manufactured
--     exposure this check was written to find, it is now a fifth of the
--     ventilated stays rather than all of them, and V10 confirms it is
--     irreducible: not one of those 5,979 stays carries a valid in-window
--     `ventendoffset` interval either. Report its share beside any result
--     that uses this covariate.
-- SELECT r.stay_id IS NOT NULL AS has_respcare,
--        APPROX_QUANTILES(f.exposure_frac, 10) AS deciles, COUNT(*) AS n
-- FROM `...v2_intervention_features_eicu` f
-- LEFT JOIN (SELECT DISTINCT patientunitstayid AS stay_id
--            FROM `physionet-data.eicu_crd.respiratorycare`
--            WHERE ventstartoffset IS NOT NULL AND ventstartoffset <= 1440) r
--   USING (stay_id)
-- WHERE f.intervention = 'invasive_vent' AND f.ever_active = 1
-- GROUP BY 1;
--
-- B7. DIURETIC DOSE PARSE. What fraction of matched medication rows
--     yielded a numeric dosage? `total_amount` is a modelled term for
--     event shapes, so a low parse rate is a real problem, not a
--     diagnostic one.
-- SELECT COUNT(*) AS n_rows,
--        COUNTIF(SAFE_CAST(REGEXP_EXTRACT(dosage, r'([0-9]+\.?[0-9]*)') AS FLOAT64) IS NULL)
--          AS n_unparsed
-- FROM `physionet-data.eicu_crd.medication`
-- WHERE drugstartoffset BETWEEN 0 AND 1439
--   AND (LOWER(drugname) LIKE '%furosemide%' OR LOWER(drugname) LIKE '%lasix%');
--
-- B8. Cross-check invasive_vent against the GCS verbal coverage from step
--     04. eICU has no gcs_unable flag, so ventilated stays with a full
--     verbal score are the contamination reconciliation section L warns
--     about.
