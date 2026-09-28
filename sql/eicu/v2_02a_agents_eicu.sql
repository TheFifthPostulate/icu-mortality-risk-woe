-- #####################################################################
-- v2 STEP 2a: eICU vasopressor / inotrope AGENT IDENTITY
--
-- PORT OF sql/mimiciv/v2_02a_agents_mimiciv.sql. Output grain is
-- identical: one row per (stay_id, hour_bin, intervention, agent).
--
-- WHY THIS IS SEPARATE, unchanged from the MIMIC file: molecule identity
-- must survive to step 5 so `n_agents` can be counted. It is destroyed by
-- any dose-collapsing step.
--
-- THE BUG THIS CODE IS WRITTEN TO AVOID, also unchanged:
--   LIKE '%epinephrine%' matches NOREPINEPHRINE. The CASE tests
--   norepinephrine FIRST and the later branches are therefore unreachable
--   for it. DO NOT REORDER THESE BRANCHES. eICU makes this worse than
--   MIMIC, not better: the drugname vocabulary includes 'Levophed',
--   'levophed', 'NeoSynephrine', 'Neosynephrine' and case variants of
--   every molecule (audit 07, 104 distinct strings), so the match is on
--   LOWER(drugname) throughout.
--
-- THE eICU-SPECIFIC TRAP, and the reason check_agent_pool() exists.
--   canonical_variable_spec.md section 7 records it: counting DISTINCT
--   `drugname` rather than distinct MOLECULE inflates eICU's n_agents by
--   2-4x, because norepinephrine alone appears under at least 20 strings
--   that differ only in unit or premix concentration. That would arrive
--   in the model as a case-mix difference and would be invisible.
--   config/config.yml declares intervention_agent_pool: vasopressor 5,
--   inotrope 2, and R/04_features.R holds the extraction to it in BOTH
--   directions. This file is what must be right for that check to pass.
--
-- INFUSIONS ONLY, matching the MIMIC rule (`ie.rate IS NOT NULL AND > 0`).
--   `infusiondrug` is the infusion table; `medication` boluses are
--   excluded. Counting boluses on one site and not the other would
--   inflate eICU's n_agents and read as case mix. eICU's `drugrate` is
--   FREE TEXT ('OFF', 'Documentation undone', blank), so the rule is
--   SAFE_CAST(drugrate AS FLOAT64) > 0. Audit 07 counts the cost:
--   'Norepinephrine' with a non-numeric rate is 1,054 stays and
--   'Norepinephrine ()' with no unit is 1,339. They lose PRESENCE, not
--   just dose. This is the same rule MIMIC applies and the loss is
--   reported, not silently absorbed - see check A5.
--
-- ---------------------------------------------------------------------
-- INTERVAL RECONSTRUCTION - the one genuinely new mechanic
-- ---------------------------------------------------------------------
-- MIMIC's `inputevents` rows are INTERVALS (starttime, endtime), so hour
-- occupancy is exact. eICU's `infusiondrug` rows are POINT observations:
-- one row per charted rate, at `infusionoffset`. Taking only the hour
-- containing each row would undercount every hour in which nothing
-- happened to be charted, and `exposure_frac` is a modelled term.
--
-- RULE, one line and pre-specified: a charted rate holds until the next
-- charted row for the SAME molecule, or for `infusion_carry_forward_min`
-- minutes, whichever comes first. At eICU's roughly hourly infusion
-- charting this reproduces continuous coverage during an infusion and
-- stops within an hour of the last chart. It never bridges a gap longer
-- than the carry-forward, so a genuine stop costs at most one hour of
-- false exposure.
--
-- This is an APPROXIMATION and it is the largest single threat to the
-- comparability of `vasopressor__exposure_frac` between the two sites.
-- Check A6 measures it directly; do not skip it.
-- #####################################################################

DECLARE infusion_carry_forward_min INT64 DEFAULT 60;

CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_agents_eicu` AS

WITH cohort AS (
  SELECT stay_id
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

hours AS (
  SELECT
    c.stay_id,
    h                AS hour_bin,
    h * 60           AS hr_start_min,
    (h + 1) * 60     AS hr_end_min
  FROM cohort c, UNNEST(GENERATE_ARRAY(0, 23)) AS h
),

-- One row per charted rate, already reduced to a MOLECULE. Rows that fail
-- the numeric-rate test are dropped here, before any interval is built.
charted AS (
  SELECT
    i.patientunitstayid                        AS stay_id,
    i.infusionoffset                           AS offset_min,
    CASE
      WHEN LOWER(i.drugname) LIKE '%dobutamine%'
        OR LOWER(i.drugname) LIKE '%milrinone%'  THEN 'inotrope'
      ELSE 'vasopressor'
    END                                        AS intervention,
    CASE
      -- ORDER MATTERS. norepinephrine must be tested before epinephrine.
      WHEN LOWER(i.drugname) LIKE '%norepinephrine%'
        OR LOWER(i.drugname) LIKE '%levophed%'      THEN 'norepinephrine'
      WHEN LOWER(i.drugname) LIKE '%epinephrine%'
        OR LOWER(i.drugname) LIKE '%adrenalin%'     THEN 'epinephrine'
      WHEN LOWER(i.drugname) LIKE '%phenylephrine%'
        OR LOWER(i.drugname) LIKE '%neosynephrine%' THEN 'phenylephrine'
      WHEN LOWER(i.drugname) LIKE '%vasopressin%'   THEN 'vasopressin'
      WHEN LOWER(i.drugname) LIKE '%dopamine%'      THEN 'dopamine'
      WHEN LOWER(i.drugname) LIKE '%dobutamine%'    THEN 'dobutamine'
      WHEN LOWER(i.drugname) LIKE '%milrinone%'     THEN 'milrinone'
    END                                        AS agent
  FROM `physionet-data.eicu_crd.infusiondrug` i
  JOIN cohort c
    ON i.patientunitstayid = c.stay_id
  WHERE i.infusionoffset BETWEEN 0 AND 1439
    AND SAFE_CAST(i.drugrate AS FLOAT64) > 0     -- infusions only, never boluses
    AND (LOWER(i.drugname) LIKE '%norepinephrine%' OR LOWER(i.drugname) LIKE '%levophed%'
      OR LOWER(i.drugname) LIKE '%epinephrine%'    OR LOWER(i.drugname) LIKE '%adrenalin%'
      OR LOWER(i.drugname) LIKE '%phenylephrine%'  OR LOWER(i.drugname) LIKE '%neosynephrine%'
      OR LOWER(i.drugname) LIKE '%vasopressin%'
      OR LOWER(i.drugname) LIKE '%dopamine%'
      OR LOWER(i.drugname) LIKE '%dobutamine%'
      OR LOWER(i.drugname) LIKE '%milrinone%')
),

-- Collapse duplicate charting of the same molecule at the same minute
-- (different drugname strings for one bag) BEFORE building intervals,
-- otherwise LEAD produces zero-length spans.
distinct_points AS (
  SELECT DISTINCT stay_id, intervention, agent, offset_min
  FROM charted
  WHERE agent IS NOT NULL
),

intervals AS (
  SELECT
    stay_id,
    intervention,
    agent,
    offset_min                                              AS start_min,
    LEAST(
      COALESCE(
        LEAD(offset_min) OVER (PARTITION BY stay_id, agent ORDER BY offset_min),
        offset_min + infusion_carry_forward_min),
      offset_min + infusion_carry_forward_min)               AS end_min
  FROM distinct_points
)

SELECT
  h.stay_id,
  h.hour_bin,
  v.intervention,
  v.agent
FROM hours h
JOIN intervals v
  ON h.stay_id  = v.stay_id
 AND v.start_min <  h.hr_end_min
 AND v.end_min   >  h.hr_start_min
GROUP BY 1, 2, 3, 4;


-- ---------------------------------------------------------------------
-- CHECKS FOR STEP 2a - run before trusting anything downstream.
-- ---------------------------------------------------------------------
-- A1. Agent vocabulary actually matched. Any NULL agent means a drugname
--     slipped past the CASE while passing the WHERE - impossible unless
--     the two lists drift apart, so this must return zero.
-- SELECT COUNTIF(agent IS NULL) AS n_unclassified FROM `...v2_agents_eicu`;
--
-- A2. THE AGENT POOL. This is the check config/intervention_agent_pool
--     exists for. It must return EXACTLY 5 for vasopressor and 2 for
--     inotrope. Anything else and R/04_features.R check_agent_pool()
--     stops the run - which is the intended behaviour, not a bug to work
--     around by editing config.
-- SELECT intervention, COUNT(DISTINCT agent) AS pool_size
-- FROM `...v2_agents_eicu` GROUP BY 1;
--
--     MILRINONE IS THE ONE AT RISK. Audit 07 swept vasopressor names only
--     and does not report it. If milrinone returns zero stays, the
--     inotrope pool is 1 at eICU and 2 at MIMIC, `n_agents` is not
--     comparable for inotrope, and that must be reported rather than
--     patched. Run this before anything else:
-- SELECT COUNT(DISTINCT patientunitstayid) FROM `physionet-data.eicu_crd.infusiondrug`
-- WHERE LOWER(drugname) LIKE '%milrinone%';
--
-- A3. Norepinephrine must be the most common agent by a wide margin. If
--     epinephrine outnumbers it, the CASE branches got reordered.
-- SELECT agent, COUNT(DISTINCT stay_id) AS n_stays
-- FROM `...v2_agents_eicu` GROUP BY 1 ORDER BY n_stays DESC;
--
-- A4. Drugname audit - confirm nothing was silently dropped. Compare the
--     matched set against the full sweep in audit 07.
-- SELECT DISTINCT drugname FROM `physionet-data.eicu_crd.infusiondrug`
-- WHERE LOWER(drugname) LIKE '%pressin%' OR LOWER(drugname) LIKE '%amine%'
--    OR LOWER(drugname) LIKE '%phrine%'  OR LOWER(drugname) LIKE '%phed%';
--
-- A5. COST OF THE NUMERIC-RATE RULE. How many stays have a matching
--     drugname in the window but NO row with a parseable positive rate?
--     Audit 07 predicts roughly 2,400 for norepinephrine alone. These
--     stays lose PRESENCE, not merely dose, so the number belongs in the
--     limitations beside the ml/hr dose problem.
-- WITH any_row AS (
--   SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.infusiondrug`
--   WHERE infusionoffset BETWEEN 0 AND 1439
--     AND (LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%levophed%')),
-- kept AS (SELECT DISTINCT stay_id FROM `...v2_agents_eicu` WHERE agent = 'norepinephrine')
-- SELECT COUNT(*) AS n_any, COUNTIF(kept.stay_id IS NULL) AS n_lost
-- FROM any_row LEFT JOIN kept ON any_row.patientunitstayid = kept.stay_id;
--
-- A6. THE CARRY-FORWARD SENSITIVITY. Re-run this file at
--     infusion_carry_forward_min = 30, 60 and 120 and compare the
--     distribution of vasopressor exposure_frac from step 2b. If the
--     median moves materially between 60 and 120, exposure_frac is being
--     set by the reconstruction rule rather than by the treatment, and
--     the cross-site comparison of that term is not supportable at face
--     value. Report the sensitivity either way; it costs one re-run of
--     two small files.
--
-- A7. Prevalence gate against MIMIC. Vasopressor prevalence at MIMIC is
--     ~25-35%. A figure far below that at eICU is under-ascertainment
--     from A5 plus the carry-forward, not case mix.
-- SELECT COUNT(DISTINCT stay_id) FROM `...v2_agents_eicu`
-- WHERE intervention = 'vasopressor';
