-- #####################################################################
-- PART A — NEW FILE: v2_02b_agents_mimiciv.sql
--
-- WHY THIS IS SEPARATE.
--   The existing `vaso` CTE in 02 reads
--   `mimiciv_3_1_derived.norepinephrine_equivalent_dose`, which is already
--   agent-collapsed: it returns one NEE per stay-time interval with no
--   record of which molecules produced it. Molecule identity is destroyed
--   before it ever reaches 05, so n_agents cannot be recovered downstream.
--   This script goes back to `inputevents` for identity only, at a
--   different grain — one row per (stay_id, hour_bin, intervention, agent).
--   Keeping it out of the main interventions table avoids fanning that
--   table out and leaves its schema untouched.
--
-- THE BUG THIS CODE IS WRITTEN TO AVOID.
--   LIKE '%epinephrine%' matches NOREPINEPHRINE as well as epinephrine.
--   A naive CASE would silently label every norepinephrine row as
--   epinephrine (or collapse the two), and n_agents would be wrong in a
--   way no prevalence check catches. The CASE below tests norepinephrine
--   FIRST and the later branches are therefore unreachable for it. Do not
--   reorder these branches.
--
-- INFUSIONS ONLY.
--   `ie.rate IS NOT NULL` selects rate-based (infusion) administrations and
--   excludes boluses. This must match the eICU side, where `infusiondrug`
--   is used and `medication` boluses are excluded. Counting boluses on one
--   site and not the other would inflate eICU's n_agents and read as case
--   mix.
-- #####################################################################

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_agents_mimiciv` AS

WITH cohort AS (
  SELECT stay_id, intime
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

hours AS (
  SELECT
    c.stay_id,
    h AS hour_bin,
    TIMESTAMP_ADD(c.intime, INTERVAL h     HOUR) AS hr_start,
    TIMESTAMP_ADD(c.intime, INTERVAL h + 1 HOUR) AS hr_end
  FROM cohort c, UNNEST(GENERATE_ARRAY(0, 23)) AS h
),

items AS (
  SELECT itemid, label, LOWER(label) AS l
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
)

SELECT
  h.stay_id,
  h.hour_bin,
  CASE
    WHEN i.l LIKE '%dobutamine%' OR i.l LIKE '%milrinone%' THEN 'inotrope'
    ELSE 'vasopressor'
  END                                                  AS intervention,
  CASE
    -- ORDER MATTERS. norepinephrine must be tested before epinephrine.
    WHEN i.l LIKE '%norepinephrine%' OR i.l LIKE '%levophed%'
      THEN 'norepinephrine'
    WHEN i.l LIKE '%epinephrine%'    OR i.l LIKE '%adrenaline%'
      THEN 'epinephrine'
    WHEN i.l LIKE '%phenylephrine%'  OR i.l LIKE '%neosynephrine%'
      THEN 'phenylephrine'
    WHEN i.l LIKE '%vasopressin%'
      THEN 'vasopressin'
    WHEN i.l LIKE '%dopamine%'
      THEN 'dopamine'
    WHEN i.l LIKE '%dobutamine%'
      THEN 'dobutamine'
    WHEN i.l LIKE '%milrinone%'
      THEN 'milrinone'
  END                                                  AS agent
FROM hours h
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie
  ON h.stay_id   = ie.stay_id
 AND ie.starttime < h.hr_end
 AND ie.endtime   > h.hr_start
JOIN items i
  ON ie.itemid = i.itemid
WHERE ie.rate IS NOT NULL           -- infusions only, never boluses
  AND ie.rate > 0
  AND (i.l LIKE '%norepinephrine%' OR i.l LIKE '%levophed%'
    OR i.l LIKE '%epinephrine%'    OR i.l LIKE '%adrenaline%'
    OR i.l LIKE '%phenylephrine%'  OR i.l LIKE '%neosynephrine%'
    OR i.l LIKE '%vasopressin%'
    OR i.l LIKE '%dopamine%'
    OR i.l LIKE '%dobutamine%'
    OR i.l LIKE '%milrinone%')
GROUP BY 1, 2, 3, 4;

-- ---------------------------------------------------------------------
-- CHECKS FOR PART A — run before trusting anything downstream.
-- ---------------------------------------------------------------------
-- A1. Agent vocabulary actually matched. Any NULL agent means a label
--     slipped past the CASE while passing the WHERE — impossible unless
--     the two lists drift apart, so this must return zero.
-- SELECT COUNTIF(agent IS NULL) AS n_unclassified FROM `...v2_agents_mimiciv`;
--
-- A2. Norepinephrine must be the most common agent by a wide margin. If
--     epinephrine outnumbers it, the CASE branches got reordered.
-- SELECT agent, COUNT(DISTINCT stay_id) AS n_stays
-- FROM `...v2_agents_mimiciv` GROUP BY 1 ORDER BY n_stays DESC;
--
-- A3. Label audit — confirm nothing was silently dropped:
-- SELECT DISTINCT label FROM `physionet-data.mimiciv_3_1_icu.d_items`
-- WHERE LOWER(label) LIKE '%pressin%' OR LOWER(label) LIKE '%amine%'
--    OR LOWER(label) LIKE '%phrine%';
--
-- A4. Coverage vs the NEE route. Stays with any agent row should closely
--     track stays with any vasopressor row in v2_interventions_mimiciv.
--     A large gap means the derived NEE table includes something the
--     label match misses, or vice versa.




