-- =====================================================================
-- v2 STEP 6: MIMIC-IV severity-score comparators (APACHE and SOFA)
--
-- ONE TABLE, TWO SCORES, because they answer two different halves of one
-- question and because splitting them would duplicate the seven columns
-- they share. Read section A for APACHE and section B for SOFA. The file
-- was called `v2_06_apache_*` until SOFA was added on 2026-08-30; the
-- rename is not cosmetic, since a table named for one score that carries
-- two is precisely the kind of quiet misnaming this project refuses
-- elsewhere (`mbp`, not `map`).
--
-- =====================================================================
-- SECTION B: WHY SOFA IS HERE AND WHY IT IS NOT A THIRD BASELINE
-- =====================================================================
-- APACHE II is a pure-physiology score, so the honest comparison is
-- against `L_meas`. SOFA folds vasopressor dose into its cardiovascular
-- component, so the honest comparison is against `L_full`, which carries
-- `L_intv`. The two scores therefore externalise the project's own
-- meas/full split rather than duplicating each other, which is the
-- specific reason `docs/v2_analysis_tiering.md` item 5 ("one severity
-- score, not both") does not settle this. That rule is about
-- DISCRIMINATION baselines and it is correct about those. It does not
-- cover the two structural reasons SOFA is here:
--
--   1. SOFA IS THE THING THE DESIGN GENERALISES.
--      `docs/v2_analytical_design_plan.md` line 12 names it as lineage:
--      "SOFA already bundles measurement with intervention (vasopressor
--      dose in the cardiovascular component), ad hoc, for one organ
--      system". The representational claim of the paper is that this is
--      generalised from one organ to twelve signal-intervention pairs.
--      Not extracting SOFA means the one prior construction the design
--      explicitly extends never appears as a number anywhere.
--      `L_mbp^full` against SOFA's cardiovascular component is the
--      sharpest single comparison available in the project: same signal,
--      same intervention, both bundled, one hand-designed and one fitted.
--
--   2. IT IS THE ONLY PUBLISHED SCORE THAT DECOMPOSES ONTO
--      `config/domains.csv`. That file already carries a `sofa_organ`
--      column mapping six of the eleven domains onto SOFA organs, frozen
--      before any result existed. Layer 2 emits one number per domain,
--      and SOFA gives each of those six an externally defined referent on
--      organ definitions the project already committed to. APACHE II
--      cannot do this: its twelve variables correspond to nothing in
--      layer 2.
--
-- SOFA COSTS ALMOST NO NEW EXTRACTION. Every input is already in the
-- pipeline. Respiration needs PaO2/FiO2, which section A already pulls.
-- Cardiovascular needs MAP and vasopressor dose, and `dx_nee_peak` is
-- already harmonised across both sites. Coagulation is platelet, liver is
-- bilirubin, CNS is GCS, renal is creatinine and urine output. What is
-- added below is two more signals off the hourly lattice, three columns
-- off the intervention table, and one join to the derived SOFA concept.
--
-- THE ONE REAL APPROXIMATION, and it is in the component that matters
-- most. SOFA's cardiovascular bands are AGENT-SPECIFIC: dopamine at or
-- below 5 or dobutamine at any dose scores 2, dopamine above 5 or
-- adrenaline or noradrenaline at or below 0.1 scores 3, and the same
-- agents above 0.1 score 4. `dx_nee_peak` is a norepinephrine-EQUIVALENT
-- dose, so it collapses the agent distinction. The consequence is
-- specific and must be stated rather than absorbed: a patient on dopamine
-- at 5 mcg/kg/min converts to roughly 0.05 NEE and will score 3 under the
-- collapse where the original scores 2. `inotrope` is carried separately
-- below so that the dobutamine tier is recovered exactly. The collapse is
-- arguably an improvement on a 1996 agent list, but it is a deviation
-- from the published definition and the paper must say so.
--
-- =====================================================================
-- SECTION A: THE APACHE COMPARATOR
-- =====================================================================
-- WHY THIS TABLE EXISTS. Every comparison arm the project has so far
-- (`xgb_raw`, `xgb_feat`, `xgb_l`) changes the LEARNER or the
-- AGGREGATION while holding the measurements fixed. None of them changes
-- the FEATURE CONSTRUCTION for a construction someone else designed.
-- APACHE is that arm: it takes essentially our 19 signals, reduces each
-- to a worst-value-in-24h, maps that value through a published point
-- table, and adds the points up. Same measurements, same 24h window, a
-- completely different and much older answer to the question "how do I
-- turn a physiologic time series into evidence".
--
-- TWO SCORES ARE EXTRACTED, AND THEY DO DIFFERENT JOBS.
--
--   1. THE NATIVE SCORE. MIMIC ships `mimiciv_derived.apsiii`, the
--      APACHE III Acute Physiology Score, computed over intime -> +24h,
--      which is our landmark window exactly. eICU ships APACHE IVa's
--      `acutephysiologyscore` in `apachePatientResult`. Both are
--      physiology-only by construction: APACHE's age and chronic-health
--      points are added on top of the APS, not inside it. So each is
--      already the "physiology-only" baseline that
--      docs/v2_state_20260828.md section 5.3 asks for on the SAPS-II
--      side. They are the strongest WITHIN-SITE baseline available and
--      they cost one join. They are NOT comparable to each other: APS III
--      and APACHE IVa are different scores with different variable sets
--      and different point tables. Never put them on one axis.
--
--   2. THE RECOMPUTED APACHE II APS. This file also extracts the raw
--      worst-in-first-24h inputs for all twelve APACHE II physiologic
--      variables, and R/09c_apache.R turns them into points with ONE
--      scoring function applied at both sites. That is the score that
--      carries the cross-site claim, and it is why the extraction stops
--      at raw values rather than emitting a total: a point table
--      duplicated in two SQL files is a point table that will diverge,
--      and hard rule 5 says no code branches on site.
--
-- APACHE II IS THE RIGHT RECOMPUTED SCORE, not APACHE III or IV.
--   * Its variable set is almost exactly ours. Nine of its twelve
--     physiologic variables are among our 19 signals.
--   * docs/mimiciv_extraction.md section 4.6 records that this project's
--     own vital-sign reference ranges ARE the APACHE II zero-point bands.
--     So the contrast is unusually clean: both constructions start from
--     the same normal ranges and diverge only in what they do with a
--     deviation. APACHE counts a band; we shrink a proportion and fit a
--     smooth.
--   * It is fully published, so it is recomputable at both sites.
--     APACHE IV's coefficients are not, which is exactly why eICU's
--     native score cannot travel.
--   * docs/v2_analysis_tiering.md item 5 names APACHE-II specifically.
--
-- WHERE THE INPUTS COME FROM, and why it is not all one source.
--   * The nine variables we already extract (temperature, MAP, heart
--     rate, respiratory rate, sodium, creatinine, WBC, bicarbonate and
--     GCS) are read from OUR OWN hourly tables, `v2_hourly_mimiciv` and
--     `v2_hourly_gcs_mimiciv`. This is deliberate. Reading them from
--     source instead would let extraction differences -- unit
--     harmonisation, artifact rejection, hourly binning -- leak into a
--     comparison that is supposed to isolate feature CONSTRUCTION. Held
--     this way, the nine variables are literally the same numbers our
--     GAMs see, and the only thing that differs is what is computed from
--     them.
--   * The variables we do not carry (potassium, haematocrit, arterial pH
--     and the oxygenation triple PaO2 / PaCO2 / FiO2) have to come from
--     source. Handicapping the baseline by dropping them would be the
--     same error as denying `xgb_raw` its `n_obs` column.
--
-- OUTPUT CONTRACT. This table must emit EXACTLY the columns of
-- R/02_validate.R CONTRACT$severity (minus `site`, which the loader
-- stamps), in the same types, and the eICU port must emit the same set.
-- Check 9 compares types across sites.
--
-- HARD RULE 8. Nothing here is fitted. This is an extraction; the point
-- table lives in R and is frozen there, not re-derived per site.
-- =====================================================================

-- =====================================================================
-- VERIFY BEFORE RUNNING. None of this can be checked from outside
-- BigQuery, and every one of them silently degrades the comparison if it
-- is wrong. Run these four, read the answers, then run the file.
-- =====================================================================
-- V1. The derived APS III table exists and is keyed on stay_id.
--     SELECT column_name, data_type
--     FROM `physionet-data.mimiciv_3_1_derived.INFORMATION_SCHEMA.COLUMNS`
--     WHERE table_name = 'apsiii' ORDER BY ordinal_position;
--     Expect: stay_id, subject_id, hadm_id, apsiii, apsiii_prob, and one
--     _score column per component. If the table is absent, the
--     `aps_native` block below returns all NULL and ONLY the recomputed
--     APACHE II arm is available. That is survivable, because the
--     recomputed arm is the one carrying the cross-site claim. Say so
--     rather than substituting a different score.
--
-- V2. Its window is intime -> intime + 24h. mimic-code's apsiii.sql uses
--     the first 24 hours of the ICU stay, which matches `landmark_time`.
--     CONFIRM IT rather than assuming: a score computed over a different
--     window is a different score, and the paper must say which.
--
-- V3. The blood gas table exists and carries a specimen label.
--     SELECT column_name FROM
--     `physionet-data.mimiciv_3_1_derived.INFORMATION_SCHEMA.COLUMNS`
--     WHERE table_name = 'bg';
--     Expect at least: subject_id, hadm_id, charttime, specimen, ph, po2,
--     pco2, fio2. `bg` is keyed on subject_id / hadm_id and NOT on
--     stay_id, which is why the join below goes through the cohort's time
--     window rather than through stay_id.
--
-- V3b. IS `bg.fio2` A PERCENTAGE OR A FRACTION? Everything downstream
--     assumes 21-100 and divides by 100 once, in R. Check before trusting
--     any oxygenation number:
--     SELECT MIN(fio2), APPROX_QUANTILES(fio2, 4), MAX(fio2)
--     FROM `physionet-data.mimiciv_3_1_derived.bg` WHERE specimen = 'ART.'
--       AND fio2 IS NOT NULL;
--     Expect a minimum near 21 and a maximum near 100. A maximum near 1.0
--     means it is a FRACTION, and both the `bg_pf` CTE below and
--     `.ap2_aado2()` in R/09c must drop their /100.
--
-- V4. Join coverage, AFTER running. The APS III table is computed over
--     all icustays and our cohort is a filtered subset, so a LEFT JOIN
--     should be near-complete. But "near" needs a number, and that number
--     belongs in the paper.
--     SELECT COUNT(*) AS n,
--            COUNTIF(aps_native IS NULL)     AS n_no_native,
--            AVG(ap2_n_vars_present)         AS mean_vars_present,
--            COUNTIF(ap2_n_vars_present < 8) AS n_thin
--     FROM `...v2_severity_mimiciv`;
--
-- V5. THE DERIVED SOFA CONCEPT. Its exact column names are the largest
--     schema risk in this file, because mimic-code has emitted them under
--     more than one naming convention across releases.
--     SELECT column_name, data_type
--     FROM `physionet-data.mimiciv_3_1_derived.INFORMATION_SCHEMA.COLUMNS`
--     WHERE table_name = 'first_day_sofa' ORDER BY ordinal_position;
--     This file assumes `stay_id, sofa, respiration, coagulation, liver,
--     cardiovascular, cns, renal`. If the release suffixes them
--     (`respiration_24hours` and so on), EDIT THE `sofa_native` CTE
--     below; do not drop the subscores to make the query run. The
--     per-organ native values are the whole point of the domain-level
--     comparison, and a total-only fallback silently reduces the SOFA arm
--     to a second APACHE arm.
--     If `first_day_sofa` is absent entirely, the hourly `sofa` concept
--     aggregated with MAX over hr in [0, 23] is the equivalent; that is
--     what `first_day_sofa` is.
-- =====================================================================

CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_severity_mimiciv`
CLUSTER BY stay_id AS

WITH cohort AS (
  SELECT stay_id, subject_id, hadm_id, intime, landmark_time, weight_kg
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

-- ---------------------------------------------------------------------
-- BLOCK 1. The native score: APACHE III's Acute Physiology Score.
--
-- Taken whole from the derived table, never re-derived. mimic-code's
-- apsiii.sql is a published implementation and the entire point of a
-- baseline is that someone else defined it. `apsiii_prob` is the
-- concept's own logistic mapping and is carried so the arm has a
-- calibrated probability and not only a rank.
-- ---------------------------------------------------------------------
native AS (
  SELECT
    stay_id,
    CAST(apsiii AS INT64)        AS aps_native,
    CAST(apsiii_prob AS FLOAT64) AS aps_native_prob
  FROM `physionet-data.mimiciv_3_1_derived.apsiii`
),

-- ---------------------------------------------------------------------
-- BLOCK 2. The nine APACHE II variables we already extract, read back
-- out of our own hourly lattice.
--
-- MIN and MAX over the 24h window, never a pre-resolved "worst". Which
-- end is worst is a property of the POINT TABLE, and the point table
-- lives in R. Resolving it here would put half the scoring rule in SQL
-- and half in R, and the SQL half would have to be duplicated in the eICU
-- file, which is exactly the divergence hard rule 5 exists to prevent.
--
-- `hr_min` / `hr_max` are the within-hour extremes, so MIN(hr_min) and
-- MAX(hr_max) are the true 24h extremes of the binned series: the same
-- values `value_min` / `value_max` carry in the feature table.
-- ---------------------------------------------------------------------
ours AS (
  SELECT
    stay_id,
    MIN(IF(signal = 'temperature', hr_min, NULL)) AS temp_min,
    MAX(IF(signal = 'temperature', hr_max, NULL)) AS temp_max,
    MIN(IF(signal = 'mbp',         hr_min, NULL)) AS mbp_min,
    MAX(IF(signal = 'mbp',         hr_max, NULL)) AS mbp_max,
    MIN(IF(signal = 'heart_rate',  hr_min, NULL)) AS hr_min_v,
    MAX(IF(signal = 'heart_rate',  hr_max, NULL)) AS hr_max_v,
    MIN(IF(signal = 'resp_rate',   hr_min, NULL)) AS rr_min,
    MAX(IF(signal = 'resp_rate',   hr_max, NULL)) AS rr_max,
    MIN(IF(signal = 'sodium',      hr_min, NULL)) AS na_min,
    MAX(IF(signal = 'sodium',      hr_max, NULL)) AS na_max,
    MIN(IF(signal = 'bicarbonate', hr_min, NULL)) AS hco3_min,
    MAX(IF(signal = 'bicarbonate', hr_max, NULL)) AS hco3_max,
    MIN(IF(signal = 'creatinine',  hr_min, NULL)) AS creat_min,
    MAX(IF(signal = 'creatinine',  hr_max, NULL)) AS creat_max,
    MIN(IF(signal = 'wbc',         hr_min, NULL)) AS wbc_min,
    MAX(IF(signal = 'wbc',         hr_max, NULL)) AS wbc_max,
    -- SOFA-only. Platelet and bilirubin are two of our nineteen signals
    -- but score no APACHE II points, so they appear here and nowhere in
    -- the `ap2_` block. Same lattice, same window, same numbers the GAMs
    -- see — which is what makes SOFA's coagulation and liver components
    -- like-for-like with `L_platelet` and `L_bilirubin_total`.
    MIN(IF(signal = 'platelet',        hr_min, NULL)) AS plt_min,
    MAX(IF(signal = 'bilirubin_total', hr_max, NULL)) AS bili_max,
    -- Urine output is a RATE in our tables (mL/kg/hr, one value per
    -- hour). APACHE II's acute-renal-failure doubling rule is stated in
    -- mL/day, so reconstitute a volume: rate x weight, summed over the
    -- covered hours. Hours with no reading contribute nothing, so this
    -- UNDER-states the daily volume whenever coverage is partial, which
    -- biases the ARF flag towards firing. That is one of the two reasons
    -- the doubling is OFF by default in config (`apache.arf_doubling`).
    SUM(IF(signal = 'urine_output_rate', hr_med, NULL)) AS uo_rate_sum
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_mimiciv`
  WHERE hour_bin BETWEEN 0 AND 23
  GROUP BY stay_id
),

-- ---------------------------------------------------------------------
-- BLOCK 3. GCS, twice, because the two sources answer different
-- questions and the difference between them is not small.
--
--   ap2_gcs_min_native  the site's canonical total GCS, from
--                       `mimiciv_derived.gcs`, which already applies the
--                       standard intubated-patient convention. This is
--                       what an APACHE II implementation would use and it
--                       is the PRIMARY input.
--   ap2_gcs_min_ours    the total rebuilt from OUR hourly components,
--                       which drop verbal at every timepoint a clinician
--                       marked the patient unassessable. Verbal is
--                       imputed as 1 for an hour where it is absent.
--
-- Our masking rule is a design choice of THIS pipeline. Importing it into
-- the baseline would handicap the baseline on exactly the sickest
-- patients, the ones who are intubated, so the primary arm does not. The
-- second column exists so that choice is a measured sensitivity rather
-- than an assertion. eICU has no `gcs_unable` analogue at all
-- (sql/eicu/v2_04_hourly_gcs_eicu.sql header), so the two columns are
-- expected to coincide there and to differ here. That asymmetry is a
-- finding about the two datasets, not a bug.
-- ---------------------------------------------------------------------
gcs_native AS (
  SELECT
    c.stay_id,
    MIN(SAFE_CAST(g.gcs AS INT64)) AS gcs_min_native
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_derived.gcs` g
    ON c.stay_id = g.stay_id
   AND g.charttime >= c.intime
   AND g.charttime <  c.landmark_time
  WHERE SAFE_CAST(g.gcs AS INT64) BETWEEN 3 AND 15
  GROUP BY c.stay_id
),

-- MEASURED 2026-08-30, audit A4: the single-constant version of this column
-- was not measuring what it claimed. At MIMIC it differs from the native
-- total on 41.4% of stays with a mean gap of 3.80 points, and eICU -- which
-- masks nothing -- shows a gap of 0.12. So the gap is the VERBAL CONVENTION
-- and not hourly binning, and its magnitude is consistent with the derived
-- concept treating unassessable verbal as NORMAL where we treated it as
-- WORST. A 3.8-point shift is 3.8 APACHE II points and up to 2 SOFA CNS
-- points on two fifths of the cohort, which is not a minor implementation
-- detail.
--
-- TWO COLUMNS, NOT A TUNED CONSTANT. The obvious repair is to pick an
-- imputation value that closes the gap against the native score. That is
-- backwards: this column exists to show what OUR masking rule does
-- differently, so fitting its constant until it agrees with the thing it is
-- a sensitivity against would leave it measuring nothing. Instead both
-- defensible conventions are emitted and the paper reports the bracket:
--
--   gcs_min_ours        verbal imputed 1 -- the PESSIMISTIC bound. Asserts an
--                       unassessable patient is unresponsive.
--   gcs_min_ours_vnorm  verbal imputed 5 -- the OPTIMISTIC bound, and the
--                       standard prospective convention ("assume normal if
--                       intubated").
--
-- The truth for any given stay is inside that interval and the two bound the
-- baseline's strength. Dropping the hour entirely was considered and
-- rejected: it would leave deeply sedated intubated patients with no GCS at
-- all, and an unmeasured variable scores zero points, so the baseline would
-- collapse on exactly the sickest subgroup -- the same asymmetry in a new
-- place, and a worse one.
gcs_ours AS (
  SELECT
    stay_id,
    MIN(gcs_total_v1)    AS gcs_min_ours,
    MIN(gcs_total_vnorm) AS gcs_min_ours_vnorm
  FROM (
    SELECT
      stay_id,
      hour_bin,
      MIN(IF(signal = 'gcs_motor', hr_min, NULL)) +
      MIN(IF(signal = 'gcs_eyes',  hr_min, NULL)) +
      COALESCE(MIN(IF(signal = 'gcs_verbal', hr_min, NULL)), 1) AS gcs_total_v1,
      MIN(IF(signal = 'gcs_motor', hr_min, NULL)) +
      MIN(IF(signal = 'gcs_eyes',  hr_min, NULL)) +
      COALESCE(MIN(IF(signal = 'gcs_verbal', hr_min, NULL)), 5) AS gcs_total_vnorm
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_gcs_mimiciv`
    WHERE hour_bin BETWEEN 0 AND 23
    GROUP BY stay_id, hour_bin
  )
  -- Bound on the pessimistic total; motor+eyes lies in [2, 10], so v1 in
  -- [3, 11] and vnorm in [7, 15] follow and neither can leave the scale.
  WHERE gcs_total_v1 BETWEEN 3 AND 15
  GROUP BY stay_id
),

-- ---------------------------------------------------------------------
-- BLOCK 4. Potassium and haematocrit, which our design does not carry.
--
-- The itemid lists are short and explicit. Audit A2 at the foot of this
-- file prints what they matched and MUST be read before the table is
-- trusted: a missed variant itemid collapses coverage in a way that looks
-- like missing data rather than like a bad whitelist.
-- ---------------------------------------------------------------------
labs_extra AS (
  SELECT
    c.stay_id,
    MIN(IF(l.itemid IN (50971, 50822), l.valuenum, NULL)) AS k_min,
    MAX(IF(l.itemid IN (50971, 50822), l.valuenum, NULL)) AS k_max,
    MIN(IF(l.itemid IN (51221, 50810), l.valuenum, NULL)) AS hct_min,
    MAX(IF(l.itemid IN (51221, 50810), l.valuenum, NULL)) AS hct_max
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON c.hadm_id = l.hadm_id
   AND l.charttime >= c.intime
   AND l.charttime <  c.landmark_time
  WHERE l.itemid IN (50971, 50822, 51221, 50810)
    AND l.valuenum IS NOT NULL
    -- Plausibility bounds, in the same spirit as signal_spec. Serum K
    -- outside 1-10 mmol/L and haematocrit outside 5-80 percent are assay
    -- or transcription artifacts, and APACHE II's outer band is 4 points
    -- either way, so a single artifact would max out the variable.
    AND ((l.itemid IN (50971, 50822) AND l.valuenum BETWEEN 1.0 AND 10.0)
      OR (l.itemid IN (51221, 50810) AND l.valuenum BETWEEN 5.0 AND 80.0))
  GROUP BY c.stay_id
),

-- ---------------------------------------------------------------------
-- BLOCK 5. Arterial pH and the oxygenation triple.
--
-- APACHE II's fifth variable is a BRANCH, not a value: score A-aDO2 when
-- FiO2 >= 0.5 and PaO2 otherwise. A-aDO2 needs PaO2, PaCO2 and FiO2 from
-- THE SAME blood gas. Taking the worst of each independently would
-- assemble a gas that never existed and would systematically overstate
-- the gradient, so the single worst arterial gas is selected by lowest
-- PaO2 and its own PaCO2 and FiO2 travel with it.
--
-- pH is taken as a min/max over all arterial gases rather than from the
-- selected one, because pH has its own two-sided point table and the
-- worst pH need not sit on the worst-oxygenation gas.
--
-- `bg` is keyed on subject_id / hadm_id, so the window join is by time.
-- ---------------------------------------------------------------------
bg_arterial AS (
  SELECT
    c.stay_id,
    SAFE_CAST(b.ph   AS FLOAT64) AS ph,
    SAFE_CAST(b.po2  AS FLOAT64) AS po2,
    SAFE_CAST(b.pco2 AS FLOAT64) AS pco2,
    SAFE_CAST(b.fio2 AS FLOAT64) AS fio2
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_derived.bg` b
    ON c.hadm_id = b.hadm_id
   AND b.charttime >= c.intime
   AND b.charttime <  c.landmark_time
  WHERE b.specimen = 'ART.'
),

bg_ph AS (
  SELECT stay_id, MIN(ph) AS ph_min, MAX(ph) AS ph_max
  FROM bg_arterial
  WHERE ph BETWEEN 6.5 AND 8.0
  GROUP BY stay_id
),

bg_oxy AS (
  SELECT
    stay_id,
    worst.po2  AS pao2,
    worst.pco2 AS paco2,
    worst.fio2 AS fio2
  FROM (
    SELECT
      stay_id,
      ARRAY_AGG(STRUCT(po2, pco2, fio2) ORDER BY po2 ASC LIMIT 1)[OFFSET(0)] AS worst
    FROM bg_arterial
    WHERE po2 BETWEEN 20 AND 700
    GROUP BY stay_id
  )
),

-- ---------------------------------------------------------------------
-- BLOCK 5b. THE WORST P/F RATIO, which is a DIFFERENT GAS from the worst
-- PaO2 and needs its own selection.
--
-- ADDED 2026-08-31 after the recomputed SOFA respiration component failed
-- its validation gate: 33.8% exact agreement against the derived concept
-- and a Spearman of 0.086, meaning we were not even ORDERING patients the
-- same way. Our mean respiration score was 0.285 against a native 1.916.
--
-- THE DIAGNOSIS, and the reason it is worth writing down. A units error
-- would have PRESERVED RANK: if `bg.fio2` were a fraction rather than a
-- percentage, every ratio would be inflated by exactly 100x, every stay
-- would collapse to tier 0, and the rank correlation would have stayed
-- high. Observing near-zero rank agreement proves that DIFFERENT GASES
-- were being selected for different patients.
--
-- And they were. `bg_oxy` above selects by lowest PaO2, which is correct
-- for APACHE II -- its oxygenation variable is PaO2-based when FiO2 < 0.5
-- -- and wrong for SOFA, whose respiration component is a RATIO. The
-- lowest-PaO2 gas is frequently a room-air sample: PaO2 95 on FiO2 21%
-- gives P/F 452 and scores 0, while a different gas at PaO2 120 on FiO2
-- 100% gives P/F 120 and scores 3 or 4. The old query systematically
-- picked the gas that scores BEST rather than worst.
--
-- Two corroborations from the same run: our respiration was non-zero on
-- only 6.0% of stays, and APACHE II's own oxygenation term scored
-- non-zero on 9.8% despite 48.2% coverage -- about one stay in five WITH
-- a gas earning any points, which is not credible for an ICU population.
--
-- So the ratio gets its own minimum, over gases where BOTH values are
-- present and plausible. `bg_oxy` is untouched and still feeds APACHE II.
-- Both columns are emitted; neither serves both scores.
bg_pf AS (
  SELECT
    stay_id,
    MIN(po2 / (fio2 / 100)) AS pf_min
  FROM bg_arterial
  WHERE po2  BETWEEN 20 AND 700
    AND fio2 BETWEEN 21 AND 100
  GROUP BY stay_id
),

-- ---------------------------------------------------------------------
-- BLOCK 6. Chronic health and admission type.
--
-- SECONDARY, AND SITE-APPROXIMATE. These are the two APACHE II terms that
-- do NOT transport cleanly. MIMIC has ICD-derived comorbidity; eICU has
-- APACHE's own prospectively collected chronic-health flags. They are
-- different constructs with different ascertainment, and no amount of
-- careful mapping makes them the same variable. The physiology comparison
-- is the headline precisely because it does not depend on them.
--
-- These columns exist for the third comparison of
-- docs/v2_state_20260828.md section 5.3, "what do our excluded covariates
-- buy", and every number computed from them must carry the caveat.
--
-- APACHE II's chronic-health criterion is biopsy-proven cirrhosis, portal
-- hypertension, NYHA class IV heart failure, severe chronic respiratory
-- disease, chronic dialysis or immunocompromise. Charlson's components
-- are the closest available and are not the same definition: Charlson
-- ascertains "has the diagnosis", APACHE ascertains "has it severely,
-- before this admission".
--
-- ADMISSION CLASS is the weaker of the two. MIMIC has no postoperative
-- flag, so surgical status is inferred from the first ICU care unit,
-- which is a unit-assignment proxy and will misclassify surgical patients
-- admitted to a medical unit and the reverse. Stated, not fixed.
-- ---------------------------------------------------------------------
chronic AS (
  SELECT
    c.stay_id,
    CAST(GREATEST(
      COALESCE(ch.aids, 0),
      COALESCE(ch.malignant_cancer, 0),
      COALESCE(ch.metastatic_solid_tumor, 0)
    ) AS INT64) AS chronic_immunocompromised,
    CAST(GREATEST(
      COALESCE(ch.severe_liver_disease, 0),
      COALESCE(ch.mild_liver_disease, 0),
      COALESCE(ch.congestive_heart_failure, 0),
      COALESCE(ch.chronic_pulmonary_disease, 0),
      COALESCE(ch.renal_disease, 0)
    ) AS INT64) AS chronic_severe_organ
  FROM cohort c
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.charlson` ch
    ON c.hadm_id = ch.hadm_id
),

admission AS (
  SELECT
    c.stay_id,
    CASE
      WHEN icu.first_careunit LIKE '%Surgical%'
        OR icu.first_careunit LIKE '%SICU%'
        OR icu.first_careunit LIKE '%CVICU%'
        OR icu.first_careunit LIKE '%Trauma%'
      THEN CASE WHEN adm.admission_type IN ('ELECTIVE', 'SURGICAL SAME DAY ADMISSION')
                THEN 'elective_postop' ELSE 'emergency_postop' END
      ELSE 'nonoperative'
    END AS admission_class
  FROM cohort c
  JOIN `physionet-data.mimiciv_3_1_icu.icustays` icu   ON c.stay_id = icu.stay_id
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm ON c.hadm_id = adm.hadm_id
),

-- Chronic dialysis is the other half of the ARF doubling rule: the
-- doubling applies to ACUTE failure only, so a chronic dialysis patient
-- must be excluded from it. Read from our own intervention table, which
-- is where `rrt` already lives, so that no second definition of renal
-- replacement enters the project.
dialysis AS (
  SELECT stay_id, 1 AS dialysis
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_intervention_features_mimiciv`
  WHERE intervention = 'rrt' AND ever_active = 1
),

-- ---------------------------------------------------------------------
-- BLOCK 7 (SOFA). The native score, taken whole from the derived concept.
--
-- Same role `apsiii` plays for APACHE: a published implementation to
-- validate our recomputation against at MIMIC, so that the recomputation
-- can then be transported to eICU, which ships no SOFA at all. The six
-- subscores matter more here than the total, because they are what the
-- domain-level comparison against `config/domains.csv` uses.
--
-- SEE VERIFY V5. These column names are the largest schema risk in the
-- file.
-- ---------------------------------------------------------------------
sofa_native AS (
  SELECT
    stay_id,
    CAST(sofa           AS INT64) AS sofa_total,
    CAST(respiration    AS INT64) AS sofa_respiration,
    CAST(coagulation    AS INT64) AS sofa_coagulation,
    CAST(liver          AS INT64) AS sofa_liver,
    CAST(cardiovascular AS INT64) AS sofa_cardiovascular,
    CAST(cns            AS INT64) AS sofa_cns,
    CAST(renal          AS INT64) AS sofa_renal
  FROM `physionet-data.mimiciv_3_1_derived.first_day_sofa`
),

-- ---------------------------------------------------------------------
-- BLOCK 8 (SOFA). The three intervention inputs SOFA needs and APACHE
-- does not.
--
-- All three come from our own intervention feature table rather than
-- from source, for the same reason the nine shared physiologic variables
-- do: the comparison is supposed to isolate construction, so the
-- vasopressor dose SOFA sees must be the same number `L_mbp^full` sees.
--
-- `dx_nee_peak` carries the `dx_` prefix and is therefore REFUSED by
-- R/05_formula.R's guard, correctly — it is a diagnostic column and may
-- never enter a model frame. That guard is about formulas. This table
-- never reaches a model frame at all, so reading it here is not a
-- workaround; it is the column being used for the one thing it was
-- extracted for.
--
-- Invasive and non-invasive ventilation are kept as SEPARATE flags. SOFA
-- scores of 3 and 4 require respiratory support, and the original 1996
-- definition says mechanical ventilation while mimic-code's reference
-- implementation counts non-invasive support too. Emitting both lets
-- R/09d_sofa.R match whichever definition the validation run shows the
-- derived concept used, and lets the other be reported as a sensitivity,
-- instead of freezing the choice here where it could not be revisited.
-- ---------------------------------------------------------------------
sofa_intv AS (
  SELECT
    stay_id,
    CAST(MAX(IF(intervention = 'invasive_vent',    ever_active, 0)) AS INT64) AS vent_inv,
    CAST(MAX(IF(intervention = 'noninvasive_vent', ever_active, 0)) AS INT64) AS vent_niv,
    CAST(MAX(IF(intervention = 'inotrope',         ever_active, 0)) AS INT64) AS inotrope,
    -- MEASURED 2026-08-30: 5,218 of eICU's 13,430 vasopressor stays (38.9%)
    -- have `ever_active = 1` and a NULL `dx_nee_peak`, because eICU's dose
    -- lives in free-text `infusiondrug.drugrate` and does not always parse.
    -- Without this flag those stays are indistinguishable from stays that
    -- never received a pressor, and R/09d would score them on the MAP tier
    -- alone -- silently weakening eICU's SOFA cardiovascular component, which
    -- is the exact component the SOFA arm exists to compare. The flag lets
    -- .sofa_cardio() floor them at tier 3 instead.
    CAST(MAX(IF(intervention = 'vasopressor',      ever_active, 0)) AS INT64) AS vaso_active,
    -- PLAUSIBILITY GUARD, in the same spirit as signal_spec and added for the
    -- same reason. MEASURED 2026-08-30: the maximum norepinephrine equivalent
    -- is 94.7 at MIMIC and 123.1 at eICU, where anything above roughly 2
    -- mcg/kg/min is already extreme. The medians are sensible (0.285 in
    -- SOFA's top native tier), so this is an outlier tail rather than a unit
    -- error -- but SOFA's cardiovascular tiers are dose THRESHOLDS, so a
    -- single artifact maxes the component out. 5.0 is generous: it is well
    -- above any defensible clinical ceiling and still removes the impossible
    -- values. Out-of-range doses become NULL, and `vaso_active` then floors
    -- the stay at tier 3 rather than dropping it to the MAP tier.
    MAX(IF(intervention = 'vasopressor'
           AND dx_nee_peak > 0 AND dx_nee_peak <= 5.0, dx_nee_peak, NULL))    AS nee_peak
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_intervention_features_mimiciv`
  GROUP BY stay_id
)

SELECT
  c.stay_id,

  -- --- native score --------------------------------------------------
  n.aps_native                                 AS aps_native,
  'apsiii'                                     AS aps_native_version,
  n.aps_native_prob                            AS aps_native_prob,
  -- APACHE III's TOTAL (APS + age + chronic health) is not shipped by the
  -- derived concept. Typed NULL so the column set matches eICU, where
  -- `apachescore` is the total and is present. Never modelled; carried so
  -- the two sites present one schema.
  CAST(NULL AS INT64)                          AS severity_total_native,

  -- --- APACHE II physiologic inputs, worst-in-first-24h -------------
  o.temp_min                                   AS ap2_temp_min,
  o.temp_max                                   AS ap2_temp_max,
  o.mbp_min                                    AS ap2_mbp_min,
  o.mbp_max                                    AS ap2_mbp_max,
  o.hr_min_v                                   AS ap2_hr_min,
  o.hr_max_v                                   AS ap2_hr_max,
  o.rr_min                                     AS ap2_rr_min,
  o.rr_max                                     AS ap2_rr_max,
  x.pao2                                       AS ap2_pao2,
  x.paco2                                      AS ap2_paco2,
  x.fio2                                       AS ap2_fio2,
  -- SOFA's respiration input. NOT derivable from the three columns above:
  -- they come from the worst-PaO2 gas and this is the worst RATIO, which
  -- is a different specimen. See block 5b.
  pf.pf_min                                    AS ap2_pf_min,
  p.ph_min                                     AS ap2_ph_min,
  p.ph_max                                     AS ap2_ph_max,
  o.hco3_min                                   AS ap2_hco3_min,
  o.hco3_max                                   AS ap2_hco3_max,
  o.na_min                                     AS ap2_sodium_min,
  o.na_max                                     AS ap2_sodium_max,
  e.k_min                                      AS ap2_potassium_min,
  e.k_max                                      AS ap2_potassium_max,
  o.creat_min                                  AS ap2_creatinine_min,
  o.creat_max                                  AS ap2_creatinine_max,
  e.hct_min                                    AS ap2_hematocrit_min,
  e.hct_max                                    AS ap2_hematocrit_max,
  o.wbc_min                                    AS ap2_wbc_min,
  o.wbc_max                                    AS ap2_wbc_max,
  gn.gcs_min_native                            AS ap2_gcs_min_native,
  go.gcs_min_ours                              AS ap2_gcs_min_ours,
  go.gcs_min_ours_vnorm                        AS ap2_gcs_min_ours_vnorm,

  -- --- acute renal failure inputs ------------------------------------
  -- mL over the covered hours, reconstituted from the rate. NULL when the
  -- stay has no usable weight, which is correct: unavailable, not zero.
  CAST(o.uo_rate_sum * c.weight_kg AS FLOAT64) AS ap2_urine_ml_24h,
  CAST(COALESCE(d.dialysis, 0) AS INT64)       AS ap2_dialysis,

  -- --- chronic health and admission type (SECONDARY) -----------------
  ch.chronic_immunocompromised                 AS chronic_immunocompromised,
  ch.chronic_severe_organ                      AS chronic_severe_organ,
  a.admission_class                            AS admission_class,

  -- --- provenance and coverage ---------------------------------------
  -- A VALUE, not a branch. Records which source supplied the variables
  -- our design does not carry, so the eICU approximation (a single worst
  -- value from apacheapsvar rather than a min/max pair) is visible in the
  -- data and not only in a document.
  'labevents+bg'                               AS ap2_gap_source,

  -- How many of the twelve APACHE II physiologic variables have at least
  -- one usable value. The validator holds this to a floor: a stay scored
  -- on four variables is not an APACHE II score, and averaging it in
  -- silently would flatter our arm by degrading the baseline.
  CAST(
    (CASE WHEN o.temp_min  IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.mbp_min   IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.hr_min_v  IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.rr_min    IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN x.pao2      IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN p.ph_min IS NOT NULL OR o.hco3_min IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.na_min    IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN e.k_min     IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.creat_min IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN e.hct_min   IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.wbc_min   IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN gn.gcs_min_native IS NOT NULL THEN 1 ELSE 0 END)
  AS INT64)                                    AS ap2_n_vars_present,

  -- --- SOFA: the native score, per organ ----------------------------
  -- The subscores are the point. `sofa_native` alone would make this a
  -- second APACHE arm; the six organs are what compare against the six
  -- domains of config/domains.csv that carry a `sofa_organ`.
  sn.sofa_total                                AS sofa_native,
  sn.sofa_respiration                          AS sofa_native_respiration,
  sn.sofa_coagulation                          AS sofa_native_coagulation,
  sn.sofa_liver                                AS sofa_native_liver,
  sn.sofa_cardiovascular                       AS sofa_native_cardiovascular,
  sn.sofa_cns                                  AS sofa_native_cns,
  sn.sofa_renal                                AS sofa_native_renal,

  -- --- SOFA: the recomputation's inputs -----------------------------
  -- Only the five SOFA needs that the APACHE block does not already
  -- carry. Respiration reuses `ap2_pao2` / `ap2_fio2`, cardiovascular
  -- reuses `ap2_mbp_min`, CNS reuses `ap2_gcs_min_native`, and renal
  -- reuses `ap2_creatinine_max` and `ap2_urine_ml_24h`. Sharing them is
  -- not a shortcut: it is what guarantees the two scores are computed
  -- from identical numbers, so any difference between the arms is the
  -- construction and never the input.
  o.plt_min                                    AS sofa_platelet_min,
  o.bili_max                                   AS sofa_bilirubin_max,
  COALESCE(si.vent_inv, 0)                     AS sofa_vent_invasive,
  COALESCE(si.vent_niv, 0)                     AS sofa_vent_noninvasive,
  COALESCE(si.vaso_active, 0)                  AS sofa_vasopressor,
  si.nee_peak                                  AS sofa_nee_peak,
  COALESCE(si.inotrope, 0)                     AS sofa_inotrope,

  -- How many of SOFA's six organs have a usable input. Same role as
  -- `ap2_n_vars_present`, and needed for the same reason: a missing SOFA
  -- component scores zero, which reads as a healthy organ rather than an
  -- unmeasured one. Cardiovascular is not counted as absent when the
  -- vasopressor columns are null, because "no vasopressor" is a real
  -- observation rather than a missing one; it needs MAP only.
  CAST(
    (CASE WHEN pf.pf_min IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.plt_min         IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.bili_max        IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.mbp_min IS NOT NULL OR si.vaso_active = 1
               OR si.inotrope = 1 THEN 1 ELSE 0 END) +
    (CASE WHEN gn.gcs_min_native IS NOT NULL THEN 1 ELSE 0 END) +
    (CASE WHEN o.creat_max IS NOT NULL OR o.uo_rate_sum IS NOT NULL THEN 1 ELSE 0 END)
  AS INT64)                                    AS sofa_n_organs_present

FROM cohort c
LEFT JOIN native      n  ON c.stay_id = n.stay_id
LEFT JOIN ours        o  ON c.stay_id = o.stay_id
LEFT JOIN gcs_native  gn ON c.stay_id = gn.stay_id
LEFT JOIN gcs_ours    go ON c.stay_id = go.stay_id
LEFT JOIN labs_extra  e  ON c.stay_id = e.stay_id
LEFT JOIN bg_ph       p  ON c.stay_id = p.stay_id
LEFT JOIN bg_oxy      x  ON c.stay_id = x.stay_id
LEFT JOIN bg_pf       pf ON c.stay_id = pf.stay_id
LEFT JOIN chronic     ch ON c.stay_id = ch.stay_id
LEFT JOIN admission   a  ON c.stay_id = a.stay_id
LEFT JOIN dialysis    d  ON c.stay_id = d.stay_id
LEFT JOIN sofa_native sn ON c.stay_id = sn.stay_id
LEFT JOIN sofa_intv   si ON c.stay_id = si.stay_id;


-- =====================================================================
-- AUDITS - run these, do not skip. Aggregates only.
-- =====================================================================
-- A1. Headline coverage. One row per cohort stay, no more and no fewer.
-- SELECT COUNT(*) AS n, COUNT(DISTINCT stay_id) AS n_distinct,
--        COUNTIF(aps_native IS NULL)       AS n_no_apsiii,
--        AVG(ap2_n_vars_present)           AS mean_vars,
--        COUNTIF(ap2_n_vars_present >= 10) AS n_ge_10,
--        COUNTIF(ap2_n_vars_present < 8)   AS n_thin
-- FROM `...v2_severity_mimiciv`;
--
-- A2. THE ITEMID AUDIT. Read this before trusting potassium or
--     haematocrit.
-- SELECT COUNTIF(ap2_potassium_min IS NOT NULL)  AS n_k,
--        COUNTIF(ap2_hematocrit_min IS NOT NULL) AS n_hct,
--        COUNT(*) AS n
-- FROM `...v2_severity_mimiciv`;
--     Expect both near-complete: potassium and a full blood count are
--     ordered on essentially every ICU admission. Anything below roughly
--     90 percent means the itemid list is wrong, NOT that the labs were
--     not drawn.
--
-- A3. Arterial gas coverage. This one is genuinely incomplete and the
--     fraction belongs in the paper: a patient with no ABG scores zero on
--     APACHE II's oxygenation variable by default, which UNDER-states the
--     baseline for exactly the least sick patients.
-- SELECT COUNTIF(ap2_pao2 IS NOT NULL) / COUNT(*)   AS frac_with_abg,
--        COUNTIF(ap2_ph_min IS NOT NULL) / COUNT(*) AS frac_with_ph,
--        COUNTIF(ap2_ph_min IS NULL AND ap2_hco3_min IS NOT NULL) / COUNT(*)
--          AS frac_hco3_substitution
-- FROM `...v2_severity_mimiciv`;
--
-- A4. The GCS disagreement. This is a RESULT, not a check: it measures
--     what our verbal masking does to a total GCS, on the site that has
--     the flag. eICU cannot produce this number.
-- SELECT COUNTIF(ap2_gcs_min_native != ap2_gcs_min_ours) AS n_differ,
--        AVG(ap2_gcs_min_native - ap2_gcs_min_ours)      AS mean_gap,
--        COUNT(*) AS n
-- FROM `...v2_severity_mimiciv`
-- WHERE ap2_gcs_min_native IS NOT NULL AND ap2_gcs_min_ours IS NOT NULL;
--
-- A3b. THE FIX'S OWN CHECK. Run this before re-running the R arm. The two
--      oxygenation selections must DISAGREE on a substantial minority of
--      stays -- if `ap2_pf_min` equals `ap2_pao2 / (ap2_fio2/100)` on
--      nearly every stay, the two CTEs are picking the same gas and the
--      fix has not done anything.
-- SELECT COUNT(*) AS n,
--        COUNTIF(ap2_pf_min IS NOT NULL) AS n_pf,
--        COUNTIF(ABS(ap2_pf_min - ap2_pao2 / (ap2_fio2/100)) < 1) AS n_same_gas,
--        APPROX_QUANTILES(ap2_pf_min, 4) AS pf_quartiles
-- FROM `...v2_severity_mimiciv` WHERE ap2_pf_min IS NOT NULL;
--      Expect the P/F quartiles to straddle 200-400. A median above 400
--      means most stays still score 0 and something else is wrong.
--
-- A5. Sanity against the native score. The recomputed APACHE II APS and
--     the derived APS III are different scores, so they must correlate
--     strongly WITHOUT agreeing. Computed in R by tests/metrics_severity.R,
--     because the recomputation lives there; noted here so the check is
--     not forgotten. A Spearman correlation below roughly 0.7 means the R
--     point table or one of the unit assumptions is wrong.
--
-- A6. Admission class distribution. The care-unit proxy is the weakest
--     construct in this file; look at the split before using it.
-- SELECT admission_class, COUNT(*), AVG(chronic_severe_organ)
-- FROM `...v2_severity_mimiciv` GROUP BY 1;
--
-- --- SOFA audits -----------------------------------------------------
--
-- A7. Native SOFA coverage and the subscore ranges. Every organ is
--     bounded 0-4 and the total 0-24; anything outside means the column
--     mapping in the `sofa_native` CTE is wrong (VERIFY V5), most likely
--     because a release renamed the subscores and the query silently
--     picked up a different column.
-- SELECT COUNTIF(sofa_native IS NULL) AS n_no_sofa,
--        MIN(sofa_native), MAX(sofa_native),
--        MAX(sofa_native_respiration), MAX(sofa_native_coagulation),
--        MAX(sofa_native_liver), MAX(sofa_native_cardiovascular),
--        MAX(sofa_native_cns), MAX(sofa_native_renal)
-- FROM `...v2_severity_mimiciv`;
--
-- A8. THE SUBSCORE SUM CHECK. The six organs must add to the total on
--     every row. If they do not, the derived concept's total is computed
--     over a different window from its components and the domain-level
--     comparison is built on sand.
-- SELECT COUNTIF(sofa_native_respiration + sofa_native_coagulation +
--                sofa_native_liver + sofa_native_cardiovascular +
--                sofa_native_cns + sofa_native_renal != sofa_native) AS n_mismatch
-- FROM `...v2_severity_mimiciv` WHERE sofa_native IS NOT NULL;
--
-- A9. SOFA input coverage, against A1's APACHE figures. Respiration is
--     the thin one at both sites because it needs an arterial gas.
-- SELECT AVG(sofa_n_organs_present) AS mean_organs,
--        COUNTIF(sofa_n_organs_present = 6) AS n_complete,
--        COUNTIF(sofa_platelet_min IS NOT NULL) / COUNT(*)  AS frac_plt,
--        COUNTIF(sofa_bilirubin_max IS NOT NULL) / COUNT(*) AS frac_bili,
--        AVG(sofa_vent_invasive) AS frac_vent_inv,
--        AVG(sofa_vent_noninvasive) AS frac_vent_niv,
--        COUNTIF(sofa_nee_peak IS NOT NULL) / COUNT(*) AS frac_on_pressor,
--        AVG(sofa_inotrope) AS frac_inotrope
-- FROM `...v2_severity_mimiciv`;
--
-- A10. WHICH VENTILATION DEFINITION THE DERIVED CONCEPT USED. This is
--      not a check, it is the measurement that decides
--      `config/sofa.resp_support` before the arm is run. Stays scoring 3
--      or 4 on native respiration must be ventilated under whichever
--      definition the concept applied; compare against each flag.
-- SELECT sofa_native_respiration,
--        COUNT(*) AS n,
--        AVG(sofa_vent_invasive) AS frac_inv,
--        AVG(GREATEST(sofa_vent_invasive, sofa_vent_noninvasive)) AS frac_any
-- FROM `...v2_severity_mimiciv`
-- WHERE sofa_native_respiration IS NOT NULL GROUP BY 1 ORDER BY 1;
--
-- A11. The cardiovascular approximation, sized. Our NEE collapse cannot
--      reproduce SOFA's dopamine-at-or-below-5 tier, so the recomputed
--      cardiovascular score will disagree with the native one on some
--      stays. Cross-tabulate before reporting the arm, and put the
--      disagreement rate in the paper rather than the point estimate
--      alone. R/09d_sofa.R writes this table too; the SQL version is
--      here so it can be checked before any parquet is exported.
-- SELECT sofa_native_cardiovascular, COUNT(*) AS n,
--        COUNTIF(sofa_nee_peak IS NULL) AS n_no_pressor,
--        AVG(sofa_inotrope) AS frac_inotrope,
--        APPROX_QUANTILES(sofa_nee_peak, 4) AS nee_quartiles
-- FROM `...v2_severity_mimiciv`
-- WHERE sofa_native_cardiovascular IS NOT NULL GROUP BY 1 ORDER BY 1;
