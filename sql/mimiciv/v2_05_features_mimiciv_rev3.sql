-- =====================================================================
-- v2 STEP 5 (rev 3): aggregation
--
-- Reads the materialised hourly lattice, never physionet-data. Re-running
-- this file costs a few hundred MB, not the ~12 GB chartevents scan in 3a.
-- Nothing below requires 3a, 3b or 04 to be re-run.
--
-- CHANGES FROM REV 2
--
--   1. SECTION 5c DELETED. The ordering flag (`o_flag`) is gone from the
--      design (2026-08-25). It stated when an intervention began relative
--      to a measurement excursion, which is a claim about treatment
--      sequencing in a framework that makes no causal claim -- a reviewer
--      reads `started_pre_excursion` as a causal assertion because that is
--      what the name says. It also belonged to neither the measurement nor
--      the intervention term set, so it was the one term preventing those
--      two sets from partitioning the joint model exactly. Removing it is
--      what makes L_full - L_intv the conditional LLR by arithmetic rather
--      than by argument.
--
--      `v2_features_ordering_mimiciv` is therefore NO LONGER PRODUCED. The
--      previously extracted copy is audit-only: it still encodes the
--      signal x intervention grid independently of pairing.csv, which the
--      R validator cross-checks when present.
--
--   2. TREND REDEFINED as an ordinary least-squares slope, scaled to
--      change-per-24h, replacing `median(late half) - median(early half)`.
--
--        slope = COVAR_POP(hour_bin, hr_med) / VAR_POP(hour_bin) * 24
--
--      Three reasons, and the third is the one that matters:
--
--      (a) No arbitrary cut point. The old rule split the window at hour
--          12 and compared two medians; a patient deteriorating from hour
--          14 onward looked identical to one deteriorating from hour 2.
--      (b) It uses every covered hour rather than two summary values, so
--          it neither discards data nor inherits the double-median
--          attenuation rev 2 documented.
--      (c) IT HAS THE SAME MEANING AT EVERY MEASUREMENT DENSITY. For a
--          dense vital with 23 covered hours it is a proper regression;
--          for a lab drawn twice it reduces EXACTLY to
--          (last - first) / (hours between them) * 24. One estimand, one
--          unit, one interpretation, whether the signal was charted 24
--          times or 2. That is what the old definition could not do, and
--          it is why sparse signals previously had no trend at all.
--
--      This is deliberately NOT true of the level terms, and the contrast
--      is worth understanding before anyone tries to unify those too: a
--      slope is a well-defined estimand at any n >= 2 with time spread,
--      whereas a tail quantile's INFORMATION CONTENT is itself a function
--      of n. No reparameterisation fixes the second, and any term that
--      adapts to n is by construction a term whose meaning depends on n.
--
--   3. TREND GUARDS: a span guard AND the half-presence guard.
--      `trend_min_span` (6h between first and last covered hour) is new --
--      without it two draws 1 hour apart extrapolate a 24x slope from a
--      single difference. `trend_min_per_half` (>= 1 covered hour in each
--      half) carries over rev 2's protection, which a span guard alone
--      does NOT provide: hours 0 and 7 satisfy a 6-hour span with nothing
--      after hour 12, and the slope would then be fitted to one end of the
--      window and read as if it described all of it.
--      Neither guard addresses density IMBALANCE -- 2 readings early and
--      10 late passes both, in rev 2 and rev 3 alike. That is an audit-6
--      question, not something a threshold fixes.
--
--      The rev 2 rationale for a coverage guard still applies verbatim and
--      is worth restating: glucose is the clearest case, where a patient
--      started on an insulin protocol gets hourly fingersticks late and
--      two serum draws early -- and the density shift was CAUSED by the
--      intervention being conditioned on.
--
--   4. TREND NO LONGER RESTRICTED TO dense/rate. Sparse labs now get a
--      slope wherever the span guard is satisfied. "Creatinine is rising"
--      is among the most clinically meaningful facts about a sparse lab --
--      AKI is DEFINED by a rise -- and rev 2 excluded it only because the
--      half-split rule could not be computed from two draws.
--      Which classes actually enter a formula stays an R-side config
--      decision (`trend_classes`), so this is extract-wide/model-narrow,
--      not a modelling change made in SQL. RUN AUDIT 5 BEFORE ENABLING IT
--      for a sparse signal: if trend is defined on under ~5% of stays it
--      is a missingness indicator wearing a physiology name.
--
--   5. TREND = 0 WHEN UNDEFINED, unchanged from rev 2 and still correct.
--      Zero is the uninformative value, exactly parallel to pi-hat
--      collapsing to the prior mean at n = 0. An indicator would
--      reintroduce missingness-as-feature, the v1 vulnerability this
--      rebuild exists to remove.
--      COST, also unchanged: it conflates "stable" with "unknown", and
--      stability is genuinely reassuring information. One limitations
--      sentence.
--
--   6. `dx_value_min` / `dx_value_max` RENAMED to `value_min` / `value_max`.
--      They are now modelled for the signals whose measurement density
--      cannot support a quantile, so the `dx_` prefix -- which the R
--      formula builder refuses by assertion -- is no longer correct.
--
--      MEASURED 2026-08-25, MIMIC-IV training set: for every sparse signal
--      `q05` equals the minimum on 100.0% of stays, and `q95` the maximum
--      on 95.6-99.9%; for the three GCS components, 99.6%+ both ways.
--      With 1-2 covered hours the fifth percentile IS the smaller value.
--      Both sets are emitted so the choice stays an R-side declaration,
--      but calling that number a percentile in the methods would be
--      describing a statistic the data cannot support.
--
--      Dense signals keep the quantiles and should: there `q05` is the 5th
--      percentile of hourly MEDIANS -- protected twice -- while
--      `value_min` is MIN(hr_min), the single most extreme raw reading in
--      the window. Step 3b already documents that a lone artifact sets the
--      hourly minimum. A damped arterial line reading MAP 35 passes the
--      physiological bounds in 3a and lands directly in `value_min`.
--
--   7. NAMING CORRECTED to match the R validator CONTRACT. Rev 2 on disk
--      emitted `dx_value_median` and `first_hour_low`, whereas the loader
--      expects `value_median` and `dx_first_hour_low` -- so the file in
--      the repository was NOT the file that produced the parquet, and
--      re-running it would have failed validation on load. Fixed here;
--      the schema below is now the authority.
--
--   8. `dx_trend_span` added, so audit 5 can check the span guard rather
--      than infer it.
--
-- UNCHANGED AND STILL LOAD-BEARING
--   Reference ranges (decision table 1), masking (decision table 2), the
--   Dirichlet-multinomial counts k_low/k_mid/k_high, and the rule that
--   magnitude and persistence are time-symmetric while trend deliberately
--   is not -- coarse directional asymmetry, ONE bit of time ordering, not
--   autocorrelation or lag structure. State that distinction in the
--   methods rather than letting a reviewer find it.
-- =====================================================================


-- =====================================================================
-- TUNABLE
-- =====================================================================
-- Minimum hours between the first and last covered hour before a slope is
-- computed. 6 on a 24-hour window means a slope is never extrapolated more
-- than 4x beyond the interval that produced it. Fixed BEFORE fitting.
DECLARE trend_min_span INT64 DEFAULT 6;

-- Minimum covered hours in EACH half of the window before a slope is computed.
-- Restores the protection rev 2 had and rev 3 initially lost: `trend_min_span`
-- alone is satisfied by hours 0 and 7, with nothing after hour 12, so the slope
-- would be fitted to one end of the window and read as if it described all of
-- it. Be clear about what this does NOT do: neither this nor rev 2 guards
-- against density IMBALANCE -- 2 readings early and 10 late passes both. That
-- remains an audit-6 question, not something a threshold fixes.
DECLARE trend_min_per_half INT64 DEFAULT 1;


-- =====================================================================
-- 5a. SIGNAL FEATURES
-- =====================================================================
CREATE OR REPLACE TABLE
  `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_features_signals_mimiciv` AS

WITH cohort AS (
  SELECT stay_id
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_cohort_mimiciv`
),

hourly AS (
  SELECT stay_id, signal, signal_class, hour_bin, hr_min, hr_med, hr_max,
         n_raw_in_hour
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_mimiciv`
  UNION ALL
  SELECT stay_id, signal, signal_class, hour_bin, hr_min, hr_med, hr_max,
         n_raw_in_hour
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_hourly_gcs_mimiciv`
),

interventions AS (
  SELECT stay_id, intervention, hour_bin
  FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.v2_interventions_mimiciv`
),

-- ---------------------------------------------------------------------
-- DECISION TABLE 1: reference ranges.
--   Labs   -> Kratz et al., N Engl J Med 2004;351:1548-1563 (MGH intervals)
--   Vitals -> APACHE II zero-point bands, Knaus et al., Crit Care Med
--             1985;13:818-829.
--   UNIFIED, not sex-specific, so pi+/- does not silently encode sex.
--   PRE-REGISTERED SENSITIVITY ANALYSES: actionable thresholds and
--   sex-specific ranges - both a one-table edit and a re-run.
-- ---------------------------------------------------------------------
signal_ranges AS (
  SELECT 'heart_rate' AS signal, 'dense' AS signal_class,  70.0 AS ref_low, 109.0 AS ref_high UNION ALL
  SELECT 'mbp',                  'dense',                  70.0,           109.0 UNION ALL
  SELECT 'resp_rate',            'dense',                  12.0,            24.0 UNION ALL
  SELECT 'spo2',                 'dense',                  95.0,           100.0 UNION ALL
  SELECT 'temperature',          'dense',                  36.0,            38.4 UNION ALL
  SELECT 'glucose',              'dense',                  70.0,           110.0 UNION ALL
  SELECT 'gcs_motor',            'dense',                   5.0,             6.0 UNION ALL
  SELECT 'gcs_eyes',             'dense',                   3.0,             4.0 UNION ALL
  SELECT 'gcs_verbal',           'dense',                   3.0,             5.0 UNION ALL
  SELECT 'urine_output_rate',    'rate',                    0.5,             3.0 UNION ALL
  SELECT 'sodium',               'sparse',                135.0,           145.0 UNION ALL
  SELECT 'bicarbonate',          'sparse',                 21.0,            28.0 UNION ALL
  SELECT 'creatinine',           'sparse',                  0.6,             1.5 UNION ALL
  SELECT 'bun',                  'sparse',                 10.0,            20.0 UNION ALL
  SELECT 'wbc',                  'sparse',                  4.5,            11.0 UNION ALL
  SELECT 'hemoglobin',           'sparse',                 12.0,            17.5 UNION ALL
  SELECT 'platelet',             'sparse',                150.0,           400.0 UNION ALL
  SELECT 'lactate',              'sparse',                  0.5,             2.2 UNION ALL
  SELECT 'bilirubin_total',      'sparse',                  0.3,             1.0
),

-- ---------------------------------------------------------------------
-- DECISION TABLE 2: measurement invalidation.
--   'from_onset'  censor hour >= first active hour. RRT: creatinine stays
--                 distorted after dialysis stops within the window - it has
--                 already been cleared by machine.
--   'concurrent'  censor only active hours. GCS recovers once sedation stops.
--   Dexmedetomidine deliberately absent: arousable sedation, GCS assessable.
-- ---------------------------------------------------------------------
mask_spec AS (
  SELECT 'creatinine'        AS signal, 'rrt'      AS masking_intervention, 'from_onset' AS mask_mode UNION ALL
  SELECT 'urine_output_rate',           'rrt',                              'from_onset' UNION ALL
  SELECT 'gcs_motor',                   'paralytic',                        'concurrent' UNION ALL
  SELECT 'gcs_eyes',                    'paralytic',                        'concurrent' UNION ALL
  SELECT 'gcs_verbal',                  'paralytic',                        'concurrent' UNION ALL
  SELECT 'gcs_motor',                   'sedation_benzo',                   'concurrent' UNION ALL
  SELECT 'gcs_eyes',                    'sedation_benzo',                   'concurrent' UNION ALL
  SELECT 'gcs_verbal',                  'sedation_benzo',                   'concurrent' UNION ALL
  SELECT 'gcs_motor',                   'sedation_propofol',                'concurrent' UNION ALL
  SELECT 'gcs_eyes',                    'sedation_propofol',                'concurrent' UNION ALL
  SELECT 'gcs_verbal',                  'sedation_propofol',                'concurrent'
),

intervention_onset AS (
  SELECT stay_id, intervention, MIN(hour_bin) AS first_hour
  FROM interventions
  GROUP BY stay_id, intervention
),

masked_hours AS (
  SELECT DISTINCT h.stay_id, h.signal, h.hour_bin
  FROM hourly h
  JOIN mask_spec m
    ON h.signal = m.signal AND m.mask_mode = 'from_onset'
  JOIN intervention_onset o
    ON h.stay_id = o.stay_id AND o.intervention = m.masking_intervention
  WHERE h.hour_bin >= o.first_hour

  UNION DISTINCT

  SELECT DISTINCT h.stay_id, h.signal, h.hour_bin
  FROM hourly h
  JOIN mask_spec m
    ON h.signal = m.signal AND m.mask_mode = 'concurrent'
  JOIN interventions i
    ON h.stay_id = i.stay_id
   AND i.intervention = m.masking_intervention
   AND i.hour_bin = h.hour_bin
),

surviving AS (
  SELECT h.*
  FROM hourly h
  LEFT JOIN masked_hours mh
    ON h.stay_id = mh.stay_id AND h.signal = mh.signal AND h.hour_bin = mh.hour_bin
  WHERE mh.stay_id IS NULL
),

mask_counts AS (
  SELECT stay_id, signal, COUNT(*) AS n_masked
  FROM masked_hours
  GROUP BY stay_id, signal
),

grid AS (
  SELECT c.stay_id, r.signal, r.signal_class, r.ref_low, r.ref_high
  FROM cohort c CROSS JOIN signal_ranges r
),

agg AS (
  SELECT
    g.stay_id,
    g.signal,
    g.signal_class,
    g.ref_low,
    g.ref_high,

    -- ---- denominator (post-masking) ---------------------------------
    COUNT(s.hour_bin)                                         AS n_obs,
    COALESCE(ANY_VALUE(mc.n_masked), 0)                       AS dx_n_masked,

    -- ---- counts for Dirichlet-multinomial shrinkage in R -------------
    COUNTIF(s.hr_med <  g.ref_low)                            AS k_low,
    COUNTIF(s.hr_med >= g.ref_low AND s.hr_med <= g.ref_high) AS k_mid,
    COUNTIF(s.hr_med >  g.ref_high)                           AS k_high,

    -- ---- magnitude: BOTH parameterisations, always -------------------
    -- Quantiles of hourly MEDIANS. Protected twice, and therefore the
    -- right choice wherever the density supports them.
    APPROX_QUANTILES(s.hr_med, 20)[OFFSET(1)]                 AS q05,
    APPROX_QUANTILES(s.hr_med, 20)[OFFSET(19)]                AS q95,
    APPROX_QUANTILES(s.hr_med,  2)[OFFSET(1)]                 AS value_median,

    -- Extremes of raw readings. Honest for sparse signals and for the
    -- bounded GCS scales, where the quantiles above ARE these numbers.
    -- Artifact-exposed for dense signals: see header note 6.
    MIN(s.hr_min)                                             AS value_min,
    MAX(s.hr_max)                                             AS value_max,

    -- ---- replicate structure for the s_e variance component -----------
    -- ADDED 2026-09-02. THESE ARE NOT MODEL COVARIATES and none of them
    -- may ever enter a formula. `^se_` is refused by the formula builder
    -- for the same reason `^dx_` and `^qc_` are: they describe how a stay
    -- was MEASURED, not what was measured.
    --
    -- WHY THEY EXIST. `delta` shrinks its residual by
    -- s_u / (s_u + s_e / n), where s_e is the within-stay noise
    -- component. Until now s_e was INFERRED from how residual variance
    -- changes with n_obs, using one residual per stay -- and that is not
    -- identified. MEASURED 2026-09-02 by re-running
    -- `.fit_var_components()` on the real residuals for all 38 delta
    -- rows: on 19 of them the two-component fit returned s_e between 1e-9
    -- and 1e-6 with a likelihood IDENTICAL to the one-component model
    -- (aic2 - aic1 = +2.000 exactly, the parameter penalty and nothing
    -- else). The shrinkage weight was therefore 1.000 at every n and the
    -- machinery was inert -- including on all three Glasgow components,
    -- which is precisely where a between-site coverage difference bites.
    -- No selection criterion could have fixed that: the likelihood is
    -- flat, so BIC, an LRT and cross-validation all tie as well.
    --
    -- The repair is to give the estimator a REPLICATE STRUCTURE. The
    -- hourly lattice already carries replicates at two levels and both
    -- are summarised here as SUFFICIENT STATISTICS rather than as
    -- estimates, so R can pool them exactly across stays:
    --
    --     s_e = SUM(se_within_ss) / SUM(se_within_n - 1)
    --
    -- over stays with se_within_n >= 2. That is the within-group mean
    -- square of a one-way random-effects model, which is what s_e was
    -- always meant to be, and it is estimable whether or not the residual
    -- variance happens to vary with n.
    --
    -- TRAINING-SITE ONLY, AND THAT IS THE POINT. s_e is fitted at MIMIC
    -- and frozen into the bundle; an apply site computes its weight from
    -- the FROZEN s_e and its own n_obs. So these columns add no transport
    -- dependency. They are extracted at both sites because hard rule 5
    -- requires byte-identical schemas, not because eICU needs them in
    -- order to be scored.
    --
    -- VAR_POP(x) * COUNT(x) is the exact sum of squared deviations and is
    -- numerically better behaved than SUM(x*x) - n*AVG(x)^2. VAR_POP,
    -- AVG and COUNTIF all ignore a NULL hr_med, so all three agree on the
    -- same rows. Both are NULL/0 at n_obs = 1, which is correct: a single
    -- hour carries no within-stay information and contributes no degrees
    -- of freedom to the pool.
    COUNTIF(s.hr_med IS NOT NULL)                             AS se_within_n,
    AVG(s.hr_med)                                             AS se_within_mean,
    VAR_POP(s.hr_med) * COUNTIF(s.hr_med IS NOT NULL)         AS se_within_ss,

    -- THE SECOND REPLICATE LEVEL, and the one that matters for the
    -- extreme statistics. `value_min` is MIN(hr_min) over RAW readings,
    -- so its sampling variance depends on within-HOUR spread as well as
    -- on between-hour spread. This is the mean within-hour range; it is
    -- exactly 0 wherever every hour holds a single reading, which is
    -- itself the diagnostic for whether the second level exists at all.
    AVG(s.hr_max - s.hr_min)                                  AS se_hour_range_mean,

    -- THE TRUE REPLICATE COUNT. `n_obs` counts HOURS; this counts
    -- READINGS, and the two differ by about 1.2x at MIMIC and far more on
    -- eICU's 5-minute vitals feed. Carried as a diagnostic and as the
    -- input to any later decision about whether the weight in s_e / w
    -- should be hours or readings. IT IS NOT THAT WEIGHT TODAY: changing
    -- w is a design change and would not be made in an extraction file.
    SUM(s.n_raw_in_hour)                                      AS se_n_raw,

    -- ---- trajectory: OLS slope inputs --------------------------------
    -- COVAR_POP and VAR_POP ignore rows where either argument is NULL, so
    -- the LEFT JOIN's unmatched rows drop out on their own. Both computed
    -- on POST-MASKING surviving hours.
    COVAR_POP(s.hour_bin, s.hr_med)                           AS trend_cov,
    VAR_POP(s.hour_bin)                                       AS trend_var,

    -- ---- coverage diagnostics, and the half-balance guard -------------
    -- These are trend inputs again as of the `trend_min_per_half` guard, and
    -- remain the way to detect the coverage asymmetry of header note 3.
    COUNTIF(s.hour_bin <  12)                                 AS dx_n_hours_early,
    COUNTIF(s.hour_bin >= 12)                                 AS dx_n_hours_late,

    -- ---- timing diagnostics ------------------------------------------
    -- 5c consumed first_hour_low/high; 5c is gone. Retained as `dx_`
    -- because they remain the cheapest way to audit excursion timing, and
    -- the `dx_` prefix is what stops them reaching a formula.
    MIN(IF(s.hr_med <  g.ref_low,  s.hour_bin, NULL))         AS dx_first_hour_low,
    MIN(IF(s.hr_med >  g.ref_high, s.hour_bin, NULL))         AS dx_first_hour_high,
    MIN(s.hour_bin)                                           AS dx_first_hour_obs,
    MAX(s.hour_bin)                                           AS dx_last_hour_obs

  FROM grid g
  LEFT JOIN surviving s
    ON g.stay_id = s.stay_id AND g.signal = s.signal
  LEFT JOIN mask_counts mc
    ON g.stay_id = mc.stay_id AND g.signal = mc.signal
  GROUP BY g.stay_id, g.signal, g.signal_class, g.ref_low, g.ref_high
),

slope AS (
  SELECT
    a.*,
    COALESCE(a.dx_last_hour_obs - a.dx_first_hour_obs, 0)     AS dx_trend_span,
    -- Defined only where a slope is estimable AND the span is long enough
    -- that it is not extrapolating a single difference across the window.
    (a.n_obs >= 2
      AND a.trend_var IS NOT NULL AND a.trend_var > 0
      AND (a.dx_last_hour_obs - a.dx_first_hour_obs) >= trend_min_span
      AND a.dx_n_hours_early >= trend_min_per_half
      AND a.dx_n_hours_late  >= trend_min_per_half) AS trend_ok
  FROM agg a
)

SELECT
  s.* EXCEPT (trend_cov, trend_var, trend_ok),

  -- TREND: change per 24 hours, positive = rising across the window.
  -- Zero when undefined (header note 5), NOT a flagged special case.
  CASE
    WHEN s.trend_ok THEN s.trend_cov / s.trend_var * 24.0
    ELSE 0.0
  END                                                         AS trend,

  -- DIAGNOSTIC ONLY. Do NOT feed this to the GAM - it is a missingness
  -- indicator, and scoring it would reintroduce exactly the v1
  -- vulnerability this design removes. Use it to check whether trend
  -- availability correlates with outcome; if it does, that is a
  -- limitation to name, not a feature to add.
  CAST(s.trend_ok AS INT64)                                   AS dx_trend_defined,

  CASE WHEN s.n_obs = 0 THEN 1 ELSE 0 END                     AS dx_value_missing

FROM slope s;


-- =====================================================================
-- 5c. ORDERING FLAG -- DELETED
--
-- `v2_features_ordering_mimiciv` is no longer produced. See header note 1.
-- Do not restore it without reading CLAUDE.md's frozen decision on
-- `o_flag` first: the objection is that it states treatment sequencing in
-- a design that makes no causal claim, and that objection is not fixed by
-- renaming the levels.
--
-- The idea considered and rejected alongside it -- extracting
-- `t(first intervention) - t(first observed deviation)` as a continuous
-- covariate -- is undefined for patients who never deviate, requires
-- reconciling per-signal and per-intervention time resolutions, and would
-- inject a time-to-event structure into a joint-probability design.
-- =====================================================================


-- =====================================================================
-- AUDITS
-- =====================================================================
-- 1. MASKING IMPACT. How much did invalidation remove? Report this.
-- SELECT signal, AVG(dx_n_masked) AS mean_hours_masked,
--        COUNTIF(dx_n_masked > 0) AS n_stays_affected,
--        COUNTIF(n_obs = 0 AND dx_n_masked > 0) AS n_fully_masked
-- FROM `...v2_features_signals_mimiciv`
-- WHERE dx_n_masked > 0 GROUP BY signal;
--
-- 2. EXCURSION BASE RATES. If k_mid/n_obs is ~0 or ~1 for a signal, the
--    reference range is not discriminating. Glucose 70-110 is the most
--    likely offender - fed ICU patients sit at 110-140 without pathology.
-- SELECT signal, SAFE_DIVIDE(SUM(k_low), SUM(n_obs))  AS rate_low,
--                SAFE_DIVIDE(SUM(k_mid), SUM(n_obs))  AS rate_mid,
--                SAFE_DIVIDE(SUM(k_high), SUM(n_obs)) AS rate_high
-- FROM `...v2_features_signals_mimiciv` WHERE n_obs > 0 GROUP BY signal;
--
-- 3. QUANTILE vs EXTREME AGREEMENT. The measurement behind header note 6,
--    re-run at each site. Where these are ~100%, the quantile IS the
--    extreme and must be named as one.
-- SELECT signal, ANY_VALUE(signal_class) AS signal_class,
--        APPROX_QUANTILES(n_obs, 2)[OFFSET(1)]        AS median_n_obs,
--        AVG(CAST(q05 = value_min AS INT64))          AS frac_q05_is_min,
--        AVG(CAST(q95 = value_max AS INT64))          AS frac_q95_is_max
-- FROM `...v2_features_signals_mimiciv` WHERE n_obs > 0 GROUP BY signal;
--
-- 4. INTERVENTION PREVALENCE sanity: vasopressor ~25-35%, invasive vent
--    ~40%, RRT ~5-10%. Large deviations mean label matching missed variants.
-- SELECT intervention, AVG(ever_active) AS prevalence,
--        AVG(IF(ever_active = 1, exposure_frac, NULL)) AS mean_exposure_if_active
-- FROM `...v2_intervention_features_mimiciv` GROUP BY intervention ORDER BY 2 DESC;
--
-- 5. TREND AVAILABILITY -- RUN THIS BEFORE ENABLING `trend` FOR A SPARSE
--    SIGNAL (header note 4). Under ~5% defined means the term is a
--    missingness indicator wearing a physiology name; leave that signal
--    out of `trend_classes` in config/config.yml.
-- SELECT signal, ANY_VALUE(signal_class) AS signal_class,
--        AVG(dx_trend_defined)                        AS frac_defined,
--        APPROX_QUANTILES(dx_trend_span, 2)[OFFSET(1)] AS median_span,
--        COUNTIF(n_obs >= 2 AND dx_trend_defined = 0) AS n_blocked_by_span
-- FROM `...v2_features_signals_mimiciv` WHERE n_obs > 0 GROUP BY signal ORDER BY 2, 1;
--
-- 6. INFORMATIVE TREND MISSINGNESS. Does trend availability track outcome?
--    If mortality differs sharply between dx_trend_defined 0 and 1 for a
--    dense signal, coverage asymmetry is severity-related -> limitation.
-- SELECT f.signal, f.dx_trend_defined, COUNT(*) AS n, AVG(c.mortality) AS mort
-- FROM `...v2_features_signals_mimiciv` f
-- JOIN `...v2_cohort_mimiciv` c USING (stay_id)
-- WHERE f.n_obs > 0 GROUP BY 1,2 ORDER BY 1,2;
--
-- 7. TREND DISTRIBUTION. Should be centred near 0 and roughly symmetric.
--    The spike at exactly 0 is the undefined cases - cross-check its
--    height against (1 - frac_defined) from audit 5. A slope has a longer
--    tail than the old half-difference, so check for implausible values:
--    a creatinine slope of 20 mg/dL/24h is a span-guard failure, not AKI.
-- SELECT signal, APPROX_QUANTILES(trend, 20) AS ventiles,
--        COUNTIF(trend = 0) AS n_exactly_zero,
--        MIN(trend) AS min_trend, MAX(trend) AS max_trend
-- FROM `...v2_features_signals_mimiciv` WHERE n_obs > 0 GROUP BY signal;
--
-- 8. TREND vs PERSISTENCE overlap. Correlate trend with k_low/n_obs and
--    k_high/n_obs per signal. Partial overlap is expected (a deteriorating
--    patient tends to have both); near-collinearity is not.
-- SELECT signal,
--        CORR(trend, SAFE_DIVIDE(k_low,  n_obs)) AS r_trend_klow,
--        CORR(trend, SAFE_DIVIDE(k_high, n_obs)) AS r_trend_khigh
-- FROM `...v2_features_signals_mimiciv` WHERE n_obs > 0 GROUP BY signal;
