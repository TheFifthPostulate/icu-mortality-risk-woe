-- =====================================================================
-- v2 STEP 5: aggregation, eICU
--
-- PORT OF sql/mimiciv/v2_05_features_mimiciv_rev3.sql.
--
-- THIS FILE IS DELIBERATELY A NEAR-VERBATIM COPY. Every decision it
-- encodes - reference ranges, masking rules and modes, the
-- Dirichlet-multinomial counts, the OLS trend definition, both trend
-- guards, trend = 0 when undefined, the dual magnitude parameterisation,
-- the dx_ prefixing - is a MODELLING DECISION frozen from MIMIC. Hard
-- rule 8 says eICU applies frozen quantities and re-derives nothing;
-- this is the SQL-side expression of that rule.
--
-- The ONLY differences from the MIMIC file are the four table names and
-- the header comments. If a diff between the two files ever shows a
-- change to `signal_ranges`, `mask_spec`, `trend_min_span`,
-- `trend_min_per_half` or the aggregation expressions, that is a bug,
-- not a port. Diff them before running.
--
-- Read the MIMIC file's header for the rev-3 rationale (OLS trend, span
-- and half-presence guards, value_min/value_max naming, the deleted
-- ordering flag). None of it is restated here, because restating it
-- invites the two copies to drift.
--
-- ---------------------------------------------------------------------
-- WHAT TO EXPECT TO BE DIFFERENT IN THE OUTPUT, and why none of it is a
-- reason to change this file
-- ---------------------------------------------------------------------
--   n_obs        Dense signals will sit near 24 rather than MIMIC's ~22,
--                because eICU charts vitals on a 5-minute interface. That
--                is the site difference the design handles inside pi_hat
--                rather than by adjusting the extraction; `n_obs` never
--                enters a formula precisely so that this cannot transport
--                as an artifact.
--   k_low/mid/high  Counted on `hr_med`, so the 5-minute interface does
--                not inflate excursion counts through single artifacts.
--                It DOES mean an excursion lasting ten minutes is
--                invisible at both sites, which is the intended symmetry.
--   dx_n_masked  Higher at eICU for creatinine and urine output than the
--                MIMIC-equivalent, because `rrt` is reconstructed
--                onset-to-window-end (step 2b, rule R3) and the mask mode
--                for those two signals is 'from_onset' anyway. The two
--                rules agree exactly here, which is why R3 is acceptable
--                for rrt and questionable elsewhere.
--   trend        The definable fraction differs sharply for GCS (0.497 vs
--                0.85, see step 04's header). `trend_classes` in config
--                admits dense and rate signals, so the GCS components DO
--                carry a trend term, and this is a live transport threat
--                that is documented rather than patched.
--   frac q05 = value_min  Re-run audit 3 here. At MIMIC it is 100% for
--                every sparse signal, which is what justified calling
--                those numbers extremes rather than percentiles. eICU's
--                lab density is similar, so it should hold - but "should"
--                is not a measurement.
-- =====================================================================


-- =====================================================================
-- TUNABLE - VERBATIM FROM THE MIMIC FILE. Frozen before fitting; these
-- are not site parameters.
-- =====================================================================
DECLARE trend_min_span INT64 DEFAULT 6;
DECLARE trend_min_per_half INT64 DEFAULT 1;


-- =====================================================================
-- 5a. SIGNAL FEATURES
-- =====================================================================
CREATE OR REPLACE TABLE
  `eicu-ext.eicu_ext_data.v2_features_signals_eicu` AS

WITH cohort AS (
  SELECT stay_id
  FROM `eicu-ext.eicu_ext_data.v2_cohort_eicu`
),

hourly AS (
  SELECT stay_id, signal, signal_class, hour_bin, hr_min, hr_med, hr_max,
         n_raw_in_hour
  FROM `eicu-ext.eicu_ext_data.v2_hourly_eicu`
  UNION ALL
  SELECT stay_id, signal, signal_class, hour_bin, hr_min, hr_med, hr_max,
         n_raw_in_hour
  FROM `eicu-ext.eicu_ext_data.v2_hourly_gcs_eicu`
),

interventions AS (
  SELECT stay_id, intervention, hour_bin
  FROM `eicu-ext.eicu_ext_data.v2_interventions_eicu`
),

-- ---------------------------------------------------------------------
-- DECISION TABLE 1: reference ranges. VERBATIM.
--   Labs   -> Kratz et al., N Engl J Med 2004;351:1548-1563 (MGH intervals)
--   Vitals -> APACHE II zero-point bands, Knaus et al., Crit Care Med
--             1985;13:818-829.
--   UNIFIED, not sex-specific, so pi+/- does not silently encode sex.
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
-- DECISION TABLE 2: measurement invalidation. VERBATIM.
--   'from_onset'  censor hour >= first active hour.
--   'concurrent'  censor only active hours.
--   Dexmedetomidine deliberately absent: arousable sedation.
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
    APPROX_QUANTILES(s.hr_med, 20)[OFFSET(1)]                 AS q05,
    APPROX_QUANTILES(s.hr_med, 20)[OFFSET(19)]                AS q95,
    APPROX_QUANTILES(s.hr_med,  2)[OFFSET(1)]                 AS value_median,

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
    COVAR_POP(s.hour_bin, s.hr_med)                           AS trend_cov,
    VAR_POP(s.hour_bin)                                       AS trend_var,

    -- ---- coverage diagnostics, and the half-balance guard -------------
    COUNTIF(s.hour_bin <  12)                                 AS dx_n_hours_early,
    COUNTIF(s.hour_bin >= 12)                                 AS dx_n_hours_late,

    -- ---- timing diagnostics ------------------------------------------
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
    (a.n_obs >= 2
      AND a.trend_var IS NOT NULL AND a.trend_var > 0
      AND (a.dx_last_hour_obs - a.dx_first_hour_obs) >= trend_min_span
      AND a.dx_n_hours_early >= trend_min_per_half
      AND a.dx_n_hours_late  >= trend_min_per_half) AS trend_ok
  FROM agg a
)

SELECT
  s.* EXCEPT (trend_cov, trend_var, trend_ok),

  CASE
    WHEN s.trend_ok THEN s.trend_cov / s.trend_var * 24.0
    ELSE 0.0
  END                                                         AS trend,

  CAST(s.trend_ok AS INT64)                                   AS dx_trend_defined,

  CASE WHEN s.n_obs = 0 THEN 1 ELSE 0 END                     AS dx_value_missing

FROM slope s;


-- =====================================================================
-- AUDITS
--
-- Every one of these is the MIMIC audit re-run at eICU. The point is not
-- the eICU number on its own; it is the PAIR. Put both columns in the
-- same table in the supplement.
-- =====================================================================
-- 1. MASKING IMPACT. Report this. eICU's rrt reconstruction is
--    onset-to-window-end, so creatinine and urine_output_rate masking
--    should be at least as aggressive as MIMIC's.
-- SELECT signal, AVG(dx_n_masked) AS mean_hours_masked,
--        COUNTIF(dx_n_masked > 0) AS n_stays_affected,
--        COUNTIF(n_obs = 0 AND dx_n_masked > 0) AS n_fully_masked
-- FROM `...v2_features_signals_eicu` WHERE dx_n_masked > 0 GROUP BY signal;
--
-- 2. EXCURSION BASE RATES. If k_mid/n_obs is ~0 or ~1 for a signal, the
--    reference range is not discriminating AT THIS SITE. Glucose 70-110
--    is the known offender; the eICU number decides whether the frozen
--    range transports or whether the pre-registered sensitivity analysis
--    on the glucose range becomes load-bearing.
-- SELECT signal, SAFE_DIVIDE(SUM(k_low), SUM(n_obs))  AS rate_low,
--                SAFE_DIVIDE(SUM(k_mid), SUM(n_obs))  AS rate_mid,
--                SAFE_DIVIDE(SUM(k_high), SUM(n_obs)) AS rate_high
-- FROM `...v2_features_signals_eicu` WHERE n_obs > 0 GROUP BY signal;
--
-- 3. QUANTILE vs EXTREME AGREEMENT. At MIMIC, q05 = value_min on 100.0%
--    of stays for every sparse signal and 99.6%+ for the GCS triple,
--    which is what justifies config/level_terms_by_class giving those
--    signals `extreme`. If eICU disagrees materially, the same column
--    name means two different statistics at the two sites.
-- SELECT signal, ANY_VALUE(signal_class) AS signal_class,
--        APPROX_QUANTILES(n_obs, 2)[OFFSET(1)]        AS median_n_obs,
--        AVG(CAST(q05 = value_min AS INT64))          AS frac_q05_is_min,
--        AVG(CAST(q95 = value_max AS INT64))          AS frac_q95_is_max
-- FROM `...v2_features_signals_eicu` WHERE n_obs > 0 GROUP BY signal;
--
-- 4. UNMEASURED FRACTION. MIMIC's is 8.7% over the signal grid, and the
--    design assigns L = 0 there. If eICU's is much larger the assignment
--    is doing more work at the validation site than at development, which
--    is a statement about the result and belongs in the abstract, not the
--    limitations.
-- SELECT signal, AVG(dx_value_missing) AS frac_unmeasured
-- FROM `...v2_features_signals_eicu` GROUP BY signal ORDER BY 2 DESC;
--
-- 5. TREND AVAILABILITY, per signal, against MIMIC (audit 02_02b).
-- SELECT signal, ANY_VALUE(signal_class) AS signal_class,
--        AVG(dx_trend_defined)                         AS frac_defined,
--        APPROX_QUANTILES(dx_trend_span, 2)[OFFSET(1)] AS median_span,
--        COUNTIF(n_obs >= 2 AND dx_trend_defined = 0)  AS n_blocked_by_span
-- FROM `...v2_features_signals_eicu` WHERE n_obs > 0 GROUP BY signal ORDER BY 2, 1;
--
-- 6. INFORMATIVE TREND MISSINGNESS. Does trend availability track
--    outcome? config/trend_classes excludes sparse signals at MIMIC for
--    exactly this reason; re-run the measurement here rather than
--    assuming the exclusion still covers it.
-- SELECT f.signal, f.dx_trend_defined, COUNT(*) AS n, AVG(c.mortality) AS mort
-- FROM `...v2_features_signals_eicu` f
-- JOIN `...v2_cohort_eicu` c USING (stay_id)
-- WHERE f.n_obs > 0 GROUP BY 1,2 ORDER BY 1,2;
--
-- 7. TREND DISTRIBUTION. Centred near 0, roughly symmetric, with a spike
--    at exactly 0 whose height matches (1 - frac_defined) from audit 5.
-- SELECT signal, APPROX_QUANTILES(trend, 20) AS ventiles,
--        COUNTIF(trend = 0) AS n_exactly_zero,
--        MIN(trend) AS min_trend, MAX(trend) AS max_trend
-- FROM `...v2_features_signals_eicu` WHERE n_obs > 0 GROUP BY signal;
--
-- 8. TREND vs PERSISTENCE overlap.
-- SELECT signal,
--        CORR(trend, SAFE_DIVIDE(k_low,  n_obs)) AS r_trend_klow,
--        CORR(trend, SAFE_DIVIDE(k_high, n_obs)) AS r_trend_khigh
-- FROM `...v2_features_signals_eicu` WHERE n_obs > 0 GROUP BY signal;
--
-- 9. SCHEMA PARITY, the cheapest and most load-bearing check in the file.
--    The column set and types must match v2_features_signals_mimiciv
--    exactly, or R/02_validate.R check 9 stops the external run.
-- SELECT column_name, data_type
-- FROM `eicu-ext.eicu_ext_data.INFORMATION_SCHEMA.COLUMNS`
-- WHERE table_name = 'v2_features_signals_eicu' ORDER BY ordinal_position;
--    Diff against the same query on the MIMIC dataset.
