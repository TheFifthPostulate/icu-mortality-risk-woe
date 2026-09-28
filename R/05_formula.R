# R/05_formula.R -------------------------------------------------------------
# Formula builder, plus the ^dx_ / ^qc_ / ^se_ refusal guard.
#
# Three models per signal (v2_analytical_design_plan.md, v2_handoff.md SS7):
#   meas - measurement variables only        -> physiological evidence
#   intv - the paired interventions only     -> treatment-context evidence
#   full - meas + intv                       -> the primary joint model
#
# The term sets of `meas` and `intv` are DISJOINT and their union is exactly
# `full`. That is what makes the subtraction well defined:
#
#   L_cond = L_full - L_intv   the measurement-deviation LLR conditional on
#                              class AND intervention context
#   L_full - L_meas            the incremental treatment evidence
#
# The partition is exact only because `o_flag` is gone (dropped 2026-08-25). It
# was the one term belonging to neither set -- its levels are defined relative
# to the first excursion, so it encoded measurement information under an
# intervention name, and it made the decomposition arguable rather than
# arithmetic. See CLAUDE.md, frozen decisions.
#
# For the 7 unpaired signals `meas` and `full` are identical by construction and
# `intv` does not exist at all: with no cognate intervention its formula would
# be `mortality ~ 1`, so L_intv is 0 by construction and is ASSIGNED, never
# fitted. build_formula() refuses it rather than emitting an intercept-only
# model that would look like a fit.
#
# No paths, no clock, no data access (hard rule 9). Everything comes from cfg
# and pairing.csv, so a formula can be inspected before a single row is read.
# ----------------------------------------------------------------------------

# `se_` joined 2026-09-02 with the replicate-structure columns. They are
# sufficient statistics for the `delta` within-stay variance component and
# describe HOW a stay was measured rather than WHAT was measured, so they
# are exactly the class this guard exists for. Closing the class rather than
# the instance is the point: a sixth `se_` column added later is refused
# without anyone remembering to add it here.
FORBIDDEN_PREFIX <- "^(dx_|qc_|se_)"

#' Variables that may never enter a formula, by NAME rather than by prefix.
#'
#' The prefix guard cannot reach these: they arrive as `{intervention}__first_hour`,
#' so `^dx_` never matches. Each was previously excluded by omission or by a
#' comment, which is not a guard — a later edit to intervention_terms() would
#' have re-admitted any of them silently.
#'
#'   first_hour             TIME-TO-EVENT. Hours from admission to first
#'                          exposure. It is exactly the class of variable
#'                          rejected with `o_flag`: it injects a time-to-event
#'                          structure into a joint-probability design, it is
#'                          undefined for the unexposed, and MIMIC and eICU
#'                          resolve intervention times differently, so it would
#'                          transport as an artifact.
#'   peak_intensity         Nullable and NULL for most rows (spec §7 defers it).
#'                          Under na.fail it would delete every unexposed patient
#'                          from the fit.
#'   max_concurrent_agents  Equal to `n_agents` on 99.03% of stays, r = 0.9916.
#'                          Dropped 2026-08-25.
#'   n_obs                  Monitoring intensity, site-specific. Enters only
#'                          inside pi_hat, as the confidence weighting.
#'   pi_mid                 Determined by pi_minus and pi_plus; entering all
#'                          three is exactly collinear.
#'
#'   ever_active            Exactly `(intensity > 0)` on 100.0000% of stays, so
#'                          never an independent covariate. Adds <= 0.0019
#'                          deviance explained on top of the intensity smooth.
#'                          Dropped 2026-08-25.
#'
#' These are assertions, like FORBIDDEN_PREFIX. Adding one back is a spec change
#' to discuss, not a line to delete.
FORBIDDEN_VARS <- c("first_hour", "peak_intensity", "max_concurrent_agents",
                    "ever_active", "n_obs", "pi_mid", "dx_nee_peak")

#' Variables DISPLACED by an active conditional-prior construct.
#'
#' Config-dependent, and therefore separate from FORBIDDEN_VARS, which is
#' static. When `magnitude_conditional` is on, `q05` has been replaced by
#' `q05_delta`; when `intensity_conditional` is on, `n_agents` has been replaced
#' by `{iv}__lambda`. In both cases the raw term must not ALSO appear -- the
#' delta and the lambda are standardised against exactly the quantity the raw
#' term carries, so entering both would re-admit the collinearity the construct
#' exists to remove, and it would do so quietly.
#'
#' This is the guard, not the comment above intervention_terms(). Turning a
#' construct off restores its raw term automatically, in both places, because
#' both read the same config flag.
#' Evaluated PER SIGNAL rather than globally, so a signal held out of the
#' magnitude construct by `magnitude_conditional_override` keeps its raw term
#' without weakening the guard for every other signal.
.displaced_vars <- function(signal, cfg) {
  if (is.null(cfg) || is.null(signal)) return(character(0))
  out <- character(0)
  if (magnitude_conditional_for(signal, cfg)) {
    # BOTH members of the pair, not just the one `excursion_side` selects: the
    # delta standardises the same underlying quantity either way.
    out <- c(out, level_vars_of(signal, cfg))
  }
  for (iv in interventions_of(signal, cfg)) {
    ls <- lambda_spec_of(iv, cfg)
    if (!is.null(ls)) out <- c(out, paste0(iv, "__", ls$accumulation))
  }
  unique(out)
}

#' The layer-1 model set. A design constant, not a config knob: adding a model
#' changes what layer 2 aggregates, which is a spec change to discuss.
#'
#' THE TWO INTERACTION MODELS JOINED IT ON 2026-09-07, and that reopens a frozen
#' decision through its own escape clause rather than against it. CLAUDE.md's
#' "Layer 1 carries no interaction terms at all" rejected `s(trend, by = o_flag)`
#' on two grounds: it was built on EXCURSION-RELATIVE TIMING, which is the part
#' that was objectionable, and it drove worst-case concurvity to about 1 in every
#' paired dense model. The clause that closed it reads "if an interaction is ever
#' wanted back, it must be built on INTERVENTION INTENSITY rather than on
#' excursion-relative timing", and that is exactly what these are:
#' `ti(measurement covariate, intervention intensity)`, with `o_flag` nowhere in
#' sight. The concurvity objection is separate and is MEASURED rather than
#' argued -- `ti()` excludes the main effects a `by=` smooth re-carries, so it
#' should not reproduce that number, and R/08 reports `concurvity_max` for every
#' one of these fits beside every other fit's.
#'
#'   full_ti_trend  `full` plus `ti(trend, <each intervention intensity>)`.
#'                  Treatment response as an explicit interaction, and the one
#'                  three independent instruments converged on.
#'   full_ti_all    `full` plus every measurement x intervention cross term.
#'                  The saturated version, kept so the trend-only result has an
#'                  upper reference rather than only a lower one.
#'
#' `full` REMAINS PRIMARY. These are additional models, not replacements: every
#' claim the paper makes about `L_full` and `L_cond` is unchanged, and the
#' interaction arms exist to bound what those claims cost.
LAYER1_MODELS <- c("meas", "full", "intv", "full_ti_trend", "full_ti_all")

#' The interaction models, as a set. Used wherever a caller has to say "the ti
#' arms" without restating the two names, which is how a third one gets missed.
LAYER1_TI_MODELS <- c("full_ti_trend", "full_ti_all")

#' Which models exist for a signal. `intv` only where there is an intervention.
models_of <- function(signal, cfg) {
  if (length(interventions_of(signal, cfg))) LAYER1_MODELS else c("meas", "full")
}

#' Where a (signal, model) cell's L COMES FROM: a fit, another model, or zero.
#'
#' THE ONE PLACE THE ALIAS AND ASSIGNMENT RULES ARE DECIDED. `layer1_jobs()`
#' enumerated them and `l_matrix()` re-implemented them, which was survivable
#' with two rules and three models and is not with five: a ti model aliases
#' `meas` on an unpaired signal, aliases `full` on a paired signal with no cross
#' terms, and is fitted otherwise, and getting that wrong in one of the two
#' places pivots one model's column under another model's name.
#'
#' @return NA_character_ where the cell is FITTED; otherwise the model whose
#'   predictions it takes, or "zero" for an assigned cell.
spec_source <- function(signal, model, cfg) {
  paired <- length(interventions_of(signal, cfg)) > 0L
  if (!paired) {
    if (model == "intv") return("zero")                # assigned, never fitted
    if (model != "meas") return("meas")                # `full` and both ti arms
    return(NA_character_)
  }
  # Paired, but the cross-term set can still be empty: `full_ti_trend` on a
  # signal whose class carries no `trend` covariate (creatinine, platelet,
  # hemoglobin) has nothing to cross, so its formula IS `full`'s.
  if (model %in% LAYER1_TI_MODELS &&
      !length(interaction_terms(signal, cfg, .ti_scope_of(model)))) return("full")
  NA_character_
}

#' Follow the alias chain from a requested model to the spec that is FITTED.
#'
#' Returns the fitted model name (`model` itself when it is fitted, `meas` for
#' an unpaired signal's `full`, ...) or the string "zero" for an assigned cell.
#' An alias can chain, and nothing in `spec_source()` promises a one-step hop,
#' so the walk is bounded and a cycle is an error rather than a hang.
#'
#' ONE RESOLVER (added 2026-09-09). `l_matrix()` carried this loop inline; the
#' nested cross-fit for the stacked `xgb_l` cell (`layer1_nested_l()`) needs
#' the same answer, and two copies of an alias walk are two places for the
#' rule to drift.
resolve_spec_source <- function(signal, model, cfg, max_hops = 3L) {
  src <- model
  for (hop in 0:max_hops) {
    nx <- spec_source(signal, src, cfg)
    if (is.na(nx)) return(src)
    if (identical(nx, "zero")) return("zero")
    src <- nx
  }
  stop(sprintf("resolve_spec_source: the alias chain for '%s'/'%s' did not ",
               "resolve in %d hops. spec_source() has a cycle.",
               signal, model, max_hops), call. = FALSE)
}

#' The guard. An assertion, not a convention (hard rule 4).
#'
#' `dx_` marks diagnostics that must never be modelled. `qc_` is in the pattern
#' because `qc_discharge_hospice` keeps its extracted name and is structurally
#' unavailable in eICU -- modelling it would produce a transportability result
#' that is really a plumbing artifact.
#'
#' Never work around this. If such a column is genuinely needed, that is a spec
#' change to discuss first.
assert_no_forbidden <- function(terms, cfg = NULL, signal = NULL) {
  v <- .term_vars(terms)
  bad <- v[grepl(FORBIDDEN_PREFIX, v)]
  if (length(bad)) {
    abort_values(paste0("formula builder refuses terms matching ", FORBIDDEN_PREFIX,
                        " (hard rule 4)"), bad)
  }

  # Checked on the BARE name and on whatever follows `{intervention}__`, because
  # a forbidden variable reaches a formula wearing an intervention prefix:
  # `vasopressor__first_hour`, not `first_hour`.
  base <- sub("^.*__", "", v)
  hit <- v[base %in% FORBIDDEN_VARS]
  if (length(hit)) {
    abort_values("formula builder refuses a forbidden variable (see FORBIDDEN_VARS)", hit)
  }

  # Config-dependent half: a raw term whose conditional replacement is active.
  disp <- .displaced_vars(signal, cfg)
  hit <- v[v %in% disp]
  if (length(hit)) {
    abort_values(paste0("formula builder refuses a raw term displaced by an active ",
                        "conditional-prior construct (see .displaced_vars); entering ",
                        "both would re-admit the collinearity the construct removes"),
                 hit)
  }
  invisible(TRUE)
}

#' Extract bare variable names from term strings like `s(trend, bs = "ts", k = 10)`.
.term_vars <- function(terms) {
  x <- paste(terms, collapse = " + ")
  x <- gsub("\\bby\\s*=\\s*", " ", x)
  x <- gsub("\\b(bs|k|m|id|sp)\\s*=\\s*[^,)]+", " ", x)   # drop smooth arguments
  v <- unlist(strsplit(x, "[^A-Za-z0-9_.]+"))
  v <- v[nzchar(v)]
  v <- setdiff(v, c("s", "te", "ti", "t2", "ts", "cr", "tp", "offset", "log", "I"))
  unique(v[!grepl("^[0-9.]+$", v)])
}

# --- term construction ------------------------------------------------------

#' Measurement terms for one signal.
#'
#' LEVEL terms answer two independent questions.
#'
#' WHICH PAIR comes from config `level_terms` — `q05`/`q95` where the
#' measurement density supports a quantile, `value_min`/`value_max` where it
#' does not (every sparse signal, plus the three bounded GCS scales).
#'
#' WHICH OF THE PAIR comes from `excursion_side` in pairing.csv, which names the
#' deviation that is clinically meaningful for this signal:
#'
#'   low       -> value_median + <low>
#'   high      -> value_median + <high>
#'   undefined -> value_median + both              (the 7 unpaired controls)
#'
#' PROPENSITY terms do NOT. pi_minus / pi_plus are the posterior means from
#' R/04_features.R, (k_j + alpha_j)/(n_obs + alpha_0), and they are two FREE
#' coordinates of one 3-simplex. Dropping one is a projection, not a variable
#' removal: two patients with the same k_low and the same n_obs but different
#' k_high have identical pi_minus, and only pi_plus separates a purely
#' hypotensive patient from a labile one. MAP is paired on the low side and its
#' high tail still carries information. So both enter by default. pi_mid never
#' does — it is determined by the other two, and entering all three is exactly
#' collinear.
#'
#' The one exception is a coordinate that is not free, named in config's
#' `signal_tails`: where a tail is structurally empty (spo2 and the GCS
#' components cannot exceed their reference ceiling), alpha_j sits at the
#' estimator floor and pi_j reduces to alpha_j/(n_obs + alpha_0) — a monotone
#' function of n_obs and nothing else. It separates no patients, so entering it
#' would smuggle measurement frequency back in under a physiology label.
#'
#' There is deliberately NO `n_obs` term. How many hours a signal was covered is
#' a clinician's ordering decision and is site-specific: monitoring intensity
#' differs between MIMIC and eICU for reasons that have nothing to do with
#' physiology, so an n_obs smooth would transport as an artifact. n enters only
#' where it belongs — inside pi_hat, as the confidence weighting that makes 0/2
#' and 0/24 different evidence. That is the entire reason for shrinking counts
#' rather than modelling k/n.
#' `trend` enters as a plain main-effect smooth and NEVER as `s(trend, by = ...)`.
#' The only interaction layer 1 ever carried was `by = o_flag`, and o_flag is
#' gone. So treatment response as an explicit interaction is not representable
#' in layer 1 as it now stands -- say that in the methods rather than let a
#' reader assume otherwise. `trend` and intervention intensity enter as additive
#' main effects, which is the conservative reading and avoids the near-1
#' concurvity the `by=` version drove in every paired dense model.
measurement_terms <- function(signal, cfg) {
  bs <- cfg$bam$smooth_basis %||% "ts"
  # k is per covariate, not per model: a few covariates are discrete-valued and
  # cannot carry the default basis. See config/smooth_k.
  sm <- function(v) sprintf("s(%s, bs = \"%s\", k = %d)", v, bs,
                            smooth_k_of(signal, v, cfg))

  side  <- excursion_side_of(signal, cfg)
  class <- signal_class_of(signal, cfg)

  # WHICH PAIR is config's `level_terms` (quantile vs extreme); WHICH OF THE
  # PAIR is `excursion_side`. Two independent questions, and conflating them is
  # how a sparse lab ends up carrying a percentile it does not have the
  # measurements to define. See config/level_terms_by_class.
  lv <- level_vars_of(signal, cfg)
  med <- if (isTRUE(cfg$level_terms_include_median)) "value_median" else character(0)
  ext <- if (isTRUE(cfg$level_terms_include_extreme %||% TRUE)) {
    raw <- if (is.na(side)) lv else if (side == "low") lv[1] else lv[2]
    # `magnitude_conditional` REPLACES the raw magnitude term, it does not add
    # one: `{var}_delta` is the same quantity standardised against the deviation
    # count, so the term count is identical either way. See R/04b_conditional.R.
    if (magnitude_conditional_for(signal, cfg)) unname(vapply(raw, delta_name_of, character(1)))
    else raw
  } else character(0)
  level <- c(med, ext)
  tails <- occupiable_tails_of(signal, cfg)
  pis <- c(if ("low" %in% tails) "pi_minus", if ("high" %in% tails) "pi_plus")

  terms <- c(vapply(level, sm, character(1)), vapply(pis, sm, character(1)))

  # trend is only defined for dense and rate signals; sparse ones carry a
  # constant 0 column, which would be a null smooth.
  # Class-gated by config's `trend_classes`, not hard-coded: `trend = 0` means
  # both "flat" and "not estimable", and for sparse signals those two groups
  # have sharply different mortality (SQL audit 6). Declared so the exclusion is
  # a reviewable decision rather than a line of code.
  if (trend_enabled_for(signal, cfg)) terms <- c(terms, sm("trend"))

  unname(terms)
}

#' Intervention terms for one signal, built from its paired interventions.
#'
#' A signal may pair with up to three interventions (the GCS components pair
#' with all three sedatives). Each contributes one block of feature columns,
#' named `{intervention}__{variable}` per spec §2.
#'
#' Every term here is a property of the treatment record ALONE. Nothing in this
#' block refers to the measurement stream, which is exactly what makes `intv` a
#' clean estimate of log[p(I|Y=1)/p(I|Y=0)] and the subtraction meaningful.
#' `present_at_admission` carries the only timing distinction that survives the
#' o_flag removal, and it is site-comparable because it is defined against
#' admission rather than against an excursion the two sites would detect at
#' different charting resolutions.
#'
#' `peak_intensity` is deliberately excluded: spec §7 defers it as nullable and
#' secondary, and it is NULL for the large majority of rows. Including it under
#' na.fail would delete every unexposed patient from the fit.
intervention_terms <- function(signal, cfg) {
  ivs <- interventions_of(signal, cfg)
  if (!length(ivs)) return(character(0))

  bs <- cfg$bam$smooth_basis %||% "ts"
  sm <- function(v) sprintf("s(%s, bs = \"%s\", k = %d)", v, bs,
                            smooth_k_of(signal, v, cfg))

  shape_of <- unlist(cfg$intervention_shape)
  with_agents <- unlist(cfg$interventions_with_agent_counts) %||% character(0)

  # NO `ever_active`. Dropped 2026-08-25. It is exactly `(intensity > 0)` on
  # 100.0000% of stays for every intervention, state- and event-shaped alike, so
  # it was never an independent covariate — only a second encoding of the point
  # mass at zero that the intensity smooth already sits on.
  #
  # MEASURED (tests/term_redundancy.R, section C): what the indicator adds ON TOP
  # OF the smooth is <= 0.0019 deviance explained across all 13 interventions and
  # usually < 0.0006; three are negative. What the smooth adds on top of the
  # indicator is 0.02-0.04. This holds for `invasive_vent` too, the case where
  # the extensive margin looks most primary: 0.0235 for the indicator alone
  # against 0.0603 for the smooth alone. In an ICU cohort 41.6% are ventilated,
  # so the indicator is a coarse split of an already-sick population and the
  # duration is what separates outcomes. For `sedation_propofol` (0.00007) and
  # `diuretic` (0.00002) the indicator alone is worth nothing at all.
  #
  # The earlier argument for keeping it — that a continuous smooth cannot
  # represent the jump at zero — did not survive measurement. `bs = "ts"` has
  # abundant data at exactly zero and gets steep enough there on its own.
  terms <- character(0)
  for (iv in ivs) {
    p <- function(v) paste0(iv, "__", v)
    shape <- shape_of[[iv]]

    # `intensity_conditional` REPLACES the second intensity coordinate with the
    # exposure-standardised one; the exposure term itself is untouched and the
    # term count is unchanged. Only interventions with TWO intensity covariates
    # get a lambda -- a state intervention without agent counts has nothing to
    # condition on. See lambda_spec_of() in R/04b_conditional.R.
    ls <- lambda_spec_of(iv, cfg)
    if (!is.null(ls)) {
      terms <- c(terms, sm(p(ls$exposure)), sm(p("lambda")), p("present_at_admission"))
      next
    }

    terms <- c(terms, if (shape == "state") sm(p("exposure_frac"))
                      else c(sm(p("n_hours")), sm(p("total_amount"))))
    # `n_agents` only. `max_concurrent_agents` was dropped 2026-08-25: MEASURED,
    # the two are equal on 99.03% of training stays and correlate at r = 0.9916,
    # giving observed concurvity 1.00 wherever both entered. The 0.97% that
    # differ are stays given two pressors at different times but never together
    # — a real distinction, but not one 400 stays can identify against a
    # covariate this collinear. Dropping the derived one and keeping the count.
    if (iv %in% with_agents) terms <- c(terms, sm(p("n_agents")))
    terms <- c(terms, p("present_at_admission"))
  }

  unname(terms)
}

#' The interaction scope a ti model name carries.
.ti_scope_of <- function(model) {
  if (!model %in% LAYER1_TI_MODELS) {
    stop(".ti_scope_of: `", model, "` is not an interaction model", call. = FALSE)
  }
  sub("^full_ti_", "", model)
}

#' The smooth-term variables in a term vector.
smooth_term_vars <- function(terms) .term_vars(grep("^s\\(", terms, value = TRUE))

#' Every measurement x intervention cross term for one signal, as `ti()` strings.
#'
#' MOVED HERE FROM `R/14_attribution_eval.R` ON 2026-09-07, VERBATIM in what it
#' emits. It was `attr_ti_terms()` there, which was itself a verbatim move from
#' `tests/coupling_attribution.R` a day earlier, for the same reason both times:
#' the interaction term set is now part of what the pipeline FITS, and a term set
#' defined in an analysis script and re-derived by the pipeline is two
#' definitions of one design. `attr_design_key()`'s `ti_all` / `ti_trend` fields
#' are hashes of these strings, so the strings themselves must not move -- and
#' they have not, which is what lets the stored replicates keep their keys.
#'
#'   scope = "trend"  cross ONLY `trend` with each intervention intensity. The
#'                    treatment-response interaction, and the one three
#'                    independent instruments converged on.
#'   scope = "all"    every measurement x intervention pair.
#'
#' `present_at_admission` is parametric and has no smooth, so it is never
#' crossed. `k` per margin is capped at `k_ti` and floored at 3, and it can only
#' ever be SMALLER than the margin's own `s()` basis -- which is what makes a ti
#' term safe wherever the main-effect smooth already passed `check_model_frame()`.
#'
#' @param k_ti the per-margin basis cap. Read from config by default (it is a
#'   design parameter and belongs in `config/config.yml`, not beside a literal);
#'   the argument stays so `attr_design_key()` can key a stored ladder to the
#'   value it was built under.
interaction_terms <- function(signal, cfg, scope = c("all", "trend"),
                              k_ti = cfg_req(cfg, "k_ti")) {
  scope <- match.arg(scope)
  if (!length(interventions_of(signal, cfg))) return(character(0))
  mv <- smooth_term_vars(measurement_terms(signal, cfg))
  iv <- smooth_term_vars(intervention_terms(signal, cfg))
  if (scope == "trend") mv <- intersect(mv, "trend")
  if (!length(mv) || !length(iv)) return(character(0))
  kt <- function(v) max(3L, min(as.integer(k_ti), smooth_k_of(signal, v, cfg)))
  unlist(lapply(mv, function(m) vapply(iv, function(i)
    sprintf("ti(%s, %s, bs = c(\"ts\", \"ts\"), k = c(%d, %d))",
            m, i, kt(m), kt(i)), character(1))))
}

# --- assembly ---------------------------------------------------------------

#' Build the layer-1 formula for one signal.
#'
#' @param model one of LAYER1_MODELS
#' @return a formula, guaranteed free of dx_/qc_ terms
build_formula <- function(signal, model = LAYER1_MODELS, cfg,
                          response = "mortality") {
  model <- match.arg(model)
  if (!signal %in% cfg$signals) abort_values("unknown signal", signal)

  if (model %in% LAYER1_TI_MODELS && !length(interventions_of(signal, cfg))) {
    stop("build_formula: '", signal, "' is unpaired, so `", model, "` would be ",
         "its `meas` formula with no cross terms. It aliases rather than fits ",
         "there -- use spec_source() / models_of() to enumerate.", call. = FALSE)
  }
  if (model == "intv" && !length(interventions_of(signal, cfg))) {
    stop("build_formula: '", signal, "' is unpaired, so an `intv` formula would ",
         "be `mortality ~ 1`. L_intv is 0 there by construction and is assigned, ",
         "never fitted. Use models_of() to enumerate.", call. = FALSE)
  }

  # A ti model IS `full` plus cross terms, so it takes `full`'s two blocks. The
  # cross terms go LAST, which is not cosmetic: `tests/attr_replicates.R` built
  # these formulas as `update(build_formula(sg, "full", cfg), . ~ . + <ti>)`,
  # which appends, and the 371 replicates on disk were fitted from that. Same
  # term order, same model matrix, same L.
  base <- if (model %in% LAYER1_TI_MODELS) "full" else model
  terms <- c(
    if (base %in% c("meas", "full")) measurement_terms(signal, cfg),
    if (base %in% c("intv", "full")) intervention_terms(signal, cfg),
    if (model %in% LAYER1_TI_MODELS)
      interaction_terms(signal, cfg, .ti_scope_of(model))
  )

  # The guard fires here, on every formula, before it can reach bam().
  assert_no_forbidden(terms, cfg = cfg, signal = signal)

  f <- stats::as.formula(
    paste(response, "~", paste(terms, collapse = " + ")),
    env = globalenv()
  )
  attr(f, "signal") <- signal
  attr(f, "model")  <- model
  f
}

#' Every formula the pipeline will fit, as a table. Cheap, data-free, and worth
#' inspecting before committing to 155 fits.
formula_table <- function(cfg) {
  rows <- list()
  for (sg in cfg$signals) {
    for (md in models_of(sg, cfg)) {
      f <- build_formula(sg, md, cfg)
      rows[[length(rows) + 1L]] <- data.frame(
        signal    = sg,
        model     = md,
        class     = signal_class_of(sg, cfg),
        side      = excursion_side_of(sg, cfg) %||% NA_character_,
        n_paired  = length(interventions_of(sg, cfg)),
        n_terms   = length(attr(stats::terms(f), "term.labels")),
        formula   = paste(deparse(f), collapse = " "),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' One row per (signal, model, term) — the machine-readable audit of every
#' formula the pipeline will fit.
#'
#' Term metadata is re-derived by PARSING the emitted term string, not carried
#' along from measurement_terms() / intervention_terms(). That is deliberate: it
#' makes this an independent check of the builder rather than a restatement of
#' it. A term built wrong shows up here as a term classified wrong, and anything
#' the classifier does not recognise is reported as `UNCLASSIFIED` rather than
#' quietly bucketed.
#'
#' Data-free, like everything else in this file: the whole table can be reviewed
#' before a single row is read.
formula_terms <- function(cfg) {
  rows <- list()
  for (sg in cfg$signals) {
    for (md in models_of(sg, cfg)) {
      f <- build_formula(sg, md, cfg)
      lab <- attr(stats::terms(f), "term.labels")
      for (i in seq_along(lab)) {
        rows[[length(rows) + 1L]] <- .classify_term(lab[i], sg, md, i, cfg)
      }
    }
  }
  do.call(rbind, rows)
}

#' Classify one term string. See formula_terms() for why this parses rather
#' than trusting the builder.
.classify_term <- function(term, signal, model, idx, cfg) {
  # A `ti()` term is classified whole rather than split into its margins: the
  # interaction is the object, and reporting it as two rows would make the
  # review table say the model carries a second `s(trend)`. `variable` is the
  # margin pair, `kind` says it is a tensor interaction, and `source` is fixed
  # rather than derived because a cross term is by construction one measurement
  # variable and one intervention variable -- if it ever is not, that is a
  # builder defect and `UNCLASSIFIED` would hide it behind a plausible label.
  if (grepl("^ti\\(", term)) {
    a <- .call_args(term)
    return(data.frame(signal = signal, model = model, idx = idx, term = term,
                      kind = "tensor_interaction",
                      variable = paste(a$vars, collapse = " x "),
                      source = "interaction",
                      intervention = {
                        iv <- a$vars[grepl("__", a$vars, fixed = TRUE)]
                        if (length(iv)) sub("__.*$", "", iv[1]) else NA_character_
                      },
                      basis = paste(a$bs, collapse = ","),
                      k = paste(a$k, collapse = ","), by = NA_character_,
                      stringsAsFactors = FALSE))
  }
  smooth <- grepl("^s\\(", term)
  bs <- k <- by <- NA_character_
  if (smooth) {
    parts <- trimws(strsplit(sub("^s\\((.*)\\)$", "\\1", term), ",")[[1]])
    v <- parts[1]
    arg <- function(nm) {
      p <- grep(paste0("^", nm, "\\s*="), parts, value = TRUE)
      if (!length(p)) NA_character_
      else trimws(gsub("\"", "", sub(paste0("^", nm, "\\s*=\\s*"), "", p)))
    }
    bs <- arg("bs"); k <- arg("k"); by <- arg("by")
  } else {
    v <- term
  }

  iv <- if (grepl("__", v, fixed = TRUE)) sub("__.*$", "", v) else NA_character_
  src <-
    if (!is.na(iv))                             "intervention"
    # `magnitude` rather than `level`: a delta is a different object from the
    # raw magnitude, not a rescaling of it, and the review table should say so.
    else if (grepl("_delta$", v) &&
             delta_base_of(v) %in% c("q05", "q95", "value_min", "value_max")) "magnitude"
    else if (v %in% c("value_median", "q05", "q95", "value_min", "value_max")) "level"
    else if (v %in% c("pi_minus", "pi_plus"))   "propensity"
    else if (v == "pi_mid")                     "ILLEGAL(pi_mid)"
    else if (v == "trend")                      "trajectory"
    else if (v == "o_flag")                     "ILLEGAL(o_flag)"
    else if (v == "n_obs")                      "ILLEGAL(n_obs)"
    else                                        "UNCLASSIFIED"

  data.frame(signal = signal, model = model, idx = idx, term = term,
             kind = if (smooth) "smooth" else "linear",
             variable = v, source = src,
             intervention = iv %||% NA_character_,
             basis = bs, k = k, by = by,
             stringsAsFactors = FALSE)
}

#' Every smooth in a formula, as (variable, k).
#'
#' Read off the FORMULA rather than recomputed from config, so that whatever
#' check_model_frame() enforces is exactly what bam() will be handed. A k that
#' drifted between the two would otherwise pass the check and fail in mgcv.
smooth_specs <- function(f) {
  lab <- grep("^(s|ti)\\(", attr(stats::terms(f), "term.labels"), value = TRUE)
  if (!length(lab)) {
    return(data.frame(variable = character(0), k = integer(0),
                      stringsAsFactors = FALSE))
  }
  # ONE ROW PER MARGIN, so a `ti()` contributes two. That is what makes the
  # basis guard cover the interaction terms at all: `check_model_frame()` reads
  # this table, and until 2026-09-07 it grepped `^s\\(` alone, so a cross term's
  # bases were the only bases in a formula nothing checked. In practice a ti
  # margin's k is capped at `k_ti` and can only be SMALLER than the same
  # variable's main-effect basis, so this can never refuse a fit the old guard
  # admitted -- which is exactly why it costs nothing to close.
  #
  # PARSED AS AN R CALL rather than by splitting on commas. The string split was
  # correct for `s(v, bs = "ts", k = 10)` and is wrong for
  # `ti(a, b, bs = c("ts", "ts"), k = c(5, 5))`, where the commas INSIDE `c()`
  # are not argument separators. Same output for every `s()` term.
  do.call(rbind, lapply(lab, function(t) {
    a <- .call_args(t)
    kk <- if (length(a$k)) as.integer(a$k) else NA_integer_
    if (length(kk) == 1L && length(a$vars) > 1L) kk <- rep(kk, length(a$vars))
    if (length(kk) != length(a$vars)) kk <- rep(NA_integer_, length(a$vars))
    data.frame(variable = a$vars, k = kk, stringsAsFactors = FALSE)
  }))
}

#' Split a smooth term string into its margin variables and its named arguments.
#'
#' `str2lang()` rather than a regex, because the arguments of a tensor term are
#' vectors: `k = c(5, 5)` has a comma in it and `bs = c("ts", "ts")` has two.
#' Unnamed arguments are the margins, in order; `bs` and `k` are evaluated in an
#' empty environment, so a term carrying anything but a literal fails loudly here
#' rather than silently returning a symbol.
.call_args <- function(term) {
  e <- tryCatch(str2lang(term), error = function(err)
    stop(".call_args: cannot parse smooth term `", term, "`", call. = FALSE))
  a <- as.list(e)[-1]
  nm <- names(a) %||% rep("", length(a))
  vars <- vapply(a[!nzchar(nm)], function(x) paste(deparse(x), collapse = ""),
                 character(1))
  lit <- function(key) {
    if (!key %in% nm) return(NULL)
    eval(a[[key]], envir = baseenv())
  }
  list(vars = unname(vars), bs = lit("bs"), k = lit("k"))
}

#' Columns a model frame must supply for a formula. R/04_features.R builds to
#' this list, so a missing column fails before bam() rather than inside it.
required_columns <- function(f) {
  v <- .term_vars(attr(stats::terms(f), "term.labels"))
  unique(c(all.vars(f)[1], v))
}
