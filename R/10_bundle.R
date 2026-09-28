# R/10_bundle.R --------------------------------------------------------------
# The object that travels, and the single code path that applies it.
#
# Hard rule 8 says run/external.R and run/survival.R fit nothing. Anything an
# external site cannot re-derive therefore has to be carried in a file, and this
# is that file's contents plus the function that consumes them.
#
# WHY THIS MATTERS MORE THAN IT LOOKS. The transportability claim is "the same
# fitted objects, evaluated on a different population". Every quantity an apply
# site re-derives from its own data silently converts that claim into a plumbing
# artifact -- and it converts it in the flattering direction, because a
# re-derived alpha, p_bar or delta is refitted to the site being scored and will
# fit it better. The failure has no symptom: the run completes, the numbers are
# plausible, and the result is about nothing.
#
# So `verify_bundle()` ASSERTS presence and never defaults. An eICU run that
# loaded a bundle without `delta` and fell back to unshrunk magnitudes would
# produce a transport number that is pure plumbing, and it is the one failure
# mode in this project that no downstream diagnostic would catch.
#
# THE DESIGN CONFIG COMES FROM THE BUNDLE, NOT FROM config/config.yml. This is
# the same argument one level up. If eICU read its formulas, whitelists,
# `smooth_k` and level-term parameterisation from the local file, then editing
# that file between the MIMIC fit and the eICU apply would change what is being
# transported without changing anything visible. `bundle_cfg()` returns the
# FROZEN design with only the local `paths` grafted on, and that is what every
# apply function is handed.
#
# NO PATHS, NO CLOCK (hard rule 9). `save_bundle()` and `load_bundle()` receive
# a path; they never build one. Nothing here calls Sys.time().
#
# AGGREGATES ONLY (hard rule 1). Every table this file returns is per (signal,
# model) or per arm. `apply_bundle()` returns row-level scores to its caller --
# a runner, which writes them to a run directory and never prints them.
# ----------------------------------------------------------------------------

BUNDLE_VERSION <- "2.0"

# The arms a bundle can score. Declared here so that MIMIC-test and eICU cannot
# be given different arm lists at their call sites -- if an arm is absent from a
# bundle, that is an error, not a shorter table.
#
# `llr_cond` JOINED THE LIST ON 2026-09-05, and it is not a new fit. CLAUDE.md's
# frozen decision on the three models per signal states that `meas` and `intv`
# partition `full` exactly, so `L_cond = L_full - L_intv` is the
# measurement-deviation LLR CONDITIONAL on class and intervention context,
# obtained by subtraction rather than by a fourth GAM. R/07 has built that
# matrix since the branch point; nothing had ever summed it into a score. The
# design's stated intent is the first term of the factorisation at the end of
# CLAUDE.md's frozen decisions, and that first term is precisely `llr_cond` --
# so the arm the paper argues for was the one arm no table reported.
#
# WHAT IT IS NOT. It is not a fourth model, it costs no fit, and it does not
# change the bundle's design hash. It is a second row-sum of matrices the apply
# path already builds, and the three LLR arms stand in an exact additive
# relation: `llr_sum = llr_cond + rowSums(L_intv)`.
# THE FOUR INTERACTION ARMS JOINED ON 2026-09-07, and unlike `llr_cond` they DO
# cost fits -- 126 of them, taking the budget from 258 to 384. They are here
# rather than in an analysis script for one reason: an arm that exists only
# inside `tests/attr_replicates.R` can be measured for attribution stability at
# MIMIC and can never be TRANSPORTED, because eICU is scored by
# `apply_bundle()` off frozen GAMs and a GAM that is not in the bundle does not
# exist at an apply site. The attribution arm was reporting level-4
# disagreement between `llr_full` and `llr_full_ti_trend` while no table
# anywhere could say what the interaction did to discrimination, calibration or
# transport -- which is the wrong way round.
#
# `llr_sum_ti_*` is `rowSums(L_full_ti_*)` and `llr_cond_ti_*` is
# `rowSums(L_full_ti_* - L_intv)`, so the additive relation holds arm for arm:
# `llr_sum_ti_x = llr_cond_ti_x + rowSums(L_intv)`, with the SAME `L_intv` in
# every one of the three families. That is what makes
# `llr_sum_ti_trend - llr_sum` a clean reading of the interaction alone.
BUNDLE_ARMS <- c("llr_sum", "llr_cond", "llr_meas",
                 "llr_sum_ti_trend", "llr_sum_ti_all",
                 "llr_cond_ti_trend", "llr_cond_ti_all",
                 "xgb_l", "xgb_feat", "xgb_raw")
XGB_DESIGNS <- c("xgb_l", "xgb_feat", "xgb_raw")

# The map from an arm name to the L matrix it sums. Declared rather than derived
# by string surgery, so an arm cannot be added to `BUNDLE_ARMS` without saying
# which matrix it reads -- `apply_bundle()` and `_targets.R` both key off this
# and would otherwise each invent their own parse of the name.
LLR_ARM_MATRIX <- c(
  llr_sum           = "full",
  llr_cond          = "cond",
  llr_meas          = "meas",
  llr_sum_ti_trend  = "full_ti_trend",
  llr_sum_ti_all    = "full_ti_all",
  llr_cond_ti_trend = "cond_ti_trend",
  llr_cond_ti_all   = "cond_ti_all")

# The contrasts every site reports, declared ONCE so MIMIC-test and eICU cannot
# be given different comparison tables. Each is (a, b) and is reported as a - b.
#
# The middle three are the ladder from docs/v2_internal_validation_20260831: the
# gap between the proposed method and a strong learner, decomposed into the
# covariate construction, the per-signal collapse to 19 scalars, and the linear
# aggregation. They sum to `llr_sum - xgb_raw` by construction, which is what
# makes the decomposition checkable rather than merely plausible.
#
# THE TRANSPORT READING. At MIMIC these three say where the gap goes. Run at
# eICU with the SAME frozen models, they say which rung survives the move --
# whether additive evidence transports, and whether the tree's interactions do.
# A rung that is large at MIMIC and near zero at eICU is a rung that was fitting
# MIMIC.
LADDER_CONTRASTS <- list(
  c("llr_sum",  "xgb_raw"),    # does the method hold up against a strong learner
  c("xgb_feat", "xgb_raw"),    # rung 1: the covariate construction
  c("xgb_l",    "xgb_feat"),   # rung 2: the per-signal collapse to 19 scalars
  c("llr_sum",  "xgb_l"),      # rung 3: linear aggregation (the Sigma-inverse case)
  c("xgb_l",    "llr_sum"),    # what a nonlinear aggregator adds over summing
  c("llr_sum",  "llr_meas"),   # what the intervention block adds over physiology
  # The two `llr_cond` contrasts, added 2026-09-05 with the arm itself. They
  # separate the two terms of the factorisation CLAUDE.md states:
  #   llr_sum - llr_cond   what the PURE INTERVENTION-PROPENSITY term adds on
  #                        top of the conditional measurement evidence. It is
  #                        exactly rowSums(L_intv), so a large positive value
  #                        says the score's discrimination leans on who got
  #                        treated rather than on what was measured.
  #   llr_cond - llr_meas  what CONDITIONING on intervention context buys over
  #                        the physiology-only models. Both arms carry
  #                        measurement evidence and neither carries a bare
  #                        propensity term, so this is the cleaner of the two
  #                        readings of "what interventions contribute".
  c("llr_sum",  "llr_cond"),
  c("llr_cond", "llr_meas"),
  # The interaction contrasts, added 2026-09-07 with the arms. Each is read
  # WITHIN a family so that the intervention-propensity block, which is
  # identical across all three, differences out and the number is the
  # interaction alone.
  #
  #   *_ti_trend - *          what treatment response as an explicit
  #                           interaction buys over the additive model. This is
  #                           the contrast the frozen no-interaction decision
  #                           declared unmeasurable in layer 1, and the reason
  #                           it says so is that its stated limitation --
  #                           "treatment response as an explicit interaction is
  #                           not expressible" -- is now a measured cost rather
  #                           than an assumed one.
  #   *_ti_all - *_ti_trend   what the remaining measurement x intervention
  #                           cross terms add on top. Near zero says the trend
  #                           interaction is the whole of it, which is what
  #                           three independent instruments predicted.
  #
  # Reported at BOTH sites. A rung that is positive at MIMIC and negative at
  # eICU is an interaction that was fitting MIMIC, and that is precisely the
  # risk a saturated tensor model carries.
  c("llr_sum_ti_trend",  "llr_sum"),
  c("llr_sum_ti_all",    "llr_sum_ti_trend"),
  c("llr_cond_ti_trend", "llr_cond"),
  c("llr_cond_ti_all",   "llr_cond_ti_trend")
)

#' The contrasts a site can actually report, given the arms it scored.
#'
#' EXTERNAL RUNNER REVIEW E6 (2026-09-09). The runners accept an `arms`
#' subset and then requested every `LADDER_CONTRASTS` pair, so a legitimate
#' reduced arm list failed inside `arm_contrasts()` after the models had been
#' applied and the scores written. Deriving the reportable pairs from the arm
#' list is what `severity_arm()` already does for `SEVERITY_CONTRASTS`; this
#' is the same rule, named, so both runners and their console ladders read
#' one definition. A pair whose arm was not scored is omitted, never zero.
contrasts_available <- function(pairs, arms) {
  arms <- as.character(arms)
  Filter(function(p) all(p %in% arms), pairs)
}

# Config keys that DEFINE the model. A difference in any of these between the
# bundle and a local config means the two are not the same design, and the
# transport comparison is not the comparison it claims to be. `paths` and the
# reporting settings are deliberately absent: those are allowed to differ, and
# in fact must.
#
# DERIVED, NOT CURATED. This list is exactly the set of `cfg$<key>` references
# in R/, minus `paths` and `pairing` which are handled separately. Anything R/
# reads at an apply site must be frozen, or `bundle_cfg()` hands back a config
# with a hole in it -- which is how `intervention_shapes` was found missing:
# the loader validates the shape level set against it, and its absence failed
# at eICU load time rather than at bundle build time.
#
# `xgboost` and `diagnostics` are here even though the apply path fits nothing
# and computes no triage. Both are DECLARED-BEFORE-THE-RUN quantities
# (CLAUDE.md), and freezing them means a hyperparameter or threshold changed
# between the fit and the apply is reported by verify_bundle() rather than
# passing unnoticed.
BUNDLE_DESIGN_KEYS <- c(
  "seed", "split", "n_folds", "fold_stratify_on", "fold_group_by", "guards",
  "signals", "signal_classes", "signal_tails",
  "interventions_extracted", "interventions_modelled",
  "intervention_shapes", "intervention_shape",
  "interventions_with_agent_counts", "intervention_agent_pool",
  "level_terms_include_median", "level_terms_include_extreme",
  "level_terms_by_class", "level_terms_override",
  "magnitude_conditional", "intensity_conditional",
  "magnitude_conditional_override",
  "magnitude_form", "magnitude_form_override", "magnitude_ordinal_scale",
  "trend_classes", "smooth_k", "bam", "xgboost", "diagnostics", "pairing",
  # `k_ti` decides the term set of the 21 interaction specs, so it is design
  # and not reporting: an apply site must not be able to choose a different
  # basis cap than the site that fitted the smooths. Joined 2026-09-07 with the
  # arms themselves, for the same reason `apache` and `sofa` joined in 2026-09.
  "k_ti",
  # The severity comparators are part of the DESIGN, not of the reporting:
  # the coverage floors, the GCS source and the scored organ set decide
  # which stays enter the comparison and how each score is built. Frozen
  # here so `bundle_cfg()` carries them to an apply site, and an edit made
  # between the fit and the apply cannot change them invisibly.
  "apache", "sofa", "severity_symmetric")

# --- stripping ---------------------------------------------------------------

#' Drop the fitted model frame from a gam object (hard rule 6).
#'
#' `predict.gam(newdata = )` does not read `$model`, so a stripped object
#' predicts identically to an unstripped one. It does read `$xlevels`,
#' `$var.summary`, `$smooth` and the coefficients, all of which survive.
#'
#' `deep = TRUE` additionally drops the length-n vectors -- residuals, fitted
#' values, weights, the linear predictor, the response. Those cost a few hundred
#' KB per model per vector across 43 models, and NOTHING in the apply path reads
#' them, because every diagnostic was already extracted by R/08 while the object
#' was alive (hard rule 6). It is off by default because CLAUDE.md mandates
#' exactly `$model <- NULL` and a bundle should carry what the frozen decision
#' says it carries; turn it on only if size becomes a real constraint, and
#' record that it was turned on.
strip_gam <- function(b, deep = FALSE) {
  b$model <- NULL
  if (isTRUE(deep)) {
    for (nm in c("residuals", "fitted.values", "linear.predictors", "weights",
                 "prior.weights", "y", "offset")) b[[nm]] <- NULL
  }
  b
}

# --- the final-only priors container -----------------------------------------

#' Restrict an `llr_priors` container to its final rows.
#'
#' CLAUDE.md: the bundle carries final-stage quantities only -- never fold-level
#' alpha, never fold GAMs. This is that restriction, applied once, at the point
#' the bundle is built, rather than trusted to every call site downstream.
#'
#' The class is preserved, so `priors_for(bundle$priors, sg, "final")` is the
#' identical call the fitting path makes. Asking this object for an `oof` row
#' then fails loudly inside `priors_of()` rather than returning a final row
#' under an out-of-fold name.
priors_final <- function(priors) {
  if (!inherits(priors, "llr_priors")) {
    stop("priors_final: expected the container from layer1_priors()", call. = FALSE)
  }
  keep <- function(d) d[!is.na(d$role) & d$role == "final", , drop = FALSE]
  structure(list(signal       = keep(priors$signal),
                 magnitude    = keep(priors$magnitude),
                 intervention = keep(priors$intervention)),
            class = "llr_priors")
}

# --- construction ------------------------------------------------------------

#' Collect everything a non-fitting site needs into one object.
#'
#' Every argument is a fitted or derived quantity produced by the internal run.
#' Nothing is computed here except the stripping and the restriction to final
#' rows: this function assembles, it does not estimate. That separation is what
#' lets `verify_bundle()` be a pure check.
#'
#' @param cfg        the resolved config. `paths` is REMOVED before storage --
#'                   a bundle carrying MIMIC's data paths and then being loaded
#'                   at eICU is an accident waiting to happen, and the paths are
#'                   not part of the design.
#' @param priors     the full `llr_priors` container; restricted here.
#' @param models     named list of final GAMs, keyed "<signal>/<model>".
#' @param sigma      named list of correlation matrices from L_oof, by model.
#' @param eigen      named list of eigenspectrum tables, by model.
#' @param cutpoints  named list of reporting-bin cut points, by arm.
#' @param xgb        named list of `xgb_fit_full()` outputs, by design.
#' @param layer2     layer-2 weights and whatever derives them. May be NULL --
#'                   the weights do not exist yet (docs/v2_next_steps_20260831
#'                   section 4) and a bundle built before they do is legitimate.
#'                   It is NULL rather than a vector of ones, so an apply site
#'                   that needs them fails instead of silently weighting
#'                   everything equally and calling the result Sigma-inverse
#'                   weighted.
#' @param severity   the frozen severity slot from `make_severity_bundle()`:
#'                   the coverage floors, the MIMIC-fitted recalibration
#'                   intercept and slope per cell, and the restricted-train
#'                   prior they are centred on. NULL is legitimate and means
#'                   the APACHE II / SOFA arms are unavailable at an apply
#'                   site, which `apply_severity()` reports as an error
#'                   rather than as a shorter table.
#' @param train_ref  reference numbers from the training site: cohort event
#'                   rate, n, and the out-of-fold AUROC per arm. Carried so a
#'                   transport table can be assembled without re-reading the
#'                   internal run, and so a bundle can be sanity-checked against
#'                   the run that produced it.
#' @param schema     the post-loader schema signature of the training tables,
#'                   from `schema_signature()`. Frozen so an apply site checks
#'                   its columns, classes and factor levels against WHAT THE
#'                   MODELS WERE FITTED ON rather than against whatever MIMIC
#'                   extraction happens to be on disk at apply time (external
#'                   runner review E7, 2026-09-09). NULL is legitimate for a
#'                   bundle built before this slot existed; the runners then
#'                   fall back to the live comparison and record that they did.
build_bundle <- function(cfg, priors, models,
                         sigma = list(), eigen = list(), cutpoints = list(),
                         xgb = list(), layer2 = NULL, train_ref = list(),
                         severity = NULL, schema = NULL,
                         domains = NULL, site = "mimic", deep_strip = FALSE,
                         source_hashes = list()) {
  design <- cfg[intersect(BUNDLE_DESIGN_KEYS, names(cfg))]
  missing_keys <- setdiff(BUNDLE_DESIGN_KEYS, names(design))
  if (length(missing_keys)) {
    abort_values("build_bundle: config is missing design key(s) the bundle must freeze",
                 missing_keys)
  }
  pf <- priors_final(priors)

  b <- list(
    version     = BUNDLE_VERSION,
    site        = site,
    cfg         = design,
    pairing     = cfg$pairing,
    domains     = domains,
    priors      = pf,
    models      = lapply(models, strip_gam, deep = deep_strip),
    p_bar_train = stats::setNames(pf$signal$p_bar, pf$signal$signal),
    sigma       = sigma,
    eigen       = eigen,
    cutpoints   = cutpoints,
    xgb         = xgb,
    layer2      = layer2,
    train_ref   = train_ref,
    severity    = severity,
    schema      = schema,
    deep_strip  = isTRUE(deep_strip),
    # PASSED IN, NEVER READ FROM DISK HERE. This was `.source_hashes()`, a call
    # into R/11_run.R that listed and hashed `R/` from inside the bundle
    # target. Two things were wrong with it (audit finding F3): a library-layer
    # file was reaching the one file allowed to build paths (hard rule 9), and
    # the `bundle` target acquired an undeclared filesystem input, so its value
    # was not a function of its declared dependencies. `_targets.R` now tracks
    # `R/` as a `format = "file"` target and hands the hashes in here, which
    # makes the same provenance a DECLARED dependency and keeps it inside the
    # object that travels to eICU.
    source_hashes = source_hashes
  )
  class(b) <- "llr_bundle"
  b
}

#' @export
print.llr_bundle <- function(x, ...) {
  cat("<llr_bundle> v", x$version, "  fitted at: ", x$site, "\n", sep = "")
  cat("  models      ", length(x$models), " final GAM(s)\n", sep = "")
  cat("  priors      signal ", nrow(x$priors$signal),
      " | magnitude ", nrow(x$priors$magnitude),
      " | intervention ", nrow(x$priors$intervention), "\n", sep = "")
  cat("  sigma       ", paste(names(x$sigma), collapse = ", "), "\n", sep = "")
  cat("  xgb         ", paste(names(x$xgb), collapse = ", "), "\n", sep = "")
  cat("  cutpoints   ", paste(names(x$cutpoints), collapse = ", "), "\n", sep = "")
  cat("  layer2      ",
      if (is.null(x$layer2)) "absent (weights not yet estimated)"
      else paste(names(x$layer2), collapse = ", "), "\n", sep = "")
  cat("  severity    ",
      if (is.null(x$severity)) "absent"
      else sprintf("%d recalibrated cell(s), p_bar %.4f",
                   nrow(x$severity$recal), x$severity$p_bar), "\n", sep = "")
  invisible(x)
}

# --- verification ------------------------------------------------------------

#' Assert a bundle is complete. Never defaults, never repairs.
#'
#' The four fitted-parameter sets (alpha, p_bar, delta, lambda) are the
#' migration risk and are checked against the grid the FORMULA BUILDER asks for,
#' reconstructed by calling the same data-free accessors layer 1 calls --
#' `magnitude_conditional_for()`, `level_vars_of()`, `lambda_spec_of()`. Checking
#' against a hard-coded list would let the two drift, which is the failure this
#' function exists to prevent.
#'
#' @param cfg  optional local config, compared against the frozen design. A
#'             difference is reported and, under `strict`, fatal.
#' @return invisible data frame of check results, one row per check
verify_bundle <- function(bundle, cfg = NULL, strict = TRUE) {
  if (!inherits(bundle, "llr_bundle")) {
    stop("verify_bundle: not an llr_bundle (got ", class(bundle)[1], ")", call. = FALSE)
  }
  res <- list()
  note <- function(check, ok, detail = "") {
    res[[length(res) + 1L]] <<- data.frame(check = check, ok = isTRUE(ok),
                                           detail = detail, stringsAsFactors = FALSE)
  }
  bc <- bundle$cfg

  # 1. version
  note("version", identical(bundle$version, BUNDLE_VERSION),
       sprintf("bundle %s, code %s", bundle$version %||% "NA", BUNDLE_VERSION))

  # 1b. the severity slot. Delegated to R/10b so the checks live beside the
  # constructor they check, and absent-is-legitimate is decided in one place.
  sev_rows <- verify_severity(bundle, cfg = cfg)
  for (i in seq_len(nrow(sev_rows))) {
    note(sev_rows$check[i], sev_rows$ok[i], sev_rows$detail[i])
  }

  # 2. every fitted final GAM present, and actually a gam
  jobs <- layer1_jobs(bc)
  want <- unique(paste(jobs$signal, jobs$model, sep = "/")[jobs$fit & jobs$role == "final"])
  miss <- setdiff(want, names(bundle$models))
  bad_class <- names(bundle$models)[!vapply(bundle$models, inherits, logical(1), "gam")]
  note("models_present", !length(miss),
       sprintf("%d/%d present%s", length(want) - length(miss), length(want),
               if (length(miss)) paste0("; missing: ", paste(miss, collapse = ", ")) else ""))
  note("models_are_gam", !length(bad_class),
       if (length(bad_class)) paste("not gam:", paste(bad_class, collapse = ", ")) else "all gam")

  # 3. no fold-level anything survived the restriction
  roles <- unique(c(bundle$priors$signal$role, bundle$priors$magnitude$role,
                    bundle$priors$intervention$role))
  note("priors_final_only", !length(roles) || all(roles == "final"),
       sprintf("roles present: %s", paste(roles, collapse = ", ")))

  # 4. alpha and p_bar, one row per signal, all finite and in range
  sp <- bundle$priors$signal
  sig_miss <- setdiff(unlist(bc$signals), sp$signal)
  a_all <- c(sp$alpha_low, sp$alpha_mid, sp$alpha_high)
  a_ok <- nrow(sp) > 0L && all(is.finite(a_all)) && all(a_all > 0)
  p_ok <- nrow(sp) > 0L && all(is.finite(sp$p_bar)) && all(sp$p_bar > 0 & sp$p_bar < 1)
  note("alpha_per_signal", !length(sig_miss) && a_ok,
       sprintf("%d signal(s)%s%s", nrow(sp),
               if (length(sig_miss)) paste0("; missing: ", paste(sig_miss, collapse = ", ")) else "",
               if (!a_ok) "; a non-positive or non-finite alpha" else ""))
  note("p_bar_per_signal", !length(sig_miss) && p_ok,
       sprintf("range [%s, %s]",
               if (nrow(sp)) round(min(sp$p_bar), 4) else NA,
               if (nrow(sp)) round(max(sp$p_bar), 4) else NA))

  # 5. delta, on exactly the grid the formula builder asks for
  want_mag <- do.call(rbind, lapply(unlist(bc$signals), function(sg) {
    if (!magnitude_conditional_for(sg, bc)) return(NULL)
    # The SAME accessor magnitude_priors() calls, unfiltered. Filtering to a
    # hard-coded variable list here would let a new magnitude variable be added
    # to the design and silently escape verification.
    v <- level_vars_of(sg, bc)
    if (!length(v)) return(NULL)
    data.frame(signal = sg, variable = v, stringsAsFactors = FALSE)
  }))
  mg <- bundle$priors$magnitude
  if (is.null(want_mag)) {
    note("delta_grid", nrow(mg) == 0L,
         "magnitude_conditional is off for every signal; no delta rows expected")
  } else {
    have <- paste(mg$signal, mg$variable)
    wants <- paste(want_mag$signal, want_mag$variable)
    d_miss <- setdiff(wants, have); d_extra <- setdiff(have, wants)

    # FORM-AWARE. `a0` is NA by construction under the ordinal form -- the
    # intercept is absorbed into the cut points -- so a blanket finiteness test
    # over every numeric column would fail a correct bundle. Each form is
    # checked against the parameters it actually has.
    fm <- mg$form %||% rep("linear", nrow(mg))
    lin <- which(fm == "linear"); ord <- which(fm == "ordinal")
    finite <- nrow(mg) > 0L &&
      all(is.finite(c(mg$a1, mg$a2, mg$s_e, mg$s_u))) &&
      (!length(lin) || all(is.finite(mg$a0[lin]))) &&
      (!length(ord) || all(is.finite(mg$a3[ord])))

    # An ordinal row is useless without its cut points, and a cut-point vector
    # of the wrong length would silently mis-bin every category at the apply
    # site. Both are checked against the DECLARED scale rather than against
    # whatever the row happens to carry.
    ord_bad <- character(0)
    for (i in ord) {
      th <- .unpack_num(mg$theta[i]); lv <- .unpack_num(mg$levels[i])
      dec <- bc$magnitude_ordinal_scale[[mg$signal[i]]]
      want_lv <- if (is.null(dec)) integer(0) else .ordinal_levels(dec)
      ok <- length(th) == length(lv) - 1L && length(th) > 0L &&
        !anyNA(th) && all(diff(th) > 0) && identical(as.numeric(lv), as.numeric(want_lv))
      if (!isTRUE(ok)) ord_bad <- c(ord_bad, paste(mg$signal[i], mg$variable[i]))
    }

    note("delta_grid", !length(d_miss) && !length(d_extra) && finite && !length(ord_bad),
         sprintf("%d/%d pair(s), %d ordinal%s%s%s%s", nrow(mg), nrow(want_mag), length(ord),
                 if (length(d_miss))  paste0("; MISSING: ", paste(d_miss, collapse = ", ")) else "",
                 if (length(d_extra)) paste0("; unexpected: ", paste(d_extra, collapse = ", ")) else "",
                 if (!finite) "; a non-finite parameter" else "",
                 if (length(ord_bad)) paste0("; BAD CUT POINTS: ",
                                             paste(ord_bad, collapse = ", ")) else ""))
  }

  # 6. lambda, likewise
  want_lam <- Filter(function(iv) !is.null(lambda_spec_of(iv, bc)),
                     unique(unlist(lapply(unlist(bc$signals), interventions_of, cfg = bc))))
  ip <- bundle$priors$intervention
  l_miss <- setdiff(want_lam, ip$intervention); l_extra <- setdiff(ip$intervention, want_lam)
  note("lambda_grid", !length(l_miss) && !length(l_extra),
       sprintf("%d/%d intervention(s)%s%s", nrow(ip), length(want_lam),
               if (length(l_miss))  paste0("; MISSING: ", paste(l_miss, collapse = ", ")) else "",
               if (length(l_extra)) paste0("; unexpected: ", paste(l_extra, collapse = ", ")) else ""))

  # 7. sigma: square, symmetric, no NA
  s_ok <- length(bundle$sigma) > 0L && all(vapply(bundle$sigma, function(S)
    is.matrix(S) && nrow(S) == ncol(S) && !anyNA(S) &&
      isTRUE(all.equal(unname(unclass(S)), unname(t(unclass(S))))), logical(1)))
  note("sigma", s_ok,
       sprintf("%d matrix/matrices: %s", length(bundle$sigma),
               paste(vapply(bundle$sigma, function(S) paste(dim(S), collapse = "x"),
                            character(1)), collapse = ", ")))

  # 8. cut points: strictly increasing, long enough to define bins
  cp_bad <- names(bundle$cutpoints)[!vapply(bundle$cutpoints, function(v)
    is.numeric(v) && length(v) > 2L && !anyNA(v) && all(diff(v) > 0), logical(1))]
  note("cutpoints", length(bundle$cutpoints) > 0L && !length(cp_bad),
       sprintf("%d arm(s)%s", length(bundle$cutpoints),
               if (length(cp_bad)) paste0("; not strictly increasing: ",
                                          paste(cp_bad, collapse = ", ")) else ""))

  # 8b. the cut points and the training reference cover EXACTLY the declared
  # arm set. Added 2026-09-08 (plumbing review F2): check 8 verified whatever
  # arms were supplied and said nothing about completeness, so a bundle built
  # from nine out-of-fold arms passed and the tenth was discovered at an apply
  # site, or -- before `score_arms()` refused a missing frozen binning -- was
  # quietly self-binned. Both sets are derived from `oof_scores`, so a
  # difference between them is itself a defect.
  arm_set <- function(have, what) {
    miss <- setdiff(BUNDLE_ARMS, have); extra <- setdiff(have, BUNDLE_ARMS)
    note(paste0(what, "_cover_arms"), !length(miss) && !length(extra),
         sprintf("%d/%d arm(s)%s%s", length(intersect(have, BUNDLE_ARMS)), length(BUNDLE_ARMS),
                 if (length(miss))  paste0("; MISSING: ", paste(miss, collapse = ", ")) else "",
                 if (length(extra)) paste0("; unexpected: ", paste(extra, collapse = ", ")) else ""))
  }
  arm_set(names(bundle$cutpoints), "cutpoints")
  arm_set(as.character(bundle$train_ref$arms$label), "train_ref")

  # 9. the boosters, with their feature sets and their training prior
  x_bad <- character(0)
  for (nm in XGB_DESIGNS) {
    m <- bundle$xgb[[nm]]
    if (is.null(m) || is.null(m$booster) || !length(m$feature_names) ||
        !is.finite(m$p_bar %||% NA_real_)) x_bad <- c(x_bad, nm)
  }
  note("xgb_models", !length(x_bad),
       sprintf("%s%s",
               paste(vapply(intersect(XGB_DESIGNS, names(bundle$xgb)), function(nm)
                 sprintf("%s:%d cols", nm, length(bundle$xgb[[nm]]$feature_names)),
                 character(1)), collapse = ", "),
               if (length(x_bad)) paste0("; INCOMPLETE: ", paste(x_bad, collapse = ", ")) else ""))

  # 10. the frozen design against a local config, when one is offered
  if (!is.null(cfg)) {
    diffs <- unlist(Filter(Negate(is.null), lapply(BUNDLE_DESIGN_KEYS, function(k) {
      if (isTRUE(all.equal(bc[[k]], cfg[[k]]))) NULL else k
    })))
    note("design_matches_local_config", !length(diffs),
         if (length(diffs)) paste("differs in:", paste(diffs, collapse = ", "))
         else "identical on all design keys")
  }

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  if (any(!out$ok)) {
    msg <- paste0("verify_bundle: ", sum(!out$ok), " check(s) failed.\n",
                  paste(sprintf("  [FAIL] %-28s %s", out$check[!out$ok], out$detail[!out$ok]),
                        collapse = "\n"),
                  "\nThe bundle is the sole input to every non-fitting run. A missing ",
                  "fitted quantity would be silently defaulted downstream and the ",
                  "transport result would be a plumbing artifact (hard rule 8).")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(out)
}

# --- round trip --------------------------------------------------------------

#' Read a bundle from disk and verify it before returning.
#'
#' Verification is not optional and not a separate step the caller may forget: a
#' bundle that fails its checks must never reach an apply function.
#' `save_bundle()` lives in R/11_run.R, because writing needs a run directory and
#' this file builds no paths (hard rule 9).
load_bundle <- function(path, cfg = NULL, strict = TRUE, verbose = TRUE) {
  if (!file.exists(path)) stop("no bundle at ", path, call. = FALSE)
  b <- qs2::qs_read(path)
  v <- verify_bundle(b, cfg = cfg, strict = strict)
  if (verbose) {
    message(sprintf("bundle loaded: v%s fitted at %s; %d checks, %d passed",
                    b$version, b$site, nrow(v), sum(v$ok)))
  }
  b
}

#' The frozen design config, with local paths grafted on.
#'
#' THE ONLY CONFIG AN APPLY SITE MAY USE. Reading `config/config.yml` at eICU
#' would let an edit between the fit and the apply change the design silently;
#' reading the bundle alone would leave the runner with no data paths. This
#' returns exactly the first plus exactly the second.
#'
#' @param paths the site's `paths` block, from its own config
bundle_cfg <- function(bundle, paths) {
  cfg <- bundle$cfg
  cfg$paths   <- paths
  cfg$pairing <- bundle$pairing
  # THE STAMP. `apply_bundle()` requires it and checks it against the bundle it
  # was handed. Comparing a few design VALUES instead would not work: a local
  # config that happens to agree with the frozen design today would pass, and
  # the whole point is to catch the day it stops agreeing. A hash of the frozen
  # design says "this config was derived from THIS bundle" and nothing else can
  # forge it by coincidence.
  cfg$.bundle_design <- .hash(bundle$cfg)
  cfg
}

#' One row per stored quantity: what a bundle contains, in counts.
#'
#' AGGREGATES ONLY (hard rule 1). Goes into the run manifest so a result can be
#' traced to the exact contents of the object that produced it.
bundle_summary <- function(bundle) {
  row <- function(item, n, detail) data.frame(item = item, n = n, detail = detail,
                                              stringsAsFactors = FALSE)
  do.call(rbind, list(
    row("version",        NA_integer_, bundle$version),
    row("site",           NA_integer_, bundle$site),
    row("models",         length(bundle$models), ""),
    row("signals",        length(unlist(bundle$cfg$signals)), ""),
    row("alpha_rows",     nrow(bundle$priors$signal), ""),
    row("delta_rows",     nrow(bundle$priors$magnitude), ""),
    row("lambda_rows",    nrow(bundle$priors$intervention), ""),
    row("sigma",          length(bundle$sigma), paste(names(bundle$sigma), collapse = ", ")),
    row("cutpoints",      length(bundle$cutpoints), paste(names(bundle$cutpoints), collapse = ", ")),
    row("xgb",            length(bundle$xgb), paste(names(bundle$xgb), collapse = ", ")),
    row("layer2_weights", length(bundle$layer2$weights),
        if (is.null(bundle$layer2)) "absent" else "present"),
    row("schema_signature", length(bundle$schema),
        if (is.null(bundle$schema))
          "absent (built before the frozen signature; apply sites compare against the live MIMIC extraction)"
        else paste(names(bundle$schema), collapse = ", ")),
    row("severity_cells",
        if (is.null(bundle$severity)) 0L else nrow(bundle$severity$recal),
        if (is.null(bundle$severity)) "absent"
        else sprintf("min_vars %s / min_organs %s / gcs %s / resp %s",
                     bundle$severity$settings$min_vars,
                     bundle$severity$settings$min_organs,
                     bundle$severity$settings$gcs_source,
                     bundle$severity$settings$resp_support))))
}

# --- the apply path ----------------------------------------------------------

#' Score a site with a frozen bundle. FITS NOTHING.
#'
#' THIS IS THE FUNCTION THAT MAKES THE TRANSPORT CLAIM CHECKABLE. MIMIC-test and
#' eICU both go through it, so "the same code path" is a property of the code
#' rather than a claim about two runner scripts that happen to look alike. The
#' runners differ only in which parquet they load and which extra tables they
#' write; every number in the comparison is produced here.
#'
#' Five arms, all on the log-odds scale and all centred on a TRAINING prior, so
#' every one of them can go straight into `score_report()`:
#'
#'   llr_sum   rowSums of L_full  -- the proposed method
#'   llr_cond  rowSums of L_cond  -- the CONDITIONAL measurement-deviation LLR,
#'             L_full - L_intv, which is the first term of the design's stated
#'             factorisation. Derived by subtraction, never fitted.
#'   llr_meas  rowSums of L_meas  -- physiology only; the APACHE II counterpart
#'   xgb_l     the frozen booster on the L matrix
#'   xgb_feat  the frozen booster on the layer-1 covariates
#'   xgb_raw   the frozen booster on the raw feature matrix
#'
#' `xgb_feat` uses `role = "final"` priors, which is the whole reason the
#' fold-dependence of `delta` and `lambda` was made explicit in
#' `xgb_design_feat()`: at an apply site there are no folds, and the frozen final
#' parameters are exactly what must be used.
#'
#' The three boosters are realigned onto their own stored `feature_names`, so a
#' site whose column set differs is scored on the MODEL's columns rather than on
#' its own. `arm_table` records how many columns were absent, which is the number
#' that says whether a weak transport result is physiology or plumbing.
#'
#' @param bundle   from load_bundle()
#' @param tabs     the site's tables, from load_tables()
#' @param cfg      MUST be bundle_cfg(bundle, <site paths>). Passed rather than
#'                 derived so the caller's intent is visible, and checked.
#' @param stay_ids stays to score
#' @param arms     which arms to compute. Absent from the bundle is an error.
#' @return list(scores, l_long, l_mats, coverage, arm_table). `scores` is
#'         row-level and must go to a run directory, never to the console.
apply_bundle <- function(bundle, tabs, cfg, stay_ids,
                         arms = BUNDLE_ARMS, verbose = TRUE) {
  bad <- setdiff(arms, BUNDLE_ARMS)
  if (length(bad)) abort_values("apply_bundle: unknown arm(s)", bad)
  stamp <- cfg$.bundle_design
  if (is.null(stamp) || !identical(stamp, .hash(bundle$cfg))) {
    stop("apply_bundle: `cfg` did not come from bundle_cfg(bundle, paths)",
         if (is.null(stamp)) " (no stamp)" else " (stamp is from a different bundle)",
         ". Reading config/config.yml at an apply site lets an edit made ",
         "between the fit and the apply change the design without changing ",
         "anything visible (hard rule 8).", call. = FALSE)
  }
  ids <- as.character(stay_ids)
  pri <- bundle$priors
  out <- list(scores = list(), arm_table = list())

  # --- layer 1, evaluated ---------------------------------------------------
  if (verbose) message("apply_bundle: layer 1 over ", length(stay_ids), " stays")
  a1 <- apply_layer1(bundle$models, tabs, cfg, pri, stay_ids, verbose = verbose)
  out$l_long   <- a1$l
  out$coverage <- a1$coverage
  out$l_mats   <- l_matrices(a1$l, tabs, cfg, stay_ids, fill = "zero")

  # ONE LOOP OVER `LLR_ARM_MATRIX`, not one `if` per arm. With three arms the
  # repeated block was readable; with seven it is seven places to forget a
  # `setNames()` or to sum the wrong matrix, and the two that would be easiest
  # to transpose -- `cond_ti_trend` and `cond_ti_all` -- differ by four
  # characters. `l_matrices()` always returns every one of these: `cond` is
  # `full - intv`, and on the 7 unpaired signals `intv` is assigned 0 rather
  # than fitted, so those columns equal `full`.
  for (nm in intersect(names(LLR_ARM_MATRIX), arms)) {
    mt <- LLR_ARM_MATRIX[[nm]]
    if (is.null(out$l_mats[[mt]])) {
      stop("apply_bundle: arm `", nm, "` was requested but no `", mt, "` L ",
           "matrix could be built. The bundle must carry every fitted spec of ",
           "that model -- an absent arm must fail rather than silently shorten ",
           "the comparison table.", call. = FALSE)
    }
    out$scores[[nm]] <- stats::setNames(rowSums(out$l_mats[[mt]]), ids)
  }

  # --- the frozen boosters --------------------------------------------------
  for (nm in intersect(XGB_DESIGNS, arms)) {
    m <- bundle$xgb[[nm]]
    if (is.null(m)) {
      stop("apply_bundle: the bundle carries no `", nm, "` booster. An absent ",
           "arm must fail rather than silently shorten the comparison table.",
           call. = FALSE)
    }
    if (verbose) message("apply_bundle: ", nm)
    X <- switch(nm,
      xgb_l    = xgb_design_L(a1$l, tabs, cfg, stay_ids, model = "full", fill = "zero"),
      xgb_feat = xgb_design_feat(tabs, cfg, pri, stay_ids, role = "final",
                                 fold = NA_integer_, feature_names = m$feature_names),
      xgb_raw  = xgb_design_raw(tabs, cfg, stay_ids, feature_names = m$feature_names))
    r <- xgb_apply(m, X)
    out$scores[[nm]] <- stats::setNames(as.numeric(r$score), ids)
    out$arm_table[[nm]] <- data.frame(
      arm = nm,
      n_cols_model = length(m$feature_names),
      n_cols_all_na_here = sum(vapply(seq_len(ncol(X)), function(j) all(is.na(X[, j])),
                                      logical(1))),
      frac_na = round(mean(is.na(X)), 5),
      p_bar_train = round(m$p_bar, 6),
      stringsAsFactors = FALSE)
  }

  # Row alignment is asserted, not assumed: five arms scored on five different
  # id orders would produce a comparison table where every paired test is wrong
  # and nothing looks wrong.
  for (nm in names(out$scores)) {
    if (length(out$scores[[nm]]) != length(ids) ||
        !identical(names(out$scores[[nm]]), ids)) {
      stop("apply_bundle: arm `", nm, "` is not aligned to `stay_ids`.", call. = FALSE)
    }
  }
  out$arm_table <- if (length(out$arm_table)) do.call(rbind, out$arm_table) else NULL
  out
}
