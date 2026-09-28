# run/external.R -------------------------------------------------------------
# eICU. Loads a bundle, scores 95,507 stays across 166 hospitals, FITS NOTHING.
#
# THE TWO QUESTIONS THIS RUN EXISTS TO ANSWER.
#
#   1. Does ADDITIVE EVIDENCE transport? `llr_sum` is a sum of 19 marginal
#      log-likelihood ratios, every one of them a frozen MIMIC smooth evaluated
#      on an eICU patient. If the sum discriminates at eICU as it does on
#      held-out MIMIC, the additive evidence geometry is a property of critical
#      illness rather than of one hospital's charting.
#
#   2. Do XGBOOST'S INTERACTIONS transport? The boosters were fitted at MIMIC
#      and carry whatever interactions MIMIC supported. The ladder --
#      `xgb_raw -> xgb_feat -> xgb_l -> llr_sum` -- is computed at both sites
#      from the SAME frozen models. A rung that is large at MIMIC and small at
#      eICU is CONSISTENT WITH that rung having fitted MIMIC; it is a hypothesis
#      the two tables raise, not a conclusion they prove, because coverage,
#      case mix and charting differ between the sites as well.
#
# HARD RULE 8, ENFORCED RATHER THAN INTENDED. Nothing here fits. alpha, p_bar,
# the delta and lambda parameters, the 43 GAMs, the three boosters and the
# reporting cut points are all read from the bundle. `assign_folds()` is NOT
# called: there is no split at eICU, every stay is scored, and there is nothing
# to hold out from a model that was never fitted here.
#
# THE DESIGN COMES FROM THE BUNDLE. `bundle_cfg()` returns the frozen MIMIC
# design with eICU's paths grafted on. config/config.yml is read ONLY for
# reporting settings and to be compared against the frozen design by
# verify_bundle(), which reports any divergence as a failed check.
#
# A FAILED DESIGN CHECK HERE IS A RESULT, NOT A CRASH -- UNLESS STRICT IS SET.
# `check_agent_pool()`, `check_signal_tails()` and the rest hold static
# declarations to the data in both directions. A declaration that holds at
# MIMIC and fails at eICU is exactly what external validation is for, so under
# `design_checks_strict: false` they warn, are recorded, and the run continues.
# Under `true` the table is written and the run stops before any model is
# applied (external runner review E2). An exception a check could not even
# evaluate is fatal in both modes.
#
# THE MANIFEST DESCRIBES THIS RUN, NOT THE TRAINING RUN (review E3). Its
# `config` block is the frozen bundle design that was applied, the full
# resolved `config/external.yml`, the CLI overrides, the reporting settings,
# the bundle's own identity (file hash, design hash, fit-time source hashes,
# producing run) and a size-plus-MD5 fingerprint of every input file. A stage
# that raises leaves status `failed` with the stage and the reason (review E1).
#
# AGGREGATES ONLY (hard rule 1). eICU is under the same PhysioNet DUA. Scores,
# L's and the per-hospital table go to the run directory as .rds; the console
# gets counts, AUROCs and a distribution. HOSPITAL IDENTIFIERS ARE NOT PRINTED.
#
#   Rscript run/external.R
#   Rscript run/external.R out/runs/internal_.../bundle.qs2
#   Rscript run/external.R - --no-hospital     # pooled only, faster
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(xgboost); library(qs2)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
ec   <- yaml::read_yaml("config/external.yml")
no_hosp <- "--no-hospital" %in% args

bundle_from_cli <- length(args) && args[1] != "-" && nzchar(args[1])
bundle_path <- if (bundle_from_cli) args[1] else ec$bundle
if (is.null(bundle_path) || !nzchar(bundle_path)) {
  stop("no bundle. Set `bundle:` in config/external.yml, or pass a path.\n",
       "It is an EXPLICIT path on purpose: a published external result must ",
       "not change because someone re-ran the internal pipeline.", call. = FALSE)
}

# --- validate the runner configuration BEFORE any input is loaded -----------
# External runner review E6. Two incompatible settings used to fail late: an
# `arms` subset passed `apply_bundle()` and then failed inside
# `arm_contrasts()` after every model had been applied and every score
# written, and `frozen_bins: false` with `self_bins: false` left no summary for
# the transport table to read. Both are refused here, in milliseconds, with
# the rest of the runner's switches checked alongside. The function returns
# the contrasts the selected arms can support, derived by
# `contrasts_available()` -- the same rule `severity_arm()` applies to its own
# contrast list -- so a reduced arm set is reported on what it scored.
#
# A reduced arm set is NOT a computation-saving control: `apply_bundle()`
# applies every primary GAM regardless, because the L matrices the LLR arms
# sum are built once for all of them. The option narrows the report, not the
# work; making application dependency-aware would be a separate change.
validate_external_config <- function(ec, no_hospital = FALSE) {
  bad <- function(...) stop("config/external.yml: ", ..., call. = FALSE)
  is_flag <- function(x) is.null(x) || (is.logical(x) && length(x) == 1L && !is.na(x))

  arms <- as.character(unlist(ec$arms %||% BUNDLE_ARMS))
  if (!length(arms)) bad("`arms` is empty")
  if (anyDuplicated(arms)) abort_values("config/external.yml: duplicated arm(s)", arms[duplicated(arms)])
  unk <- setdiff(arms, BUNDLE_ARMS)
  if (length(unk)) abort_values("config/external.yml: `arms` names arm(s) no bundle can score", unk)

  frozen <- isTRUE(ec$frozen_bins %||% TRUE)
  self   <- isTRUE(ec$self_bins   %||% TRUE)
  if (!frozen && !self) {
    bad("`frozen_bins` and `self_bins` are both false, so no score summary ",
        "would exist for the transport table. Enable at least one binning.")
  }

  for (k in c("schema_equality", "design_checks_strict", "prior_fit_diagnostics")) {
    if (!is_flag(ec$checks[[k]])) bad("`checks.", k, "` must be true or false")
  }
  if (!is_flag(ec$severity$enabled)) bad("`severity.enabled` must be true or false")
  if (!is_flag(ec$hospital$enabled)) bad("`hospital.enabled` must be true or false")

  hp <- ec$hospital %||% list()
  if (isTRUE(hp$enabled) && !no_hospital) {
    if (is.null(ec$hospital_table) || !nzchar(ec$hospital_table)) {
      bad("`hospital.enabled` is true but `hospital_table` names no file")
    }
    if (!is.character(hp$group_col %||% "hospitalid")) bad("`hospital.group_col` must be a string")
    for (k in c("min_stays", "min_events")) {
      v <- hp[[k]]
      if (is.null(v) || !is.numeric(v) || length(v) != 1L || v < 1 || v != round(v)) {
        bad("`hospital.", k, "` must be a positive integer")
      }
    }
    if (is.null(hp$unmatched_policy) ||
        !hp$unmatched_policy %in% c("exclude", "error")) {
      bad("`hospital.unmatched_policy` must be declared as `exclude` or `error` ",
          "(review E4): it decides what population the heterogeneity table describes")
    }
  }

  need <- c("cohort", "signal_features", "intervention_features")
  miss <- setdiff(need, names(ec$paths %||% list()))
  if (length(miss)) abort_values("config/external.yml: `paths` is missing", miss)

  ladder <- contrasts_available(LADDER_CONTRASTS, arms)
  dropped <- Filter(function(p) !all(p %in% arms), LADDER_CONTRASTS)
  list(arms = arms, frozen = frozen, self = self, ladder = ladder,
       dropped_contrasts = vapply(dropped, paste, character(1), collapse = " - "))
}
rc <- validate_external_config(ec, no_hospital = no_hosp)
arms <- rc$arms

cfg_local <- load_config("config/config.yml")
bundle    <- load_bundle(bundle_path, cfg = cfg_local, strict = TRUE)
cfg       <- bundle_cfg(bundle, ec$paths)

nb  <- cfg_local$metrics$n_bins %||% 20L
nbt <- cfg_local$metrics$n_boot %||% 200L

# --- the run and its provenance ---------------------------------------------
# External runner review E3. `new_run()` used to receive `cfg_local`, the MIMIC
# config, so the manifest's `config` and `config_hash` described inputs this
# run never read and omitted the ones it did. What is snapshotted now:
#
#   design    the frozen bundle design actually applied (`bundle$cfg`)
#   external  the complete resolved config/external.yml
#   cli       the command line, including the bundle override and --no-hospital
#   reporting bin and bootstrap counts, seed, and the local config's hash
#   bundle    file MD5 and size, design hash, version, fit-time source hashes,
#             and the producing internal run's id, status and config hash
#   inputs    size and MD5 of every input file, the hospital table included
#
# None of it is row-level: hashes, counts, paths and declarations only.
bundle_manifest <- tryCatch(read_manifest(dirname(bundle_path)), error = function(e) NULL)
input_paths <- c(ec$paths, if (isTRUE(ec$hospital$enabled) && !no_hosp)
                              list(hospital_table = ec$hospital_table))
run_cfg <- list(
  runner   = "run/external.R",
  site     = ec$site %||% "eicu",
  design   = bundle$cfg,
  pairing_hash = .hash(bundle$pairing),
  external = ec,
  cli      = list(args = as.list(args), bundle_from_cli = isTRUE(bundle_from_cli),
                  no_hospital = no_hosp),
  reporting = list(n_bins = nb, n_boot = nbt, seed = cfg$seed,
                   local_config = "config/config.yml",
                   local_config_hash = .hash(cfg_local)),
  bundle   = list(path = bundle_path,
                  bytes = as.numeric(file.size(bundle_path)),
                  md5 = unname(as.character(tools::md5sum(bundle_path))),
                  version = bundle$version, fitted_at = bundle$site,
                  design_hash = .hash(bundle$cfg),
                  has_schema_signature = !is.null(bundle$schema),
                  has_severity_slot = !is.null(bundle$severity),
                  source_hashes = bundle$source_hashes,
                  producing_run = if (is.null(bundle_manifest)) "no manifest beside the bundle"
                                  else list(run_id = bundle_manifest$run_id,
                                            status = bundle_manifest$status,
                                            config_hash = bundle_manifest$config_hash,
                                            started_at = bundle_manifest$started_at)),
  inputs   = .input_fingerprint(input_paths)
)

run <- new_run("external", run_cfg, note = sprintf(
  "eICU external validation. bundle %s. Fits nothing (hard rule 8).",
  basename(dirname(bundle_path))))
log_msg(run, "bundle: ", bundle_path, "  md5 ", run_cfg$bundle$md5)
if (length(rc$dropped_contrasts)) {
  log_msg(run, "arms subset: ", length(rc$dropped_contrasts), " ladder contrast(s) not ",
          "reportable: ", paste(rc$dropped_contrasts, collapse = "; "))
}
save_table(run, verify_bundle(bundle, cfg = cfg_local, strict = FALSE),
           "bundle_checks", subdir = "diagnostics")
save_table(run, bundle_summary(bundle), "bundle_contents", subdir = "diagnostics")

# The manifest content known so far, for a `failed` manifest as much as for
# the final one. Extended as stages complete.
arm_status <- list(severity = "pending", hospital = "pending",
                   schema_reference = "pending", site_checks = "pending")
manifest_core <- function() list(
  bundle_path = bundle_path,
  compare_to  = ec$compare_to %||% "",
  arms        = as.list(arms),
  contrasts_reported = vapply(rc$ladder, paste, character(1), collapse = " - "),
  contrasts_dropped  = as.list(rc$dropped_contrasts),
  arm_status  = arm_status)

# --- load -------------------------------------------------------------------
tabs <- run_stage(run, "load_tables",
                  load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE),
                  extra = manifest_core())
log_msg(run, sprintf("eICU: %d stays, %d signal rows, %d intervention rows",
                     nrow(tabs$cohort), nrow(tabs$signal_features),
                     nrow(tabs$intervention_features)))

val <- run_stage(run, "validate_tables", validate_tables(tabs, cfg, strict = TRUE),
                 extra = manifest_core())
save_table(run, val, "validation", subdir = "diagnostics")

# Column-by-column against the TRAINING schema. This is what turns "the schemas
# are supposed to be identical" (hard rule 5) into something that is checked at
# the moment it would otherwise silently stop being true.
#
# External runner review E7. The reference is the post-loader schema signature
# FROZEN IN THE BUNDLE when it carries one: names, classes and factor levels
# of the tables the models were fitted on, compared without loading a single
# MIMIC row. A bundle built before that slot existed falls back to loading the
# live MIMIC extraction -- which is a comparison against whatever is on disk
# today, not against the fit -- and the manifest says which reference was
# used.
arm_status$schema_reference <- "disabled"
if (isTRUE(ec$checks$schema_equality)) {
  if (!is.null(bundle$schema)) {
    se <- compare_schema_signatures(bundle$schema, schema_signature(tabs), strict = FALSE,
                                    label_a = "bundle (mimic, fit time)", label_b = "eicu")
    arm_status$schema_reference <- "bundle_signature"
  } else {
    log_msg(run, "schema check: the bundle carries no frozen schema signature; ",
            "comparing against the LIVE MIMIC extraction (pre-E7 fallback)")
    tabs_mimic <- load_tables(cfg_local$paths$mimiciv, cfg_local, site = "mimic",
                              verbose = FALSE)
    se <- validate_schema_equality(tabs_mimic, tabs, strict = FALSE)
    rm(tabs_mimic); invisible(gc())
    arm_status$schema_reference <- "live_mimic_extraction"
  }
  save_table(run, as.data.frame(se), "schema_equality", subdir = "diagnostics")
}

# --- the two known site risks ----------------------------------------------
# Both already mitigated in code and unverified in eICU data until now
# (docs/v2_next_steps_20260831 section 3). Non-strict: a failure is a finding.
strict_checks <- isTRUE(ec$checks$design_checks_strict)

# Capture the outcome without letting a non-strict failure vanish. Under
# `strict = FALSE` these checks WARN rather than stop, and a bare tryCatch would
# record `ok = TRUE` for a warning -- turning the one class of eICU finding this
# section exists to surface into a silent pass.
#
# THREE OUTCOMES, NOT TWO (review E2). `warning` is a declaration discrepancy
# the validator chose to report rather than stop on; `error` is either the
# strict-mode stop or an exception the validator could not get past. The two
# used to collapse into `ok = FALSE`, and the run applied the models either
# way, so `design_checks_strict: true` enforced nothing.
.capture_check <- function(fn) {
  msgs <- character(0); kind <- "pass"
  val <- withCallingHandlers(
    tryCatch(fn(), error = function(e) {
      msgs <<- c(msgs, conditionMessage(e)); kind <<- "error"; NULL
    }),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      if (kind != "error") kind <<- "warning"
      invokeRestart("muffleWarning")
    })
  list(ok = kind == "pass", outcome = kind, value = val,
       msg = paste(msgs, collapse = " | "))
}

site_checks <- list(
  # The four resampling declarations, held to eICU's cohort. eICU has 166
  # hospitals and repeat stays, so `subject_id` is NOT one-per-row here the way
  # it is at MIMIC, and this is the check that says so out loud rather than
  # leaving it to be inferred from a fold table nobody builds at an apply site.
  resampling_cols = .capture_check(function() check_resampling_cols(tabs, cfg, strict = strict_checks)),
  excursion_sides = .capture_check(function() check_excursion_sides(tabs, cfg, strict = strict_checks)),
  signal_tails    = .capture_check(function() check_signal_tails(tabs, cfg, strict = strict_checks)),
  ordinal_scales  = .capture_check(function() check_ordinal_scales(tabs, cfg, strict = strict_checks)),
  agent_pool      = .capture_check(function() check_agent_pool(tabs, cfg, strict = strict_checks))
)
chk_tab <- do.call(rbind, lapply(names(site_checks), function(nm) data.frame(
  check = nm, ok = isTRUE(site_checks[[nm]]$ok),
  outcome = site_checks[[nm]]$outcome, strict = strict_checks,
  detail = substr(site_checks[[nm]]$msg, 1, 400), stringsAsFactors = FALSE)))
save_table(run, chk_tab, "site_checks", subdir = "diagnostics")
arm_status$site_checks <- sprintf("%d/%d passed; %d warning(s); %d error(s); strict=%s",
                                  sum(chk_tab$ok), nrow(chk_tab),
                                  sum(chk_tab$outcome == "warning"),
                                  sum(chk_tab$outcome == "error"), strict_checks)

# THE TABLE IS WRITTEN, THEN THE POLICY IS APPLIED. Strict mode: any failure
# stops here, before a model is applied. Either mode: an exception a check
# could not evaluate past is not a finding about eICU and stops here too.
if (any(chk_tab$outcome == "error")) {
  bad <- chk_tab[chk_tab$outcome == "error", , drop = FALSE]
  reason <- paste0(bad$check, ": ", bad$detail, collapse = " || ")
  fail_run(run, "site_checks", reason, extra = manifest_core())
  stop("site checks: ", nrow(bad), " check(s) raised", if (strict_checks) " under design_checks_strict: true",
       ". The table is in diagnostics/site_checks; no model was applied.\n  ",
       reason, call. = FALSE)
}

# --- prior-fit diagnostics, at eICU -----------------------------------------
# `config/external.yml` has declared `checks.prior_fit_diagnostics: true` since
# it was written, and until 2026-09-03 nothing read it and nothing computed the
# check (audit finding F1). This is that check.
#
# WHAT IT ASKS. The bundle carries alpha, and the fitted delta and lambda
# conditional means, all frozen at MIMIC. Those parameters are now being
# EVALUATED ON eICU's COVERAGE DISTRIBUTION, which is not the distribution they
# were fitted against: `a2` is the log(n) coverage adjustment and is the least
# fold-stable parameter at MIMIC, so it is the first thing that would fail
# here. Differencing these tables between the MIMIC run and this one is the
# whole of the prior-transport reading, and the `lambda_fit` table in
# particular is the cleanest cross-site quantity the project has.
#
# NOTHING HERE IS FITTED IN THE SENSE HARD RULE 8 FORBIDS. The parameters come
# out of the bundle and are not re-estimated. `delta_fit_diagnostics()` does
# fit one small smooth of a residual on log1p(k) -- that is the diagnostic
# itself, its basis is declared in `diagnostics.decoupling_smooth`, and its
# output reaches no score, no L and no bundle. Same standing as
# `check_smooth_k()`, which also computes at an apply site in order to report.
#
# The measured-row set is built by `measured_rows()`, the SAME function
# `.measured_train_rows()` is built from, so the column scoping cannot diverge
# between the two sites.
if (isTRUE(ec$checks$prior_fit_diagnostics %||% TRUE)) {
  s_site <- measured_rows(tabs, cfg)
  pf <- run_stage(run, "prior_fit_diagnostics", list(
    dm_shrinkage      = dm_shrinkage_table(bundle$priors$signal, s_site, role = "final"),
    # The DM's goodness of fit and its drift across coverage, evaluated on
    # eICU's own counts under the FROZEN MIMIC alphas. Wired 2026-09-05: both
    # existed and neither runner called them, so the surviving construct's fit
    # had never been checked at the validation site. Differencing these two
    # tables against the internal run is the DM half of the prior-transport
    # reading, exactly as it is for delta and lambda.
    dm_ppc            = dm_ppc_all(bundle$priors$signal, s_site,
                                   role = "final", seed = cfg$seed),
    dm_drift          = dm_dispersion_drift(bundle$priors$signal, s_site, cfg,
                                            role = "final"),
    delta_fit         = delta_fit_diagnostics(bundle$priors$magnitude, s_site,
                                              cfg, role = "final"),
    delta_ordinal_fit = delta_ordinal_diagnostics(bundle$priors$magnitude, s_site,
                                                  cfg, role = "final",
                                                  seed = cfg$seed),
    lambda_fit        = lambda_fit_diagnostics(bundle$priors$intervention,
                                               role = "final",
                                               ivf = tabs$intervention_features,
                                               cfg = cfg)),
    extra = manifest_core())
  for (nm in names(pf)) {
    if (!is.null(pf[[nm]]) && nrow(pf[[nm]])) {
      save_table(run, pf[[nm]], paste0("prior_", nm), subdir = "diagnostics")
    }
  }
  log_msg(run, sprintf(paste0("prior-fit diagnostics: %d shrinkage row(s), %d PPC ",
                             "cell(s), %d drift row(s), %d delta row(s), %d lambda row(s)"),
    NROW(pf$dm_shrinkage), NROW(pf$dm_ppc), NROW(pf$dm_drift),
    NROW(pf$delta_fit), NROW(pf$lambda_fit)))
} else {
  pf <- NULL
  log_msg(run, "prior-fit diagnostics: SKIPPED by config (checks.prior_fit_diagnostics)")
}

# --- apply ------------------------------------------------------------------
stay_ids <- sort(unique(tabs$cohort$stay_id))
y <- as.integer(tabs$cohort$mortality[match(as.character(stay_ids),
                                            as.character(tabs$cohort$stay_id))])
if (anyNA(y)) stop("external: a stay has no mortality value", call. = FALSE)
p_bar_site <- mean(y)
# The bootstrap unit (statistical review S4): the patient of every eICU stay,
# resolved through the same accessor the training folds use, so repeat stays
# of one patient are drawn together. Never printed. Hospital-level inference is
# a separate question and is answered by the heterogeneity table below, not by
# this interval.
grp_site <- patient_group_of(tabs$cohort, cfg, stay_ids)

manifest_cohort <- function() list(n_stays = length(stay_ids), n_events = sum(y),
                                   p_bar_site = p_bar_site)

log_msg(run, sprintf("scoring %d stays, %d deaths, event rate %.4f (MIMIC train cohort rate %.4f)",
                     length(stay_ids), sum(y), p_bar_site,
                     bundle$train_ref$p_bar_cohort %||% NA_real_))

t0 <- Sys.time()
ap <- run_stage(run, "apply_bundle",
                apply_bundle(bundle, tabs, cfg, stay_ids, arms = arms, verbose = FALSE),
                extra = c(manifest_core(), manifest_cohort()))
log_msg(run, sprintf("apply_bundle: %.1f min",
                     as.numeric(difftime(Sys.time(), t0, units = "mins"))))

save_object(run, ap$l_long, "l_eicu")
save_object(run, ap$scores, "scores_eicu")
save_table(run, ap$coverage, "layer1_coverage", subdir = "diagnostics")
if (!is.null(ap$arm_table)) save_table(run, ap$arm_table, "arm_design", subdir = "diagnostics")

# --- score ------------------------------------------------------------------
# `p_bar` is eICU's OWN event rate, and only for the calibration reference
# intercept: it is the prevalence the slope's intercept is read against, not a
# quantity anything is centred on. Every score arriving here is already
# centred on a MIMIC TRAINING prior, carried in the bundle and never re-derived
# (hard rule 8). The gap between the two rates is itself a finding and is
# printed below.
#
# WHAT THE DASHED REFERENCE IS, THEN. For an LLR arm, `logit(p_bar_site) +
# score` is the Bayes line a perfectly calibrated LLR would sit on at THIS
# prevalence. For a booster arm, whose score is `logit(p_hat) -
# logit(p_bar_train)`, the same line is a TARGET-PREVALENCE-ADJUSTED reference,
# not the frozen booster's own probability: a centred score of zero maps to
# the training prevalence in the booster and to the site prevalence on the
# line. The slope does not depend on `p_bar`; only `expected_intercept` does,
# and the calibration table records `p_bar_reference` so the two cannot be
# confused.
sc_f <- if (rc$frozen)
  run_stage(run, "score_arms_frozen",
            score_arms(run, ap$scores, y, p_bar_site, breaks = bundle$cutpoints,
                       suffix = "_frozen", n_bins = nb, n_boot = nbt, seed = cfg$seed,
                       group = grp_site),
            extra = c(manifest_core(), manifest_cohort())) else NULL
sc_s <- if (rc$self)
  run_stage(run, "score_arms_self",
            score_arms(run, ap$scores, y, p_bar_site, breaks = NULL, suffix = "_self",
                       n_bins = nb, n_boot = nbt, seed = cfg$seed, group = grp_site),
            extra = c(manifest_core(), manifest_cohort())) else NULL

summ <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(sc_f)) cbind(binning = "frozen", sc_f$summary),
  if (!is.null(sc_s)) cbind(binning = "self",   sc_s$summary))))
save_table(run, summ, "score_summary")

ct <- NULL
if (length(rc$ladder)) {
  ct <- run_stage(run, "arm_contrasts",
                  arm_contrasts(ap$scores, y, rc$ladder, n_boot = nbt, seed = cfg$seed,
                                group = grp_site),
                  extra = c(manifest_core(), manifest_cohort()))
  save_table(run, ct, "arm_contrasts")
} else {
  log_msg(run, "arm contrasts: the selected arms support no ladder contrast; none written")
}

# --- transport --------------------------------------------------------------
# eICU against the training site's own out-of-fold numbers, which the bundle
# carries. The like-for-like held-out comparator is MIMIC TEST, produced by
# run/test_look.R; `compare_to` in config/external.yml names that run and the
# join is done in analysis, not here, so this runner stays free of any
# dependency on a second run directory existing.
#
# THE PRIMARY COLUMNS ARE THE TWO AUROCS AND THEIR DIFFERENCE. `auroc_ratio`
# is the raw quotient the table used to call `retained`; it credits chance,
# so a model at 0.5 against a training 0.8 reads as 0.625.
# `auroc_retained_above_chance` measures both sides from 0.5 and is NA when
# the training AUROC is within `RETENTION_MIN_EXCESS` of chance (review E8).
base <- if (!is.null(sc_s)) sc_s$summary else sc_f$summary
tr <- bundle$train_ref
transport <- merge(
  data.frame(arm = sub("_(self|frozen)$", "", base$label),
             auroc_eicu = base$auroc, auroc_lo_eicu = base$auroc_lo,
             auroc_hi_eicu = base$auroc_hi, auprc_eicu = base$auprc,
             cal_slope_eicu = base$cal_slope, stringsAsFactors = FALSE),
  data.frame(arm = tr$arms$label, auroc_mimic_oof = tr$arms$auroc,
             auprc_mimic_oof = tr$arms$auprc, stringsAsFactors = FALSE),
  by = "arm", all.x = TRUE)
transport$d_auroc <- round(transport$auroc_eicu - transport$auroc_mimic_oof, 5)
transport$d_auprc <- round(transport$auprc_eicu - transport$auprc_mimic_oof, 5)
transport$auroc_ratio <- round(transport$auroc_eicu / transport$auroc_mimic_oof, 4)
transport$auroc_retained_above_chance <-
  retention_above_chance(transport$auroc_eicu, transport$auroc_mimic_oof)
transport$binning_of_eicu_row <- if (!is.null(sc_s)) "self" else "frozen"
save_table(run, transport, "transport")

# --- the severity arm -------------------------------------------------------
# APACHE II and SOFA, through `severity_arm()` -- the SAME function
# run/test_look.R calls, so the two sites differ only in which rows they hand it.
# Nothing is fitted: the point tables are frozen in code and the recalibration
# intercept and slope come out of the bundle (hard rule 8).
#
# The coverage floors restrict every cell, ours included, so `llr_sum` appears
# twice in this run on two different row sets. The gap between them is the size
# of the selection the floors introduce and is reported rather than assumed away.
#
# FOUR STATES, NOT TWO (review E1). The arm is `disabled` by config;
# `unavailable` when the bundle has no severity slot or the site configures no
# severity table, both of which are declared conditions and are recorded; or
# it runs. If it runs and raises -- a recalibration fault, a row misalignment,
# a plotting error, a defect -- the error is NOT swallowed: `run_stage()`
# stamps the manifest `failed` with the stage and reason and re-raises. A
# manifest that says `complete` therefore means every requested arm finished.
# A configured severity path whose file is missing already fails in
# `load_tables()`, which is the consistent policy: a declared input that is
# absent is an error, an undeclared one is an unavailable arm.
sa <- NULL
arm_status$severity <- {
  if (!isTRUE(ec$severity$enabled %||% TRUE)) "disabled"
  else if (is.null(bundle$severity)) "unavailable: bundle carries no severity slot"
  else if (is.null(tabs$severity)) "unavailable: no severity table configured for this site"
  else "requested"
}
if (arm_status$severity == "requested") {
  sa <- run_stage(run, "severity_arm",
    severity_arm(run, bundle, tabs, cfg, stay_ids, y, ap$scores,
                 l_full = ap$l_mats$full, domains = bundle$domains,
                 n_bins = nb, n_boot = nbt, seed = cfg$seed, verbose = FALSE,
                 group = grp_site),
    extra = c(manifest_core(), manifest_cohort(),
              list(auroc = as.list(stats::setNames(base$auroc, base$label)))))
  arm_status$severity <- "complete"
} else {
  log_msg(run, "severity arm: ", arm_status$severity)
}

# --- between-hospital heterogeneity ----------------------------------------
# The result eICU is uniquely able to supply. The pooled AUROC and the median
# within-hospital AUROC answer different questions: a score can discriminate
# well overall partly because hospitals differ in case mix while discriminating
# less well inside any one of them. `pooled_minus_median` puts the two side by
# side. It is a DESCRIPTIVE CONTRAST between two weightings of the eligible
# cohort, not a decomposition into within- and between-hospital parts (review
# E8); see `group_metrics()`'s header for the counterexample.
#
# THE JOIN IS CHECKED, NOT ASSUMED (review E4). The hospital table must carry
# `stay_id` and the grouping column; a stay key that appears twice with two
# different hospitals is a conflict and stops the run; a duplicate with one
# hospital is de-duplicated and counted; a missing hospital value is an error.
# A scored stay the table does not cover is handled by the DECLARED
# `unmatched_policy`, and the summary reports the scored, matched and eligible
# denominators separately so a coverage fraction cannot hide an exclusion.
# Zero eligible hospitals is a result with a `status`, not an infinite
# extremum. An enabled analysis whose table is missing is an error: the config
# declared the file.
hosp_summ <- NULL; hosp_meta <- NULL
arm_status$hospital <- {
  if (!isTRUE(ec$hospital$enabled)) "disabled"
  else if (no_hosp) "skipped: --no-hospital" else "requested"
}
if (arm_status$hospital == "requested") {
  hp <- ec$hospital
  hosp <- run_stage(run, "hospital_heterogeneity", {
    if (!file.exists(ec$hospital_table)) {
      stop("hospital analysis is enabled but `hospital_table` does not exist: ",
           ec$hospital_table, call. = FALSE)
    }
    hz <- as.data.frame(arrow::read_parquet(ec$hospital_table))
    gcol <- hp$group_col %||% "hospitalid"
    miss <- setdiff(c("stay_id", gcol), names(hz))
    if (length(miss)) abort_values("external: hospital table lacks column(s)", miss)

    ids_h <- as.character(hz$stay_id)
    n_na_grp <- sum(is.na(hz[[gcol]]))
    if (n_na_grp) {
      stop(sprintf("hospital table: %d row(s) have a missing `%s`", n_na_grp, gcol),
           call. = FALSE)
    }
    dup <- duplicated(ids_h) | duplicated(ids_h, fromLast = TRUE)
    n_dup_keys <- length(unique(ids_h[dup]))
    if (n_dup_keys) {
      conflict <- tapply(as.character(hz[[gcol]][dup]), ids_h[dup],
                         function(v) length(unique(v)) > 1L)
      n_conf <- sum(conflict)
      if (n_conf) {
        stop(sprintf(paste0("hospital table: %d stay key(s) appear more than once and %d ",
                            "of them carry conflicting `%s` assignments"),
                     n_dup_keys, n_conf, gcol), call. = FALSE)
      }
      log_msg(run, sprintf("hospital table: %d stay key(s) duplicated with one consistent ",
                           "%s each; de-duplicated", n_dup_keys, gcol))
      hz <- hz[!duplicated(ids_h), , drop = FALSE]; ids_h <- as.character(hz$stay_id)
    }

    grp <- hz[[gcol]][match(as.character(stay_ids), ids_h)]
    n_unmatched <- sum(is.na(grp))
    n_not_scored <- sum(!ids_h %in% as.character(stay_ids))
    log_msg(run, sprintf(paste0("hospital table: %d row(s), %d distinct %s; %d scored stay(s) ",
                                "unmatched (policy: %s); %d table stay(s) not in the scored cohort"),
                         nrow(hz), length(unique(hz[[gcol]])), gcol, n_unmatched,
                         hp$unmatched_policy, n_not_scored))
    if (n_unmatched && identical(hp$unmatched_policy, "error")) {
      stop(sprintf("hospital table: %d scored stay(s) have no hospital and ",
                   "`unmatched_policy` is `error`", n_unmatched), call. = FALSE)
    }

    gm <- lapply(names(ap$scores), function(nm)
      group_metrics(ap$scores[[nm]], y, grp, label = nm,
                    min_n = as.integer(hp$min_stays), min_events = as.integer(hp$min_events)))
    names(gm) <- names(ap$scores)
    list(summary = do.call(rbind, lapply(gm, function(z) z$summary)),
         per_group = lapply(gm, function(z) z$per_group),
         meta = list(n_table_rows = nrow(hz), n_unmatched = n_unmatched,
                     n_table_stays_not_scored = n_not_scored,
                     n_duplicate_keys = n_dup_keys,
                     unmatched_policy = hp$unmatched_policy,
                     min_stays = as.integer(hp$min_stays),
                     min_events = as.integer(hp$min_events), group_col = gcol))
  }, extra = c(manifest_core(), manifest_cohort()))
  hosp_summ <- hosp$summary; hosp_meta <- hosp$meta
  save_table(run, hosp_summ, "hospital_summary")
  # Row-level in the sense that matters: a hospital id is an identifier under
  # the DUA. Written as an object, never printed, never a CSV.
  save_object(run, hosp$per_group, "hospital_per_group")
  arm_status$hospital <- if (all(hosp_summ$status == "ok")) "complete"
                         else paste0("complete: ", paste(unique(hosp_summ$status), collapse = "/"))
}

# --- report -----------------------------------------------------------------
cat("\n=== eICU EXTERNAL VALIDATION: does the evidence geometry transport? ===\n\n")
cat(sprintf("  bundle        %s (fitted at %s; md5 %s)\n", bundle_path, bundle$site,
            substr(run_cfg$bundle$md5, 1, 12)))
cat(sprintf("  eICU          %d stays, %d deaths (%.2f%%)\n",
            length(stay_ids), sum(y), 100 * p_bar_site))
cat(sprintf("  MIMIC train   %d stays, cohort event rate %.2f%%\n",
            bundle$train_ref$n %||% NA_integer_,
            100 * (bundle$train_ref$p_bar_cohort %||% NA_real_)))
cat("  Every score below is centred on a MIMIC TRAINING prior carried in the\n")
cat("  bundle. None of it was re-derived here (hard rule 8).\n\n")

cat("=== site checks: declarations that held at MIMIC, tested here ===\n\n")
print(chk_tab[, c("check", "ok", "outcome", "detail")], row.names = FALSE)
cat(sprintf("\n  strict = %s. A `warning` is a declaration that did not hold here and is a\n",
            strict_checks))
cat("  FINDING; an `error` would have stopped this run before any model was\n")
cat("  applied. `agent_pool` is the one to read first: counting drug names\n")
cat("  rather than molecules inflates eICU n_agents two- to fourfold and would\n")
cat("  read as case mix.\n")

if (!is.null(pf)) {
  cat("\n=== PRIOR TRANSPORT: the frozen priors, evaluated on eICU coverage ===\n")
  cat("  Nothing below is re-fitted. These are the MIMIC-fitted parameters read\n")
  cat("  against eICU's own coverage distribution, which is the distribution\n")
  cat("  they were NOT fitted against. Difference each row against the same row\n")
  cat("  in the internal run's diagnostics to get the transport reading.\n\n")
  if (!is.null(pf$dm_shrinkage) && nrow(pf$dm_shrinkage)) {
    cat("  Dirichlet shrinkage: the posterior weight on a stay's OWN counts at\n")
    cat("  each quartile of coverage. A weight that moves between sites means\n")
    cat("  the same alpha is shrinking differently because coverage differs.\n\n")
    print(pf$dm_shrinkage[, c("signal", "n_stays", "alpha0", "n_med",
                              "w_q25", "w_med", "w_q75")], row.names = FALSE)
  }
  if (!is.null(pf$lambda_fit) && nrow(pf$lambda_fit)) {
    cat("\n  lambda, the intensity conditional mean. The cleanest cross-site\n")
    cat("  quantity the project has: same parameters, different exposure mix.\n\n")
    print(pf$lambda_fit[, c("intervention", "family", "n_stays", "n_exposed_site",
                            "c0", "c1", "b0", "b1", "resid_mean", "resid_sd",
                            "noise_identified", "spread_rho_iqr")], row.names = FALSE)
  }
  if (!is.null(pf$delta_fit) && nrow(pf$delta_fit)) {
    cat("\n  delta, the magnitude conditional mean. `resid_smooth_r2` is the\n")
    cat("  leftover count dependence: near zero means the frozen conditional\n")
    cat("  mean still absorbs the count relationship here. `a2` is the log(n)\n")
    cat("  coverage adjustment and is the parameter most exposed to a coverage\n")
    cat("  shift, so read it first. Worst 10 rows by leftover:\n\n")
    print(utils::head(pf$delta_fit[, c("signal", "variable", "form", "n_stays",
                                       "a1", "a2", "resid_smooth_r2",
                                       "resid_step_r2")], 10), row.names = FALSE)
  }
}

cat("\n=== discrimination at eICU, all arms ===\n")
cat("  `auroc_lo`/`auroc_hi` resample PATIENTS (`boot_unit`). `cal_slope` is a\n")
cat("  GLM slope whose interval, in the saved calibration tables, is a\n")
cat("  stay-level profile interval; the bin intervals are stay-level Wilson.\n")
cat("  Those two are labelled `ci_method` in their tables and are NOT\n")
cat("  patient-clustered.\n\n")
print(base[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi", "boot_unit",
               "auprc", "auprc_lift", "spearman", "cal_slope")], row.names = FALSE)

cat("\n=== TRANSPORT: eICU against MIMIC out-of-fold ===\n")
cat("  Read `d_auroc` first. `auroc_ratio` is the raw quotient and credits\n")
cat("  chance (0.5 against 0.8 would read 0.625); `auroc_retained_above_chance`\n")
cat(sprintf("  measures both from 0.5 and is NA when MIMIC is within %.2f of chance.\n",
            RETENTION_MIN_EXCESS))
cat("  Read beside MIMIC TEST from run/test_look.R -- out-of-fold train and a\n")
cat("  fully external cohort are two different kinds of held-out.\n\n")
print(transport[, c("arm", "auroc_mimic_oof", "auroc_eicu", "d_auroc", "auroc_ratio",
                    "auroc_retained_above_chance", "d_auprc", "cal_slope_eicu")],
      row.names = FALSE)

if (!is.null(ct)) {
  cat("\n=== paired contrasts at eICU, identical rows ===\n")
  cat("  Interval and `auroc_p` are the paired bootstrap resampled by `boot_unit`\n")
  cat("  (patients). DeLong is the observation-level iid reference and is printed\n")
  cat("  separately below; the two rest on different assumptions.\n\n")
  print(ct[, c("a", "b", "d_auroc", "auroc_lo", "auroc_hi", "auroc_p", "boot_unit",
               "d_auprc", "auprc_lo", "auprc_hi", "auprc_p")], row.names = FALSE)
  cat("\n  DeLong, iid sensitivity (not cluster-adjusted):\n\n")
  print(ct[, c("a", "b", "d_auroc", "delong_se", "delong_z", "delong_p")], row.names = FALSE)

  cat("\n=== DO THE INTERACTIONS TRANSPORT? the ladder at eICU ===\n\n")
  rung <- function(a, b, label) {
    v <- ct$d_auroc[ct$a == a & ct$b == b]
    if (length(v) == 1L) cat(sprintf("  %-24s %+.4f   (%s - %s)\n", label, v, a, b))
    else cat(sprintf("  %-24s    n/a    (%s or %s not scored)\n", label, a, b))
  }
  rung("xgb_feat", "xgb_raw",  "covariate construction")
  rung("xgb_l",    "xgb_feat", "per-signal collapse")
  rung("llr_sum",  "xgb_l",    "linear aggregation")
  rung("llr_sum",  "xgb_raw",  "method vs strong learner")
  cat("\n  Compare each rung with the same rung at MIMIC. A rung that is large\n")
  cat("  there and near zero here is CONSISTENT WITH that rung having fitted\n")
  cat("  MIMIC-specific structure. It is a hypothesis the two tables raise, not\n")
  cat("  a proof: coverage, case mix and charting also differ between sites.\n")
} else {
  cat("\n=== paired contrasts: none reportable for the selected arms ===\n")
}

if (!is.null(sa)) report_severity_arm(sa, "eICU")

if (!is.null(hosp_summ)) {
  cat(sprintf(paste0("\n=== BETWEEN-HOSPITAL HETEROGENEITY (%d hospitals, %d reported; ",
                     "%d of %d scored stays matched, policy %s) ===\n\n"),
              hosp_summ[["n_groups"]][1], hosp_summ[["n_reported"]][1],
              hosp_summ[["n_matched"]][1], hosp_summ[["n_scored"]][1],
              hosp_meta$unmatched_policy))
  print(hosp_summ[, c("label", "status", "n_reported", "frac_kept_of_scored",
                      "frac_kept_of_matched", "pooled_auroc_matched",
                      "pooled_auroc_eligible", "auroc_median", "auroc_q1", "auroc_q3",
                      "auroc_iqr", "pooled_minus_median")], row.names = FALSE)
  cat("\n  `pooled_minus_median` is eligible-pooled minus the median over the same\n")
  cat("  eligible hospitals: a DESCRIPTIVE CONTRAST between two weightings of\n")
  cat("  one cohort. It is not a case-mix decomposition -- selection and unequal\n")
  cat("  within-hospital discrimination move it too. `frac_kept_of_scored` is\n")
  cat("  the coverage of the scored cohort; `_of_matched` is eligibility among\n")
  cat("  matched stays. Hospital identifiers are not printed.\n")
  # THE SAME DISTRIBUTION ON THE PRECISION-RECALL SIDE, added 2026-09-07. AUPRC
  # is summarised as a LIFT over each hospital's own event rate because the
  # event rate is exactly what differs between them -- a median raw AUPRC over
  # 166 hospitals is a median over 166 different floors, and a hospital would
  # rank high purely by being sicker. See group_metrics()'s header.
  cat("\n=== the same, on the precision-recall side (lift over each ",
      "hospital's own event rate) ===\n\n", sep = "")
  print(hosp_summ[, c("label", "pooled_event_rate_eligible", "pooled_auprc_eligible",
                      "pooled_auprc_lift_eligible", "auprc_lift_median", "auprc_lift_q1",
                      "auprc_lift_q3", "auprc_lift_iqr",
                      "pooled_minus_median_lift")], row.names = FALSE)
  cat("\n  A lift of 1.0 is no better than the event rate at that hospital.\n")
  cat("  `pooled_auprc_lift_eligible` is over the POOLED eligible rate, which is\n")
  cat("  not the mean of the per-hospital floors, so `pooled_minus_median_lift`\n")
  cat("  mixes two things and is a description rather than a decomposition.\n")
}

cat("\n=== layer-1 coverage at eICU ===\n")
cat("  `frac_scored` is the measured fraction per signal. A signal much thinner\n")
cat("  here than at MIMIC contributes less evidence, and that is a coverage\n")
cat("  difference rather than a transport failure. Read the two side by side.\n\n")
print(ap$coverage[, c("signal", "model", "n_scored", "frac_scored",
                      "l_mean", "l_sd")], row.names = FALSE)

finalize_run(run, extra = c(manifest_core(), manifest_cohort(), list(
  auroc = as.list(stats::setNames(base$auroc, base$label)),
  hospital = if (is.null(hosp_summ)) list(status = arm_status$hospital) else c(
    list(status = arm_status$hospital,
         n_groups = hosp_summ$n_groups[1], n_reported = hosp_summ$n_reported[1],
         n_scored = hosp_summ$n_scored[1], n_matched = hosp_summ$n_matched[1]),
    hosp_meta),
  severity_scored = !is.null(sa),
  severity_n_kept = if (is.null(sa)) NA_integer_ else sum(sa$keep),
  site_checks_passed = sum(chk_tab$ok), site_checks_total = nrow(chk_tab),
  site_checks_strict = strict_checks)))

cat(sprintf("\n  run directory: %s\n\n", run$path))
