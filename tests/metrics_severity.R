# tests/metrics_severity.R ---------------------------------------------------
# The severity-score comparison arms, APACHE and SOFA. Reads a finished
# branch_point run, fits no layer-1 model, and scores every cell through the
# SAME code path as the proposed method (R/09's score_report) so the numbers are
# commensurable.
#
# WHAT THESE ARMS ARE FOR, and how they differ from the XGBoost one. R/09b
# varies the learner and the aggregation while holding the measurements fixed.
# This file varies WHO DESIGNED THE FEATURES. APACHE II and SOFA see the same
# measurements over the same 24-hour window and reduce each signal to a worst
# value mapped through a published point table. We shrink a proportion, fit a
# smooth and take a log-likelihood ratio. The gap is attributable to the
# construction rather than to the data, the window or the estimator, which is
# the one thing no other arm in the project isolates.
#
# TWO SCORES, NOT ONE, AND THEY ARE NOT REDUNDANT. APACHE II is pure
# physiology, so its honest counterpart is `L_meas`. SOFA folds vasopressor
# dose into its cardiovascular component, so its honest counterpart is
# `L_full`, which carries `L_intv`. The two externalise the project's own
# meas/full split rather than duplicating each other. That is the specific
# reason docs/v2_analysis_tiering.md item 5 -- one severity score, not both --
# does not settle this: it is a rule about discrimination baselines and it is
# correct about those.
#
# THE CELLS.
#   llr_sum        the naive summed L_full             <- the proposed method
#   llr_meas       the naive summed L_meas             <- physiology only
#   apache2_aps    recomputed APACHE II acute physiology score, physiology only
#   apache2_total  APACHE II total: APS + age + chronic health
#   aps_native     the site's own published APS (APS III here, APACHE IVa at
#                  eICU), taken whole from the derived table
#   sofa           recomputed SOFA total               <- bundles intervention
#   sofa_native    the derived SOFA concept (MIMIC only; eICU ships none)
#   llr_plus_demo  llr_sum with age, chronic health and admission class added
#
# READING THEM. The first two are the claims; the rest are the context a
# reviewer will ask for.
#   llr_meas vs apache2_aps     THE APACHE CLAIM. Both are physiology-only,
#                               both are built on the same measurements, and
#                               neither carries age or comorbidity.
#   llr_sum  vs sofa            THE SOFA CLAIM. Both bundle measurement with
#                               intervention. SOFA does it for one organ, ad
#                               hoc; the design does it for twelve pairs,
#                               systematically. This is the comparison the
#                               representational argument actually rests on.
#   llr_sum vs apache2_total    the honest external-comparability number. We
#                               lose the age and comorbidity APACHE carries,
#                               and a reviewer will look for this one first.
#   llr_plus_demo vs apache2_total  what our excluded covariates are worth, and
#                               whether excluding them costs discrimination.
#   llr_sum vs aps_native       against the site's canonical implementation
#                               rather than our recomputation. Guards against
#                               the whole arm resting on our point table.
#   sofa vs sofa_native         VALIDATION, not a result. Our SOFA against a
#                               published SOFA, per organ. This gates the arm.
#
# AND THE DOMAIN TABLE, which is what SOFA buys that nothing else does. Layer 2
# emits one number per domain and until now that number could only be compared
# with itself. config/domains.csv already carries a `sofa_organ` column, frozen
# before any result existed, so six of the eleven domains get an externally
# defined referent on organ definitions the project already committed to. The
# five domains SOFA has no organ for are reported alongside, because "our
# domains beat SOFA's organs" is a much weaker claim than "and there are five
# more domains SOFA cannot express at all".
#
# THREE ASYMMETRIES, and all three favour the proposed method. Say so.
#   1. A missing APACHE II variable or SOFA organ scores ZERO, so thin coverage
#      makes both baselines look healthy rather than unmeasured. That is the
#      standard retrospective convention and it weakens them. The primary
#      analysis therefore RESTRICTS to stays with at least
#      `config/apache.min_vars_present` of the twelve variables AND
#      `config/sofa.min_organs_present` of the six organs; the unrestricted
#      numbers are printed beside it, never instead of it.
#   2. APACHE II is scored once, on the whole cohort, with no fitting at all.
#      Our L's are cross-fitted. The baseline gets no chance to overfit, which
#      is a real advantage to it, but the point scale is fixed at 1985 values
#      and never adapted to this cohort, which is a real disadvantage.
#      `apache2_aps` is therefore ALSO reported after a cross-fitted logistic
#      recalibration onto the log-odds scale, which is what makes the
#      calibration column comparable. AUROC is a rank statistic, so
#      recalibration cannot move it: any AUROC difference is the score, not the
#      rescaling.
#   3. The oxygenation variable is missing for every stay with no arterial gas,
#      which is disproportionately the less sick. See the coverage table.
#
# PAIRED TESTS, NOT OVERLAPPING INTERVALS. Every cell is scored on identical
# rows, so the AUROCs are strongly correlated and marginal bootstrap intervals
# overlap long after the difference is reliable. DeLong on AUROC and a paired
# bootstrap on AUPRC settle it; overlap does not.
#
# AGGREGATES ONLY (hard rule 1). Scores are row-level and never printed.
#
#   Rscript tests/metrics_severity.R                        # newest branch run
#   Rscript tests/metrics_severity.R out/runs/branch_...    # a specific one
#   Rscript tests/metrics_severity.R - --all-stays          # skip the coverage
#                                                             restriction
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args     <- commandArgs(trailingOnly = TRUE)
src_dir  <- if (length(args) >= 1L && args[1] != "-") args[1] else latest_run("branch")
all_stays <- "--all-stays" %in% args
if (is.null(src_dir) || !dir.exists(src_dir)) {
  stop("no completed branch run found. Run tests/branch_point.R first, or pass a run directory.",
       call. = FALSE)
}
l_path <- file.path(src_dir, "tables", "l_oof.rds")
if (!file.exists(l_path)) stop("no l_oof.rds in ", src_dir, call. = FALSE)

cfg <- load_config("config/config.yml")
if (is.null(cfg$paths$mimiciv$severity)) {
  stop("config/config.yml has no paths.mimiciv.severity. Run ",
       "sql/mimiciv/v2_06_severity_mimiciv.sql, export the parquet, then ",
       "uncomment the path.", call. = FALSE)
}
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
if (is.null(tabs$severity)) stop("severity table did not load", call. = FALSE)
folds <- assign_folds(tabs$cohort, cfg)

l_long   <- readRDS(l_path)
stay_ids <- sort(unique(l_long$stay_id))
stay_ids <- intersect(folds$stay_id[folds$split == "train"], stay_ids)

# One alignment, used by every cell. Everything below indexes off `stay_ids`,
# so no cell can silently be scored on a different row set from another - which
# is the assumption every paired test in this file depends on.
ci <- match(as.character(stay_ids), as.character(tabs$cohort$stay_id))
ai <- match(as.character(stay_ids), as.character(tabs$severity$stay_id))
if (anyNA(ci)) stop("metrics_severity: a scored stay has no cohort row", call. = FALSE)
if (anyNA(ai)) stop("metrics_severity: a scored stay has no severity row", call. = FALSE)

y       <- tabs$cohort$mortality[ci]
age     <- tabs$cohort$age[ci]
ap      <- tabs$severity[ai, , drop = FALSE]
fold_of <- folds$fold[match(stay_ids, folds$stay_id)]
if (anyNA(fold_of)) stop("metrics_severity: a scored stay has no fold", call. = FALSE)

min_vars   <- cfg$apache$min_vars_present %||% 10L
arf        <- isTRUE(cfg$apache$arf_doubling)
gcs_src    <- cfg$apache$gcs_source %||% "native"
min_organs <- cfg$sofa$min_organs_present %||% 4L
resp_sup   <- cfg$sofa$resp_support %||% "invasive"
sofa_organs <- as.character(unlist(cfg$sofa$organs %||% SOFA_ORGANS))

# --- score both severity scores ---------------------------------------------
# Same table, same rows, and every shared input read from the same column. Any
# difference between the two arms is therefore construction and never input.
ap2 <- apache2_score(ap, age = age, arf_doubling = arf, gcs_source = gcs_src)
sf  <- sofa_score(ap, resp_support = resp_sup, organs = sofa_organs)

# The restriction. Applied to EVERY cell, not only the baseline ones: comparing
# a restricted baseline against an unrestricted proposed method would be a
# different and much weaker experiment. Both coverage rules apply, so the
# retained set is scoreable under BOTH constructions -- which is what makes the
# APACHE and SOFA cells comparable to each other and not only to llr_sum.
# `sf` rather than `ap` for the SOFA mask: the extraction's
# `sofa_n_organs_present` counts six, while `sofa_score()` recomputes it over
# the declared organ set. Using the raw column would apply a 4-of-6 floor where
# the score is five organs.
keep <- if (all_stays) rep(TRUE, length(stay_ids)) else
  (apache_complete_mask(ap, min_vars = min_vars) &
   sofa_complete_mask(sf, min_organs = min_organs))

# THE MASK MUST BE THE LENGTH OF THE COHORT. A shorter one recycles and a
# zero-length one empties everything while reporting "0 dropped" -- and because
# `all(is.na(numeric(0)))` is TRUE, every cell then declares itself entirely NA
# and the run dies much later on an unrelated rbind. Assert it here, once,
# where the diagnosis is one line rather than a bisect.
if (length(keep) != length(stay_ids)) {
  stop(sprintf(paste0("coverage mask is length %d against %d stays. A mask ",
                      "built from a missing column is zero-length and fails ",
                      "silently downstream."),
               length(keep), length(stay_ids)), call. = FALSE)
}
if (!any(keep)) stop("coverage mask retained no stays; check the floors", call. = FALSE)
n_drop <- sum(!keep)

run <- new_run("severity", cfg, note = sprintf(
  "APACHE + SOFA comparison arms from %s; %d of %d stays kept (min_vars %s / min_organs %s); gcs_source = %s; resp_support = %s",
  basename(src_dir), sum(keep), length(keep),
  if (all_stays) "off" else min_vars, if (all_stays) "off" else min_organs,
  gcs_src, resp_sup))
log_msg(run, sprintf("source run: %s | %d stays | %d dropped for thin APACHE coverage",
                     basename(src_dir), length(stay_ids), n_drop))
log_msg(run, sprintf("arf_doubling = %s | gcs_source = %s | resp_support = %s | native version(s): %s",
                     arf, gcs_src, resp_sup,
                     paste(levels(droplevels(ap$aps_native_version[keep])), collapse = ",")))

yk       <- y[keep]
fold_k   <- fold_of[keep]
p_bar    <- mean(yk)
nb       <- cfg$metrics$n_bins %||% 20L
nboot    <- cfg$metrics$n_boot %||% 200L

# --- the cells --------------------------------------------------------------
mats     <- l_matrices(l_long, tabs, cfg, stay_ids, fill = "zero")
llr_full <- rowSums(mats$full)[keep]
# L_meas is the physiology-only aggregate and is the honest counterpart to
# APACHE II, which carries no intervention term at all. It exists for the 19
# signals exactly as L_full does; for the 7 unpaired ones the two are identical
# by construction (CLAUDE.md, frozen decisions), so llr_meas differs from
# llr_sum only through the 12 paired signals -- which is precisely the
# difference the SOFA comparison is about.
llr_meas <- if (!is.null(mats$meas)) rowSums(mats$meas)[keep] else NULL

# Age and chronic health as a log-odds block, so that `llr_plus_demo` is
# llr_sum PLUS something rather than a refit of llr_sum. Cross-fitted on the
# same folds, so the added covariates get no in-sample advantage over the L's
# they are being added to.
demo <- data.frame(
  llr        = llr_full,
  age        = age[keep],
  chronic    = as.integer(ap$chronic_immunocompromised[keep] == 1L |
                          ap$chronic_severe_organ[keep] == 1L),
  admission  = droplevels(ap$admission_class[keep]))
llr_demo <- rep(NA_real_, nrow(demo))
for (f in sort(unique(fold_k))) {
  te <- fold_k == f
  fit <- stats::glm(yk[!te] ~ llr + age + chronic + admission,
                    data = demo[!te, , drop = FALSE], family = stats::binomial())
  llr_demo[te] <- stats::predict(fit, newdata = demo[te, , drop = FALSE], type = "link")
}
llr_demo <- llr_demo - logit(p_bar)

# The three APACHE scores, each recalibrated out-of-fold onto the log-odds
# scale so the calibration column and the risk-curve bins mean something.
# The recalibration is monotone WITHIN a fold but not globally, so it can move
# AUROC very slightly; that gap is measured below rather than assumed away.
aps2_raw   <- ap2$aps[keep]
tot2_raw   <- ap2$total[keep]
natv_raw   <- as.numeric(ap$aps_native[keep])

sofa_raw   <- sf$sofa[keep]
sofan_raw  <- as.numeric(ap$sofa_native[keep])

cells <- list(
  llr_sum       = llr_full,
  apache2_aps   = recalibrate_oof(aps2_raw, yk, fold_k, p_bar),
  apache2_total = recalibrate_oof(tot2_raw, yk, fold_k, p_bar),
  sofa          = recalibrate_oof(sofa_raw, yk, fold_k, p_bar),
  llr_plus_demo = llr_demo)
if (!is.null(llr_meas)) cells$llr_meas <- llr_meas
if (any(!is.na(natv_raw))) {
  cells$aps_native <- recalibrate_oof(natv_raw, yk, fold_k, p_bar)
} else {
  log_msg(run, "aps_native is entirely NA; the native-APS cell is skipped")
}
if (any(!is.na(sofan_raw))) {
  cells$sofa_native <- recalibrate_oof(sofan_raw, yk, fold_k, p_bar)
} else {
  # The normal state at eICU, which ships no SOFA. NOT a pass: it means the
  # recomputation has nothing to be checked against at this site, and it must
  # therefore have been validated at MIMIC before this run is believed.
  log_msg(run, "sofa_native is entirely NA (expected at eICU); validation table skipped")
}
cells <- cells[intersect(c("llr_sum", "llr_meas", "apache2_aps", "apache2_total",
                           "aps_native", "sofa", "sofa_native", "llr_plus_demo"),
                         names(cells))]

res <- list()
for (nm in names(cells)) {
  s <- cells[[nm]]
  if (all(is.na(s))) { log_msg(run, sprintf("cell %s all NA, skipped", nm)); next }
  res[[nm]] <- score_report(run, s, yk, p_bar, label = nm,
                            n_bins = nb, n_boot = nboot, seed = cfg$seed,
                            title = sprintf("Out-of-fold risk ordering - %s", nm))
}
summ <- do.call(rbind, lapply(res, function(z) z$summary))
save_table(run, summ, "score_summary")

# --- the tables that have to sit beside the numbers -------------------------
cov_tab <- apache_coverage(ap)
save_table(run, cov_tab, "apache_coverage", subdir = "diagnostics")

pts <- data.frame(
  variable = sub("^pt_", "", grep("^pt_", names(ap2), value = TRUE)),
  mean_points = round(vapply(grep("^pt_", names(ap2), value = TRUE),
                             function(v) mean(ap2[[v]][keep]), numeric(1)), 4),
  frac_nonzero = round(vapply(grep("^pt_", names(ap2), value = TRUE),
                              function(v) mean(ap2[[v]][keep] > 0), numeric(1)), 4),
  stringsAsFactors = FALSE, row.names = NULL)
save_table(run, pts, "apache2_points_by_variable", subdir = "diagnostics")

agree <- if (any(!is.na(natv_raw)))
  apache_native_agreement(tot2_raw, natv_raw) else NULL
if (!is.null(agree)) save_table(run, agree, "apache2_vs_native", subdir = "diagnostics")

# THE GATE ON THE SOFA ARM. Unlike the APACHE sanity check, which compares two
# different scores and can only ask for correlation, this compares our SOFA
# against a published SOFA and can ask for near-identity. Cardiovascular is
# expected to disagree -- the NEE collapse in .sofa_cardio() guarantees it --
# and any OTHER organ disagreeing materially means a band, a unit or a
# missing-value rule is wrong, which eICU would then inherit invisibly.
sofa_agree <- sofa_native_agreement(sf[keep, , drop = FALSE], ap[keep, , drop = FALSE],
                                    scored = sofa_organs)
if (!is.null(sofa_agree)) save_table(run, sofa_agree, "sofa_vs_native", subdir = "diagnostics")

sofa_pts <- data.frame(
  organ = sub("^sofa_", "", grep("^sofa_(resp|coag|liver|cardio|cns|renal)",
                                 names(sf), value = TRUE)),
  mean_score = round(vapply(grep("^sofa_(resp|coag|liver|cardio|cns|renal)",
                                 names(sf), value = TRUE),
                            function(v) mean(sf[[v]][keep]), numeric(1)), 4),
  frac_nonzero = round(vapply(grep("^sofa_(resp|coag|liver|cardio|cns|renal)",
                                   names(sf), value = TRUE),
                              function(v) mean(sf[[v]][keep] > 0), numeric(1)), 4),
  stringsAsFactors = FALSE, row.names = NULL)
save_table(run, sofa_pts, "sofa_score_by_organ", subdir = "diagnostics")

# The domain-level comparison. Equal weights, because layer 2's weight vector
# does not exist yet (v2_state_20260828.md section 3.4); D_k reduces to the
# unweighted sum at w = 1 and the comparison sharpens rather than changes shape
# when the weights land.
domains  <- load_domains(cfg$paths$domains %||% "config/domains.csv")
dom_tab  <- sofa_domain_table(mats$full[keep, , drop = FALSE], sf[keep, , drop = FALSE],
                              yk, domains, organs = sofa_organs)
uncov    <- sofa_uncovered_domains(mats$full[keep, , drop = FALSE], yk, domains)
save_table(run, dom_tab, "domain_vs_sofa_organ")
save_table(run, uncov,   "domains_sofa_cannot_express")

# What the out-of-fold recalibration cost in discrimination. It is monotone
# within a fold but each fold gets its own mapping, so stays in different folds
# can swap and AUROC can move. Expect a difference in the fourth decimal. If it
# is ever large enough to matter against the effect being reported, the RAW
# column is the one to quote and the recalibrated score is doing work it should
# not be.
raw_of <- list(apache2_aps = aps2_raw, apache2_total = tot2_raw, sofa = sofa_raw,
               aps_native = natv_raw, sofa_native = sofan_raw)
raw_of <- raw_of[vapply(names(raw_of), function(k)
  k %in% names(cells) && any(!is.na(raw_of[[k]])), logical(1))]
recal <- data.frame(
  score = names(raw_of),
  auroc_raw   = round(vapply(names(raw_of), function(k) .auroc(raw_of[[k]], yk), numeric(1)), 5),
  auroc_recal = round(vapply(names(raw_of), function(k) .auroc(cells[[k]], yk), numeric(1)), 5),
  stringsAsFactors = FALSE, row.names = NULL)
recal$delta <- round(recal$auroc_recal - recal$auroc_raw, 5)
save_table(run, recal, "recalibration_cost", subdir = "diagnostics")

# --- paired tests -----------------------------------------------------------
# These are the numbers that decide anything. Marginal intervals do not.
pairs <- list(
  # the two claims, first
  c("llr_meas", "apache2_aps"),
  c("llr_sum",  "sofa"),
  # context
  c("llr_sum", "apache2_aps"),
  c("llr_sum", "apache2_total"),
  c("llr_plus_demo", "apache2_total"),
  c("llr_sum", "aps_native"),
  # our own internal contrast, for free: what the intervention terms add
  c("llr_sum", "llr_meas"),
  # validation
  c("apache2_aps", "aps_native"),
  c("sofa", "sofa_native"))
pairs <- Filter(function(p) all(p %in% names(res)), pairs)

dl <- do.call(rbind, lapply(pairs, function(p) {
  cbind(a = p[1], b = p[2], delong_test(cells[[p[1]]], cells[[p[2]]], yk))
}))
save_table(run, dl, "delong_auroc")

pb <- do.call(rbind, lapply(pairs, function(p) {
  paired_boot_diff(cells[[p[1]]], cells[[p[2]]], yk, metric = .auprc,
                   n_boot = nboot, seed = cfg$seed,
                   label = paste(p[1], "-", p[2]))
}))
save_table(run, pb, "paired_boot_auprc")

# --- report -----------------------------------------------------------------
cat("\n=== cohort for this arm ===\n\n")
cat(sprintf("  %d stays scored, %d dropped for thin APACHE coverage (< %s of 12 variables)\n",
            sum(keep), n_drop, if (all_stays) "restriction off" else min_vars))
cat(sprintf("  event rate %.4f | arf_doubling = %s | gcs_source = %s\n",
            p_bar, arf, gcs_src))
cat("  Every cell is scored on THESE rows. The llr_sum AUROC here is therefore\n")
cat("  not the headline from tests/metrics.R, which uses all scored stays.\n")

cat("\n=== discrimination: the cells ===\n\n")
print(summ[, c("label", "n", "n_events", "auroc", "auroc_lo", "auroc_hi",
               "auprc", "auprc_lo", "auprc_hi", "auprc_lift")], row.names = FALSE)
cat("\n  These marginal intervals OVERLAP by construction: the cells are scored\n")
cat("  on identical rows and are strongly correlated. Read the paired tests.\n")

cat("\n=== paired AUROC differences (DeLong) ===\n\n")
print(dl[, c("a", "b", "auroc_1", "auroc_2", "delta", "ci_lo", "ci_hi", "p_value")],
      row.names = FALSE)

cat("\n=== paired AUPRC differences (bootstrap on identical resamples) ===\n\n")
print(pb, row.names = FALSE)

cat("\n=== calibration slope (1.000 = the score is in log-odds units) ===\n\n")
cal <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$cal)))
print(cal, row.names = FALSE)
cat("\n  The APACHE cells were recalibrated out-of-fold, so a slope near 1 is\n")
cat("  expected and is NOT evidence about APACHE. It is the recalibration\n")
cat("  working.\n")

cat("\n=== what the recalibration cost in discrimination ===\n\n")
print(recal, row.names = FALSE)
cat("\n  The mapping is monotone WITHIN a fold, but each fold gets its own, so\n")
cat("  stays in different folds can swap and AUROC can move a little. Expect\n")
cat("  the fourth decimal. If delta is ever large against the effect being\n")
cat("  reported, quote auroc_raw and say the recalibration was cosmetic.\n")

cat("\n=== monotonicity of the 20-tile curve ===\n\n")
mono <- do.call(rbind, lapply(names(res), function(k) cbind(label = k, res[[k]]$mono)))
print(mono, row.names = FALSE)

if (!is.null(sofa_agree)) {
  cat("\n=== GATE: our recomputed SOFA against the derived concept, per organ ===\n\n")
  print(sofa_agree, row.names = FALSE)
  cat(sprintf("\n  The total is over %d organ(s): %s\n", length(sofa_organs),
              paste(sofa_organs, collapse = ", ")))
  ex <- setdiff(SOFA_ALL_ORGANS, sofa_organs)
  if (length(ex)) {
    cat(sprintf("  EXCLUDED but still reported above: %s. `frac_input` is why.\n",
                paste(ex, collapse = ", ")))
    cat("  For those rows read `exact_computable`, not `exact_agree`: it is the\n")
    cat("  agreement where our input actually EXISTS, and it separates a coverage\n")
    cat("  failure from a scoring one. See config/sofa.organs.\n")
  }
  cat("\n  Cardiovascular IS EXPECTED TO DISAGREE: dx_nee_peak is a norepinephrine\n")
  cat("  equivalent and collapses SOFA's agent-specific tiers, so a patient on\n")
  cat("  dopamine at 5 lands at 3 here where the original gives 2. Any OTHER\n")
  cat("  organ with low exact_agree means a band, a unit or a missing-value rule\n")
  cat("  is wrong -- and eICU, which ships no SOFA, would inherit that error with\n")
  cat("  nothing to catch it. Fix before reporting either site.\n")
} else {
  cat("\n=== no native SOFA at this site: the recomputation is UNVALIDATED here ===\n\n")
  cat("  Expected at eICU. It means this run cannot check its own SOFA. Do not\n")
  cat("  report it unless the MIMIC validation table above passed first.\n")
}

cat("\n=== SOFA by organ: where the score comes from ===\n\n")
print(sofa_pts[order(-sofa_pts$mean_score), ], row.names = FALSE)

cat("\n=== DOMAIN LEVEL: each D_k against its SOFA organ ===\n\n")
print(dom_tab, row.names = FALSE)
cat("\n  D_k is the EQUAL-WEIGHT sum of the L's in the domain; layer 2's weights\n")
cat("  do not exist yet. auroc_D and auroc_organ say which discriminates better;\n")
cat("  `spearman` says whether they ORDER THE SAME PATIENTS as sick, which is a\n")
cat("  different question and the more interesting one. Two constructions can\n")
cat("  discriminate equally well while disagreeing about who is ill.\n")

cat("\n=== the five domains SOFA has no organ for ===\n\n")
print(uncov, row.names = FALSE)
cat("\n  domains.csv is disjoint and complete over the 19 signals, so these are\n")
cat("  exactly the evidence SOFA cannot represent at all. Report them beside the\n")
cat("  six-organ table: 'our domains beat SOFA's organs' is a much weaker claim\n")
cat("  than that plus 'and there are five more SOFA has no way to express'.\n")

cat("\n=== APACHE II input coverage (this is a caveat, not a diagnostic) ===\n\n")
print(cov_tab, row.names = FALSE)
cat("\n  A variable with low coverage contributes ZERO points on the stays that\n")
cat("  lack it, which makes those stays look healthy rather than unmeasured.\n")
cat("  Every such variable weakens the baseline and therefore flatters llr_sum.\n")

cat("\n=== where APACHE II's points actually come from ===\n\n")
print(pts[order(-pts$mean_points), ], row.names = FALSE)
cat("\n  Compare the ordering against signal_auroc from tests/metrics.R. Two\n")
cat("  hand-built and machine-built constructions agreeing on which signals\n")
cat("  carry the evidence is a result; disagreeing is a bigger one.\n")

if (!is.null(agree)) {
  cat("\n=== sanity: recomputed APACHE II total vs the site's native APS ===\n\n")
  print(agree, row.names = FALSE)
  cat("\n  These are DIFFERENT scores (APACHE II here, APS III natively), so they\n")
  cat("  must correlate strongly WITHOUT agreeing. Spearman below about 0.7\n")
  cat("  means the point table, a unit assumption or a sentinel guard is wrong,\n")
  cat("  and the arm must not be reported until it is found. This is audit A5\n")
  cat("  of sql/mimiciv/v2_06_severity_mimiciv.sql.\n")
}

g <- function(a, b, f) summ[[f]][summ$label == a] - summ[[f]][summ$label == b]
cat("\n=== the readings ===\n\n")
if (all(c("llr_meas", "apache2_aps") %in% summ$label)) {
  cat(sprintf("  llr_meas vs apache2_aps     dAUROC %+.4f  dAUPRC %+.4f   THE APACHE CLAIM (physiology only)\n",
              g("llr_meas", "apache2_aps", "auroc"), g("llr_meas", "apache2_aps", "auprc")))
}
if (all(c("llr_sum", "sofa") %in% summ$label)) {
  cat(sprintf("  llr_sum  vs sofa            dAUROC %+.4f  dAUPRC %+.4f   THE SOFA CLAIM (both bundle intervention)\n",
              g("llr_sum", "sofa", "auroc"), g("llr_sum", "sofa", "auprc")))
}
if (all(c("llr_sum", "llr_meas") %in% summ$label)) {
  cat(sprintf("  llr_sum  vs llr_meas        dAUROC %+.4f  dAUPRC %+.4f   what our intervention terms add\n",
              g("llr_sum", "llr_meas", "auroc"), g("llr_sum", "llr_meas", "auprc")))
}
if (all(c("llr_sum", "apache2_aps") %in% summ$label)) {
  cat(sprintf("  llr_sum  vs apache2_aps     dAUROC %+.4f  dAUPRC %+.4f   feature construction, full aggregate\n",
              g("llr_sum", "apache2_aps", "auroc"), g("llr_sum", "apache2_aps", "auprc")))
}
if (all(c("llr_sum", "apache2_total") %in% summ$label)) {
  cat(sprintf("  llr_sum vs apache2_total    dAUROC %+.4f  dAUPRC %+.4f   external comparability\n",
              g("llr_sum", "apache2_total", "auroc"), g("llr_sum", "apache2_total", "auprc")))
}
if (all(c("llr_plus_demo", "apache2_total") %in% summ$label)) {
  cat(sprintf("  llr_plus_demo vs total      dAUROC %+.4f  dAUPRC %+.4f   what age+comorbidity buy\n",
              g("llr_plus_demo", "apache2_total", "auroc"),
              g("llr_plus_demo", "apache2_total", "auprc")))
}
if (all(c("llr_sum", "aps_native") %in% summ$label)) {
  cat(sprintf("  llr_sum vs aps_native       dAUROC %+.4f  dAUPRC %+.4f   vs the canonical implementation\n",
              g("llr_sum", "aps_native", "auroc"), g("llr_sum", "aps_native", "auprc")))
}

a_of <- function(k) if (k %in% summ$label) summ$auroc[summ$label == k] else NA_real_
finalize_run(run, extra = list(
  source_run         = basename(src_dir),
  n_stays_scored     = sum(keep),
  n_stays_dropped    = n_drop,
  min_vars_present   = if (all_stays) NA_integer_ else min_vars,
  min_organs_present = if (all_stays) NA_integer_ else min_organs,
  arf_doubling       = arf,
  gcs_source         = gcs_src,
  resp_support       = resp_sup,
  sofa_validated     = !is.null(sofa_agree),
  auroc_llr_sum       = a_of("llr_sum"),
  auroc_llr_meas      = a_of("llr_meas"),
  auroc_apache2_aps   = a_of("apache2_aps"),
  auroc_apache2_total = a_of("apache2_total"),
  auroc_aps_native    = a_of("aps_native"),
  auroc_sofa          = a_of("sofa"),
  auroc_sofa_native   = a_of("sofa_native")))
cat(sprintf("\n  run directory: %s\n", run$path))
