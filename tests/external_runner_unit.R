# tests/external_runner_unit.R ------------------------------------------------
# SYNTHETIC UNIT TESTS FOR THE EXTERNAL RUNNER'S REPAIRS. Reads no data, no
# bundle, no targets store; every input is built in this file. Runs in seconds.
#
# Written 2026-09-09 against `docs/external_runner_review_20260909.md`. Each
# block reproduces one of that review's synthetic checks and asserts the
# repaired behaviour:
#
#   E1  a failing stage stamps the manifest `failed` with stage and reason
#   E2  `.capture_check()` distinguishes pass / warning / error
#   E3  the orchestration hash covers config/; input fingerprints are MD5
#   E4  group_metrics reports three denominators and handles empty cases
#   E5  bin and calibration tables label their interval assumptions;
#       contrasts carry the clustered p and the resampling unit
#   E6  incompatible runner settings are refused before loading; contrasts
#       are derived from the selected arms
#   E7  the schema signature compares what the live comparison compared, and
#       carries no row value
#   E8  above-chance retention and the declared near-chance policy
#
# Runner-local functions (`.capture_check`, `validate_external_config`) are
# taken from the PARSED runner, never by sourcing it, so nothing is loaded.
#
#   Rscript tests/external_runner_unit.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(yaml); library(digest)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

FAILS <- 0L
check <- function(what, ok) {
  ok <- isTRUE(ok)
  cat(sprintf("  %-72s %s\n", what, if (ok) "PASS" else "*** FAIL ***"))
  if (!ok) FAILS <<- FAILS + 1L
  invisible(ok)
}
expect_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
n_warnings <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") })
  n
}

# Runner-local functions, lifted from the parsed script.
runner_fn <- function(name) {
  ex <- parse("run/external.R", keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(as.character(e[[2]]), name)) return(eval(e[[3]], envir = globalenv()))
  }
  stop("run/external.R defines no `", name, "`")
}
.capture_check <- runner_fn(".capture_check")
validate_external_config <- runner_fn("validate_external_config")

cat("\n=== external runner: synthetic unit tests ===\n\n")

# --- E2: three outcomes -------------------------------------------------------
cat("E2  site-check capture\n")
r <- .capture_check(function() 1)
check("clean check: outcome pass, ok TRUE, value kept", r$outcome == "pass" && r$ok && r$value == 1)
r <- .capture_check(function() { warning("declaration differs"); 2 })
check("warning: outcome warning, ok FALSE, value kept, message recorded",
      r$outcome == "warning" && !r$ok && r$value == 2 && grepl("declaration", r$msg))
r <- .capture_check(function() stop("cannot evaluate"))
check("error: outcome error, ok FALSE, value NULL, message recorded",
      r$outcome == "error" && !r$ok && is.null(r$value) && grepl("cannot", r$msg))
r <- .capture_check(function() { warning("first"); stop("then") })
check("warning then error: outcome error, both messages recorded",
      r$outcome == "error" && grepl("first", r$msg) && grepl("then", r$msg))

# --- E6: config validation before load --------------------------------------
cat("E6  runner configuration validation\n")
ec <- yaml::read_yaml("config/external.yml")
rc <- validate_external_config(ec)
check("shipped config validates; full ladder reportable",
      length(rc$ladder) == length(LADDER_CONTRASTS) && !length(rc$dropped_contrasts))
ec2 <- ec; ec2$frozen_bins <- FALSE; ec2$self_bins <- FALSE
check("both binnings off is refused", expect_error(validate_external_config(ec2)))
ec3 <- ec; ec3$arms <- c("llr_sum", "not_an_arm")
check("unknown arm is refused", expect_error(validate_external_config(ec3)))
ec4 <- ec; ec4$arms <- c("llr_sum", "xgb_raw", "llr_meas")
rc4 <- validate_external_config(ec4)
want4 <- list(c("llr_sum", "xgb_raw"), c("llr_sum", "llr_meas"))
check("reduced arm set: contrasts derived from the arms, dropped ones named",
      identical(unname(rc4$ladder), want4) &&
        "xgb_feat - xgb_raw" %in% rc4$dropped_contrasts)
ec5 <- ec; ec5$hospital$unmatched_policy <- NULL
check("hospital enabled without an unmatched policy is refused",
      expect_error(validate_external_config(ec5)))
check("... unless the hospital analysis is skipped from the CLI",
      !expect_error(validate_external_config(ec5, no_hospital = TRUE)))
ec6 <- ec; ec6$hospital$unmatched_policy <- "ignore"
check("an undeclared unmatched policy value is refused", expect_error(validate_external_config(ec6)))
ec7 <- ec; ec7$hospital$min_stays <- 0
check("a non-positive hospital floor is refused", expect_error(validate_external_config(ec7)))
ec8 <- ec; ec8$checks$design_checks_strict <- "yes"
check("a non-logical check switch is refused", expect_error(validate_external_config(ec8)))
check("contrasts_available() keeps only pairs whose arms were scored",
      length(contrasts_available(LADDER_CONTRASTS, c("llr_sum", "llr_cond"))) == 1L)

# --- E4: denominators and empty cases ----------------------------------------
cat("E4  group_metrics denominators\n")
set.seed(4)
n <- 400
y <- rbinom(n, 1, 0.3)
s <- y * 1.2 + rnorm(n)
g <- rep(c("A", "B", "C"), c(200, 150, 50))
g[1:40] <- NA                                     # 40 scored stays with no hospital
gm <- group_metrics(s, y, g, label = "t", min_n = 100L, min_events = 10L)
sm <- gm$summary
check("per-group columns unchanged for tests/coupling_strata.R",
      identical(names(gm$per_group),
                c("group", "n", "deaths", "event_rate", "reported", "auroc", "auprc", "auprc_lift")))
check("n_scored / n_matched / n_unmatched are 400 / 360 / 40",
      sm$n_scored == 400 && sm$n_matched == 360 && sm$n_unmatched == 40)
kept <- sum(gm$per_group$n[gm$per_group$reported])
check("frac_kept_of_scored and frac_kept_of_matched use their own denominators",
      isTRUE(all.equal(sm$frac_kept_of_scored, round(kept / 400, 4))) &&
        isTRUE(all.equal(sm$frac_kept_of_matched, round(kept / 360, 4))) &&
        sm$frac_kept_of_scored < sm$frac_kept_of_matched)
ok <- !is.na(g)
elig <- ok & g %in% gm$per_group$group[gm$per_group$reported]
check("pooled_auroc_matched is the AUROC over matched rows",
      isTRUE(all.equal(sm$pooled_auroc_matched, round(.auroc(s[ok], y[ok]), 5))))
check("pooled_auroc_eligible is the AUROC over rows of reported groups",
      isTRUE(all.equal(sm$pooled_auroc_eligible, round(.auroc(s[elig], y[elig]), 5))))
check("pooled_minus_median is eligible-pooled minus the eligible median",
      abs(sm$pooled_minus_median - (sm$pooled_auroc_eligible - sm$auroc_median)) < 2e-5)
check("status ok when at least one group is reported", sm$status == "ok")

# The review's four-row case: two unmatched, the rest one eligible group.
gm2 <- group_metrics(c(0.1, 0.9, 0.4, 0.6), c(0L, 1L, 0L, 1L), c("H", "H", NA, NA),
                     min_n = 2L, min_events = 1L)
check("review case: frac_kept_of_matched 1 but frac_kept_of_scored 0.5",
      gm2$summary$frac_kept_of_matched == 1 && gm2$summary$frac_kept_of_scored == 0.5 &&
        gm2$summary$n_unmatched == 2)

nw <- n_warnings(gm3 <- group_metrics(s, y, g, min_n = 10000L, min_events = 10L))
check("no eligible group: status set, extrema NA, no warning",
      gm3$summary$status == "no_eligible_group" && is.na(gm3$summary$auroc_min) &&
        is.na(gm3$summary$auroc_max) && is.na(gm3$summary$pooled_minus_median) &&
        gm3$summary$n_reported == 0 && nw == 0L)
check("no eligible group: matched pooled AUROC still reported",
      is.finite(gm3$summary$pooled_auroc_matched))
nw <- n_warnings(gm4 <- try(group_metrics(s, y, rep(NA_character_, n)), silent = TRUE))
check("all rows unmatched: no error, status no_matched_rows, zero groups",
      !inherits(gm4, "try-error") && gm4$summary$status == "no_matched_rows" &&
        gm4$summary$n_groups == 0 && gm4$summary$n_matched == 0 && nw == 0L)

# --- E8: retention and the pooled-minus-median contrast ----------------------
cat("E8  retention ratios and the heterogeneity contrast\n")
check("chance-level external against 0.8 training retains 0 above chance (ratio 0.625)",
      retention_above_chance(0.5, 0.8) == 0 && round(0.5 / 0.8, 4) == 0.625)
check("half the excess retained reads 0.5", retention_above_chance(0.65, 0.8) == 0.5)
check("training within 0.05 of chance gives NA by policy",
      is.na(retention_above_chance(0.6, 0.53)) && !is.na(retention_above_chance(0.6, 0.55)))
check("non-finite inputs give NA, vectorised",
      identical(is.na(retention_above_chance(c(NA, 0.7, 0.7), c(0.8, NA, 0.8))), c(TRUE, TRUE, FALSE)))

# Three hospitals, identical prevalence (50%) and identical marginal binary
# score distribution (half the stays score 1), unequal sizes and discrimination.
# Case mix cannot explain any pooled-versus-median gap here.
yh <- c(1, 0,  1, 0,  1, 1, 1, 1, 0, 0, 0, 0)
sh <- c(0, 1,  0, 1,  1, 1, 1, 0, 0, 0, 0, 1)
hh <- c("h1", "h1", "h2", "h2", rep("h3", 8))
gm5 <- group_metrics(sh, yh, hh, min_n = 2L, min_events = 1L)
check("identical-prevalence hospitals: every event rate is 0.5",
      all(gm5$per_group$event_rate == 0.5))
check("... yet the pooled-minus-median contrast is non-zero (not a case-mix decomposition)",
      gm5$summary$auroc_median == 0 && gm5$summary$pooled_minus_median > 0 &&
        isTRUE(all.equal(gm5$summary$pooled_auroc_eligible, round(.auroc(sh, yh), 5))))

# --- E5: interval labels and the clustered contrast columns ------------------
cat("E5  interval assumptions are labelled\n")
rb <- risk_bins(s, y, n_bins = 5L)
check("risk_bins carries ci_method = wilson_stay_iid", all(rb$ci_method == "wilson_stay_iid"))
ca <- llr_calibration(s, y, 0.3)
check("llr_calibration carries ci_method and the p_bar it was given",
      ca$ci_method == "glm_profile_stay_iid" && ca$p_bar_reference == 0.3)
grp <- as.character(rep(seq_len(n / 2), each = 2))
ct <- arm_contrasts(list(a = s, b = s + rnorm(n, sd = 0.3)), y,
                    list(c("a", "b")), n_boot = 20L, seed = 1L, group = grp)
check("arm_contrasts carries auroc_p, delong_p and boot_unit = patient",
      all(c("auroc_p", "delong_p", "delong_se", "boot_unit") %in% names(ct)) &&
        ct$boot_unit == "patient")

# --- E7: schema signature -----------------------------------------------------
cat("E7  frozen schema signature\n")
mk <- function(site, a_class = as.integer, lv = c("x", "y"), with_ord = TRUE) {
  tb <- list(
    cohort = data.frame(site = factor(site), a = a_class(1:3),
                        f = factor(c("x", "y", "x"), levels = lv),
                        txt = c("SECRETROWVALUE1", "SECRETROWVALUE2", "SECRETROWVALUE3"),
                        stringsAsFactors = FALSE),
    signal_features = data.frame(site = factor(site), v = c(1.5, 2.5, 3.5)))
  if (with_ord) tb$ordering <- data.frame(site = factor(site), o = 1:2)
  tb
}
ta <- mk("mimic"); tb <- mk("eicu", with_ord = FALSE)
sig <- schema_signature(ta)
check("signature carries names, classes and non-site factor levels only",
      identical(names(sig$cohort), c("names", "classes", "levels")) &&
        identical(names(sig$cohort$levels), "f") && !"site" %in% names(sig$cohort$levels))
check("signature carries no row value (hard rule 1)",
      !any(grepl("SECRETROWVALUE", deparse(sig))) && !any(grepl("1.5", deparse(sig$signal_features), fixed = TRUE)))
se <- suppressMessages(compare_schema_signatures(sig, schema_signature(tb), strict = FALSE))
check("matching schemas: no FAIL; the one-sided table is a SKIP row",
      sum(se$status == "FAIL") == 0 && any(se$check == "9 tables" & se$status == "SKIP"))
se2 <- suppressMessages(compare_schema_signatures(sig, schema_signature(mk("eicu", a_class = as.numeric)), strict = FALSE))
check("a class change is a `9 types` FAIL",
      any(se2$check == "9 types" & se2$status == "FAIL" & grepl("a \\(integer vs numeric\\)", se2$detail)))
se3 <- suppressMessages(compare_schema_signatures(sig, schema_signature(mk("eicu", lv = c("x", "y", "z"))), strict = FALSE))
check("a level-set change is a `9 levels` FAIL",
      any(se3$check == "9 levels" & se3$status == "FAIL" & se3$detail == "f"))
check("strict mode stops on a FAIL",
      expect_error(suppressMessages(compare_schema_signatures(sig, schema_signature(mk("eicu", lv = c("x", "y", "z"))), strict = TRUE))))
live <- suppressMessages(validate_schema_equality(ta, mk("eicu"), strict = FALSE))
froz <- suppressMessages(compare_schema_signatures(sig, schema_signature(mk("eicu")), strict = FALSE))
check("the live comparison and the frozen comparison produce identical tables", identical(live, froz))
check("build_bundle() accepts a `schema` slot, NULL by default",
      "schema" %in% names(formals(build_bundle)) && is.null(formals(build_bundle)$schema))

# --- E1 / E3: manifest status and provenance ---------------------------------
cat("E1  failed-stage manifest\n")
root <- file.path(tempdir(), "llr_unit_runs")
run <- new_run("unittest", list(a = 1), root = root, note = "unit")
check("run_stage returns the stage's value on success", run_stage(run, "ok_stage", 41 + 1) == 42)
err <- try(run_stage(run, "boom_stage", stop("synthetic boom"), extra = list(k = 7L)), silent = TRUE)
m <- read_manifest(run$path)
check("a raising stage re-raises after stamping the manifest",
      inherits(err, "try-error") && grepl("synthetic boom", as.character(err)))
check("manifest status failed, with stage, reason and the extra content",
      identical(m$status, "failed") && identical(m$failed_stage, "boom_stage") &&
        grepl("synthetic boom", m$failure_reason) && identical(m$k, 7L))
finalize_run(run)
check("finalize_run still flips a run to complete", identical(read_manifest(run$path)$status, "complete"))
unlink(root, recursive = TRUE)

cat("E3  provenance\n")
oh <- .orchestration_hashes()
check("orchestration hashes cover config/ as well as run/ and _targets.R",
      all(c("external.yml", "config.yml", "pairing.csv", "external.R") %in% names(oh)))
tmp <- tempfile(fileext = ".bin"); writeBin(as.raw(1:64), tmp)
fp <- .input_fingerprint(list(present = tmp, absent = file.path(tempdir(), "no_such_file.parquet")))
check("input fingerprint: MD5 and size for a present file, NA and exists=FALSE for an absent one",
      fp$present$exists && fp$present$bytes == 64 &&
        identical(fp$present$md5, unname(as.character(tools::md5sum(tmp)))) &&
        !fp$absent$exists && is.na(fp$absent$md5))
unlink(tmp)

cat(sprintf("\n%s\n", if (FAILS) sprintf("*** %d FAILURE(S) ***", FAILS) else "all checks passed"))
if (FAILS) quit(status = 1L)
