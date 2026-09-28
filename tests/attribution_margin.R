# tests/attribution_margin.R -------------------------------------------------
# IS THE ATTRIBUTION BRITTLE, OR IS THE METRIC?
#
# tests/coupling_attribution.R found that two defensible specifications name a
# different top-contributing DOMAIN for about 19% of patients, and a different
# top SIGNAL for about 24%. Even `full` against `full_ti_trend` -- two models
# that differ by 0.0014 AUROC -- disagree about the top domain for 8% of
# patients. That looks alarming and it may not be.
#
# THE QUESTION THIS SCRIPT ANSWERS. "Top-1 agreement" is a DISCRETE statistic
# computed on what is often a near-tie. If a patient's largest and second
# largest domain contributions are 0.81 and 0.79 nats, then an arbitrarily small
# displacement reorders them, and the two models have not disagreed about
# anything real -- they have both declined to separate two domains that are not
# separated. If instead the disagreement persists where the leader is well
# clear, the attribution genuinely is brittle.
#
# So every agreement statistic here is reported CONDITIONAL ON THE MARGIN: the
# gap in nats between the top contributor and the runner-up, under the reference
# arm. Two readings follow and they are different claims:
#
#   the margin DISTRIBUTION   how often is there a defensible winner at all?
#                             This is a property of the design, not of any
#                             comparison, and it bounds what ANY model can claim.
#   agreement BY margin bin   given that there is a clear winner, do two
#                             specifications agree on it?
#
# If agreement rises to near 1 in the wide-margin bins, the honest reporting
# rule is "name a top contributor only when the margin clears X nats", and that
# is a publishable rule rather than a caveat.
#
# ALSO HERE: where flips live relative to zero, and what age and sex are worth.
# The second is not a coupling question but it is the same kind of question --
# what is the score missing -- and it costs nothing to answer alongside.
#
# FITS NOTHING. It reads the cached ladder written by
# tests/coupling_attribution.R, which holds all eight L matrices.
#
# Aggregates only (hard rule 1): counts, quantiles and agreement rates. The L
# matrices and the demographics are row-level and are never printed.
#
#   Rscript tests/attribution_margin.R out/runs/coupattr_20260905T184528
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

args    <- commandArgs(trailingOnly = TRUE)
ladder_d <- if (length(args) >= 1L && nzchar(args[1])) args[1] else
              latest_run("coupattr", require_complete = FALSE)
lad <- readRDS(file.path(ladder_d, "tables", "l_oof_ladder.rds"))
A   <- lad$arms

cfg     <- tar_read(cfg)
tabs    <- tar_read(tabs)
tr      <- tar_read(train_ids)
y       <- as.integer(tar_read(y_train))
domains <- load_domains("config/domains.csv")
sigs    <- unlist(cfg$signals)
stopifnot(identical(as.character(lad$stay_id), as.character(tr)))

run <- new_run("attrmargin", cfg, note = sprintf(
  "margin-conditioned attribution agreement, ladder from %s, fits nothing",
  basename(ladder_d)))

# Bin edges in nats, declared before anything is computed.
BINS <- c(0, 0.05, 0.10, 0.25, 0.50, 1.00, Inf)
BINL <- c("<0.05", "0.05-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.00", ">1.00")

CMP <- list(
  c("cond",          "cond_ti_all"),
  c("cond",          "cond_ti_trend"),
  c("cond_ti_trend", "cond_ti_all"),
  c("meas",          "cond"),
  c("meas",          "cond_ti_all"),
  c("full",          "full_ti_all"),
  c("full",          "full_ti_trend"))

# --- helpers ----------------------------------------------------------------

#' Top-1 index and the margin to the runner-up, per row.
#'
#' The margin is in NATS on the same scale as the L's, so it is comparable
#' across patients only to the extent their total evidence is. `rel` divides it
#' by the row's total absolute evidence, which is the scale-free version.
top_margin <- function(M) {
  ab <- abs(M)
  i1 <- max.col(ab, ties.method = "first")
  v1 <- ab[cbind(seq_len(nrow(ab)), i1)]
  ab2 <- ab; ab2[cbind(seq_len(nrow(ab)), i1)] <- -Inf
  v2 <- ab2[cbind(seq_len(nrow(ab2)), max.col(ab2, ties.method = "first"))]
  tot <- rowSums(ab)
  list(top = i1, margin = v1 - v2, rel = (v1 - v2) / pmax(tot, .Machine$double.eps),
       lead = v1, total = tot)
}

dm  <- domains$domain[match(sigs, domains$signal)]
dnm <- sort(unique(dm))
to_dom <- function(M) {
  D <- matrix(0, nrow(M), length(dnm), dimnames = list(rownames(M), dnm))
  for (k in dnm) {
    j <- which(dm == k)
    D[, k] <- if (length(j) == 1L) M[, j] else rowSums(M[, j, drop = FALSE])
  }
  D
}

#' Agreement on the top contributor, within bins of the reference arm's margin.
agree_by_margin <- function(Ma, Mb, label_a, label_b, level) {
  ta <- top_margin(Ma); tb <- top_margin(Mb)
  ok <- ta$top == tb$top
  b  <- cut(ta$margin, breaks = BINS, labels = BINL, right = FALSE)
  do.call(rbind, lapply(seq_along(BINL), function(i) {
    s <- !is.na(b) & b == BINL[i]
    data.frame(level = level, arm_a = label_a, arm_b = label_b,
               margin_bin = BINL[i], n_patients = sum(s),
               share_of_cohort = round(mean(s), 5),
               top1_agree = round(if (any(s)) mean(ok[s]) else NA_real_, 5),
               stringsAsFactors = FALSE)
  }))
}

# --- A. how often is there a defensible winner at all? ----------------------
#
# A property of the DESIGN, before any comparison. If most patients have no
# clear leader, then "the top domain for this patient" is not a claim any model
# can support, and the specification disagreement is a symptom rather than the
# disease.
cat("\n=== A. the margin distribution: is there a winner to disagree about? ===\n\n")
cat("    margin = |L| of the top contributor minus |L| of the runner-up, nats.\n")
cat("    rel    = the same divided by the patient's total absolute evidence.\n\n")
cat(sprintf("%-16s %-8s %9s %9s %9s %9s %10s\n", "arm", "level",
            "med", "p25", "p75", "med rel", "frac<0.10"))
MD <- list()
for (a in names(A)) {
  for (lv in c("signal", "domain")) {
    M <- if (lv == "signal") A[[a]] else to_dom(A[[a]])
    t <- top_margin(M)
    MD[[length(MD) + 1L]] <- data.frame(arm = a, level = lv,
      margin_median = round(stats::median(t$margin), 5),
      margin_p25 = round(unname(stats::quantile(t$margin, 0.25)), 5),
      margin_p75 = round(unname(stats::quantile(t$margin, 0.75)), 5),
      rel_median = round(stats::median(t$rel), 5),
      frac_below_0.10 = round(mean(t$margin < 0.10), 5),
      frac_below_0.25 = round(mean(t$margin < 0.25), 5),
      stringsAsFactors = FALSE)
    z <- MD[[length(MD)]]
    cat(sprintf("%-16s %-8s %9.4f %9.4f %9.4f %9.4f %10.4f\n", a, lv,
                z$margin_median, z$margin_p25, z$margin_p75, z$rel_median,
                z$frac_below_0.10))
  }
}
save_table(run, do.call(rbind, MD), "margin_distribution", subdir = "diagnostics")

# --- B. agreement conditional on the margin ---------------------------------
cat("\n=== B. top-1 agreement WITHIN bins of the reference arm's margin ===\n\n")
cat("    If agreement approaches 1 in the wide bins, the attribution is sound\n")
cat("    and the pooled number was measuring near-ties. If it stays low in the\n")
cat("    wide bins, the attribution is genuinely specification-dependent.\n")
G <- list()
for (lv in c("domain", "signal")) {
  cat(sprintf("\n  --- %s level ---\n\n", lv))
  cat(sprintf("  %-15s %-15s", "arm a", "arm b"))
  for (bl in BINL) cat(sprintf(" %10s", bl))
  cat("\n")
  for (p in CMP) {
    Ma <- if (lv == "signal") A[[p[1]]] else to_dom(A[[p[1]]])
    Mb <- if (lv == "signal") A[[p[2]]] else to_dom(A[[p[2]]])
    g <- agree_by_margin(Ma, Mb, p[1], p[2], lv)
    G[[length(G) + 1L]] <- g
    cat(sprintf("  %-15s %-15s", p[1], p[2]))
    for (i in seq_along(BINL)) cat(sprintf(" %10.4f", g$top1_agree[i]))
    cat("\n")
  }
  cat(sprintf("  %-15s %-15s", "", "share of cohort"))
  gg <- G[[length(G)]]
  for (i in seq_along(BINL)) cat(sprintf(" %10.4f", gg$share_of_cohort[i]))
  cat("\n")
}
save_table(run, do.call(rbind, G), "agreement_by_margin", subdir = "diagnostics")

# --- C. a reporting rule ----------------------------------------------------
#
# The constructive form of section B. For each threshold, how much of the cohort
# survives it and what is the agreement among those who do.
cat("\n=== C. a candidate reporting rule: name a top domain only above a margin ===\n\n")
cat(sprintf("%-15s %-15s %8s %10s %10s %10s\n", "arm a", "arm b",
            "thresh", "kept", "agree", "agree_all"))
RR <- list()
for (p in list(c("cond", "cond_ti_all"), c("meas", "cond"), c("meas", "cond_ti_all"))) {
  Da <- to_dom(A[[p[1]]]); Db <- to_dom(A[[p[2]]])
  ta <- top_margin(Da); tb <- top_margin(Db)
  ok <- ta$top == tb$top
  base <- mean(ok)
  for (th in c(0, 0.10, 0.25, 0.50, 1.00)) {
    s <- ta$margin >= th
    RR[[length(RR) + 1L]] <- data.frame(arm_a = p[1], arm_b = p[2], level = "domain",
      threshold = th, frac_kept = round(mean(s), 5),
      top1_agree_kept = round(mean(ok[s]), 5), top1_agree_all = round(base, 5),
      stringsAsFactors = FALSE)
    cat(sprintf("%-15s %-15s %8.2f %10.4f %10.4f %10.4f\n", p[1], p[2], th,
                mean(s), mean(ok[s]), base))
  }
}
save_table(run, do.call(rbind, RR), "reporting_rule", subdir = "diagnostics")

# --- D. where do the sign flips live relative to zero? ----------------------
#
# A flip of a cell that was already near zero is not a disagreement about
# anything. This asks what fraction of flips are of that kind, by reporting the
# magnitude the LARGER of the two arms assigned.
cat("\n=== D. flipped cells: how much evidence was actually at stake? ===\n\n")
meas_ok <- measured_matrix(tabs, cfg, tr)
FD <- list()
cat(sprintf("%-15s %-15s %10s %10s %10s %10s %10s\n", "arm a", "arm b",
            "n_flip", "med mag", "p90 mag", "frac<0.05", "frac>0.25"))
for (p in CMP) {
  a <- A[[p[1]]]; b <- A[[p[2]]]
  fl <- (sign(a) != sign(b)) & meas_ok
  mg <- pmax(abs(a), abs(b))[fl]
  FD[[length(FD) + 1L]] <- data.frame(arm_a = p[1], arm_b = p[2],
    n_flips = sum(fl), med_magnitude = round(stats::median(mg), 5),
    p90_magnitude = round(unname(stats::quantile(mg, 0.90)), 5),
    frac_below_0.05 = round(mean(mg < 0.05), 5),
    frac_above_0.25 = round(mean(mg > 0.25), 5), stringsAsFactors = FALSE)
  z <- FD[[length(FD)]]
  cat(sprintf("%-15s %-15s %10d %10.4f %10.4f %10.4f %10.4f\n", p[1], p[2],
              z$n_flips, z$med_magnitude, z$p90_magnitude, z$frac_below_0.05,
              z$frac_above_0.25))
}
save_table(run, do.call(rbind, FD), "flip_magnitude", subdir = "diagnostics")

# --- E. what are age and sex worth? -----------------------------------------
#
# NOT a coupling question, and asked here only because it costs nothing. The
# design excludes both from layer 1, and `tests/metrics_severity.R` already
# measures age plus chronic health plus admission class as a SCORE-LEVEL block
# (`llr_plus_demo`). What is NOT measured anywhere is SEX, which appears in no
# arm of this project at all.
#
# Fitted as smooths, because section 6 of the coupling document established that
# a linear increment test systematically under-detects here.
cat("\n=== E. what age and sex are worth on top of the LLR ===\n\n")
co  <- tabs$cohort
m   <- match(tr, co$stay_id)
age <- as.numeric(co$age[m]); sex <- droplevels(co$gender[m])
cat(sprintf("  cohort: %d stays, %d with age, sex levels {%s}, share F = %.4f\n\n",
            length(tr), sum(!is.na(age)), paste(levels(sex), collapse = ","),
            mean(sex == "F", na.rm = TRUE)))
E <- list()
for (arm in c("meas", "cond", "cond_ti_all", "full")) {
  s <- rowSums(A[[arm]])
  d <- data.frame(y = y, s = s, age = age, sex = sex)
  d <- d[stats::complete.cases(d), ]
  g0 <- mgcv::gam(y ~ s(s), family = binomial(), method = "REML", data = d)
  g1 <- mgcv::gam(y ~ s(s) + s(age), family = binomial(), method = "REML", data = d)
  g2 <- mgcv::gam(y ~ s(s) + sex, family = binomial(), method = "REML", data = d)
  g3 <- mgcv::gam(y ~ s(s) + s(age) + sex, family = binomial(), method = "REML", data = d)
  au <- function(g) .auroc(stats::predict(g, type = "link"), d$y)
  E[[arm]] <- data.frame(arm = arm, n = nrow(d),
    auroc = round(au(g0), 5), auroc_age = round(au(g1), 5),
    auroc_sex = round(au(g2), 5), auroc_age_sex = round(au(g3), 5),
    d_auroc_age = round(au(g1) - au(g0), 5),
    d_auroc_sex = round(au(g2) - au(g0), 5),
    d_dev_age = round(summary(g1)$dev.expl - summary(g0)$dev.expl, 5),
    d_dev_sex = round(summary(g2)$dev.expl - summary(g0)$dev.expl, 5),
    d_aic_age = round(AIC(g1) - AIC(g0), 1), d_aic_sex = round(AIC(g2) - AIC(g0), 1),
    stringsAsFactors = FALSE)
  z <- E[[arm]]
  cat(sprintf("  %-12s AUROC %.5f  +age %+.5f (dAIC %+.0f)  +sex %+.5f (dAIC %+.0f)  +both %.5f\n",
              arm, z$auroc, z$d_auroc_age, z$d_aic_age, z$d_auroc_sex, z$d_aic_sex,
              z$auroc_age_sex))
}
save_table(run, do.call(rbind, E), "demographics_increment", subdir = "diagnostics")

cat("\n  Compare against tests/metrics_severity.R's `llr_plus_demo`, which adds\n")
cat("  age + chronic health + admission class to llr_sum and measured\n")
cat("  0.8457 -> 0.8601 AUROC on 2026-08-31. Sex is in no arm of that cell.\n")

finalize_run(run, extra = list(ladder = basename(ladder_d),
                               bins = paste(BINS, collapse = ",")))
cat(sprintf("\nwritten: %s\n", run$path))
