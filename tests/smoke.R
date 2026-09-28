# tests/smoke.R --------------------------------------------------------------
# End-to-end check of everything built so far, stopping at the fit boundary.
# Nothing here fits a GAM; it answers "would the 155 fits start cleanly, and are
# the quantities feeding them scoped correctly".
#
#   Rscript tests/smoke.R
#
# Sections A-B are data-free and run in milliseconds. C onward touch data/ and
# take ~15 s, dominated by the parquet read.
#
# OUTPUT IS AGGREGATES ONLY — counts, rates, level sets, dimensions. No rows, no
# identifiers, no individual measurements (hard rule 1). Anything added here
# must keep that property.
# ----------------------------------------------------------------------------

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

FAIL <- 0L
ok <- function(label, cond, note = "") {
  cond <- isTRUE(cond)
  if (!cond) FAIL <<- FAIL + 1L
  # `note` IS COERCED TO LENGTH ONE, and that is a bug fix rather than defensive
  # coding. `sprintf()` is vectorised over its arguments and returns
  # `character(0)` if any of them is zero-length, so `cat()` printed NOTHING --
  # a failure was counted in `FAIL` and its line never appeared. FOUND
  # 2026-09-07, when `layer1_budget()`'s `aliased` column was split into
  # `alias_meas` and `alias_full` and this file's `sprintf("%d", bd$aliased)`
  # became `sprintf("%d", NULL)`. The run reported four failures and printed
  # two, which is the worst possible behaviour for a check script: it knew and
  # would not say.
  if (!length(note)) note <- "(no note: the value was NULL or zero-length)"
  cat(sprintf("  [%s] %-58s %s\n", if (cond) "ok" else "FAIL", label,
              paste(as.character(note)[1], collapse = "")))
}
sect <- function(x) cat(sprintf("\n== %s %s\n", x, strrep("-", max(0, 60 - nchar(x)))))

# --- A. Dirichlet-multinomial estimator, on synthetic data ------------------
# The only part of the pipeline whose correctness is not visible by inspection.
# Simulate from a known alpha and check it comes back. No project data involved,
# so a failure here is unambiguously an estimator bug.
sect("A. fit_dm_alpha recovers a known alpha")

rdirichlet <- function(n, a) {
  g <- matrix(stats::rgamma(n * length(a), rep(a, each = n)), nrow = n)
  g / rowSums(g)
}
sim_dm <- function(N, n_i, alpha) {
  p <- rdirichlet(N, alpha)
  t(vapply(seq_len(N), function(i) stats::rmultinom(1, n_i[i], p[i, ])[, 1], numeric(3)))
}

# Cases span the operating regime and then some. MEASURED 2026-08-25: the 19
# signals sit at alpha0 = 0.28-2.96 on N = 17,861-41,195 stays, so cases 1-3
# bracket it with an order of magnitude of headroom.
#
# Tolerances are set from the identification argument, not from observed error.
# alpha0 is identified only through overdispersion relative to a multinomial,
# and that vanishes as alpha0 grows past max(n_obs) = 24 — so alpha0's sampling
# variance climbs while the prior MEAN alpha/alpha0 stays sharp. Case 4 is
# deliberately in that flat-likelihood regime and is therefore asserted on the
# mean and only loosely on alpha0. Seeds are set per case so a pass is a fact
# about the estimator and not about the draw.
run_case <- function(seed, a, N = 20000L) {
  set.seed(seed)
  n_i <- sample(1:24, N, replace = TRUE)
  fit_dm_alpha(sim_dm(N, n_i, a), max_iter = 50000L)
}
for (tc in list(list(s = 11L, a = c(0.3, 0.8, 0.2), lab = "weak prior     (a0 = 1.3)"),
                list(s = 12L, a = c(2.0, 5.0, 1.0), lab = "operating range (a0 = 8)"),
                list(s = 13L, a = c(6.0, 10.0, 4.0), lab = "above range    (a0 = 20)"))) {
  fit <- run_case(tc$s, tc$a)
  rel <- max(abs(fit$alpha - tc$a) / tc$a)
  ok(tc$lab, fit$converged && rel < 0.10, sprintf("rel.err %.3f in %d iters", rel, fit$iter))
}
fit <- run_case(14L, c(20, 40, 10))
rel_mean <- max(abs(fit$alpha / fit$alpha0 - c(20, 40, 10) / 70) / (c(20, 40, 10) / 70))
ok("flat-likelihood regime (a0 = 70)", fit$converged && rel_mean < 0.05 &&
     fit$alpha0 > 35 && fit$alpha0 < 140,
   sprintf("prior-mean rel.err %.3f | a0 %.1f vs 70 | %d iters", rel_mean, fit$alpha0, fit$iter))

# A category that never occurs must be reported, not silently returned as a
# plausible small number, because it makes pi_hat constant and mgcv's error for
# that reads like a basis-dimension problem.
k0 <- sim_dm(5000, sample(1:24, 5000, TRUE), c(2, 5, 1)); k0[, 3] <- 0
ok("degenerate category flagged", any(fit_dm_alpha(k0)$degenerate))

# Posterior means are a probability vector, and shrinkage moves toward the prior
# mean by an amount set by n_obs. Both are load-bearing for pi_mid being
# determined by the other two.
a <- c(2, 5, 1)
p1 <- shrink_pi(matrix(c(1, 0, 0), 1), a)
p24 <- shrink_pi(matrix(c(24, 0, 0), 1), a)
ok("shrink_pi rows sum to 1", isTRUE(all.equal(sum(p1), 1)) && isTRUE(all.equal(sum(p24), 1)))
ok("shrinkage decreases with n_obs", p1[1, "pi_minus"] < p24[1, "pi_minus"],
   sprintf("n=1: %.3f   n=24: %.3f", p1[1, "pi_minus"], p24[1, "pi_minus"]))

# --- B. Design files and formulas -------------------------------------------
sect("B. config, pairing, formulas (data-free)")

cfg <- load_config("config/config.yml")
ok("config + pairing agree", TRUE, sprintf("%d signals, %d modelled interventions",
                                           length(cfg$signals), length(cfg$interventions_modelled)))

ft <- formula_table(cfg)
ok("74 formulas build (19 meas + 19 full + 12 intv + 12 ti_all + 12 ti_trend)",
   nrow(ft) == sum(vapply(cfg$signals, function(sg) length(models_of(sg, cfg)), numeric(1))),
   sprintf("terms %d-%d", min(ft$n_terms), max(ft$n_terms)))
ok("meas == full for the 7 unpaired",
   all(vapply(split(ft, ft$signal), function(z)
     (z$n_paired[1] > 0) || identical(z$formula[1], z$formula[2]), logical(1))))
ok("full >= meas terms for the 12 paired",
   all(vapply(split(ft, ft$signal), function(z)
     z$n_paired[1] == 0 || z$n_terms[z$model == "full"] > z$n_terms[z$model == "meas"],
     logical(1))))

# `o_flag` was dropped from the design on 2026-08-25. It was the only source of
# a `by=` smooth and the only term belonging to neither the measurement nor the
# intervention set, so these three checks stand where the old interaction checks
# did: the term must be gone, no by= smooth may reappear, and the meas/intv term
# sets must partition `full` exactly — which is what makes L_full - L_intv the
# conditional term rather than an approximation to it.
ok("o_flag appears in no formula",
   !any(vapply(cfg$signals, function(sg) any(vapply(models_of(sg, cfg), function(md)
     "o_flag" %in% required_columns(build_formula(sg, md, cfg)), logical(1))), logical(1))))
ok("no formula carries a by= smooth",
   !any(vapply(cfg$signals, function(sg) any(vapply(models_of(sg, cfg), function(md)
     any(grepl("by\\s*=", attr(stats::terms(build_formula(sg, md, cfg)), "term.labels"))),
     logical(1))), logical(1))))
ok("at most one trend smooth per formula",
   all(vapply(cfg$signals, function(sg) max(vapply(models_of(sg, cfg), function(md) {
     lab <- attr(stats::terms(build_formula(sg, md, cfg)), "term.labels")
     sum(grepl("^s\\(trend[,)]", lab))
   }, numeric(1))) <= 1, logical(1))))

# The partition, checked as sets rather than as counts: meas and intv disjoint,
# and their union exactly full. A term that drifted into both would make the
# subtraction double-count it, silently.
tset <- function(sg, md) attr(stats::terms(build_formula(sg, md, cfg)), "term.labels")
paired_sigs <- cfg$signals[vapply(cfg$signals, function(sg)
  length(interventions_of(sg, cfg)) > 0L, logical(1))]
ok("meas and intv term sets are disjoint",
   all(vapply(paired_sigs, function(sg)
     !length(intersect(tset(sg, "meas"), tset(sg, "intv"))), logical(1))))
ok("meas union intv == full, exactly",
   all(vapply(paired_sigs, function(sg)
     setequal(c(tset(sg, "meas"), tset(sg, "intv")), tset(sg, "full")), logical(1))))
ok("intv is refused for the 7 unpaired signals",
   all(vapply(setdiff(cfg$signals, paired_sigs), function(sg)
     inherits(try(build_formula(sg, "intv", cfg), silent = TRUE), "try-error"), logical(1))))

# n_obs is deliberately not a term anywhere: monitoring frequency is a
# site-specific ordering behaviour and would transport as an artifact. n enters
# only inside pi_hat, as the confidence weighting.
ok("n_obs is in no formula",
   !any(vapply(cfg$signals, function(sg) any(vapply(models_of(sg, cfg), function(md)
     "n_obs" %in% required_columns(build_formula(sg, md, cfg)), logical(1))), logical(1))))

# The LEVEL term follows excursion_side; the PROPENSITY terms do not. pi is a
# simplex and both coordinates are free, so both enter except where config's
# signal_tails names one as structurally pinned.
# TWO independent rules now. `level_terms` picks the PAIR (quantile vs extreme);
# `excursion_side` picks WHICH OF THE PAIR. Conflating them is how a sparse lab
# ends up carrying a percentile it has no measurements to define.
ALL_LEVEL <- c("q05", "q95", "value_min", "value_max")
# Which magnitude quantities a formula carries, whatever they are NAMED. With
# `magnitude_conditional` on they arrive as `{var}_delta`, so the checks below
# resolve back to the base variable and keep testing the design rule rather
# than the spelling.
mag_bases <- function(sg) {
  req <- required_columns(build_formula(sg, "meas", cfg))
  d <- req[grepl("_delta$", req)]
  intersect(c(intersect(req, ALL_LEVEL), if (length(d)) delta_base_of(d)), ALL_LEVEL)
}
ok("level terms follow excursion_side", all(vapply(cfg$signals, function(sg) {
  side <- excursion_side_of(sg, cfg)
  lv <- level_vars_of(sg, cfg)
  setequal(mag_bases(sg),
           if (is.na(side)) lv else if (side == "low") lv[1] else lv[2])
}, logical(1))))
ok("level terms follow the quantile/extreme declaration",
   all(vapply(cfg$signals, function(sg)
     !length(intersect(mag_bases(sg), setdiff(ALL_LEVEL, level_vars_of(sg, cfg)))),
     logical(1))))
# The raw magnitude term must never appear beside its own conditional version.
ok("no raw magnitude term survives an active conditional construct",
   all(vapply(cfg$signals, function(sg) {
     if (!magnitude_conditional_for(sg, cfg)) return(TRUE)
     req <- required_columns(build_formula(sg, "meas", cfg))
     !length(intersect(req, ALL_LEVEL))
   }, logical(1))))
ok("every sparse signal and the GCS triple use extremes, not quantiles",
   all(vapply(c(cfg$signals[vapply(cfg$signals, function(s)
                  signal_class_of(s, cfg) == "sparse", logical(1))],
                "gcs_motor", "gcs_eyes", "gcs_verbal"),
              function(sg) level_scale_of(sg, cfg) == "extreme", logical(1))))
ok("the dense vitals keep quantiles",
   all(vapply(c("mbp", "heart_rate", "spo2", "resp_rate", "temperature"),
              function(sg) level_scale_of(sg, cfg) == "quantile", logical(1))))
# trend is class-gated by config, and audit 6 is why no sparse signal is in it.
ok("trend follows trend_classes", all(vapply(cfg$signals, function(sg) {
  has <- "trend" %in% required_columns(build_formula(sg, "meas", cfg))
  has == trend_enabled_for(sg, cfg)
}, logical(1))))
ok("no sparse signal carries trend (SQL audit 6)",
   !any(vapply(cfg$signals, function(sg)
     signal_class_of(sg, cfg) == "sparse" && trend_enabled_for(sg, cfg), logical(1))))
ok("pi terms follow signal_tails, not excursion_side", all(vapply(cfg$signals, function(sg) {
  req <- required_columns(build_formula(sg, "meas", cfg))
  tl <- occupiable_tails_of(sg, cfg)
  setequal(intersect(req, c("pi_minus", "pi_plus")),
           c(if ("low" %in% tl) "pi_minus", if ("high" %in% tl) "pi_plus"))
}, logical(1))))
# The case that separates the two rules: MAP is paired low, but keeps pi_plus.
mbp_req <- required_columns(build_formula("mbp", "meas", cfg))
mbp_lo <- if (magnitude_conditional_for("mbp", cfg)) delta_name_of("q05") else "q05"
mbp_hi <- if (magnitude_conditional_for("mbp", cfg)) delta_name_of("q95") else "q95"
ok("mbp keeps pi_plus despite side = low",
   all(c(mbp_lo, "pi_minus", "pi_plus") %in% mbp_req) && !mbp_hi %in% mbp_req,
   sprintf("magnitude term: %s", mbp_lo))
ok("pi_mid never enters",
   !any(vapply(cfg$signals, function(sg) any(vapply(c("meas", "full"), function(md)
     "pi_mid" %in% required_columns(build_formula(sg, md, cfg)), logical(1))), logical(1))))

# Hard rule 4 is an assertion, so test that it actually fires.
ok("guard refuses dx_", inherits(try(assert_no_forbidden("s(dx_value_min)"), silent = TRUE), "try-error"))
# FORBIDDEN_VARS: excluded by NAME, because they arrive wearing an intervention
# prefix and `^dx_` cannot reach them. Each of these was previously excluded by
# omission or by a comment, which is not a guard.
for (fv in c("vasopressor__first_hour", "fio2__peak_intensity",
             "vasopressor__max_concurrent_agents", "invasive_vent__ever_active",
             "n_obs", "pi_mid")) {
  ok(sprintf("guard refuses %s", fv),
     inherits(try(assert_no_forbidden(sprintf("s(%s, bs = 'ts', k = 10)", fv)),
                  silent = TRUE), "try-error"))
}
ok("guard refuses qc_", inherits(try(assert_no_forbidden("qc_discharge_hospice"), silent = TRUE), "try-error"))
ok("guard admits legal terms", isTRUE(assert_no_forbidden(c("s(value_median, bs = \"ts\", k = 10)", "vasopressor__present_at_admission"))))

# --- C. Load and validate ---------------------------------------------------
sect("C. load + nine validator checks")

t0 <- Sys.time()
tabs <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
rep <- validate_tables(tabs, cfg, strict = FALSE)
ok("no validator FAILs", sum(rep$status == "FAIL") == 0L,
   sprintf("%d pass, %d skip, %.1fs", sum(rep$status == "PASS"),
           sum(rep$status == "SKIP"), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- D. Folds ---------------------------------------------------------------
sect("D. split and folds")

folds <- assign_folds(tabs$cohort, cfg)
ok("check_folds passes", isTRUE(check_folds(folds, cfg)))
ok("reproducible", identical(folds, assign_folds(tabs$cohort, cfg)))
ok("test is untouched by folds", all(is.na(folds$fold[folds$split == "test"])))
print(fold_summary(folds), row.names = FALSE)

# --- E. Training priors -----------------------------------------------------
sect("E. Dirichlet alpha and p_bar, training rows only")

pri <- layer1_priors(tabs, folds, cfg, verbose = TRUE)
sp  <- pri$signal
ok("one row per signal x (5 folds + final)", nrow(sp) == length(cfg$signals) * 6L)
ok("all converged", all(sp$converged),
   sprintf("%d-%d iterations", min(sp$iter), max(sp$iter)))
ok("alpha strictly positive", all(sp[, c("alpha_low", "alpha_mid", "alpha_high")] > 0))

# --- conditional-prior constructs (R/04b) -----------------------------------
# Both constructs are switchable; when off the tables are empty and the checks
# below are vacuously true, which is correct for an ablation run.
mp <- pri$magnitude; ip <- pri$intervention
n_mag_var <- sum(vapply(cfg$signals, function(sg)
  if (magnitude_conditional_for(sg, cfg)) length(level_vars_of(sg, cfg)) else 0L, integer(1)))
ok("magnitude priors: one row per (signal, magnitude var) x 6",
   nrow(mp) == n_mag_var * 6L, sprintf("%d rows", nrow(mp)))
ok("magnitude: none degenerate", !nrow(mp) || !any(mp$degenerate),
   if (nrow(mp)) sprintf("%d/%d variance fits converged", sum(mp$converged), nrow(mp)) else "off")
ok("magnitude: variance components non-negative",
   !nrow(mp) || all(mp$s_u >= 0 & mp$s_e >= 0),
   sprintf("%d of %d shrink; the rest are standardised only", sum(mp$shrinks), nrow(mp)))
# The failure this guards against is silent: an unidentified variance split can
# zero a column that check_model_frame still accepts as non-constant.
ok("no magnitude column is numerically null",
   !nrow(mp) || all(vapply(cfg$signals, function(sg) {
     if (!magnitude_conditional_for(sg, cfg)) return(TRUE)
     d <- signal_frame(sg, "meas", tabs, cfg, priors_for(pri, sg, "final"),
                       stay_ids = folds$stay_id[folds$split == "train"])
     all(vapply(grep("_delta$", names(d), value = TRUE),
                function(v) stats::sd(d[[v]]) > 1e-6, logical(1)))
   }, logical(1))))

ok("intensity priors: one row per intervention x 6",
   nrow(ip) == length(unique(ip$intervention)) * 6L,
   sprintf("%d rows: %s", nrow(ip), paste(unique(ip$intervention), collapse = ", ")))
ok("intensity: none degenerate", !nrow(ip) || !any(ip$degenerate))
ok("intensity: every family is known",
   !nrow(ip) || all(ip$family %in% c("binomial", "lognormal")))
if (nrow(ip)) {
  cat("
  intensity models, final fit:
")
  fin <- ip[ip$role == "final", , drop = FALSE]
  for (i in seq_len(nrow(fin))) {
    cat(sprintf("   %-22s %-14s n=%6d  %s
", fin$intervention[i], fin$family[i],
                fin$n_stays[i],
                if (fin$family[i] == "binomial")
                  sprintf("M=%d  c1=%+.3f", fin$M[i], fin$c1[i])
                else sprintf("b1=%+.3f  shrinks=%s", fin$b1[i], fin$shrinks[i])))
  }
}

# The molecule pool is DECLARED in config; this holds it to the counts, which is
# what turns the eICU drugname-vs-molecule trap (spec SS7) into a loud failure.
ap <- check_agent_pool(tabs, cfg, strict = FALSE)
ok("intervention_agent_pool covers the observed molecule counts",
   is.null(ap) || all(ap$ok),
   if (is.null(ap)) "not declared" else
     paste(sprintf("%s %d/%d", ap$intervention, ap$max_observed, ap$declared_pool),
           collapse = "; "))

# The tail each signal actually models must be occupied. Occupancy in the
# unmodelled tail is expected, not an error: MAP has hypertensive hours and is
# still paired to vasopressor on the low side.
es <- check_excursion_sides(tabs, cfg, strict = FALSE)
ok("excursion_side points at an occupied side", all(es$ok),
   sprintf("thinnest: %s (%d stays)", es$signal[which.min(es$n_thin)], min(es$n_thin)))

# signal_tails must match occupancy in BOTH directions: declaring a pinned tail
# occupiable smuggles n_obs in as physiology, and declaring a free one empty
# discards a simplex coordinate.
tl <- check_signal_tails(tabs, cfg, strict = FALSE)
ok("signal_tails matches the counts", all(tl$ok),
   sprintf("one-sided: %s", paste(tl$signal[tl$declared != "low+high"], collapse = ", ")))
print(merge(tl[, c("signal", "declared", "stays_low", "stays_high")],
            es[, c("signal", "modelled")], by = "signal", sort = FALSE), row.names = FALSE)

# alpha_j hits the estimator floor exactly where a tail is structurally empty,
# and no formula may ask for such a coordinate.
deg <- unique(sp$signal[sp$degenerate])
ok("degenerate alpha only on tails no formula uses",
   setequal(deg, tl$signal[tl$stays_low == 0 | tl$stays_high == 0]) &&
     all(vapply(deg, function(sg) {
       z <- tl[tl$signal == sg, ]
       (z$stays_low  > 0 || !z$declared_low) && (z$stays_high > 0 || !z$declared_high)
     }, logical(1))), paste(deg, collapse = ", "))

# THE leakage check. Fold f's alpha must be fitted on the other four folds, so
# its n_stays must equal (final n_stays) - (fold f's own measured train stays).
sf <- tabs$signal_features
tr <- folds$stay_id[folds$split == "train"]
fold_of <- folds$fold[match(sf$stay_id, folds$stay_id)]
meas_tr <- sf$n_obs > 0 & sf$stay_id %in% tr
per_fold <- table(sf$signal[meas_tr], fold_of[meas_tr])
bad <- 0L
for (sg in cfg$signals) {
  fin <- sp$n_stays[sp$signal == sg & sp$role == "final"]
  for (f in seq_len(cfg$n_folds)) {
    got <- sp$n_stays[sp$signal == sg & sp$role == "oof" & sp$fold == f]
    if (got != fin - per_fold[sg, as.character(f)]) bad <- bad + 1L
  }
}
ok("oof alpha excludes its own fold", bad == 0L, sprintf("%d mismatches over 95 fits", bad))

# p_bar must be the MEASURED subpopulation's prior, not the cohort's. If it were
# the cohort rate it would be constant across signals, which is the exact bug
# spec §5.5 point 3 exists to prevent.
cohort_rate <- mean(tabs$cohort$mortality)
fin <- sp[sp$role == "final", ]
ok("p_bar varies by signal", stats::sd(fin$p_bar) > 0.01,
   sprintf("cohort %.4f | p_bar range %.4f-%.4f", cohort_rate, min(fin$p_bar), max(fin$p_bar)))
ok("p_bar never equals the cohort rate", all(abs(fin$p_bar - cohort_rate) > 1e-6))

cat("\n  alpha0 (prior strength, hours) and measured-subpopulation prior:\n")
print(data.frame(signal = fin$signal, n = fin$n_stays,
                 alpha0 = round(fin$alpha0, 2), p_bar = round(fin$p_bar, 4)),
      row.names = FALSE)

# --- F. Model frames --------------------------------------------------------
sect("F. model frames — all 74, on the training set")

t0 <- Sys.time()
fr <- frame_report(tabs, cfg, pri, stay_ids = tr, role = "final")
ok("74 frames attempted", nrow(fr) == 74L,
   sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ok("all frames pass check_model_frame", all(fr$status == "ok"),
   sprintf("%d rejected", sum(fr$status != "ok")))
ok("no frame is empty", all(is.na(fr$n_rows) | fr$n_rows > 0))
print(fr[, setdiff(names(fr), "status")], row.names = FALSE)
if (any(fr$status != "ok")) {
  cat("\n  REJECTED FRAMES:\n")
  for (i in which(fr$status != "ok")) cat(sprintf("   %s [%s]  %s\n", fr$signal[i], fr$model[i], fr$status[i]))
}

# config/smooth_k against the measured distinct-value counts, in BOTH
# directions. A thin-plate basis needs k < the number of distinct values, and
# 13 of the 116 smooth covariates are discrete enough to fall below the default
# k = 10 — the three GCS scales, the vasopressor agent counts, and platelet
# transfusion hours. Undeclared and needed means bam() dies at that fit;
# declared and unnecessary means basis flexibility is being thrown away.
sk <- check_smooth_k(tabs, cfg, pri, stay_ids = tr, strict = FALSE)
ok("smooth_k matches the distinct-value counts", all(sk$ok),
   sprintf("%d of %d covariates need an override", sum(sk$needed), nrow(sk)))
cat("\n  smooth covariates carrying a reduced basis:\n")
print(sk[sk$declared, c("signal", "variable", "k", "n_unique")], row.names = FALSE)

# The same check on every FOLD's fitting subset, which is the stricter one: a
# fold holds 4/5 of train and can thin a covariate the full set supports. This
# is the check that would otherwise surface as an mgcv error at fit 173 of 258.
fold_ok <- TRUE
for (fd in seq_len(cfg$n_folds)) {
  ids <- job_ids("oof", fd, folds)$fit_ids
  z <- try(check_smooth_k(tabs, cfg, pri, stay_ids = ids, strict = FALSE), silent = TRUE)
  if (inherits(z, "try-error") || any(z$k >= z$n_unique)) fold_ok <- FALSE
}
ok("every fold's fitting subset supports its declared k", fold_ok)

# The intv block must be a property of the treatment record ALONE, so two
# signals sharing an intervention must receive byte-identical lambda columns.
# If this fails, L_intv has stopped being log p(I|Y=1)/p(I|Y=0).
li <- try(check_lambda_invariance(tabs, cfg, pri, stay_ids = tr, strict = TRUE), silent = TRUE)
ok("lambda is invariant across signals sharing an intervention",
   !inherits(li, "try-error"),
   if (inherits(li, "try-error")) conditionMessage(attr(li, "condition")) else "")

# The layer-1 job table. Data-free, but checked here because it is what the
# targets graph will map over and a wrong budget is expensive to discover late.
jb <- layer1_jobs(cfg)
bd <- layer1_budget(cfg)
cat("\n  layer-1 fit budget:\n")
print(bd, row.names = FALSE)
# THE BUDGET, AS LITERALS ON PURPOSE. Everywhere else in this project a literal
# beside a config-derived value is the F9-to-F11 defect; here it is the whole
# point. This is an INDEPENDENT statement of what the design should cost, and
# deriving it from `layer1_budget()` would make it assert that a function equals
# itself. A change to any of these four numbers is a spec change and must be
# noticed here.
#
# 2026-09-07: 43 -> 64 specs and 258 -> 384 fits, when `full_ti_trend` and
# `full_ti_all` joined LAYER1_MODELS. The old numbers are kept in this comment
# so a document quoting them can be dated.
ok("64 distinct specs (12 paired x 3 + 9 ti_trend + 12 ti_all + 7 unpaired x 1)",
   bd$specs_distinct == 64L, sprintf("%d", bd$specs_distinct))
ok("384 fits total (320 oof + 64 final)", bd$fits_total == 384L,
   sprintf("%d = %d + %d", bd$fits_total, bd$fits_oof, bd$fits_final))
# COUNTED SEPARATELY because the two aliases mean different things. `alias_meas`
# is a signal with no intervention at all: 7 signals x 3 models (`full` and both
# ti arms) x 6 roles. `alias_full` is a PAIRED signal whose interaction set is
# empty -- `full_ti_trend` on creatinine, platelet and hemoglobin, whose class
# carries no `trend` covariate: 3 x 1 x 6.
ok("unpaired full and both ti arms alias meas rather than refitting",
   bd$alias_meas == 126L, sprintf("%d", bd$alias_meas))
ok("full_ti_trend aliases full where there is no trend covariate",
   bd$alias_full == 18L, sprintf("%d", bd$alias_full))
ok("unpaired intv is assigned L = 0, never fitted", bd$assigned_zero == 42L,
   sprintf("%d", bd$assigned_zero))

# --- the interaction models, checked as a PARTITION ---------------------------
# The three primary models partition exactly (checked in section B). The two
# interaction models must be `full` PLUS cross terms and nothing else: no term
# removed, no term duplicated, and every added term a `ti()`. If a ti model ever
# dropped or altered one of `full`'s terms, `L_full_ti_* - L_intv` would stop
# being the conditional term and the whole `cond_ti_*` family would be measuring
# something nobody declared.
.tl <- function(sg, md) attr(stats::terms(build_formula(sg, md, cfg)), "term.labels")
.pti <- Filter(function(sg) length(interventions_of(sg, cfg)) > 0L, cfg$signals)
ok("every ti model is `full` plus cross terms, exactly",
   all(vapply(.pti, function(sg) all(vapply(LAYER1_TI_MODELS, function(md) {
     a <- .tl(sg, "full"); b <- .tl(sg, md)
     identical(b[seq_along(a)], a) && !anyDuplicated(b) &&
       all(grepl("^ti\\(", b[-seq_along(a)]))
   }, logical(1))), logical(1))),
   sprintf("%d paired signals", length(.pti)))
ok("full_ti_trend crosses ONLY trend",
   all(vapply(.pti, function(sg) {
     x <- interaction_terms(sg, cfg, "trend")
     !length(x) || all(grepl("^ti\\(trend,", x))
   }, logical(1))))
ok("full_ti_trend's cross terms are a subset of full_ti_all's",
   all(vapply(.pti, function(sg)
     all(interaction_terms(sg, cfg, "trend") %in% interaction_terms(sg, cfg, "all")),
     logical(1))))
ok("no ti margin carries a larger basis than its own main-effect smooth",
   all(vapply(.pti, function(sg) all(vapply(LAYER1_TI_MODELS, function(md) {
     sp <- smooth_specs(build_formula(sg, md, cfg))
     mx <- tapply(sp$k, sp$variable, max)
     all(vapply(names(mx), function(v) mx[[v]] <= smooth_k_of(sg, v, cfg), logical(1)))
   }, logical(1))), logical(1))))
# `spec_source()` is read by BOTH `layer1_jobs()` and `l_matrix()`, so an alias
# it names must resolve to something that is actually fitted -- otherwise the
# pivot reads an empty slice of `l_long` and stops after 320 fits.
ok("every alias resolves to a fitted spec",
   all(vapply(cfg$signals, function(sg) all(vapply(LAYER1_MODELS, function(md) {
     src <- spec_source(sg, md, cfg); hops <- 0L
     while (!is.na(src) && !identical(src, "zero") && hops < 5L) {
       md <- src; src <- spec_source(sg, md, cfg); hops <- hops + 1L
     }
     is.na(src) || identical(src, "zero")
   }, logical(1))), logical(1))))
ok("every fitted job has a buildable formula",
   all(vapply(which(jb$fit), function(i)
     inherits(try(build_formula(jb$signal[i], jb$model[i], cfg), silent = TRUE),
              "formula"), logical(1))))
ok("no job is both fitted and sourced elsewhere",
   all(is.na(jb$source[jb$fit])) && all(!is.na(jb$source[!jb$fit])))

# Out-of-fold scoping: the rows a fold is fitted on and the rows it is scored on
# must be disjoint, and must together be the whole training set. This is the
# property the entire out-of-fold design rests on, so it is checked directly.
for (fd in seq_len(cfg$n_folds)) {
  ids <- job_ids("oof", fd, folds)
  ok(sprintf("fold %d: fit and predict ids are disjoint", fd),
     !length(intersect(ids$fit_ids, ids$predict_ids)))
  ok(sprintf("fold %d: fit + predict == train", fd),
     setequal(c(ids$fit_ids, ids$predict_ids), folds$stay_id[folds$split == "train"]))
  ok(sprintf("fold %d: no test stay in scope", fd),
     !length(intersect(c(ids$fit_ids, ids$predict_ids), folds$stay_id[folds$split == "test"])))
}
ok("a final job predicts on nothing (test is touched once, elsewhere)",
   is.null(job_ids("final", NA_integer_, folds)$predict_ids))

# The measured-subset policy, checked rather than assumed: a frame must hold
# exactly the measured stays of the population it was asked for, and never a
# test stay.
te <- folds$stay_id[folds$split == "test"]
d <- signal_frame("lactate", "full", tabs, cfg,
                  priors_for(pri, "lactate", "final"), stay_ids = tr)
n_expect <- sum(sf$signal == "lactate" & sf$n_obs > 0 & sf$stay_id %in% tr)
ok("frame == measured subset of the requested stays", nrow(d) == n_expect,
   sprintf("%d rows", nrow(d)))
ok("no test stay leaks into a train frame", !any(d$stay_id %in% te))

# A fold frame built with that fold's oof alpha — the shape 06_layer1 will use.
f1 <- folds$stay_id[folds$split == "train" & folds$fold == 1L]
d1 <- signal_frame("lactate", "full", tabs, cfg,
                   priors_for(pri, "lactate", "oof", 1L), stay_ids = f1)
ok("held-out fold frame builds", nrow(d1) > 0 && !anyNA(d1), sprintf("%d rows", nrow(d1)))
ok("same columns as the training frame", identical(names(d), names(d1)))

# --- G. Run machinery -------------------------------------------------------
sect("G. run directory conventions")

r <- new_run("smoke", cfg, root = file.path(tempdir(), "smoke_runs"), note = "tests/smoke.R")
save_table(r, ft, "formula_table")
save_table(r, fr, "frame_report", subdir = "diagnostics")
save_table(r, pri$signal, "training_priors", subdir = "diagnostics")
if (nrow(pri$magnitude)) save_table(r, pri$magnitude, "magnitude_priors", subdir = "diagnostics")
if (nrow(pri$intervention)) save_table(r, pri$intervention, "intensity_priors", subdir = "diagnostics")
finalize_run(r)
ok("run dir written", file.exists(file.path(r$path, "manifest.yml")))
ok("diagnostics land in diagnostics/", file.exists(file.path(r$path, "diagnostics", "frame_report.rds")))
ok("tables land in tables/", file.exists(file.path(r$path, "tables", "formula_table.rds")))
ok("external run refuses a bundle",
   inherits(try(save_bundle(new_run("external", cfg, root = file.path(tempdir(), "smoke_runs")),
                            list(), "b"), silent = TRUE), "try-error"))

# ----------------------------------------------------------------------------
cat(sprintf("\n%s\nsmoke: %s (%d failure%s)\n", strrep("=", 64),
            if (FAIL == 0L) "ALL CHECKS PASSED" else "FAILURES PRESENT",
            FAIL, if (FAIL == 1L) "" else "s"))
if (FAIL > 0L) quit(status = 1L)
