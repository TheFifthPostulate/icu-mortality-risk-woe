# tests/coupling_strata.R ----------------------------------------------------
# DOES CONDITIONING ORDER PATIENTS BETTER WHERE TREATMENT IS HELD FIXED?
#
# The companion to tests/coupling_interaction.R, and the cheap half of the
# question. That script asks whether the additive design is LOSING an
# interaction. This one asks a question that does not depend on the answer:
# within a set of patients whose treatment is similar, does `llr_cond` order
# them better than `llr_meas`?
#
# WHY POOLED AUROC CANNOT ANSWER IT. A pooled AUROC averages over treatment
# strata. It is identical whether or not the measurement-risk map depends on the
# stratum, because between-stratum separation dominates the pairs: most
# discordant pairs in this cohort put a heavily treated patient against a lightly
# treated one, and any score carrying treatment information wins those pairs
# without saying anything about physiology. The within-stratum comparison throws
# those pairs away and keeps the ones a clinician actually faces -- two patients
# on the same support, which one is deteriorating.
#
# WHAT THE ARMS ARE, and why `llr_full` is not among them. Within a stratum
# where the intervention covariates barely vary, `L_intv` is nearly constant, so
# `llr_sum` and `llr_cond` differ by nearly a constant and AUROC cannot separate
# them -- that is arithmetic, not a finding. The comparison that carries content
# is `llr_cond` against `llr_meas`, which differ in what the MEASUREMENT smooths
# were fitted alongside. `llr_sum` is reported to show the between-stratum effect
# it gets credit for disappearing.
#
# FITS NOTHING. Everything here is post-processing of `oof_scores` and
# `l_mats_zero`, both already in the targets cache. It reuses `group_metrics()`
# and `delong_test()` from R/09 so a stratum is scored on the identical code
# path as a hospital at eICU.
#
# Aggregates only, never a row (hard rule 1): counts, AUROC, AUPRC, and paired
# test statistics. Scores and L columns are row-level and are never printed.
#
#   Rscript tests/coupling_strata.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(targets); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

cfg    <- tar_read(cfg)
tabs   <- tar_read(tabs)
tr     <- tar_read(train_ids)
y      <- as.integer(tar_read(y_train))
scores <- tar_read(oof_scores)
Lz     <- tar_read(l_mats_zero)

run <- new_run("coupstrat", cfg, note =
  "within-treatment-stratum discrimination, llr_cond vs llr_meas, fits nothing")

MIN_N <- 300L; MIN_EV <- 25L   # declared here, before any stratum is seen

# --- the stratifier ---------------------------------------------------------
#
# TREATMENT BURDEN is the number of DISTINCT interventions a stay ever received,
# over the 13 extracted interventions. It is chosen over a total-hours measure
# for one reason: hours are on incomparable scales across interventions (a
# vasopressor hour and a transfusion hour are not the same unit), so summing
# them would build a stratifier whose bins mean different things at different
# points of the range. A count of distinct interventions is unit-free.
#
# `ever_active` is used here and is NOT a violation of its exclusion from the
# models. It was dropped as a MODEL TERM because it is exactly
# `(intensity > 0)` and therefore redundant with the intensity smooth. As a
# STRATIFIER redundancy is irrelevant -- the bins only have to partition the
# cohort in a way that means something clinically.
iv <- tabs$intervention_features
iv <- iv[iv$stay_id %in% tr, , drop = FALSE]
burden_tab <- tapply(as.integer(iv$ever_active), factor(iv$stay_id, levels = as.character(tr)),
                     sum, default = 0L)
burden <- as.integer(burden_tab[as.character(tr)])
burden[is.na(burden)] <- 0L
burden_bin <- cut(burden, breaks = c(-1, 0, 1, 2, 3, Inf),
                  labels = c("0", "1", "2", "3", "4+"))

cat("\n=== stratifier: distinct interventions ever received ===\n\n")
cat(sprintf("%-8s %8s %8s %10s\n", "burden", "stays", "deaths", "rate"))
for (g in levels(burden_bin)) {
  s <- burden_bin == g
  cat(sprintf("%-8s %8d %8d %9.2f%%\n", g, sum(s), sum(y[s]), 100 * mean(y[s])))
}

# --- A. arm discrimination within each burden stratum -----------------------
cat("\n=== A. discrimination within a treatment-burden stratum ===\n\n")
cat("    A score that wins pooled but not within stratum is winning on\n")
cat("    between-stratum separation, which is treatment assignment.\n\n")

arms <- c("llr_meas", "llr_cond", "llr_sum")
GA <- list()
for (a in arms) {
  gm <- group_metrics(scores[[a]], y, as.character(burden_bin), label = a,
                      min_n = MIN_N, min_events = MIN_EV)
  z <- gm$per_group; z$arm <- a
  GA[[a]] <- z
}
GA <- do.call(rbind, GA)
save_table(run, GA, "strata_arm_auroc", subdir = "diagnostics")

cat(sprintf("%-8s %8s %8s", "burden", "stays", "deaths"))
for (a in arms) cat(sprintf(" %12s", a))
cat("   cond-meas\n")
for (g in levels(burden_bin)) {
  s <- burden_bin == g
  r <- GA[GA$group == g, , drop = FALSE]
  if (!nrow(r) || !all(r$reported)) {
    cat(sprintf("%-8s %8d %8d   (below the declared floor: n>=%d, deaths>=%d)\n",
                g, sum(s), sum(y[s]), MIN_N, MIN_EV))
    next
  }
  cat(sprintf("%-8s %8d %8d", g, r$n[1], r$deaths[1]))
  for (a in arms) cat(sprintf(" %12.4f", r$auroc[r$arm == a]))
  dt <- delong_test(scores$llr_cond[s], scores$llr_meas[s], y[s])
  cat(sprintf("   %+.4f (p=%.3g)\n", dt$delta, dt$p_value))
}

# --- B. the same, per signal, in that signal's own treatment stratum ---------
#
# THE UNEXPOSED CELL IS THE CLEAN TEST, AND IT IS THE UNEXPOSED ONE. That is the
# opposite of the intuitive reading, so here is the algebra.
#
# VERIFIED 2026-09-05 by counting distinct values: among stays unexposed to
# EVERY intervention paired with a signal, all of that signal's intervention
# covariates take exactly one value -- one distinct value for each of the 21
# covariates checked across `mbp`, `creatinine`, `spo2` and `gcs_motor`,
# including each `__lambda` and each `__present_at_admission`. So on that
# subset the whole intervention block is a constant, and
#
#   L_full = b0 + f(M) + g(0)   ->   L_cond = L_full - L_intv = f(M) + const
#   L_meas = b0'+ f~(M)
#
# where `f` is the measurement block fitted ALONGSIDE the interventions and `f~`
# is the same block fitted alone. A constant does not change an ordering, so
# within the unexposed stratum the AUROC difference between `L_cond` and
# `L_meas` is EXACTLY the effect of having fitted the measurement smooths in the
# presence of the intervention block, and nothing else. That is the coupling
# effect, isolated, at the score level.
#
# THE EXPOSED CELL IS CONTAMINATED and is reported for contrast rather than as a
# test. There `g(I) - g~(I)` varies from stay to stay, so the difference mixes
# the change in `f` with a residual intervention term. A large negative value
# there means the two blocks were explaining the same thing and the subtraction
# removed it -- evidence that M and I are strongly dependent, which is not the
# same claim as treatment context changing what a measurement means.
cat("\n=== B. per signal, exposed vs unexposed to its own paired intervention ===\n\n")
cat("    d = AUROC(L_cond column) - AUROC(L_meas column), within stratum.\n")
cat("    UNEXPOSED is the clean cell: the intervention block is constant there,\n")
cat("    so d is exactly the effect of fitting the measurement smooths beside\n")
cat("    the interventions. This is the coupling test.\n")
cat("    EXPOSED is contaminated: g(I) - g~(I) varies, so d there mixes that\n")
cat("    change with a residual intervention term. Read it as how much evidence\n")
cat("    the two blocks SHARE, not as effect modification.\n\n")

paired <- Filter(function(s) length(interventions_of(s, cfg)) > 0L, cfg$signals)
Lm <- Lz$meas; Lc <- Lz$cond
if (is.null(Lm)) stop("l_mats_zero has no `meas` matrix; nothing to compare against")

cat(sprintf("%-18s | %8s %9s %9s %9s %9s | %8s %9s %9s %9s %5s\n",
            "signal", "n_unexp", "d_UNEXP", "auroc_m", "auroc_c", "p",
            "n_expo", "d_expo", "auroc_m", "auroc_c", "n_iv0"))
GB <- list()
for (sg in paired) {
  ivs <- interventions_of(sg, cfg)
  z <- iv[iv$intervention %in% ivs, , drop = FALSE]
  hit <- tapply(as.integer(z$ever_active), factor(z$stay_id, levels = as.character(tr)),
                max, default = 0L)
  expo <- as.integer(hit[as.character(tr)]); expo[is.na(expo)] <- 0L
  expo <- expo > 0L

  # The algebra above asserted rather than assumed: every intervention covariate
  # that REACHES A FORMULA must be constant across the unexposed subset.
  #
  # THE SET IS TAKEN FROM `build_formula()` AND NOT FROM `names(z)`. The
  # intervention feature table is wider than the model: it also carries metadata
  # (`shape`, which labels an intervention state- or event-shaped) and the
  # columns the formula builder explicitly refuses (`first_hour`,
  # `peak_intensity`, `max_concurrent_agents`). Checking the table's columns is
  # a SECOND DECLARATION of the model's covariate set, and the two disagree --
  # `shape` varies among the unexposed for any signal paired with both a
  # state-shaped and an event-shaped intervention, which is `creatinine` and
  # `urine_output_rate`, and says nothing whatever about the algebra.
  #
  # `__lambda` has no raw column: it is a deterministic function of the
  # accumulation column beside it and the frozen parameters, so it is constant
  # whenever that column is, and checking it does not need the priors loaded.
  form_iv <- setdiff(.term_vars(build_formula(sg, "full", cfg)),
                     c("mortality", .term_vars(measurement_terms(sg, cfg))))
  raw_cols <- intersect(unique(sub("^.*__", "", form_iv)), names(z))
  zu <- z[!(as.character(z$stay_id) %in% as.character(tr[expo])), , drop = FALSE]
  n_iv_const <- sum(vapply(raw_cols, function(v)
    length(unique(zu[[v]][!is.na(zu[[v]])])) <= 1L, logical(1)))
  if (n_iv_const < length(raw_cols)) {
    cat(sprintf("  WARNING %s: %d of %d MODELLED intervention covariates VARY\n",
                sg, length(raw_cols) - n_iv_const, length(raw_cols)))
    cat("          among the unexposed, so its unexposed cell is not the clean test.\n")
  }

  one <- function(s) {
    if (sum(s) < MIN_N || sum(y[s]) < MIN_EV || sum(1 - y[s]) < MIN_EV) return(NULL)
    dt <- delong_test(Lc[s, sg], Lm[s, sg], y[s])
    list(n = sum(s), k = sum(y[s]), am = dt$auroc_2, ac = dt$auroc_1,
         d = dt$delta, p = dt$p_value)
  }
  e <- one(expo); u <- one(!expo)
  if (is.null(e) && is.null(u)) next
  GB[[sg]] <- data.frame(
    signal = sg, n_exposed = e$n %||% NA_integer_, deaths_exposed = e$k %||% NA_integer_,
    auroc_meas_exposed = round(e$am %||% NA_real_, 5),
    auroc_cond_exposed = round(e$ac %||% NA_real_, 5),
    d_exposed = round(e$d %||% NA_real_, 5), p_exposed = e$p %||% NA_real_,
    n_unexposed = u$n %||% NA_integer_, deaths_unexposed = u$k %||% NA_integer_,
    auroc_meas_unexposed = round(u$am %||% NA_real_, 5),
    auroc_cond_unexposed = round(u$ac %||% NA_real_, 5),
    d_unexposed = round(u$d %||% NA_real_, 5), p_unexposed = u$p %||% NA_real_,
    n_iv_cols = length(raw_cols), n_iv_const_unexposed = n_iv_const,
    unexposed_cell_clean = n_iv_const == length(raw_cols),
    stringsAsFactors = FALSE)
  cat(sprintf("%-18s | %8s %+9s %9s %9s %9s | %8s %+9s %9s %9s %5s\n", sg,
              if (is.null(u)) "-" else format(u$n),
              if (is.null(u)) "  -" else sprintf("%.4f", u$d),
              if (is.null(u)) "-" else sprintf("%.4f", u$am),
              if (is.null(u)) "-" else sprintf("%.4f", u$ac),
              if (is.null(u)) "-" else sprintf("%.2g", u$p),
              if (is.null(e)) "-" else format(e$n),
              if (is.null(e)) "  -" else sprintf("%.4f", e$d),
              if (is.null(e)) "-" else sprintf("%.4f", e$am),
              if (is.null(e)) "-" else sprintf("%.4f", e$ac),
              sprintf("%d/%d", n_iv_const, length(raw_cols))))
}
GB <- do.call(rbind, GB)
save_table(run, GB, "strata_signal_exposed", subdir = "diagnostics")

# --- C. the top of the risk distribution ------------------------------------
#
# AUROC weights every discordant pair equally and clinical use does not. This is
# the tail the score would be read at.
cat("\n=== C. precision among the highest-scoring stays ===\n\n")
fracs <- c(0.10, 0.05, 0.02, 0.01)
tp <- function(x, f) { k <- ceiling(f * length(x)); mean(y[order(-x)][seq_len(k)]) }
cat(sprintf("%-12s %8s", "arm", "AUROC"))
for (f in fracs) cat(sprintf(" %9s", sprintf("top %.0f%%", 100 * f)))
cat("\n")
TC <- list()
for (a in names(scores)) {
  v <- vapply(fracs, function(f) tp(scores[[a]], f), numeric(1))
  TC[[a]] <- data.frame(arm = a, auroc = round(.auroc(scores[[a]], y), 5),
                        frac = fracs, precision = round(v, 5),
                        n_top = as.integer(ceiling(fracs * length(y))),
                        stringsAsFactors = FALSE)
  cat(sprintf("%-12s %8.4f", a, .auroc(scores[[a]], y)))
  for (x in v) cat(sprintf(" %9.4f", x))
  cat("\n")
}
save_table(run, do.call(rbind, TC), "top_tail_precision", subdir = "diagnostics")
cat(sprintf("\n  event rate %.4f, so 1.0 would be a perfect tail.\n", mean(y)))

finalize_run(run, extra = list(min_n = MIN_N, min_events = MIN_EV,
                               n_signals = length(paired)))
cat(sprintf("\nwritten: %s\n", run$path))
