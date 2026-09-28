# R/08_diagnostics.R ---------------------------------------------------------
# One row per fitted model, extracted while the gam object is still alive.
#
# THIS FILE EXISTS BECAUSE OF HARD RULE 6. Fold fits are transient: fit_one()
# returns predictions and this row, never the gam object. k.check(), edf and
# concurvity() all need the object, so if they are not pulled out here they
# cannot be recovered later without re-running every fold fit. Diagnostics are
# therefore extracted INSIDE the fit, not downstream of it.
#
# NO CLOCK (hard rules 7 and 9). Nothing in this file records elapsed time, and
# nothing here may. A diagnostics row becomes a target VALUE, and the moment a
# target's value contains a wall-clock number its hash changes on every run —
# so every downstream target rebuilds and the whole cache is gone. Timing is a
# property of a run, not of a result: benchmark_fit() in R/06_layer1.R measures
# it at the call site and never returns it into the graph.
#
# Aggregates only, never a row (hard rule 1). Every number here is a count, a
# rate, or a fitted statistic.
# ----------------------------------------------------------------------------

# --- the row ----------------------------------------------------------------

#' Diagnostics for one fitted layer-1 model.
#'
#' @param b       a fitted mgcv bam/gam object
#' @param signal,model,role,fold  what was fitted. `fold` is NA for final fits.
#' @param cfg     from load_config(); supplies the k used, for edf_ratio
#' @param extras  named list merged into the row (n_fit, n_predict, ...)
#' @return a one-row data frame
model_diagnostics <- function(b, signal, model, role, fold = NA_integer_,
                              cfg = NULL, extras = list()) {
  k_cfg <- cfg$bam$k %||% 10

  kc <- .k_check(b)
  cc <- .concurvity(b)

  # bam() reports convergence differently depending on method. Take whichever
  # is present and record NA rather than guessing TRUE — an unconverged fit
  # that reports itself as converged is the one failure mode worth a false
  # alarm over.
  conv <- b$converged
  if (is.null(conv)) conv <- b$outer.info$conv
  conv <- if (is.null(conv)) NA else isTRUE(conv) || identical(conv, "full convergence")

  edf_total <- tryCatch(sum(b$edf), error = function(e) NA_real_)
  n_smooth  <- length(b$smooth %||% list())

  row <- data.frame(
    signal      = signal,
    model       = model,
    role        = role,
    fold        = as.integer(fold),
    n_rows      = as.integer(stats::nobs(b)),
    n_events    = as.integer(sum(b$y)),
    n_terms     = length(attr(stats::terms(stats::formula(b)), "term.labels")),
    n_smooth    = as.integer(n_smooth),
    edf_total   = round(edf_total, 3),
    dev_expl    = round(.dev_expl(b), 5),
    aic         = round(tryCatch(stats::AIC(b), error = function(e) NA_real_), 2),
    scale_est   = round(b$scale %||% NA_real_, 5),
    converged   = conv,
    k_index_min = kc$k_index_min,
    k_index_p   = kc$p_min,
    k_worst     = kc$worst_term,
    edf_ratio_max = kc$edf_ratio_max,
    edf_worst   = kc$edf_worst_term,
    concurvity_max = cc$worst,
    concurvity_obs = cc$observed,
    concurvity_est = cc$estimate,
    concurvity_term = cc$term,
    notes       = paste(c(kc$note, cc$note), collapse = "; "),
    stringsAsFactors = FALSE
  )

  for (nm in names(extras)) row[[nm]] <- extras[[nm]]
  # edf_ratio is edf/(k-1), and k comes from config rather than from the object
  # so that a mismatch between the two shows up as a wrong ratio rather than
  # being silently absorbed.
  attr(row, "k_config") <- k_cfg
  row
}

#' Deviance explained, the same quantity summary.gam() prints.
.dev_expl <- function(b) {
  tryCatch((b$null.deviance - b$deviance) / b$null.deviance,
           error = function(e) NA_real_)
}

#' mgcv::k.check(), reduced to the four numbers worth thresholding.
#'
#' Returns the WORST smooth on each axis rather than a per-smooth table: with
#' 258 fits x up to 11 smooths there is no reviewing a full table by eye, and
#' the triage question is only ever "does this model have a problem". The term
#' name is carried so a flagged row says where to look.
.k_check <- function(b) {
  na <- list(k_index_min = NA_real_, p_min = NA_real_, worst_term = NA_character_,
             edf_ratio_max = NA_real_, edf_worst_term = NA_character_, note = NULL)
  if (!length(b$smooth %||% list())) {
    return(utils::modifyList(na, list(note = "no smooths")))
  }
  kc <- tryCatch(mgcv::k.check(b), error = function(e) e)
  if (inherits(kc, "error")) {
    return(utils::modifyList(na, list(note = paste0("k.check failed: ", conditionMessage(kc)))))
  }
  if (is.null(kc) || !nrow(kc)) return(utils::modifyList(na, list(note = "k.check empty")))

  kp   <- kc[, "k'"]
  edf  <- kc[, "edf"]
  kidx <- kc[, "k-index"]
  pval <- kc[, "p-value"]
  # k' is the effective basis dimension after identifiability constraints, so
  # edf/k' is the right saturation measure; the config threshold is written as
  # edf/(k-1), which is the same number for a centred smooth.
  ratio <- ifelse(kp > 0, edf / kp, NA_real_)

  i <- which.min(kidx)
  j <- which.max(ratio)
  list(
    k_index_min    = round(unname(kidx[i]), 4),
    p_min          = round(min(pval, na.rm = TRUE), 4),
    worst_term     = rownames(kc)[i],
    edf_ratio_max  = round(unname(ratio[j]), 4),
    edf_worst_term = rownames(kc)[j],
    note           = NULL
  )
}

#' Pairwise concurvity, mgcv::concurvity(full = FALSE).
#'
#' full = FALSE is the pairwise version and is the one that names a culprit
#' pair; full = TRUE only says "this term is explained by the rest". The
#' diagonal is dropped because a term is trivially concurve with itself.
#'
#' ALL THREE measures are recorded, not only the one the threshold reads.
#' `worst` is an upper bound and is deliberately pessimistic: it asks how bad
#' concurvity could be for the worst function in the span, not how bad it is for
#' the function actually fitted. `observed` and `estimate` are the realised
#' versions. Recording `worst` alone makes a flagged row unreadable, because the
#' reader cannot separate a genuine collinearity from the bound being loose.
#'
#' The THRESHOLD still reads `worst`, exactly as config declares. Thresholds are
#' fixed before fitting and are never retuned after seeing which models trip
#' them (CLAUDE.md, frozen decisions) — this adds information for triage and
#' changes no decision rule.
.concurvity <- function(b) {
  na3 <- list(worst = NA_real_, observed = NA_real_, estimate = NA_real_,
              term = NA_character_, note = NULL)
  if (length(b$smooth %||% list()) < 2L) return(na3)

  cc <- tryCatch(mgcv::concurvity(b, full = FALSE), error = function(e) e)
  if (inherits(cc, "error")) {
    return(utils::modifyList(na3, list(
      note = paste0("concurvity failed: ", conditionMessage(cc)))))
  }
  w <- cc[["worst"]]
  if (is.null(w)) return(utils::modifyList(na3, list(note = "no worst matrix")))

  offdiag_max <- function(m) {
    if (is.null(m)) return(NA_real_)
    diag(m) <- NA_real_
    if (all(is.na(m))) NA_real_ else max(m, na.rm = TRUE)
  }
  diag(w) <- NA_real_
  if (all(is.na(w))) return(na3)
  idx <- which(w == max(w, na.rm = TRUE), arr.ind = TRUE)[1, , drop = TRUE]
  list(worst    = round(max(w, na.rm = TRUE), 4),
       observed = round(offdiag_max(cc[["observed"]]), 4),
       estimate = round(offdiag_max(cc[["estimate"]]), 4),
       term     = paste0(rownames(w)[idx[["row"]]], " ~ ", colnames(w)[idx[["col"]]]),
       note     = NULL)
}

# --- triage -----------------------------------------------------------------

#' Which fits tripped a threshold, and which one.
#'
#' Thresholds are fixed in config BEFORE fitting and are never chosen after
#' seeing which models trip them (CLAUDE.md, frozen decisions). This function
#' only reads them.
#'
#' A flagged row is not a failed model. It is a row that has to be looked at
#' before its L is used — which is the entire point of emitting one row per fit
#' rather than a summary.
#'
#' @param diag  rbind of model_diagnostics() rows
#' @param cfg   from load_config()
#' @return the flagged subset, with a `flags` column, most-flagged first
triage <- function(diag, cfg) {
  d <- cfg$diagnostics
  f <- vector("list", nrow(diag))

  for (i in seq_len(nrow(diag))) {
    r <- diag[i, ]
    hit <- character(0)
    # k-index below threshold AND significantly so. Either alone is weak
    # evidence; mgcv's own advice is to read them together.
    if (!is.na(r$k_index_min) && !is.na(r$k_index_p) &&
        r$k_index_min < (d$k_index_min %||% 0.9) &&
        r$k_index_p   < (d$k_index_p_max %||% 0.05)) {
      hit <- c(hit, sprintf("k_index=%.3f (p=%.3f) on %s", r$k_index_min, r$k_index_p, r$k_worst))
    }
    if (!is.na(r$edf_ratio_max) && r$edf_ratio_max > (d$edf_ratio_max %||% 0.8)) {
      hit <- c(hit, sprintf("edf_ratio=%.2f on %s", r$edf_ratio_max, r$edf_worst))
    }
    if (!is.na(r$concurvity_max) && r$concurvity_max > (d$concurvity_max %||% 0.8)) {
      hit <- c(hit, sprintf("concurvity worst=%.2f obs=%.2f (%s)",
                            r$concurvity_max, r$concurvity_obs, r$concurvity_term))
    }
    if (!is.na(r$dev_expl) && r$dev_expl < (d$dev_expl_min %||% 0.005)) {
      hit <- c(hit, sprintf("dev_expl=%.4f", r$dev_expl))
    }
    if (isTRUE(d$require_convergence) && !isTRUE(r$converged)) {
      hit <- c(hit, if (is.na(r$converged)) "convergence unknown" else "did not converge")
    }
    if (nzchar(r$notes)) hit <- c(hit, r$notes)
    f[[i]] <- hit
  }

  n <- lengths(f)
  out <- diag[n > 0, , drop = FALSE]
  if (!nrow(out)) return(out[0, , drop = FALSE])
  out$n_flags <- n[n > 0]
  out$flags   <- vapply(f[n > 0], paste, character(1), collapse = "; ")
  out <- out[order(-out$n_flags, out$signal, out$model), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Console summary of a diagnostics table. Counts only.
diagnostics_summary <- function(diag, cfg) {
  tr <- triage(diag, cfg)
  message(sprintf("diagnostics: %d fits, %d flagged (%s)",
                  nrow(diag), nrow(tr),
                  if (nrow(tr)) paste(sort(unique(paste0(tr$signal, "/", tr$model))),
                                      collapse = ", ") else "none"))
  invisible(tr)
}

# --- the permutation null for concurvity ------------------------------------
#
# WHY A NULL IS NEEDED AT ALL.
#
# mgcv::concurvity(full = FALSE) does not return 0 when covariates are
# unrelated. On this model class it largely reads the ATOM STRUCTURE of the
# covariates -- above all the point mass at zero that every intervention column
# carries on unexposed stays, and every measurement column carries on stays with
# no excursion.
#
# MEASURED 2026-08-28, and the demonstration is in tests/concurvity_null.R:
#
#   * one smooth, one covariate of PURE NOISE with a share of rows zeroed, and
#     nothing else in the model -- so there is no second covariate to concurve
#     with -- returns the zero fraction itself, to three decimals:
#     10% -> 0.082, 50% -> 0.489, 68% -> 0.670, 87% -> 0.866, 95% -> 0.948.
#   * two STATISTICALLY INDEPENDENT noise covariates sharing one zero block
#     reach worst 0.974 at 68% zeros and 0.991 at 87%.
#
# So "205 of 215 fits exceeded 0.80" is largely a statement about how many stays
# were never treated, not about collinearity between physiological quantities.
# A raw concurvity number is uninterpretable here without a matched null, and
# the fixed 0.80 threshold in config was mis-specified for this model class.
#
# WHAT THE NULL PRESERVES, AND WHY IT IS THE RIGHT ONE.
#
# The synthetic curve above is parameterised only by the zero fraction, and that
# is not enough: `fio2__exposure_frac` has a SECOND atom at 1.0 (continuously
# ventilated stays) and reads 0.994 against a zero-fraction prediction of 0.46,
# while the event-shaped `n_hours` columns pile up at 1 and 2. So the null
# permutes the covariate's NON-ZERO VALUES AMONG THE ROWS THAT HOLD THEM,
# leaving intact:
#
#   - the zero pattern, exactly -- so the shared support is unchanged
#   - the marginal distribution, exactly -- so every atom survives
#   - every other column, untouched
#
# and destroying only the row-wise association with the other covariates. That
# is precisely the null hypothesis "this pair overlaps no more than two columns
# with these marginals and this shared support must".
#
# THE NULL IS VALID FOR `worst` AND NOT FOR `observed`. This is the one thing to
# get right before quoting a number from it.
#
# `worst` is computed from the MODEL MATRIX alone and ignores the fitted
# coefficients, so permuting a covariate changes exactly the quantity the null
# is about -- basis overlap -- and nothing else. The null is exact.
#
# `observed` uses the fitted coefficients. Permuting a covariate destroys its
# relationship to the OUTCOME as well as to the other covariates, so its fitted
# smooth collapses to a flat function that cannot overlap anything. MEASURED: the
# `observed` null mean falls to ~0.015, and the "excess" becomes almost the whole
# observed value -- a number that says nothing. There is no permutation that
# preserves a covariate's association with the outcome while destroying its
# association with the other covariates, so this cannot be patched.
#
# For `observed`, use the SYNTHETIC null instead -- sections A and B of
# tests/concurvity_null.R, matched on the covariate's zero fraction. That null is
# approximate, because it is parameterised by the zero fraction alone and does not
# reproduce a second atom, but it does not flatten the smooth. Both are reported
# and each is labelled with which null it came from.
#
# COST. One refit per permutation, so this is a CALIBRATION run to be done once
# and saved, never a per-fit diagnostic: at ~2 s/fit and B = 20, a single
# (signal, model) costs ~40 s. tests/concurvity_null.R writes the table; anything
# downstream joins against it rather than recomputing.

#' Permute a covariate's non-zero values among the rows that hold them.
#'
#' `zero_at` names the value treated as the point mass. It is 0 for every column
#' in this design; the argument exists so the choice is visible at the call site
#' rather than assumed.
.permute_preserving_zeros <- function(x, zero_at = 0) {
  nz <- which(x != zero_at)
  if (length(nz) > 1L) x[nz] <- x[sample(nz)]
  x
}

#' The null distribution of pairwise concurvity for one fitted model.
#'
#' For each smooth covariate in turn, permute it as above, refit, and record the
#' whole pairwise matrix. Every pair involving the permuted covariate is then a
#' draw from the null; pairs not involving it are unchanged and are ignored.
#'
#' @param f,d      the formula and model frame that were actually fitted
#' @param cfg      from load_config(); supplies the bam settings
#' @param B        permutations per covariate
#' @param seed     base seed; each replicate uses seed + b so the run is
#'                 reproducible without reseeding the global stream
#' @param measure  "worst" or "observed". Both are returned; this names the one
#'                 the summary columns are built from.
#' @return data frame, one row per (pair, measure) with the observed value, the
#'   null mean and 2.5/97.5 percentiles, and the excess. Aggregates only.
concurvity_null <- function(f, d, cfg, B = 20L, seed = 1L) {
  fit <- function(dat) .bam_fit(f, dat, cfg)
  b0 <- fit(d)
  if (length(b0$smooth %||% list()) < 2L) return(.empty_concurvity_null())

  pair_vec <- function(b) {
    cc <- tryCatch(mgcv::concurvity(b, full = FALSE), error = function(e) NULL)
    if (is.null(cc)) return(NULL)
    out <- list()
    for (m in c("worst", "observed")) {
      x <- cc[[m]]
      if (is.null(x)) next
      for (i in seq_len(nrow(x))) for (j in seq_len(ncol(x))) {
        if (i == j) next
        out[[length(out) + 1L]] <- data.frame(
          measure = m, row = rownames(x)[i], col = colnames(x)[j],
          value = x[i, j], stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  }

  obs <- pair_vec(b0)
  if (is.null(obs)) return(.empty_concurvity_null())
  obs$pair <- paste(obs$row, "~", obs$col)

  # Only smooth covariates are permuted: `para` is the parametric block and has
  # no single column to shuffle.
  vars <- vapply(b0$smooth, function(s) s$term[1], character(1))
  vars <- intersect(unique(vars), names(d))

  draws <- list()
  for (v in vars) {
    for (b in seq_len(B)) {
      dd <- d
      dd[[v]] <- with_seed(seed + b, .permute_preserving_zeros(d[[v]]))
      bb <- try(fit(dd), silent = TRUE)
      if (inherits(bb, "try-error")) next
      pv <- pair_vec(bb)
      if (is.null(pv)) next
      pv$pair <- paste(pv$row, "~", pv$col)
      # A pair is a null draw only if the permuted covariate is in it.
      pv <- pv[grepl(v, pv$pair, fixed = TRUE), , drop = FALSE]
      if (!nrow(pv)) next
      pv$permuted <- v
      draws[[length(draws) + 1L]] <- pv
    }
  }
  if (!length(draws)) return(.empty_concurvity_null())
  D <- do.call(rbind, draws)

  agg <- do.call(rbind, lapply(split(D, list(D$measure, D$pair), drop = TRUE), function(z) {
    data.frame(measure = z$measure[1], pair = z$pair[1], n_draw = nrow(z),
               null_mean = mean(z$value),
               null_lo = unname(stats::quantile(z$value, 0.025)),
               null_hi = unname(stats::quantile(z$value, 0.975)),
               stringsAsFactors = FALSE)
  }))
  m <- match(paste(agg$measure, agg$pair), paste(obs$measure, obs$pair))
  agg$observed <- obs$value[m]
  agg$excess <- agg$observed - agg$null_mean
  # A pair is "real" only if the observed value clears the null's upper tail.
  agg$above_null <- agg$observed > agg$null_hi
  # Carried on every row so a table read out of context cannot be misquoted.
  agg$null_valid <- agg$measure == "worst"
  rownames(agg) <- NULL
  agg[order(-agg$excess), , drop = FALSE]
}

.empty_concurvity_null <- function() {
  data.frame(measure = character(0), pair = character(0), n_draw = integer(0),
             null_mean = numeric(0), null_lo = numeric(0), null_hi = numeric(0),
             observed = numeric(0), excess = numeric(0), above_null = logical(0),
             null_valid = logical(0), stringsAsFactors = FALSE)
}

#' Collapse a concurvity_null() table to one row per (signal, model, measure).
#'
#' What triage wants is not the whole pair matrix but the worst pair and how far
#' above its own null it sits. Reported rather than thresholded: the fixed 0.80
#' in config reads `worst` and is left exactly as declared.
concurvity_excess <- function(nulltab, signal, model) {
  if (!nrow(nulltab)) return(NULL)
  do.call(rbind, lapply(split(nulltab, nulltab$measure), function(z) {
    i <- which.max(z$excess)
    data.frame(signal = signal, model = model, measure = z$measure[1],
               worst_pair = z$pair[i], observed = z$observed[i],
               null_mean = z$null_mean[i], null_hi = z$null_hi[i],
               excess = z$excess[i], above_null = z$above_null[i],
               null_valid = z$null_valid[1],
               n_pairs_above_null = sum(z$above_null),
               n_pairs = nrow(z), stringsAsFactors = FALSE)
  }))
}

#' The dominant atom: the largest share of rows sitting on one value, taken over
#' a model's smooth covariates.
#'
#' This is the parameter the SYNTHETIC null is indexed by, so it is what a model
#' has to be matched on to read that null. Almost always the point mass at zero
#' -- unexposed stays, or stays with no excursion -- but `fio2__exposure_frac`
#' has its dominant atom at 1.0 (continuously ventilated), which is exactly the
#' case a zero-fraction-only null gets wrong.
#'
#' @return list(frac, variable, value)
dominant_atom <- function(f, d) {
  v <- .term_vars(grep("^s[(]", attr(stats::terms(f), "term.labels"), value = TRUE))
  v <- intersect(v, names(d))
  if (!length(v)) return(list(frac = 0, variable = NA_character_, value = NA_real_))
  best <- list(frac = 0, variable = NA_character_, value = NA_real_)
  for (x in v) {
    tb <- table(d[[x]])
    fr <- max(tb) / nrow(d)
    if (fr > best$frac) {
      best <- list(frac = unname(fr), variable = x,
                   value = as.numeric(names(tb)[which.max(tb)]))
    }
  }
  best
}

#' Read the synthetic null (section B of tests/concurvity_null.R) at an
#' arbitrary atom fraction, by linear interpolation between the tabulated points.
#'
#' Approximate by construction: the synthetic curve is indexed by one number and
#' cannot reproduce a second atom. Use it for `observed`, where the permutation
#' null is invalid, and say which null a quoted number came from.
synthetic_null_at <- function(synth, frac, measure = c("observed", "worst")) {
  measure <- match.arg(measure)
  z <- synth[synth$section == "B", , drop = FALSE]
  if (!nrow(z)) return(NA_real_)
  z <- z[order(z$zero_frac), , drop = FALSE]
  stats::approx(z$zero_frac, z[[measure]], xout = frac, rule = 2)$y
}
