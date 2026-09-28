# run/patient_card.R ---------------------------------------------------------
# ONE PATIENT CARD from the frozen bundle: the composite score (llr_sum, the sum of channel evidence)
# with a posterior uncertainty band, the MIMIC-IV reference-mortality stratum
# the score falls in, and the 19 per-signal evidence values with bands, ordered
# by evidence.
#
# APPLICATION, NOT REFITTING (hard rule 8). Every quantity is a function of the
# bundle's frozen parameters and one stay's covariates: the final GAMs, their
# smoothing-corrected posterior covariance `Vc`, the final priors, `p_bar_train`,
# and the reporting cut points. Nothing is estimated from the stay being scored
# and nothing is estimated from the site.
#
# THE BAND. For each fitted final GAM the linear predictor at the stay is
# X_p beta with X_p the lpmatrix row; beta ~ N(beta_hat, Vc) is the Bayesian
# posterior of the smooth coefficients conditional on the fitted smoothing
# parameters (mgcv's `Vc` carries the smoothing-parameter correction, and it
# must be read OFF THE OBJECT: predict.bam() drops it). B draws of beta give B
# values of L = X_p beta - logit(p_bar_train), and the band is their 2.5th and
# 97.5th percentile. The 19 GAMs are separate fits, so their posteriors are
# independent given the data and the band on the sum is the band of the sum of
# per-signal draws. This is the `posterior` route of tests/attr_replicates.R
# (contrast L3P) applied at a single stay; the bootstrap-refit route (L3) is
# the more conservative estimand and is not computed here. The GAM bootstrap
# of 2026-09-12 measured the analytic band as adequate for resolved terms and
# as reporting false certainty for terms the shrinkage basis toggles, which is
# the stated limitation of this band.
#
# An unmeasured signal has L = 0 by assignment with a zero-width band: absence
# of a measurement is never evidence (CLAUDE.md).
#
# THE STRATUM. The score is placed among the 20 equal-count bins of the MIMIC-IV
# training out-of-fold score (bundle$cutpoints$llr_sum, frozen) and the observed
# training mortality of that bin, its decile and its quartile are reported from
# the reference run's exported bin table. The band's end points are placed the
# same way, so the card says which strata the patient could plausibly occupy.
#
# PROFILES (--profiles). Three stays chosen by the frozen score alone, one per
# reference-risk profile (quartile 1; quartiles 2-3; quartile 4 of the frozen
# MIMIC-IV cut points), each drawn at random among stays with the terms signal
# measured. The outcome never enters the choice. See `score_pool()` below.
#
# TERMS (--terms-signal, default resp_rate). For that signal the card also
# splits L into its terms: each smooth block and parametric column of the
# linear predictor, with a band from the same posterior draws as the channel.
# The terms plus the intercept minus logit(p_bar) reproduce L exactly; the
# script stops if they do not.
#
# AGGREGATES ONLY IN THE CONSOLE (hard rule 1). The card itself, which carries
# one stay's identifier and derived values, is WRITTEN TO A FILE under the run
# directory for the analyst to read; the console prints the path and structural
# checks only.
#
#   Rscript run/patient_card.R                             # random MIMIC-IV test stay
#   Rscript run/patient_card.R --profiles                  # three stays by frozen quartile
#   Rscript run/patient_card.R --site eicu --seed 7        # random eICU stay
#   Rscript run/patient_card.R --stay <stay_id>            # a named MIMIC-IV test stay
#   Rscript run/patient_card.R --site eicu --stay <id> --draws 1000
#   Rscript run/patient_card.R --ref-run out/runs/internal_...   # bin table source
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(arrow); library(yaml); library(qs2)
})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
SITE    <- match.arg(.opt("--site", "mimic_test"), c("mimic_test", "eicu", "eicu_demo"))
STAY    <- .opt("--stay")
SEED    <- as.integer(.opt("--seed", 1L))
B       <- as.integer(.opt("--draws", 400L))
REF_RUN <- .opt("--ref-run")
# The QC run whose support grid flags extrapolated terms (tables/curves_1d.rds:
# per smooth, a grid with a `supported` flag for >= 50 stays and 5 deaths in a
# 5% window). Model-level aggregate, read for flags only.
QC_RUN  <- .opt("--qc-run", "out/runs/gamqc_20260909T180107")
# Frozen charting density for the propensity reference values
# (paper/make_anchor_n.R): median observed hours among measured MIMIC-IV
# training stays, per channel.
ANCHOR_N_PATH <- .opt("--anchor-n", "paper/figs/coefs/anchor_n.csv")
PROFILE_MODE <- "--profiles" %in% args
# One or more signals, comma-separated, e.g. --terms-signal resp_rate,heart_rate
TERMS_SIGNALS <- strsplit(.opt("--terms-signal", "resp_rate"), ",", fixed = TRUE)[[1]]
if (B < 50L) stop("--draws below 50 gives an unusable 95% band", call. = FALSE)

# --- bundle and site, exactly as the apply runners do ------------------------
rc <- yaml::read_yaml("config/internal.yml")
bundle_path <- .opt("--bundle", rc$test_look$bundle)
if (is.null(bundle_path) || !nzchar(bundle_path)) stop("no bundle path", call. = FALSE)
cfg_local <- load_config(rc$config %||% "config/config.yml")
bundle    <- load_bundle(bundle_path, cfg = cfg_local, strict = TRUE, verbose = FALSE)

if (SITE == "mimic_test") {
  cfg  <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
  tabs <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
  folds    <- assign_folds(tabs$cohort, cfg)            # deterministic partition, not a fit
  pool_ids <- folds$stay_id[folds$split == "test"]
  site_lab <- "MIMIC-IV held-out test"
} else {
  # eicu: the credentialed external cohort. eicu_demo: the openly licensed
  # eICU-CRD demo (ODbL v1.0), extracted by the same SQL (sql/eicu_demo/),
  # scored with the same frozen bundle; the only site whose cards may be
  # published.
  ec   <- yaml::read_yaml(if (SITE == "eicu_demo") "config/external_demo.yml" else "config/external.yml")
  if (normalizePath(ec$bundle, mustWork = FALSE) != normalizePath(bundle_path, mustWork = FALSE)) {
    stop("the site config names a different bundle from the card's bundle", call. = FALSE)
  }
  cfg  <- bundle_cfg(bundle, ec$paths)
  tabs <- load_tables(cfg$paths, cfg, site = "eicu", verbose = FALSE)
  pool_ids <- sort(unique(tabs$cohort$stay_id))
  site_lab <- if (SITE == "eicu_demo") "eICU-CRD demo" else "eICU"
}
if (is.null(REF_RUN)) REF_RUN <- dirname(bundle_path)
ref_bins_path <- file.path(REF_RUN, "tables", "risk_bins_llr_sum_oof.csv")
if (!file.exists(ref_bins_path)) {
  stop("reference bin table not found: ", ref_bins_path,
       ". Pass --ref-run <internal run directory that produced the bundle>.", call. = FALSE)
}
ref_bins <- read.csv(ref_bins_path, stringsAsFactors = FALSE)
if (!all(c("bin", "n", "deaths", "obs_rate") %in% names(ref_bins)) || nrow(ref_bins) != 20L) {
  stop("reference bin table has an unexpected shape", call. = FALSE)
}

# --- reference distributions for percentiles ----------------------------------
# WHERE A VALUE SITS WITHIN ITS OWN CHANNEL. tests/channel_semantics.R measured
# that every channel's L is on the common log-odds unit (calibration slope 1 at
# MIMIC), so +0.5 is the same marginal risk from any signal; but it is the 80th
# percentile of BUN evidence and the 95th of platelet evidence, so it is not the
# same degree of abnormality. The percentile answers "how unusual is this value
# for this channel" and is a DISPLAY quantity: it is never summed and never
# enters a score. The reference is the MIMIC-IV training out-of-fold L_full of
# MEASURED stays, read from the reference run's exported long table, so it is
# the same frozen population as the reference strata. The alias rule is the one
# `l_matrix()` applies: an unpaired signal's `full` is its `meas`.
ref_l_path <- file.path(REF_RUN, "tables", "l_oof.rds")
ref_s_path <- file.path(REF_RUN, "tables", "oof_scores.rds")
if (!file.exists(ref_l_path) || !file.exists(ref_s_path)) {
  stop("reference L table or OOF scores not found under ", REF_RUN, "/tables", call. = FALSE)
}
ref_l <- readRDS(ref_l_path)
sup_path <- file.path(QC_RUN, "tables", "curves_1d.rds")
if (!file.exists(sup_path)) stop("support grid not found: ", sup_path, ". Pass --qc-run.", call. = FALSE)
sup_grid <- readRDS(sup_path)[, c("key", "term", "grid_x", "supported")]
if (!file.exists(ANCHOR_N_PATH)) stop("anchor table not found: ", ANCHOR_N_PATH, call. = FALSE)
ANCHOR_N <- read.csv(ANCHOR_N_PATH, stringsAsFactors = FALSE)
#' Is a stay's value of one smooth outside the term's training support? TRUE
#' when it lies beyond the fitted grid, or its nearest grid point is not
#' supported. NA when the QC run has no grid for the term.
outside_support <- function(key, term, x) {
  z <- sup_grid[sup_grid$key == key & sup_grid$term == term, ]
  if (!nrow(z) || is.na(x)) return(NA)
  if (x < min(z$grid_x) || x > max(z$grid_x)) return(TRUE)
  !z$supported[which.min(abs(z$grid_x - x))]
}
ref_dist <- lapply(cfg$signals, function(sg) {
  z <- ref_l$l[ref_l$signal == sg & ref_l$model == "full"]
  if (!length(z)) z <- ref_l$l[ref_l$signal == sg & ref_l$model == "meas"]
  if (!length(z)) stop("no reference L for signal ", sg, call. = FALSE)
  sort(z)
})
names(ref_dist) <- cfg$signals
ref_sum <- sort(as.numeric(readRDS(ref_s_path)[["llr_sum"]]))
if (!length(ref_sum)) stop("reference OOF llr_sum is empty", call. = FALSE)
pctl <- function(x, ref) 100 * mean(ref <= x)

# --- posterior draws ------------------------------------------------------------
#' Draw from N(mu, V) with a clamped spectrum. Same construction as the
#' generator's `rmvn_clamped()` inside tests/attr_replicates.R (the owner);
#' the number of clamped eigenvalues is reported, never swallowed.
rmvn_clamped <- function(nd, mu, V) {
  e <- eigen(V, symmetric = TRUE)
  neg <- sum(e$values < 0)
  d <- sqrt(pmax(e$values, 0))
  Z <- matrix(stats::rnorm(nd * length(mu)), nd, length(mu))
  list(draws = sweep(Z %*% (t(e$vectors) * d), 2L, mu, `+`), n_clamped = neg)
}

sf   <- tabs$signal_features
sigs <- cfg$signals
jobs <- layer1_jobs(cfg)
jobs <- jobs[jobs$fit & jobs$role == "final" & jobs$model %in% c("meas", "full", "intv"), , drop = FALSE]
# The model whose L is a channel's evidence in llr_sum: `full`, or `meas` for an
# unpaired signal (the alias rule of l_matrix()).
full_model_of <- function(sg) if (!is.null(bundle$models[[paste(sg, "full", sep = "/")]])) "full" else "meas"
if (!all(TERMS_SIGNALS %in% sigs)) stop("--terms-signal names an undeclared signal", call. = FALSE)

# --- reference stratum on frozen MIMIC-IV cut points -------------------------
cp <- bundle$cutpoints[["llr_sum"]]
if (length(cp) != 21L) stop("expected 21 llr_sum cut points in the bundle", call. = FALSE)
bin_of <- function(x) findInterval(x, cp, rightmost.closed = TRUE, all.inside = TRUE)
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c0 - h, c0 + h)
}
stratum <- function(bins) {
  s <- ref_bins[ref_bins$bin %in% bins, ]
  k <- sum(s$deaths); n <- sum(s$n); w <- wilson(k, n)
  sprintf("%.1f%% (95%% CI %.1f%%-%.1f%%; %d of %d training stays)", 100 * k / n, 100 * w[1], 100 * w[2], k, n)
}
strat_row <- function(bins) {
  s <- ref_bins[ref_bins$bin %in% bins, ]; k <- sum(s$deaths); n <- sum(s$n); w <- wilson(k, n)
  c(rate = k / n, lo = w[1], hi = w[2])
}
cohort_rate <- sum(ref_bins$deaths) / sum(ref_bins$n)

# --- one card -----------------------------------------------------------------
#' Build one stay's card. Returns the per-signal table, the summary row and,
#' for each of TERMS_SIGNALS, the term-by-term decomposition of its evidence. Nothing is
#' printed: the caller writes the returned tables into the run directory.
make_card <- function(STAY, how, seed) {
  measured <- vapply(sigs, function(sg) {
    any(sf$signal == sg & as.character(sf$stay_id) == STAY & sf$n_obs > 0)
  }, logical(1))

  fitted <- list()      # key "<signal>/<model>" -> list(point, draws, n_clamped, se)
  terms  <- NULL
  extrap <- list()      # signal -> list(n_terms, outside = term labels)
  ctr <- 0L
  for (i in seq_len(nrow(jobs))) {
    sg <- jobs$signal[i]; md <- jobs$model[i]; key <- paste(sg, md, sep = "/")
    if (!measured[[sg]]) next
    b   <- bundle$models[[key]]
    if (is.null(b)) stop("bundle is missing final GAM ", key, call. = FALSE)
    pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
    nd  <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = STAY, stage = "predict")
    if (nrow(nd) != 1L) stop("expected one predict row for ", key, call. = FALSE)
    Xp  <- stats::predict(b, newdata = nd, type = "lpmatrix", discrete = FALSE)
    V   <- if (!is.null(b$Vc)) b$Vc else b$Vp
    ctr <- ctr + 1L
    dr  <- with_seed(seed * 1000L + ctr, rmvn_clamped(B, stats::coef(b), V))
    eta <- as.numeric(Xp %*% t(dr$draws)) - logit(pri$p_bar)
    pt  <- as.numeric(Xp %*% stats::coef(b)) - logit(pri$p_bar)
    fitted[[key]] <- list(point = pt, draws = eta, n_clamped = dr$n_clamped,
                          se_analytic = sqrt(as.numeric(Xp %*% V %*% t(Xp))),
                          used_vc = !is.null(b$Vc))
    # EXTRAPOLATION FLAG, on the model whose L enters the card.
    if (md == full_model_of(sg)) {
      flags <- vapply(b$smooth, function(s) outside_support(key, s$label, as.numeric(nd[[s$term[1]]])), logical(1))
      extrap[[sg]] <- list(n_terms = length(flags),
                           outside = vapply(b$smooth, `[[`, "", "label")[!is.na(flags) & flags])
    }
    # TERM DECOMPOSITION, ANCHORED (R/15_anchors.R). The linear predictor is a
    # sum over coefficient blocks: one per smooth, one column per parametric
    # term, and the intercept. Each term is reported relative to its reference
    # value x0 (a quiet day on the channel): X_p[, block] beta - X_0[, block]
    # beta. The constants go to the REFERENCE EVIDENCE of the model, L_ref =
    # intercept + sum f(x0) - logit(p_bar), so L = L_ref + sum of anchored
    # terms, exactly. Every band comes from the SAME posterior draws as the
    # channel's band, so terms, L_ref and the channel agree draw by draw.
    if (sg %in% TERMS_SIGNALS && md == full_model_of(sg)) {
      cf <- stats::coef(b)
      n_anc <- ANCHOR_N$n_median[match(sg, ANCHOR_N$signal)]
      blocks <- lapply(b$smooth, function(s) {
        x0 <- term_anchor(s$term[1], pri, n_anc)
        list(label = s$label, var = s$term[1], idx = s$first.para:s$last.para,
             x0 = x0, row0 = as.numeric(.smooth_rows(s, x0)))
      })
      in_sm <- unlist(lapply(blocks, `[[`, "idx"))
      para  <- setdiff(seq_along(cf), c(in_sm, match("(Intercept)", names(cf))))
      blocks <- c(blocks, lapply(para, function(j) list(label = names(cf)[j], var = names(cf)[j], idx = j, x0 = 0, row0 = 0)))
      tt <- do.call(rbind, lapply(blocks, function(k) {
        d  <- Xp[, k$idx, drop = FALSE] - matrix(k$row0, nrow = 1, ncol = length(k$idx))
        tp <- as.numeric(d %*% cf[k$idx])
        td <- as.numeric(d %*% t(dr$draws[, k$idx, drop = FALSE]))
        q  <- stats::quantile(td, c(0.025, 0.975), names = FALSE)
        data.frame(signal = sg, model = md, term = k$label, variable = k$var,
                   x = as.numeric(nd[[k$var]]), anchor = k$x0, contribution = tp, lo = q[1], hi = q[2],
                   stringsAsFactors = FALSE)
      }))
      re   <- reference_evidence(b, pri, n_anc)
      re_d <- as.numeric(dr$draws %*% re$row) - logit(pri$p_bar)
      recon <- sum(tt$contribution) + re$value
      if (abs(recon - pt) > 1e-8) stop("anchored decomposition does not reproduce L for ", key, call. = FALSE)
      tt$reference_evidence <- re$value
      tt$ref_lo <- stats::quantile(re_d, 0.025, names = FALSE)
      tt$ref_hi <- stats::quantile(re_d, 0.975, names = FALSE)
      tt$L_channel <- pt
      terms <- rbind(terms, tt)
    }
  }
  for (ts in TERMS_SIGNALS) if (measured[[ts]] && !any(terms$signal == ts)) stop("no term decomposition for ", ts, call. = FALSE)

  # Alias and assignment rules, as l_matrix() applies them: an unpaired
  # signal's `full` IS its `meas`; its `intv` is 0. `cond = full - intv`.
  pick <- function(sg, md) {
    key <- paste(sg, md, sep = "/")
    if (!measured[[sg]]) return(list(point = 0, draws = rep(0, B), n_clamped = 0L, se_analytic = 0, used_vc = NA))
    if (!is.null(fitted[[key]])) return(fitted[[key]])
    if (md == "full") return(pick(sg, "meas"))
    if (md == "intv") return(list(point = 0, draws = rep(0, B), n_clamped = 0L, se_analytic = 0, used_vc = NA))
    stop("no fitted `", md, "` for measured signal ", sg, call. = FALSE)
  }
  per <- do.call(rbind, lapply(sigs, function(sg) {
    fu <- pick(sg, "full"); iv <- pick(sg, "intv"); me <- pick(sg, "meas")
    q  <- stats::quantile(fu$draws, c(0.025, 0.975), names = FALSE)
    data.frame(signal = sg, measured = measured[[sg]],
               L_full = fu$point, lo = q[1], hi = q[2],
               sd_post = stats::sd(fu$draws), se_analytic = fu$se_analytic,
               L_intv = iv$point, L_cond = fu$point - iv$point, L_meas = me$point,
               pct_train = if (measured[[sg]]) pctl(fu$point, ref_dist[[sg]]) else NA_real_,
               pct_lo    = if (measured[[sg]]) pctl(q[1], ref_dist[[sg]]) else NA_real_,
               pct_hi    = if (measured[[sg]]) pctl(q[2], ref_dist[[sg]]) else NA_real_,
               n_clamped = fu$n_clamped,
               n_terms = if (measured[[sg]]) extrap[[sg]]$n_terms else 0L,
               n_terms_outside = if (measured[[sg]]) length(extrap[[sg]]$outside) else 0L,
               terms_outside = if (measured[[sg]]) paste(extrap[[sg]]$outside, collapse = "; ") else "",
               stringsAsFactors = FALSE)
  }))
  full_draws <- sapply(sigs, function(sg) pick(sg, "full")$draws)   # B x 19
  sum_draws  <- rowSums(full_draws)
  llr_sum    <- sum(per$L_full)
  sum_q      <- stats::quantile(sum_draws, c(0.025, 0.975), names = FALSE)
  llr_meas   <- sum(per$L_meas); llr_cond <- sum(per$L_cond)
  sum_pct    <- pctl(llr_sum, ref_sum)
  sum_pct_q  <- c(pctl(sum_q[1], ref_sum), pctl(sum_q[2], ref_sum))

  b_pt <- bin_of(llr_sum); b_lo <- bin_of(sum_q[1]); b_hi <- bin_of(sum_q[2])
  dec  <- ceiling(b_pt / 2); qua <- ceiling(b_pt / 5)
  dec_bins <- (2 * dec - 1):(2 * dec); qua_bins <- (5 * qua - 4):(5 * qua)

  per_o <- per[order(-per$L_full), ]
  fmt <- function(x) formatC(x, format = "f", digits = 2, width = 6)
  y_stay <- tabs$cohort$mortality[match(STAY, as.character(tabs$cohort$stay_id))]
  lines <- c(
    "PATIENT EVIDENCE CARD  (joint evidence model, llr_sum)",
    sprintf("site: %s   stay: %s   selection: %s", site_lab, STAY, how),
    sprintf("bundle: %s   posterior draws: %d   seed: %d", basename(dirname(bundle_path)), B, seed),
    sprintf("observed outcome: %s", ifelse(is.na(y_stay), "unknown", ifelse(y_stay == 1, "died in hospital", "survived"))),
    "",
    sprintf("COMPOSITE SCORE (sum of channel evidence)   llr_sum = %s nats   95%% band [%s, %s]   (llr_meas %s, llr_cond %s)",
            fmt(llr_sum), fmt(sum_q[1]), fmt(sum_q[2]), fmt(llr_meas), fmt(llr_cond)),
    sprintf("  positive = more evidence for death than the MIMIC-IV training cohort rate of %.1f%%", 100 * cohort_rate),
    sprintf("  percentile among MIMIC-IV training stays: %.0f (band end points at percentiles %.0f and %.0f)",
            sum_pct, sum_pct_q[1], sum_pct_q[2]),
    "",
    "REFERENCE MORTALITY STRATUM  (MIMIC-IV training, out-of-fold, frozen cut points)",
    sprintf("  QUARTILE       %2d of 4  : %s", qua, stratum(qua_bins)),
    sprintf("  decile         %2d of 10 : %s", dec, stratum(dec_bins)),
    sprintf("  20-bin stratum %2d of 20 : %s", b_pt, stratum(b_pt)),
    sprintf("  band end points fall in bins %d to %d of 20", b_lo, b_hi),
    "",
    "PER-SIGNAL EVIDENCE  (L_full with 95% posterior band; ordered by evidence, descending)",
    sprintf("  %-18s %8s %8s %8s %12s   %8s %8s %8s  %s", "signal", "L", "lo", "hi", "pctl [lo-hi]", "L_meas", "L_intv", "L_cond", "note"),
    vapply(seq_len(nrow(per_o)), function(i) {
      r <- per_o[i, ]
      note <- if (!r$measured) "unmeasured: L = 0 by assignment"
              else paste(c(if (r$n_terms_outside > 0) sprintf("EXTRAPOLATED: %d of %d terms outside training support (%s)",
                                                             r$n_terms_outside, r$n_terms, r$terms_outside),
                           if (r$n_clamped > 0) sprintf("%d covariance eigenvalue(s) clamped", r$n_clamped)), collapse = "; ")
      pc <- if (is.na(r$pct_train)) "           -" else sprintf("%3.0f [%3.0f-%3.0f]", r$pct_train, r$pct_lo, r$pct_hi)
      sprintf("  %-18s %s %s %s %s   %s %s %s  %s", r$signal, fmt(r$L_full), fmt(r$lo), fmt(r$hi), pc,
              fmt(r$L_meas), fmt(r$L_intv), fmt(r$L_cond), note)
    }, character(1)),
    "",
    if (!is.null(terms)) unlist(lapply(unique(terms$signal), function(ts) {
      tt <- terms[terms$signal == ts, ]
      c(sprintf("TERMS OF %s  (%s model; contribution relative to the reference value, 95%% posterior band)",
                toupper(ts), tt$model[1]),
        vapply(seq_len(nrow(tt)), function(i) {
          r <- tt[i, ]
          sprintf("  %-34s x = %9.3f (ref %7.3f)   %s  [%s, %s]", r$term, r$x, r$anchor, fmt(r$contribution), fmt(r$lo), fmt(r$hi))
        }, character(1)),
        sprintf("  reference evidence %s  [%s, %s]   + terms = L %s", fmt(tt$reference_evidence[1]),
                fmt(tt$ref_lo[1]), fmt(tt$ref_hi[1]), fmt(tt$L_channel[1])),
        "")
    })) else character(0),
    "READING",
    "  L is the log-likelihood ratio contribution of the signal in nats: the shift in log-odds of death",
    "  relative to the training cohort, given that signal's 24-hour profile and intervention context.",
    "  The band is the 2.5th-97.5th percentile of draws from the frozen model's posterior (estimation",
    "  uncertainty of the fitted smooths conditional on the training sample); it does not include",
    "  measurement noise, resampling of the training set, or the unmeasured-equals-zero assignment.",
    "  pctl is where the value sits among MIMIC-IV training stays that had the signal measured, with the",
    "  band's end points mapped the same way: equal L is equal marginal risk from any signal, but not",
    "  equal rarity. Percentiles are never summed.",
    "  llr_sum is an evidence sum, not a calibrated probability; read risk from the reference stratum.",
    "  EXTRAPOLATED marks a channel with a term evaluated outside its training support (fewer than 50",
    "  stays or 5 deaths in a 5% window, or beyond the fitted range): its value and band rest on",
    "  extrapolation of the fitted function and should be read with that caveat.",
    "  TERMS are each relative to a reference value (a quiet day on the channel: no deviation at the",
    "  channel's typical number of observed hours, a tail as extreme as expected, a flat trend, no",
    "  intervention). The reference evidence is the channel's L at every reference value; it plus the",
    "  terms equals L.")

  sb <- strat_row(b_pt); sd_ <- strat_row(dec_bins); sq <- strat_row(qua_bins)
  summary <- data.frame(site = site_lab, llr_sum = llr_sum, sum_lo = sum_q[1], sum_hi = sum_q[2],
    llr_sum_pct = sum_pct, sum_pct_lo = sum_pct_q[1], sum_pct_hi = sum_pct_q[2],
    llr_meas = llr_meas, llr_cond = llr_cond, bin20 = b_pt, bin20_lo = b_lo, bin20_hi = b_hi, decile = dec, quartile = qua,
    bin20_rate = sb[1], bin20_rate_lo = sb[2], bin20_rate_hi = sb[3],
    decile_rate = sd_[1], decile_rate_lo = sd_[2], decile_rate_hi = sd_[3],
    quartile_rate = sq[1], quartile_rate_lo = sq[2], quartile_rate_hi = sq[3],
    cohort_rate = cohort_rate, n_measured = sum(measured), n_draws = B)
  list(lines = lines, per = per_o, summary = summary, terms = terms, fitted = fitted, measured = measured)
}

# --- which stays ----------------------------------------------------------------
# PROFILES. Three stays chosen by the frozen score alone, never by outcome: the
# card is meant to be read without knowing the outcome, so the outcome must not
# pick it. The test pool is scored with the frozen bundle (application: the
# final GAMs and priors at each stay, no fitting), each stay is placed in a
# frozen MIMIC-IV quartile of llr_sum, and one stay is drawn at random within
# each profile among stays with every TERMS_SIGNALS signal measured, so every card carries a
# term decomposition. Low reference risk = quartile 1 (bins 1-5), intermediate
# = quartiles 2-3 (bins 6-15), high = quartile 4 (bins 16-20).
PROFILES <- list(low = 1:5, intermediate = 6:15, high = 16:20)
PROFILE_LAB <- c(low = "low reference risk (quartile 1)", intermediate = "intermediate reference risk (quartiles 2-3)",
                 high = "high reference risk (quartile 4)")

score_pool <- function(ids) {
  ids <- as.character(ids)
  L <- matrix(0, length(ids), length(sigs), dimnames = list(NULL, sigs))
  for (sg in sigs) {
    md  <- full_model_of(sg)
    b   <- bundle$models[[paste(sg, md, sep = "/")]]
    pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
    nd  <- signal_frame(sg, md, tabs, cfg, pri, stay_ids = ids, stage = "predict")
    eta <- as.numeric(stats::predict(b, newdata = nd, type = "link", discrete = FALSE))
    L[match(as.character(nd$stay_id), ids), sg] <- eta - logit(pri$p_bar)
  }
  stats::setNames(rowSums(L), ids)
}

run <- new_run("card", cfg_local, note = sprintf("patient card, %s, %d posterior draws, %s", site_lab, B,
                                                  if (PROFILE_MODE) "three score profiles" else "single stay"))
save_table(run, do.call(rbind, lapply(cfg$signals, function(sg) {
  # The channel reference distribution in nats, so the exhibit can draw each
  # channel's training spread behind its point (aggregate: five quantiles per
  # signal over measured MIMIC-IV training stays).
  q <- stats::quantile(ref_dist[[sg]], c(0.01, 0.10, 0.50, 0.90, 0.99), names = FALSE)
  data.frame(signal = sg, n_ref = length(ref_dist[[sg]]), p01 = q[1], p10 = q[2], p50 = q[3], p90 = q[4], p99 = q[5])
})), "signal_reference", subdir = "tables")

write_card <- function(card, tag) {
  suffix <- if (is.null(tag)) "" else paste0("_", tag)
  writeLines(card$lines, file.path(run$path, paste0("card", suffix, ".txt")))
  save_table(run, card$per, paste0("per_signal", suffix), subdir = "tables")
  save_table(run, cbind(profile = if (is.null(tag)) NA_character_ else tag, card$summary),
             paste0("summary", suffix), subdir = "tables")
  if (!is.null(card$terms)) save_table(run, card$terms, paste0("terms", suffix), subdir = "tables")
}

if (PROFILE_MODE) {
  if (!is.null(STAY)) stop("--profiles and --stay are exclusive", call. = FALSE)
  pool_scores <- score_pool(pool_ids)
  has_term_sig <- Reduce(`&`, lapply(TERMS_SIGNALS, function(ts)
    as.character(pool_ids) %in% as.character(sf$stay_id[sf$signal == ts & sf$n_obs > 0])))
  pool_bins <- bin_of(pool_scores)
  cards <- list(); pool_counts <- list()
  for (p in seq_along(PROFILES)) {
    tag <- names(PROFILES)[p]
    eligible <- names(pool_scores)[pool_bins %in% PROFILES[[tag]] & has_term_sig]
    pool_counts[[tag]] <- data.frame(profile = tag, n_in_profile = sum(pool_bins %in% PROFILES[[tag]]),
                                     n_eligible = length(eligible))
    if (!length(eligible)) stop("no eligible stay in profile ", tag, call. = FALSE)
    stay <- with_seed(SEED + p, sample(eligible, 1L))
    how  <- sprintf("random draw (seed %d) among %d %s stays in the %s with %s measured; chosen by score, not outcome",
                    SEED + p, length(eligible), site_lab, PROFILE_LAB[[tag]], paste(TERMS_SIGNALS, collapse = " and "))
    card <- make_card(stay, how, SEED + p)
    # The pool score and the card score are two computations of one number:
    # the vectorised link prediction and the single-stay lpmatrix product.
    if (abs(card$summary$llr_sum - pool_scores[[stay]]) > 1e-6) {
      stop("pool score and card score disagree for the ", tag, " profile", call. = FALSE)
    }
    if (!card$summary$bin20 %in% PROFILES[[tag]]) stop("card stay falls outside its profile", call. = FALSE)
    write_card(card, tag)
    cards[[tag]] <- card
  }
  save_table(run, do.call(rbind, pool_counts), "profile_pool", subdir = "tables")
} else {
  if (is.null(STAY)) {
    STAY <- with_seed(SEED, sample(as.character(pool_ids), 1L))
    how  <- sprintf("random draw from %d %s stays (seed %d)", length(pool_ids), site_lab, SEED)
  } else {
    if (!STAY %in% as.character(pool_ids)) {
      stop("the requested stay is not in the ", site_lab, " scoring population; ",
           "identifier not echoed (hard rule 1).", call. = FALSE)
    }
    how <- "named on the command line"
  }
  cards <- list(single = make_card(STAY, how, SEED))
  write_card(cards$single, NULL)
}
finalize_run(run, extra = list(site = SITE, n_draws = B, seed = SEED, bundle_path = bundle_path,
                                profiles = PROFILE_MODE, terms_signal = paste(TERMS_SIGNALS, collapse = ","),
                                n_cards = length(cards)))

# --- console: structure only --------------------------------------------------
cat("\ncards written under:", run$path, "\n")
if (PROFILE_MODE) {
  pc <- do.call(rbind, pool_counts)
  cat(sprintf("  pool %d %s stays scored; per profile: %s\n", length(pool_ids), site_lab,
              paste(sprintf("%s %d in profile, %d eligible", pc$profile, pc$n_in_profile, pc$n_eligible), collapse = "; ")))
}
for (tag in names(cards)) {
  cd <- cards[[tag]]; per <- cd$per
  cat(sprintf("  %s: %d of 19 signals measured; %d GAMs evaluated; term rows %d; clamped GAMs %d\n",
              tag, sum(cd$measured), length(cd$fitted), if (is.null(cd$terms)) 0L else nrow(cd$terms),
              sum(per$n_clamped > 0)))
}
