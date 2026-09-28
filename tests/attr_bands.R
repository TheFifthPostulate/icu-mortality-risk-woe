# tests/attr_bands.R ----------------------------------------------------------
# PER-PATIENT UNCERTAINTY BANDS ON THE EVIDENCE VALUES, AS A POPULATION SUMMARY.
#
# A NON-FITTING CONSUMER of the attribution replicate store. It reads, for each
# LLR method, the anchor (out-of-fold L from the fold fits), the B posterior
# draws (route `posterior`: beta ~ N(beta_hat, Vc) of each frozen out-of-fold
# fit; contrast L3P, estimation uncertainty conditional on the sample) and the
# bootstrap refits (route `bootstrap`: the 38 shared bags; contrast L3, sampling
# variability), and turns each into a 95% band per patient per signal and per
# patient on the sum. The two routes are different estimands and are reported
# side by side, never pooled (config/attribution_eval.yml).
#
# What is reported, aggregates only (hard rule 1):
#   band_width_signal   per signal, over MEASURED patients: median and 90th
#                       percentile of the 95% band width in nats, both routes
#   band_width_sum      the band on llr_sum / llr_meas / llr_cond: quantiles of
#                       width over patients; share of bands containing zero
#   leader_resolution   share of patients whose anchor leader's lower band end
#                       exceeds the runner-up's upper band end (leader resolved
#                       at 95%), and the same at the top-3 boundary
#   stratum_stability   share of patients whose llr_sum band end points fall in
#                       the same 20-bin stratum, decile and quartile of the
#                       frozen MIMIC-IV cut points as the anchor score
#   ratio_post_vs_boot  per signal, median over patients of the posterior band
#                       width divided by the bootstrap band width, so the
#                       reader knows how much the cheaper band understates the
#                       sampling one (the subsample inflates the bootstrap by
#                       about 1.26, stated not corrected)
#
#   Rscript tests/attr_bands.R                          # store from config/attribution_eval.yml or latest attrgen
#   Rscript tests/attr_bands.R --store out/runs/attrgen_... --bundle out/runs/internal_.../bundle.qs2
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(qs2); library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
STORE  <- .opt("--store", "out/runs/attrgen_20260909T115047")
BUNDLE <- .opt("--bundle", yaml::read_yaml("config/internal.yml")$test_look$bundle)
METHODS <- c("llr_full", "llr_meas", "llr_cond")
Z <- c(0.025, 0.975)

cfg_local <- load_config("config/config.yml")
mf <- read.csv(file.path(STORE, "manifest_replicates.csv"), stringsAsFactors = FALSE)
bundle <- load_bundle(BUNDLE, cfg = cfg_local, strict = TRUE, verbose = FALSE)
measured <- qs2::qs_read(file.path(STORE, "measured.qs2"))
rd <- function(key) qs2::qs_read(file.path(STORE, "replicates", paste0(key, ".qs2")))

run <- new_run("attrbands", cfg_local, note = sprintf("per-patient 95%% bands from %s; fits nothing", basename(STORE)))
log_msg(run, "store: ", STORE, "; bundle: ", BUNDLE)

# Quantile bands over a list of matrices, one patient x signal cell at a time,
# done column-wise to keep memory at one 41,250 x B block.
band_of <- function(mats, anchor) {
  n <- nrow(anchor); p <- ncol(anchor); B <- length(mats)
  lo <- hi <- matrix(NA_real_, n, p, dimnames = dimnames(anchor))
  for (j in seq_len(p)) {
    M <- vapply(mats, function(m) m[, j], numeric(n))          # n x B
    q <- apply(M, 1L, stats::quantile, probs = Z, names = FALSE)
    lo[, j] <- q[1, ]; hi[, j] <- q[2, ]
  }
  sums <- vapply(mats, rowSums, numeric(n))                     # n x B
  qs <- apply(sums, 1L, stats::quantile, probs = Z, names = FALSE)
  list(lo = lo, hi = hi, sum_lo = qs[1, ], sum_hi = qs[2, ])
}

cp <- bundle$cutpoints[["llr_sum"]]
bin_of <- function(x) findInterval(x, cp, rightmost.closed = TRUE, all.inside = TRUE)

W_SIG <- list(); W_SUM <- list(); RES <- list(); STR <- list(); RATIO <- list()
for (m in METHODS) {
  keys <- function(route) mf$key[mf$method == m & mf$route == route]
  anchor <- rd(keys("ladder")[1])
  stopifnot(identical(dimnames(anchor), dimnames(measured)))
  for (route in c("posterior", "bootstrap")) {
    ks <- keys(route)
    log_msg(run, sprintf("%s / %s: %d replicates", m, route, length(ks)))
    mats <- lapply(ks, rd)
    for (x in mats) stopifnot(identical(dimnames(x), dimnames(anchor)))
    bd <- band_of(mats, anchor); rm(mats); gc()
    width <- bd$hi - bd$lo
    # per signal over measured patients
    for (j in seq_len(ncol(anchor))) {
      w <- width[measured[, j], j]
      W_SIG[[length(W_SIG) + 1L]] <- data.frame(method = m, route = route, signal = colnames(anchor)[j],
        n_measured = length(w), width_median = stats::median(w), width_p90 = stats::quantile(w, 0.9, names = FALSE),
        width_max = max(w), abs_l_median = stats::median(abs(anchor[measured[, j], j])),
        frac_band_excludes_zero = mean(bd$lo[measured[, j], j] > 0 | bd$hi[measured[, j], j] < 0))
    }
    # the sum
    sw <- bd$sum_hi - bd$sum_lo; s_pt <- rowSums(anchor)
    W_SUM[[length(W_SUM) + 1L]] <- data.frame(method = m, route = route, n = length(sw),
      width_median = stats::median(sw), width_p10 = stats::quantile(sw, 0.1, names = FALSE),
      width_p90 = stats::quantile(sw, 0.9, names = FALSE), width_max = max(sw),
      abs_sum_median = stats::median(abs(s_pt)),
      frac_band_contains_zero = mean(bd$sum_lo <= 0 & bd$sum_hi >= 0),
      frac_width_below_1nat = mean(sw < 1), frac_width_below_2nat = mean(sw < 2))
    # leader resolution at 95%: leader by anchor among measured signals
    A <- anchor; A[!measured] <- -Inf
    ord <- t(apply(A, 1L, order, decreasing = TRUE))
    n <- nrow(A); i1 <- cbind(seq_len(n), ord[, 1]); i2 <- cbind(seq_len(n), ord[, 2]); i3 <- cbind(seq_len(n), ord[, 3]); i4 <- cbind(seq_len(n), ord[, 4])
    ok2 <- A[i2] > -Inf; ok4 <- A[i4] > -Inf
    RES[[length(RES) + 1L]] <- data.frame(method = m, route = route,
      n_with_runner_up = sum(ok2),
      frac_leader_resolved = mean((bd$lo[i1] > bd$hi[i2])[ok2]),
      frac_leader_band_overlaps_runner_up = mean((bd$lo[i1] <= bd$hi[i2])[ok2]),
      n_with_fourth = sum(ok4),
      frac_top3_boundary_resolved = mean((bd$lo[i3] > bd$hi[i4])[ok4]))
    # stratum stability of the sum band on frozen cut points (llr_sum strata for every method, as a common yardstick)
    b_pt <- bin_of(s_pt); b_lo <- bin_of(bd$sum_lo); b_hi <- bin_of(bd$sum_hi)
    STR[[length(STR) + 1L]] <- data.frame(method = m, route = route,
      frac_same_bin20 = mean(b_lo == b_pt & b_hi == b_pt),
      frac_same_decile = mean(ceiling(b_lo / 2) == ceiling(b_pt / 2) & ceiling(b_hi / 2) == ceiling(b_pt / 2)),
      frac_same_quartile = mean(ceiling(b_lo / 5) == ceiling(b_pt / 5) & ceiling(b_hi / 5) == ceiling(b_pt / 5)),
      frac_within_one_bin20 = mean(b_hi - b_lo <= 1),
      median_bins_spanned = stats::median(b_hi - b_lo + 1))
    assign(paste0("width_", route), width)
  }
  for (j in seq_len(ncol(anchor))) {
    r <- width_posterior[measured[, j], j] / width_bootstrap[measured[, j], j]
    r <- r[is.finite(r)]
    RATIO[[length(RATIO) + 1L]] <- data.frame(method = m, signal = colnames(anchor)[j],
      ratio_median = stats::median(r), ratio_q1 = stats::quantile(r, 0.25, names = FALSE),
      ratio_q3 = stats::quantile(r, 0.75, names = FALSE))
  }
}
save_table(run, do.call(rbind, W_SIG), "band_width_signal")
save_table(run, do.call(rbind, W_SUM), "band_width_sum")
save_table(run, do.call(rbind, RES), "leader_resolution")
save_table(run, do.call(rbind, STR), "stratum_stability")
save_table(run, do.call(rbind, RATIO), "ratio_post_vs_boot")
finalize_run(run, extra = list(store = STORE, bundle = BUNDLE, methods = as.list(METHODS),
                                n_posterior = sum(mf$method == "llr_full" & mf$route == "posterior"),
                                n_bootstrap = sum(mf$method == "llr_full" & mf$route == "bootstrap")))
cat("\n=== band on the sum ===\n"); print(do.call(rbind, W_SUM), row.names = FALSE, digits = 3)
cat("\n=== leader resolution ===\n"); print(do.call(rbind, RES), row.names = FALSE, digits = 3)
cat("\n=== stratum stability ===\n"); print(do.call(rbind, STR), row.names = FALSE, digits = 3)
cat("\nrun directory:", run$path, "\n")
