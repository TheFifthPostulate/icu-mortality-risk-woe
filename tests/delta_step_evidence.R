# tests/delta_step_evidence.R ------------------------------------------------
# THE DISCONTINUITY AT k = 0, shown rather than asserted.
#
# `delta`'s conditional mean is `a0 + a1 log(1+k) + a2 log(n)`. That is one
# smooth curve through every stay, and it cannot represent a JUMP. This script
# shows that the relationship between magnitude and deviation count genuinely
# has one, at the single point k = 0, and that no reshaping of a smooth
# function of k could absorb it.
#
# HOW THE COVERAGE CONFOUND IS REMOVED. Plotting mean(v | k) directly would mix
# the count relationship with the coverage relationship, because stays with more
# deviant hours also tend to have been watched longer. So everything below is
# computed on the residual of `v ~ log(n)` -- the magnitude with the coverage
# component already taken out. What remains is the k relationship alone, which
# is the thing the conditional mean gets wrong.
#
# WHAT TO LOOK FOR. On each panel the x axis is log(1+k), so k = 0 sits at
# exactly 0 and every deviant stay sits to its right. The dashed line is what
# the current model fits: ONE line through all the points, dragged toward the
# k = 0 mass. The solid line is the same model plus `I(k > 0)`. The k = 0 point
# is drawn as an open square because it is a single point on the covariate axis,
# not part of a continuum -- that is the whole argument.
#
# AGGREGATES ONLY (hard rule 1). Every plotted point is a mean over at least
# `min_bin` stays, and the table is per (signal, variable). No row reaches the
# console, the figure or this file's output.
#
#   Rscript tests/delta_step_evidence.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({library(mgcv); library(arrow); library(yaml)})
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

min_bin <- 30L      # a plotted point needs at least this many stays
max_k   <- 24L      # k above this is pooled into one bin

cfg   <- load_config("config/config.yml")
tabs  <- load_tables(cfg$paths$mimiciv, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
s     <- .measured_train_rows(tabs, folds, cfg)

run <- new_run("deltastep", cfg, note = paste(
  "Evidence for the k = 0 discontinuity in delta's conditional mean.",
  "Aggregates only; every plotted point is a mean over >=", min_bin, "stays."))

# --- assemble, per (signal, magnitude variable) -----------------------------
specs <- list()
for (sg in unlist(cfg$signals)) {
  if (!magnitude_conditional_for(sg, cfg)) next
  used <- unique(unlist(lapply(models_of(sg, cfg), function(m) {
    rc <- required_columns(build_formula(sg, m, cfg))
    delta_base_of(rc[grepl("_delta$", rc)])
  })))
  for (v in used) specs[[length(specs) + 1L]] <- c(sg, v)
}

one <- function(sg, vr) {
  cv <- delta_count_of(vr)
  z <- s[s$signal == sg, , drop = FALSE]
  v <- z[[vr]]; k <- z[[cv]]; n <- z$n_obs
  keep <- !is.na(v) & !is.na(k) & !is.na(n) & n > 0
  v <- v[keep]; k <- k[keep]; n <- n[keep]
  if (length(v) < 500L || stats::sd(v) == 0 || sum(k > 0) < 100L) return(NULL)

  # coverage removed once, so the k relationship is what is left
  rn <- stats::residuals(stats::lm(v ~ log(n)))
  lk <- log1p(k)
  pos <- k > 0

  fitA <- stats::lm(rn ~ lk)                       # current: one line
  fitB <- stats::lm(rn ~ lk + pos)                 # ANCOVA: two intercepts
  a1A <- unname(stats::coef(fitA)[2])
  a1B <- unname(stats::coef(fitB)[2])
  a3B <- unname(stats::coef(fitB)[3])

  kb <- pmin(k, max_k)
  agg <- do.call(rbind, lapply(sort(unique(kb)), function(j) {
    i <- kb == j
    if (sum(i) < min_bin) return(NULL)
    data.frame(k = j, n_stays = sum(i), mean_r = mean(rn[i]),
               se = stats::sd(rn[i]) / sqrt(sum(i)), stringsAsFactors = FALSE)
  }))
  if (is.null(agg) || nrow(agg) < 4L) return(NULL)

  m0 <- mean(rn[!pos]); m1 <- if (any(k == 1)) mean(rn[k == 1]) else NA_real_
  sdr <- stats::sd(rn)

  list(
    agg = agg, a0A = unname(stats::coef(fitA)[1]), a1A = a1A,
    a0B = unname(stats::coef(fitB)[1]), a1B = a1B, a3B = a3B, sdr = sdr,
    row = data.frame(
      signal = sg, variable = vr, count_var = cv,
      n_stays = length(v), n_k0 = sum(!pos), frac_k0 = round(mean(!pos), 3),
      # what the data says the step is
      emp_mean_k0 = round(m0, 4), emp_mean_k1 = round(m1, 4),
      emp_jump = round(m1 - m0, 4),
      emp_jump_sd = round((m1 - m0) / sdr, 3),
      # what the CURRENT model can express as a step from k=0 to k=1
      model_jump = round(a1A * log(2), 4),
      model_jump_sd = round(a1A * log(2) / sdr, 3),
      # the fitted step, once it is allowed one
      a3_fitted = round(a3B, 4), a3_sd = round(a3B / sdr, 3),
      # how much the slope changes once the k=0 anchor is released
      slope_before = round(a1A, 4), slope_after = round(a1B, 4),
      slope_ratio = round(a1B / a1A, 2),
      stringsAsFactors = FALSE))
}

res <- Filter(Negate(is.null), lapply(specs, function(p) one(p[1], p[2])))
tab <- do.call(rbind, lapply(res, `[[`, "row"))
tab <- tab[order(-abs(tab$emp_jump_sd)), ]
rownames(tab) <- NULL
save_table(run, tab, "delta_step_evidence")

# --- the figure -------------------------------------------------------------
ord <- order(-abs(vapply(res, function(z) z$row$emp_jump_sd, numeric(1))))
pick <- res[ord][seq_len(min(8L, length(res)))]

p <- save_fig(run, "delta_step_discontinuity", width = 11, height = 6.5, dpi = 150)
op <- graphics::par(mfrow = c(2, 4), mar = c(4.2, 4.2, 3.0, 0.8), mgp = c(2.4, 0.8, 0))
for (z in pick) {
  a <- z$agg
  xs <- log1p(a$k)
  lo <- a$mean_r - 1.96 * a$se; hi <- a$mean_r + 1.96 * a$se
  gx <- seq(0, max(xs), length.out = 200)
  yA <- z$a0A + z$a1A * gx
  ylim <- range(c(lo, hi, yA), finite = TRUE)
  plot(xs, a$mean_r, type = "n", xlim = c(-0.15, max(xs) + 0.1), ylim = ylim,
       xlab = "log(1 + k)", ylab = "mean magnitude | coverage removed",
       main = sprintf("%s / %s", z$row$signal, z$row$variable), cex.main = 1.0)
  graphics::abline(h = 0, col = "grey85")
  graphics::segments(xs, lo, xs, hi, col = "grey55")

  # THE CURRENT MODEL: one line, exactly as fitted on the rows. It must pass
  # through the k = 0 mass because 86% of stays are there, and it is then wrong
  # across the whole k > 0 arm.
  graphics::lines(gx, yA, lty = 2, lwd = 2, col = "grey30")
  # THE ANCOVA: the k > 0 arm gets its own intercept, so its slope is free.
  gx1 <- seq(log(2), max(xs), length.out = 200)
  graphics::lines(gx1, (z$a0B + z$a3B) + z$a1B * gx1, lty = 1, lwd = 2.4, col = "black")
  graphics::points(0, z$a0B, pch = 22, cex = 1.6, lwd = 2, bg = "white")
  graphics::points(xs[a$k > 0], a$mean_r[a$k > 0], pch = 19, cex = 0.75)
  k0 <- a$mean_r[a$k == 0]
  if (length(k0)) graphics::points(0, k0, pch = 22, cex = 2.1, lwd = 2.2, bg = "white")
  graphics::legend(if (z$a1A > 0) "topleft" else "bottomleft", bty = "n", cex = 0.68,
                   legend = c("one line (current)", "two intercepts (ANCOVA)",
                              "k = 0 (single point)"),
                   lty = c(2, 1, NA), pch = c(NA, NA, 22), lwd = c(2, 2.4, 2.2),
                   col = c("grey30", "black", "black"))
}
graphics::par(op)
grDevices::dev.off()
log_msg(run, "figure written: ", basename(p))

# --- report -----------------------------------------------------------------
cat("\n=== the jump from k = 0 to k = 1, against what the current model can express ===\n")
cat("  All columns are on the coverage-removed magnitude. `_sd` columns are in\n")
cat("  units of that residual's SD, so they are comparable across signals.\n\n")
print(tab[, c("signal", "variable", "frac_k0", "emp_jump_sd", "model_jump_sd",
              "a3_sd", "slope_before", "slope_after", "slope_ratio")],
      row.names = FALSE)

cat("\n  emp_jump_sd    what the DATA says the k=0 -> k=1 step is\n")
cat("  model_jump_sd  the largest step the CURRENT model can express there,\n")
cat("                 which is a1 * log(2) because log1p is smooth at 0\n")
cat("  a3_sd          the step once the model is allowed one\n")
cat("  slope_ratio    how much the log(1+k) slope changes once the k = 0\n")
cat("                 anchor is released. Far from 1 means the single line was\n")
cat("                 being held by the k = 0 mass rather than fitting k > 0.\n")

cat(sprintf("\n  rows where the empirical jump exceeds what the model can express: %d of %d\n",
            sum(abs(tab$emp_jump_sd) > abs(tab$model_jump_sd), na.rm = TRUE), nrow(tab)))
cat(sprintf("  median |slope_ratio - 1|: %.2f\n",
            stats::median(abs(tab$slope_ratio - 1), na.rm = TRUE)))

finalize_run(run, extra = list(
  n_specs = nrow(tab), min_bin = min_bin,
  median_abs_emp_jump_sd = round(stats::median(abs(tab$emp_jump_sd), na.rm = TRUE), 3)))
cat(sprintf("\n  run directory: %s\n", run$path))
cat(sprintf("  figure:        %s\n\n", file.path(run$path, "figs",
                                                 "delta_step_discontinuity.png")))
