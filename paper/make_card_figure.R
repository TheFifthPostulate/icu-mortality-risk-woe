# paper/make_card_figure.R -----------------------------------------------------
# Draw the patient-card exhibits from a card run's exported tables. The stay
# identifier is never read; the figures show derived evidence values only.
#   Rscript paper/make_card_figure.R out/runs/card_<id>
#
# A single-stay run (tables/per_signal.csv) gives F14_patient_card. A profile
# run (run/patient_card.R --profiles; tables/per_signal_<profile>.csv) gives
# one card per reference-risk profile, F14_patient_card_<profile>, and the term
# figures F17_card_terms_<signal>: the fitted terms of each decomposed signal's final model
# (the QC curve table, as F15) with each profile's stay marked on every term.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
args <- commandArgs(trailingOnly = TRUE)
RUN <- if (length(args)) args[1] else stop("pass the card run directory")
G <- "out/runs/gamqc_20260909T180107/tables"
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; BLUE <- "#2a78d6"
thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        plot.title = element_text(colour = INK, face = "bold", size = 11),
        plot.subtitle = element_text(colour = INK2, size = 9), legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA))
PROFILE_LAB <- c(low = "Low reference risk", intermediate = "Intermediate reference risk", high = "High reference risk")
PROFILE_COL <- c("Low reference risk" = "#1baf7a", "Intermediate reference risk" = "#eda100",
                 "High reference risk" = "#eb6834")
CH_LAB <- c(mbp = "Mean blood pressure", heart_rate = "Heart rate", spo2 = "Oxygen saturation",
            resp_rate = "Respiratory rate", glucose = "Glucose", gcs_motor = "Glasgow motor",
            gcs_eyes = "Glasgow eyes", gcs_verbal = "Glasgow verbal", urine_output_rate = "Urine output rate",
            creatinine = "Creatinine", platelet = "Platelets", hemoglobin = "Hemoglobin",
            temperature = "Temperature", sodium = "Sodium", bicarbonate = "Bicarbonate", bun = "Urea nitrogen",
            wbc = "White cell count", lactate = "Lactate", bilirubin_total = "Bilirubin")
save_both <- function(p, name, w, h) {
  for (ext in c("png", "pdf")) {
    f <- file.path("paper/figs", paste0(name, ".", ext))
    if (ext == "png") ragg::agg_png(f, width = w, height = h, units = "in", res = 300) else grDevices::cairo_pdf(f, width = w, height = h)
    print(p); dev.off()
  }
  cat("wrote", name, "\n")
}

ref_path <- file.path(RUN, "tables", "signal_reference.csv")
has_ref <- file.exists(ref_path)
rq_all <- if (has_ref) read.csv(ref_path, stringsAsFactors = FALSE) else NULL

# --- one card -------------------------------------------------------------------
draw_card <- function(per, sm, title) {
  # Signed descending: most evidence for death at the top, most for survival at
  # the bottom; unmeasured signals sit at zero between them.
  # Paper names; a dagger marks a channel with a term outside training support.
  lab <- unname(ifelse(per$signal %in% names(CH_LAB), CH_LAB[per$signal], per$signal))
  flagged <- if ("n_terms_outside" %in% names(per)) per$n_terms_outside > 0 else rep(FALSE, nrow(per))
  lab[flagged] <- paste0(lab[flagged], " †")
  lab_of <- stats::setNames(lab, per$signal)
  per$signal <- factor(lab, rev(lab[order(-per$L_full)]))
  per$kind <- ifelse(per$measured, "measured", "unmeasured (L = 0 by assignment)")
  # The channel's training distribution, drawn behind each point on the same
  # nats axis: a light strip for the 1st to 99th percentile and a darker one
  # for the 10th to 90th, with a tick at the median.
  STRIP <- "#d9d6cf"; STRIP2 <- "#bfbbb1"
  pA <- ggplot(per, aes(L_full, signal))
  if (has_ref) {
    rq <- rq_all; rq$signal <- factor(unname(lab_of[rq$signal]), levels(per$signal))
    pA <- pA +
      geom_tile(data = rq, aes(x = (p01 + p99) / 2, y = signal, width = p99 - p01, height = 0.64),
                inherit.aes = FALSE, fill = STRIP, colour = NA) +
      geom_tile(data = rq, aes(x = (p10 + p90) / 2, y = signal, width = p90 - p10, height = 0.64),
                inherit.aes = FALSE, fill = STRIP2, colour = NA) +
      geom_tile(data = rq, aes(x = p50, y = signal, width = 0.012, height = 0.64),
                inherit.aes = FALSE, fill = "white", colour = NA)
  }
  pA <- pA +
    geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 1.2, colour = BLUE, alpha = 0.6, orientation = "y") +
    geom_point(aes(shape = kind), colour = BLUE, size = 2.4, fill = "white") +
    scale_shape_manual(values = c("measured" = 16, "unmeasured (L = 0 by assignment)" = 21)) +
    labs(x = "Evidence L (nats) with 95% posterior band", y = NULL,
         title = "Evidence per channel, ordered by evidence",
         subtitle = if (has_ref) paste0("Grey strips: the channel's MIMIC-IV training distribution\n(1st-99th, 10th-90th percentile; white tick = median). Open circles: unmeasured, L = 0.",
                                        if (any(flagged)) "\n† A term of this channel lies outside its training support; the value is extrapolated." else "")
                    else "Open circles: unmeasured channels, L = 0 by assignment") + thm
  # Where the value and its band end points sit within the channel's MIMIC-IV
  # training distribution (measured stays only).
  pp <- per[per$measured, ]
  pP <- ggplot(pp, aes(pct_train, signal)) +
    geom_vline(xintercept = 50, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0, linewidth = 1.2, colour = "#c9722f", alpha = 0.55, orientation = "y") +
    geom_point(colour = "#c9722f", size = 2.4) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Percentile within channel (MIMIC-IV training)", y = NULL,
         title = "Rarity of the value within its channel",
         subtitle = "Band end points mapped through the same reference") +
    thm + theme(axis.text.y = element_blank())
  # Header text, one line per annotation so the quartile line can be bold.
  # Order: quartile first (bold), then decile, then the 20-bin stratum.
  hdr <- data.frame(
    txt = c(sprintf("Composite score (sum of channel evidence)  %.2f nats  [95%% band %.2f to %.2f]", sm$llr_sum, sm$sum_lo, sm$sum_hi),
            sprintf("  percentile among MIMIC-IV training stays: %.0f (band: %.0f to %.0f)", sm$llr_sum_pct, sm$sum_pct_lo, sm$sum_pct_hi),
            "Reference stratum (MIMIC-IV training, frozen cut points):",
            sprintf("  quartile %d of 4: mortality %.1f%% (%.1f-%.1f%%)", sm$quartile, 100 * sm$quartile_rate, 100 * sm$quartile_rate_lo, 100 * sm$quartile_rate_hi),
            sprintf("  decile %d of 10: %.1f%% (%.1f-%.1f%%)", sm$decile, 100 * sm$decile_rate, 100 * sm$decile_rate_lo, 100 * sm$decile_rate_hi),
            sprintf("  bin %d of 20: %.1f%% (%.1f-%.1f%%)   band spans bins %d to %d", sm$bin20, 100 * sm$bin20_rate, 100 * sm$bin20_rate_lo, 100 * sm$bin20_rate_hi, sm$bin20_lo, sm$bin20_hi),
            sprintf("Cohort rate %.1f%%.  %d of 19 channels measured.  %d draws.", 100 * sm$cohort_rate, sm$n_measured, sm$n_draws)),
    face = c("plain", "plain", "plain", "bold", "plain", "plain", "plain"))
  hdr$y <- 1 - (seq_len(nrow(hdr)) - 1) * 0.14
  pB <- ggplot(hdr) + geom_text(aes(x = 0, y = y, label = txt, fontface = face), hjust = 0, vjust = 1,
                                size = 3.2, colour = INK, family = "mono") +
    xlim(0, 1) + ylim(0, 1) + theme_void() + labs(title = title) +
    theme(plot.title = element_text(colour = INK, face = "bold", size = 11))
  pB / (pA | pP) + plot_layout(heights = c(1, 3.2), widths = c(2.2, 1))
}

tabs_dir <- file.path(RUN, "tables")
if (file.exists(file.path(tabs_dir, "per_signal.csv"))) {
  per <- read.csv(file.path(tabs_dir, "per_signal.csv"), stringsAsFactors = FALSE)
  sm  <- read.csv(file.path(tabs_dir, "summary.csv"), stringsAsFactors = FALSE)
  save_both(draw_card(per, sm, sprintf("Patient card, %s", sm$site)), "F14_patient_card", 9, 7.5)
  quit(save = "no")
}

profiles <- names(PROFILE_LAB)[file.exists(file.path(tabs_dir, paste0("per_signal_", names(PROFILE_LAB), ".csv")))]
if (length(profiles) != 3L) stop("expected three profile cards under ", tabs_dir)
for (tag in profiles) {
  per <- read.csv(file.path(tabs_dir, paste0("per_signal_", tag, ".csv")), stringsAsFactors = FALSE)
  sm  <- read.csv(file.path(tabs_dir, paste0("summary_", tag, ".csv")), stringsAsFactors = FALSE)
  save_both(draw_card(per, sm, sprintf("Patient card, %s: %s", sm$site, PROFILE_LAB[[tag]])),
            paste0("F14_patient_card_", tag), 9, 7.5)
}

# --- F17: the three stays on the fitted terms, one figure per decomposed signal --
tm_all <- do.call(rbind, lapply(profiles, function(tag) {
  z <- read.csv(file.path(tabs_dir, paste0("terms_", tag, ".csv")), stringsAsFactors = FALSE)
  z$profile <- factor(PROFILE_LAB[[tag]], PROFILE_LAB); z
}))
# Anchored evidence functions (paper/make_anchored_curves.R): each term relative
# to its reference value, matching the anchored term values of the card.
cv_all <- readRDS("paper/figs/coefs/anchored_curves.rds")
ANCH   <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
pc_all <- read.csv("paper/figs/coefs/parametric_coefs.csv", stringsAsFactors = FALSE)
TERM_LAB <- c("s(pi_minus)" = "Deviation propensity, low side", "s(pi_plus)" = "Deviation propensity, high side",
              "s(q05_delta)" = "Frequency-adjusted tail extremity, low side", "s(q95_delta)" = "Frequency-adjusted tail extremity, high side",
              "s(trend)" = "Trend", "s(invasive_vent__exposure_frac)" = "Invasive ventilation, exposure fraction",
              "invasive_vent__present_at_admission" = "Invasive ventilation, admission indicator",
              "s(vasopressor__exposure_frac)" = "Vasopressor, exposure fraction",
              "s(vasopressor__lambda)" = "Vasopressor, exposure-adjusted intensity",
              "vasopressor__present_at_admission" = "Vasopressor, admission indicator")
Y_LAB <- "Relative to reference (nats)"
UNIT <- c(resp_rate = "breaths/min", heart_rate = "beats/min")
x_label <- function(sg, t) {
  u <- if (sg %in% names(UNIT)) UNIT[[sg]] else "units"
  if (t == "s(pi_minus)") return("Propensity to fall below the reference range")
  if (t == "s(pi_plus)")  return("Propensity to rise above the reference range")
  if (endsWith(t, "_delta)")) return(sprintf("Extremity relative to expected (%s)", u))
  if (t == "s(trend)") return(sprintf("Slope of hourly medians (%s per 24 h)", u))
  if (endsWith(t, "__exposure_frac)")) return("Fraction of the 24-hour window exposed")
  if (endsWith(t, "__lambda)")) return("Intensity relative to what exposure predicts")
  NULL
}
# A term shrunk to (near) zero must look flat: every panel spans at least
# this many nats on its y axis, centred on its own data.
MIN_SPAN <- 0.5
# Propensities first, then tail extremities, then trend, then the intervention terms:
# a tail extremity is read after the propensity on its side.
rank_of <- function(t) ifelse(startsWith(t, "s(pi_"), 1, ifelse(endsWith(t, "_delta)"), 2, ifelse(t == "s(trend)", 3, 4)))

for (sg in unique(tm_all$signal)) {
  tm  <- tm_all[tm_all$signal == sg, ]
  key <- paste(sg, tm$model[1], sep = "/")
  cv  <- cv_all[cv_all$key == key, ]
  pc  <- pc_all[pc_all$key == key, ]
  # Consistency check between the two sources: the QC curve evaluated at a
  # stay's value must equal the stay's term contribution from the bundle. One
  # aggregate number is printed (the largest disagreement), no stay values.
  sm_rows <- tm[tm$term %in% cv$term, ]
  dev <- vapply(seq_len(nrow(sm_rows)), function(i) {
    cz <- cv[cv$term == sm_rows$term[i], ]; cz <- cz[order(cz$grid_x), ]
    if (sm_rows$x[i] < min(cz$grid_x) || sm_rows$x[i] > max(cz$grid_x)) return(NA_real_)
    abs(stats::approx(cz$grid_x, cz$fit, xout = sm_rows$x[i])$y - sm_rows$contribution[i])
  }, numeric(1))
  cat(sprintf("%s: curve vs card term contribution: max |difference| %.4f nats over %d in-grid term values (%d off grid)\n",
              sg, max(dev, na.rm = TRUE), sum(!is.na(dev)), sum(is.na(dev))))

  terms_order <- unique(tm$term); terms_order <- terms_order[order(rank_of(terms_order), seq_along(terms_order))]
  panels <- lapply(terms_order, function(t) {
    pts <- tm[tm$term == t, ]
    if (t %in% cv$term) {
      cz <- cv[cv$term == t, ]; cz <- cz[order(cz$grid_x), ]
      cz$lo <- cz$fit - 1.96 * cz$se; cz$hi <- cz$fit + 1.96 * cz$se
      x0 <- ANCH$anchor[ANCH$key == key & ANCH$term == t][1]
      g <- ggplot() + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
        geom_vline(xintercept = x0, colour = INK2, linewidth = 0.3, linetype = "22") +
        geom_ribbon(data = cz, aes(grid_x, ymin = lo, ymax = hi), fill = BLUE, alpha = 0.12)
      g <- if (cz$kind[1] == "atoms") g + geom_point(data = cz, aes(grid_x, fit), colour = BLUE, size = 1.1, alpha = 0.6)
           else g + geom_line(data = cz, aes(grid_x, fit), colour = BLUE, linewidth = 0.6, alpha = 0.8)
    } else {
      cf <- pc[pc$term == t, ]
      if (nrow(cf) != 1L) stop("no parametric coefficient for ", key, " ", t)
      base <- data.frame(grid_x = c(0, 1), fit = c(0, cf$estimate), lo = c(0, cf$estimate - 1.96 * cf$se),
                         hi = c(0, cf$estimate + 1.96 * cf$se))
      g <- ggplot() + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
        geom_errorbar(data = base, aes(grid_x, ymin = lo, ymax = hi), width = 0.06, colour = BLUE, alpha = 0.5) +
        scale_x_continuous(breaks = c(0, 1), labels = c("not at admission", "at admission"), limits = c(-0.3, 1.3))
    }
    g + geom_errorbar(data = pts, aes(x, ymin = lo, ymax = hi, colour = profile), width = 0, linewidth = 1) +
      geom_point(data = pts, aes(x, contribution, colour = profile), size = 2.4) +
      scale_colour_manual(values = PROFILE_COL, drop = FALSE) +
      expand_limits(y = c(-MIN_SPAN / 2, MIN_SPAN / 2)) +
      labs(title = if (t %in% names(TERM_LAB)) TERM_LAB[[t]] else t, x = x_label(sg, t), y = Y_LAB) + thm +
      theme(plot.title = element_text(size = 9))
  })
  refe <- tm$reference_evidence[1]
  p17 <- wrap_plots(panels, ncol = 3, guides = "collect") +
    plot_annotation(subtitle = sprintf(paste0("Each term relative to its reference value (dotted line). Reference evidence of the ",
                                              "channel: %+.2f nats; a stay's evidence for the channel is this plus its terms."), refe),
                    theme = theme(plot.subtitle = element_text(size = 9, colour = INK2))) &
    theme(legend.position = "bottom", legend.title = element_blank())
  save_both(p17, paste0("F17_card_terms_", sg), 9.5, 3.2 * ceiling(length(panels) / 3) + 0.6)
}
