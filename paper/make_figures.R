# paper/make_figures.R ---------------------------------------------------------
# Preprint figures, drawn ONLY from aggregate run tables (hard rule 1). No data/
# access, no model objects, no row-level scores. Every input is a CSV or RDS of
# aggregates under out/runs/. Run IDs are pinned to the runs of record named in
# docs/preprint_results_ledger_20260912.md.
#
#   Rscript paper/make_figures.R
#
# Writes paper/figs/F<k>_<name>.png (300 dpi) and .pdf.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

I <- "out/runs/internal_20260909T112643/tables"
T <- "out/runs/test_20260912T215546/tables"
E <- "out/runs/external_20260909T112824/tables"
A <- "out/runs/attrmetrics_20260909T165755/diagnostics"
X <- "out/runs/attrext_20260909T175854/diagnostics"            # eICU, frozen bundle only (old F12/F13 blocks)
Z <- "out/runs/attrmetricsext_20261006T121114/diagnostics"     # eICU, whole-bag refits on the 38 internal bags
ATTR_SITE <- c(mimic = A, eicu = Z)
ATTR_SITE_LAB <- c(mimic = "MIMIC-IV train", eicu = "eICU")
one_row <- function(d, what) { if (nrow(d) != 1L) stop(what, ": expected one row, found ", nrow(d)); d }
G <- "out/runs/gamqc_20260909T180107/tables"
OUT <- "paper/figs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(dir, name) read.csv(file.path(dir, paste0(name, ".csv")), stringsAsFactors = FALSE)

# --- palette: fixed slot per entity, never cycled ------------------------------
ARM_LAB <- c(llr_sum = "WoE paired", llr_meas = "WoE measurement", llr_cond = "WoE conditional",
             xgb_raw = "XGBoost raw", xgb_feat = "XGBoost constructed", xgb_l = "XGBoost stacked",
             shap_xgb_feat = "SHAP (XGBoost constructed)", llr_full = "WoE paired")
ARM_COL <- c("WoE paired" = "#2a78d6", "XGBoost raw" = "#eb6834", "WoE measurement" = "#1baf7a",
             "WoE conditional" = "#eda100", "XGBoost stacked" = "#e87ba4", "XGBoost constructed" = "#008300",
             "SHAP (XGBoost constructed)" = "#4a3aa7")
SITE_COL <- c("MIMIC-IV train (out-of-fold)" = "#2a78d6", "MIMIC-IV test (held out)" = "#1baf7a", "eICU" = "#eb6834")
SITE_LEV <- names(SITE_COL)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"

thm <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
        strip.text = element_text(colour = INK, face = "bold", hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(colour = INK, face = "bold", size = 11),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.background = element_rect(fill = "white", colour = NA))
save_fig <- function(p, name, w, h) {
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h)
  print(p); dev.off()
  cat("wrote", name, "\n")
}
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  data.frame(rate = p, lo = c0 - h, hi = c0 + h)
}

# ==== F1: 20-bin mortality curves on frozen MIMIC cut points ===================
bins <- do.call(rbind, lapply(c("llr_sum", "xgb_feat", "xgb_raw"), function(arm) rbind(
  cbind(site = SITE_LEV[1], arm = ARM_LAB[[arm]], rd(I, paste0("risk_bins_", arm, "_oof"))[, c("bin", "n", "deaths", "obs_rate", "lo", "hi")]),
  cbind(site = SITE_LEV[2], arm = ARM_LAB[[arm]], rd(T, paste0("risk_bins_", arm, "_frozen"))[, c("bin", "n", "deaths", "obs_rate", "lo", "hi")]),
  cbind(site = SITE_LEV[3], arm = ARM_LAB[[arm]], rd(E, paste0("risk_bins_", arm, "_frozen"))[, c("bin", "n", "deaths", "obs_rate", "lo", "hi")]))))
bins$site <- factor(bins$site, SITE_LEV)
bins$arm <- factor(bins$arm, c("WoE paired", "XGBoost constructed", "XGBoost raw"))
p1 <- ggplot(bins, aes(bin, obs_rate, colour = site)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, alpha = 0.7,
                position = position_dodge(width = 0.55)) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 0.55)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.55)) +
  facet_wrap(~arm, ncol = 3) +
  scale_colour_manual(values = SITE_COL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Risk bin on MIMIC-IV training cut points (20 equal-count bins of the out-of-fold score)",
       y = "Observed in-hospital mortality (Wilson 95% CI)") + thm
save_fig(p1, "F1_risk_curves_frozen_cuts", 9.5, 3.6)

# ==== F11: reference-risk strata (quartiles) ===================================
strata_of <- function(dir, suffix, arm, site) {
  r <- rd(dir, paste0("risk_bins_", arm, suffix))
  q <- list(Q1 = 1:5, Q2 = 6:10, Q3 = 11:15, Q4 = 16:20)
  do.call(rbind, lapply(names(q), function(nm) {
    s <- r[r$bin %in% q[[nm]], ]; w <- wilson(sum(s$deaths), sum(s$n))
    data.frame(site = site, arm = ARM_LAB[[arm]], stratum = nm, n = sum(s$n), w)
  }))
}
arms5 <- c("llr_sum", "llr_meas", "xgb_l", "xgb_feat", "xgb_raw")   # conditional arm dropped 2026-10-03
st <- do.call(rbind, lapply(arms5, function(a) rbind(
  strata_of(I, "_oof", a, SITE_LEV[1]), strata_of(T, "_frozen", a, SITE_LEV[2]), strata_of(E, "_frozen", a, SITE_LEV[3]))))
st$site <- factor(st$site, SITE_LEV); st$arm <- factor(st$arm, ARM_LAB[arms5])
p11 <- ggplot(st, aes(arm, rate, colour = site)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  facet_wrap(~stratum, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = SITE_COL) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = NULL, y = "Observed mortality (Wilson 95% CI)") + thm +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig(p11, "F11_reference_risk_strata", 9.5, 4)

# ==== F2: within-hospital AUROC and AUPRC distributions ===========================
hp <- readRDS(file.path(E, "hospital_per_group.rds"))
hs <- rd(E, "hospital_summary")
keep_arms <- ARM_LAB[c("llr_sum", "llr_meas", "xgb_l", "xgb_feat", "xgb_raw")]
hosp <- do.call(rbind, lapply(names(keep_arms), function(a) {
  d <- hp[[a]]; d <- d[d$reported %in% TRUE, ]
  rbind(data.frame(arm = ARM_LAB[[a]], metric = "AUROC", v = d$auroc),   # the group column is never read
        data.frame(arm = ARM_LAB[[a]], metric = "AUPRC", v = d$auprc))
}))
n_hosp <- unique(as.vector(table(hosp$arm[hosp$metric == "AUROC"])))
if (length(n_hosp) != 1L) stop("arms differ in their number of reported hospitals")
hosp$arm <- factor(hosp$arm, rev(keep_arms))
hs <- hs[hs$label %in% names(keep_arms), ]
pooled <- rbind(data.frame(arm = ARM_LAB[hs$label], metric = "AUROC", v = hs$pooled_auroc_eligible),
                data.frame(arm = ARM_LAB[hs$label], metric = "AUPRC", v = hs$pooled_auprc_eligible))
pooled$arm <- factor(pooled$arm, rev(keep_arms))
hosp_panel <- function(m, xlab, show_y) {
  ggplot(hosp[hosp$metric == m, ], aes(v, arm, colour = arm)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.4, fill = NA) +
    geom_point(position = position_jitter(width = 0, height = 0.18, seed = 1), size = 0.9, alpha = 0.55) +
    geom_point(data = pooled[pooled$metric == m, ], shape = 23, size = 2.6, fill = "white", stroke = 0.8) +
    scale_colour_manual(values = ARM_COL, guide = "none") +
    labs(x = xlab, y = NULL) + thm +
    (if (show_y) theme() else theme(axis.text.y = element_blank()))
}
p2 <- (hosp_panel("AUROC", "Within-hospital AUROC", TRUE) | hosp_panel("AUPRC", "Within-hospital AUPRC", FALSE)) +
  plot_annotation(subtitle = sprintf(paste0("eICU, %d hospitals with at least 100 stays and 10 deaths. Points: hospitals. ",
                                            "Box: median and IQR. Open diamond: pooled value on the same stays."), n_hosp),
                  theme = theme(plot.subtitle = element_text(colour = INK2, size = 9)))
save_fig(p2, "F2_hospital_auroc_auprc", 9.5, 3.6)

# ==== F3: transport ladder =====================================================
tr <- rd(I, "train_ref_arms"); ts <- rd(T, "score_summary"); es <- rd(E, "score_summary")
ts <- ts[ts$binning == "self", ]; ts$label <- sub("_self$", "", ts$label)
es <- es[es$binning == "frozen", ]; es$label <- sub("_frozen$", "", es$label)
lad <- rbind(data.frame(site = SITE_LEV[1], label = tr$label, auroc = tr$auroc, lo = tr$auroc_lo, hi = tr$auroc_hi),
             data.frame(site = SITE_LEV[2], label = ts$label, auroc = ts$auroc, lo = ts$auroc_lo, hi = ts$auroc_hi),
             data.frame(site = SITE_LEV[3], label = es$label, auroc = es$auroc, lo = es$auroc_lo, hi = es$auroc_hi))
main_arms <- c("llr_sum", "llr_meas", "xgb_l", "xgb_feat", "xgb_raw")
lad <- lad[lad$label %in% main_arms, ]
lad$arm <- factor(ARM_LAB[lad$label], rev(ARM_LAB[main_arms])); lad$site <- factor(lad$site, SITE_LEV)
p3 <- ggplot(lad, aes(auroc, arm, colour = site)) +
  geom_line(aes(group = arm), colour = GRID, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = SITE_COL) +
  labs(x = "AUROC (patient-resampled bootstrap 95% CI)", y = NULL) + thm
save_fig(p3, "F3_transport_ladder", 6.5, 3.2)

# ==== F4: severity comparators =================================================
sev_tab <- function(dir, site, cell, suffix, file) {
  s <- rd(dir, file); s$label <- sub(suffix, "", s$label)
  data.frame(site = site, cell = cell, label = s$label, auroc = s$auroc, lo = s$auroc_lo, hi = s$auroc_hi, n = s$n)
}
sv <- rbind(sev_tab(E, "eICU", "A: baseline-restricted", "_sev$", "severity_score_summary"),
            sev_tab(E, "eICU", "B: symmetric", "_symsev$", "severity_score_summary_symmetric"),
            sev_tab(T, "MIMIC-IV test", "A: baseline-restricted", "_sev$", "severity_score_summary"),
            sev_tab(T, "MIMIC-IV test", "B: symmetric", "_symsev$", "severity_score_summary_symmetric"))
SEV_LAB <- c(llr_meas = "WoE measurement", apache2_aps = "APACHE II APS", apache2_total = "APACHE II total",
             llr_sum = "WoE paired", sofa = "SOFA")
sv <- sv[sv$label %in% names(SEV_LAB), ]
sv$score <- factor(SEV_LAB[sv$label], rev(c("WoE measurement", "APACHE II APS", "APACHE II total", "WoE paired", "SOFA")))
sv$family <- ifelse(grepl("^llr", sv$label), "Proposed", "Classical (recalibrated on MIMIC)")
sv$panel <- paste(sv$site, sv$cell, sep = " · ")
p4 <- ggplot(sv, aes(auroc, score, colour = family)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(values = c("Proposed" = "#2a78d6", "Classical (recalibrated on MIMIC)" = "#e34948")) +
  labs(x = "AUROC (patient-resampled bootstrap 95% CI)", y = NULL) + thm
save_fig(p4, "F4_severity_comparators", 7.2, 4.6)

# ==== F5: eigenspectrum ========================================================
eg <- rbind(cbind(matrix = "Unmeasured as zero", rd(I, "eigen_full_zero")),
            cbind(matrix = "Unmeasured as missing", rd(I, "eigen_full_na")))
eg$matrix <- factor(eg$matrix, c("Unmeasured as zero", "Unmeasured as missing"))
p5 <- ggplot(eg, aes(component, prop, colour = matrix)) +
  geom_hline(yintercept = 1 / 19, colour = INK2, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 19)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Principal component of the 19-signal evidence correlation matrix", y = "Share of variance",
       subtitle = "Dashed line: 1/19, the share under no correlation") + thm
save_fig(p5, "F5_eigenspectrum", 6, 3.4)

# ==== F6: per-signal single-column AUROC, joint vs conditional ================
sf <- rd(I, "signal_auroc_full"); sc <- rd(I, "signal_auroc_cond")
ps <- rbind(data.frame(model = "WoE paired", signal = sf$signal, auroc = sf$auroc),
            data.frame(model = "WoE conditional", signal = sc$signal, auroc = sc$auroc))
ord <- sf$signal[order(sf$auroc)]
ps$signal <- factor(ps$signal, ord); ps$model <- factor(ps$model, c("WoE paired", "WoE conditional"))
p6 <- ggplot(ps, aes(auroc, signal, colour = model)) +
  geom_line(aes(group = signal), colour = GRID, linewidth = 1.2) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 0.5, colour = INK2, linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = ARM_COL[c("WoE paired", "WoE conditional")]) +
  labs(x = "Out-of-fold AUROC of the single signal's evidence column", y = NULL) + thm
save_fig(p6, "F6_signal_auroc", 6, 4.4)

# ==== F7: partial-effect curves with the support mask ==========================
# Anchored evidence functions (paper/make_anchored_curves.R, R/15_anchors.R):
# the QC grid and support flags, each term re-anchored at its reference value
# (a quiet day on the channel), exact bands from the frozen bundle.
cv <- readRDS("paper/figs/coefs/anchored_curves.rds")
ANCH <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
anchor_at <- function(key, term) ANCH$anchor[ANCH$key == key & ANCH$term == term][1]
# Paper vocabulary for the model terms (methods: Random Variables). Evidence is
# in nats; a smooth is centered, so its value is the term's contribution to the
# channel's weight of evidence relative to the term's average contribution.
TERM_LAB <- c(
  "s(pi_minus)"  = "Deviation propensity, low side",
  "s(pi_plus)"   = "Deviation propensity, high side",
  "s(q05_delta)" = "Frequency-adjusted tail extremity, low side",
  "s(q95_delta)" = "Frequency-adjusted tail extremity, high side",
  "s(trend)"     = "Trend",
  "s(invasive_vent__exposure_frac)"   = "Invasive ventilation, exposure fraction",
  "s(sedation_dexmed__exposure_frac)" = "Dexmedetomidine, exposure fraction")
X_LAB <- c(
  "s(pi_minus)"  = "Propensity to fall below the reference range",
  "s(pi_plus)"   = "Propensity to rise above the reference range",
  "s(invasive_vent__exposure_frac)"   = "Fraction of the 24-hour window ventilated",
  "s(sedation_dexmed__exposure_frac)" = "Fraction of the 24-hour window exposed",
  "s(vasopressor__exposure_frac)"     = "Fraction of the 24-hour window on a vasopressor")
UNIT <- c(spo2 = "percentage points", resp_rate = "breaths/min", mbp = "mmHg", heart_rate = "beats/min")
x_label <- function(signal, term) {
  if (term %in% names(X_LAB)) return(X_LAB[[term]])
  u <- UNIT[[signal]]
  if (term %in% c("s(q05_delta)", "s(q95_delta)")) return(sprintf("Extremity relative to expected (%s)", u))
  if (term == "s(trend)") return(sprintf("Slope of hourly medians (%s per 24 h)", u))
  term
}
Y_LAB <- "Contribution relative to reference (nats)"
MODEL_COL <- c("Paired model" = "#2a78d6", "Measurement model" = "#1baf7a")

# One panel: one or more curves of the same term, each with its analytic 95%
# band, unsupported grid shaded, dotted (or open points) outside support.
curve_panel <- function(specs, title, xlab, show_legend = FALSE) {
  d <- do.call(rbind, lapply(specs, function(s) {
    z <- cv[cv$key == s$key & cv$term == s$term, ]
    if (!nrow(z)) stop("no curve for ", s$key, " ", s$term)
    z <- z[order(z$grid_x), ]; z$lo <- z$fit - 1.96 * z$se; z$hi <- z$fit + 1.96 * z$se
    # one id per contiguous supported / unsupported run, so a line never joins
    # two separate unsupported stretches across a supported one
    z$run <- cumsum(c(1L, diff(as.integer(z$supported)) != 0L))
    # repeat the boundary point so solid and dotted runs meet without a gap
    nxt <- which(diff(z$run) != 0L)
    if (length(nxt) && z$kind[1] != "atoms") {
      b <- z[nxt + 1L, ]; b$run <- z$run[nxt]; b$supported <- z$supported[nxt]
      z <- rbind(z, b); z <- z[order(z$grid_x, z$run), ]
    }
    z$model_lab <- s$lab; z
  }))
  d$model_lab <- factor(d$model_lab, names(MODEL_COL))
  ref <- d[d$model_lab == levels(droplevels(d$model_lab))[1], ]
  ref <- ref[!duplicated(ref$grid_x), ]
  half <- if (nrow(ref) > 1L) stats::median(diff(ref$grid_x)) / 2 else 0
  r <- rle(!ref$supported); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  rects <- data.frame(xmin = ref$grid_x[starts[r$values]] - half, xmax = ref$grid_x[ends[r$values]] + half)
  atoms <- d$kind[1] == "atoms"
  g <- ggplot(d, aes(grid_x, fit, colour = model_lab, fill = model_lab))
  if (nrow(rects)) g <- g + geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                       inherit.aes = FALSE, fill = "#f0efec")
  x0 <- anchor_at(specs[[1]]$key, specs[[1]]$term)
  g <- g + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_vline(xintercept = x0, colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = model_lab), colour = NA, alpha = 0.16)
  g <- if (atoms) {
    g + geom_point(aes(shape = supported), size = 1.6) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none")
  } else {
    g + geom_line(aes(linetype = supported, group = interaction(model_lab, run)), linewidth = 0.7) +
      scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = "none")
  }
  g + scale_colour_manual(values = MODEL_COL, drop = TRUE, guide = if (show_legend) "legend" else "none") +
    scale_fill_manual(values = MODEL_COL, drop = TRUE, guide = if (show_legend) "legend" else "none") +
    labs(title = title, x = xlab, y = Y_LAB) + thm + theme(plot.title = element_text(size = 9))
}
one <- function(key, term, title) {
  sg <- sub("/.*", "", key)
  lab <- if (grepl("/meas$", key)) "Measurement model" else "Paired model"
  curve_panel(list(list(key = key, term = term, lab = lab)), title, x_label(sg, term))
}

p7 <- (one("spo2/meas", "s(q05_delta)", "Oxygen saturation, measurement model:\nfrequency-adjusted tail extremity, low side") |
       one("gcs_motor/full", "s(pi_minus)", "Glasgow motor, paired model:\ndeviation propensity, low side")) /
      (one("resp_rate/full", "s(pi_plus)", "Respiratory rate, paired model:\ndeviation propensity, high side") |
       one("gcs_motor/full", "s(sedation_dexmed__exposure_frac)", "Glasgow motor, paired model:\ndexmedetomidine exposure"))
save_fig(p7, "F7_partial_effects_support", 7.5, 5.8)

# ==== F15: one channel's evidence, term by term (respiratory rate, joint) =======
# The fitted parts of Equation eq:gam_resp_rate: five smooths from the QC curve
# table and the admission coefficient from the frozen bundle
# (paper/make_figure_coefs.R; model parameters, no patient rows).
# Propensities first: a tail extremity is read after the propensity on its side,
# since for a stay with no deviations it is the tail position within range.
rr_terms <- c("s(pi_minus)", "s(pi_plus)", "s(q95_delta)", "s(trend)", "s(invasive_vent__exposure_frac)")
rr_panels <- lapply(rr_terms, function(tm) one("resp_rate/full", tm, TERM_LAB[[tm]]))
pc <- read.csv("paper/figs/coefs/parametric_coefs.csv", stringsAsFactors = FALSE)
g_adm <- pc[pc$key == "resp_rate/full" & pc$term == "invasive_vent__present_at_admission", ]
if (nrow(g_adm) != 1L) stop("admission coefficient for resp_rate/full not found")
adm <- data.frame(x = factor(c("Not ventilated\nat admission", "Ventilated\nat admission"),
                             c("Not ventilated\nat admission", "Ventilated\nat admission")),
                  fit = c(0, g_adm$estimate), lo = c(0, g_adm$estimate - 1.96 * g_adm$se),
                  hi = c(0, g_adm$estimate + 1.96 * g_adm$se))
rr_panels[[6]] <- ggplot(adm, aes(x, fit)) + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, colour = "#2a78d6") +
  geom_point(size = 2, colour = "#2a78d6") +
  labs(title = "Invasive ventilation, admission indicator", x = NULL, y = Y_LAB) + thm +
  theme(plot.title = element_text(size = 9))
p15 <- wrap_plots(rr_panels, ncol = 3)
save_fig(p15, "F15_channel_terms_resp_rate", 9.5, 5.8)

# ==== F16: the same measurement term with and without the intervention block ===
# For the four channels whose evidence is most affected by conditioning, the
# largest-range measurement term under the measurement model and the joint model.
ov <- list(
  list(sg = "mbp",        tm = "s(pi_minus)",  title = "Mean blood pressure:\ndeviation propensity, low side"),
  list(sg = "heart_rate", tm = "s(pi_plus)",   title = "Heart rate:\ndeviation propensity, high side"),
  list(sg = "spo2",       tm = "s(q05_delta)", title = "Oxygen saturation:\nfrequency-adjusted tail extremity, low side"),
  list(sg = "resp_rate",  tm = "s(pi_plus)",   title = "Respiratory rate:\ndeviation propensity, high side"))
ov_panels <- lapply(seq_along(ov), function(i) {
  o <- ov[[i]]
  curve_panel(list(list(key = paste0(o$sg, "/meas"), term = o$tm, lab = "Measurement model"),
                   list(key = paste0(o$sg, "/full"), term = o$tm, lab = "Paired model")),
              o$title, x_label(o$sg, o$tm), show_legend = (i == 1))
})
p16 <- wrap_plots(ov_panels, ncol = 2, guides = "collect") & theme(legend.position = "bottom")
save_fig(p16, "F16_measurement_terms_meas_vs_joint", 7.5, 6.2)

# ==== F19: mean blood pressure, the measurement term and the vasopressor term ====
# Results 4.5: the low-side deviation propensity under the measurement and the
# paired model, beside the paired model's vasopressor exposure term.
p19 <- (curve_panel(list(list(key = "mbp/meas", term = "s(pi_minus)", lab = "Measurement model"),
                         list(key = "mbp/full", term = "s(pi_minus)", lab = "Paired model")),
                    "Deviation propensity, low side", x_label("mbp", "s(pi_minus)"), show_legend = TRUE) |
        one("mbp/full", "s(vasopressor__exposure_frac)", "Vasopressor, exposure fraction (paired model)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p19, "F19_mbp_measurement_and_vasopressor", 7.5, 3.5)

# ==== F9: attribution stability ladder, both sites ===============================
# MIMIC-IV: the internal replicate store; eICU: whole-bag refits on the same bags,
# applied at eICU (tests/attr_external_bags.R, tests/attr_metrics.R --site eicu).
lab <- function(m, contrast) {
  base <- ARM_LAB[[m]]
  if (contrast == "L3")  paste0(base, ": resampled")
  else if (contrast == "L2") paste0(base, ": re-seeded only")
  else if (contrast == "L3T") paste0(base, ": resampled and re-seeded")
  else NA
}
lev <- c("WoE measurement: resampled", "WoE paired: resampled",
         "SHAP (XGBoost constructed): re-seeded only", "SHAP (XGBoost constructed): resampled",
         "SHAP (XGBoost constructed): resampled and re-seeded")
THR_LAB <- c(noise_calibrated_q0.5 = "At the noise threshold", noise_calibrated_q0.9 = "At the 90th-percentile threshold")
KLAB <- function(k) factor(paste0("Top-", k, " set differs"), c("Top-1 set differs", "Top-3 set differs"))
dd <- do.call(rbind, lapply(names(ATTR_SITE), function(st) {
  d <- rd(ATTR_SITE[[st]], "disagreement_distributions")
  d <- d[d$agg == "signal" & d$delta == 0 & d$bag_population == "common", ]
  keep <- (d$contrast == "L3" & d$held == "seed=0" & d$method_a %in% c("llr_full", "llr_meas", "shap_xgb_feat")) |
          (d$contrast == "L2" & d$held == "resampled") | (d$contrast == "L3T")
  d <- d[keep, ]
  if (anyDuplicated(paste(d$method_a, d$contrast, d$k))) stop("F9: duplicate rows at ", st)
  data.frame(site = ATTR_SITE_LAB[[st]], label = mapply(lab, d$method_a, d$contrast), method = ARM_LAB[d$method_a],
             k = KLAB(d$k), min = d$min, p05 = d$p05, median = d$median, p95 = d$p95, max = d$max)
}))
if (nrow(dd) != 2 * 5 * 2) stop("F9: expected 20 rows, found ", nrow(dd))
dd$label <- factor(dd$label, rev(lev)); dd$site <- factor(dd$site, ATTR_SITE_LAB)
# Thresholded readings: the median top-k disagreement at the method's noise
# threshold and at the 90th percentile of its pooled differences, from the
# 16-replicate magnitude table, on the resampled rows.
thr <- do.call(rbind, lapply(names(ATTR_SITE), function(st) {
  mg <- rd(ATTR_SITE[[st]], "magnitude_metric_distributions")
  mg <- mg[mg$agg == "signal" & mg$contrast == "L3" & mg$held == "seed=0" & mg$metric == "disagree" &
           mg$method_a %in% c("llr_full", "llr_meas", "shap_xgb_feat") &
           mg$delta_kind %in% names(THR_LAB), ]
  if (nrow(mg) != 3 * 2 * 2) stop("F9: expected 12 thresholded rows at ", st, ", found ", nrow(mg))
  data.frame(site = ATTR_SITE_LAB[[st]], label = paste0(ARM_LAB[mg$method_a], ": resampled"),
             method = ARM_LAB[mg$method_a], k = KLAB(mg$k), reading = factor(THR_LAB[mg$delta_kind], THR_LAB),
             median = mg$median)
}))
thr$label <- factor(thr$label, rev(lev)); thr$site <- factor(thr$site, ATTR_SITE_LAB)
p9 <- ggplot(dd, aes(y = label, colour = method)) +
  geom_errorbarh(aes(xmin = min, xmax = max), height = 0, linewidth = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = p05, xmax = p95), height = 0, linewidth = 1.4) +
  geom_point(aes(x = median), size = 2.4) +
  geom_point(data = thr, aes(x = median, shape = reading), size = 2.3, fill = "white", stroke = 0.8) +
  facet_grid(site ~ k, scales = "free_x") +
  scale_colour_manual(values = ARM_COL, guide = "none") +
  scale_shape_manual(values = c(21, 23), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Share of patients whose leading set differs between replicates",
       y = NULL, subtitle = paste0("Filled dot: median, no threshold. Thick: 5th–95th percentile. Thin: range.\n",
                                   "Open markers: median with a threshold (16 replicates per method).")) + thm +
  theme(strip.text.y = element_text(angle = 0))
save_fig(p9, "F9_attribution_stability", 7.5, 5.6)

# ==== F10: leader collapse by threshold, both sites (bag distributions) ==========
pairs <- list(c("llr_meas", "llr_full"), c("llr_full", "shap_xgb_feat"))
dir_lab <- function(p, ab) if (ab) sprintf("%s leader in %s", ARM_LAB[[p[1]]], ARM_LAB[[p[2]]]) else sprintf("%s leader in %s", ARM_LAB[[p[2]]], ARM_LAB[[p[1]]])
S_SITE <- c(mimic = "MIMIC-IV train", eicu = "eICU")
pm_l4 <- lapply(ATTR_SITE, function(dir) {
  pm <- rd(dir, "pairwise_metric_distributions")
  pm[pm$agg == "signal" & pm$contrast == "L4" & pm$bag_population == "common", ]
})
lc <- do.call(rbind, lapply(names(ATTR_SITE), function(st) do.call(rbind, lapply(pairs, function(p)
  do.call(rbind, lapply(c(TRUE, FALSE), function(ab) {
    d <- if (ab) "ab" else "ba"
    do.call(rbind, lapply(c("0.01", "0.02", "0.05"), function(t) {
      r <- one_row(pm_l4[[st]][pm_l4[[st]]$method_a == p[1] & pm_l4[[st]]$method_b == p[2] &
                                 pm_l4[[st]]$metric == paste0("collapse_", d, "_", t), ], "F10")
      data.frame(site = S_SITE[[st]], direction = dir_lab(p, ab), threshold = paste0(as.numeric(t) * 100, "%"),
                 collapse = r$median, lo = r$p05, hi = r$p95)
    }))
  }))))))
lc$site <- factor(lc$site, S_SITE); lc$direction <- factor(lc$direction, unique(lc$direction))
p10 <- ggplot(lc, aes(threshold, collapse, colour = site, group = site)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.5) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~direction, ncol = 2) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Threshold: leader of one model carries less than this share of the other model's total evidence",
       y = "Share of patients (leader collapse)",
       subtitle = "Median and 5th–95th percentile over the 38 bags at each site") + thm
save_fig(p10, "F10_leader_collapse", 7.5, 4.6)

# ==== F12: rank displacement of one model's leader in the other's ranking ====
# (Superseded by the redrawn F12 below; not in the paper. Kept reproducible.)
pm <- pm_l4$mimic
S_M <- "MIMIC-IV train (median, 5th-95th pct over 38 bags)"; S_E <- "eICU (frozen bundle, single application)"
rd_e <- rd(X, "leader_collapse_distributions")
rd_e <- rd_e[rd_e$agg == "signal" & rd_e$metric == "rank_displacement" & rd_e$population == "all_scored", ]
disp <- do.call(rbind, lapply(pairs, function(p) do.call(rbind, lapply(c(TRUE, FALSE), function(ab) {
  d <- if (ab) "ab" else "ba"
  m <- do.call(rbind, lapply(c("median", "p90", "p95"), function(q) {
    r <- pm[pm$method_a == p[1] & pm$method_b == p[2] & pm$metric == paste0("disp_", d, "_", q), ]
    data.frame(site = S_M, direction = dir_lab(p, ab), quantile = q, disp = r$median, lo = r$p05, hi = r$p95)
  }))
  r <- rd_e[rd_e$method_a == p[1] & rd_e$method_b == p[2] & rd_e$direction == ifelse(ab, "a_leader_in_b", "b_leader_in_a"), ][1, ]
  e <- data.frame(site = S_E, direction = dir_lab(p, ab), quantile = c("median", "p90", "p95"),
                  disp = c(r$median, r$p90, r$p95), lo = NA, hi = NA)
  rbind(m, e)
}))))
disp$quantile <- factor(c(median = "Median patient", p90 = "90th percentile", p95 = "95th percentile")[disp$quantile],
                        c("Median patient", "90th percentile", "95th percentile"))
disp$site <- factor(disp$site, c(S_M, S_E)); disp$direction <- factor(disp$direction, unique(disp$direction))
p12 <- ggplot(disp, aes(quantile, disp, fill = site)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.5, na.rm = TRUE, position = position_dodge(width = 0.7)) +
  facet_wrap(~direction, ncol = 2) +
  scale_fill_manual(values = c("#2a78d6", "#eb6834")) +
  scale_y_continuous(breaks = seq(0, 18, 3), limits = c(0, 18), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Quantile over patients", y = "Rank of the leader in the other model (0 = same leader)") + thm
save_fig(p12, "F12_rank_displacement", 7.5, 6)

# ==== F13: level-3 rank displacement (one method, two resamples) ============
lc_m <- rd(A, "leader_collapse_distributions")
l3 <- lc_m[lc_m$site == "mimic" & lc_m$agg == "signal" & lc_m$contrast == "L3" & lc_m$held == "seed=0" &
           lc_m$metric == "rank_displacement" & lc_m$bag_population == "common" & lc_m$direction == "a_leader_in_b" &
           lc_m$method_a %in% c("llr_meas", "llr_full", "llr_cond", "shap_xgb_feat"), ]
l3 <- do.call(rbind, lapply(seq_len(nrow(l3)), function(i) data.frame(
  method = ARM_LAB[[l3$method_a[i]]],
  population = ifelse(l3$population[i] == "all_scored", "All patients", "Patients whose leader changed"),
  quantile = c("Median patient", "90th percentile", "95th percentile"),
  disp = c(l3$median[i], l3$p90[i], l3$p95[i]))))
l3$method <- factor(l3$method, ARM_LAB[c("llr_meas", "llr_full", "llr_cond", "shap_xgb_feat")])
l3$quantile <- factor(l3$quantile, c("Median patient", "90th percentile", "95th percentile"))
l3$population <- factor(l3$population, c("All patients", "Patients whose leader changed"))
p13 <- ggplot(l3, aes(quantile, disp, fill = method)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  facet_wrap(~population, ncol = 2) +
  scale_fill_manual(values = ARM_COL) +
  scale_y_continuous(breaks = 0:6, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Quantile over patients", y = "Rank of the leader in the other replicate",
       subtitle = "Two refits of the same method on two shared bags; pooled over 150 bag pairs") + thm
save_fig(p13, "F13_rank_displacement_L3", 7.5, 3.8)

# ==== Revision 2026-10-01: figures for the README and the paper ================
# F18 discrimination and transport (new), F3b AUPRC transport ladder (new), and the
# redrawn F4 (severity bars), F6 (joint vs measurement-only), F12 (joint vs SHAP)
# and F13 (changed leaders only). The earlier F4/F6/F12/F13 blocks above are kept
# so their files stay reproducible; the paper now includes the files written here.
ARMS6 <- c("llr_sum", "llr_meas", "xgb_l", "xgb_feat", "xgb_raw")   # conditional arm dropped 2026-10-03

# ---- shared: six arms, AUROC and AUPRC with intervals, three sites ------------
SITE3 <- c("MIMIC-IV\ntrain (out-of-fold)", "MIMIC-IV\ntest", "eICU\n(pooled)")
ev_of <- function(d, site) data.frame(site = site, label = d$label, n = d$n, n_events = d$n_events,
  auroc = d$auroc, auroc_lo = d$auroc_lo, auroc_hi = d$auroc_hi,
  auprc = d$auprc, auprc_lo = d$auprc_lo, auprc_hi = d$auprc_hi, stringsAsFactors = FALSE)
tr6 <- rd(I, "train_ref_arms")
ts6 <- rd(T, "score_summary"); ts6 <- ts6[ts6$binning == "frozen", ]; ts6$label <- sub("_frozen$", "", ts6$label)
es6 <- rd(E, "score_summary"); es6 <- es6[es6$binning == "frozen", ]; es6$label <- sub("_frozen$", "", es6$label)
dt <- rbind(ev_of(tr6, SITE3[1]), ev_of(ts6, SITE3[2]), ev_of(es6, SITE3[3]))
dt <- dt[dt$label %in% ARMS6, ]

# ==== F18: discrimination and transport =========================================
F18_COL <- c("WoE paired" = "#1d4f9c", "WoE measurement" = "#6f9fdc",
             "XGBoost raw" = "#b5401a", "XGBoost constructed" = "#e2723f", "XGBoost stacked" = "#f0a982")
HEAVY <- c("WoE paired", "XGBoost raw")
long <- rbind(
  data.frame(dt[, c("site", "label")], metric = "AUROC", v = dt$auroc, lo = dt$auroc_lo, hi = dt$auroc_hi),
  data.frame(dt[, c("site", "label")], metric = "AUPRC", v = dt$auprc, lo = dt$auprc_lo, hi = dt$auprc_hi))
long$arm <- factor(ARM_LAB[long$label], names(F18_COL))
long$x <- match(long$site, SITE3)
# manual dodge (the spacing position_dodge(width = 0.32) gives five arms), so the
# value labels can point at the dodged points
long$xd <- long$x + (as.integer(long$arm) - (nlevels(long$arm) + 1) / 2) * 0.32 / nlevels(long$arm)
long$heavy <- factor(long$arm %in% HEAVY, c(FALSE, TRUE))
rate_of <- function(s) { r <- dt[dt$site == s, ][1, ]; r$n_events / r$n }
chance_prc <- data.frame(x = 1:3, v = vapply(SITE3, rate_of, numeric(1)))
slope_panel <- function(m, ttl, ylab, chance = NULL) {
  d <- long[long$metric == m, ]
  ends <- d[d$x == 3 & d$heavy == TRUE, ]
  mids <- d[d$x < 3 & d$heavy == TRUE, ]          # MIMIC-IV out-of-fold and test values
  p <- ggplot(d, aes(xd, v, colour = arm, group = arm))
  if (!is.null(chance)) p <- p +
    geom_segment(data = chance, aes(x = x - 0.28, xend = x + 0.28, y = v, yend = v), inherit.aes = FALSE,
                 colour = INK2, linetype = "dashed", linewidth = 0.4)
  p + geom_line(aes(linewidth = heavy)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45) +
    geom_point(size = 2.2) +
    ggrepel::geom_text_repel(data = ends, aes(label = sprintf("%s  %.3f", arm, v)),
                             nudge_x = 0.5, direction = "y", hjust = 0, size = 3, colour = INK,
                             segment.colour = INK2, segment.size = 0.3, min.segment.length = 0, seed = 1) +
    ggrepel::geom_label_repel(data = mids, aes(label = sprintf("%.3f", v)), fill = "white", label.size = 0,
                              label.padding = 0.1, point.padding = 0.5, nudge_x = -0.32, direction = "y", hjust = 1, size = 2.8, colour = INK,
                             segment.colour = INK2, segment.size = 0.3, min.segment.length = 0, seed = 1) +
    scale_colour_manual(values = F18_COL, drop = FALSE) +
    scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.35), guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = SITE3, limits = c(0.45, 4.1), expand = c(0, 0)) +
    labs(title = ttl, x = NULL, y = ylab) + thm + theme(panel.grid.major.x = element_blank())
}
pA <- slope_panel("AUROC", "A  Discrimination (AUROC)", "AUROC (bootstrap 95% CI)")
pB <- slope_panel("AUPRC", "B  Precision-recall (AUPRC)", "AUPRC (bootstrap 95% CI)") +
  labs(subtitle = sprintf("Chance level is the mortality: %.1f%% in MIMIC-IV, %.1f%% at eICU",
                          100 * chance_prc$v[1], 100 * chance_prc$v[3]))
ret <- do.call(rbind, lapply(ARMS6, function(a) {
  t1 <- dt[dt$site == SITE3[1] & dt$label == a, ]; e1 <- dt[dt$site == SITE3[3] & dt$label == a, ]
  pt <- t1$n_events / t1$n; pe <- e1$n_events / e1$n
  data.frame(arm = ARM_LAB[[a]], metric = c("AUROC", "AUPRC"),
             r = c((e1$auroc - 0.5) / (t1$auroc - 0.5), ((e1$auprc - pe) / (1 - pe)) / ((t1$auprc - pt) / (1 - pt))))
}))
ord_ret <- ret$arm[ret$metric == "AUPRC"][order(ret$r[ret$metric == "AUPRC"])]
ret$arm <- factor(ret$arm, ord_ret); ret$metric <- factor(ret$metric, c("AUPRC", "AUROC"))
pC <- ggplot(ret, aes(r, arm, fill = arm)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * r)), hjust = -0.15, size = 3, colour = INK) +
  facet_wrap(~metric, nrow = 1) +
  scale_fill_manual(values = F18_COL, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(title = "C  Above-chance discrimination retained at eICU", x = NULL, y = NULL) + thm +
  theme(panel.grid.major.y = element_blank())
p18 <- (pA | pB) / pC + plot_layout(heights = c(1.25, 0.8), guides = "collect") &
  theme(legend.position = "bottom") & guides(colour = guide_legend(nrow = 1))
save_fig(p18, "F18_discrimination_transport", 10, 8)

# ==== F3b: transport ladder in AUPRC =============================================
SL <- setNames(SITE_LEV, SITE3)
lad_prc <- data.frame(site = factor(SL[dt$site], SITE_LEV), arm = factor(ARM_LAB[dt$label], rev(ARM_LAB[ARMS6])),
                      v = dt$auprc, lo = dt$auprc_lo, hi = dt$auprc_hi)
p3b <- ggplot(lad_prc, aes(v, arm, colour = site)) +
  geom_line(aes(group = arm), colour = GRID, linewidth = 1.2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = SITE_COL) +
  labs(x = "AUPRC (patient-resampled bootstrap 95% CI)", y = NULL) + thm
save_fig(p3b, "F3_transport_ladder_auprc", 6.5, 3.2)

# ==== F4 (redrawn): severity comparison as bars, baseline-restricted cell ========
SEV5 <- c(llr_sum = "WoE paired", llr_meas = "WoE measurement", apache2_aps = "APACHE II APS",
          apache2_total = "APACHE II total", sofa = "SOFA")
sev_ab <- function(dir, site) {
  s <- rd(dir, "severity_score_summary"); s$label <- sub("_sev$", "", s$label); s <- s[s$label %in% names(SEV5), ]
  rbind(data.frame(site = site, metric = "AUPRC", label = s$label, v = s$auprc, lo = s$auprc_lo, hi = s$auprc_hi,
                   chance = s$n_events / s$n),
        data.frame(site = site, metric = "AUROC", label = s$label, v = s$auroc, lo = s$auroc_lo, hi = s$auroc_hi,
                   chance = 0.5))
}
sb <- rbind(sev_ab(T, "MIMIC-IV test"), sev_ab(E, "eICU"))
sb$score <- SEV5[sb$label]
sb$family <- ifelse(grepl("^llr", sb$label), "Proposed", "Classical (recalibrated on MIMIC-IV)")
sb$panel <- factor(paste(sb$site, sb$metric, sep = " · "),
                   c("MIMIC-IV test · AUPRC", "eICU · AUPRC", "MIMIC-IV test · AUROC", "eICU · AUROC"))
sb$key <- paste(sb$panel, sb$score, sep = "__")
sb$key <- factor(sb$key, sb$key[order(sb$panel, sb$v)])
ch4 <- unique(sb[, c("panel", "chance")]); ch4 <- ch4[!duplicated(ch4$panel), ]
p4b <- ggplot(sb, aes(v, key, fill = family)) +
  geom_col(width = 0.64) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.4, colour = INK2) +
  geom_text(aes(x = hi, label = sprintf("%.3f", v)), hjust = -0.2, size = 2.9, colour = INK) +
  geom_vline(data = ch4, aes(xintercept = chance), colour = INK2, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~panel, ncol = 2, scales = "free") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c("Proposed" = "#2a78d6", "Classical (recalibrated on MIMIC-IV)" = "#e34948")) +
  labs(x = "Patient-resampled bootstrap 95% CI; dashed line: chance level", y = NULL) + thm +
  theme(panel.grid.major.y = element_blank())
save_fig(p4b, "F4_severity_bars", 8, 5.4)

# ==== F6 (redrawn): per-channel AUROC, joint vs measurement-only ================
sfj <- rd(I, "signal_auroc_full")
sfm <- read.csv("paper/figs/coefs/signal_auroc_meas_oof.csv", stringsAsFactors = FALSE)   # make_signal_auroc_meas.R
ps2 <- rbind(data.frame(model = "WoE paired", signal = sfj$signal, auroc = sfj$auroc),
             data.frame(model = "WoE measurement", signal = sfm$signal, auroc = sfm$auroc))
# display names as in the atlas (make_atlas.R)
CH_LAB <- c(mbp = "Mean blood pressure", heart_rate = "Heart rate", spo2 = "Oxygen saturation",
            resp_rate = "Respiratory rate", glucose = "Glucose", gcs_motor = "Glasgow motor",
            gcs_eyes = "Glasgow eyes", gcs_verbal = "Glasgow verbal", urine_output_rate = "Urine output rate",
            creatinine = "Creatinine", platelet = "Platelets", hemoglobin = "Hemoglobin",
            temperature = "Temperature", sodium = "Sodium", bicarbonate = "Bicarbonate", bun = "Urea nitrogen",
            wbc = "White cell count", lactate = "Lactate", bilirubin_total = "Bilirubin")
if (!all(ps2$signal %in% names(CH_LAB))) stop("signal without a display name")
ps2$signal <- factor(CH_LAB[ps2$signal], CH_LAB[sfj$signal[order(sfj$auroc)]])
ps2$model <- factor(ps2$model, c("WoE paired", "WoE measurement"))
p6b <- ggplot(ps2, aes(auroc, signal, colour = model)) +
  geom_line(aes(group = signal), colour = GRID, linewidth = 1.2) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 0.5, colour = INK2, linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = ARM_COL[c("WoE paired", "WoE measurement")]) +
  labs(x = "Out-of-fold AUROC of the channel's evidence column alone", y = NULL) + thm
save_fig(p6b, "F6_signal_auroc_meas", 6, 4.4)

# ==== F12 (redrawn): leader displacement between the paired WoE and SHAP, both sites ==
pp <- c("llr_full", "shap_xgb_feat")
QLAB <- c(median = "Median patient", p90 = "90th percentile", p95 = "95th percentile")
num_lab <- function(x) formatC(x, format = "fg", digits = 3)
dsh <- do.call(rbind, lapply(names(ATTR_SITE), function(st) do.call(rbind, lapply(c(TRUE, FALSE), function(ab) {
  d <- if (ab) "ab" else "ba"
  lab <- if (ab) "WoE paired leader, ranked in SHAP" else "SHAP leader, ranked in WoE paired"
  do.call(rbind, lapply(names(QLAB), function(q) {
    r <- one_row(pm_l4[[st]][pm_l4[[st]]$method_a == pp[1] & pm_l4[[st]]$method_b == pp[2] &
                               pm_l4[[st]]$metric == paste0("disp_", d, "_", q), ], "F12")
    data.frame(site = S_SITE[[st]], direction = lab, quantile = q, disp = r$median, lo = r$p05, hi = r$p95)
  }))
}))))
dsh$quantile <- factor(QLAB[dsh$quantile], QLAB)
dsh$site <- factor(dsh$site, S_SITE); dsh$direction <- factor(dsh$direction, unique(dsh$direction))
p12b <- ggplot(dsh, aes(quantile, disp, fill = site)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.5,
                position = position_dodge(width = 0.7)) +
  geom_text(aes(y = disp / 2, label = num_lab(disp)), position = position_dodge(width = 0.7), size = 2.9,
            colour = "white", fontface = "bold") +
  facet_wrap(~direction, ncol = 2, strip.position = "bottom") +
  scale_fill_manual(values = c("#2a78d6", "#eb6834")) +
  scale_y_reverse(breaks = seq(0, 18, 3), expand = expansion(mult = c(0.04, 0.04))) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "Rank in the other model (0 = same leader)",
       subtitle = "Bars: median over the 38 bags at each site; whiskers: 5th–95th percentile") + thm +
  theme(strip.placement = "outside")
save_fig(p12b, "F12_rank_displacement_shap", 7.5, 4)

# ==== F13 (redrawn): leader displacement under resampling, changed leaders, both sites ==
c3 <- do.call(rbind, lapply(names(ATTR_SITE), function(st) {
  lm3 <- rd(ATTR_SITE[[st]], "leader_collapse_distributions")
  x <- lm3[lm3$site == st & lm3$agg == "signal" & lm3$contrast == "L3" & lm3$held == "seed=0" &
           lm3$metric == "rank_displacement" & lm3$bag_population == "common" & lm3$direction == "a_leader_in_b" &
           lm3$population == "leader_differs" & lm3$method_a %in% c("llr_meas", "llr_full", "shap_xgb_feat"), ]
  if (nrow(x) != 3L || anyDuplicated(x$method_a)) stop("F13: expected one row per method at ", st)
  np <- unique(x$n_pairs); if (length(np) != 1L) stop("F13: methods differ in pair count at ", st)
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) data.frame(
    site = sprintf("%s (pooled over %d bag pairs)", S_SITE[[st]], np),
    method = ARM_LAB[[x$method_a[i]]], quantile = QLAB, disp = c(x$median[i], x$p90[i], x$p95[i]))))
}))
c3$site <- factor(c3$site, unique(c3$site))
c3$method <- factor(c3$method, ARM_LAB[c("llr_meas", "llr_full", "shap_xgb_feat")])
c3$quantile <- factor(c3$quantile, QLAB)
p13b <- ggplot(c3, aes(quantile, disp, fill = method)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = num_lab(disp)), position = position_dodge(width = 0.78), vjust = 1.6, size = 2.9, colour = INK) +
  facet_wrap(~site, ncol = 2) +
  scale_fill_manual(values = ARM_COL) +
  scale_y_reverse(breaks = 0:6, expand = expansion(mult = c(0.2, 0.04))) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "Rank of the leader in the other replicate",
       subtitle = "Patients whose leading channel changed between two refits on two shared bags") + thm
save_fig(p13b, "F13_rank_displacement_L3_changed", 9.5, 3.8)

cat("done\n")
