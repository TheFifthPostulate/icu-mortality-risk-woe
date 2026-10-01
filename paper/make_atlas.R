# paper/make_atlas.R ------------------------------------------------------------
# Appendix D: the evidence-function atlas. One page per channel with every
# fitted term of its final evidence models, drawn from model-level aggregates
# only (hard rule 1): the QC curve table (fit, analytic SE, support flag on a
# grid), the QC curve summary (shrunk-out flag) and the parametric
# coefficients read from the frozen bundle (paper/make_figure_coefs.R). The
# window counts stored with the curves (n, events) are NOT drawn or exported:
# in sparse regions they are small cells.
#
# Each panel is one random variable of the channel:
#   measurement terms: measurement model (green) and joint model (blue);
#   intervention terms: intervention model (amber) and joint model (blue);
#   unpaired channels: one model, which is both the measurement and the joint.
# Panel order is the reading order: propensities, tail extremities, trend,
# then each paired intervention (exposure, exposure-adjusted intensity,
# admission indicator).
#
#   Rscript paper/make_atlas.R
#
# Writes paper/figs/atlas/atlas_<signal>.pdf and .png, and paper/appendix_atlas.tex.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

G   <- "out/runs/gamqc_20260909T180107/tables"
OUT <- "paper/figs/atlas"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
# Anchored evidence functions (paper/make_anchored_curves.R, R/15_anchors.R):
# each term relative to its reference value, a quiet day on the channel.
cv  <- readRDS("paper/figs/coefs/anchored_curves.rds")
cv  <- cv[, c("key", "signal", "model", "term", "kind", "grid_x", "fit", "se", "supported")]
ANCH <- read.csv("paper/figs/coefs/anchors.csv", stringsAsFactors = FALSE)
REFE <- read.csv("paper/figs/coefs/reference_evidence.csv", stringsAsFactors = FALSE)
anchor_at <- function(key, term) ANCH$anchor[ANCH$key == key & ANCH$term == term][1]
cs  <- read.csv(file.path(G, "curve_summary.csv"), stringsAsFactors = FALSE)
pc  <- read.csv("paper/figs/coefs/parametric_coefs.csv", stringsAsFactors = FALSE)

INK2 <- "#52514e"; GRID <- "#e6e5e1"
MODEL_COL <- c("Measurement model" = "#1baf7a", "Intervention model" = "#eda100", "Joint model" = "#2a78d6",
               "Measurement model (unpaired: also the joint)" = "#1baf7a")
thm <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2, size = 8),
        plot.title = element_text(size = 8.5, face = "bold"), legend.position = "bottom",
        legend.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA))

CH <- c(mbp = "Mean blood pressure", heart_rate = "Heart rate", spo2 = "Oxygen saturation",
        resp_rate = "Respiratory rate", glucose = "Glucose", gcs_motor = "Glasgow motor",
        gcs_eyes = "Glasgow eyes", gcs_verbal = "Glasgow verbal", urine_output_rate = "Urine output rate",
        creatinine = "Creatinine", platelet = "Platelets", hemoglobin = "Hemoglobin",
        temperature = "Temperature", sodium = "Sodium", bicarbonate = "Bicarbonate", bun = "Urea nitrogen",
        wbc = "White cell count", lactate = "Lactate", bilirubin_total = "Bilirubin")
UNIT <- c(mbp = "mmHg", heart_rate = "beats/min", spo2 = "percentage points", resp_rate = "breaths/min",
          glucose = "mg/dL", gcs_motor = "logit scale", gcs_eyes = "logit scale", gcs_verbal = "logit scale",
          urine_output_rate = "mL/kg/h", creatinine = "mg/dL", platelet = "10^9/L", hemoglobin = "g/dL",
          temperature = "°C", sodium = "mmol/L", bicarbonate = "mmol/L", bun = "mg/dL",
          wbc = "10^9/L", lactate = "mmol/L", bilirubin_total = "mg/dL")
IV <- c(vasopressor = "Vasopressor", inotrope = "Inotrope", invasive_vent = "Invasive ventilation",
        fio2 = "Inspired oxygen", insulin = "Insulin infusion", sedation_benzo = "Benzodiazepine",
        sedation_propofol = "Propofol", sedation_dexmed = "Dexmedetomidine", diuretic = "Diuretic",
        rrt = "Renal replacement", transfusion_prbc = "Red-cell transfusion",
        transfusion_platelet = "Platelet transfusion")

# --- vocabulary for one term ------------------------------------------------------
var_of <- function(term) sub("^s\\((.*)\\)$", "\\1", term)
describe <- function(sg, term) {
  v <- var_of(term); u <- UNIT[[sg]]
  if (grepl("__", v)) {
    iv <- IV[[sub("__.*", "", v)]]; what <- sub(".*__", "", v)
    return(switch(what,
      exposure_frac = list(title = paste0(iv, ", exposure fraction"), x = "Fraction of the 24-hour window exposed"),
      n_hours = list(title = paste0(iv, ", hours with an administration"), x = "Hours with an administration"),
      lambda = list(title = paste0(iv, ", exposure-adjusted intensity"), x = "Intensity relative to what exposure predicts"),
      present_at_admission = list(title = paste0(iv, ", admission indicator"), x = NULL),
      list(title = v, x = v)))
  }
  switch(v,
    pi_minus = list(title = "Deviation propensity, low side", x = "Propensity to fall below the reference range"),
    pi_plus  = list(title = "Deviation propensity, high side", x = "Propensity to rise above the reference range"),
    q05_delta = , value_min_delta = list(title = "Frequency-adjusted tail extremity, low side",
                                         x = sprintf("Extremity relative to expected (%s)", u)),
    q95_delta = , value_max_delta = list(title = "Frequency-adjusted tail extremity, high side",
                                         x = sprintf("Extremity relative to expected (%s)", u)),
    trend = list(title = "Trend", x = sprintf("Slope of hourly medians (%s per 24 h)",
                                             if (u == "logit scale") "score points" else u)),
    list(title = v, x = v))
}
rank_of <- function(term) {
  v <- var_of(term)
  if (startsWith(v, "pi_")) return(1)
  if (endsWith(v, "_delta")) return(2)
  if (v == "trend") return(3)
  4 + match(sub("__.*", "", v), names(IV)) / 100 +
    match(sub(".*__", "", v), c("exposure_frac", "n_hours", "lambda", "present_at_admission")) / 1000
}

# Two-line titles: break after the comma so no title is cut at the panel edge.
wrap_title <- function(x) sub(", ", ",\n", x, fixed = TRUE)

# --- one smooth panel, one or two models -------------------------------------------
smooth_panel <- function(sg, term, models, shrunk) {
  d <- do.call(rbind, lapply(names(models), function(lab) {
    z <- cv[cv$key == paste(sg, models[[lab]], sep = "/") & cv$term == term, ]
    if (!nrow(z)) return(NULL)
    z <- z[order(z$grid_x), ]; z$lo <- z$fit - 1.96 * z$se; z$hi <- z$fit + 1.96 * z$se
    z$run <- cumsum(c(1L, diff(as.integer(z$supported)) != 0L)); z$model_lab <- lab; z
  }))
  d$model_lab <- factor(d$model_lab, names(MODEL_COL))
  ref <- d[d$model_lab == levels(droplevels(d$model_lab))[length(levels(droplevels(d$model_lab)))], ]
  half <- if (nrow(ref) > 1L) stats::median(diff(ref$grid_x)) / 2 else 0
  r <- rle(!ref$supported); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  rects <- data.frame(xmin = ref$grid_x[starts[r$values]] - half, xmax = ref$grid_x[ends[r$values]] + half)
  lab <- describe(sg, term)
  g <- ggplot(d, aes(grid_x, fit, colour = model_lab, fill = model_lab))
  if (nrow(rects)) g <- g + geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                       inherit.aes = FALSE, fill = "#f0efec")
  g <- g + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_vline(xintercept = anchor_at(paste(sg, models[[length(models)]], sep = "/"), term),
               colour = INK2, linewidth = 0.3, linetype = "22") +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = model_lab), colour = NA, alpha = 0.14)
  g <- if (d$kind[1] == "atoms") g + geom_point(aes(shape = supported), size = 1.2) +
         scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none")
       else g + geom_line(aes(linetype = supported, group = interaction(model_lab, run)), linewidth = 0.6) +
         scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = "none")
  g + scale_colour_manual(values = MODEL_COL, guide = "none") + scale_fill_manual(values = MODEL_COL, guide = "none") +
    expand_limits(y = c(-0.25, 0.25)) +
    labs(title = wrap_title(paste0(lab$title, if (shrunk) " (shrunk)" else "")), x = lab$x,
         y = "Relative to reference (nats)") + thm
}
coef_panel <- function(sg, term, models) {
  d <- do.call(rbind, lapply(names(models), function(lab) {
    z <- pc[pc$key == paste(sg, models[[lab]], sep = "/") & pc$term == term, ]
    if (!nrow(z)) return(NULL)
    data.frame(model_lab = lab, fit = z$estimate, lo = z$estimate - 1.96 * z$se, hi = z$estimate + 1.96 * z$se)
  }))
  d$model_lab <- factor(d$model_lab, names(MODEL_COL))
  ggplot(d, aes(model_lab, fit, colour = model_lab)) + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) + geom_point(size = 1.8) +
    scale_colour_manual(values = MODEL_COL, guide = "none") + expand_limits(y = c(-0.25, 0.25)) +
    labs(title = wrap_title(describe(sg, term)$title), x = "At admission, against not", y = "Relative to reference (nats)") +
    thm + theme(axis.text.x = element_blank())
}

# --- one page per channel ---------------------------------------------------------
pages <- character(0)
for (sg in names(CH)) {
  keys  <- unique(cv$model[cv$signal == sg])
  paired <- "full" %in% keys
  jm    <- if (paired) "full" else "meas"
  terms <- unique(c(cv$term[cv$signal == sg & cv$model == jm], if (paired) cv$term[cv$signal == sg & cv$model == "intv"]))
  paras <- unique(pc$term[pc$key %in% paste(sg, c(jm, "intv"), sep = "/")])
  all_t <- c(terms, paras); all_t <- all_t[order(vapply(all_t, rank_of, numeric(1)))]
  panels <- lapply(all_t, function(t) {
    is_iv <- grepl("__", var_of(t))
    models <- if (!paired) list("Measurement model (unpaired: also the joint)" = "meas")
              else if (is_iv) list("Intervention model" = "intv", "Joint model" = "full")
              else list("Measurement model" = "meas", "Joint model" = "full")
    if (t %in% paras) return(coef_panel(sg, t, models))
    sh <- cs$shrunk_out[cs$key == paste(sg, jm, sep = "/") & cs$term == t]
    smooth_panel(sg, t, models, isTRUE(as.logical(sh)))
  })
  nr <- ceiling(length(panels) / 3)
  re <- REFE[REFE$signal == sg, ]; rv <- function(m) re$reference_evidence[re$model == m]
  key_txt <- if (paired) sprintf(paste0("Green: measurement model.  Amber: intervention model.  Blue: joint model.
",
                                        "Reference evidence (all terms at their reference): measurement %+.2f, intervention %+.2f, joint %+.2f nats."),
                                 rv("meas"), rv("intv"), rv("full"))
             else sprintf(paste0("Green: measurement model, which is also the joint model (no paired intervention).
",
                                 "Reference evidence (all terms at their reference): %+.2f nats."), rv("meas"))
  p <- wrap_plots(panels, ncol = 3) +
    plot_annotation(title = CH[[sg]], subtitle = key_txt,
                    theme = theme(plot.title = element_text(face = "bold", size = 12),
                                  plot.subtitle = element_text(colour = INK2, size = 9)))
  h <- min(9.2, 2.35 * nr + 0.9)
  f <- file.path(OUT, paste0("atlas_", sg))
  ragg::agg_png(paste0(f, ".png"), width = 8, height = h, units = "in", res = 200); print(p); dev.off()
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = 8, height = h); print(p); dev.off()
  pages <- c(pages, sg)
  cat(sprintf("%-18s %2d panels\n", sg, length(panels)))
}

# --- the appendix -------------------------------------------------------------------
tex <- c(
  "% paper/appendix_atlas.tex -- GENERATED by paper/make_atlas.R; do not edit by hand.",
  "\\clearpage",
  "\\setcounter{figure}{0}",
  "\\renewcommand{\\thefigure}{D\\arabic{figure}}",
  "\\section*{Appendix D. Atlas of the evidence functions}",
  "\\label{sm:atlas}",
  "",
  paste("This appendix draws every fitted term of the final evidence models, one channel per page. Each panel is",
        "one random variable of the channel. A measurement term is drawn under the measurement model and the joint",
        "model, and an intervention term under the intervention model and the joint model, so the effect of adding",
        "the other block can be read directly. A channel without a paired intervention has one model, which is both",
        "its measurement model and its joint model.",
        "

",
        "Each term is drawn relative to a reference value of its random variable, marked by a dotted vertical",
        "line (Section~\\ref{sec:evidence_models}). The reference values describe a quiet day on the channel: no",
        "deviation on either side at the channel's typical number of observed hours, a tail as extreme as the",
        "deviation count predicts, a flat trend, and no intervention. A point on a curve is therefore the evidence,",
        "in nats, that the value adds to or removes from the channel's weight of evidence compared with the reference",
        "value, and the band is the 95\\% interval of that difference. The weight of evidence of a stay whose",
        "variables all sit at their reference values is the reference evidence, given for each model at the top of",
        "its page. A stay's weight of evidence is the reference evidence plus the value of each term at the stay.",
        "

",
        "Shaded regions are outside empirical support, with fewer than 50 stays or 5 deaths in a 5\\% window, and",
        "are drawn dotted or as open points. An admission indicator is a single coefficient against not being",
        "present at admission, with its 95\\% interval. The panels follow the reading order of the card: the",
        "deviation propensities, then the frequency-adjusted tail extremity, which is read after the propensity on",
        "its side, then the trend, then the terms of each paired intervention. A term marked as shrunk was removed",
        "by the shrinkage basis in the joint model. The tail extremity of the Glasgow components is on the logit",
        "scale of the proportional-odds residual (Equation~\\ref{eq:midpit}), and their trend is in score points.",
        "Their reference value of zero falls between the few values that the ordinal residual takes, so for the",
        "motor and verbal components the curve at the reference is interpolated between supported values."),
  "")
for (sg in pages) {
  tex <- c(tex, "\\begin{figure}[p]", "\\centering",
           sprintf("\\includegraphics[width=\\linewidth,height=0.92\\textheight,keepaspectratio]{figs/atlas/atlas_%s.pdf}", sg),
           sprintf("\\caption{Evidence functions of %s.}", tolower(CH[[sg]])),
           sprintf("\\label{fig:atlas_%s}", sg), "\\end{figure}", "")
}
writeLines(tex, "paper/appendix_atlas.tex")
cat("wrote paper/appendix_atlas.tex with", length(pages), "figures\n")
