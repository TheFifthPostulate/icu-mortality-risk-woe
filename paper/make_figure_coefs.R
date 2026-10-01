# paper/make_figure_coefs.R ----------------------------------------------------
# Parametric coefficients of the final evidence models, for the channel
# decomposition figure. The QC run stores smooth curves only (curves_1d.rds);
# the admission indicator enters each model as a single parametric coefficient,
# which lives only in the frozen bundle. This script reads the bundle's final
# GAMs and writes their parametric table: model parameters, no patient rows
# (hard rule 1). Application-side reading of frozen parameters, no fitting.
#
#   Rscript paper/make_figure_coefs.R
#
# Writes paper/figs/coefs/parametric_coefs.csv.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(mgcv); library(qs2); library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

rc <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config %||% "config/config.yml")
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg_local, strict = TRUE, verbose = FALSE)

keys <- grep("/(meas|full|intv)$", names(bundle$models), value = TRUE)
rows <- lapply(keys, function(k) {
  b  <- bundle$models[[k]]
  V  <- if (!is.null(b$Vc)) b$Vc else b$Vp
  cf <- stats::coef(b)
  para <- names(cf)[!grepl("^s\\(|^te\\(|^ti\\(", names(cf))]
  para <- setdiff(para, "(Intercept)")
  if (!length(para)) return(NULL)
  data.frame(key = k, term = para, estimate = unname(cf[para]),
             se = sqrt(diag(V)[match(para, names(cf))]), stringsAsFactors = FALSE)
})
out <- do.call(rbind, rows)
dir.create("paper/figs/coefs", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "paper/figs/coefs/parametric_coefs.csv", row.names = FALSE)
cat("parametric terms written:", nrow(out), "from", length(keys), "final models\n")
