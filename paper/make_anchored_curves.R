# paper/make_anchored_curves.R ---------------------------------------------------
# The anchored evidence functions that every term figure reads (F7, F15, F16,
# F17, the atlas). Same grid and support flags as the QC curve table
# (curves_1d.rds); values re-anchored at the reference anchors of
# R/15_anchors.R, with exact bands from the frozen bundle. Model parameters and
# aggregate support flags only; no patient rows, no window counts (hard rule 1).
#
#   Rscript paper/make_anchored_curves.R
#
# Writes paper/figs/coefs/anchored_curves.rds (+ .csv), anchors.csv and
# reference_evidence.csv. Stops unless the centered curves recomputed from the
# bundle reproduce the QC table, which proves both describe the same models.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(mgcv); library(qs2); library(yaml) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
rc <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config %||% "config/config.yml")
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg_local, strict = TRUE, verbose = FALSE)
cv <- readRDS("out/runs/gamqc_20260909T180107/tables/curves_1d.rds")
cv <- cv[cv$model %in% c("meas", "full", "intv"), ]
an <- read.csv("paper/figs/coefs/anchor_n.csv", stringsAsFactors = FALSE)
n_of <- function(sg) an$n_median[match(sg, an$signal)]

out <- list(); anc <- list(); ref <- list(); worst <- 0
for (key in unique(cv$key)) {
  sg <- sub("/.*", "", key); md <- sub(".*/", "", key)
  b <- bundle$models[[key]]; if (is.null(b)) stop("bundle is missing ", key)
  pri <- priors_for(bundle$priors, sg, "final", NA_integer_)
  for (s in b$smooth) {
    z <- cv[cv$key == key & cv$term == s$label, ]; z <- z[order(z$grid_x), ]
    if (!nrow(z)) stop("no QC grid for ", key, " ", s$label)
    idx <- s$first.para:s$last.para
    centered <- as.numeric(.smooth_rows(s, z$grid_x) %*% coef(b)[idx])
    worst <- max(worst, abs(centered - z$fit))
    x0 <- term_anchor(s$term[1], pri, n_of(sg))
    ac <- anchored_curve(b, s, z$grid_x, x0)
    sup_x0 <- if (x0 < min(z$grid_x) || x0 > max(z$grid_x)) FALSE else z$supported[which.min(abs(z$grid_x - x0))]
    out[[length(out) + 1]] <- data.frame(key = key, signal = sg, model = md, term = s$label, variable = s$term[1],
      kind = z$kind, grid_x = z$grid_x, fit = ac$fit, se = ac$se, supported = z$supported, fit_centered = centered)
    anc[[length(anc) + 1]] <- data.frame(key = key, signal = sg, model = md, term = s$label, variable = s$term[1],
      anchor = x0, n_anchor = n_of(sg), anchor_in_grid = x0 >= min(z$grid_x) && x0 <= max(z$grid_x), anchor_supported = sup_x0)
  }
  re <- reference_evidence(b, pri, n_of(sg))
  ref[[length(ref) + 1]] <- data.frame(key = key, signal = sg, model = md, reference_evidence = re$value, se = re$se)
}
if (worst > 1e-5) stop(sprintf("recomputed centered curves differ from the QC table by %.2e", worst))
cur <- do.call(rbind, out); anc <- do.call(rbind, anc); ref <- do.call(rbind, ref)
dir.create("paper/figs/coefs", showWarnings = FALSE, recursive = TRUE)
saveRDS(cur, "paper/figs/coefs/anchored_curves.rds")
write.csv(anc, "paper/figs/coefs/anchors.csv", row.names = FALSE)
write.csv(ref, "paper/figs/coefs/reference_evidence.csv", row.names = FALSE)
cat(sprintf("anchored %d terms of %d models; centered check max |diff| %.1e\n", nrow(anc), nrow(ref), worst))
bad <- anc[!anc$anchor_supported, ]
cat(sprintf("anchors outside training support: %d of %d\n", nrow(bad), nrow(anc)))
if (nrow(bad)) print(bad[, c("key", "term", "anchor", "anchor_in_grid")], row.names = FALSE)
