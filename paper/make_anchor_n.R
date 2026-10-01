# paper/make_anchor_n.R ------------------------------------------------------------
# The charting density at which the deviation-propensity anchor is set: for each
# channel, the median number of observed hours (n_obs, after masking) among
# MIMIC-IV TRAINING stays with the measurement present. Training split only, by
# the pipeline's own deterministic partition (assign_folds fits nothing); test
# stays and eICU do not enter. Output is one integer per channel (aggregate,
# hard rule 1), frozen in paper/figs/coefs/anchor_n.csv and read by
# R/15_anchors.R::term_anchor().
#
#   Rscript paper/make_anchor_n.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(arrow); library(yaml); library(qs2); library(mgcv) })
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
rc <- yaml::read_yaml("config/internal.yml")
cfg_local <- load_config(rc$config %||% "config/config.yml")
bundle <- load_bundle(rc$test_look$bundle, cfg = cfg_local, strict = TRUE, verbose = FALSE)
cfg  <- bundle_cfg(bundle, cfg_local$paths$mimiciv)
tabs <- load_tables(cfg$paths, cfg, site = "mimic", verbose = FALSE)
folds <- assign_folds(tabs$cohort, cfg)
train_ids <- as.character(folds[["stay_id"]][folds[["split"]] == "train"])

sf <- tabs$signal_features
res <- do.call(rbind, lapply(cfg$signals, function(sg) {
  sel <- sf[["signal"]] == sg & sf[["n_obs"]] > 0 & as.character(sf[["stay_id"]]) %in% train_ids
  nob <- sf[["n_obs"]][sel]
  data.frame(signal = sg, n_median = as.integer(round(stats::median(nob))), n_stays = length(nob))
}))
dir.create("paper/figs/coefs", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "paper/figs/coefs/anchor_n.csv", row.names = FALSE)
print(res, row.names = FALSE)
