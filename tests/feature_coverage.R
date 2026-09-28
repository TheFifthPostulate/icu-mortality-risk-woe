# tests/feature_coverage.R ----------------------------------------------------
# COVERAGE OF THE DELIVERED FEATURES AT BOTH SITES. Non-fitting consumer.
#
# Per measurement: share of stays with the signal measured (n_obs > 0), mean
# and quartiles of n_obs among measured stays, share of measured stays with a
# defined trend. Per intervention: prevalence (ever_active) and mean exposure
# fraction or hours among exposed stays. Both sites side by side, from the
# loaded tables through the bundle's frozen design, so the eICU rows describe
# exactly the tables the external run scored.
#
# Aggregates only (hard rule 1). No row is printed.
#
#   Rscript tests/feature_coverage.R
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(arrow); library(yaml); library(qs2) })
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

xcfg   <- yaml::read_yaml("config/external.yml")
bpath  <- as.character(cfg_req(xcfg, "bundle"))
cfg_m  <- load_config("config/config.yml")
bundle <- load_bundle(bpath, cfg = cfg_m, strict = TRUE, verbose = FALSE)
cfg_e  <- bundle_cfg(bundle, paths = cfg_req(xcfg, "paths"))
tabs   <- list(mimic = load_tables(cfg_m$paths$mimiciv, cfg_m, site = "mimic", verbose = FALSE),
               eicu  = load_tables(cfg_e$paths, cfg_e, site = "eicu", verbose = FALSE))

sig_cov <- function(sf, site) {
  do.call(rbind, lapply(as.character(unlist(cfg_m$signals)), function(sg) {
    z <- sf[sf$signal == sg, ]
    m <- z$n_obs > 0
    q <- stats::quantile(z$n_obs[m], c(.25, .5, .75), names = FALSE)
    data.frame(site = site, signal = sg, n_stays = nrow(z),
               share_measured = round(mean(m), 4),
               mean_hours = round(mean(z$n_obs[m]), 2),
               q25_hours = q[1], median_hours = q[2], q75_hours = q[3],
               share_trend_defined = round(mean(z$dx_trend_defined[m] == 1), 4),
               share_masked_any = round(mean(z$dx_n_masked > 0), 4),
               stringsAsFactors = FALSE)
  }))
}
iv_cov <- function(ivf, site) {
  do.call(rbind, lapply(sort(unique(as.character(ivf$intervention))), function(iv) {
    z <- ivf[ivf$intervention == iv, ]
    ex <- z$ever_active == 1
    st <- unique(as.character(z$shape))
    data.frame(site = site, intervention = iv, shape = st[1], n_stays = nrow(z),
               prevalence = round(mean(ex), 4),
               mean_exposure_frac = if (st[1] == "state") round(mean(z$exposure_frac[ex]), 3) else NA_real_,
               share_exposure_one = if (st[1] == "state") round(mean(z$exposure_frac[ex] == 1), 3) else NA_real_,
               mean_hours = if (st[1] == "event") round(mean(z$n_hours[ex]), 2) else NA_real_,
               present_at_admission = round(mean(z$present_at_admission[ex] == 1), 3),
               stringsAsFactors = FALSE)
  }))
}
sc <- rbind(sig_cov(tabs$mimic$signal_features, "mimic"), sig_cov(tabs$eicu$signal_features, "eicu"))
ic <- rbind(iv_cov(tabs$mimic$intervention_features, "mimic"), iv_cov(tabs$eicu$intervention_features, "eicu"))

run <- new_run("featcov", cfg_m, note = "coverage of delivered features at both sites")
save_table(run, sc, "signal_coverage", subdir = "tables")
save_table(run, ic, "intervention_coverage", subdir = "tables")
cat("\nSIGNAL COVERAGE\n"); print(sc, row.names = FALSE)
cat("\nINTERVENTION COVERAGE\n"); print(ic, row.names = FALSE)
finalize_run(run)
