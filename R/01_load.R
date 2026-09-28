# R/01_load.R ----------------------------------------------------------------
# The adapter. Reads parquet, coerces types, sets factor levels, and performs
# the only data repairs in the project. Everything it does is declared here and
# logged; nothing is repaired silently.
#
# It does NOT construct paths (hard rule 9) — the caller passes them in from
# config. It does NOT branch on site (hard rule 5) — `site` is a constant
# argument that gets stamped, not a condition that changes behaviour.
#
# It does NOT print rows. Tables are 50k-980k rows and a single print would cost
# a large slice of the context window (hard rule 1). Logs are counts only.
# ----------------------------------------------------------------------------

# --- config -----------------------------------------------------------------

#' Read config.yml and pairing.csv, and check they agree with each other.
#'
#' This is a static-consistency check between two design files. It has nothing
#' to do with the data and runs before any parquet is opened, so a typo in a
#' whitelist fails in milliseconds rather than after a 980k-row read.
load_config <- function(config_path, pairing_path = NULL) {
  cfg <- yaml::read_yaml(config_path)
  cfg$pairing <- utils::read.csv(
    pairing_path %||% cfg$paths$pairing,
    colClasses = "character", na.strings = character()
  )

  p <- cfg$pairing
  stopifnot(identical(names(p), c("signal", "intervention", "signal_class", "excursion_side")))

  # Signals: pairing and whitelist must be the same set, no exceptions.
  if (!same_set(p$signal, cfg$signals)) {
    abort_values("pairing.csv signals disagree with config whitelist",
                 c(setdiff(p$signal, cfg$signals), setdiff(cfg$signals, p$signal)))
  }
  # Interventions: modelled must be a subset of extracted, and pairing may only
  # reference modelled ones.
  bad <- setdiff(cfg$interventions_modelled, cfg$interventions_extracted)
  if (length(bad)) abort_values("modelled interventions absent from extraction", bad)

  ivs <- setdiff(unique(p$intervention), "")
  bad <- setdiff(ivs, cfg$interventions_modelled)
  if (length(bad)) abort_values("pairing.csv references non-modelled interventions", bad)

  # Every extracted intervention needs a shape, and shapes must be legal.
  bad <- setdiff(cfg$interventions_extracted, names(cfg$intervention_shape))
  if (length(bad)) abort_values("interventions missing from intervention_shape map", bad)
  bad <- setdiff(unlist(cfg$intervention_shape), cfg$intervention_shapes)
  if (length(bad)) abort_values("illegal shape in intervention_shape map", bad)

  # signal_class must be single-valued per signal, and legal.
  by_sig <- tapply(p$signal_class, p$signal, function(z) length(unique(z)))
  if (any(by_sig != 1L)) abort_values("signal_class varies within a signal", names(by_sig)[by_sig != 1L])
  bad <- setdiff(p$signal_class, cfg$signal_classes)
  if (length(bad)) abort_values("illegal signal_class in pairing.csv", bad)

  # excursion_side must be low/high for paired rows and empty for unpaired ones.
  paired <- p$intervention != ""
  bad <- setdiff(p$excursion_side[paired], c("low", "high"))
  if (length(bad)) abort_values("illegal excursion_side on a paired row", bad)
  if (any(p$excursion_side[!paired] != "")) {
    abort_values("unpaired signal carries an excursion_side",
                 p$signal[!paired & p$excursion_side != ""])
  }

  # signal_tails: which pi coordinates enter, defaulting to both.
  if (length(cfg$signal_tails)) {
    bad <- setdiff(names(cfg$signal_tails), cfg$signals)
    if (length(bad)) abort_values("signal_tails names a signal outside the whitelist", bad)
    bad <- setdiff(unlist(cfg$signal_tails), c("low", "high"))
    if (length(bad)) abort_values("illegal tail in signal_tails (want low/high)", bad)
    empty <- names(cfg$signal_tails)[lengths(cfg$signal_tails) == 0L]
    if (length(empty)) {
      abort_values("signal_tails leaves a signal with no pi coordinate at all", empty)
    }
  }

  # level_terms: every class needs a parameterisation, and only two are legal.
  legal_scale <- c("quantile", "extreme")
  bad <- setdiff(unlist(cfg$level_terms_by_class), legal_scale)
  if (length(bad)) abort_values("illegal level_terms_by_class value (want quantile/extreme)", bad)
  bad <- setdiff(cfg$signal_classes, names(cfg$level_terms_by_class))
  if (length(bad)) abort_values("signal class missing from level_terms_by_class", bad)
  if (length(cfg$level_terms_override)) {
    bad <- setdiff(names(cfg$level_terms_override), cfg$signals)
    if (length(bad)) abort_values("level_terms_override names a signal outside the whitelist", bad)
    bad <- setdiff(unlist(cfg$level_terms_override), legal_scale)
    if (length(bad)) abort_values("illegal level_terms_override value (want quantile/extreme)", bad)
  }

  # The four resampling declarations. STATIC ONLY here: each must be a single
  # non-empty column name. Whether the column exists and is usable is a
  # question about DATA and belongs to check_resampling_cols() (R/03_folds.R),
  # which runs at both sites -- same division as smooth_k, declared here and
  # held to the counts there.
  for (k in c("split.group_by", "split.stratify_on",
              "fold_group_by", "fold_stratify_on")) {
    v <- if (startsWith(k, "split.")) cfg$split[[sub("^split[.]", "", k)]] else cfg[[k]]
    # Absent is legal: the fold-level keys default to the split-level ones, and
    # the split-level ones default in R/03_folds.R. Present-but-malformed is not.
    if (is.null(v)) next
    if (!is.character(v) || length(v) != 1L || !nzchar(v)) {
      abort_values(paste0("config `", k, "` must be a single non-empty column name"),
                   as.character(v))
    }
  }

  # trend_classes: a subset of the declared classes, and never empty.
  tc <- as.character(unlist(cfg$trend_classes %||% c("dense", "rate")))
  bad <- setdiff(tc, cfg$signal_classes)
  if (length(bad)) abort_values("trend_classes names a class outside signal_classes", bad)
  if (!length(tc)) stop("trend_classes is empty; no signal would carry a trend term", call. = FALSE)

  # smooth_k: per-(signal, covariate) basis dimension overrides.
  if (length(cfg$smooth_k)) {
    bad <- setdiff(names(cfg$smooth_k), cfg$signals)
    if (length(bad)) abort_values("smooth_k names a signal outside the whitelist", bad)
    kdef <- cfg$bam$k %||% 10
    for (sg in names(cfg$smooth_k)) {
      m <- cfg$smooth_k[[sg]]
      v <- unlist(m)
      if (!length(v) || !is.numeric(v)) {
        abort_values(paste0("smooth_k[", sg, "] values must be numeric"), names(m))
      }
      # An override may only ever REDUCE the basis. A larger k here would be a
      # silent per-signal model change wearing the name of a compatibility fix.
      bad <- names(m)[v >= kdef | v < 3]
      if (length(bad)) {
        abort_values(paste0("smooth_k[", sg, "] must be in [3, ", kdef,
                            "); an override may only reduce the basis"), bad)
      }
    }
  }

  cfg
}

#' Which magnitude parameterisation a signal uses: "quantile" or "extreme".
#'
#' Class default, overridden per signal. Design information, not a data-derived
#' fact, so R/05_formula.R stays data-free and the choice is frozen from MIMIC
#' rather than re-derived at eICU (hard rule 8) — which matters here more than
#' most, because the quantile-versus-extreme question is decided by measurement
#' DENSITY, and density is exactly what differs between the two sites.
level_scale_of <- function(signal, cfg) {
  v <- cfg$level_terms_override[[signal]]
  if (!is.null(v)) return(as.character(v))
  cls <- signal_class_of(signal, cfg)
  v <- cfg$level_terms_by_class[[cls]]
  if (is.null(v)) "quantile" else as.character(v)
}

#' The (low, high) level variable names for a signal, given its parameterisation.
level_vars_of <- function(signal, cfg) {
  if (identical(level_scale_of(signal, cfg), "extreme")) c("value_min", "value_max")
  else c("q05", "q95")
}

#' Does this signal carry a `trend` term? Class-gated, declared in config.
trend_enabled_for <- function(signal, cfg) {
  cls <- cfg$trend_classes %||% c("dense", "rate")
  signal_class_of(signal, cfg) %in% as.character(unlist(cls))
}

#' The basis dimension for one smooth covariate: the config override if there is
#' one, otherwise `bam.k`.
#'
#' Design information, not a data-derived fact — which is what keeps
#' R/05_formula.R data-free and freezes the basis from MIMIC instead of letting
#' eICU re-derive it (hard rule 8). R/04_features.R's check_smooth_k() holds the
#' declaration to the measured counts in both directions.
smooth_k_of <- function(signal, var, cfg) {
  v <- cfg$smooth_k[[signal]][[var]]
  if (is.null(v)) cfg$bam$k %||% 10 else as.integer(v)
}

#' Lookup: signal -> signal_class, and signal -> excursion_side.
#' Unpaired signals get excursion_side NA, which the formula builder reads as
#' "undefined" and answers with both tails (CLAUDE.md, frozen decisions).
signal_class_of <- function(signal, cfg) {
  m <- unique(cfg$pairing[, c("signal", "signal_class")])
  m$signal_class[match(signal, m$signal)]
}

excursion_side_of <- function(signal, cfg) {
  p <- cfg$pairing[cfg$pairing$signal == signal, ]
  s <- unique(p$excursion_side[p$excursion_side != ""])
  if (!length(s)) return(NA_character_)
  if (length(s) > 1L) {
    stop("signal '", signal, "' has conflicting excursion_side values: ",
         paste(s, collapse = ", "), call. = FALSE)
  }
  s
}

#' Which pi coordinates enter this signal's formula: "low", "high", or both.
#' Both unless config overrides it.
#'
#' NOT the same question as `excursion_side`, which selects the tail quantile.
#' pi_hat is a simplex and both coordinates are normally free;
#' config names only those signals where one is structurally pinned and would
#' therefore carry `n_obs` rather than physiology. See config/signal_tails.
occupiable_tails_of <- function(signal, cfg) {
  t <- cfg$signal_tails[[signal]]
  if (is.null(t)) c("low", "high") else as.character(unlist(t))
}

#' Interventions paired with a signal, in pairing.csv order. Empty for the 7.
interventions_of <- function(signal, cfg) {
  iv <- cfg$pairing$intervention[cfg$pairing$signal == signal]
  iv[iv != ""]
}

# --- tables -----------------------------------------------------------------

#' Load the modelled tables for one site.
#'
#' THREE are required. `ordering` is optional and audit-only: `o_flag` was
#' dropped from the design on 2026-08-25 and was the only modelled column this
#' table carried. What remains is an independent encoding of the signal x
#' intervention grid, which validator check 4 cross-checks against pairing.csv
#' when the table is present. Removing `paths$ordering` from config turns that
#' check off and changes nothing else — no model frame reads it.
#'
#' @param paths named list with cohort/signal_features/intervention_features,
#'   and optionally ordering
#' @param site  the constant stamped into every table. Not a branch.
load_tables <- function(paths, cfg, site, verbose = TRUE) {
  need <- c("cohort", "signal_features", "intervention_features")
  miss <- setdiff(need, names(paths))
  if (length(miss)) abort_values("missing path entries", miss)
  want <- c(need,
            if (!is.null(paths$ordering)) "ordering",
            if (!is.null(paths$severity)) "severity")
  missing_files <- vapply(paths[want], function(p) !file.exists(p), logical(1))
  if (any(missing_files)) abort_values("parquet file not found", unlist(paths[want])[missing_files])

  tabs <- list(
    cohort                = load_cohort(paths$cohort, cfg, site),
    signal_features       = load_signal_features(paths$signal_features, cfg, site),
    intervention_features = load_intervention_features(paths$intervention_features, cfg, site)
  )
  if ("ordering" %in% want) tabs$ordering <- load_ordering(paths$ordering, cfg, site)
  if ("severity" %in% want) tabs$severity <- load_severity(paths$severity, cfg, site)
  if (verbose) {
    for (nm in names(tabs)) {
      message(sprintf("  loaded %-22s %8d x %2d", nm, nrow(tabs[[nm]]), ncol(tabs[[nm]])))
    }
  }
  tabs
}

#' Stamp `site` if absent; verify it if already present.
#'
#' Three of the four tables ship without `site` and one ships with it. Stamping
#' a constant is not branching on site — the value is an argument, and every
#' table takes the same code path (hard rule 5, spec §5.1).
.stamp_site <- function(x, site, table_name) {
  if ("site" %in% names(x)) {
    got <- unique(as.character(x$site))
    if (length(got) != 1L || !identical(got, site)) {
      abort_values(paste0("`site` already present in ", table_name,
                          " but disagrees with the requested site '", site, "'"), got)
    }
    x$site <- factor(site, levels = site)
  } else {
    x$site <- factor(site, levels = site)
  }
  # Keep `site` first so the four tables present a consistent leading key.
  x[, c("site", setdiff(names(x), "site")), drop = FALSE]
}

load_cohort <- function(path, cfg, site) {
  x <- arrow::read_parquet(path)
  x <- as.data.frame(x)
  x <- .stamp_site(x, site, "cohort")

  x$gender <- factor(x$gender, levels = c("F", "M"))
  if (anyNA(x$gender)) stop("cohort: `gender` has values outside {F, M}", call. = FALSE)

  # Repair 1 of 1 in this table: the 30-300 weight guard. Applied in SQL as of
  # the 2026-08-22 re-extraction, so this is normally a no-op assertion. It
  # stays because urine_output_rate is mL/kg/hr — a 1 kg weight inflates a
  # paired signal roughly seventyfold, and that must never depend on which
  # extraction produced the file.
  lo <- cfg$guards$weight_kg_min %||% 30
  hi <- cfg$guards$weight_kg_max %||% 300
  bad <- !is.na(x$weight_kg) & (x$weight_kg < lo | x$weight_kg > hi)
  if (any(bad)) {
    message(sprintf("  cohort: weight_kg guard [%g, %g] set %d value(s) to NA", lo, hi, sum(bad)))
    x$weight_kg[bad] <- NA_real_
  }

  # Assertion, never a repair: the age cap belongs in SQL (LEAST(age, 90)).
  # Silently capping here would hide a re-extraction that forgot it, and an
  # uncapped MIMIC scored against a capped eICU is a transportability artifact.
  cap <- cfg$guards$age_max %||% 90
  if (any(x$age > cap, na.rm = TRUE)) {
    stop("cohort: `age` exceeds ", cap, " in ", sum(x$age > cap, na.rm = TRUE),
         " row(s). The cap belongs in SQL (LEAST(age, 90)); re-extract.", call. = FALSE)
  }
  x
}

load_signal_features <- function(path, cfg, site) {
  x <- as.data.frame(arrow::read_parquet(path))
  x <- .stamp_site(x, site, "signal_features")
  x$signal       <- .as_checked_factor(x$signal, cfg$signals, "signal_features$signal")
  x$signal_class <- .as_checked_factor(x$signal_class, cfg$signal_classes, "signal_features$signal_class")
  x
}

load_intervention_features <- function(path, cfg, site) {
  x <- as.data.frame(arrow::read_parquet(path))
  x <- .stamp_site(x, site, "intervention_features")
  x$intervention <- .as_checked_factor(x$intervention, cfg$interventions_extracted,
                                       "intervention_features$intervention")
  x$shape <- .as_checked_factor(x$shape, cfg$intervention_shapes,
                                "intervention_features$shape")
  x
}

#' Audit-only since 2026-08-25 (see load_tables). `o_flag` is read as-is and
#' never coerced to a factor: it is not a modelled variable any more, so
#' declaring its levels here would assert a design contract that no longer
#' exists. Nothing but validator check 4 looks at this table.
load_ordering <- function(path, cfg, site) {
  x <- as.data.frame(arrow::read_parquet(path))
  x <- .stamp_site(x, site, "ordering")
  x$signal       <- .as_checked_factor(x$signal, cfg$signals, "ordering$signal")
  x$intervention <- .as_checked_factor(x$intervention, cfg$interventions_extracted,
                                       "ordering$intervention")
  x$excursion_side <- .as_checked_factor(x$excursion_side, c("low", "high"),
                                         "ordering$excursion_side")
  x
}

#' The severity-comparator table (APACHE and SOFA). OPTIONAL, like `ordering`.
#'
#' Nothing here is a modelling input. It is the comparison arms' raw material,
#' scored by R/09c_apache.R and R/09d_sofa.R, and no model frame reads it. It is loaded through
#' the same adapter as everything else so that it gets the same column contract,
#' the same declared factor levels and the same schema-equality check across
#' sites — which is the whole reason the two SQL files were written to one
#' output contract in the first place.
#'
#' TWO REPAIRS, and both are assertions in the normal case.
#'
#' The `-1` sentinel. eICU codes "not measured" as -1 throughout `apacheApsVar`
#' and `apachePatientResult`, and every guard is already applied in SQL. It is
#' re-applied here for the same reason the weight guard is: a -1 potassium
#' scores four points on APACHE II's low band, so a sentinel that slipped
#' through would inflate the baseline for every unmeasured stay and nothing
#' downstream would show it. MIMIC has no sentinel convention, so this is a
#' no-op there — which is the point. It does not branch on site; it applies one
#' rule to both, and one of them never trips it.
#'
#' `admission_class` and `aps_native_version` become factors with declared
#' levels, so an unexpected value is an error rather than a silent NA.
load_severity <- function(path, cfg, site) {
  x <- as.data.frame(arrow::read_parquet(path))
  x <- .stamp_site(x, site, "severity")

  sentinel_cols <- intersect(
    c(grep("^(ap2_|sofa_)", names(x), value = TRUE),
      "aps_native", "aps_native_prob", "severity_total_native"),
    names(x))
  sentinel_cols <- sentinel_cols[vapply(x[sentinel_cols], is.numeric, logical(1))]
  n_sent <- 0L
  for (v in sentinel_cols) {
    bad <- !is.na(x[[v]]) & x[[v]] < 0
    n_sent <- n_sent + sum(bad)
    x[[v]][bad] <- NA
  }
  if (n_sent > 0L) {
    message(sprintf("  severity: %d negative sentinel value(s) set to NA across %d column(s)",
                    n_sent, length(sentinel_cols)))
  }

  x$admission_class <- .as_checked_factor(
    x$admission_class, c("nonoperative", "emergency_postop", "elective_postop"),
    "severity$admission_class")
  x$aps_native_version <- .as_checked_factor(
    x$aps_native_version, c("apsiii", "apache_iv", "apache_iva", "apache_none"),
    "severity$aps_native_version")
  x
}

#' Factor with declared levels, erroring on any value outside them.
#'
#' `factor()` would turn an unknown value into NA and carry on. That is the
#' failure mode spec §8 check 2 exists to prevent, so it is an error here.
.as_checked_factor <- function(x, levels, what) {
  x <- as.character(x)
  bad <- setdiff(unique(x[!is.na(x)]), as.character(levels))
  if (length(bad)) abort_values(paste0(what, ": value(s) outside the declared level set"), bad)
  factor(x, levels = as.character(levels))
}
