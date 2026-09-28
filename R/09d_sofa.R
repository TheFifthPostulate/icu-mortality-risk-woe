# R/09d_sofa.R ---------------------------------------------------------------
# The SOFA comparison arm. Sibling of R/09c_apache.R and it reuses that file's
# `.ap2_band()`, because a band table is a band table.
#
# WHY SOFA IS A DIFFERENT ARM FROM APACHE II, and not a second one. APACHE II is
# pure physiology, so the honest comparison is against `L_meas`. SOFA folds
# vasopressor dose into its cardiovascular component, so the honest comparison
# is against `L_full`, which carries `L_intv`. The two scores externalise the
# project's own meas/full split instead of duplicating each other. That is the
# specific reason `docs/v2_analysis_tiering.md` item 5 -- one severity score,
# not both -- does not settle the question: it is a rule about discrimination
# baselines and it is right about those.
#
# SOFA IS ALSO THE THING THE DESIGN GENERALISES.
# `docs/v2_analytical_design_plan.md` line 12 names it as lineage: SOFA already
# bundles measurement with intervention, ad hoc, for ONE organ system. The
# paper's representational claim is that this is generalised to twelve
# signal-intervention pairs. `L_mbp^full` against SOFA's cardiovascular
# component is therefore the sharpest single comparison in the project: same
# signal, same intervention, both bundled, one hand-designed and one fitted.
#
# AND IT IS THE ONLY PUBLISHED SCORE THAT DECOMPOSES ONTO `config/domains.csv`,
# which already carries a `sofa_organ` column, frozen before any result existed.
# `sofa_domain_table()` at the foot of this file is what that buys.
#
# NO PATHS, NO CLOCK (hard rules 7 and 9). AGGREGATES ONLY (hard rule 1).
#
# SOURCE. Vincent JL, Moreno R, Takala J, et al. The SOFA (Sepsis-related Organ
# Failure Assessment) score to describe organ dysfunction/failure. Intensive
# Care Med 1996;22(7):707-710. Bands transcribed from that paper's Table 1 and
# neither tuned nor re-fitted.
# ----------------------------------------------------------------------------

# --- the bands --------------------------------------------------------------

#' SOFA's five band-scored organs. Cardiovascular is not here: it is a
#' branching rule over dose and agent rather than a band table, and it is
#' handled explicitly in `sofa_score()`.
#'
#' Every band is ONE-SIDED, unlike APACHE II's, because SOFA measures organ
#' dysfunction in a single direction. That is why this file takes one value per
#' organ where R/09c takes a min/max pair: there is no worse end to choose.
SOFA_BANDS <- list(
  # PaO2 / FiO2, mmHg. Scores 3 and 4 additionally require respiratory support;
  # `sofa_score()` applies that gate after the band.
  pf         = list(breaks = c(100, 200, 300, 400),
                    points = c(4, 3, 2, 1, 0)),
  # platelets, x10^3/uL
  platelet   = list(breaks = c(20, 50, 100, 150),
                    points = c(4, 3, 2, 1, 0)),
  # bilirubin, mg/dL
  bilirubin  = list(breaks = c(1.2, 2.0, 6.0, 12.0),
                    points = c(0, 1, 2, 3, 4)),
  # Glasgow Coma Scale
  gcs        = list(breaks = c(6, 10, 13, 15),
                    points = c(4, 3, 2, 1, 0)),
  # creatinine, mg/dL
  creatinine = list(breaks = c(1.2, 2.0, 3.5, 5.0),
                    points = c(0, 1, 2, 3, 4)),
  # urine output, mL/day. The renal score is the worse of this and creatinine.
  urine      = list(breaks = c(200, 500),
                    points = c(4, 3, 0))
)

#' Structural assertions, same role as `check_apache_tables()`.
#'
#' Cannot check the numbers against Vincent 1996 -- nothing in R can. Checks
#' what a transcription slip breaks: lengths, ordering, the 0-4 range that
#' every SOFA organ is bounded by, and the presence of a zero band. A table
#' with no zero band would score every patient as dysfunctional and the arm
#' would look absurdly strong for a reason that has nothing to do with SOFA.
check_sofa_tables <- function() {
  for (nm in names(SOFA_BANDS)) {
    b <- SOFA_BANDS[[nm]]
    if (length(b$points) != length(b$breaks) + 1L) {
      stop("SOFA_BANDS[", nm, "]: points must be one longer than breaks", call. = FALSE)
    }
    if (is.unsorted(b$breaks)) {
      stop("SOFA_BANDS[", nm, "]: breaks must be ascending", call. = FALSE)
    }
    if (any(b$points < 0) || any(b$points > 4)) {
      stop("SOFA_BANDS[", nm, "]: points outside [0, 4]", call. = FALSE)
    }
    if (!any(b$points == 0)) {
      stop("SOFA_BANDS[", nm, "]: no zero band; every value would score", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' SOFA's six organs, and the five the RECOMPUTED score actually sums.
#'
#' MEASURED 2026-08-31, and this is why respiration is not in `SOFA_ORGANS`.
#' `mimiciv_derived.bg.fio2` is NULL on 93% of stays, so a P/F ratio is
#' computable on only 12.3% of the cohort against 48.1% coverage for PaO2
#' alone. Of the stays where we scored 0 and the derived concept scored
#' something, 96% were missing the input entirely rather than computing a ratio
#' above 400. eICU is thinner still: arterial gas coverage there is 26.4%.
#'
#' THE SCORING LOGIC WAS NEVER WRONG, and that distinction is the whole reason
#' this is an exclusion rather than a bug fix. Restricted to stays where the
#' input exists, respiration agrees with the derived concept at 0.507 exact and
#' Spearman 0.665, and our median P/F falls monotonically from 341 at native
#' tier 0 to 130 at tier 4. It is a COVERAGE failure, not a construction one.
#'
#' Dropping it takes the recomputed total from Spearman 0.875 and 58.0% exact
#' agreement to **0.956 and 79.7%** — from a score that fails its own validation
#' gate to one that passes it convincingly. A five-organ SOFA that is right is
#' worth more than a six-organ SOFA whose sixth organ is noise.
#'
#' What is lost is stated rather than hidden: the respiration row of
#' `sofa_domain_table()`, and any claim about SOFA's respiratory organ. The
#' NATIVE six-organ SOFA is unaffected and remains the MIMIC within-site
#' comparator; it is only the score that TRAVELS that drops to five.
#'
#' Restoring respiration needs FiO2 built from `chartevents` (itemid 223835)
#' plus ventilator settings, carried forward and joined to each arterial gas —
#' which is what mimic-code's own SOFA does — with an eICU analogue. Deferred,
#' not abandoned. `config/sofa.organs` is the one place to change it.
SOFA_ALL_ORGANS <- c("respiration", "coagulation", "liver",
                     "cardiovascular", "cns", "renal")
SOFA_ORGANS <- c("coagulation", "liver", "cardiovascular", "cns", "renal")

# --- the score --------------------------------------------------------------

#' SOFA's cardiovascular component: the one organ that is a branching rule.
#'
#' THIS IS THE COMPONENT THE WHOLE ARM IS FOR, and it is also the one carrying
#' the arm's only real approximation. Vincent 1996 scores by AGENT:
#'
#'   MAP >= 70                                                       -> 0
#'   MAP  < 70                                                       -> 1
#'   dopamine <= 5, or dobutamine at any dose                        -> 2
#'   dopamine  > 5, or adrenaline <= 0.1, or noradrenaline <= 0.1    -> 3
#'   dopamine  > 15, or adrenaline > 0.1, or noradrenaline > 0.1     -> 4
#'
#' `dx_nee_peak` is a norepinephrine-EQUIVALENT dose, so it collapses the agent
#' distinction. The consequence is specific and is stated rather than absorbed:
#' a patient on dopamine at 5 mcg/kg/min converts to roughly 0.05 NEE and lands
#' at 3 here where the original gives 2. `inotrope` is carried separately so the
#' dobutamine tier is recovered exactly, which is why tier 2 below is an
#' inotrope test rather than a dose test.
#'
#' The collapse is arguably an improvement on a 1996 agent list -- NEE is how
#' vasopressor burden is compared in modern practice -- but it is a deviation
#' from the published definition, and audit A11 in the SQL sizes the resulting
#' disagreement against the native score. Report that disagreement rate; do not
#' report the point estimate alone.
#' A DOSE THAT DID NOT PARSE IS NOT AN ABSENT PRESSOR. MEASURED 2026-08-30:
#' 5,218 of eICU's 13,430 vasopressor stays -- 38.9% -- carry `ever_active = 1`
#' with a NULL `dx_nee_peak`, because eICU's dose lives in free-text
#' `infusiondrug.drugrate`. An earlier version of this function treated a NULL
#' dose as zero dose, which scored those stays on the MAP tier alone and
#' silently weakened eICU's cardiovascular component -- the exact component
#' this arm exists to compare, weakened at the exact site where the comparison
#' is made. `vaso_active` now floors such a stay at tier 3, the lowest
#' catecholamine tier: they are demonstrably on a pressor, so 0 and 1 are
#' excluded, and 4 cannot be justified without a dose.
.sofa_cardio <- function(mbp_min, nee_peak, inotrope, vaso_active) {
  s <- rep(NA_real_, length(mbp_min))
  have <- !is.na(mbp_min)
  s[have] <- ifelse(mbp_min[have] < 70, 1, 0)

  nee <- ifelse(is.na(nee_peak), 0, nee_peak)
  ino <- ifelse(is.na(inotrope), 0L, inotrope)
  vas <- ifelse(is.na(vaso_active), 0L, vaso_active)
  # On a pressor, dose unknown or out of the plausibility guard.
  dose_unknown <- vas == 1L & is.na(nee_peak)

  # Tiers are applied upward, so the highest reached wins. `pmax` against the
  # MAP tier rather than assignment, because a patient can be on a pressor and
  # still have MAP >= 70 -- which is the entire point of the component.
  s <- pmax(s, ifelse(ino == 1L, 2, 0), na.rm = FALSE)
  s <- pmax(s, ifelse((nee > 0 & nee <= 0.1) | dose_unknown, 3, 0), na.rm = FALSE)
  s <- pmax(s, ifelse(nee > 0.1, 4, 0), na.rm = FALSE)

  # If MAP was never measured but support was given, the support tiers still
  # apply: the organ is observed to be supported. Only a stay with neither is NA.
  no_map <- is.na(mbp_min)
  s[no_map] <- ifelse(nee[no_map] > 0.1, 4,
               ifelse(nee[no_map] > 0 | dose_unknown[no_map], 3,
               ifelse(ino[no_map] == 1L, 2, NA_real_)))
  s
}

#' SOFA, from the extracted inputs.
#'
#' MISSING ORGANS SCORE ZERO, with the same caveat and the same direction as the
#' APACHE arm: a stay with no arterial gas scores 0 on respiration exactly as a
#' stay with a normal gas does, so thin coverage makes SOFA look weaker than it
#' is and therefore flatters our method. `n_organs_present` travels with every
#' row and `tests/metrics_severity.R` restricts on it by default.
#'
#' @param sv           a data frame with the `sofa_*` and shared `ap2_*` columns
#'   of CONTRACT$severity
#' @param resp_support "any" counts non-invasive support towards SOFA's
#'   respiration tiers 3 and 4, "invasive" counts only invasive ventilation.
#'   Declared in `config/sofa.resp_support`. Vincent 1996 says mechanical
#'   ventilation; mimic-code's reference implementation counts non-invasive
#'   support too, and SQL audit A10 measures which the derived concept actually
#'   used. Set this from that measurement, before running the arm.
#' @return data frame: one `sofa_*` column per organ, plus `sofa` (the total)
#'   and `n_organs_present`.
#' @param organs which organs enter the TOTAL. Every organ is still computed
#'   and returned; this governs `sofa` and `n_organs_present` only, so
#'   respiration stays visible as a diagnostic while being excluded from the
#'   score. Declared in `config/sofa.organs`; see `SOFA_ORGANS` for why
#'   respiration is out.
sofa_score <- function(sv, resp_support = c("any", "invasive"),
                       organs = SOFA_ORGANS) {
  check_sofa_tables()
  resp_support <- match.arg(resp_support)
  bad <- setdiff(organs, SOFA_ALL_ORGANS)
  if (length(bad)) abort_values("sofa_score: unknown organ(s)", bad)
  if (!length(organs)) stop("sofa_score: `organs` is empty", call. = FALSE)

  need <- c("ap2_pf_min", "ap2_mbp_min", "ap2_creatinine_max",
            "ap2_gcs_min_native", "ap2_urine_ml_24h",
            "sofa_platelet_min", "sofa_bilirubin_max",
            "sofa_vent_invasive", "sofa_vent_noninvasive",
            "sofa_vasopressor", "sofa_nee_peak", "sofa_inotrope")
  miss <- setdiff(need, names(sv))
  if (length(miss)) abort_values("sofa_score: missing input columns", miss)

  b <- function(v, x) .ap2_band(x, SOFA_BANDS[[v]]$breaks, SOFA_BANDS[[v]]$points)

  # Respiration. `ap2_pf_min` is the worst P/F RATIO, extracted as its own
  # column because it comes from a DIFFERENT arterial gas than the worst
  # PaO2 that feeds APACHE II's oxygenation term.
  #
  # FIXED 2026-08-31. The first version computed the ratio here as
  # `ap2_pao2 / (ap2_fio2/100)` from the worst-PaO2 gas, and it failed the
  # validation gate badly: 33.8% exact agreement against the derived SOFA
  # concept and a Spearman of 0.086, with our mean respiration at 0.285
  # against a native 1.916. The near-zero RANK correlation is what
  # identified the cause -- a units error would have preserved rank, so
  # different gases were being selected for different patients. They were:
  # the lowest-PaO2 gas is often a room-air sample scoring P/F 452 and tier
  # 0, while the gas that actually determines the score sits at high FiO2.
  # The selection now happens in SQL, once, on the right criterion.
  resp <- b("pf", sv$ap2_pf_min)
  supported <- if (resp_support == "invasive") {
    sv$sofa_vent_invasive == 1L
  } else {
    sv$sofa_vent_invasive == 1L | sv$sofa_vent_noninvasive == 1L
  }
  supported[is.na(supported)] <- FALSE
  # Tiers 3 and 4 require respiratory support. Without it the component is
  # capped at 2 -- the published rule, and the reason the ventilation flag is an
  # input to a physiology score at all.
  resp[!supported & !is.na(resp)] <- pmin(resp[!supported & !is.na(resp)], 2)

  # Renal is the WORSE of the creatinine and urine tiers, not the creatinine
  # tier with urine as a tie-break. A stay can be anuric at a normal creatinine
  # early in an injury, which is exactly the case the urine limb exists for.
  renal_cr <- b("creatinine", sv$ap2_creatinine_max)
  renal_uo <- b("urine",      sv$ap2_urine_ml_24h)
  renal    <- pmax(renal_cr, renal_uo, na.rm = TRUE)

  p <- data.frame(
    sofa_respiration    = resp,
    sofa_coagulation    = b("platelet",  sv$sofa_platelet_min),
    sofa_liver          = b("bilirubin", sv$sofa_bilirubin_max),
    sofa_cardiovascular = .sofa_cardio(sv$ap2_mbp_min, sv$sofa_nee_peak,
                                       sv$sofa_inotrope, sv$sofa_vasopressor),
    sofa_cns            = b("gcs",       sv$ap2_gcs_min_native),
    sofa_renal          = renal,
    stringsAsFactors = FALSE)

  P <- as.matrix(p)
  present <- !is.na(P)
  P[is.na(P)] <- 0
  p[colnames(P)] <- P
  # The total and the coverage count run over the DECLARED organs only. Every
  # organ is still returned, so an excluded one remains inspectable rather than
  # disappearing -- which is what keeps the exclusion honest.
  scored <- paste0("sofa_", organs)
  p$n_organs_present <- rowSums(present[, scored, drop = FALSE])
  p$sofa <- rowSums(P[, scored, drop = FALSE])
  # WHICH organs were actually scored. Without this every downstream comparison
  # silently conflates "scored 0 because normal" with "scored 0 because the
  # input was missing" -- exactly how a 93%-missing FiO2 hid behind a
  # plausible-looking respiration score on 2026-08-31.
  #
  # CARRIED AS COLUMNS, NOT AN ATTRIBUTE. The first version used
  # `attr(p, "organ_present")` to leave the frame's shape untouched. That was
  # wrong: attributes are for scalar metadata, and `[.data.frame` carries this
  # one through row subsetting WITHOUT subsetting it, so `sf[keep, ]` returned a
  # frame of one length beside a presence matrix of another. The next comparison
  # recycled them against each other, `sd()` went NA, and an `if (NA)` blew up
  # three functions away. Row-aligned data belongs in rows.
  #
  # `organs` stays an attribute because it is genuinely scalar metadata -- the
  # declared organ set, not a per-stay fact -- and every reader defaults
  # sensibly when it is dropped.
  for (o in SOFA_ALL_ORGANS) p[[paste0("has_", o)]] <- unname(present[, paste0("sofa_", o)])
  attr(p, "organs") <- organs
  p
}

#' Restrict to stays whose SOFA is actually a SOFA. Mirrors
#' `apache_complete_mask()`; `min_organs` is declared in config before a run.
#'
#' NOTE the denominator moved with `SOFA_ORGANS`: `sofa_n_organs_present` in the
#' extraction still counts SIX, but `sofa_score()` recomputes it over the
#' declared set. Pass the scored frame, not the raw table, or the floor means
#' something different from what it says.
#' TWO COLUMNS CAN CARRY THIS COUNT AND THEY MEAN DIFFERENT THINGS.
#' `sofa_score()` emits `n_organs_present` over the DECLARED organ set;
#' the extraction emits `sofa_n_organs_present` over all SIX. The scored frame
#' is preferred, because a 4-of-6 floor applied to a five-organ score is not the
#' floor the config asked for.
#'
#' THE `stop()` IS THE POINT. Reading a column that does not exist returns NULL,
#' `is.na(NULL)` is `logical(0)`, and a zero-length mask then propagates in
#' total silence: `sum(!keep)` reports 0 dropped, `y[keep]` is empty, and
#' `all(is.na(numeric(0)))` is TRUE — so every downstream cell reports itself as
#' entirely NA and the run dies much later on an unrelated `rbind`. That is
#' exactly what happened on 2026-08-31 when this function was handed the scored
#' frame while still looking for the extraction's column name. Fail here, where
#' the name is wrong, not fifty lines downstream.
sofa_complete_mask <- function(sv, min_organs = 4L) {
  v <- if (!is.null(sv[["n_organs_present"]])) {
    sv[["n_organs_present"]]
  } else if (!is.null(sv[["sofa_n_organs_present"]])) {
    message("sofa_complete_mask: using the extraction's six-organ count; pass ",
            "the sofa_score() frame to count over the declared organ set")
    sv[["sofa_n_organs_present"]]
  } else {
    stop("sofa_complete_mask: neither `n_organs_present` (from sofa_score()) ",
         "nor `sofa_n_organs_present` (from the extraction) is present. ",
         "A missing column here returns a ZERO-LENGTH mask that silently ",
         "empties every downstream cell.", call. = FALSE)
  }
  !is.na(v) & v >= min_organs
}

# --- validation against the native concept ----------------------------------

#' Per-organ agreement between our recomputation and the derived concept.
#'
#' THE GATE ON THE WHOLE ARM, and a stricter one than the APACHE equivalent.
#' `apache_native_agreement()` compares two DIFFERENT scores and can only ask
#' for strong correlation. This compares our SOFA against a published SOFA, so
#' it can ask for near-identity: exact agreement should be high and the mean
#' absolute difference near zero on every organ except cardiovascular, where the
#' NEE collapse documented in `.sofa_cardio()` guarantees some disagreement.
#'
#' An organ other than cardiovascular disagreeing materially means a band, a
#' unit, or a missing-value rule is wrong, and the eICU SOFA -- which has no
#' native score to check against -- would inherit that error invisibly. Fix it
#' before the arm is reported at either site.
#'
#' Returns NULL when no native score is present, which is the normal state at
#' eICU and must not be treated as a pass.
#' @param scored the organs that enter the total. Taken from the `organs`
#'   attribute when present, but accepted explicitly because an attribute can be
#'   dropped by subsetting and a silently wrong organ set would compare our
#'   five-organ sum against a six-organ native total.
sofa_native_agreement <- function(ours, sv,
                                  scored = attr(ours, "organs") %||% SOFA_ORGANS) {
  # EVERY organ is reported, including any excluded from the total, because an
  # exclusion that stops being visible stops being defensible.
  organs <- SOFA_ALL_ORGANS
  if (!("sofa_native" %in% names(sv)) || all(is.na(sv$sofa_native))) return(NULL)

  # `ours` is zero-filled, so a raw comparison cannot tell a genuine 0 from an
  # unmeasured organ. The `has_` columns separate them, and `frac_input` /
  # `exact_computable` are what make a coverage failure visible as a coverage
  # failure rather than as a scoring disagreement.
  rows <- lapply(organs, function(o) {
    a <- ours[[paste0("sofa_", o)]]
    b <- sv[[paste0("sofa_native_", o)]]
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) return(NULL)
    hv <- ours[[paste0("has_", o)]]
    have <- if (is.null(hv)) rep(TRUE, length(a)) else hv
    stopifnot(length(have) == length(a))
    ok2 <- ok & have
    data.frame(organ = o, n = sum(ok),
               frac_input = round(mean(have[ok]), 4),
               exact_agree = round(mean(a[ok] == b[ok]), 4),
               exact_computable = if (any(ok2)) round(mean(a[ok2] == b[ok2]), 4) else NA_real_,
               n_computable = sum(ok2),
               within_1    = round(mean(abs(a[ok] - b[ok]) <= 1), 4),
               mean_abs_diff = round(mean(abs(a[ok] - b[ok])), 4),
               mean_ours   = round(mean(a[ok]), 3),
               mean_native = round(mean(b[ok]), 3),
               # Same isTRUE() guard as the computable column below. A
               # constant organ is a real state at a validation site, and it
               # should return NA rather than a warning and an NA.
               spearman = if (isTRUE(stats::sd(a[ok]) > 0) &&
                              isTRUE(stats::sd(b[ok]) > 0))
                 round(stats::cor(a[ok], b[ok], method = "spearman"), 4) else NA_real_,
               # isTRUE(), because sd() of a constant-or-empty vector is NA or 0
               # and `if (NA)` is an error rather than a FALSE.
               spearman_computable = if (isTRUE(sum(ok2) > 2) &&
                                         isTRUE(stats::sd(a[ok2]) > 0) &&
                                         isTRUE(stats::sd(b[ok2]) > 0))
                 round(stats::cor(a[ok2], b[ok2], method = "spearman"), 4) else NA_real_,
               stringsAsFactors = FALSE)
  })
  # The TOTAL compares like with like: our declared-organ sum against the NATIVE
  # sum over the SAME organs, not against the native six-organ total.
  nat_cols <- paste0("sofa_native_", scored)
  nat_tot  <- rowSums(sv[, nat_cols, drop = FALSE])
  ok <- !is.na(ours$sofa) & !is.na(nat_tot)
  tot <- data.frame(organ = paste0("TOTAL (", length(scored), " organs)"),
                    n = sum(ok),
                    frac_input = round(mean(ours$n_organs_present[ok] / length(scored)), 4),
                    exact_agree = round(mean(ours$sofa[ok] == nat_tot[ok]), 4),
                    exact_computable = NA_real_,
                    # NA rather than a count: "computable" is a per-organ idea
                    # and a total has no single input to be present or absent.
                    # It must still be here, or rbind fails on a column
                    # mismatch that says nothing about what is wrong.
                    n_computable = NA_integer_,
                    within_1 = round(mean(abs(ours$sofa[ok] - nat_tot[ok]) <= 1), 4),
                    mean_abs_diff = round(mean(abs(ours$sofa[ok] - nat_tot[ok])), 4),
                    mean_ours = round(mean(ours$sofa[ok]), 3),
                    mean_native = round(mean(nat_tot[ok]), 3),
                    spearman = round(stats::cor(ours$sofa[ok], nat_tot[ok],
                                                method = "spearman"), 4),
                    spearman_computable = NA_real_,
                    stringsAsFactors = FALSE)
  out <- do.call(rbind, c(Filter(Negate(is.null), rows), list(tot)))
  rownames(out) <- NULL
  out
}

# --- the domain-level comparison --------------------------------------------

#' Read the frozen domain partition and its SOFA organ mapping.
#'
#' `config/domains.csv` is an aggregation guide for layer 2 and changes no
#' model. Six of its eleven domains carry a `sofa_organ`; the other five are the
#' leftovers that SOFA has no organ for, and they are reported as such rather
#' than stretched onto one.
load_domains <- function(path = "config/domains.csv") {
  d <- utils::read.csv(path, colClasses = "character", na.strings = "")
  stopifnot(all(c("signal", "domain", "sofa_organ") %in% names(d)))
  # domains.csv says "respiratory"; SOFA says "respiration". One rename, here,
  # so neither file has to be edited to match the other's vocabulary.
  d$sofa_organ[d$sofa_organ %in% "respiratory"] <- "respiration"
  d
}

#' Domain-level comparison: each `D_k` against its SOFA organ.
#'
#' THIS IS WHAT THE SOFA ARM BUYS THAT NOTHING ELSE DOES. Layer 2 emits one
#' number per domain and, until now, that number could only be compared with
#' itself. Six domains carry a SOFA organ, so six of them get an externally
#' defined referent on organ definitions frozen before any result existed.
#'
#' `D_k` is the sum of the L's in domain k. Layer 2's weights do not exist yet
#' (Sigma and the weight vector are still unbuilt, docs/v2_state_20260828.md
#' section 3.4), so this is the equal-weight case, which is what `D_k` reduces
#' to at w = 1. When the weights land, pass them in and the comparison sharpens
#' rather than changing shape.
#'
#' Three numbers per domain, and the third is the interesting one:
#'   auroc_D       what the domain's summed L discriminates on its own
#'   auroc_organ   what SOFA's organ score discriminates on its own
#'   spearman      whether the two ORDER patients the same way, which is a
#'                 question about agreement rather than about accuracy. Two
#'                 constructions can discriminate equally well while ranking
#'                 different patients as sick, and that difference is the
#'                 finding.
#'
#' EACH AUROC HAS AN AUPRC BESIDE IT AS OF 2026-09-07, and this table is the one
#' where the difference is most likely to bite. A SOFA ORGAN SCORE IS A SMALL
#' INTEGER -- 0 to 4 -- so it partitions the cohort into five blocks and orders
#' patients within a block not at all. AUROC handles ties gracefully and reports
#' the average over them; AUPRC does not, because precision at the top of the
#' ranking is exactly what a five-level score cannot resolve. A domain whose
#' `auroc_D` barely beats `auroc_organ` while its `auprc_D` beats it clearly is
#' a domain where the continuous evidence is doing real work where the deaths
#' are, and the AUROC-only table said nothing about that either way.
#'
#' `auprc_lift_*` divides by the event rate on the SAME rows the AUROC is
#' computed on (`ok`), so the two columns describe one population.
#'
#' TWO VERSIONS OF `D_k`, AND THE DIFFERENCE BETWEEN THEM IS ITSELF A NUMBER.
#' `domains.csv` tags the SIGNAL, not the domain, and its "leftover rule" puts
#' untagged signals into tagged domains: `heart_rate` joins `mbp` in
#' hemodynamic_support, `resp_rate` joins `spo2` in respiratory_support, `bun`
#' joins the renal pair. So there are two defensible aggregates and reporting
#' only one would hide a choice:
#'
#'   auroc_D         `D_k` over the WHOLE domain. This is what layer 2 actually
#'                   emits and what the paper reports, so it is primary.
#'   auroc_D_tagged  the sum over only the signals carrying that SOFA organ.
#'                   The like-for-like aggregate, over exactly the measurements
#'                   SOFA itself uses.
#'
#' Their difference is what the leftover rule contributes, which is a property
#' of the frozen partition worth knowing on its own.
#'
#' @param M       an L matrix (stays x signals), rows aligned to `sofa_tab`
#' @param sofa_tab output of `sofa_score()`, same rows
#' @param w       optional named weight vector over signals; defaults to 1
sofa_domain_table <- function(M, sofa_tab, y, domains, w = NULL,
                              organs = attr(sofa_tab, "organs") %||% SOFA_ORGANS) {
  agg <- function(sigs) {
    sigs <- intersect(sigs, colnames(M))
    if (!length(sigs)) return(NULL)
    ww <- if (is.null(w)) rep(1, length(sigs)) else unname(w[sigs])
    list(sigs = sigs, D = as.numeric(M[, sigs, drop = FALSE] %*% ww))
  }
  # A zero-variance column makes cor() warn and return NA. That is a real state
  # at a validation site where an organ can be degenerate, so it is handled
  # rather than warned about.
  sp <- function(a, b) {
    if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    round(stats::cor(a, b, method = "spearman"), 4)
  }

  # Only domains whose SOFA organ is actually SCORED. A domain whose organ was
  # excluded from the total has no comparator, and reporting it against an
  # organ the recomputed score does not use would be quietly comparing against
  # a number we have declared unreliable.
  domains <- domains[is.na(domains$sofa_organ) | domains$sofa_organ %in% organs, ,
                     drop = FALSE]
  covered <- unique(domains$domain[!is.na(domains$sofa_organ)])
  rows <- lapply(covered, function(dm) {
    z <- domains[domains$domain == dm, , drop = FALSE]
    org <- unique(stats::na.omit(z$sofa_organ))
    # One domain must map to at most one organ, or `D_k` has no single referent
    # and the comparison is undefined. Assert rather than silently take the
    # first: domains.csv is frozen, so this can only fire if it is edited.
    if (length(org) > 1L) {
      stop("domain '", dm, "' maps to more than one SOFA organ: ",
           paste(org, collapse = ", "), call. = FALSE)
    }
    a_all <- agg(z$signal)
    a_tag <- agg(z$signal[!is.na(z$sofa_organ)])
    if (is.null(a_all) || is.null(a_tag)) return(NULL)
    S  <- sofa_tab[[paste0("sofa_", org)]]
    ok <- !is.na(a_all$D) & !is.na(S)
    data.frame(
      sofa_organ = org,
      domain     = dm,
      n_signals  = length(a_all$sigs),
      n_tagged   = length(a_tag$sigs),
      signals    = paste(a_all$sigs, collapse = ","),
      auroc_D        = round(.auroc(a_all$D[ok], y[ok]), 5),
      auroc_D_tagged = round(.auroc(a_tag$D[ok], y[ok]), 5),
      auroc_organ    = round(.auroc(S[ok], y[ok]), 5),
      delta      = round(.auroc(a_all$D[ok], y[ok]) - .auroc(S[ok], y[ok]), 5),
      auprc_D        = round(.auprc(a_all$D[ok], y[ok]), 5),
      auprc_D_tagged = round(.auprc(a_tag$D[ok], y[ok]), 5),
      auprc_organ    = round(.auprc(S[ok], y[ok]), 5),
      delta_auprc = round(.auprc(a_all$D[ok], y[ok]) - .auprc(S[ok], y[ok]), 5),
      event_rate  = round(mean(y[ok]), 5),
      auprc_lift_D     = round(.auprc(a_all$D[ok], y[ok]) / mean(y[ok]), 4),
      auprc_lift_organ = round(.auprc(S[ok], y[ok]) / mean(y[ok]), 4),
      spearman   = sp(a_all$D[ok], S[ok]),
      mean_organ = round(mean(S[ok]), 3),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out[order(-out$auroc_organ), , drop = FALSE]
}

#' The domains SOFA has no organ for at all. Reported, not hidden.
#'
#' A domain counts as uncovered only when NONE of its signals carries a SOFA
#' organ. Splitting on the signal-level tag instead would put
#' hemodynamic_support in both tables, because `mbp` is tagged and `heart_rate`
#' is not, and the two tables would then double-count.
#'
#' `domains.csv` is disjoint and complete over the 19 signals, so these are
#' exactly the evidence SOFA cannot represent in any form. That belongs beside
#' the organ comparison, because "our domains beat SOFA's organs" is a much
#' weaker claim than that plus "and there are more domains SOFA has no way to
#' express".
sofa_uncovered_domains <- function(M, y, domains, w = NULL) {
  covered <- unique(domains$domain[!is.na(domains$sofa_organ)])
  dd <- domains[!domains$domain %in% covered, , drop = FALSE]
  rows <- lapply(split(dd, dd$domain), function(z) {
    sigs <- intersect(z$signal, colnames(M))
    if (!length(sigs)) return(NULL)
    ww <- if (is.null(w)) rep(1, length(sigs)) else unname(w[sigs])
    D  <- as.numeric(M[, sigs, drop = FALSE] %*% ww)
    data.frame(domain = z$domain[1], n_signals = length(sigs),
               signals = paste(sigs, collapse = ","),
               auroc_D = round(.auroc(D, y), 5),
               # These domains have NO SOFA comparator at all, so the only
               # available floor is the event rate itself -- which is what
               # `auprc_lift_D` is against. A lift near 1 on an uncovered
               # domain would mean the evidence SOFA cannot express is also
               # evidence that does not separate the deaths, and the AUROC
               # column alone could not say so.
               auprc_D = round(.auprc(D, y), 5),
               auprc_lift_D = round(.auprc(D, y) / mean(y), 4),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out[order(-out$auroc_D), , drop = FALSE]
}
