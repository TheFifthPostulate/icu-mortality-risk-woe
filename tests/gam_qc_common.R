# tests/gam_qc_common.R --------------------------------------------------------
# Helpers shared by tests/gam_qc.R and tests/gam_qc_bootstrap.R. ONE definition
# each, sourced by both, for the reason tests/audit_common.R exists: the first
# draft carried `.neutral_row()` in both scripts and an extremum counter in
# each, and two counters that can drift apart make `extrema_agree` a comparison
# between two conventions rather than between two fits.
#
# Nothing here reads data, builds a path or reads the clock. Every function
# takes what it needs as an argument; the config values come from the caller.
# ------------------------------------------------------------------------------

#' A one-row frame with every model variable at a neutral value.
#'
#' `type = "terms"` evaluates a smooth from its own covariate only, so the
#' other columns never reach a reported number; they exist because predict()
#' needs the full frame and na.fail refuses a gap. Factor columns take their
#' first level so the prediction frame carries a level the model has seen.
qc_neutral_row <- function(d, vars) {
  out <- lapply(vars, function(v) {
    x <- d[[v]]
    if (is.factor(x)) factor(levels(x)[1], levels = levels(x))
    else if (is.logical(x)) FALSE
    else stats::median(x, na.rm = TRUE)
  })
  names(out) <- vars
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' The evaluation grid for one covariate, with the bin edges around each point.
#'
#' A covariate with at most `max_distinct` distinct values is evaluated AT
#' those values (`kind = "atoms"`); otherwise on `n_points` uniform points over
#' the [q_lo, q_hi] quantiles of its training distribution (`kind = "uniform"`).
qc_grid_of <- function(x, n_points, q_lo, q_hi, max_distinct) {
  ux <- sort(unique(x))
  if (length(ux) <= as.integer(max_distinct)) {
    edges <- if (length(ux) > 1L) c(-Inf, (ux[-1] + ux[-length(ux)]) / 2, Inf) else c(-Inf, Inf)
    return(list(grid = ux, edges = edges, kind = "atoms"))
  }
  q <- stats::quantile(x, c(as.numeric(q_lo), as.numeric(q_hi)), names = FALSE, type = 7)
  grid <- seq(q[1], q[2], length.out = as.integer(n_points))
  step <- grid[2] - grid[1]
  list(grid = grid, edges = c(grid - step / 2, grid[length(grid)] + step / 2), kind = "uniform")
}

#' Sum a count vector over a window of `h` steps either side of each index.
qc_window_sum <- function(n, h) {
  if (h < 1L) return(n)
  cs <- c(0, cumsum(n)); len <- length(n); i <- seq_len(len)
  cs[pmin(i + h, len) + 1L] - cs[pmax(i - h, 1L)]
}

#' Interior extrema of a curve, ignoring steps smaller than `rel_tol` of the
#' curve's range. Without the tolerance a curve shrunk to a flat line reports
#' every round-off wobble as an extremum.
qc_n_extrema <- function(f, rel_tol = 1e-3) {
  rg <- max(f) - min(f)
  if (!is.finite(rg) || rg <= 0) return(0L)
  s <- sign(diff(f)); s[abs(diff(f)) < rel_tol * rg] <- 0
  s <- s[s != 0]
  if (length(s) > 1L) sum(diff(s) != 0) else 0L
}

#' Monotone direction of a curve over a covariate: +1, -1, or 0 when the
#' Spearman correlation is below `rho_min` in magnitude (non-monotone) or
#' undefined (flat).
qc_dir_of <- function(x, f, rho_min) {
  if (!(stats::sd(f) > 0) || !(stats::sd(x) > 0)) return(0)
  r <- stats::cor(x, f, method = "spearman")
  if (is.na(r) || abs(r) < rho_min) 0 else sign(r)
}

#' A memoising frame builder. The three primary models of a paired signal have
#' three column sets, and the two interaction models have exactly `full`'s
#' (a `ti()` term introduces no variable), so the five frames of one signal
#' are three distinct objects. Built once per (signal, column set) and handed
#' back with THIS model's own formula attached, so `attr(d, "formula")` is
#' always the formula that would have been fitted.
#'
#' The cache is an environment passed by the caller, so its lifetime is the
#' caller's choice: one pass over the training set, or one bag.
qc_frame_for <- function(cache, signal, model, tabs, cfg, pri, stay_ids) {
  f   <- build_formula(signal, model, cfg)
  key <- paste(signal, paste(sort(required_columns(f)), collapse = ","), sep = "|")
  d <- cache[[key]]
  if (is.null(d)) {
    d <- signal_frame(signal, model, tabs, cfg, pri, stay_ids = stay_ids)
    cache[[key]] <- d
  }
  attr(d, "model")   <- model
  attr(d, "formula") <- f
  d
}
