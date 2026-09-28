# tests/attr_eval_unit.R -------------------------------------------------------
# SYNTHETIC UNIT TESTS FOR THE ATTRIBUTION-EVALUATION LIBRARY. Reads no data,
# no targets store and no replicate; every input is built in this file.
#
# Written 2026-09-09 against `docs/attribution_eval_review_20260909.md`. Each
# block reproduces one of that review's synthetic counterexamples and asserts
# the repaired behaviour, so a regression of any of them fails here in seconds
# rather than in a seventeen-minute consumer run or a five-hour generator run.
#
#   A1  the store fingerprint moves when the design key does not
#   A2  the posterior draw layout is nd-dependent (documented, versioned)
#   A3  the bitmask fast path equals the tie-aware definition, ties included
#   A6  domain-level displacement is not the signal-level displacement
#   A7  the one-pass delta grid equals `attr_agree_k()` at every delta
#   A8  share metrics are oriented before dominance; both populations exist
#   A9  duplicate coordinates and reordered rows are refused
#   A11 a refit of the same bam formula is bitwise identical, so a memoised
#       alias fit reproduces the refit path exactly
#
#   Rscript tests/attr_eval_unit.R
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv); library(yaml); library(digest)
})
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)

FAILS <- 0L
check <- function(what, ok) {
  ok <- isTRUE(ok)
  cat(sprintf("  %-70s %s\n", what, if (ok) "PASS" else "*** FAIL ***"))
  if (!ok) FAILS <<- FAILS + 1L
  invisible(ok)
}
expect_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

cat("\n=== attribution-evaluation library: synthetic unit tests ===\n\n")

# --- A3: tie-aware top-k on the bitmask path ---------------------------------
cat("A3  top-k agreement semantics\n")
A <- rbind(c(3, 2, 1, 1), c(5, 4, 3, 2), c(0, 0, 2, 0))
B <- rbind(c(1, 1, 3, 2), c(5, 4, 3, 2), c(2, 0, 0, 0))
ta <- attr_topk(A, c(1L, 3L)); tb <- attr_topk(B, c(1L, 3L))
pa <- attr_prep(A); pb <- attr_prep(B)
check("review counterexample: bitmask top-3 == tie-aware (row 1 agrees)",
      attr_topk_agree(ta, tb, 3L) == mean(attr_agree_k(pa, pb, 3L, 0)))
check("review counterexample: row 1 top-3 is an agreement under ties",
      isTRUE(attr_agree_k(pa, pb, 3L, 0)[1]))
set.seed(11)
for (trial in 1:20) {
  n <- 400; p <- sample(c(5L, 11L, 19L), 1)
  # Coarse values so boundary ties are common, plus assigned zeros.
  X <- matrix(sample(0:3, n * p, replace = TRUE), n, p) * sample(c(-1, 1), n * p, TRUE)
  Y <- X; sw <- sample(length(X), n); Y[sw] <- sample(0:3, n, TRUE)
  tx <- attr_topk(X, c(1L, 3L)); ty <- attr_topk(Y, c(1L, 3L))
  px <- attr_prep(X); py <- attr_prep(Y)
  for (k in c(1L, 3L)) {
    if (abs(attr_topk_agree(tx, ty, k) - mean(attr_agree_k(px, py, k, 0))) > 1e-12) {
      check(sprintf("random tie-heavy matrices, trial %d, k = %d", trial, k), FALSE)
    }
  }
}
check("20 random tie-heavy trials: bitmask == tie-aware at k = 1 and 3", TRUE)

# --- A7: the one-pass delta grid -----------------------------------------------
cat("A7  delta grid\n")
set.seed(5)
X <- matrix(rnorm(1500 * 19), 1500); Y <- X + matrix(rnorm(1500 * 19, sd = 0.4), 1500)
px <- attr_prep(X); py <- attr_prep(Y)
grid <- c(0, 0.02, 0.05, 0.1, 0.25, 0.5)
ok <- TRUE
for (k in c(1L, 3L)) {
  g <- attr_agree_k_grid(px, py, k, grid)
  h <- vapply(grid, function(d) mean(attr_agree_k(px, py, k, d)), numeric(1))
  ok <- ok && max(abs(g - h)) < 1e-12
}
check("absolute grid equals attr_agree_k at every delta, k = 1 and 3", ok)
Z <- X * 2   # same ordering, doubled totals
pz <- attr_prep(Z)
gr <- attr_agree_k_grid(px, pz, 3L, c(0.01, 0.05), relative = TRUE)
check("relative grid is symmetric in the arms", isTRUE(all.equal(
  gr, attr_agree_k_grid(pz, px, 3L, c(0.01, 0.05), relative = TRUE))))

# --- A6: domain calibration is not signal calibration -----------------------
cat("A6  domain against signal displacement\n")
m1 <- matrix(c(1, 1), 1); m2 <- matrix(c(2, 0), 1)
dom <- data.frame(signal = c("a", "b"), domain = c("d", "d"))
colnames(m1) <- colnames(m2) <- c("a", "b")
ds <- attr_displacement(list(m1, m2), max_pairs = 10, cell_frac = 1, seed = 1)
dd <- attr_displacement(list(attr_to_domain(m1, dom), attr_to_domain(m2, dom)),
                        max_pairs = 10, cell_frac = 1, seed = 1)
check("signal displacement median 1, domain displacement median 0",
      stats::median(ds) == 1 && stats::median(dd) == 0)

# --- A8: orientation and populations ------------------------------------------
cat("A8  orientation and leader populations\n")
x <- c(0.8, 0.85, 0.9); y <- c(0.1, 0.15, 0.2)
raw <- attr_dominance(x, y)$p_dominates
ox <- attr_orient(x, "share_ab_median"); oy <- attr_orient(y, "share_ab_median")
check("raw share dominance reads 1; oriented share-loss dominance reads 0",
      raw == 1 && attr_dominance(ox$value, oy$value)$p_dominates == 0)
check("orientation map: share_* is agreement, everything else disagreement",
      attr_metric_orientation("share_ba_p05") == "agreement" &&
      attr_metric_orientation("cosine_dissim") == "disagreement" &&
      attr_metric_orientation("top1") == "disagreement")
S <- matrix(abs(rnorm(200 * 6)), 200); T <- S; T[1:50, ] <- T[1:50, 6:1]
cl <- attr_leader_cells(S, T)
ld <- attr_leader_distribution(cl, shares = 0.02)
check("identical leaders are excluded from the `leader_differs` population",
      all(ld$n_scored[ld$population == "leader_differs"] ==
          sum(cl$ok & !cl$same_leader)) &&
      all(ld$n_scored[ld$population == "all_scored"] == sum(cl$ok)))
pool <- .pool_leader_new()
pool <- .pool_leader_add(pool, cl); pool <- .pool_leader_add(pool, cl)
ld2 <- attr_leader_distribution(pool, shares = 0.02)
check("a pooled cell is labelled patient_pairs and counts both pairs",
      all(ld2$unit == "patient_pairs") && all(ld2$n_pairs == 2L) &&
      all(ld2$n_patients == 400L))
pool_old <- .pool_leader_flatten(.pool_leader_cells(list(cl, cl)))
pool_fl  <- .pool_leader_flatten(pool)
check("incremental pooling equals list pooling",
      identical(pool_old$disp_a, pool_fl$disp_a) && identical(pool_old$sh_a_in_b, pool_fl$sh_a_in_b) &&
      identical(pool_fl$disp_a, c(cl$disp_a[cl$ok], cl$disp_a[cl$ok])))
# A pool built from scratch extracts (what the bag-outer loop reads back).
pool_x <- .pool_leader_add(.pool_leader_add(.pool_leader_new(), .pool_leader_extract(cl)), .pool_leader_extract(cl))
check("pooling compact extracts equals pooling full cells",
      identical(attr_leader_distribution(pool_x, 0.02), ld2))

# --- A9: manifest and matrix validation ----------------------------------------
cat("A9  manifest and replicate validation\n")
man <- data.frame(method = "llr_full", route = "bootstrap", boot_id = c(1L, 1L),
                  seed_id = 0L, draw_id = 0L, key = c("k1", "k2"),
                  n_rows = 3L, n_cols = 2L, stringsAsFactors = FALSE)
check("duplicate coordinate is refused", expect_error(attr_validate_manifest(man, "t")))
man$boot_id <- c(1L, 2L)
check("distinct coordinates pass", isTRUE(attr_validate_manifest(man, "t")))
man$n_rows <- c(3L, 4L)
check("mixed shapes are refused", expect_error(attr_validate_manifest(man, "t")))
M <- matrix(1:6 + 0, 3, 2, dimnames = list(c("s1", "s2", "s3"), c("a", "b")))
check("a conforming matrix passes",
      isTRUE(attr_check_replicate(M, c("s1", "s2", "s3"), c("a", "b"))))
check("a reordered matrix is refused",
      expect_error(attr_check_replicate(M[c(2, 1, 3), ], c("s1", "s2", "s3"), c("a", "b"))))
M2 <- M; M2[1, 1] <- NA
check("a non-finite cell is refused",
      expect_error(attr_check_replicate(M2, c("s1", "s2", "s3"), c("a", "b"))))
plan <- data.frame(method = "llr_full", route = "bootstrap", index = 501:503,
                   boot_id = 1:3, seed_id = 0L, draw_id = 0L, fits = NA_integer_,
                   note = "", stringsAsFactors = FALSE)
man2 <- data.frame(method = "llr_full", route = "bootstrap", boot_id = 1L,
                   seed_id = 0L, draw_id = 0L, key = "k1", stringsAsFactors = FALSE)
cov <- attr_replicate_coverage(plan, man2, tombstoned = 2L)
check("coverage: 3 planned, 1 present, 1 tombstoned, 1 missing, not complete",
      cov$planned == 3L && cov$present == 1L && cov$tombstoned == 1L &&
      cov$missing == 1L && !cov$complete)

# --- A1: the fingerprint moves when the key does not ---------------------------
cat("A1  store fingerprint\n")
cfg0 <- list(seed = 1L, signals = list("a"), xgboost = list(max_depth = 4L, eta = 0.05),
             intervention_agent_pool = list(vasopressor = c("x", "y")))
ecfg0 <- list(bootstrap = list(seed_base = 1L),
              levels = list(seed = list(seed_base = 1000L),
                            sample = list(draw_seed_base = 1L, b = 40L,
                                          llr_route = "posterior")))
dk <- list(bam = "x", folds = "f")
f0 <- attr_design_fingerprint(cfg0, NULL, 5L, priors = list(p = 1), ecfg0, design_key = dk)
cfg1 <- cfg0; cfg1$xgboost$max_depth <- 6L
f1 <- attr_design_fingerprint(cfg1, NULL, 5L, priors = list(p = 1), ecfg0, design_key = dk)
d1 <- attr_fingerprint_diff(f0, f1)
check("xgboost.max_depth moves the fingerprint, names xgboost and cfg_design",
      setequal(d1$field, c("xgboost", "cfg_design")))
cfg2 <- cfg0; cfg2$intervention_agent_pool$vasopressor <- c("x", "y", "z")
d2 <- attr_fingerprint_diff(f0, attr_design_fingerprint(cfg2, NULL, 5L, list(p = 1), ecfg0, design_key = dk))
check("a wider vasopressor agent pool moves cfg_design", "cfg_design" %in% d2$field)
d3 <- attr_fingerprint_diff(f0, attr_design_fingerprint(cfg0, NULL, 5L, list(p = 2), ecfg0, design_key = dk))
check("a changed priors object moves priors", identical(d3$field, "priors"))
ecfg1 <- ecfg0; ecfg1$levels$sample$b <- 60L
d4 <- attr_fingerprint_diff(f0, attr_design_fingerprint(cfg0, NULL, 5L, list(p = 1), ecfg1, design_key = dk))
check("extending posterior B moves eval_sample_b and names the posterior route",
      identical(d4$field, "eval_sample_b") && d4$routes_affected == "posterior")
check("an unchanged design gives an empty diff",
      nrow(attr_fingerprint_diff(f0, f0)) == 0L)

# --- A2: the draw layout is nd-dependent, and is versioned --------------------
cat("A2  posterior draw layout\n")
set.seed(1); Z6 <- matrix(rnorm(6 * 2), 6, 2)
set.seed(1); Z4 <- matrix(rnorm(4 * 2), 4, 2)
check("column-major layout: draw 1..4 differ between B = 6 and B = 4 (documented)",
      !identical(Z6[1:4, ], Z4))
check("the layout is named in ATTR_GENERATOR_VERSION",
      identical(ATTR_GENERATOR_VERSION$posterior_draw_layout, "rnorm_colmajor_nd_dependent_v1"))

# --- A11: bam determinism, so a memoised alias fit is the refit ---------------
cat("A11 bam refit identity\n")
set.seed(7)
d <- data.frame(x = runif(400), z = runif(400))
d$y <- rbinom(400, 1, plogis(-1 + 2 * d$x))
f <- y ~ s(x, bs = "ts", k = 6) + s(z, bs = "ts", k = 6)
b1 <- mgcv::bam(f, data = d, family = binomial(), method = "fREML", discrete = TRUE,
                nthreads = 2, gamma = 1.5, na.action = na.fail)
b2 <- mgcv::bam(f, data = d, family = binomial(), method = "fREML", discrete = TRUE,
                nthreads = 2, gamma = 1.5, na.action = na.fail)
check("two bam fits of one formula: identical coefficients",
      identical(coef(b1), coef(b2)))
check("two bam fits of one formula: identical Vc (or Vp when Vc is absent)",
      identical(b1$Vc %||% b1$Vp, b2$Vc %||% b2$Vp))
dr1 <- with_seed(99, matrix(rnorm(5 * length(coef(b1))), 5))
dr2 <- with_seed(99, matrix(rnorm(5 * length(coef(b2))), 5))
check("with_seed reproduces the standard-normal draw", identical(dr1, dr2))

cat(sprintf("\n=== %s ===\n\n", if (FAILS == 0L) "ALL PASS" else paste(FAILS, "FAILURE(S)")))
if (FAILS > 0L) quit(save = "no", status = 1L)
