## Profiled phase LOD (M02, submission_notes/M02_phase_lod_audit.md).
##
## For the single-family two-point model the repulsion table is the coupling
## table with r -> 1 - r, each phase-specific log-likelihood is concave on (0,1),
## and therefore the losing phase's constrained maximum on [1e-6, 0.5] is the
## no-linkage null whenever the winner has r_hat < 0.5. The profiled phase LOD
## then equals the linkage LOD exactly. These tests pin that behaviour down on
## the production entry point and on the single-pair harnesses, and check the
## general (multi-dam) profiling against an independent R reference.

.pl_pw <- function(G_list, M_list, ...) {
  HalfSibMap:::pairwise_rf_estimation_multi_parallel_cpp(G_list, M_list, ...)
}
## exact mirror of joint_child_probs_3x3(); phase 0 = coupling, 1 = repulsion
.pl_joint <- function(phase, r, qi, qj) {
  F0 <- function(q) c(1 - q, q, 0); F1 <- function(q) c(0, 1 - q, q)
  A <- outer(F0(qi), F0(qj)) + outer(F1(qi), F1(qj))
  B <- outer(F0(qi), F1(qj)) + outer(F1(qi), F0(qj))
  if (phase == 0) (1 - r) / 2 * A + r / 2 * B else (1 - r) / 2 * B + r / 2 * A
}
.pl_ll <- function(C3, P, tiny = 1e-12) sum(ifelse(C3 > 0, C3 * log(pmax(P, tiny)), 0))
.pl_counts <- function(n, r, phase, qi, qj, seed) {
  set.seed(seed)
  nonrec <- if (phase == "C") list(c(0L, 0L), c(1L, 1L)) else list(c(1L, 0L), c(0L, 1L))
  rec    <- if (phase == "C") list(c(0L, 1L), c(1L, 0L)) else list(c(0L, 0L), c(1L, 1L))
  C3 <- matrix(0L, 3, 3)
  for (k in seq_len(n)) {
    m <- if (stats::runif(1) < 1 - r) nonrec[[sample.int(2, 1)]] else rec[[sample.int(2, 1)]]
    yi <- m[1] + stats::rbinom(1, 1, qi); yj <- m[2] + stats::rbinom(1, 1, qj)
    C3[yi + 1L, yj + 1L] <- C3[yi + 1L, yj + 1L] + 1L
  }
  C3
}

# 1 ---------------------------------------------------------------------------
test_that("P_R(r) = P_C(1 - r) cell by cell for engineered single-dam tables", {
  set.seed(202)
  for (k in 1:20) {
    qi <- stats::runif(1, 1e-6, 1 - 1e-6); qj <- stats::runif(1, 1e-6, 1 - 1e-6)
    r  <- stats::runif(1, 1e-6, 0.5)
    expect_lt(max(abs(.pl_joint(1, r, qi, qj) - .pl_joint(0, 1 - r, qi, qj))), 1e-15)
  }
  ## and the identity is what the C++ builders implement: the Newton harness
  ## reports the per-phase maxima; the losing phase must sit at r = 0.5 and its
  ## maximum must equal the null likelihood whenever the winner has r_hat < 0.5
  for (k in 1:12) {
    qi <- stats::runif(1, 0.1, 0.9); qj <- stats::runif(1, 0.1, 0.9)
    C3 <- .pl_counts(200L, stats::runif(1, 0.01, 0.35), sample(c("C", "R"), 1), qi, qj, 300 + k)
    n <- HalfSibMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), qi, qj)
    if (n$r_hat >= 0.5 - 1e-6) next
    r_lose  <- if (n$phase == 1L) n$r_R  else n$r_C
    ll_lose <- if (n$phase == 1L) n$ll_R else n$ll_C
    expect_identical(r_lose, 0.5)
    expect_lt(abs(ll_lose - n$ll_half), 1e-9)
  }
})

# 2 ---------------------------------------------------------------------------
test_that("single family: profiled phase support equals the linkage LOD exactly", {
  RcppParallel::setThreadOptions(numThreads = 1)
  sim <- sim_multi_pop(T_markers = 30L, n_pops = 1L, n_ind_per_pop = 250L,
                       r_vec = c(rep(0.02, 10), rep(0.1, 9), rep(0.3, 10)),
                       phase_mode = "random", repulsion_rate = 0.4,
                       maternal_geno_mode = "all_het", paternal_pA_base = 0.35,
                       error_rate = 0.01, seed = 91)
  res <- .pl_pw(sim$G_list, sim$M_list, lambda = 0, q0 = 0.5, optimizer = "auto")
  expect_identical(res$optimizer, "newton-concave")
  ut <- upper.tri(res$r)
  called <- ut & !is.na(res$mom_phase_list[[1]])
  expect_gt(sum(called), 300L)
  ## assigned, hence bit-identical, wherever a phase is called ...
  expect_identical(res$lod_ph[called], res$lod_r[called])
  expect_identical(res$lod_ph_list[[1]][called], res$lod_r[called])
  ## ... and never above the linkage LOD anywhere
  expect_true(all(res$lod_ph[ut] <= pmax(res$lod_r[ut], 0) + 1e-12, na.rm = TRUE))
  expect_true(all(res$lod_ph[ut] >= 0, na.rm = TRUE))
  ## the same holds on the grid reference optimizer
  rg <- .pl_pw(sim$G_list, sim$M_list, lambda = 0, q0 = 0.5, optimizer = "grid")
  cg <- ut & !is.na(rg$mom_phase_list[[1]])
  expect_identical(rg$lod_ph[cg], rg$lod_r[cg])
  ## and it agrees with an explicit per-phase profile computed in R
  G <- sim$G_list[[1]]; q <- res$q_list[[1]]
  set.seed(5); pick <- which(called, arr.ind = TRUE)
  pick <- pick[sample.int(nrow(pick), 25L), , drop = FALSE]
  for (p in seq_len(nrow(pick))) {
    i <- pick[p, 1]; j <- pick[p, 2]
    yi <- G[, i]; yj <- G[, j]; ok <- !is.na(yi) & !is.na(yj)
    C3 <- matrix(tabulate(3L * yi[ok] + yj[ok] + 1L, 9L), 3, 3, byrow = TRUE)
    prof <- function(ph) stats::optimize(function(r) .pl_ll(C3, .pl_joint(ph, r, q[i], q[j])),
                                         c(1e-6, 0.5), maximum = TRUE, tol = 1e-12)$objective
    ## absolute tolerance: near-null pairs have LOD ~ 0.01, below which a relative
    ## tolerance would be tighter than the optimizer tolerance (1e-6 in r)
    expect_lt(abs(res$lod_ph[i, j] - abs(prof(0) - prof(1)) / log(10)), 1e-6)
  }
})

# 3 ---------------------------------------------------------------------------
test_that("r_hat = 0.5 gives phase NA and phase support 0; exact ties likewise", {
  ## bit-exact null (same construction as test-two-point-optimizer-switch.R)
  G <- rbind(matrix(rep(c(0L, 0L), each = 25), 25, 2), matrix(rep(c(2L, 2L), each = 25), 25, 2),
             matrix(rep(c(0L, 2L), each = 25), 25, 2), matrix(rep(c(2L, 0L), each = 25), 25, 2))
  colnames(G) <- c("m1", "m2")
  res <- .pl_pw(list(P1 = G), list(P1 = c(m1 = 1L, m2 = 1L)), lambda = 0, q0 = 0.5)
  expect_identical(res$r["m1", "m2"], 0.5)
  expect_true(is.na(res$mom_phase_list$P1["m1", "m2"]))
  expect_identical(res$lod_ph["m1", "m2"], 0)
  expect_identical(res$lod_ph_list$P1["m1", "m2"], 0)
  ## flat likelihood (only double-het offspring, q = 0.5): tie at r_hat < 0.5
  C3 <- matrix(0L, 3, 3); C3[2, 2] <- 200L
  n <- HalfSibMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), 0.5, 0.5)
  expect_true(is.na(n$phase)); expect_identical(n$lod_ph, 0)
})

# 4 ---------------------------------------------------------------------------
test_that("the preferred phase sign is the phase with the larger profiled maximum", {
  set.seed(77)
  for (k in 1:40) {
    qi <- stats::runif(1, 0.05, 0.95); qj <- stats::runif(1, 0.05, 0.95)
    C3 <- .pl_counts(sample(c(15L, 60L, 300L), 1), stats::runif(1, 0.005, 0.45),
                     sample(c("C", "R"), 1), qi, qj, 700 + k)
    n <- HalfSibMap:::two_point_pair_newton_cpp(C3, integer(3), integer(3), qi, qj)
    if (is.na(n$phase)) next
    ## sign from the shared-r_hat comparison (as implemented) ...
    ## ... must agree with the sign of the per-phase profiled maxima
    expect_identical(n$phase, if (n$ll_C > n$ll_R) 1L else 0L)
    ## and the reported support is the profiled contrast = lod_r
    expect_lt(abs(n$lod_ph - abs(n$ll_C - n$ll_R) / log(10)), 1e-9)
    expect_identical(n$lod_ph, max(n$lod_r, 0))
  }
})

# 5 ---------------------------------------------------------------------------
test_that("several dams: per-dam support is the profiled ratio against an R reference", {
  RcppParallel::setThreadOptions(numThreads = 1)
  M <- list(A = c(m1 = 1L, m2 = 1L), B = c(m1 = 1L, m2 = 1L))
  mk <- function(C3) {          # expand a 3x3 count table into a genotype matrix
    idx <- which(C3 > 0, arr.ind = TRUE)
    g <- do.call(rbind, lapply(seq_len(nrow(idx)), function(k)
      matrix(rep(c(idx[k, 1] - 1L, idx[k, 2] - 1L), C3[idx[k, 1], idx[k, 2]]), ncol = 2, byrow = TRUE)))
    colnames(g) <- c("m1", "m2"); storage.mode(g) <- "integer"; g
  }
  CA <- .pl_counts(120L, 0.08, "C", 0.4, 0.6, 901)     # dam A coupling, tight
  CB <- .pl_counts(60L,  0.25, "R", 0.5, 0.5, 902)     # dam B repulsion, loose
  res <- .pl_pw(list(A = mk(CA), B = mk(CB)), M, lambda = 0, q0 = 0.5)
  expect_identical(res$optimizer, "grid+local-refine")
  qA <- res$q_list$A; qB <- res$q_list$B
  llA <- function(ph, r) .pl_ll(CA, .pl_joint(ph, r, qA[1], qA[2]))
  llB <- function(ph, r) .pl_ll(CB, .pl_joint(ph, r, qB[1], qB[2]))
  pooled <- function(r) max(llA(0, r), llA(1, r)) + max(llB(0, r), llB(1, r))
  rs <- seq(1e-6, 0.5, length.out = 20001L)
  ll_hat <- max(vapply(rs, pooled, 1))
  expect_lt(abs(res$logLik[1, 2] - ll_hat), 1e-4)
  callA <- res$mom_phase_list$A[1, 2]; callB <- res$mom_phase_list$B[1, 2]
  expect_identical(callA, 1L); expect_identical(callB, 0L)
  ## alternative for dam A: A forced to the other phase, B and r free
  altA <- function(r) llA(if (callA == 1L) 1 else 0, r) + max(llB(0, r), llB(1, r))
  altB <- function(r) max(llA(0, r), llA(1, r)) + llB(if (callB == 1L) 1 else 0, r)
  refA <- (ll_hat - max(vapply(rs, altA, 1))) / log(10)
  refB <- (ll_hat - max(vapply(rs, altB, 1))) / log(10)
  expect_equal(res$lod_ph_list$A[1, 2], refA, tolerance = 1e-3)
  expect_equal(res$lod_ph_list$B[1, 2], refB, tolerance = 1e-3)
  expect_equal(res$lod_ph[1, 2], refA + refB, tolerance = 2e-3)
  ## each is bounded by the pooled linkage LOD and by the old shared-r contrast
  expect_lte(res$lod_ph_list$A[1, 2], res$lod_r[1, 2] + 1e-9)
  expect_lte(res$lod_ph_list$B[1, 2], res$lod_r[1, 2] + 1e-9)
  r_hat <- res$r[1, 2]
  expect_lt(res$lod_ph_list$A[1, 2], abs(llA(0, r_hat) - llA(1, r_hat)) / log(10))
})
