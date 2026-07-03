# Exact harmonic energy-error metric and the step-adapted minimum-error schemes.

# Worst-case relative energy error of a named/adaptive scheme over (0, nu_max].
band_error <- function(ops, nu_max, n = 64) {
  nu <- nu_max * seq_len(n) / n
  max(vapply(nu, function(x) SIMP:::simp_harmonic_error_cpp(ops$kind, ops$coef, x),
             numeric(1)))
}

test_that("leapfrog's harmonic stability limit is nu = 2", {
  lf <- SIMP:::simp_scheme_ops_cpp("leapfrog")
  expect_equal(SIMP:::simp_stability_limit_cpp(lf$kind, lf$coef, 10), 2,
               tolerance = 1e-6)
})

test_that("energy error vanishes as the step shrinks and grows near the cliff", {
  lf <- SIMP:::simp_scheme_ops_cpp("leapfrog")
  expect_lt(SIMP:::simp_harmonic_error_cpp(lf$kind, lf$coef, 0.05), 1e-2)
  expect_gt(SIMP:::simp_harmonic_error_cpp(lf$kind, lf$coef, 1.9),
            SIMP:::simp_harmonic_error_cpp(lf$kind, lf$coef, 1.0))
  expect_true(is.infinite(
    SIMP:::simp_harmonic_error_cpp(lf$kind, lf$coef, 2.5)))  # beyond stability
})

test_that("minerror2 beats leapfrog on the harmonic energy error", {
  lf <- SIMP:::simp_scheme_ops_cpp("leapfrog")
  me <- SIMP:::simp_scheme_ops_cpp("minerror2")
  # Match per-gradient cost: leapfrog at half the step is one gradient/step like
  # minerror2 is two/step, so compare minerror2(nu) with leapfrog(nu/2).
  for (nu in c(0.5, 1.0, 1.5)) {
    e_me <- SIMP:::simp_harmonic_error_cpp(me$kind, me$coef, nu)
    e_lf <- SIMP:::simp_harmonic_error_cpp(lf$kind, lf$coef, nu / 2)
    expect_lt(e_me, e_lf)
  }
})

test_that("the adaptive coefficient minimizes the band error it targets", {
  nu_max <- 1.4
  b <- simp_b_optimal(nu_max, stages = 2)
  expect_true(b > 0 && b < 0.5)

  ad <- SIMP:::simp_adaptive_ops_cpp(nu_max, 2L)
  me <- SIMP:::simp_scheme_ops_cpp("minerror2")
  lf <- SIMP:::simp_scheme_ops_cpp("leapfrog")

  # Tuned to this band, adaptive2 is no worse than the fixed minerror2 and both
  # beat leapfrog there.
  expect_lte(band_error(ad, nu_max), band_error(me, nu_max) + 1e-8)
  expect_lt(band_error(ad, nu_max), band_error(lf, nu_max))
})

test_that("adaptive3 holds a small error over a wider band than adaptive2", {
  nu_max <- 2.5  # a demanding band where two stages struggle
  ad2 <- SIMP:::simp_adaptive_ops_cpp(nu_max, 2L)
  ad3 <- SIMP:::simp_adaptive_ops_cpp(nu_max, 3L)
  expect_lt(band_error(ad3, nu_max), band_error(ad2, nu_max))

  co <- simp_b_optimal(nu_max, stages = 3)
  expect_length(co, 2)
  expect_true(all(co > 0 & co < 0.5))
})

test_that("adaptive integrators run through simp_integrate and simp_hmc", {
  grad_U <- function(q) q
  out <- simp_integrate(c(1, 0), c(0, 1), grad_U, eps = 0.4, n_steps = 30,
                        integrator = "adaptive2", nu_max = 1.2)
  expect_length(out$q, 2)

  Sigma <- matrix(c(1, 0.6, 0.6, 1), 2)
  P <- solve(Sigma)
  U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
  gU <- function(q) as.numeric(P %*% q)
  set.seed(1)
  fit <- simp_hmc(U, gU, init = c(0, 0), n_iter = 300, eps = 0.3, n_steps = 15,
                  integrator = "adaptive2", nu_max = 1.5, warmup = 150)
  expect_gt(fit$accept_rate, 0.6)
  expect_lt(max(abs(colMeans(fit$draws))), 0.4)
})

test_that("adaptive integrators require nu_max", {
  grad_U <- function(q) q
  expect_error(
    simp_integrate(1, 0, grad_U, eps = 0.1, n_steps = 1, integrator = "adaptive2"),
    "nu_max must be supplied")
})
