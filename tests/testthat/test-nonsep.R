# Generalized (implicit) leapfrog for non-separable Hamiltonians.

test_that("on a separable Hamiltonian it reduces to leapfrog", {
  # H = 0.5 p'p + 0.5 q'q. dH/dq = q (no p), dH/dp = p (no q).
  q0 <- c(1, -0.4); p0 <- c(0.2, 0.7)
  imp <- simp_integrate_nonsep(
    q0, p0,
    dH_dq = function(q, p) q, dH_dp = function(q, p) p,
    eps = 0.1, n_steps = 50)
  lf <- simp_integrate(q0, p0, grad_U = function(q) q,
                       eps = 0.1, n_steps = 50, integrator = "leapfrog")
  expect_true(imp$converged)
  expect_equal(imp$q, lf$q, tolerance = 1e-8)
  expect_equal(imp$p, lf$p, tolerance = 1e-8)
})

test_that("it conserves energy on a non-separable Hamiltonian", {
  # Position-dependent mass m(q) = 1 + 0.5 q^2:
  #   H(q, p) = 0.5 p^2 / m(q) + 0.5 q^2.
  m  <- function(q) 1 + 0.5 * q^2
  H  <- function(q, p) 0.5 * p^2 / m(q) + 0.5 * q^2
  dq <- function(q, p) -0.5 * p^2 * q / m(q)^2 + q   # dH/dq (depends on p)
  dp <- function(q, p) p / m(q)                       # dH/dp (depends on q)

  out <- simp_integrate_nonsep(
    q0 = 1.0, p0 = 0.5, dH_dq = dq, dH_dp = dp,
    eps = 0.02, n_steps = 3000, trajectory = TRUE, H = H)

  expect_true(out$converged)
  rng <- max(out$energy) - min(out$energy)
  expect_lt(rng, 1e-2)  # bounded oscillation, no secular drift
  half <- floor(length(out$energy) / 2)
  drift <- abs(mean(out$energy[(half + 1):length(out$energy)]) -
                 mean(out$energy[1:half]))
  expect_lt(drift, rng)
})

test_that("it is time reversible on a non-separable Hamiltonian", {
  m  <- function(q) 1 + 0.5 * q^2
  dq <- function(q, p) -0.5 * p^2 * q / m(q)^2 + q
  dp <- function(q, p) p / m(q)

  fwd <- simp_integrate_nonsep(1.0, 0.5, dq, dp, eps = 0.05, n_steps = 40)
  back <- simp_integrate_nonsep(fwd$q, -fwd$p, dq, dp, eps = 0.05, n_steps = 40)
  expect_equal(back$q, 1.0, tolerance = 1e-7)
  expect_equal(back$p, -0.5, tolerance = 1e-7)
})
