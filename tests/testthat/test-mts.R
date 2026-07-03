# Multiple-time-stepping integrator.

test_that("with no slow force, MTS reduces to fine-grained leapfrog", {
  # step_mts with slow force 0 is m inner leapfrog substeps of size eps/m, so
  # n_steps outer steps equal n_steps*m leapfrog steps at eps/m.
  grad_fast <- function(q) q          # unit oscillator
  grad_zero <- function(q) rep(0, length(q))
  q0 <- c(1.3, -0.4)
  p0 <- c(0.2, 0.9)
  eps <- 0.3
  n_steps <- 25
  m <- 8

  mts <- simp_integrate_mts(q0, p0, grad_slow = grad_zero, grad_fast = grad_fast,
                            eps = eps, n_steps = n_steps, m = m)
  lf <- simp_integrate(q0, p0, grad_U = grad_fast, eps = eps / m,
                       n_steps = n_steps * m, integrator = "leapfrog")

  expect_equal(mts$q, lf$q, tolerance = 1e-10)
  expect_equal(mts$p, lf$p, tolerance = 1e-10)
})

test_that("MTS is time reversible", {
  grad_slow <- function(q) 0.1 * q
  grad_fast <- function(q) 40 * q
  q0 <- c(0.7, -1.1)
  p0 <- c(-0.3, 0.5)

  fwd <- simp_integrate_mts(q0, p0, grad_slow, grad_fast, eps = 0.2,
                            n_steps = 40, m = 12)
  back <- simp_integrate_mts(fwd$q, -fwd$p, grad_slow, grad_fast, eps = 0.2,
                             n_steps = 40, m = 12)

  expect_equal(back$q, q0, tolerance = 1e-9)
  expect_equal(-back$p, p0, tolerance = 1e-9)
})

test_that("MTS keeps the energy bounded on a stiff split Gaussian", {
  # U = 0.5 * (k_fast q1^2 + k_slow q2^2); the stiff coordinate rides the inner
  # loop, the soft one the outer step.
  k_fast <- 100
  k_slow <- 1
  grad_fast <- function(q) c(k_fast * q[1], 0)
  grad_slow <- function(q) c(0, k_slow * q[2])
  U <- function(q) 0.5 * (k_fast * q[1]^2 + k_slow * q[2]^2)

  out <- simp_integrate_mts(q0 = c(1, 1), p0 = c(0, 0), grad_slow, grad_fast,
                            eps = 0.25, n_steps = 200, m = 20,
                            trajectory = TRUE, U = U)
  drift <- max(out$energy) - min(out$energy)
  expect_lt(drift / out$energy[1], 0.05)
})

test_that("MTS validates its arguments", {
  g <- function(q) q
  expect_error(simp_integrate_mts(1, 0, g, g, eps = 0.1, n_steps = 1, m = 0),
               "m must be a positive integer")
  expect_error(
    simp_integrate_mts(c(1, 0), 0, g, g, eps = 0.1, n_steps = 1, m = 2),
    "same length")
})
