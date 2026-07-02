# Harmonic oscillator U(q) = 0.5 q^2, exact flow q(t) = q0 cos t + p0 sin t,
# p(t) = -q0 sin t + p0 cos t. A fixed horizon T = 2*pi returns to the start.
osc_grad <- function(q) q

osc_error <- function(integrator, eps, horizon = 2 * pi, q0 = 1, p0 = 0) {
  n <- max(1L, round(horizon / eps))
  eps <- horizon / n
  out <- simp_integrate(q0, p0, osc_grad, eps, n, integrator = integrator)
  q_true <- q0 * cos(horizon) + p0 * sin(horizon)
  p_true <- -q0 * sin(horizon) + p0 * cos(horizon)
  sqrt((out$q - q_true)^2 + (out$p - p_true)^2)
}

test_that("each scheme converges at its declared order", {
  # Halving eps must divide the global error by about 2^order. Base steps are
  # placed in the asymptotic regime and stopped before the error reaches the
  # rounding floor, where the ratio would stop being informative. A wrong
  # coefficient lowers the observed order, so the ratio catches it.
  cases <- list(
    leapfrog = list(order = 2, base = 0.08),
    yoshida4 = list(order = 4, base = 0.125),
    yoshida6 = list(order = 6, base = 0.18),
    yoshida8 = list(order = 8, base = 0.30)
  )
  for (nm in names(cases)) {
    ord <- cases[[nm]]$order
    epss <- cases[[nm]]$base / 2^(0:2)
    errs <- vapply(epss, function(e) osc_error(nm, e), numeric(1))
    expect_true(all(diff(errs) < 0),
                info = paste(nm, "error must fall as eps shrinks"))
    expect_true(all(errs > 1e-11),
                info = paste(nm, "errors must stay above the rounding floor"))
    log2ratio <- log2(errs[-length(errs)] / errs[-1])
    expect_true(all(abs(log2ratio - ord) < 0.6),
                info = paste(nm, "observed order:",
                             paste(round(log2ratio, 2), collapse = ", ")))
  }
})

test_that("minerror2 cancels the leading energy error on a Gaussian target", {
  U <- function(q) 0.5 * sum(q^2)
  energy_range <- function(scheme, eps, nstep) {
    out <- simp_integrate(1.0, 0.0, osc_grad, eps, nstep, integrator = scheme,
                          trajectory = TRUE, U = U)
    max(out$energy) - min(out$energy)
  }
  # On the harmonic oscillator the leading (h^2) energy error is cancelled, so
  # the residual scales like h^4 -- halving eps divides it by about 16.
  e1 <- energy_range("minerror2", 0.20, 200)
  e2 <- energy_range("minerror2", 0.10, 400)
  expect_gt(log2(e1 / e2), 3.4)
  # And it is far below leapfrog's error at the same step size.
  expect_lt(e1, energy_range("leapfrog", 0.20, 200) / 50)
})

test_that("minerror2 is order two on a nonlinear problem", {
  # Anharmonic U = q^2/2 + q^4/4, no closed form: use self-convergence, where
  # halving eps divides the endpoint change by about 2^order.
  grad <- function(q) q + q^3
  horizon <- 3.0
  endpoint <- function(eps) {
    nn <- round(horizon / eps)
    out <- simp_integrate(0.8, 0.2, grad, horizon / nn, nn,
                          integrator = "minerror2")
    c(out$q, out$p)
  }
  e <- lapply(c(0.08, 0.04, 0.02, 0.01), endpoint)
  r1 <- log2(sqrt(sum((e[[1]] - e[[2]])^2)) / sqrt(sum((e[[2]] - e[[3]])^2)))
  r2 <- log2(sqrt(sum((e[[2]] - e[[3]])^2)) / sqrt(sum((e[[3]] - e[[4]])^2)))
  expect_true(abs(r1 - 2) < 0.6, info = paste("r1 =", round(r1, 2)))
  expect_true(abs(r2 - 2) < 0.6, info = paste("r2 =", round(r2, 2)))
})

test_that("integration is time reversible", {
  grad_U <- function(q) q^3  # nonlinear, U = 0.25 q^4
  q0 <- c(0.7, -0.3)
  p0 <- c(0.1, 0.5)
  for (nm in simp_schemes()$integrator) {
    fwd <- simp_integrate(q0, p0, grad_U, eps = 0.05, n_steps = 40,
                          integrator = nm)
    back <- simp_integrate(fwd$q, -fwd$p, grad_U, eps = 0.05, n_steps = 40,
                           integrator = nm)
    expect_equal(back$q, q0, tolerance = 1e-9, info = nm)
    expect_equal(back$p, -p0, tolerance = 1e-9, info = nm)
  }
})

test_that("total energy stays bounded with no secular drift", {
  # Anharmonic 1-D system, long run.
  U <- function(q) 0.5 * q^2 + 0.25 * q^4
  grad_U <- function(q) q + q^3
  n <- 4000
  out <- simp_integrate(1.0, 0.0, grad_U, eps = 0.05, n_steps = n,
                        integrator = "leapfrog", trajectory = TRUE, U = U)
  e <- out$energy
  rng <- max(e) - min(e)
  expect_lt(rng, 1e-2)  # bounded oscillation, not growth
  half <- floor(length(e) / 2)
  drift <- abs(mean(e[(half + 1):length(e)]) - mean(e[1:half]))
  expect_lt(drift, rng)  # no linear energy drift between halves
})

test_that("one step is symplectic (Jacobian preserves the form)", {
  grad_U <- function(q) q^3
  eps <- 0.1
  one_step <- function(z) {
    out <- simp_integrate(z[1], z[2], grad_U, eps = eps, n_steps = 1,
                          integrator = "leapfrog")
    c(out$q, out$p)
  }
  z0 <- c(0.6, -0.2)
  h <- 1e-6
  J <- matrix(0, 2, 2)
  for (k in 1:2) {
    zp <- z0; zp[k] <- zp[k] + h
    zm <- z0; zm[k] <- zm[k] - h
    J[, k] <- (one_step(zp) - one_step(zm)) / (2 * h)
  }
  Omega <- matrix(c(0, -1, 1, 0), 2, 2)  # [[0,1],[-1,0]] column-major
  expect_equal(det(J), 1, tolerance = 1e-6)
  expect_equal(t(J) %*% Omega %*% J, Omega, tolerance = 1e-5)
})
