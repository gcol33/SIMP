test_that("HMC recovers a correlated Gaussian target", {
  skip_on_cran()
  set.seed(20260703)

  Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
  P <- solve(Sigma)
  U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
  grad_U <- function(q) as.numeric(P %*% q)

  fit <- simp_hmc(U, grad_U, init = c(0, 0), n_iter = 4000, eps = 0.18,
                  n_steps = 18, integrator = "leapfrog", warmup = 1000)

  expect_gt(fit$accept_rate, 0.6)
  expect_equal(fit$n_divergent, 0L)

  m <- colMeans(fit$draws)
  expect_lt(max(abs(m)), 0.1)

  S <- stats::cov(fit$draws)
  expect_lt(max(abs(S - Sigma)), 0.15)
})

test_that("a higher-order integrator also samples the target", {
  skip_on_cran()
  set.seed(11)

  U <- function(q) 0.5 * sum(q^2)
  grad_U <- function(q) q

  fit <- simp_hmc(U, grad_U, init = c(0, 0, 0), n_iter = 3000, eps = 0.3,
                  n_steps = 8, integrator = "yoshida4", warmup = 800)

  expect_gt(fit$accept_rate, 0.7)
  expect_equal(fit$n_divergent, 0L)
  expect_lt(max(abs(colMeans(fit$draws))), 0.12)
  expect_lt(max(abs(apply(fit$draws, 2, stats::sd) - 1)), 0.15)
})
