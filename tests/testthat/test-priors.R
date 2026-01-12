# Tests for prior specification

test_that("prior_normal creates valid object", {
  p <- prior_normal(0, 2.5)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "normal")
  expect_equal(p$parameters$mean, 0)
  expect_equal(p$parameters$sd, 2.5)
})

test_that("prior_normal rejects invalid sd", {
  expect_error(prior_normal(0, -1), "must be positive")
  expect_error(prior_normal(0, 0), "must be positive")
})

test_that("prior_half_normal creates valid object", {
  p <- prior_half_normal(2)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "half_normal")
  expect_equal(p$parameters$sd, 2)
})

test_that("prior_half_cauchy creates valid object", {
  p <- prior_half_cauchy(0, 2)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "half_cauchy")
  expect_equal(p$parameters$location, 0)
  expect_equal(p$parameters$scale, 2)
})

test_that("prior_gamma creates valid object", {
  p <- prior_gamma(2, 0.1)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "gamma")
  expect_equal(p$parameters$shape, 2)
  expect_equal(p$parameters$rate, 0.1)
})

test_that("prior_exponential creates valid object", {
  p <- prior_exponential(0.5)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "exponential")
  expect_equal(p$parameters$rate, 0.5)
})

test_that("prior_beta creates valid object", {
  p <- prior_beta(2, 2)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "beta")
  expect_equal(p$parameters$alpha, 2)
  expect_equal(p$parameters$beta, 2)
})

test_that("prior_pc creates valid object with lambda", {
  p <- prior_pc(u = 1, alpha = 0.01)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "pc")
  expect_equal(p$parameters$u, 1)
  expect_equal(p$parameters$alpha, 0.01)
  # lambda = -log(0.01) / 1 ≈ 4.605
  expect_equal(p$parameters$lambda, -log(0.01) / 1, tolerance = 1e-10)
})

test_that("prior_pc requires both arguments", {
  expect_error(prior_pc(u = 1), "required")
  expect_error(prior_pc(alpha = 0.01), "required")
})

test_that("prior_lkj creates valid object", {
  p <- prior_lkj(2)

  expect_s3_class(p, "simp_prior")
  expect_equal(p$distribution, "lkj")
  expect_equal(p$parameters$eta, 2)
})

test_that("simp_priors creates valid container", {
  priors <- simp_priors()

  expect_s3_class(priors, "simp_priors")
  expect_true(all(c("beta", "intercept", "sigma", "phi", "cor") %in% names(priors)))
})

test_that("simp_priors accepts custom priors", {
  priors <- simp_priors(
    beta = prior_normal(0, 10),
    sigma = prior_half_normal(2)
  )

  expect_equal(priors$beta$parameters$sd, 10)
  expect_equal(priors$sigma$distribution, "half_normal")
})

test_that("simp_priors rejects non-prior objects", {
  expect_error(simp_priors(beta = "wrong"), "must be a")
})

test_that("default_priors returns valid priors for families", {
  priors_dir <- default_priors("dirichlet")
  expect_s3_class(priors_dir, "simp_priors")

  priors_ln <- default_priors("logistic_normal")
  expect_s3_class(priors_ln, "simp_priors")

  priors_hd <- default_priors("hurdle_dirichlet")
  expect_s3_class(priors_hd, "simp_priors")
})

test_that("is_prior and is_priors work", {
  expect_true(is_prior(prior_normal(0, 1)))
  expect_false(is_prior("not a prior"))

  expect_true(is_priors(simp_priors()))
  expect_false(is_priors(list(a = 1)))
})
