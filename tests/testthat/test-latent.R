# Tests for latent factor specification

test_that("latent_factor creates valid specification", {
  lf <- latent_factor()

  expect_s3_class(lf, "simp_latent")
  expect_equal(lf$n_factors, 1L)
  expect_s3_class(lf$prior, "simp_prior")
})

test_that("latent_factor accepts custom number of factors", {
  lf <- latent_factor(n = 3)

  expect_equal(lf$n_factors, 3L)
})

test_that("latent_factor rejects invalid n", {
  expect_error(latent_factor(n = 0), "positive integer")
  expect_error(latent_factor(n = -1), "positive integer")
  expect_error(latent_factor(n = 1.5), "positive integer")
  expect_error(latent_factor(n = c(1, 2)), "scalar")
})

test_that("latent_factor accepts custom prior", {
  custom_prior <- prior_half_normal(2)
  lf <- latent_factor(prior = custom_prior)

  expect_equal(lf$prior$distribution, "half_normal")
  expect_equal(lf$prior$parameters$sd, 2)
})

test_that("latent_factor rejects invalid prior", {
  expect_error(latent_factor(prior = "invalid"), "simp_prior")
  expect_error(latent_factor(prior = list()), "simp_prior")
})

test_that("latent_factor uses PC prior by default", {
  lf <- latent_factor()

  expect_equal(lf$prior$distribution, "pc")
})

test_that("is_latent works correctly", {
  expect_true(is_latent(latent_factor()))
  expect_false(is_latent("not latent"))
  expect_false(is_latent(list(n_factors = 1)))
})

test_that("latent_factor print method works", {
  lf <- latent_factor(n = 2)

  expect_output(print(lf), "Latent factor")
  expect_output(print(lf), "2")
})
