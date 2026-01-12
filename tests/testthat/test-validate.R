# Tests for validation functions

test_that("validate_composition detects row sum issues", {
  x <- matrix(c(0.3, 0.5, 0.3,  # sums to 1.1
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  result <- validate_composition(x, verbose = FALSE)

  expect_false(result$valid)
  expect_equal(result$row_sums$n_issues, 1)
})

test_that("validate_composition detects negative values", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0.4, -0.1, 0.7), nrow = 2, byrow = TRUE)

  result <- validate_composition(x, verbose = FALSE)

  expect_false(result$valid)
  expect_true(result$issues$has_negatives)
})

test_that("validate_composition detects values > 1", {
  x <- matrix(c(1.2, 0.5, -0.7,
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  result <- validate_composition(x, verbose = FALSE)

  expect_false(result$valid)
  expect_true(result$issues$has_above_one)
})

test_that("validate_composition counts zeros correctly", {
  x <- matrix(c(0.5, 0.5, 0,
                0, 0.6, 0.4,
                0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)

  result <- validate_composition(x, verbose = FALSE)

  expect_equal(result$zeros$total, 2)
})

test_that("validate_composition accepts valid composition", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  result <- validate_composition(x, verbose = FALSE)

  expect_true(result$valid)
})

test_that("is_composition returns logical", {
  x_valid <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  x_invalid <- matrix(c(0.3, 0.5, 0.3), nrow = 1)

  expect_true(is_composition(x_valid))
  expect_false(is_composition(x_invalid))
})

test_that("validate_simp_data checks minimum observations", {
  x <- matrix(c(0.3, 0.5, 0.2), nrow = 1)

  expect_error(
    validate_simp_data(x, min_obs = 10),
    "Insufficient observations"
  )
})

test_that("validate_simp_data checks minimum components", {
  x <- matrix(rep(1, 20), nrow = 20, ncol = 1)

  expect_error(
    validate_simp_data(x, min_obs = 5, min_components = 2),
    "Insufficient components"
  )
})

test_that("validate_simp_data checks zero rows", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0, 0, 0,
                0.4, 0.4, 0.2,
                0.3, 0.4, 0.3,
                0.2, 0.5, 0.3,
                0.1, 0.6, 0.3,
                0.4, 0.3, 0.3,
                0.5, 0.2, 0.3,
                0.6, 0.1, 0.3,
                0.3, 0.3, 0.4,
                0.2, 0.4, 0.4), nrow = 11, byrow = TRUE)

  expect_error(
    validate_simp_data(x),
    "zero"
  )
})

test_that("validate_simp_data checks zero columns", {
  x <- matrix(c(0.5, 0.5, 0,
                0.4, 0.6, 0,
                0.3, 0.7, 0,
                0.6, 0.4, 0,
                0.5, 0.5, 0,
                0.4, 0.6, 0,
                0.3, 0.7, 0,
                0.6, 0.4, 0,
                0.5, 0.5, 0,
                0.4, 0.6, 0,
                0.3, 0.7, 0), nrow = 11, byrow = TRUE)

  expect_error(
    validate_simp_data(x),
    "zeros"
  )
})
