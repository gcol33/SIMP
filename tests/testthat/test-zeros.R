# Tests for zero handling functions

test_that("identify_zeros counts zeros correctly", {
  x <- simp_composition(
    matrix(c(0.4, 0.4, 0.2, 0,
             0.5, 0.3, 0.2, 0,
             0.3, 0.5, 0, 0.2), nrow = 3, byrow = TRUE),
    components = c("sp1", "sp2", "sp3", "sp4")
  )

  result <- identify_zeros(x)

  expect_equal(result$n_zeros, 3)
})

test_that("identify_zeros classifies structural zeros", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0.6, 0,
             0.3, 0.7, 0,
             0.6, 0.4, 0,
             0.5, 0.5, 0), nrow = 5, byrow = TRUE),
    components = c("sp1", "sp2", "sp3")
  )

  result <- identify_zeros(x, structural_threshold = 0.8)

  expect_true("sp3" %in% result$structural_components)
})

test_that("impute_zeros multiplicative method works", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("sp1", "sp2", "sp3")
  )

  result <- impute_zeros(x, method = "multiplicative")

  # Should have no zeros
  expect_true(all(as.matrix(result) > 0))

  # Should still sum to 1
  expect_true(all(abs(rowSums(as.matrix(result)) - 1) < 1e-10))
})

test_that("impute_zeros additive method works", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("sp1", "sp2", "sp3")
  )

  result <- impute_zeros(x, method = "additive")

  # Should have no zeros
  expect_true(all(as.matrix(result) > 0))

  # Should still sum to 1
  expect_true(all(abs(rowSums(as.matrix(result)) - 1) < 1e-10))
})

test_that("impute_zeros with custom delta works", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0), nrow = 1),
    components = c("sp1", "sp2", "sp3")
  )

  result <- impute_zeros(x, method = "additive", delta = 0.001)

  # Check that zero was replaced with something small
  expect_true(as.matrix(result)[1, 3] < 0.01)
})

test_that("impute_zeros preserve_structural works", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0.6, 0,
             0.3, 0.7, 0,
             0.6, 0.4, 0,
             0.5, 0.4, 0.1), nrow = 5, byrow = TRUE),
    components = c("sp1", "sp2", "sp3")
  )

  result <- impute_zeros(x, preserve_structural = TRUE,
                         structural_threshold = 0.8)

  # sp3 should still have zeros (structural) except row 5
  mat <- as.matrix(result)
  expect_equal(mat[1:4, 3], rep(0, 4))
  expect_true(mat[5, 3] > 0)
})

test_that("has_zeros returns correct value", {
  x_with <- simp_composition(
    matrix(c(0.5, 0.5, 0), nrow = 1),
    components = c("a", "b", "c")
  )

  x_without <- simp_composition(
    matrix(c(0.3, 0.5, 0.2), nrow = 1),
    components = c("a", "b", "c")
  )

  expect_true(has_zeros(x_with))
  expect_false(has_zeros(x_without))
})

test_that("count_zeros returns correct count", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0, 0.6), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  expect_equal(count_zeros(x), 2)
})

test_that("impute_zeros returns original if no zeros", {
  x <- simp_composition(
    matrix(c(0.3, 0.5, 0.2), nrow = 1),
    components = c("a", "b", "c")
  )

  expect_message(result <- impute_zeros(x), "No zeros")
  expect_identical(result, x)
})
