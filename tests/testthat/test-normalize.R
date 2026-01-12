# Tests for normalization functions

test_that("normalize_cover rescales rows to sum to 1", {
  x <- matrix(c(30, 50, 20,
                40, 40, 20), nrow = 2, byrow = TRUE)

  result <- normalize_cover(x)

  expect_true(all(abs(rowSums(result) - 1) < 1e-10))
})
test_that("normalize_cover preserves valid compositions", {
  x <- simp_composition(
    matrix(c(0.3, 0.5, 0.2,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  result <- normalize_cover(x)

  expect_equal(as.matrix(result), as.matrix(x), tolerance = 1e-10)
})

test_that("normalize_cover errors on negative values", {
  x <- matrix(c(-0.1, 0.6, 0.5), nrow = 1)

  expect_error(normalize_cover(x), "negative")
})

test_that("normalize_cover with truncate method works", {
  x <- matrix(c(1.5, 0.5, 0.2), nrow = 1)

  result <- normalize_cover(x, method = "truncate")

  expect_true(all(result <= 1))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("add_residual closes open composition", {
  x <- matrix(c(0.3, 0.3, 0.2,
                0.4, 0.2, 0.2), nrow = 2, byrow = TRUE)
  colnames(x) <- c("sp1", "sp2", "sp3")

  result <- add_residual(x, name = "other")

  expect_equal(ncol(result), 4)
  expect_true("other" %in% colnames(result))
  expect_true(all(abs(rowSums(result) - 1) < 1e-10))
})

test_that("add_residual errors when sum > 1", {
  x <- matrix(c(0.4, 0.4, 0.4), nrow = 1)

  expect_error(add_residual(x), "exceed")
})

test_that("add_residual position argument works", {
  x <- matrix(c(0.3, 0.5), nrow = 1)
  colnames(x) <- c("a", "b")

  result_last <- add_residual(x, name = "other", position = "last")
  result_first <- add_residual(x, name = "other", position = "first")

  expect_equal(colnames(result_last), c("a", "b", "other"))
  expect_equal(colnames(result_first), c("other", "a", "b"))
})

test_that("remove_residual works for simp_cover", {
  x <- simp_cover(
    matrix(c(0.3, 0.5, 0.2,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    species = c("grass", "shrub", "bare"),
    residual = "bare"
  )

  result <- remove_residual(x)

  expect_equal(n_components(result), 2)
  expect_false("bare" %in% components(result))
})

test_that("transform_composition CLR works", {
  x <- simp_composition(
    matrix(c(0.3, 0.5, 0.2,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  result <- transform_composition(x, to = "clr")

  # CLR should sum to 0 for each row
  expect_true(all(abs(rowSums(result)) < 1e-10))
})

test_that("transform_composition ALR works", {
  x <- simp_composition(
    matrix(c(0.3, 0.5, 0.2,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  result <- transform_composition(x, to = "alr", reference = "c")

  # ALR should have K-1 columns
  expect_equal(ncol(result), 2)
})

test_that("transform_composition ILR works", {
  x <- simp_composition(
    matrix(c(0.3, 0.5, 0.2,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  result <- transform_composition(x, to = "ilr")

  # ILR should have K-1 columns
  expect_equal(ncol(result), 2)
})

test_that("transform_composition errors on zeros", {
  x <- simp_composition(
    matrix(c(0.5, 0.5, 0,
             0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
    components = c("a", "b", "c")
  )

  expect_error(transform_composition(x, to = "clr"), "non-zero")
})
