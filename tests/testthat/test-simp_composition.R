# Tests for simp_composition S3 class

test_that("simp_composition creates valid objects from matrix", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0.4, 0.4, 0.2,
                0.5, 0.3, 0.2), nrow = 3, byrow = TRUE)

  comp <- simp_composition(x, components = c("a", "b", "c"))

  expect_s3_class(comp, "simp_composition")
  expect_equal(n_obs(comp), 3)
  expect_equal(n_components(comp), 3)
  expect_equal(components(comp), c("a", "b", "c"))
})

test_that("simp_composition creates valid objects from data.frame", {
  df <- data.frame(sp1 = c(0.3, 0.4), sp2 = c(0.5, 0.4), sp3 = c(0.2, 0.2))

  comp <- simp_composition(df)

  expect_s3_class(comp, "simp_composition")
  expect_equal(components(comp), c("sp1", "sp2", "sp3"))
})

test_that("simp_composition rejects invalid input", {
  # Rows don't sum to 1
  x <- matrix(c(0.3, 0.5, 0.3,
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  expect_error(simp_composition(x), "Invalid composition")
})

test_that("simp_composition with normalize = TRUE fixes row sums", {
  x <- matrix(c(30, 50, 20,
                40, 40, 20), nrow = 2, byrow = TRUE)

  comp <- simp_composition(x, normalize = TRUE)

  expect_true(all(abs(rowSums(as.matrix(comp)) - 1) < 1e-10))
})

test_that("simp_composition subsetting preserves class", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0.4, 0.4, 0.2,
                0.5, 0.3, 0.2), nrow = 3, byrow = TRUE)

  comp <- simp_composition(x, components = c("a", "b", "c"))

  # Row subset
  sub1 <- comp[1:2, ]
  expect_s3_class(sub1, "simp_composition")
  expect_equal(n_obs(sub1), 2)

  # Column subset
  sub2 <- comp[, 1:2]
  expect_s3_class(sub2, "simp_composition")
  expect_equal(n_components(sub2), 2)
})

test_that("simp_composition rbind works", {
  x1 <- simp_composition(
    matrix(c(0.3, 0.5, 0.2), nrow = 1),
    components = c("a", "b", "c")
  )
  x2 <- simp_composition(
    matrix(c(0.4, 0.4, 0.2), nrow = 1),
    components = c("a", "b", "c")
  )

  combined <- rbind(x1, x2)

  expect_s3_class(combined, "simp_composition")
  expect_equal(n_obs(combined), 2)
})

test_that("simp_composition detects zeros", {
  x <- matrix(c(0.5, 0.5, 0,
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  comp <- simp_composition(x, components = c("a", "b", "c"))

  expect_true(attr(comp, "has_zeros"))
})

test_that("as.matrix and as.data.frame work", {
  x <- matrix(c(0.3, 0.5, 0.2,
                0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  comp <- simp_composition(x, components = c("a", "b", "c"))

  mat <- as.matrix(comp)
  expect_true(is.matrix(mat))
  expect_false(inherits(mat, "simp_composition"))

  df <- as.data.frame(comp)
  expect_s3_class(df, "data.frame")
  expect_equal(names(df), c("a", "b", "c"))
})
