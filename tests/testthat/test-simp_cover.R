# Tests for simp_cover S3 class

test_that("simp_cover creates valid objects", {
  x <- matrix(c(0.3, 0.2, 0.1, 0.4,
                0.4, 0.3, 0.2, 0.1), nrow = 2, byrow = TRUE)

  cover <- simp_cover(x, species = c("oak", "pine", "grass", "bare"))

  expect_s3_class(cover, "simp_cover")
  expect_s3_class(cover, "simp_composition")
  expect_equal(n_species(cover), 4)
  expect_equal(species(cover), c("oak", "pine", "grass", "bare"))
})

test_that("simp_cover handles layers", {
  x <- matrix(c(0.3, 0.2, 0.1, 0.4), nrow = 1)

  cover <- simp_cover(x,
                      species = c("oak", "pine", "grass", "bare"),
                      layers = c("tree", "tree", "herb", "ground"))

  expect_equal(layers(cover), c("tree", "tree", "herb", "ground"))
  expect_equal(n_layers(cover), 3)
})

test_that("simp_cover handles residual", {
  x <- matrix(c(0.3, 0.5, 0.2), nrow = 1)

  cover <- simp_cover(x,
                      species = c("grass", "shrub", "bare"),
                      residual = "bare")

  expect_equal(attr(cover, "residual"), "bare")
})

test_that("simp_cover inherits from simp_composition", {
  x <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  cover <- simp_cover(x, species = c("a", "b", "c"))

  # Should work with simp_composition functions
  expect_equal(n_components(cover), 3)
  expect_equal(components(cover), c("a", "b", "c"))
})

test_that("simp_cover subsetting preserves class", {
  x <- matrix(c(0.3, 0.2, 0.1, 0.4,
                0.4, 0.3, 0.2, 0.1), nrow = 2, byrow = TRUE)

  cover <- simp_cover(x,
                      species = c("oak", "pine", "grass", "bare"),
                      layers = c("tree", "tree", "herb", "ground"))

  # Row subset
  sub1 <- cover[1, ]
  expect_s3_class(sub1, "simp_cover")

  # Column subset
  sub2 <- cover[, 1:2]
  expect_s3_class(sub2, "simp_cover")
  expect_equal(n_species(sub2), 2)
})

test_that("simp_cover normalize works", {
  x <- matrix(c(30, 50, 20), nrow = 1)

  cover <- simp_cover(x, species = c("a", "b", "c"), normalize = TRUE)

  expect_equal(sum(as.matrix(cover)), 1, tolerance = 1e-10)
})

test_that("is_simp_cover returns correct value", {
  x <- matrix(c(0.3, 0.5, 0.2), nrow = 1)

  cover <- simp_cover(x, species = c("a", "b", "c"))
  comp <- simp_composition(x, components = c("a", "b", "c"))

  expect_true(is_simp_cover(cover))
  expect_false(is_simp_cover(comp))
})

test_that("species function extracts names", {
  x <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  cover <- simp_cover(x, species = c("grass", "shrub", "tree"))

  expect_equal(species(cover), c("grass", "shrub", "tree"))
})
