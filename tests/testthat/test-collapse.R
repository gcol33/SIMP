# Tests for collapse functions

test_that("collapse_components aggregates correctly", {
  x <- simp_composition(
    matrix(c(0.1, 0.2, 0.3, 0.2, 0.2,
             0.2, 0.1, 0.2, 0.3, 0.2), nrow = 2, byrow = TRUE),
    components = c("grass1", "grass2", "shrub1", "shrub2", "tree")
  )

  result <- collapse_components(x, mapping = list(
    grass = c("grass1", "grass2"),
    shrub = c("shrub1", "shrub2"),
    tree = "tree"
  ))

  expect_equal(n_components(result), 3)
  expect_equal(components(result), c("grass", "shrub", "tree"))

  # Check values
  expect_equal(unname(as.matrix(result)[1, "grass"]), 0.3, tolerance = 1e-10)
  expect_equal(unname(as.matrix(result)[1, "shrub"]), 0.5, tolerance = 1e-10)
})

test_that("collapse_components errors on invalid mapping", {
  x <- simp_composition(
    matrix(c(0.5, 0.3, 0.2), nrow = 1),
    components = c("a", "b", "c")
  )

  expect_error(
    collapse_components(x, mapping = list(x = "nonexistent")),
    "not found"
  )
})

test_that("collapse_components warns about unmapped components", {
  x <- simp_composition(
    matrix(c(0.4, 0.3, 0.2, 0.1), nrow = 1),
    components = c("a", "b", "c", "d")
  )

  expect_warning(
    collapse_components(x, mapping = list(ab = c("a", "b"))),
    "not in mapping"
  )
})

test_that("collapse_layers works for simp_cover", {
  x <- simp_cover(
    matrix(c(0.1, 0.2, 0.3, 0.4,
             0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
    species = c("oak", "pine", "grass", "bare"),
    layers = c("tree", "tree", "herb", "ground")
  )

  result <- collapse_layers(x)

  expect_s3_class(result, "simp_composition")
  expect_equal(n_components(result), 3)
  expect_equal(components(result), c("tree", "herb", "ground"))

  # Check aggregation
  expect_equal(unname(as.matrix(result)[1, "tree"]), 0.3, tolerance = 1e-10)
})

test_that("split_layer extracts single layer", {
  x <- simp_cover(
    matrix(c(0.1, 0.2, 0.3, 0.4,
             0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
    species = c("oak", "pine", "grass", "bare"),
    layers = c("tree", "tree", "herb", "ground")
  )

  result <- split_layer(x, "tree")

  expect_s3_class(result, "simp_cover")
  expect_equal(n_components(result), 2)
})

test_that("identify_rare finds rare components", {
  x <- simp_composition(
    matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
             0.5, 0.3, 0.2, 0, 0,
             0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
    components = c("sp1", "sp2", "sp3", "sp4", "sp5")
  )

  # sp4 and sp5 have non-zero in 1/3 = 33% of observations
  # So they occur in < 50% and have low mean cover
  rare <- identify_rare(x, threshold = 0.5, min_cover = 0.1)

  expect_true("sp4" %in% rare)
  expect_true("sp5" %in% rare)
})

test_that("collapse_rare aggregates rare components", {
  x <- simp_composition(
    matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
             0.5, 0.3, 0.2, 0, 0,
             0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
    components = c("sp1", "sp2", "sp3", "sp4", "sp5")
  )

  result <- suppressMessages(collapse_rare(x, threshold = 0.5, min_cover = 0.1))

  expect_true("other" %in% components(result))
  expect_false("sp4" %in% components(result))
  expect_false("sp5" %in% components(result))
})
