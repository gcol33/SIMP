# Tests for formula parsing

test_that("parse_formula extracts fixed effects", {
  f <- parse_formula(y ~ x1 + x2)

  expect_equal(f$fixed$terms, c("x1", "x2"))
  expect_false(f$has_random)
})

test_that("parse_formula extracts random intercepts", {
  f <- parse_formula(y ~ x1 + (1 | group))

  expect_true(f$has_random)
  expect_length(f$random, 1)
  expect_equal(f$random[[1]]$group, "group")
  expect_true(f$random[[1]]$has_intercept)
  expect_length(f$random[[1]]$terms, 0)
})

test_that("parse_formula extracts random slopes", {
  f <- parse_formula(y ~ x1 + (1 + x1 | group))

  expect_true(f$random[[1]]$has_intercept)
  expect_equal(f$random[[1]]$terms, "x1")
  expect_true(f$random[[1]]$correlated)
})

test_that("parse_formula handles uncorrelated syntax", {
  f <- parse_formula(y ~ x1 + (1 + x1 || group))

  expect_false(f$random[[1]]$correlated)
})

test_that("parse_formula handles slope-only (0 + x)", {
  f <- parse_formula(y ~ x1 + (0 + x1 | group))

  expect_false(f$random[[1]]$has_intercept)
  expect_equal(f$random[[1]]$terms, "x1")
})

test_that("parse_formula expands nested effects", {
  f <- parse_formula(y ~ x1 + (1 | site/plot))

  # Should expand to (1|site) + (1|site:plot)
  expect_length(f$random, 2)
  expect_equal(f$random[[1]]$group, "site")
  expect_equal(f$random[[2]]$group, "site:plot")
})

test_that("parse_formula handles crossed effects", {
  f <- parse_formula(y ~ x1 + (1 | site) + (1 | year))

  expect_length(f$random, 2)
  expect_equal(f$random[[1]]$group, "site")
  expect_equal(f$random[[2]]$group, "year")
})

test_that("extract_response handles univariate response", {
  f <- parse_formula(y ~ x1)

  expect_equal(f$response$type, "univariate")
  expect_equal(f$response$variables, "y")
})

test_that("extract_response handles cbind response", {
  f <- parse_formula(cbind(sp1, sp2, sp3) ~ x1)

  expect_equal(f$response$type, "multivariate")
  expect_equal(f$response$variables, c("sp1", "sp2", "sp3"))
  expect_equal(f$response$n_components, 3)
})

test_that("remove_random_effects works", {
  f_orig <- y ~ x1 + x2 + (1 | group) + (1 + x1 | site)
  f_fixed <- remove_random_effects(f_orig)

  expect_equal(deparse(f_fixed), "y ~ x1 + x2")
})

test_that("validate_formula validates against data", {
  data <- data.frame(y = 1:5, x1 = rnorm(5), group = rep(1:2, c(3, 2)))

  expect_true(validate_formula(y ~ x1 + (1 | group), data))
})

test_that("validate_formula errors on missing variables", {
  data <- data.frame(y = 1:5, x1 = rnorm(5))

  expect_error(
    validate_formula(y ~ x1 + x2, data),
    "not found"
  )
})

test_that("validate_formula errors on missing grouping factor", {
  data <- data.frame(y = 1:5, x1 = rnorm(5))

  expect_error(
    validate_formula(y ~ x1 + (1 | group), data),
    "not found"
  )
})
