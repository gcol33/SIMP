# Tests for SIMP family functions

test_that("simp_dirichlet creates valid family object", {
  fam <- simp_dirichlet()

  expect_s3_class(fam, "simp_family")
  expect_equal(fam$family, "dirichlet")
  expect_equal(fam$link, "log")
  expect_null(fam$concentration)
  expect_false(fam$handles_zeros)
})

test_that("simp_dirichlet with fixed concentration works", {
  fam <- simp_dirichlet(concentration = 10)

  expect_equal(fam$concentration, 10)
})

test_that("simp_dirichlet rejects invalid concentration", {
  expect_error(simp_dirichlet(concentration = -1), "positive")
  expect_error(simp_dirichlet(concentration = c(1, 2)), "scalar")
})

test_that("simp_dirichlet validates link", {
  expect_error(simp_dirichlet(link = "invalid"), "Invalid link")
})

test_that("simp_logistic_normal creates valid family object", {
  fam <- simp_logistic_normal()

  expect_s3_class(fam, "simp_family")
  expect_equal(fam$family, "logistic_normal")
  expect_equal(fam$reference, "last")
})

test_that("simp_logistic_normal accepts reference options", {
  fam1 <- simp_logistic_normal(reference = "first")
  expect_equal(fam1$reference, "first")

  fam2 <- simp_logistic_normal(reference = 2)
  expect_equal(fam2$reference, 2)
})

test_that("simp_hurdle_dirichlet creates valid family object", {
  fam <- simp_hurdle_dirichlet()

  expect_s3_class(fam, "simp_family")
  expect_equal(fam$family, "hurdle_dirichlet")
  expect_true(fam$handles_zeros)
  expect_equal(fam$link_zero, "logit")
})

test_that("simp_hurdle_dirichlet accepts zero formula", {
  fam <- simp_hurdle_dirichlet(zero_formula = ~ elevation + temp)

  expect_equal(fam$zero_formula, ~ elevation + temp)
})

test_that("simp_hurdle_dirichlet rejects invalid zero_formula", {
  expect_error(simp_hurdle_dirichlet(zero_formula = "not a formula"), "formula")
})

test_that("simp_stick_breaking creates valid family object", {
  fam <- simp_stick_breaking()

  expect_s3_class(fam, "simp_family")
  expect_equal(fam$family, "stick_breaking")
  expect_null(fam$order)
})

test_that("simp_stick_breaking accepts order", {
  fam <- simp_stick_breaking(order = c("tree", "shrub", "herb"))

  expect_equal(fam$order, c("tree", "shrub", "herb"))
})

test_that("is_family works correctly", {
  expect_true(is_family(simp_dirichlet()))
  expect_false(is_family("not a family"))
  expect_false(is_family(list(family = "fake")))
})

test_that("family print method works", {
  fam <- simp_dirichlet()

  output <- capture.output(print(fam))
  output_str <- paste(output, collapse = " ")
  expect_true(grepl("dirichlet", output_str, ignore.case = TRUE))
  expect_true(grepl("log", output_str))
})
