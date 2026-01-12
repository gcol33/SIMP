# Tests for simp() model specification function

# Helper to create test data
make_test_data <- function(n = 50, n_sites = 10) {
  set.seed(42)
  site <- factor(rep(seq_len(n_sites), length.out = n))
  temp <- rnorm(n)
  precip <- rnorm(n)

  # Compositional response
  props <- matrix(runif(n * 3), ncol = 3)
  props <- props / rowSums(props)
  cover <- simp_composition(props, components = c("grass", "shrub", "bare"))

  list(
    cover = cover,
    data = data.frame(temp = temp, precip = precip, site = site)
  )
}

test_that("simp creates valid model object", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp, data = td$data)

  expect_s3_class(model, "simp_model")
  expect_false(model$fitted)
})

test_that("simp parses fixed effects correctly", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp + precip, data = td$data)

  expect_true("(Intercept)" %in% model$fixed$terms)
  expect_true("temp" %in% model$fixed$terms)
  expect_true("precip" %in% model$fixed$terms)
  expect_equal(nrow(model$fixed$X), 50)
})

test_that("simp parses random intercepts correctly", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp + (1 | site), data = td$data)

  expect_length(model$random, 1)
  expect_equal(model$random$site$group, "site")
  expect_true(model$random$site$has_intercept)
  expect_equal(model$random$site$n_levels, 10)
})

test_that("simp parses random slopes correctly", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp + (1 + temp | site), data = td$data)

  expect_equal(model$random$site$terms, "temp")
  expect_true(model$random$site$has_intercept)
  expect_true(model$random$site$correlated)
})

test_that("simp parses uncorrelated random effects", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp + (1 + temp || site), data = td$data)

  expect_false(model$random$site$correlated)
})

test_that("simp handles nested random effects", {
  td <- make_test_data()
  td$data$plot <- factor(rep(1:5, 10))

  model <- simp(td$cover ~ temp + (1 | site/plot), data = td$data)

  # Should expand to two random effects
  expect_length(model$random, 2)
  expect_true("site" %in% names(model$random))
  expect_true("site:plot" %in% names(model$random))
})

test_that("simp handles crossed random effects", {
  td <- make_test_data()
  td$data$year <- factor(rep(1:5, 10))

  model <- simp(td$cover ~ temp + (1 | site) + (1 | year), data = td$data)

  expect_length(model$random, 2)
  expect_true("site" %in% names(model$random))
  expect_true("year" %in% names(model$random))
})

test_that("simp validates formula argument", {
  td <- make_test_data()

  expect_error(simp("not a formula", data = td$data), "formula")
})
test_that("simp validates data argument", {
  td <- make_test_data()

  expect_error(simp(td$cover ~ temp, data = "not a data.frame"), "data frame")
})

test_that("simp validates family argument", {
  td <- make_test_data()

  expect_error(simp(td$cover ~ temp, data = td$data, family = "invalid"), "simp_family")
})

test_that("simp validates prior argument", {
  td <- make_test_data()

  expect_error(
    simp(td$cover ~ temp, data = td$data, prior = "invalid"),
    "simp_priors"
  )
})

test_that("simp uses default priors from family", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp, data = td$data)

  expect_s3_class(model$prior, "simp_priors")
})

test_that("simp accepts custom priors", {
  td <- make_test_data()

  custom_prior <- simp_priors(beta = prior_normal(0, 10))
  model <- simp(td$cover ~ temp, data = td$data, prior = custom_prior)

  expect_equal(model$prior$beta$parameters$sd, 10)
})

test_that("simp extracts response correctly", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp, data = td$data)

  expect_equal(model$response$n_obs, 50)
  expect_equal(model$response$n_components, 3)
  expect_equal(model$response$component_names, c("grass", "shrub", "bare"))
})

test_that("simp warns about zeros with non-hurdle family", {
  td <- make_test_data()
  # Add zeros
  mat <- as.matrix(td$cover)
  mat[1:5, 1] <- 0
  mat <- mat / rowSums(mat)
  cover_zeros <- simp_composition(mat, components = c("grass", "shrub", "bare"))

  expect_warning(
    simp(cover_zeros ~ temp, data = td$data),
    "zeros"
  )
})

test_that("simp accepts hurdle family for zero data", {
  td <- make_test_data()
  mat <- as.matrix(td$cover)
  mat[1:5, 1] <- 0
  mat <- mat / rowSums(mat)
  cover_zeros <- simp_composition(mat, components = c("grass", "shrub", "bare"))

  # Should not warn with hurdle family
  expect_no_warning(
    simp(cover_zeros ~ temp, data = td$data, family = simp_hurdle_dirichlet())
  )
})

test_that("simping returns TRUE", {
  expect_true(simping())
})

test_that("simp print method works", {
  td <- make_test_data()

  model <- simp(td$cover ~ temp + (1 | site), data = td$data)

  output <- capture.output(print(model))
  output_str <- paste(output, collapse = " ")
  # Check key content is present (CLI formatting may vary)
  expect_true(grepl("temp", output_str))
  expect_true(grepl("site", output_str))
  expect_true(grepl("dirichlet", output_str))
})
