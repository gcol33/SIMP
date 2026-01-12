# Tests for model fitting

# Helper to create test data
make_test_data <- function(n = 30, n_sites = 5) {
  set.seed(42)
  site <- factor(rep(seq_len(n_sites), length.out = n))
  temp <- rnorm(n)

  # Compositional response
  props <- matrix(runif(n * 3), ncol = 3)
  props <- props / rowSums(props)
  cover <- simp_composition(props, components = c("grass", "shrub", "bare"))

  list(
    cover = cover,
    data = data.frame(temp = temp, site = site)
  )
}

test_that("fit validates model argument", {
  expect_error(fit("not a model"), "simp_model")
})

test_that("fit validates numeric arguments", {
  td <- make_test_data()
  model <- simp(td$cover ~ temp, data = td$data)

  expect_error(fit(model, chains = 0), "positive integer")
  expect_error(fit(model, iter = -1), "positive integer")
  expect_error(fit(model, adapt_delta = 1.5), "between 0 and 1")
})

test_that("cpp_dirichlet_n_params returns correct count", {
  # 2 predictors (intercept + temp), 3 components, plus phi
  expect_equal(cpp_dirichlet_n_params(2, 3), 2 * 3 + 1)
})

test_that("cpp_dirichlet_init returns correct length", {
  init <- cpp_dirichlet_init(2, 3)
  expect_length(init, 7)  # 2*3 + 1
})

test_that("cpp_test_log_dirichlet computes correct value", {
  # Uniform Dirichlet with alpha = 1,1,1 should give constant
  y <- c(1/3, 1/3, 1/3)
  alpha <- c(1, 1, 1)

  log_prob <- cpp_test_log_dirichlet(y, alpha)

  # log(Gamma(3) / (Gamma(1)^3)) = log(2) - 0 = log(2)
  expect_equal(log_prob, log(2), tolerance = 1e-6)
})

test_that("cpp_test_log_dirichlet_softmax works", {
  y <- c(0.3, 0.5, 0.2)
  eta <- c(0, 0, 0)  # Uniform softmax
  phi <- 10

  log_prob <- cpp_test_log_dirichlet_softmax(y, eta, phi)

  # Should be finite
 expect_true(is.finite(log_prob))
})

test_that("fit runs without error on simple model", {
  skip_on_cran()

  td <- make_test_data(n = 20)
  model <- simp(td$cover ~ temp, data = td$data)

  # Quick run with minimal iterations (use zero init for stability)
  fit_obj <- fit(model, chains = 1, iter = 50, warmup = 100,
                 verbose = FALSE, seed = 123, init = "zero")

  expect_s3_class(fit_obj, "simp_fit")
  expect_equal(dim(fit_obj$draws)[1], 50)  # iterations
  expect_equal(dim(fit_obj$draws)[2], 1)   # chains
})

test_that("fit summary works", {
  skip_on_cran()

  td <- make_test_data(n = 20)
  model <- simp(td$cover ~ temp, data = td$data)

  fit_obj <- fit(model, chains = 1, iter = 50, warmup = 100,
                 verbose = FALSE, seed = 123, init = "zero")

  summ <- summary(fit_obj)
  expect_s3_class(summ, "summary.simp_fit")
  expect_true("mean" %in% names(summ$stats))
  expect_true("rhat" %in% names(summ$stats))
})

test_that("as.data.frame.simp_fit works", {
  skip_on_cran()

  td <- make_test_data(n = 20)
  model <- simp(td$cover ~ temp, data = td$data)

  fit_obj <- fit(model, chains = 1, iter = 50, warmup = 100,
                 verbose = FALSE, seed = 123, init = "zero")

  df <- as.data.frame(fit_obj)
  expect_s3_class(df, "data.frame")
  expect_true(".chain" %in% names(df))
  expect_true(".iteration" %in% names(df))
})

test_that("is_simp_fit works", {
  td <- make_test_data(n = 20)
  model <- simp(td$cover ~ temp, data = td$data)

  expect_false(is_simp_fit(model))

  skip_on_cran()
  fit_obj <- fit(model, chains = 1, iter = 50, warmup = 100,
                 verbose = FALSE, seed = 123, init = "zero")
  expect_true(is_simp_fit(fit_obj))
})
