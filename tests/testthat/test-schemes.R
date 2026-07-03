test_that("registry lists leapfrog, the 2-stage scheme, and the Yoshida family", {
  s <- simp_schemes()
  expect_true(all(c("leapfrog", "minerror2", "yoshida4", "yoshida6", "yoshida8")
                  %in% s$integrator))

  ord <- stats::setNames(s$order, s$integrator)
  expect_equal(unname(ord["leapfrog"]), 2L)
  expect_equal(unname(ord["minerror2"]), 2L)
  expect_equal(unname(ord["yoshida4"]), 4L)
  expect_equal(unname(ord["yoshida6"]), 6L)
  expect_equal(unname(ord["yoshida8"]), 8L)

  ge <- stats::setNames(s$gradient_evals, s$integrator)
  # Triple-jump composition triples the stage count at each order step; the
  # 2-stage scheme costs two force evaluations per step.
  expect_equal(unname(ge["leapfrog"]), 1L)
  expect_equal(unname(ge["minerror2"]), 2L)
  expect_equal(unname(ge["yoshida4"]), 3L)
  expect_equal(unname(ge["yoshida6"]), 9L)
  expect_equal(unname(ge["yoshida8"]), 27L)
})

test_that("an unknown integrator name errors", {
  grad_U <- function(q) q
  expect_error(
    simp_integrate(1, 0, grad_U, eps = 0.1, n_steps = 1, integrator = "nope"),
    "unknown integrator"
  )
})

test_that("built-in and adaptive schemes are valid palindromic splittings", {
  # kicks sum to one, drifts sum to one, and the op sequence is a palindrome, so
  # the step is consistent and time reversible.
  check_valid <- function(ops) {
    kind <- ops$kind
    coef <- ops$coef
    expect_equal(sum(coef[kind == 1]), 1, tolerance = 1e-12)  # kicks
    expect_equal(sum(coef[kind == 0]), 1, tolerance = 1e-12)  # drifts
    expect_equal(kind, rev(kind))
    expect_equal(coef, rev(coef), tolerance = 1e-12)
  }
  for (nm in simp_schemes()$integrator) {
    check_valid(SIMP:::simp_scheme_ops_cpp(nm))
  }
  check_valid(SIMP:::simp_adaptive_ops_cpp(1.3, 2L))
  check_valid(SIMP:::simp_adaptive_ops_cpp(1.3, 3L))
})
