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
