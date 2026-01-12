# Prior specification for SIMP models
# Ported from ratiod/R/priors.R with adaptations for compositional data

# -----------------------------------------------------------------------------
# Prior constructors
# -----------------------------------------------------------------------------

#' Create a normal prior
#'
#' @param mean Prior mean (default 0)
#' @param sd Prior standard deviation (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_normal(0, 2.5)
#' prior_normal(mean = 0, sd = 10)
prior_normal <- function(mean = 0, sd = 1) {
  check_scalar(mean, "mean")
  check_scalar(sd, "sd")

  check_positive(sd, "sd")

  new_simp_prior(
    distribution = "normal",
    parameters = list(mean = mean, sd = sd)
  )
}

#' Create a half-normal prior
#'
#' Half-normal distribution for positive parameters. The distribution is a
#' normal distribution truncated at 0.
#'
#' @param sd Prior scale parameter (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_half_normal(2.5)
prior_half_normal <- function(sd = 1) {
  check_scalar(sd, "sd")
  check_positive(sd, "sd")

  new_simp_prior(
    distribution = "half_normal",
    parameters = list(sd = sd)
  )
}

#' Create a half-Cauchy prior
#'
#' Half-Cauchy distribution for positive parameters. Recommended for variance
#' components due to heavy tails that allow for large values while maintaining
#' good behavior near zero.
#'
#' @param location Location parameter (default 0)
#' @param scale Scale parameter (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_half_cauchy(0, 2)
prior_half_cauchy <- function(location = 0, scale = 1) {
  check_scalar(location, "location")
  check_scalar(scale, "scale")
  check_non_negative(location, "location")
  check_positive(scale, "scale")

  new_simp_prior(
    distribution = "half_cauchy",
    parameters = list(location = location, scale = scale)
  )
}

#' Create a gamma prior
#'
#' Gamma distribution for positive parameters.
#'
#' @param shape Shape parameter (default 1)
#' @param rate Rate parameter (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_gamma(2, 0.1)
prior_gamma <- function(shape = 1, rate = 1) {
  check_scalar(shape, "shape")
  check_scalar(rate, "rate")
  check_positive(shape, "shape")
  check_positive(rate, "rate")

  new_simp_prior(
    distribution = "gamma",
    parameters = list(shape = shape, rate = rate)
  )
}

#' Create an exponential prior
#'
#' Exponential distribution for positive parameters.
#'
#' @param rate Rate parameter (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_exponential(0.5)
prior_exponential <- function(rate = 1) {
  check_scalar(rate, "rate")
  check_positive(rate, "rate")

  new_simp_prior(
    distribution = "exponential",
    parameters = list(rate = rate)
  )
}

#' Create a beta prior
#'
#' Beta distribution for parameters bounded in \[0, 1\].
#'
#' @param alpha First shape parameter (default 1)
#' @param beta Second shape parameter (default 1)
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_beta(2, 2)  # Symmetric, peaked at 0.5
#' prior_beta(1, 1)  # Uniform
prior_beta <- function(alpha = 1, beta = 1) {
  check_scalar(alpha, "alpha")
  check_scalar(beta, "beta")
  check_positive(alpha, "alpha")
  check_positive(beta, "beta")

  new_simp_prior(
    distribution = "beta",
    parameters = list(alpha = alpha, beta = beta)
  )
}

#' Create a penalized complexity (PC) prior
#'
#' Penalized complexity priors shrink toward a simple base model. The prior is
#' specified through P(sigma > u) = alpha, where u is a threshold and alpha is
#' the tail probability.
#'
#' @param u Upper bound threshold
#' @param alpha Tail probability P(sigma > u)
#' @return A `simp_prior` object
#' @references Simpson et al. (2017) Penalising Model Component Complexity:
#'   A Principled, Practical Approach to Constructing Priors.
#'   Statistical Science, 32(1), 1-28.
#' @export
#' @examples
#' # P(sigma > 1) = 0.01
#' prior_pc(u = 1, alpha = 0.01)
prior_pc <- function(u, alpha) {
  if (missing(u) || missing(alpha)) {
    cli::cli_abort(c(
      "Both {.arg u} and {.arg alpha} are required for PC priors",
      "i" = "Use {.code prior_pc(u = 1, alpha = 0.01)} to set P(sigma > 1) = 0.01"
    ))
  }
  check_scalar(u, "u")
  check_scalar(alpha, "alpha")
  check_positive(u, "u")
  check_probability(alpha, "alpha")

  # Compute rate parameter: lambda = -log(alpha) / u
  lambda <- -log(alpha) / u

  new_simp_prior(
    distribution = "pc",
    parameters = list(u = u, alpha = alpha, lambda = lambda)
  )
}

#' Create a LKJ prior for correlation matrices
#'
#' LKJ distribution for correlation matrices. The parameter eta controls the
#' concentration around the identity matrix.
#'
#' @param eta Concentration parameter. eta = 1 gives uniform distribution over
#'   correlation matrices. eta > 1 shrinks toward identity matrix.
#' @return A `simp_prior` object
#' @export
#' @examples
#' prior_lkj(2)  # Shrink toward uncorrelated
prior_lkj <- function(eta = 1) {
  check_scalar(eta, "eta")
  check_positive(eta, "eta")

  new_simp_prior(
    distribution = "lkj",
    parameters = list(eta = eta)
  )
}

# -----------------------------------------------------------------------------
# Prior container
# -----------------------------------------------------------------------------

#' Specify priors for a SIMP model
#'
#' Create a container of priors for different model parameters.
#'
#' @param beta Prior for fixed effect coefficients (default: normal(0, 2.5))
#' @param intercept Prior for intercept (default: normal(0, 5))
#' @param sigma Prior for random effect standard deviations
#'   (default: half_cauchy(0, 2))
#' @param phi Prior for Dirichlet concentration parameter
#'   (default: gamma(2, 0.1))
#' @param cor Prior for correlation matrices in random effects
#'   (default: lkj(2))
#' @param ... Additional named priors for custom parameters
#' @return A `simp_priors` object
#' @export
#' @examples
#' # Use defaults
#' simp_priors()
#'
#' # Custom priors
#' simp_priors(
#'   beta = prior_normal(0, 10),
#'   sigma = prior_half_normal(2)
#' )
simp_priors <- function(beta = prior_normal(0, 2.5),
                        intercept = prior_normal(0, 5),
                        sigma = prior_half_cauchy(0, 2),
                        phi = prior_gamma(2, 0.1),
                        cor = prior_lkj(2),
                        ...) {

  # Validate that all arguments are simp_prior objects
  priors <- list(
    beta = beta,
    intercept = intercept,
    sigma = sigma,
    phi = phi,
    cor = cor,
    ...
  )

  for (nm in names(priors)) {
    if (!inherits(priors[[nm]], "simp_prior")) {
      cli::cli_abort(c(
        "{.arg {nm}} must be a {.cls simp_prior} object",
        "i" = "Use {.fn prior_normal}, {.fn prior_half_cauchy}, etc."
      ))
    }
  }

  structure(
    priors,
    class = "simp_priors"
  )
}

# -----------------------------------------------------------------------------
# S3 class constructor and methods
# -----------------------------------------------------------------------------

#' Create a new simp_prior object
#' @noRd
new_simp_prior <- function(distribution, parameters) {
  structure(
    list(
      distribution = distribution,
      parameters = parameters
    ),
    class = "simp_prior"
  )
}

#' @export
print.simp_prior <- function(x, ...) {
  params <- x$parameters
  param_str <- paste(names(params), "=", format_num(unlist(params)), collapse = ", ")
  cat(x$distribution, "(", param_str, ")\n", sep = "")
  invisible(x)
}

#' @export
print.simp_priors <- function(x, ...) {
  cli::cli_h3("SIMP Priors")
  for (nm in names(x)) {
    cat("  ", nm, ": ", sep = "")
    print(x[[nm]])
  }
  invisible(x)
}

#' Get default priors for a model family
#'
#' @param family A SIMP model family (character or family object)
#' @return A `simp_priors` object with appropriate defaults
#' @export
#' @examples
#' default_priors("dirichlet")
default_priors <- function(family = "dirichlet") {
  family_name <- if (is.character(family)) family else family$family

  switch(family_name,
    dirichlet = simp_priors(
      beta = prior_normal(0, 2.5),
      intercept = prior_normal(0, 5),
      sigma = prior_half_cauchy(0, 2),
      phi = prior_gamma(2, 0.1),
      cor = prior_lkj(2)
    ),
    logistic_normal = simp_priors(
      beta = prior_normal(0, 2.5),
      intercept = prior_normal(0, 5),
      sigma = prior_half_cauchy(0, 2),
      Sigma = prior_lkj(2),
      cor = prior_lkj(2)
    ),
    hurdle_dirichlet = simp_priors(
      beta = prior_normal(0, 2.5),
      intercept = prior_normal(0, 5),
      sigma = prior_half_cauchy(0, 2),
      phi = prior_gamma(2, 0.1),
      gamma = prior_normal(0, 2.5),  # hurdle coefficients
      cor = prior_lkj(2)
    ),
    cli::cli_abort("Unknown family: {.val {family_name}}")
  )
}

#' Check if an object is a simp_prior
#'
#' @param x Object to check
#' @return Logical
#' @export
is_prior <- function(x) {
  inherits(x, "simp_prior")
}

#' Check if an object is a simp_priors container
#'
#' @param x Object to check
#' @return Logical
#' @export
is_priors <- function(x) {
  inherits(x, "simp_priors")
}
