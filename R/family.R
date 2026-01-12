# Family specification for SIMP models
# Defines likelihood families for compositional data

# -----------------------------------------------------------------------------
# Family constructors
# -----------------------------------------------------------------------------

#' Dirichlet family for compositional data
#'
#' The Dirichlet distribution is the natural choice for compositional data,
#' producing vectors that sum to 1 with all positive components.
#'
#' @param link Link function for the concentration parameter (default "log")
#' @param concentration Fixed concentration parameter, or NULL to estimate
#' @return A `simp_family` object
#' @export
#' @examples
#' # Default: estimate concentration
#' simp_dirichlet()
#'
#' # Fixed concentration
#' simp_dirichlet(concentration = 10)
simp_dirichlet <- function(link = "log", concentration = NULL) {
  if (!is.null(concentration)) {
    check_scalar(concentration, "concentration")
    check_positive(concentration, "concentration")
  }

  validate_link(link, valid = c("log", "identity"))

  new_simp_family(
    family = "dirichlet",
    link = link,
    n_params = 1L,
    param_names = "phi",
    concentration = concentration,
    handles_zeros = FALSE,
    description = "Dirichlet distribution for compositional data"
  )
}

#' Logistic-normal family for compositional data
#'
#' The logistic-normal distribution models compositional data via a
#' multivariate normal distribution on the log-ratio transformed space.
#' This allows flexible correlation structures between components.
#'
#' @param link Link function for the latent space (default "identity")
#' @param reference Reference component for the additive log-ratio transform.
#'   Either "last" (default), "first", or a component name/index.
#' @return A `simp_family` object
#' @export
#' @examples
#' simp_logistic_normal()
#' simp_logistic_normal(reference = "first")
simp_logistic_normal <- function(link = "identity", reference = "last") {
  validate_link(link, valid = c("identity", "log"))

  if (!is.character(reference) && !is.numeric(reference))

    cli::cli_abort("{.arg reference} must be character or numeric")

  new_simp_family(
    family = "logistic_normal",
    link = link,
    n_params = 2L,
    param_names = c("mu", "Sigma"),
    reference = reference,
    handles_zeros = FALSE,
    description = "Logistic-normal distribution for compositional data"
  )
}

#' Hurdle Dirichlet family for zero-inflated compositional data
#'
#' A two-part model for compositional data with structural zeros. The first
#' part models the probability of each component being zero; the second part
#' models the non-zero composition using a Dirichlet distribution.
#'
#' @param zero_formula Formula for the zero model (default `~ 1`)
#' @param link_zero Link function for zero probability (default "logit")
#' @param link_positive Link function for positive concentration (default "log")
#' @return A `simp_family` object
#' @export
#' @examples
#' # Intercept-only zero model
#' simp_hurdle_dirichlet()
#'
#' # Covariate-dependent zeros
#' simp_hurdle_dirichlet(zero_formula = ~ elevation)
simp_hurdle_dirichlet <- function(zero_formula = ~ 1,
                                   link_zero = "logit",
                                   link_positive = "log") {
  if (!inherits(zero_formula, "formula")) {
    cli::cli_abort("{.arg zero_formula} must be a formula")
  }

  validate_link(link_zero, valid = c("logit", "probit", "cloglog"))
  validate_link(link_positive, valid = c("log", "identity"))

  new_simp_family(
    family = "hurdle_dirichlet",
    link = link_positive,
    link_zero = link_zero,
    zero_formula = zero_formula,
    n_params = 2L,
    param_names = c("phi", "gamma"),
    handles_zeros = TRUE,
    description = "Hurdle Dirichlet for zero-inflated compositional data"
  )
}

#' Stick-breaking family for ordered compositional data
#'
#' The stick-breaking construction models compositional data as a sequence
#' of binary splits. Useful when components have a natural ordering.
#'
#' @param order Component ordering. Either NULL (estimate or use data order),
#'   or a character vector of component names specifying the order.
#' @param link Link function for the break probabilities (default "logit")
#' @return A `simp_family` object
#' @export
#' @examples
#' simp_stick_breaking()
#'
#' # Specified order
#' simp_stick_breaking(order = c("tree", "shrub", "herb", "bare"))
simp_stick_breaking <- function(order = NULL, link = "logit") {
  if (!is.null(order)) {
    check_character(order, "order")
  }

  validate_link(link, valid = c("logit", "probit"))

  new_simp_family(
    family = "stick_breaking",
    link = link,
    n_params = 1L,
    param_names = "theta",
    order = order,
    handles_zeros = FALSE,
    description = "Stick-breaking construction for ordered compositions"
  )
}

# -----------------------------------------------------------------------------
# Family constructor and methods
# -----------------------------------------------------------------------------

#' Low-level constructor for simp_family
#' @noRd
new_simp_family <- function(family, link, n_params, param_names, ...) {
  structure(
    list(
      family = family,
      link = link,
      n_params = n_params,
      param_names = param_names,
      ...
    ),
    class = "simp_family"
  )
}

#' @export
print.simp_family <- function(x, ...) {
  cli::cli_h3("SIMP Family: {x$family}")
  cat("Link:", x$link, "\n")

  if (!is.null(x$description)) {
    cat("Description:", x$description, "\n")
  }

  if (x$family == "hurdle_dirichlet") {
    cat("Zero link:", x$link_zero, "\n")
    cat("Zero formula:", deparse(x$zero_formula), "\n")
  }

  if (x$family == "logistic_normal") {
    cat("Reference:", x$reference, "\n")
  }

  if (x$family == "stick_breaking" && !is.null(x$order)) {
    cat("Order:", paste(x$order, collapse = " > "), "\n")
  }

  if (x$family == "dirichlet" && !is.null(x$concentration)) {
    cat("Concentration (fixed):", x$concentration, "\n")
  }

  cat("Handles zeros:", x$handles_zeros, "\n")

  invisible(x)
}

#' Check if object is a simp_family
#'
#' @param x Object to check
#' @return Logical
#' @export
is_family <- function(x) {
  inherits(x, "simp_family")
}

# -----------------------------------------------------------------------------
# Link function validation
# -----------------------------------------------------------------------------

#' Validate link function name
#' @noRd
validate_link <- function(link, valid) {
  if (!link %in% valid) {
    cli::cli_abort(c(
      "Invalid link function: {.val {link}}",
      "i" = "Valid options: {.val {valid}}"
    ))
  }
  invisible(link)
}

#' Get link function
#' @noRd
get_link_function <- function(link) {
  switch(link,
    log = list(
      linkfun = log,
      linkinv = exp,
      name = "log"
    ),
    identity = list(
      linkfun = identity,
      linkinv = identity,
      name = "identity"
    ),
    logit = list(
      linkfun = function(x) log(x / (1 - x)),
      linkinv = function(x) 1 / (1 + exp(-x)),
      name = "logit"
    ),
    probit = list(
      linkfun = qnorm,
      linkinv = pnorm,
      name = "probit"
    ),
    cloglog = list(
      linkfun = function(x) log(-log(1 - x)),
      linkinv = function(x) 1 - exp(-exp(x)),
      name = "cloglog"
    ),
    cli::cli_abort("Unknown link function: {.val {link}}")
  )
}

# -----------------------------------------------------------------------------
# Family-specific utilities
# -----------------------------------------------------------------------------

#' Check if family handles structural zeros
#'
#' @param family A simp_family object
#' @return Logical
#' @noRd
family_handles_zeros <- function(family) {
  if (!inherits(family, "simp_family")) {
    cli::cli_abort("{.arg family} must be a {.cls simp_family}")
  }
  isTRUE(family$handles_zeros)
}

#' Get default priors for a family
#'
#' @param family A simp_family object
#' @return A simp_priors object
#' @noRd
family_default_priors <- function(family) {
  if (!inherits(family, "simp_family")) {
    cli::cli_abort("{.arg family} must be a {.cls simp_family}")
  }

  default_priors(family$family)
}
