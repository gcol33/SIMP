# S3 class for SIMP model objects
# Contains model specification before and after fitting

# -----------------------------------------------------------------------------
# Constructor
# -----------------------------------------------------------------------------

#' Create a new simp_model object
#'
#' Low-level constructor for model specification objects.
#'
#' @param formula Original model formula
#' @param data Processed data
#' @param family A simp_family object
#' @param prior A simp_priors object
#' @param response Response specification
#' @param fixed Fixed effects specification
#' @param random Random effects specification
#' @param call Original function call
#' @return A simp_model object
#' @noRd
new_simp_model <- function(formula, data, family, prior,
                            response, fixed, random, call = NULL) {
  structure(
    list(
      formula = formula,
      data = data,
      family = family,
      prior = prior,
      response = response,
      fixed = fixed,
      random = random,
      call = call,
      fitted = FALSE,
      fit = NULL
    ),
    class = "simp_model"
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.simp_model <- function(x, ...) {
  cli::cli_h3("SIMP Model Specification")
  cli::cli_rule()

  # Response
  cat("Response:", x$response$name)
  if (x$response$n_components > 1) {
    cat(sprintf(" (%d components)\n", x$response$n_components))
  } else {
    cat("\n")
  }

  # Fixed effects
  cat("\nFixed effects:")
  if (length(x$fixed$terms) > 0) {
    terms_str <- paste(x$fixed$terms, collapse = " + ")
    if (nchar(terms_str) > 50) {
      cat("\n  ", terms_str, "\n")
    } else {
      cat(" ", terms_str, "\n")
    }
  } else {
    cat(" (intercept only)\n")
  }

  # Random effects
  if (length(x$random) > 0) {
    cat("\nRandom effects:\n")
    for (re in x$random) {
      terms_str <- format_re_terms(re)
      cat(sprintf("  - (%s | %s): %d levels",
                  terms_str, re$group, re$n_levels))
      if (re$correlated && length(re$terms) > 0) {
        cat(", correlated")
      }
      cat("\n")
    }
  }

  # Family
  cat("\nFamily:", x$family$family)
  if (x$family$link != "log" && x$family$link != "identity") {
    cat(sprintf(" (link: %s)", x$family$link))
  }
  cat("\n")

  # Data summary
  cat(sprintf("Observations: %d\n", x$response$n_obs))

  # Status
  cat("\nStatus: ")
  if (x$fitted) {
    cli::cli_alert_success("fitted")
  } else {
    cli::cli_alert_info("not fitted")
  }

  invisible(x)
}

#' Format random effect terms for display
#' @noRd
format_re_terms <- function(re) {
  if (re$has_intercept && length(re$terms) > 0) {
    paste("1 +", paste(re$terms, collapse = " + "))
  } else if (re$has_intercept) {
    "1"
  } else {
    paste("0 +", paste(re$terms, collapse = " + "))
  }
}

#' @export
summary.simp_model <- function(object, ...) {
  if (!object$fitted) {
    cli::cli_warn("Model has not been fitted yet")
    print(object)
    return(invisible(object))
  }

  # Summary for fitted model (v0.3.0+)
  structure(
    list(
      call = object$call,
      formula = object$formula,
      family = object$family,
      n_obs = object$response$n_obs,
      n_components = object$response$n_components,
      fixed = object$fixed,
      random = object$random,
      fit = object$fit
    ),
    class = "summary.simp_model"
  )
}

#' @export
print.summary.simp_model <- function(x, ...) {
  cli::cli_h3("SIMP Model Summary")
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  cat("Family:", x$family$family, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Components:", x$n_components, "\n")

  # Coefficients will be added in v0.3.0
  if (!is.null(x$fit)) {
    cat("\n(Coefficient summary available after fitting)\n")
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# Accessor functions
# -----------------------------------------------------------------------------

#' Extract formula from simp_model
#'
#' @param x A simp_model object
#' @param ... Ignored
#' @return The model formula
#' @export
formula.simp_model <- function(x, ...) {
  x$formula
}

#' Check if model is fitted
#'
#' @param x A simp_model object
#' @return Logical
#' @export
is_fitted <- function(x) {
  if (!inherits(x, "simp_model")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_model}")
  }
  x$fitted
}

#' Get number of observations from model
#'
#' @param object A simp_model object
#' @param ... Ignored
#' @return Integer count of observations
#' @export
nobs.simp_model <- function(object, ...) {
  object$response$n_obs
}

#' Get model family
#'
#' @param object A simp_model object
#' @param ... Ignored
#' @return A simp_family object
#' @export
family.simp_model <- function(object, ...) {
  object$family
}

# -----------------------------------------------------------------------------
# Model information
# -----------------------------------------------------------------------------

#' Get model component names
#'
#' @param model A simp_model object
#' @return Character vector of component names
#' @export
model_components <- function(model) {
  if (!inherits(model, "simp_model")) {
    cli::cli_abort("{.arg model} must be a {.cls simp_model}")
  }
  model$response$component_names
}

#' Get number of model parameters
#'
#' @param model A simp_model object
#' @return Named integer vector with parameter counts
#' @export
n_params <- function(model) {

  if (!inherits(model, "simp_model")) {
    cli::cli_abort("{.arg model} must be a {.cls simp_model}")
  }

  K <- model$response$n_components

  # Fixed effects: each predictor has K coefficients (or K-1 for some families)
  n_fixed <- length(model$fixed$terms) * K

  # Random effects
  n_random <- 0
  n_sigma <- 0
  n_cor <- 0

  for (re in model$random) {
    n_re_terms <- if (re$has_intercept) 1L else 0L
    n_re_terms <- n_re_terms + length(re$terms)

    # Random effects per group
    n_random <- n_random + re$n_levels * n_re_terms * K

    # Variance parameters
    n_sigma <- n_sigma + n_re_terms

    # Correlation parameters (if correlated)
    if (re$correlated && n_re_terms > 1) {
      n_cor <- n_cor + n_re_terms * (n_re_terms - 1) / 2
    }
  }

  # Family-specific parameters
  n_family <- switch(model$family$family,
    dirichlet = if (is.null(model$family$concentration)) 1L else 0L,
    logistic_normal = K * (K - 1) / 2,  # correlation matrix
    hurdle_dirichlet = K,  # zero probabilities
    stick_breaking = 0L,
    0L
  )

  c(
    fixed = n_fixed,
    random = n_random,
    sigma = n_sigma,
    correlation = n_cor,
    family = n_family,
    total = n_fixed + n_random + n_sigma + n_cor + n_family
  )
}

#' Check if object is a simp_model
#'
#' @param x Object to check
#' @return Logical
#' @export
is_simp_model <- function(x) {
  inherits(x, "simp_model")
}
