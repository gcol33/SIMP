# Main model specification function for SIMP
# Creates model objects without fitting (fitting comes in v0.3.0)

#' Check if currently modeling with SIMP
#'
#' @return Logical; TRUE if you are indeed simping
#' @export
#' @examples
#' simping()
simping <- function() {
 TRUE
}

#' Specify a SIMP model
#'
#' Create a model specification for compositional data. This function validates
#' the formula, data, and priors, and constructs the model object. The model
#' is not fitted until `fit()` is called (v0.3.0+).
#'
#' @param formula A model formula. The response should be a `simp_composition`
#'   or `simp_cover` object, or a matrix/data frame of compositional data.
#'   Random effects use lme4-style syntax: `(1 | group)`, `(1 + x | group)`, etc.
#' @param data A data frame containing the predictor variables
#' @param family A family object specifying the likelihood. See
#'   [simp_dirichlet()], [simp_logistic_normal()], [simp_hurdle_dirichlet()],
#'   [simp_stick_breaking()].
#' @param prior A `simp_priors` object, or NULL to use family defaults
#' @param validate Logical; validate the response data (default TRUE)
#' @return A `simp_model` object
#' @export
#' @examples
#' # Create some example data
#' set.seed(42)
#' n <- 50
#' site <- factor(rep(1:10, each = 5))
#' temp <- rnorm(n)
#'
#' # Compositional response (3 components)
#' props <- matrix(runif(n * 3), ncol = 3)
#' props <- props / rowSums(props)
#' cover <- simp_composition(props, components = c("grass", "shrub", "bare"))
#'
#' # Create data frame
#' dat <- data.frame(temp = temp, site = site)
#'
#' # Basic model specification
#' model <- simp(cover ~ temp, data = dat)
#' print(model)
#'
#' # With random effects
#' model <- simp(cover ~ temp + (1 | site), data = dat)
#'
#' # With custom family and priors
#' model <- simp(
#'   cover ~ temp + (1 | site),
#'   data = dat,
#'   family = simp_dirichlet(),
#'   prior = simp_priors(beta = prior_normal(0, 5))
#' )
simp <- function(formula, data, family = simp_dirichlet(),
                 prior = NULL, validate = TRUE) {
  call <- match.call()

  # Validate inputs
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg formula} must be a formula")
  }

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame")
  }

  if (!inherits(family, "simp_family")) {
    cli::cli_abort(c(
      "{.arg family} must be a {.cls simp_family} object",
      "i" = "Use {.fn simp_dirichlet}, {.fn simp_logistic_normal}, etc."
    ))
  }

  # Set default priors if not provided
  if (is.null(prior)) {
    prior <- family_default_priors(family)
  } else if (!inherits(prior, "simp_priors")) {
    cli::cli_abort(c(
      "{.arg prior} must be a {.cls simp_priors} object",
      "i" = "Use {.fn simp_priors} to specify priors"
    ))
  }

  # Extract and validate response
  response <- extract_model_response(formula, data, validate = validate)

  # Check zeros compatibility with family
  if (response$has_zeros && !family_handles_zeros(family)) {
    cli::cli_warn(c(
      "Response contains zeros but {.val {family$family}} family does not handle zeros",
      "i" = "Consider using {.fn simp_hurdle_dirichlet} or {.fn impute_zeros}"
    ))
  }

  # Parse formula
  parsed <- parse_formula(formula, data)

  # Build fixed effects structure
  fixed <- build_fixed_effects(parsed, data)

  # Build random effects structure
  random <- build_random_effects(parsed, data)

  # Create model object
  model <- new_simp_model(
    formula = formula,
    data = data,
    family = family,
    prior = prior,
    response = response,
    fixed = fixed,
    random = random,
    call = call
  )

  model
}

# -----------------------------------------------------------------------------
# Response extraction
# -----------------------------------------------------------------------------

#' Extract and validate model response
#' @noRd
extract_model_response <- function(formula, data, validate = TRUE) {
  if (length(formula) < 3) {
    cli::cli_abort("Formula must have a response (left-hand side)")
  }

  lhs <- formula[[2]]
  response_name <- deparse(lhs)

  # Try to evaluate the response in the data environment
  response <- tryCatch(
    eval(lhs, envir = data, enclos = parent.frame(2)),
    error = function(e) {
      cli::cli_abort(c(
        "Could not find response variable",
        "x" = "Response {.val {response_name}} not found in data",
        "i" = "The response should be a {.cls simp_composition} or matrix"
      ))
    }
  )

  # Handle different response types
  if (inherits(response, "simp_composition")) {
    components <- if (inherits(response, "simp_cover")) {
      attr(response, "species")
    } else {
      attr(response, "components")
    }

    list(
      name = response_name,
      data = response,
      n_obs = attr(response, "n_obs"),
      n_components = attr(response, "n_components"),
      component_names = components,
      has_zeros = attr(response, "has_zeros"),
      type = if (inherits(response, "simp_cover")) "cover" else "composition"
    )
  } else if (is.matrix(response) || is.data.frame(response)) {
    # Convert to simp_composition
    mat <- as_matrix(response)

    if (validate) {
      validation <- validate_composition(mat, verbose = FALSE)
      if (!validation$valid) {
        cli::cli_abort(c(
          "Response is not valid compositional data",
          "i" = "Use {.fn simp_composition} with {.code validate = FALSE} to skip checks"
        ))
      }
    }

    comp_names <- colnames(mat)
    if (is.null(comp_names)) {
      comp_names <- paste0("comp", seq_len(ncol(mat)))
    }

    has_zeros <- any(mat == 0, na.rm = TRUE)

    list(
      name = response_name,
      data = mat,
      n_obs = nrow(mat),
      n_components = ncol(mat),
      component_names = comp_names,
      has_zeros = has_zeros,
      type = "matrix"
    )
  } else {
    cli::cli_abort(c(
      "Response must be compositional data",
      "x" = "Got {.cls {class(response)}}",
      "i" = "Use {.fn simp_composition} or {.fn simp_cover} to create response"
    ))
  }
}

# -----------------------------------------------------------------------------
# Fixed effects construction
# -----------------------------------------------------------------------------

#' Build fixed effects structure
#' @noRd
build_fixed_effects <- function(parsed, data) {
  # Get terms without random effects
  fixed_formula <- parsed$fixed$formula
  terms <- parsed$fixed$terms

  # Check for special terms (varying, shared)
  special <- extract_special_terms(terms)

  # Build design matrix
  if (length(terms) == 0 || (length(terms) == 1 && terms[1] == "1")) {
    # Intercept only
    X <- matrix(1, nrow = nrow(data), ncol = 1)
    colnames(X) <- "(Intercept)"
    terms <- "(Intercept)"
  } else {
    # Use model.matrix to build design matrix
    mm_formula <- fixed_formula
    X <- tryCatch(
      stats::model.matrix(mm_formula, data = data),
      error = function(e) {
        cli::cli_abort(c(
          "Could not build design matrix",
          "x" = e$message
        ))
      }
    )
    terms <- colnames(X)
  }

  list(
    formula = fixed_formula,
    terms = terms,
    X = X,
    n_terms = ncol(X),
    special = special
  )
}

#' Extract special terms (varying, shared)
#' @noRd
extract_special_terms <- function(terms) {
  varying <- character(0)
  shared <- character(0)

  for (term in terms) {
    if (grepl("^varying\\(", term)) {
      var <- gsub("varying\\(([^)]+)\\)", "\\1", term)
      varying <- c(varying, var)
    } else if (grepl("^shared\\(", term)) {
      var <- gsub("shared\\(([^)]+)\\)", "\\1", term)
      shared <- c(shared, var)
    }
  }

  list(
    varying = varying,
    shared = shared,
    has_special = length(varying) > 0 || length(shared) > 0
  )
}

# -----------------------------------------------------------------------------
# Random effects construction
# -----------------------------------------------------------------------------

#' Build random effects structure
#' @noRd
build_random_effects <- function(parsed, data) {
  if (!parsed$has_random) {
    return(list())
  }

  re_list <- lapply(parsed$random, function(re) {
    # Get grouping variable
    group_var <- re$group

    # Handle interaction groups (from nested effects)
    if (grepl(":", group_var)) {
      parts <- strsplit(group_var, ":")[[1]]
      group_factor <- interaction(data[parts], drop = TRUE)
    } else {
      if (!group_var %in% names(data)) {
        cli::cli_abort(c(
          "Grouping variable not found in data",
          "x" = "Variable {.val {group_var}} not in data"
        ))
      }
      group_factor <- as.factor(data[[group_var]])
    }

    n_levels <- nlevels(group_factor)

    # Build Z matrix (design matrix for random effects)
    Z <- build_re_design_matrix(re, data, group_factor)

    list(
      group = re$group,
      terms = re$terms,
      has_intercept = re$has_intercept,
      correlated = re$correlated,
      n_levels = n_levels,
      levels = levels(group_factor),
      group_factor = group_factor,
      Z = Z,
      nested_in = re$nested_in
    )
  })

  # Name by group
  names(re_list) <- vapply(re_list, function(x) x$group, character(1))

  re_list
}

#' Build design matrix for random effects
#' @noRd
build_re_design_matrix <- function(re, data, group_factor) {
  n <- nrow(data)
  n_levels <- nlevels(group_factor)
  group_idx <- as.integer(group_factor)

  # Number of terms (intercept + slopes)
  n_terms <- if (re$has_intercept) 1L else 0L
  n_terms <- n_terms + length(re$terms)

  if (n_terms == 0) {
    cli::cli_abort("Random effect must have at least one term")
  }

  # Build matrix
  # For now, create a simple indicator matrix for intercepts

  # Full sparse matrix construction will be in v0.3.0 for efficiency

  if (re$has_intercept && length(re$terms) == 0) {
    # Random intercept only: indicator matrix
    Z <- matrix(0, nrow = n, ncol = n_levels)
    for (i in seq_len(n)) {
      Z[i, group_idx[i]] <- 1
    }
    colnames(Z) <- paste0(re$group, "_", levels(group_factor))
  } else {
    # Random slopes: more complex structure
    # Each level gets a column for each term
    Z <- matrix(0, nrow = n, ncol = n_levels * n_terms)

    col_names <- character(n_levels * n_terms)
    col <- 1

    for (j in seq_len(n_levels)) {
      if (re$has_intercept) {
        for (i in seq_len(n)) {
          if (group_idx[i] == j) {
            Z[i, col] <- 1
          }
        }
        col_names[col] <- paste0(re$group, "_", levels(group_factor)[j], "_int")
        col <- col + 1
      }

      for (term in re$terms) {
        if (!term %in% names(data)) {
          cli::cli_abort("Random slope variable {.val {term}} not in data")
        }
        x <- data[[term]]
        for (i in seq_len(n)) {
          if (group_idx[i] == j) {
            Z[i, col] <- x[i]
          }
        }
        col_names[col] <- paste0(re$group, "_", levels(group_factor)[j], "_", term)
        col <- col + 1
      }
    }

    colnames(Z) <- col_names
  }

  Z
}
