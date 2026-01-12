# Latent factor specification for SIMP models
# Handles unmeasured confounders via observation-level latent factors

#' Specify latent factors
#'
#' Add latent factors to a SIMP model to account for unmeasured confounders.
#' Latent factors are observation-level random effects that can capture
#' shared variation across components not explained by the fixed effects.
#'
#' @param n Number of latent factors (default 1)
#' @param prior Prior for the factor loadings. Default is a PC prior that
#'   shrinks toward zero factors.
#' @return A `simp_latent` specification object
#' @references
#' Simpson et al. (2017) Penalising Model Component Complexity.
#' Statistical Science 32(1):1-28.
#' @export
#' @examples
#' # Single latent factor with default prior
#' latent_factor()
#'
#' # Multiple factors with custom prior
#' latent_factor(n = 2, prior = prior_pc(u = 1, alpha = 0.01))
#'
#' # Usage in model formula (v0.3.0+)
#' # simp(cover ~ temp + latent_factor(2) + (1 | site), data = dat)
latent_factor <- function(n = 1, prior = NULL) {
  check_scalar(n, "n")

  if (n < 1 || n != floor(n)) {
    cli::cli_abort(c(
      "{.arg n} must be a positive integer",
      "x" = "Got {.val {n}}"
    ))
  }

  # Default prior: PC prior shrinking toward zero variance
 if (is.null(prior)) {
    prior <- prior_pc(u = 1, alpha = 0.01)
  }

  if (!inherits(prior, "simp_prior")) {
    cli::cli_abort(c(
      "{.arg prior} must be a {.cls simp_prior} object",
      "i" = "Use {.fn prior_pc}, {.fn prior_half_normal}, etc."
    ))
  }

  structure(
    list(
      n_factors = as.integer(n),
      prior = prior
    ),
    class = "simp_latent"
  )
}

#' @export
print.simp_latent <- function(x, ...) {
  cat("Latent factor specification\n")
  cat("  Factors:", x$n_factors, "\n")
  cat("  Prior: ")
  print(x$prior)
  invisible(x)
}

#' Check if object is a latent factor specification
#'
#' @param x Object to check
#' @return Logical
#' @export
is_latent <- function(x) {
  inherits(x, "simp_latent")
}

#' Extract latent factor specification from formula
#'
#' @param formula A model formula
#' @return A simp_latent object or NULL if no latent factors specified
#' @noRd
extract_latent_from_formula <- function(formula) {
  formula_str <- deparse(formula, width.cutoff = 500)
  formula_str <- paste(formula_str, collapse = "")

  # Check for latent_factor() call
  if (!grepl("latent_factor\\s*\\(", formula_str)) {
    return(NULL)
  }

  # Extract the latent_factor call
  # This is a simplified parser; full implementation in v0.3.0
  match <- regmatches(
    formula_str,
    regexpr("latent_factor\\s*\\([^)]*\\)", formula_str)
  )

  if (length(match) == 0) {
    return(NULL)
  }

  # Parse and evaluate
  latent_call <- tryCatch(
    eval(parse(text = match)),
    error = function(e) {
      cli::cli_warn("Could not parse latent_factor specification: {e$message}")
      NULL
    }
  )

  latent_call
}

#' Remove latent factor terms from formula
#'
#' @param formula A model formula
#' @return Formula with latent_factor() terms removed
#' @noRd
remove_latent_from_formula <- function(formula) {
  formula_str <- deparse(formula, width.cutoff = 500)
  formula_str <- paste(formula_str, collapse = "")

  # Remove latent_factor(...) terms
  formula_str <- gsub("\\+?\\s*latent_factor\\s*\\([^)]*\\)\\s*\\+?", " + ", formula_str)

  # Clean up
  formula_str <- gsub("~\\s*\\+", "~", formula_str)
  formula_str <- gsub("\\+\\s*$", "", formula_str)
  formula_str <- gsub("\\+\\s*\\+", "+", formula_str)
  formula_str <- gsub("\\s+", " ", formula_str)
  formula_str <- trimws(formula_str)

  # Handle intercept-only case
  if (grepl("~\\s*$", formula_str)) {
    formula_str <- paste0(formula_str, "1")
  }

  as.formula(formula_str, env = environment(formula))
}
