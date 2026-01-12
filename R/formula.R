# Formula parsing for SIMP models
# Ported from ratiod/R/formula.R with extensions for compositional responses

# -----------------------------------------------------------------------------
# Main parsing function
# -----------------------------------------------------------------------------

#' Parse a SIMP model formula
#'
#' Extracts fixed effects, random effects, and response information from a
#' model formula. Supports all standard random effect structures.
#'
#' @param formula A model formula
#' @param data Optional data frame for validation
#' @return A list with parsed formula components
#' @export
#' @examples
#' # Fixed effects only
#' parse_formula(y ~ x1 + x2)
#'
#' # Random intercepts
#' parse_formula(y ~ x1 + (1 | group))
#'
#' # Random slopes
#' parse_formula(y ~ x1 + (1 + x1 | group))
#'
#' # Nested random effects
#' parse_formula(y ~ x1 + (1 | site/plot))
#'
#' # Crossed random effects
#' parse_formula(y ~ x1 + (1 | site) + (1 | year))
#'
#' # Uncorrelated random effects
#' parse_formula(y ~ x1 + (1 + x1 || group))
parse_formula <- function(formula, data = NULL) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg formula} must be a formula")
  }

  # Extract response
  response <- extract_response(formula)

  # Extract random effects (before removing from formula)
  re_terms <- extract_random_effects(formula)

  # Get fixed effects formula (remove RE terms)
  fixed_formula <- remove_random_effects(formula)

  # Parse fixed effects
  fixed_terms <- if (length(fixed_formula) == 3) {
    all.vars(fixed_formula[[3]])
  } else {
    character(0)
  }

  # Validate against data if provided
  if (!is.null(data)) {
    validate_formula_data(formula, data, fixed_terms, re_terms)
  }

  structure(
    list(
      formula = formula,
      response = response,
      fixed = list(
        formula = fixed_formula,
        terms = fixed_terms
      ),
      random = re_terms,
      has_random = length(re_terms) > 0
    ),
    class = "simp_parsed_formula"
  )
}

# -----------------------------------------------------------------------------
# Response extraction
# -----------------------------------------------------------------------------

#' Extract response from formula
#' @noRd
extract_response <- function(formula) {
  if (length(formula) < 3) {
    return(NULL)
  }

  lhs <- formula[[2]]

  # Check for cbind() response (multivariate)
  if (is.call(lhs) && identical(lhs[[1]], as.name("cbind"))) {
    vars <- as.character(lhs[-1])
    return(list(
      type = "multivariate",
      variables = vars,
      n_components = length(vars)
    ))
  }

  # Single response
  list(
    type = "univariate",
    variables = deparse(lhs),
    n_components = 1
  )
}

# -----------------------------------------------------------------------------
# Random effects extraction
# -----------------------------------------------------------------------------

#' Extract random effect terms from formula
#'
#' Parses random effect specifications including:
#' - `(1 | group)` - random intercepts
#' - `(1 + x | group)` - correlated random intercept and slope
#' - `(0 + x | group)` - random slope only (no intercept)
#' - `(1 + x || group)` - uncorrelated random effects
#' - `(1 | g1/g2)` - nested random effects
#'
#' @param formula A model formula
#' @return A list of random effect specifications
#' @export
#' @examples
#' extract_random_effects(y ~ x + (1 | group))
#' extract_random_effects(y ~ x + (1 + x | group))
#' extract_random_effects(y ~ x + (1 | site/plot))
extract_random_effects <- function(formula) {
  formula_str <- deparse(formula, width.cutoff = 500)
  formula_str <- paste(formula_str, collapse = "")

  # Find all RE terms using balanced parenthesis matching
  re_strings <- extract_re_terms_balanced(formula_str)

  if (length(re_strings) == 0) {
    return(list())
  }

  # Parse each RE term
  re_list <- lapply(re_strings, parse_single_re)

  # Expand nested effects
  re_list <- expand_nested_effects(re_list)

  # Name by grouping factor
  names(re_list) <- vapply(re_list, function(x) x$group, character(1))

  re_list
}

#' Extract RE terms with balanced parenthesis matching
#' @noRd
extract_re_terms_balanced <- function(formula_str) {
  # Find positions of ( that might start a RE term
  # RE terms have pattern: ( ... | ... )

  results <- character(0)
  chars <- strsplit(formula_str, "")[[1]]
  n <- length(chars)
  i <- 1

  while (i <= n) {
    if (chars[i] == "(") {
      # Find matching close paren
      depth <- 1
      start <- i
      j <- i + 1

      while (j <= n && depth > 0) {
        if (chars[j] == "(") depth <- depth + 1
        if (chars[j] == ")") depth <- depth - 1
        j <- j + 1
      }

      if (depth == 0) {
        content <- substr(formula_str, start + 1, j - 2)
        # Check if it contains a bar (RE indicator)
        if (grepl("\\|", content)) {
          results <- c(results, content)
        }
      }
      i <- j
    } else {
      i <- i + 1
    }
  }

  results
}

#' Parse a single random effect term
#' @noRd
parse_single_re <- function(re_string) {
  # Handle uncorrelated syntax: ||
  uncorrelated <- grepl("\\|\\|", re_string)

  # Split on | or ||
  parts <- if (uncorrelated) {
    strsplit(re_string, "\\|\\|")[[1]]
  } else {
    strsplit(re_string, "\\|")[[1]]
  }

  if (length(parts) != 2) {
    cli::cli_abort(c(
      "Invalid random effect specification: {.val {re_string}}",
      "i" = "Expected format: (terms | group) or (terms || group)"
    ))
  }

  lhs <- trimws(parts[1])
  rhs <- trimws(parts[2])

  # Parse the left-hand side (terms)
  terms <- parse_re_terms(lhs)

  # Check for nested structure in grouping factor
  is_nested <- grepl("/", rhs)

  list(
    raw = re_string,
    group = gsub("/.*", "", rhs),  # First group for nested
    full_group = rhs,
    terms = terms$terms,
    has_intercept = terms$has_intercept,
    correlated = !uncorrelated,
    nested = is_nested,
    nested_in = if (is_nested) gsub(".*/", "", rhs) else NULL
  )
}

#' Parse terms from LHS of random effect
#' @noRd
parse_re_terms <- function(lhs) {
  # Handle "1" for intercept only
  if (lhs == "1") {
    return(list(terms = character(0), has_intercept = TRUE))
  }

  # Handle "0 + x" for slope only
  if (grepl("^0\\s*\\+", lhs)) {
    terms <- trimws(strsplit(gsub("^0\\s*\\+\\s*", "", lhs), "\\+")[[1]])
    return(list(terms = terms, has_intercept = FALSE))
  }

  # Handle "1 + x" for intercept and slope
  if (grepl("^1\\s*\\+", lhs)) {
    terms <- trimws(strsplit(gsub("^1\\s*\\+\\s*", "", lhs), "\\+")[[1]])
    return(list(terms = terms, has_intercept = TRUE))
  }

  # Just variable names
  terms <- trimws(strsplit(lhs, "\\+")[[1]])
  list(terms = terms, has_intercept = FALSE)
}

#' Expand nested random effects into separate terms
#' @noRd
expand_nested_effects <- function(re_list) {
  expanded <- list()

  for (re in re_list) {
    if (re$nested) {
      # (1 | g1/g2) expands to (1 | g1) + (1 | g1:g2)
      outer_group <- re$group
      inner_group <- re$nested_in

      # Outer group
      outer_re <- re
      outer_re$nested <- FALSE
      outer_re$nested_in <- NULL
      outer_re$full_group <- outer_group
      expanded[[length(expanded) + 1]] <- outer_re

      # Inner group (interaction)
      inner_re <- re
      inner_re$group <- paste0(outer_group, ":", inner_group)
      inner_re$nested <- FALSE
      inner_re$nested_in <- outer_group
      inner_re$full_group <- inner_re$group
      expanded[[length(expanded) + 1]] <- inner_re
    } else {
      expanded[[length(expanded) + 1]] <- re
    }
  }

  expanded
}

# -----------------------------------------------------------------------------
# Formula manipulation
# -----------------------------------------------------------------------------

#' Remove random effects from formula
#'
#' @param formula A model formula
#' @return Formula with random effect terms removed
#' @export
#' @examples
#' remove_random_effects(y ~ x + (1 | group))
#' # Returns: y ~ x
remove_random_effects <- function(formula) {
  formula_str <- deparse(formula, width.cutoff = 500)
  formula_str <- paste(formula_str, collapse = "")

  # Remove RE terms with balanced parenthesis matching
  re_strings <- extract_re_terms_balanced(formula_str)

  for (re_str in re_strings) {
    # Escape special regex characters in the RE string
    escaped <- gsub("([\\(\\)\\|\\+])", "\\\\\\1", re_str)
    # Remove the term including surrounding parentheses and plus signs
    pattern <- paste0("\\s*\\+?\\s*\\(", escaped, "\\)\\s*\\+?\\s*")
    formula_str <- gsub(pattern, " + ", formula_str)
  }

  # Clean up formula string
  formula_str <- gsub("~\\s*\\+", "~", formula_str)
  formula_str <- gsub("\\+\\s*$", "", formula_str)
  formula_str <- gsub("\\s+", " ", formula_str)
  formula_str <- trimws(formula_str)

  # Handle case where all predictors were removed
  if (grepl("~\\s*$", formula_str)) {
    formula_str <- paste0(formula_str, "1")
  }

  as.formula(formula_str, env = environment(formula))
}

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

#' Validate formula against data
#' @noRd
validate_formula_data <- function(formula, data, fixed_terms, re_terms) {
  data_names <- names(data)

  # Check fixed effects
  missing_fixed <- setdiff(fixed_terms, data_names)
  if (length(missing_fixed) > 0) {
    cli::cli_abort(c(
      "Variables not found in data",
      "x" = "Missing: {.val {missing_fixed}}"
    ))
  }

  # Check random effect grouping factors
  for (re in re_terms) {
    groups <- strsplit(re$group, ":")[[1]]
    missing_groups <- setdiff(groups, data_names)
    if (length(missing_groups) > 0) {
      cli::cli_abort(c(
        "Grouping factor not found in data",
        "x" = "Missing: {.val {missing_groups}}",
        "i" = "From random effect: ({re$raw})"
      ))
    }

    # Check random slope terms
    missing_slopes <- setdiff(re$terms, data_names)
    if (length(missing_slopes) > 0) {
      cli::cli_abort(c(
        "Random slope variable not found in data",
        "x" = "Missing: {.val {missing_slopes}}",
        "i" = "From random effect: ({re$raw})"
      ))
    }
  }

  invisible(TRUE)
}

#' Validate a SIMP formula
#'
#' Check that a formula is valid for SIMP models.
#'
#' @param formula A model formula
#' @param data Data frame
#' @return TRUE invisibly if valid, error otherwise
#' @export
validate_formula <- function(formula, data) {
  parsed <- parse_formula(formula, data)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Print methods
# -----------------------------------------------------------------------------

#' @export
print.simp_parsed_formula <- function(x, ...) {
  cli::cli_h3("Parsed SIMP Formula")

  cat("Response:", x$response$variables, "\n")
  cat("  Type:", x$response$type, "\n")

  cat("\nFixed effects:\n")
  if (length(x$fixed$terms) > 0) {
    cat(" ", paste(x$fixed$terms, collapse = " + "), "\n")
  } else {
    cat("  (intercept only)\n")
  }

  if (x$has_random) {
    cat("\nRandom effects:\n")
    for (re in x$random) {
      terms_str <- if (re$has_intercept && length(re$terms) > 0) {
        paste("1 +", paste(re$terms, collapse = " + "))
      } else if (re$has_intercept) {
        "1"
      } else {
        paste("0 +", paste(re$terms, collapse = " + "))
      }
      cor_str <- if (re$correlated) "|" else "||"
      cat(sprintf("  (%s %s %s)\n", terms_str, cor_str, re$full_group))
    }
  }

  invisible(x)
}
