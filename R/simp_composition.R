# S3 class for compositional data
# Core data structure for SIMP package

#' Create a composition object
#'
#' Constructor for the `simp_composition` class, representing data where each
#' row sums to 1 and all values are in \[0, 1\].
#'
#' @param x A matrix or data frame of compositional data
#' @param components Character vector of component names. If NULL, uses column
#'   names from x or generates default names.
#' @param validate Logical; validate the composition (default TRUE)
#' @param normalize Logical; normalize rows to sum to 1 if needed (default FALSE)
#' @return A `simp_composition` object
#' @export
#' @examples
#' # From matrix
#' x <- matrix(c(0.3, 0.5, 0.2,
#'               0.4, 0.4, 0.2,
#'               0.5, 0.3, 0.2), nrow = 3, byrow = TRUE)
#' comp <- simp_composition(x, components = c("sp1", "sp2", "sp3"))
#'
#' # From data frame
#' df <- data.frame(grass = c(0.3, 0.4), shrub = c(0.5, 0.4), tree = c(0.2, 0.2))
#' comp <- simp_composition(df)
simp_composition <- function(x, components = NULL, validate = TRUE,
                             normalize = FALSE) {
  # Convert to matrix
  mat <- as_matrix(x)

  # Set component names
  if (is.null(components)) {
    if (!is.null(colnames(mat))) {
      components <- colnames(mat)
    } else {
      components <- paste0("comp", seq_len(ncol(mat)))
    }
  }

  check_character(components, "components")

  if (length(components) != ncol(mat)) {
    cli::cli_abort(c(
      "Length of {.arg components} must match number of columns",
      "x" = "Got {length(components)} names for {ncol(mat)} columns"
    ))
  }

  colnames(mat) <- components

  # Normalize if requested
  if (normalize) {
    row_sums <- rowSums(mat, na.rm = TRUE)
    mat <- mat / row_sums
  }

  # Validate
  if (validate) {
    validation <- validate_composition(mat, verbose = FALSE)
    if (!validation$valid) {
      issues <- character(0)
      if (validation$issues$has_negatives) {
        issues <- c(issues, "negative values")
      }
      if (validation$issues$has_above_one) {
        issues <- c(issues, "values > 1")
      }
      if (validation$row_sums$n_issues > 0) {
        issues <- c(issues, sprintf("%d rows don't sum to 1", validation$row_sums$n_issues))
      }
      cli::cli_abort(c(
        "Invalid composition",
        "x" = paste(issues, collapse = ", "),
        "i" = "Use {.code validate = FALSE} to skip validation",
        "i" = "Use {.code normalize = TRUE} to auto-normalize row sums"
      ))
    }
  }

  # Compute zero info
  zero_matrix <- mat == 0
  has_zeros <- any(zero_matrix)

  new_simp_composition(
    mat,
    components = components,
    n_components = ncol(mat),
    n_obs = nrow(mat),
    has_zeros = has_zeros,
    zero_pattern = if (has_zeros) zero_matrix else NULL
  )
}

#' Low-level constructor for simp_composition
#' @noRd
new_simp_composition <- function(x, components, n_components, n_obs,
                                  has_zeros, zero_pattern) {
  structure(
    x,
    components = components,
    n_components = n_components,
    n_obs = n_obs,
    has_zeros = has_zeros,
    zero_pattern = zero_pattern,
    class = c("simp_composition", "matrix", "array")
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.simp_composition <- function(x, ..., n = 6) {
  n_obs <- attr(x, "n_obs")
  n_comp <- attr(x, "n_components")
  components <- attr(x, "components")

  cli::cli_h3("simp_composition")
  cat(sprintf("%d observations x %d components\n", n_obs, n_comp))
  cat(sprintf("Components: %s\n", paste(components, collapse = ", ")))

  if (attr(x, "has_zeros")) {
    n_zeros <- sum(attr(x, "zero_pattern"))
    pct <- 100 * n_zeros / (n_obs * n_comp)
    cat(sprintf("Zeros: %d (%.1f%%)\n", n_zeros, pct))
  }

  # Show preview
  cat("\n")
  if (n_obs > n) {
    print(unclass(x)[1:n, , drop = FALSE])
    cat(sprintf("... and %d more rows\n", n_obs - n))
  } else {
    print(unclass(x))
  }

  invisible(x)
}

#' @export
summary.simp_composition <- function(object, ...) {
  n_obs <- attr(object, "n_obs")
  n_comp <- attr(object, "n_components")
  components <- attr(object, "components")

  # Component summaries
  comp_summary <- data.frame(
    component = components,
    min = apply(object, 2, min, na.rm = TRUE),
    q25 = apply(object, 2, quantile, 0.25, na.rm = TRUE),
    median = apply(object, 2, median, na.rm = TRUE),
    mean = apply(object, 2, mean, na.rm = TRUE),
    q75 = apply(object, 2, quantile, 0.75, na.rm = TRUE),
    max = apply(object, 2, max, na.rm = TRUE),
    zeros = colSums(object == 0, na.rm = TRUE),
    row.names = NULL
  )

  structure(
    list(
      n_obs = n_obs,
      n_components = n_comp,
      components = comp_summary,
      has_zeros = attr(object, "has_zeros")
    ),
    class = "summary.simp_composition"
  )
}

#' @export
print.summary.simp_composition <- function(x, ...) {
  cli::cli_h3("simp_composition Summary")
  cat(sprintf("Observations: %d\n", x$n_obs))
  cat(sprintf("Components: %d\n\n", x$n_components))

  print(x$components, row.names = FALSE)

  invisible(x)
}

#' @export
as.data.frame.simp_composition <- function(x, ...) {
  df <- as.data.frame(unclass(x))
  names(df) <- attr(x, "components")
  df
}

#' @export
as.matrix.simp_composition <- function(x, ...) {
  mat <- unclass(x)
  attr(mat, "components") <- NULL
  attr(mat, "n_components") <- NULL
  attr(mat, "n_obs") <- NULL
  attr(mat, "has_zeros") <- NULL
  attr(mat, "zero_pattern") <- NULL
  class(mat) <- c("matrix", "array")
  mat
}

#' @export
`[.simp_composition` <- function(x, i, j, drop = FALSE) {
  components <- attr(x, "components")
  zero_pattern <- attr(x, "zero_pattern")

  # Get underlying matrix
  mat <- unclass(x)

  # Subset
  if (missing(j)) {
    result <- mat[i, , drop = drop]
    new_components <- components
    if (!is.null(zero_pattern)) {
      zero_pattern <- zero_pattern[i, , drop = drop]
    }
  } else {
    result <- mat[i, j, drop = drop]
    new_components <- components[j]
    if (!is.null(zero_pattern)) {
      zero_pattern <- zero_pattern[i, j, drop = drop]
    }
  }

  # Return matrix if dropped to vector or single value
  if (!is.matrix(result)) {
    return(result)
  }

  # Rebuild composition
  new_simp_composition(
    result,
    components = new_components,
    n_components = ncol(result),
    n_obs = nrow(result),
    has_zeros = any(zero_pattern),
    zero_pattern = zero_pattern
  )
}

#' @export
rbind.simp_composition <- function(..., deparse.level = 1) {
  objs <- list(...)

  # Check all are simp_composition
  if (!all(vapply(objs, inherits, logical(1), "simp_composition"))) {
    cli::cli_abort("All arguments must be {.cls simp_composition} objects")
  }

  # Check components match
  components <- attr(objs[[1]], "components")
  for (i in seq_along(objs)[-1]) {
    if (!identical(attr(objs[[i]], "components"), components)) {
      cli::cli_abort("All compositions must have the same components")
    }
  }

  # Combine matrices
  mats <- lapply(objs, as.matrix)
  combined <- do.call(rbind, mats)

  # Combine zero patterns
  zero_patterns <- lapply(objs, attr, "zero_pattern")
  if (all(vapply(zero_patterns, is.null, logical(1)))) {
    combined_zeros <- NULL
    has_zeros <- FALSE
  } else {
    # Replace NULL with FALSE matrix
    zero_patterns <- lapply(seq_along(objs), function(i) {
      zp <- zero_patterns[[i]]
      if (is.null(zp)) {
        matrix(FALSE, nrow = attr(objs[[i]], "n_obs"),
               ncol = attr(objs[[i]], "n_components"))
      } else {
        zp
      }
    })
    combined_zeros <- do.call(rbind, zero_patterns)
    has_zeros <- any(combined_zeros)
  }

  new_simp_composition(
    combined,
    components = components,
    n_components = length(components),
    n_obs = nrow(combined),
    has_zeros = has_zeros,
    zero_pattern = combined_zeros
  )
}

#' @export
dim.simp_composition <- function(x) {
  c(attr(x, "n_obs"), attr(x, "n_components"))
}

#' @export
dimnames.simp_composition <- function(x) {
  list(rownames(unclass(x)), attr(x, "components"))
}

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

#' Check if object is a simp_composition
#'
#' @param x Object to check
#' @return Logical
#' @export
is_simp_composition <- function(x) {
  inherits(x, "simp_composition")
}

#' Get component names
#'
#' @param x A simp_composition object
#' @return Character vector of component names
#' @export
components <- function(x) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition}")
  }
  attr(x, "components")
}

#' Get number of components
#'
#' @param x A simp_composition object
#' @return Integer
#' @export
n_components <- function(x) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition}")
  }
  attr(x, "n_components")
}

#' Get number of observations
#'
#' @param x A simp_composition object
#' @return Integer
#' @export
n_obs <- function(x) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition}")
  }
  attr(x, "n_obs")
}
