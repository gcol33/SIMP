# Normalization functions for compositional data
# Functions to ensure data satisfies compositional constraints

#' Normalize cover data to sum to 1
#'
#' Rescale rows of compositional data so they sum to 1.
#'
#' @param x A matrix, data frame, or simp_composition/simp_cover object
#' @param method Normalization method:
#'   - `"proportional"`: Scale all values proportionally (default)
#'   - `"truncate"`: Truncate values > 1 before scaling
#' @return Object of the same class as input with normalized values
#' @export
#' @examples
#' # Data that doesn't sum to 1
#' x <- matrix(c(30, 50, 20,
#'               40, 40, 25), nrow = 2, byrow = TRUE)
#' normalize_cover(x)
normalize_cover <- function(x, method = c("proportional", "truncate")) {
  method <- match.arg(method)

  # Get matrix and metadata
  is_cover <- inherits(x, "simp_cover")
  is_comp <- inherits(x, "simp_composition")

  mat <- if (is_cover || is_comp) as.matrix(x) else as_matrix(x)

  # Check for negative values
  if (any(mat < 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "Cannot normalize data with negative values",
      "i" = "Found {sum(mat < 0, na.rm = TRUE)} negative value(s)"
    ))
  }

  # Apply method
  if (method == "truncate") {
    # First truncate values > 1
    mat[mat > 1] <- 1
  }

  # Proportional scaling
  row_sums <- rowSums(mat, na.rm = TRUE)

  # Check for zero rows
  zero_rows <- row_sums == 0
  if (any(zero_rows)) {
    cli::cli_abort(c(
      "Cannot normalize rows that sum to zero",
      "x" = "{sum(zero_rows)} row(s) have zero sum"
    ))
  }

  mat_normalized <- mat / row_sums

  # Return appropriate type
  if (is_cover) {
    simp_cover(mat_normalized,
               species = attr(x, "species"),
               layers = attr(x, "layers"),
               residual = attr(x, "residual"),
               validate = FALSE)
  } else if (is_comp) {
    simp_composition(mat_normalized,
                     components = attr(x, "components"),
                     validate = FALSE)
  } else {
    mat_normalized
  }
}

#' Add residual component to close composition
#'
#' If compositional data doesn't sum to 1, add a residual component that
#' captures the remainder.
#'
#' @param x A matrix, data frame, or simp_composition/simp_cover object
#' @param name Name for the residual component (default "other")
#' @param position Where to add the residual: "last" (default) or "first"
#' @return Object with added residual component
#' @export
#' @examples
#' # Open composition (sums to 0.8)
#' x <- matrix(c(0.3, 0.3, 0.2,
#'               0.4, 0.2, 0.2), nrow = 2, byrow = TRUE)
#' colnames(x) <- c("sp1", "sp2", "sp3")
#'
#' # Close by adding residual
#' add_residual(x, name = "bare")
add_residual <- function(x, name = "other", position = c("last", "first")) {
  position <- match.arg(position)

  # Get matrix and metadata
  is_cover <- inherits(x, "simp_cover")
  is_comp <- inherits(x, "simp_composition")

  mat <- if (is_cover || is_comp) as.matrix(x) else as_matrix(x)
  components <- if (is_cover) {
    attr(x, "species")
  } else if (is_comp) {
    attr(x, "components")
  } else {
    colnames(mat)
  }

  if (is.null(components)) {
    components <- paste0("comp", seq_len(ncol(mat)))
  }

  # Check if residual name already exists
  if (name %in% components) {
    cli::cli_abort(c(
      "Residual name already exists in components",
      "x" = "{.val {name}} is already a component name",
      "i" = "Choose a different name for the residual"
    ))
  }

  # Calculate residual
  row_sums <- rowSums(mat, na.rm = TRUE)
  residual <- 1 - row_sums

  # Check for negative residuals
  if (any(residual < -1e-6)) {
    n_over <- sum(residual < -1e-6)
    cli::cli_abort(c(
      "Cannot add residual when rows already exceed 1",
      "x" = "{n_over} row(s) sum to more than 1",
      "i" = "Use {.fn normalize_cover} first to rescale"
    ))
  }

  # Clip small negative values to 0
  residual[residual < 0 & residual > -1e-6] <- 0

  # Add residual column
  if (position == "last") {
    mat_new <- cbind(mat, residual)
    new_components <- c(components, name)
  } else {
    mat_new <- cbind(residual, mat)
    new_components <- c(name, components)
  }

  colnames(mat_new) <- new_components

  # Return appropriate type
  if (is_cover) {
    layers <- attr(x, "layers")
    new_layers <- if (!is.null(layers)) {
      if (position == "last") {
        c(layers, "residual")
      } else {
        c("residual", layers)
      }
    } else {
      NULL
    }

    simp_cover(mat_new,
               species = new_components,
               layers = new_layers,
               residual = name,
               validate = FALSE)
  } else if (is_comp) {
    simp_composition(mat_new,
                     components = new_components,
                     validate = FALSE)
  } else {
    mat_new
  }
}

#' Remove residual component from composition
#'
#' Remove a specified residual/other component, returning an open composition.
#'
#' @param x A simp_composition or simp_cover object
#' @param name Name of the residual component to remove. If NULL, uses the
#'   residual attribute for simp_cover objects.
#' @return Object without the residual component
#' @export
#' @examples
#' x <- simp_cover(
#'   matrix(c(0.3, 0.5, 0.2,
#'            0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
#'   species = c("grass", "shrub", "bare"),
#'   residual = "bare"
#' )
#' remove_residual(x)
remove_residual <- function(x, name = NULL) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition} or {.cls simp_cover}")
  }

  # Get residual name
  if (is.null(name)) {
    if (inherits(x, "simp_cover")) {
      name <- attr(x, "residual")
    }
    if (is.null(name)) {
      cli::cli_abort(c(
        "No residual component specified",
        "i" = "Provide {.arg name} or use a {.cls simp_cover} with residual attribute"
      ))
    }
  }

  # Find residual column
  components <- if (inherits(x, "simp_cover")) {
    attr(x, "species")
  } else {
    attr(x, "components")
  }

  idx <- match(name, components)
  if (is.na(idx)) {
    cli::cli_abort(c(
      "Residual component not found",
      "x" = "{.val {name}} is not a component"
    ))
  }

  # Remove column
  x[, -idx, drop = FALSE]
}

#' Convert between compositional representations
#'
#' Transform compositional data between different representations.
#'
#' @param x A simp_composition or simp_cover object
#' @param to Target representation:
#'   - `"alr"`: Additive log-ratio (requires reference component)
#'   - `"clr"`: Centered log-ratio
#'   - `"ilr"`: Isometric log-ratio
#' @param reference For ALR, the index or name of the reference component
#' @return Transformed data as a matrix
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.3, 0.5, 0.2,
#'            0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
#'   components = c("a", "b", "c")
#' )
#'
#' # Centered log-ratio
#' transform_composition(x, to = "clr")
#'
#' # Additive log-ratio with "c" as reference
#' transform_composition(x, to = "alr", reference = "c")
transform_composition <- function(x, to = c("clr", "alr", "ilr"),
                                  reference = NULL) {
  to <- match.arg(to)

  mat <- if (inherits(x, "simp_composition")) as.matrix(x) else as_matrix(x)

  # Check for zeros
  if (any(mat == 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "Log-ratio transforms require non-zero values",
      "i" = "Use {.fn impute_zeros} first to handle zeros"
    ))
  }

  log_mat <- log(mat)
  K <- ncol(mat)

  switch(to,
    clr = {
      # Centered log-ratio: log(x_k) - mean(log(x))
      geometric_mean_log <- rowMeans(log_mat, na.rm = TRUE)
      result <- log_mat - geometric_mean_log
      colnames(result) <- if (inherits(x, "simp_cover")) {
        attr(x, "species")
      } else if (inherits(x, "simp_composition")) {
        attr(x, "components")
      } else {
        colnames(mat)
      }
      result
    },
    alr = {
      # Additive log-ratio: log(x_k / x_ref)
      if (is.null(reference)) {
        reference <- K  # Default to last column
      }
      if (is.character(reference)) {
        components <- if (inherits(x, "simp_cover")) {
          attr(x, "species")
        } else if (inherits(x, "simp_composition")) {
          attr(x, "components")
        } else {
          colnames(mat)
        }
        reference <- match(reference, components)
        if (is.na(reference)) {
          cli::cli_abort("Reference component not found")
        }
      }
      ref_log <- log_mat[, reference]
      result <- log_mat[, -reference, drop = FALSE] - ref_log
      result
    },
    ilr = {
      # Isometric log-ratio using Helmert subcomposition
      # Creates K-1 orthonormal coordinates
      result <- matrix(0, nrow = nrow(mat), ncol = K - 1)
      for (j in seq_len(K - 1)) {
        # Helmert contrast
        result[, j] <- sqrt(j / (j + 1)) *
          (log_mat[, j] - rowMeans(log_mat[, (j + 1):K, drop = FALSE]))
      }
      colnames(result) <- paste0("ilr", seq_len(K - 1))
      result
    }
  )
}
