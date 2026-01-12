# Zero handling for compositional data
# Functions for identifying, classifying, and imputing zeros

#' Identify zeros in compositional data
#'
#' Analyze zero patterns in compositional data and classify them as potentially
#' structural (true absence) or sampling (rounded/below detection).
#'
#' @param x A simp_composition or simp_cover object
#' @param structural_threshold Components with zeros in more than this
#'   proportion of observations are likely structural (default 0.8)
#' @return A list with zero pattern information
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.4, 0.4, 0.2, 0,
#'            0.5, 0.3, 0.2, 0,
#'            0.3, 0.5, 0, 0.2), nrow = 3, byrow = TRUE),
#'   components = c("sp1", "sp2", "sp3", "sp4")
#' )
#'
#' identify_zeros(x)
identify_zeros <- function(x, structural_threshold = 0.8) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition} or {.cls simp_cover}")
  }

  mat <- as.matrix(x)
  components <- if (inherits(x, "simp_cover")) {
    attr(x, "species")
  } else {
    attr(x, "components")
  }

  n_obs <- nrow(mat)
  n_comp <- ncol(mat)

  # Zero counts
  zero_counts <- colSums(mat == 0, na.rm = TRUE)
  zero_props <- zero_counts / n_obs

  # Classify zeros
  structural <- zero_props >= structural_threshold
  sampling <- zero_props > 0 & zero_props < structural_threshold

  # Pattern matrix
  zero_matrix <- mat == 0

  # Row-level info
  zeros_per_row <- rowSums(zero_matrix, na.rm = TRUE)

  result <- list(
    n_zeros = sum(zero_counts),
    pct_zeros = 100 * sum(zero_counts) / (n_obs * n_comp),
    per_component = data.frame(
      component = components,
      n_zeros = zero_counts,
      pct_zeros = 100 * zero_props,
      classification = ifelse(structural, "structural",
                              ifelse(sampling, "sampling", "none")),
      row.names = NULL
    ),
    per_row = data.frame(
      row = seq_len(n_obs),
      n_zeros = zeros_per_row
    ),
    structural_components = components[structural],
    sampling_zero_components = components[sampling],
    zero_matrix = zero_matrix
  )

  class(result) <- "simp_zero_analysis"
  result
}

#' @export
print.simp_zero_analysis <- function(x, ...) {
  cli::cli_h3("Zero Analysis")

  cat(sprintf("Total zeros: %d (%.1f%%)\n", x$n_zeros, x$pct_zeros))

  if (length(x$structural_components) > 0) {
    cat("\nStructural zeros (likely true absence):\n")
    cat(" ", paste(x$structural_components, collapse = ", "), "\n")
  }

  if (length(x$sampling_zero_components) > 0) {
    cat("\nSampling zeros (may need imputation):\n")
    cat(" ", paste(x$sampling_zero_components, collapse = ", "), "\n")
  }

  cat("\nPer-component summary:\n")
  print(x$per_component[x$per_component$n_zeros > 0, ], row.names = FALSE)

  invisible(x)
}

#' Impute zeros in compositional data
#'
#' Replace zeros with small positive values using various imputation methods.
#' This is necessary for log-ratio transformations which require non-zero values.
#'
#' @param x A simp_composition or simp_cover object
#' @param method Imputation method:
#'   - `"multiplicative"`: Multiplicative replacement (default)
#'   - `"additive"`: Simple additive replacement
#'   - `"dl"`: Detection limit replacement
#' @param delta Replacement value. If NULL, automatically determined based on
#'   method. For multiplicative, uses 0.65 times minimum non-zero value.
#' @param preserve_structural Logical; if TRUE, don't impute zeros that appear
#'   to be structural (default FALSE)
#' @param structural_threshold Threshold for structural zeros (default 0.8)
#' @return Object with imputed values
#' @references
#' Martin-Fernandez, J.A., Barcelo-Vidal, C., Pawlowsky-Glahn, V. (2003)
#' Dealing with zeros and missing values in compositional data sets using
#' nonparametric imputation. Mathematical Geology, 35(3), 253-278.
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.4, 0.4, 0.2, 0,
#'            0.5, 0.3, 0.2, 0,
#'            0.3, 0.5, 0, 0.2), nrow = 3, byrow = TRUE),
#'   components = c("sp1", "sp2", "sp3", "sp4")
#' )
#'
#' # Multiplicative replacement
#' impute_zeros(x)
#'
#' # Additive with custom delta
#' impute_zeros(x, method = "additive", delta = 0.001)
impute_zeros <- function(x, method = c("multiplicative", "additive", "dl"),
                         delta = NULL, preserve_structural = FALSE,
                         structural_threshold = 0.8) {
  method <- match.arg(method)

  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition} or {.cls simp_cover}")
  }

  mat <- as.matrix(x)
  is_cover <- inherits(x, "simp_cover")

  components <- if (is_cover) attr(x, "species") else attr(x, "components")

  # Identify zeros to impute
  zero_analysis <- identify_zeros(x, structural_threshold = structural_threshold)

  if (zero_analysis$n_zeros == 0) {
    cli::cli_inform("No zeros to impute")
    return(x)
  }

  # Determine which zeros to replace
  zero_mask <- zero_analysis$zero_matrix

  if (preserve_structural && length(zero_analysis$structural_components) > 0) {
    structural_idx <- match(zero_analysis$structural_components, components)
    zero_mask[, structural_idx] <- FALSE
    n_impute <- sum(zero_mask)
    if (n_impute == 0) {
      cli::cli_inform("All zeros are structural; nothing to impute")
      return(x)
    }
    cli::cli_inform("Preserving structural zeros; imputing {n_impute} sampling zeros")
  }

  # Determine delta
  if (is.null(delta)) {
    non_zero <- mat[mat > 0]
    min_nz <- min(non_zero, na.rm = TRUE)
    delta <- switch(method,
      multiplicative = 0.65 * min_nz,
      additive = min_nz / 2,
      dl = min_nz
    )
  }

  # Apply imputation
  result <- switch(method,
    multiplicative = impute_multiplicative(mat, zero_mask, delta),
    additive = impute_additive(mat, zero_mask, delta),
    dl = impute_additive(mat, zero_mask, delta)  # Same as additive
  )

  colnames(result) <- components

  # Return appropriate type
  if (is_cover) {
    simp_cover(result,
               species = attr(x, "species"),
               layers = attr(x, "layers"),
               residual = attr(x, "residual"),
               validate = FALSE)
  } else {
    simp_composition(result,
                     components = attr(x, "components"),
                     validate = FALSE)
  }
}

#' Multiplicative replacement
#' @noRd
impute_multiplicative <- function(mat, zero_mask, delta) {
  result <- mat
  n <- nrow(mat)
  K <- ncol(mat)

  for (i in seq_len(n)) {
    zeros <- zero_mask[i, ]
    n_zeros <- sum(zeros)

    if (n_zeros > 0) {
      # Replacement value for this row
      row_sum <- sum(mat[i, ], na.rm = TRUE)
      replacement <- delta * row_sum

      # Replace zeros
      result[i, zeros] <- replacement

      # Adjust non-zeros to maintain sum
      nonzeros <- !zeros
      adjustment <- 1 - n_zeros * delta
      result[i, nonzeros] <- result[i, nonzeros] * adjustment
    }
  }

  result
}

#' Additive replacement
#' @noRd
impute_additive <- function(mat, zero_mask, delta) {
  result <- mat

  # Replace zeros
  result[zero_mask] <- delta

  # Renormalize rows
  row_sums <- rowSums(result, na.rm = TRUE)
  result <- result / row_sums

  result
}

#' Check for zeros that would cause issues
#'
#' Quick check to determine if data has zeros that need handling.
#'
#' @param x A simp_composition or simp_cover object
#' @return Logical indicating if zeros are present
#' @export
has_zeros <- function(x) {
  if (inherits(x, "simp_composition")) {
    attr(x, "has_zeros")
  } else {
    any(as_matrix(x) == 0, na.rm = TRUE)
  }
}

#' Count zeros in compositional data
#'
#' @param x A simp_composition or simp_cover object
#' @return Integer count of zeros
#' @export
count_zeros <- function(x) {
  mat <- if (inherits(x, "simp_composition")) as.matrix(x) else as_matrix(x)
  sum(mat == 0, na.rm = TRUE)
}
