# Validation functions for compositional data
# Diagnostic tools for data quality assessment

#' Validate compositional data
#'
#' Check if data satisfies compositional constraints and provide diagnostic
#' information about violations.
#'
#' @param x A matrix, data frame, or simp_composition object
#' @param tolerance Tolerance for row sum check (default 1e-6)
#' @param verbose Print diagnostic output (default TRUE)
#' @return A list with validation results (invisibly)
#' @export
#' @examples
#' # Create some test data
#' x <- matrix(c(0.3, 0.5, 0.2,
#'               0.4, 0.4, 0.2,
#'               0.6, 0.3, 0.15), nrow = 3, byrow = TRUE)
#' validate_composition(x)
validate_composition <- function(x, tolerance = 1e-6, verbose = TRUE) {
  # Convert to matrix
  if (inherits(x, "simp_composition")) {
    mat <- as.matrix(x)
    species <- attr(x, "species")
  } else {
    mat <- as_matrix(x)
    species <- colnames(mat)
  }

  if (is.null(species)) {
    species <- paste0("V", seq_len(ncol(mat)))
  }

  n_obs <- nrow(mat)
  n_comp <- ncol(mat)

  # Compute diagnostics
  row_sums <- rowSums(mat, na.rm = TRUE)
  has_negatives <- any(mat < 0, na.rm = TRUE)
  has_above_one <- any(mat > 1, na.rm = TRUE)
  has_na <- any(is.na(mat))

  # Zero analysis
  n_zeros <- sum(mat == 0, na.rm = TRUE)
  pct_zeros <- 100 * n_zeros / (n_obs * n_comp)
  zeros_per_component <- colSums(mat == 0, na.rm = TRUE)
  high_zero_components <- names(zeros_per_component)[zeros_per_component > 0.5 * n_obs]

  # Row sum analysis
  row_sum_issues <- abs(row_sums - 1) > tolerance
  n_row_issues <- sum(row_sum_issues)

  # Determine validity
  is_valid <- !has_negatives && !has_above_one && n_row_issues == 0 && !has_na

  # Build result
  result <- list(
    valid = is_valid,
    n_obs = n_obs,
    n_components = n_comp,
    row_sums = list(
      min = min(row_sums),
      max = max(row_sums),
      mean = mean(row_sums),
      n_issues = n_row_issues,
      pct_issues = 100 * n_row_issues / n_obs
    ),
    zeros = list(
      total = n_zeros,
      pct = pct_zeros,
      per_component = zeros_per_component,
      high_zero_components = high_zero_components
    ),
    issues = list(
      has_negatives = has_negatives,
      has_above_one = has_above_one,
      has_na = has_na,
      row_sum_issues = n_row_issues
    )
  )

  # Generate recommendations
  result$recommendations <- generate_recommendations(result)

  class(result) <- "simp_validation"

  if (verbose) {
    print(result)
  }

  invisible(result)
}

#' Generate recommendations based on validation results
#' @noRd
generate_recommendations <- function(result) {
  recs <- character(0)

  if (result$issues$has_negatives) {
    recs <- c(recs, "Negative values found. Check data entry or measurement errors.")
  }

  if (result$issues$has_above_one) {
    recs <- c(recs, "Values > 1 found. Ensure data is in proportion scale [0, 1].")
  }

  if (result$issues$has_na) {
    recs <- c(recs, "Missing values found. Consider imputation or removal.")
  }

  if (result$row_sums$n_issues > 0) {
    if (result$row_sums$max > 1.1) {
      recs <- c(recs, "Use normalize_cover() to rescale rows to sum to 1.")
    } else if (result$row_sums$max < 0.99) {
      recs <- c(recs, "Use add_residual() to close the composition.")
    }
  }

  if (length(result$zeros$high_zero_components) > 0) {
    recs <- c(recs, sprintf(
      "Components with >50%% zeros: %s. Consider collapsing sparse components.",
      paste(result$zeros$high_zero_components, collapse = ", ")
    ))
  }

  if (result$zeros$pct > 50) {
    recs <- c(recs, "High zero proportion. Consider a hurdle model (simp_hurdle_dirichlet).")
  }

  recs
}

#' @export
print.simp_validation <- function(x, ...) {
  cli::cli_h3("Composition Diagnostics")

  cat(sprintf("Observations: %d\n", x$n_obs))
  cat(sprintf("Components: %d\n", x$n_components))

  cat("\nRow sums:\n")
  cat(sprintf("  Min: %.3f, Max: %.3f, Mean: %.3f\n",
              x$row_sums$min, x$row_sums$max, x$row_sums$mean))
  if (x$row_sums$n_issues > 0) {
    cat(sprintf("  Rows not summing to 1: %d (%.1f%%)\n",
                x$row_sums$n_issues, x$row_sums$pct_issues))
  }

  cat("\nZeros:\n")
  cat(sprintf("  Total zeros: %d (%.1f%%)\n",
              x$zeros$total, x$zeros$pct))
  if (length(x$zeros$high_zero_components) > 0) {
    cat(sprintf("  Components with >50%% zeros: %s\n",
                paste(x$zeros$high_zero_components, collapse = ", ")))
  }

  if (length(x$issues) > 0) {
    issues <- c()
    if (x$issues$has_negatives) issues <- c(issues, "negative values")
    if (x$issues$has_above_one) issues <- c(issues, "values > 1")
    if (x$issues$has_na) issues <- c(issues, "missing values")
    if (length(issues) > 0) {
      cat("\nIssues found:", paste(issues, collapse = ", "), "\n")
    }
  }

  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (rec in x$recommendations) {
      cat("  -", rec, "\n")
    }
  }

  if (x$valid) {
    cli::cli_alert_success("Data is a valid composition")
  } else {
    cli::cli_alert_warning("Data is NOT a valid composition")
  }

  invisible(x)
}

#' Check if data is a valid composition
#'
#' @param x A matrix, data frame, or simp_composition object
#' @param tolerance Tolerance for row sum check (default 1e-6)
#' @return Logical indicating if data is valid
#' @export
#' @examples
#' x <- matrix(c(0.3, 0.5, 0.2,
#'               0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)
#' is_composition(x)
is_composition <- function(x, tolerance = 1e-6) {
  result <- validate_composition(x, tolerance = tolerance, verbose = FALSE)
  result$valid
}

#' Validate that data is suitable for SIMP modeling
#'
#' Performs comprehensive checks beyond basic composition validity.
#'
#' @param x Data to validate
#' @param min_obs Minimum observations required (default 10)
#' @param min_components Minimum components required (default 2)
#' @param max_zero_pct Maximum allowed zero percentage (default 95)
#' @return TRUE invisibly if valid, error otherwise
#' @export
validate_simp_data <- function(x, min_obs = 10, min_components = 2,
                               max_zero_pct = 95) {
  # Convert to matrix
  mat <- if (inherits(x, "simp_composition")) {
    as.matrix(x)
  } else {
    as_matrix(x)
  }

  n_obs <- nrow(mat)
  n_comp <- ncol(mat)

  # Check dimensions

if (n_obs < min_obs) {
    cli::cli_abort(c(
      "Insufficient observations",
      "x" = "Found {n_obs}, need at least {min_obs}"
    ))
  }

  if (n_comp < min_components) {
    cli::cli_abort(c(
      "Insufficient components",
      "x" = "Found {n_comp}, need at least {min_components}"
    ))
  }

  # Check for all-zero rows or columns
  row_all_zero <- rowSums(mat, na.rm = TRUE) == 0
  if (any(row_all_zero)) {
    cli::cli_abort(c(
      "Found rows with all zeros",
      "x" = "{sum(row_all_zero)} row(s) sum to zero"
    ))
  }

  col_all_zero <- colSums(mat, na.rm = TRUE) == 0
  if (any(col_all_zero)) {
    cli::cli_abort(c(
      "Found components with all zeros",
      "x" = "Components: {.val {colnames(mat)[col_all_zero]}}"
    ))
  }

  # Check zero percentage
  zero_pct <- 100 * sum(mat == 0, na.rm = TRUE) / length(mat)
  if (zero_pct > max_zero_pct) {
    cli::cli_abort(c(
      "Too many zeros in data",
      "x" = "{format(zero_pct, digits = 1)}% zeros exceeds {max_zero_pct}% threshold",
      "i" = "Consider collapsing sparse components"
    ))
  }

  invisible(TRUE)
}
