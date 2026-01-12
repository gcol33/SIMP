# Internal utilities for SIMP
# These functions are not exported

#' @importFrom stats as.formula median quantile
NULL

#' Check if input is a valid matrix-like object
#' @noRd
is_matrix_like <- function(x) {

inherits(x, c("matrix", "data.frame", "Matrix"))
}

#' Convert to matrix
#' @noRd
as_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) {
    # Only convert numeric columns
    num_cols <- vapply(x, is.numeric, logical(1))
    if (!all(num_cols)) {
      cli::cli_abort(c(
        "All columns must be numeric",
        "x" = "Non-numeric columns: {.val {names(x)[!num_cols]}}"
      ))
    }
    return(as.matrix(x))
  }
  if (inherits(x, "Matrix")) {
    return(as.matrix(x))
  }
  cli::cli_abort("Cannot convert {.cls {class(x)}} to matrix")
}

#' Validate that x is positive
#' @noRd
check_positive <- function(x, name = "x") {
  if (any(x <= 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "{.arg {name}} must be positive",
      "x" = "Found {sum(x <= 0, na.rm = TRUE)} non-positive value(s)"
    ))
  }
  invisible(x)
}

#' Validate that x is non-negative
#' @noRd
check_non_negative <- function(x, name = "x") {
  if (any(x < 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "{.arg {name}} must be non-negative",
      "x" = "Found {sum(x < 0, na.rm = TRUE)} negative value(s)"
    ))
  }
  invisible(x)
}

#' Validate that x is in 0 to 1
#' @noRd
check_probability <- function(x, name = "x") {
  if (any(x < 0 | x > 1, na.rm = TRUE)) {
    cli::cli_abort(c(
      "{.arg {name}} must be in [0, 1]",
      "x" = "Found {sum(x < 0 | x > 1, na.rm = TRUE)} value(s) outside range"
    ))
  }
  invisible(x)
}

#' Validate that x is a scalar
#' @noRd
check_scalar <- function(x, name = "x") {
  if (length(x) != 1) {
    cli::cli_abort(c(
      "{.arg {name}} must be a scalar",
      "x" = "Got length {length(x)}"
    ))
  }
  invisible(x)
}

#' Validate that x is a character vector
#' @noRd
check_character <- function(x, name = "x") {
  if (!is.character(x)) {
    cli::cli_abort(c(
      "{.arg {name}} must be a character vector",
      "x" = "Got {.cls {class(x)}}"
    ))
  }
  invisible(x)
}

#' Validate that x is a logical scalar
#' @noRd
check_flag <- function(x, name = "x") {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    cli::cli_abort(c(
      "{.arg {name}} must be TRUE or FALSE",
      "x" = "Got {.val {x}}"
    ))
  }
  invisible(x)
}

#' Format a numeric vector for display
#' @noRd
format_num <- function(x, digits = 3) {
  format(round(x, digits), nsmall = digits)
}

#' Pluralize a word based on count
#' @noRd
plural <- function(n, singular, plural = paste0(singular, "s")) {
  if (n == 1) singular else plural
}
