# Special term functions for SIMP models
# varying() and shared() specify component-specific vs shared effects

#' Specify a varying effect
#'
#' Mark a predictor as having component-specific (varying) coefficients.
#' Each component of the compositional response will have its own
#' coefficient for this predictor.
#'
#' @param x A variable name or expression
#' @return A term marker for formula parsing
#' @export
#' @examples
#' # Each component gets its own temperature coefficient
#' # cover ~ varying(temp) + (1 | site)
#'
#' # Equivalent to the default behavior for fixed effects
#' # cover ~ temp + (1 | site)
varying <- function(x) {
  term <- deparse(substitute(x))
  structure(
    term,
    class = c("simp_varying_term", "simp_special_term", "character")
  )
}

#' Specify a shared effect
#'
#' Mark a predictor as having a shared coefficient across all components.
#' All components of the compositional response will use the same
#' coefficient for this predictor.
#'
#' @param x A variable name or expression
#' @return A term marker for formula parsing
#' @export
#' @examples
#' # All components share the same temperature coefficient
#' # cover ~ shared(temp) + (1 | site)
#'
#' # Mix shared and varying effects
#' # cover ~ shared(temp) + varying(precip) + (1 | site)
shared <- function(x) {

  term <- deparse(substitute(x))
  structure(
    term,
    class = c("simp_shared_term", "simp_special_term", "character")
  )
}

#' Check if term is a varying effect
#'
#' @param x Object to check
#' @return Logical
#' @noRd
is_varying <- function(x) {
  inherits(x, "simp_varying_term")
}

#' Check if term is a shared effect
#'
#' @param x Object to check
#' @return Logical
#' @noRd
is_shared <- function(x) {
  inherits(x, "simp_shared_term")
}

#' Check if term is a special SIMP term
#'
#' @param x Object to check
#' @return Logical
#' @noRd
is_special_term <- function(x) {
  inherits(x, "simp_special_term")
}
