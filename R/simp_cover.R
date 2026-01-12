# S3 class for vegetation cover data
# Subclass of simp_composition with layer support

#' Create a cover data object
#'
#' Constructor for the `simp_cover` class, a subclass of `simp_composition`
#' designed for vegetation cover data with optional layer (stratum) support.
#'
#' @param x A matrix or data frame of cover data
#' @param species Character vector of species/component names
#' @param layers Optional character vector specifying the vegetation layer
#'   for each column (e.g., "tree", "shrub", "herb")
#' @param residual Name of the residual component (e.g., "bare", "other").
#'   If NULL and rows don't sum to 1, a residual will be suggested.
#' @param validate Logical; validate the composition (default TRUE)
#' @param normalize Logical; normalize rows to sum to 1 (default FALSE)
#' @return A `simp_cover` object
#' @export
#' @examples
#' # Simple cover data
#' x <- matrix(c(0.3, 0.2, 0.1, 0.4,
#'               0.4, 0.3, 0.2, 0.1), nrow = 2, byrow = TRUE)
#' cover <- simp_cover(x, species = c("oak", "pine", "grass", "bare"))
#'
#' # With layers
#' cover <- simp_cover(x,
#'   species = c("oak", "pine", "grass", "bare"),
#'   layers = c("tree", "tree", "herb", "ground")
#' )
simp_cover <- function(x, species = NULL, layers = NULL, residual = NULL,
                       validate = TRUE, normalize = FALSE) {
  # Convert to matrix
  mat <- as_matrix(x)

  # Set species names
  if (is.null(species)) {
    if (!is.null(colnames(mat))) {
      species <- colnames(mat)
    } else {
      species <- paste0("sp", seq_len(ncol(mat)))
    }
  }

  check_character(species, "species")

  if (length(species) != ncol(mat)) {
    cli::cli_abort(c(
      "Length of {.arg species} must match number of columns",
      "x" = "Got {length(species)} names for {ncol(mat)} columns"
    ))
  }

  # Validate layers if provided
  if (!is.null(layers)) {
    check_character(layers, "layers")
    if (length(layers) != ncol(mat)) {
      cli::cli_abort(c(
        "Length of {.arg layers} must match number of columns",
        "x" = "Got {length(layers)} layers for {ncol(mat)} columns"
      ))
    }
  }

  colnames(mat) <- species

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
        "Invalid cover composition",
        "x" = paste(issues, collapse = ", "),
        "i" = "Use {.code validate = FALSE} to skip validation",
        "i" = "Use {.code normalize = TRUE} to auto-normalize row sums"
      ))
    }
  }

  # Compute zero info
  zero_matrix <- mat == 0
  has_zeros <- any(zero_matrix)

  # Identify residual component
  residual_idx <- NULL
  if (!is.null(residual)) {
    residual_idx <- match(residual, species)
    if (is.na(residual_idx)) {
      cli::cli_warn("Residual component {.val {residual}} not found in species")
      residual <- NULL
      residual_idx <- NULL
    }
  }

  # Compute layer info
  unique_layers <- if (!is.null(layers)) unique(layers) else NULL
  n_layers <- if (!is.null(layers)) length(unique_layers) else 1L

  new_simp_cover(
    mat,
    species = species,
    n_species = ncol(mat),
    n_obs = nrow(mat),
    has_zeros = has_zeros,
    zero_pattern = if (has_zeros) zero_matrix else NULL,
    layers = layers,
    unique_layers = unique_layers,
    n_layers = n_layers,
    residual = residual,
    residual_idx = residual_idx
  )
}

#' Low-level constructor for simp_cover
#' @noRd
new_simp_cover <- function(x, species, n_species, n_obs, has_zeros, zero_pattern,
                           layers, unique_layers, n_layers, residual,
                           residual_idx) {
  structure(
    x,
    components = species,  # For compatibility with simp_composition
    species = species,
    n_components = n_species,
    n_species = n_species,
    n_obs = n_obs,
    has_zeros = has_zeros,
    zero_pattern = zero_pattern,
    layers = layers,
    unique_layers = unique_layers,
    n_layers = n_layers,
    residual = residual,
    residual_idx = residual_idx,
    class = c("simp_cover", "simp_composition", "matrix", "array")
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.simp_cover <- function(x, ..., n = 6) {
  n_obs <- attr(x, "n_obs")
  n_species <- attr(x, "n_species")
  species <- attr(x, "species")
  layers <- attr(x, "unique_layers")

  cli::cli_h3("simp_cover")
  cat(sprintf("%d observations x %d species\n", n_obs, n_species))

  if (!is.null(layers)) {
    cat(sprintf("Layers: %s\n", paste(layers, collapse = ", ")))
  }

  if (!is.null(attr(x, "residual"))) {
    cat(sprintf("Residual: %s\n", attr(x, "residual")))
  }

  if (attr(x, "has_zeros")) {
    n_zeros <- sum(attr(x, "zero_pattern"))
    pct <- 100 * n_zeros / (n_obs * n_species)
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
summary.simp_cover <- function(object, ...) {
  n_obs <- attr(object, "n_obs")
  n_species <- attr(object, "n_species")
  species <- attr(object, "species")
  layers <- attr(object, "layers")

  # Species summaries
  species_summary <- data.frame(
    species = species,
    layer = if (!is.null(layers)) layers else NA_character_,
    min = apply(object, 2, min, na.rm = TRUE),
    median = apply(object, 2, median, na.rm = TRUE),
    mean = apply(object, 2, mean, na.rm = TRUE),
    max = apply(object, 2, max, na.rm = TRUE),
    zeros = colSums(object == 0, na.rm = TRUE),
    pct_zero = 100 * colSums(object == 0, na.rm = TRUE) / n_obs,
    row.names = NULL
  )

  # Layer summaries if applicable
  layer_summary <- NULL
  if (!is.null(layers)) {
    unique_layers <- attr(object, "unique_layers")
    layer_summary <- data.frame(
      layer = unique_layers,
      n_species = vapply(unique_layers, function(l) sum(layers == l), integer(1)),
      mean_cover = vapply(unique_layers, function(l) {
        mean(rowSums(object[, layers == l, drop = FALSE], na.rm = TRUE))
      }, numeric(1))
    )
  }

  structure(
    list(
      n_obs = n_obs,
      n_species = n_species,
      n_layers = attr(object, "n_layers"),
      species = species_summary,
      layers = layer_summary,
      has_zeros = attr(object, "has_zeros"),
      residual = attr(object, "residual")
    ),
    class = "summary.simp_cover"
  )
}

#' @export
print.summary.simp_cover <- function(x, ...) {
  cli::cli_h3("simp_cover Summary")
  cat(sprintf("Observations: %d\n", x$n_obs))
  cat(sprintf("Species: %d\n", x$n_species))

  if (x$n_layers > 1) {
    cat(sprintf("Layers: %d\n", x$n_layers))
  }

  if (!is.null(x$residual)) {
    cat(sprintf("Residual: %s\n", x$residual))
  }

  cat("\nSpecies summary:\n")
  print(x$species, row.names = FALSE)

  if (!is.null(x$layers)) {
    cat("\nLayer summary:\n")
    print(x$layers, row.names = FALSE)
  }

  invisible(x)
}

#' @export
`[.simp_cover` <- function(x, i, j, drop = FALSE) {
  species <- attr(x, "species")
  layers <- attr(x, "layers")
  zero_pattern <- attr(x, "zero_pattern")
  residual <- attr(x, "residual")
  residual_idx <- attr(x, "residual_idx")

  # Get underlying matrix
  mat <- unclass(x)

  # Subset
  if (missing(j)) {
    result <- mat[i, , drop = drop]
    new_species <- species
    new_layers <- layers
    if (!is.null(zero_pattern)) {
      zero_pattern <- zero_pattern[i, , drop = drop]
    }
  } else {
    result <- mat[i, j, drop = drop]
    new_species <- species[j]
    new_layers <- if (!is.null(layers)) layers[j] else NULL
    if (!is.null(zero_pattern)) {
      zero_pattern <- zero_pattern[i, j, drop = drop]
    }
    # Update residual index
    if (!is.null(residual)) {
      new_residual_idx <- match(residual, new_species)
      if (is.na(new_residual_idx)) {
        residual <- NULL
        residual_idx <- NULL
      } else {
        residual_idx <- new_residual_idx
      }
    }
  }

  # Return matrix if dropped to vector or single value
  if (!is.matrix(result)) {
    return(result)
  }

  new_simp_cover(
    result,
    species = new_species,
    n_species = ncol(result),
    n_obs = nrow(result),
    has_zeros = if (!is.null(zero_pattern)) any(zero_pattern) else FALSE,
    zero_pattern = zero_pattern,
    layers = new_layers,
    unique_layers = if (!is.null(new_layers)) unique(new_layers) else NULL,
    n_layers = if (!is.null(new_layers)) length(unique(new_layers)) else 1L,
    residual = residual,
    residual_idx = residual_idx
  )
}

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

#' Check if object is a simp_cover
#'
#' @param x Object to check
#' @return Logical
#' @export
is_simp_cover <- function(x) {
  inherits(x, "simp_cover")
}

#' Get species names
#'
#' @param x A simp_cover object
#' @return Character vector of species names
#' @export
species <- function(x) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover}")
  }
  attr(x, "species")
}

#' Get layer assignments
#'
#' @param x A simp_cover object
#' @return Character vector of layer names (or NULL if no layers)
#' @export
layers <- function(x) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover}")
  }
  attr(x, "layers")
}

#' Get number of species
#'
#' @param x A simp_cover object
#' @return Integer
#' @export
n_species <- function(x) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover}")
  }
  attr(x, "n_species")
}

#' Get number of layers
#'
#' @param x A simp_cover object
#' @return Integer
#' @export
n_layers <- function(x) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover}")
  }
  attr(x, "n_layers")
}
