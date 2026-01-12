# Functions for collapsing compositional components
# Aggregate components by grouping or layer

#' Collapse components according to a mapping
#'
#' Aggregate multiple components into groups according to a named mapping.
#'
#' @param x A simp_composition or simp_cover object
#' @param mapping A named list where names are new component names and values
#'   are character vectors of original component names to aggregate
#' @param na.rm Remove NA values when summing (default TRUE)
#' @return Object with collapsed components
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.1, 0.2, 0.3, 0.2, 0.2,
#'            0.2, 0.1, 0.2, 0.3, 0.2), nrow = 2, byrow = TRUE),
#'   components = c("grass1", "grass2", "shrub1", "shrub2", "tree")
#' )
#'
#' # Collapse grasses and shrubs
#' collapse_components(x, mapping = list(
#'   grass = c("grass1", "grass2"),
#'   shrub = c("shrub1", "shrub2"),
#'   tree = "tree"
#' ))
collapse_components <- function(x, mapping, na.rm = TRUE) {
  if (!inherits(x, "simp_composition")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_composition} or {.cls simp_cover}")
  }

  mat <- as.matrix(x)
  is_cover <- inherits(x, "simp_cover")

  components <- if (is_cover) attr(x, "species") else attr(x, "components")

  # Validate mapping
  if (!is.list(mapping) || is.null(names(mapping))) {
    cli::cli_abort("{.arg mapping} must be a named list")
  }

  all_mapped <- unlist(mapping)
  missing <- setdiff(all_mapped, components)
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Components in mapping not found in data",
      "x" = "Missing: {.val {missing}}"
    ))
  }

  # Check for duplicates in mapping
  if (anyDuplicated(all_mapped)) {
    dups <- all_mapped[duplicated(all_mapped)]
    cli::cli_abort(c(
      "Components appear in multiple groups",
      "x" = "Duplicated: {.val {unique(dups)}}"
    ))
  }

  # Warn about unmapped components
  unmapped <- setdiff(components, all_mapped)
  if (length(unmapped) > 0) {
    cli::cli_warn(c(
      "Some components not in mapping will be dropped",
      "i" = "Unmapped: {.val {unmapped}}"
    ))
  }

  # Collapse
  new_names <- names(mapping)
  n_new <- length(new_names)
  result <- matrix(0, nrow = nrow(mat), ncol = n_new)
  colnames(result) <- new_names

  for (i in seq_along(mapping)) {
    cols <- mapping[[i]]
    col_idx <- match(cols, components)
    result[, i] <- rowSums(mat[, col_idx, drop = FALSE], na.rm = na.rm)
  }

  # Handle layers for simp_cover
  if (is_cover) {
    old_layers <- attr(x, "layers")
    new_layers <- if (!is.null(old_layers)) {
      # Use layer of first component in each group
      vapply(mapping, function(cols) {
        idx <- match(cols[1], components)
        old_layers[idx]
      }, character(1))
    } else {
      NULL
    }

    simp_cover(result,
               species = new_names,
               layers = new_layers,
               validate = FALSE)
  } else {
    simp_composition(result,
                     components = new_names,
                     validate = FALSE)
  }
}

#' Collapse vegetation layers
#'
#' Aggregate cover data by vegetation layer, summing species within each layer.
#'
#' @param x A simp_cover object with layer information
#' @param fun Aggregation function (default `sum`). Can also use `max` for
#'   situations where overlapping cover exists.
#' @param na.rm Remove NA values when aggregating (default TRUE)
#' @return A simp_composition with layers as components
#' @export
#' @examples
#' x <- simp_cover(
#'   matrix(c(0.1, 0.2, 0.3, 0.4,
#'            0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
#'   species = c("oak", "pine", "grass", "bare"),
#'   layers = c("tree", "tree", "herb", "ground")
#' )
#'
#' # Sum within layers
#' collapse_layers(x)
collapse_layers <- function(x, fun = sum, na.rm = TRUE) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover} object")
  }

  layers <- attr(x, "layers")
  if (is.null(layers)) {
    cli::cli_abort(c(
      "No layer information in {.cls simp_cover} object",
      "i" = "Specify layers when creating the cover object"
    ))
  }

  mat <- as.matrix(x)
  unique_layers <- attr(x, "unique_layers")
  n_layers <- length(unique_layers)

  # Aggregate by layer
  result <- matrix(0, nrow = nrow(mat), ncol = n_layers)
  colnames(result) <- unique_layers

  for (i in seq_along(unique_layers)) {
    layer <- unique_layers[i]
    col_idx <- which(layers == layer)
    if (length(col_idx) == 1) {
      result[, i] <- mat[, col_idx]
    } else {
      result[, i] <- apply(mat[, col_idx, drop = FALSE], 1, fun, na.rm = na.rm)
    }
  }

  simp_composition(result,
                   components = unique_layers,
                   validate = FALSE)
}

#' Split composition by layer
#'
#' Extract species within a specific vegetation layer.
#'
#' @param x A simp_cover object with layer information
#' @param layer Name of the layer to extract
#' @param normalize Normalize the extracted layer to sum to 1 (default FALSE)
#' @return A simp_cover object with only species from the specified layer
#' @export
#' @examples
#' x <- simp_cover(
#'   matrix(c(0.1, 0.2, 0.3, 0.4,
#'            0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
#'   species = c("oak", "pine", "grass", "bare"),
#'   layers = c("tree", "tree", "herb", "ground")
#' )
#'
#' # Extract tree layer only
#' split_layer(x, "tree")
split_layer <- function(x, layer, normalize = FALSE) {
  if (!inherits(x, "simp_cover")) {
    cli::cli_abort("{.arg x} must be a {.cls simp_cover} object")
  }

  layers <- attr(x, "layers")
  if (is.null(layers)) {
    cli::cli_abort(c(
      "No layer information in {.cls simp_cover} object",
      "i" = "Specify layers when creating the cover object"
    ))
  }

  if (!layer %in% layers) {
    cli::cli_abort(c(
      "Layer not found",
      "x" = "{.val {layer}} is not in layers: {.val {unique(layers)}}"
    ))
  }

  # Extract columns for this layer
  col_idx <- which(layers == layer)
  result <- x[, col_idx, drop = FALSE]

  if (normalize) {
    result <- normalize_cover(result)
  }

  result
}

#' Identify rare components
#'
#' Find components that occur rarely across observations.
#'
#' @param x A simp_composition or simp_cover object
#' @param threshold Threshold for "rare" as proportion of observations
#'   with non-zero values (default 0.1, i.e., present in < 10% of obs)
#' @param min_cover Minimum mean cover to not be considered rare (default 0.01)
#' @return Character vector of rare component names
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
#'            0.5, 0.3, 0.2, 0, 0,
#'            0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
#'   components = c("sp1", "sp2", "sp3", "sp4", "sp5")
#' )
#'
#' # Find components present in < 50% of observations
#' identify_rare(x, threshold = 0.5)
identify_rare <- function(x, threshold = 0.1, min_cover = 0.01) {
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

  # Frequency of occurrence (proportion of obs with non-zero)
  freq <- colSums(mat > 0, na.rm = TRUE) / n_obs

  # Mean cover
  mean_cover <- colMeans(mat, na.rm = TRUE)

  # Rare: low frequency AND low mean cover
  rare <- (freq < threshold) & (mean_cover < min_cover)

  components[rare]
}

#' Collapse rare components into "other"
#'
#' Identify rare components and collapse them into a single "other" category.
#'
#' @param x A simp_composition or simp_cover object
#' @param threshold Threshold for "rare" (default 0.1)
#' @param min_cover Minimum mean cover (default 0.01)
#' @param other_name Name for the collapsed category (default "other")
#' @return Object with rare components collapsed
#' @export
#' @examples
#' x <- simp_composition(
#'   matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
#'            0.5, 0.3, 0.2, 0, 0,
#'            0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
#'   components = c("sp1", "sp2", "sp3", "sp4", "sp5")
#' )
#'
#' collapse_rare(x, threshold = 0.5)
collapse_rare <- function(x, threshold = 0.1, min_cover = 0.01,
                          other_name = "other") {

  rare <- identify_rare(x, threshold = threshold, min_cover = min_cover)

  if (length(rare) == 0) {
    cli::cli_inform("No rare components found with current thresholds")
    return(x)
  }

  cli::cli_inform("Collapsing {length(rare)} rare component(s): {.val {rare}}")

  # Build mapping
  components <- if (inherits(x, "simp_cover")) {
    attr(x, "species")
  } else {
    attr(x, "components")
  }

  keep <- setdiff(components, rare)

  mapping <- as.list(keep)
  names(mapping) <- keep
  mapping[[other_name]] <- rare

  collapse_components(x, mapping)
}
