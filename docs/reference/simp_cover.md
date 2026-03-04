# Create a cover data object

Constructor for the `simp_cover` class, a subclass of `simp_composition`
designed for vegetation cover data with optional layer (stratum)
support.

## Usage

``` r
simp_cover(
  x,
  species = NULL,
  layers = NULL,
  residual = NULL,
  validate = TRUE,
  normalize = FALSE
)
```

## Arguments

- x:

  A matrix or data frame of cover data

- species:

  Character vector of species/component names

- layers:

  Optional character vector specifying the vegetation layer for each
  column (e.g., "tree", "shrub", "herb")

- residual:

  Name of the residual component (e.g., "bare", "other"). If NULL and
  rows don't sum to 1, a residual will be suggested.

- validate:

  Logical; validate the composition (default TRUE)

- normalize:

  Logical; normalize rows to sum to 1 (default FALSE)

## Value

A `simp_cover` object

## Examples

``` r
# Simple cover data
x <- matrix(c(0.3, 0.2, 0.1, 0.4,
              0.4, 0.3, 0.2, 0.1), nrow = 2, byrow = TRUE)
cover <- simp_cover(x, species = c("oak", "pine", "grass", "bare"))

# With layers
cover <- simp_cover(x,
  species = c("oak", "pine", "grass", "bare"),
  layers = c("tree", "tree", "herb", "ground")
)
```
