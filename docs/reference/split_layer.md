# Split composition by layer

Extract species within a specific vegetation layer.

## Usage

``` r
split_layer(x, layer, normalize = FALSE)
```

## Arguments

- x:

  A simp_cover object with layer information

- layer:

  Name of the layer to extract

- normalize:

  Normalize the extracted layer to sum to 1 (default FALSE)

## Value

A simp_cover object with only species from the specified layer

## Examples

``` r
x <- simp_cover(
  matrix(c(0.1, 0.2, 0.3, 0.4,
           0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
  species = c("oak", "pine", "grass", "bare"),
  layers = c("tree", "tree", "herb", "ground")
)

# Extract tree layer only
split_layer(x, "tree")
```
