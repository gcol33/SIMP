# Collapse vegetation layers

Aggregate cover data by vegetation layer, summing species within each
layer.

## Usage

``` r
collapse_layers(x, fun = sum, na.rm = TRUE)
```

## Arguments

- x:

  A simp_cover object with layer information

- fun:

  Aggregation function (default `sum`). Can also use `max` for
  situations where overlapping cover exists.

- na.rm:

  Remove NA values when aggregating (default TRUE)

## Value

A simp_composition with layers as components

## Examples

``` r
x <- simp_cover(
  matrix(c(0.1, 0.2, 0.3, 0.4,
           0.2, 0.1, 0.4, 0.3), nrow = 2, byrow = TRUE),
  species = c("oak", "pine", "grass", "bare"),
  layers = c("tree", "tree", "herb", "ground")
)

# Sum within layers
collapse_layers(x)
```
