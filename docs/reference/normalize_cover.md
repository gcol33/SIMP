# Normalize cover data to sum to 1

Rescale rows of compositional data so they sum to 1.

## Usage

``` r
normalize_cover(x, method = c("proportional", "truncate"))
```

## Arguments

- x:

  A matrix, data frame, or simp_composition/simp_cover object

- method:

  Normalization method:

  - `"proportional"`: Scale all values proportionally (default)

  - `"truncate"`: Truncate values \> 1 before scaling

## Value

Object of the same class as input with normalized values

## Examples

``` r
# Data that doesn't sum to 1
x <- matrix(c(30, 50, 20,
              40, 40, 25), nrow = 2, byrow = TRUE)
normalize_cover(x)
```
