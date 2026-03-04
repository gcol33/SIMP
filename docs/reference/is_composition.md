# Check if data is a valid composition

Check if data is a valid composition

## Usage

``` r
is_composition(x, tolerance = 1e-06)
```

## Arguments

- x:

  A matrix, data frame, or simp_composition object

- tolerance:

  Tolerance for row sum check (default 1e-6)

## Value

Logical indicating if data is valid

## Examples

``` r
x <- matrix(c(0.3, 0.5, 0.2,
              0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)
is_composition(x)
```
