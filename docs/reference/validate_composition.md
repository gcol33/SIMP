# Validate compositional data

Check if data satisfies compositional constraints and provide diagnostic
information about violations.

## Usage

``` r
validate_composition(x, tolerance = 1e-06, verbose = TRUE)
```

## Arguments

- x:

  A matrix, data frame, or simp_composition object

- tolerance:

  Tolerance for row sum check (default 1e-6)

- verbose:

  Print diagnostic output (default TRUE)

## Value

A list with validation results (invisibly)

## Examples

``` r
# Create some test data
x <- matrix(c(0.3, 0.5, 0.2,
              0.4, 0.4, 0.2,
              0.6, 0.3, 0.15), nrow = 3, byrow = TRUE)
validate_composition(x)
```
