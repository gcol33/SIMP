# Add residual component to close composition

If compositional data doesn't sum to 1, add a residual component that
captures the remainder.

## Usage

``` r
add_residual(x, name = "other", position = c("last", "first"))
```

## Arguments

- x:

  A matrix, data frame, or simp_composition/simp_cover object

- name:

  Name for the residual component (default "other")

- position:

  Where to add the residual: "last" (default) or "first"

## Value

Object with added residual component

## Examples

``` r
# Open composition (sums to 0.8)
x <- matrix(c(0.3, 0.3, 0.2,
              0.4, 0.2, 0.2), nrow = 2, byrow = TRUE)
colnames(x) <- c("sp1", "sp2", "sp3")

# Close by adding residual
add_residual(x, name = "bare")
```
