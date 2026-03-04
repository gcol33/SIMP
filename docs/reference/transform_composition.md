# Convert between compositional representations

Transform compositional data between different representations.

## Usage

``` r
transform_composition(x, to = c("clr", "alr", "ilr"), reference = NULL)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- to:

  Target representation:

  - `"alr"`: Additive log-ratio (requires reference component)

  - `"clr"`: Centered log-ratio

  - `"ilr"`: Isometric log-ratio

- reference:

  For ALR, the index or name of the reference component

## Value

Transformed data as a matrix

## Examples

``` r
x <- simp_composition(
  matrix(c(0.3, 0.5, 0.2,
           0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
  components = c("a", "b", "c")
)

# Centered log-ratio
transform_composition(x, to = "clr")

# Additive log-ratio with "c" as reference
transform_composition(x, to = "alr", reference = "c")
```
