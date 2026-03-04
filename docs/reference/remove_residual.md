# Remove residual component from composition

Remove a specified residual/other component, returning an open
composition.

## Usage

``` r
remove_residual(x, name = NULL)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- name:

  Name of the residual component to remove. If NULL, uses the residual
  attribute for simp_cover objects.

## Value

Object without the residual component

## Examples

``` r
x <- simp_cover(
  matrix(c(0.3, 0.5, 0.2,
           0.4, 0.4, 0.2), nrow = 2, byrow = TRUE),
  species = c("grass", "shrub", "bare"),
  residual = "bare"
)
remove_residual(x)
```
