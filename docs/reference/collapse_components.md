# Collapse components according to a mapping

Aggregate multiple components into groups according to a named mapping.

## Usage

``` r
collapse_components(x, mapping, na.rm = TRUE)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- mapping:

  A named list where names are new component names and values are
  character vectors of original component names to aggregate

- na.rm:

  Remove NA values when summing (default TRUE)

## Value

Object with collapsed components

## Examples

``` r
x <- simp_composition(
  matrix(c(0.1, 0.2, 0.3, 0.2, 0.2,
           0.2, 0.1, 0.2, 0.3, 0.2), nrow = 2, byrow = TRUE),
  components = c("grass1", "grass2", "shrub1", "shrub2", "tree")
)

# Collapse grasses and shrubs
collapse_components(x, mapping = list(
  grass = c("grass1", "grass2"),
  shrub = c("shrub1", "shrub2"),
  tree = "tree"
))
```
