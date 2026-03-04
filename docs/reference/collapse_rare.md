# Collapse rare components into "other"

Identify rare components and collapse them into a single "other"
category.

## Usage

``` r
collapse_rare(x, threshold = 0.1, min_cover = 0.01, other_name = "other")
```

## Arguments

- x:

  A simp_composition or simp_cover object

- threshold:

  Threshold for "rare" (default 0.1)

- min_cover:

  Minimum mean cover (default 0.01)

- other_name:

  Name for the collapsed category (default "other")

## Value

Object with rare components collapsed

## Examples

``` r
x <- simp_composition(
  matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
           0.5, 0.3, 0.2, 0, 0,
           0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
  components = c("sp1", "sp2", "sp3", "sp4", "sp5")
)

collapse_rare(x, threshold = 0.5)
```
