# Identify rare components

Find components that occur rarely across observations.

## Usage

``` r
identify_rare(x, threshold = 0.1, min_cover = 0.01)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- threshold:

  Threshold for "rare" as proportion of observations with non-zero
  values (default 0.1, i.e., present in \< 10% of obs)

- min_cover:

  Minimum mean cover to not be considered rare (default 0.01)

## Value

Character vector of rare component names

## Examples

``` r
x <- simp_composition(
  matrix(c(0.4, 0.4, 0.1, 0.05, 0.05,
           0.5, 0.3, 0.2, 0, 0,
           0.3, 0.5, 0.2, 0, 0), nrow = 3, byrow = TRUE),
  components = c("sp1", "sp2", "sp3", "sp4", "sp5")
)

# Find components present in < 50% of observations
identify_rare(x, threshold = 0.5)
```
