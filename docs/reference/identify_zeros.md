# Identify zeros in compositional data

Analyze zero patterns in compositional data and classify them as
potentially structural (true absence) or sampling (rounded/below
detection).

## Usage

``` r
identify_zeros(x, structural_threshold = 0.8)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- structural_threshold:

  Components with zeros in more than this proportion of observations are
  likely structural (default 0.8)

## Value

A list with zero pattern information

## Examples

``` r
x <- simp_composition(
  matrix(c(0.4, 0.4, 0.2, 0,
           0.5, 0.3, 0.2, 0,
           0.3, 0.5, 0, 0.2), nrow = 3, byrow = TRUE),
  components = c("sp1", "sp2", "sp3", "sp4")
)

identify_zeros(x)
```
