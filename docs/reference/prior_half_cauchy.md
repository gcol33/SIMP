# Create a half-Cauchy prior

Half-Cauchy distribution for positive parameters. Recommended for
variance components due to heavy tails that allow for large values while
maintaining good behavior near zero.

## Usage

``` r
prior_half_cauchy(location = 0, scale = 1)
```

## Arguments

- location:

  Location parameter (default 0)

- scale:

  Scale parameter (default 1)

## Value

A `simp_prior` object

## Examples

``` r
prior_half_cauchy(0, 2)
```
