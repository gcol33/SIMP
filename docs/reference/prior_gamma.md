# Create a gamma prior

Gamma distribution for positive parameters.

## Usage

``` r
prior_gamma(shape = 1, rate = 1)
```

## Arguments

- shape:

  Shape parameter (default 1)

- rate:

  Rate parameter (default 1)

## Value

A `simp_prior` object

## Examples

``` r
prior_gamma(2, 0.1)
```
