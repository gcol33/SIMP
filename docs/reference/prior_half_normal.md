# Create a half-normal prior

Half-normal distribution for positive parameters. The distribution is a
normal distribution truncated at 0.

## Usage

``` r
prior_half_normal(sd = 1)
```

## Arguments

- sd:

  Prior scale parameter (default 1)

## Value

A `simp_prior` object

## Examples

``` r
prior_half_normal(2.5)
```
