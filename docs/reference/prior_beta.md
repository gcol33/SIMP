# Create a beta prior

Beta distribution for parameters bounded in \[0, 1\].

## Usage

``` r
prior_beta(alpha = 1, beta = 1)
```

## Arguments

- alpha:

  First shape parameter (default 1)

- beta:

  Second shape parameter (default 1)

## Value

A `simp_prior` object

## Examples

``` r
prior_beta(2, 2)  # Symmetric, peaked at 0.5
prior_beta(1, 1)  # Uniform
```
