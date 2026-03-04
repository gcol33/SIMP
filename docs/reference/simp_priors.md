# Specify priors for a SIMP model

Create a container of priors for different model parameters.

## Usage

``` r
simp_priors(
  beta = prior_normal(0, 2.5),
  intercept = prior_normal(0, 5),
  sigma = prior_half_cauchy(0, 2),
  phi = prior_gamma(2, 0.1),
  cor = prior_lkj(2),
  ...
)
```

## Arguments

- beta:

  Prior for fixed effect coefficients (default: normal(0, 2.5))

- intercept:

  Prior for intercept (default: normal(0, 5))

- sigma:

  Prior for random effect standard deviations (default: half_cauchy(0,
  2))

- phi:

  Prior for Dirichlet concentration parameter (default: gamma(2, 0.1))

- cor:

  Prior for correlation matrices in random effects (default: lkj(2))

- ...:

  Additional named priors for custom parameters

## Value

A `simp_priors` object

## Examples

``` r
# Use defaults
simp_priors()

# Custom priors
simp_priors(
  beta = prior_normal(0, 10),
  sigma = prior_half_normal(2)
)
```
