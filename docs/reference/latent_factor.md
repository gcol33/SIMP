# Specify latent factors

Add latent factors to a SIMP model to account for unmeasured
confounders. Latent factors are observation-level random effects that
can capture shared variation across components not explained by the
fixed effects.

## Usage

``` r
latent_factor(n = 1, prior = NULL)
```

## Arguments

- n:

  Number of latent factors (default 1)

- prior:

  Prior for the factor loadings. Default is a PC prior that shrinks
  toward zero factors.

## Value

A `simp_latent` specification object

## References

Simpson et al. (2017) Penalising Model Component Complexity. Statistical
Science 32(1):1-28.

## Examples

``` r
# Single latent factor with default prior
latent_factor()

# Multiple factors with custom prior
latent_factor(n = 2, prior = prior_pc(u = 1, alpha = 0.01))

# Usage in model formula (v0.3.0+)
# simp(cover ~ temp + latent_factor(2) + (1 | site), data = dat)
```
