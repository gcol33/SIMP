# Create a penalized complexity (PC) prior

Penalized complexity priors shrink toward a simple base model. The prior
is specified through P(sigma \> u) = alpha, where u is a threshold and
alpha is the tail probability.

## Usage

``` r
prior_pc(u, alpha)
```

## Arguments

- u:

  Upper bound threshold

- alpha:

  Tail probability P(sigma \> u)

## Value

A `simp_prior` object

## References

Simpson et al. (2017) Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors. Statistical
Science, 32(1), 1-28.

## Examples

``` r
# P(sigma > 1) = 0.01
prior_pc(u = 1, alpha = 0.01)
```
