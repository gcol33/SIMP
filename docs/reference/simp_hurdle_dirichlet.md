# Hurdle Dirichlet family for zero-inflated compositional data

A two-part model for compositional data with structural zeros. The first
part models the probability of each component being zero; the second
part models the non-zero composition using a Dirichlet distribution.

## Usage

``` r
simp_hurdle_dirichlet(
  zero_formula = ~1,
  link_zero = "logit",
  link_positive = "log"
)
```

## Arguments

- zero_formula:

  Formula for the zero model (default `~ 1`)

- link_zero:

  Link function for zero probability (default "logit")

- link_positive:

  Link function for positive concentration (default "log")

## Value

A `simp_family` object

## Examples

``` r
# Intercept-only zero model
simp_hurdle_dirichlet()

# Covariate-dependent zeros
simp_hurdle_dirichlet(zero_formula = ~ elevation)
```
