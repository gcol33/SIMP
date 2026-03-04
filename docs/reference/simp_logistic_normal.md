# Logistic-normal family for compositional data

The logistic-normal distribution models compositional data via a
multivariate normal distribution on the log-ratio transformed space.
This allows flexible correlation structures between components.

## Usage

``` r
simp_logistic_normal(link = "identity", reference = "last")
```

## Arguments

- link:

  Link function for the latent space (default "identity")

- reference:

  Reference component for the additive log-ratio transform. Either
  "last" (default), "first", or a component name/index.

## Value

A `simp_family` object

## Examples

``` r
simp_logistic_normal()
simp_logistic_normal(reference = "first")
```
