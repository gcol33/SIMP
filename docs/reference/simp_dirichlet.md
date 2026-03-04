# Dirichlet family for compositional data

The Dirichlet distribution is the natural choice for compositional data,
producing vectors that sum to 1 with all positive components.

## Usage

``` r
simp_dirichlet(link = "log", concentration = NULL)
```

## Arguments

- link:

  Link function for the concentration parameter (default "log")

- concentration:

  Fixed concentration parameter, or NULL to estimate

## Value

A `simp_family` object

## Examples

``` r
# Default: estimate concentration
simp_dirichlet()

# Fixed concentration
simp_dirichlet(concentration = 10)
```
