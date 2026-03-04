# Stick-breaking family for ordered compositional data

The stick-breaking construction models compositional data as a sequence
of binary splits. Useful when components have a natural ordering.

## Usage

``` r
simp_stick_breaking(order = NULL, link = "logit")
```

## Arguments

- order:

  Component ordering. Either NULL (estimate or use data order), or a
  character vector of component names specifying the order.

- link:

  Link function for the break probabilities (default "logit")

## Value

A `simp_family` object

## Examples

``` r
simp_stick_breaking()

# Specified order
simp_stick_breaking(order = c("tree", "shrub", "herb", "bare"))
```
