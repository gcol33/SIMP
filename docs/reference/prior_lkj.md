# Create a LKJ prior for correlation matrices

LKJ distribution for correlation matrices. The parameter eta controls
the concentration around the identity matrix.

## Usage

``` r
prior_lkj(eta = 1)
```

## Arguments

- eta:

  Concentration parameter. eta = 1 gives uniform distribution over
  correlation matrices. eta \> 1 shrinks toward identity matrix.

## Value

A `simp_prior` object

## Examples

``` r
prior_lkj(2)  # Shrink toward uncorrelated
```
