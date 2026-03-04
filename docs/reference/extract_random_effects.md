# Extract random effect terms from formula

Parses random effect specifications including:

- `(1 | group)` - random intercepts

- `(1 + x | group)` - correlated random intercept and slope

- `(0 + x | group)` - random slope only (no intercept)

- `(1 + x || group)` - uncorrelated random effects

- `(1 | g1/g2)` - nested random effects

## Usage

``` r
extract_random_effects(formula)
```

## Arguments

- formula:

  A model formula

## Value

A list of random effect specifications

## Examples

``` r
extract_random_effects(y ~ x + (1 | group))
extract_random_effects(y ~ x + (1 + x | group))
extract_random_effects(y ~ x + (1 | site/plot))
```
