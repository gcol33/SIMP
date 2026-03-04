# Remove random effects from formula

Remove random effects from formula

## Usage

``` r
remove_random_effects(formula)
```

## Arguments

- formula:

  A model formula

## Value

Formula with random effect terms removed

## Examples

``` r
remove_random_effects(y ~ x + (1 | group))
# Returns: y ~ x
```
