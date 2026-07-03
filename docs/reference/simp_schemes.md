# Available symplectic integrators

The built-in scheme family: leapfrog and the Yoshida members generated
from it by triple-jump composition.

## Usage

``` r
simp_schemes()
```

## Value

A data frame with one row per integrator giving its name, order of
accuracy, and number of gradient evaluations per step.

## Examples

``` r
simp_schemes()
#>   integrator order gradient_evals
#> 1   leapfrog     2              1
#> 2  minerror2     2              2
#> 3   yoshida4     4              3
#> 4   yoshida6     6              9
#> 5   yoshida8     8             27
```
