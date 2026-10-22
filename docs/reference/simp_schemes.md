# Available symplectic integrators

The built-in scheme family: leapfrog, the two-stage minimum-error
member, and the Yoshida members generated from leapfrog by triple-jump
composition. The step-adapted `"adaptive2"` / `"adaptive3"` schemes are
not listed here because their coefficients depend on the step band; see
[`simp_b_optimal()`](https://gcol33.github.io/SIMP/reference/simp_b_optimal.md).

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
```
