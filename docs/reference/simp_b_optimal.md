# Minimum-error multistage coefficients for a step band

The multistage kick and drift weights that minimize the worst-case
energy error a splitting integrator induces on a harmonic (Gaussian)
target across the band of dimensionless steps `(0, nu_max]`. A separable
quadratic target excites oscillator modes at `nu_k = omega_k * eps`, so
with `nu_max = omega_max * eps` (the largest such step, from the highest
target frequency) these coefficients give the minimum-error integrator
for the problem at hand rather than a fixed compromise. Feed the
coefficients back through
[`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md)
with `integrator = "adaptive2"` or `"adaptive3"`, which recompute them
internally from the same `nu_max`.

## Usage

``` r
simp_b_optimal(nu_max, stages = 2)
```

## Arguments

- nu_max:

  Largest dimensionless step `omega_max * eps`, a single positive
  number.

- stages:

  `2` for the two-stage scheme (returns the kick weight `b`) or `3` for
  the three-stage scheme (returns `c(b, a)`, the kick and drift free
  weights).

## Value

A numeric vector of the free coefficients: length 1 for two stages,
length 2 for three.

## Examples

``` r
simp_b_optimal(nu_max = 1.5, stages = 2)
```
