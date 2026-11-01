
# SIMP

*symplectic integrators for Hamiltonian Monte Carlo*

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/gcol33/SIMP/graph/badge.svg)](https://app.codecov.io/gh/gcol33/SIMP)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

SIMP (Symplectic Integrators for Monte-Carlo Proposals) integrates
Hamiltonian systems and builds the Hamiltonian Monte Carlo proposals
that run on them.

The integrators are splitting schemes: an alternating sequence of
position and momentum updates whose coefficients define the method. A
symplectic scheme conserves a nearby energy, so a trajectory keeps its
total energy within a bounded band over long horizons and the proposals
it produces are accepted at high rates. Leapfrog and the Yoshida order
4, 6 and 8 members come from a single triple-jump composition, so the
whole family shares one stepper and one coefficient source, written as
header-only C++. You can compare schemes on your target by their energy
error and their gradient cost per step, and choose the one that drives
your sampler.

SIMP supports:

- **A family of schemes.** Leapfrog, the two-stage minimum-error member
  `minerror2` and Yoshida orders 4, 6 and 8.
  [`simp_schemes()`](https://gcol33.github.io/SIMP/reference/simp_schemes.html)
  lists each with its order and gradient evaluations per step.
- **Coefficients tuned to the step band.** The `"adaptive2"` and
  `"adaptive3"` schemes take the two- or three-stage weights that
  minimise the worst-case energy error on a Gaussian target over a band
  of dimensionless steps, computed from an exact harmonic analysis.
  [`simp_b_optimal()`](https://gcol33.github.io/SIMP/reference/simp_b_optimal.html)
  returns the weights for a given band.
- **Hamiltonian Monte Carlo.** `simp_hmc()` draws each proposal by
  integrating the target’s Hamiltonian flow and accepts it by the
  Metropolis rule on the energy, with any scheme as the integrator.
- **Split potentials.** `simp_integrate_mts()` advances a potential
  split into a slow and a fast part at two rates (RESPA), so the fast
  part takes small inner steps without evaluating the slow part at the
  same rate.
- **Non-separable Hamiltonians.** `simp_integrate_nonsep()` uses the
  generalized (implicit) leapfrog, as in Riemannian dynamics with a
  position-dependent mass.
- **A C++ core for other packages.** The steppers are scalar-templated
  headers in `inst/include`. A package reaches them with
  `LinkingTo: SIMP` and drives them with its own gradient, under
  automatic differentiation as well as in double precision.

## Installation

SIMP compiles C++ code. Install the build tools for your operating
system first:

| Operating system | Build tools |
|----|----|
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version |
| macOS | Run `xcode-select --install` in Terminal |
| Debian / Ubuntu | Install `r-base-dev` |

Check that R can find the compiler:

``` r
pkgbuild::has_build_tools(debug = TRUE)
```

Then install SIMP from GitHub:

``` r
# install.packages("pak")
pak::pak("gcol33/SIMP")
```

The C++ code uses the Eigen headers that RcppEigen supplies and needs no
system library.

## Usage

1.  Write the gradient of the potential energy `U(q)` as an R function
    of the position.
2.  Call `simp_integrate()` with a starting position and momentum, the
    gradient, the step size and the number of steps.
3.  Sample a target with `simp_hmc()`, which runs the same integrator
    inside each proposal.

The potential of a unit Gaussian target is `U(q) = 0.5 q'q`, so its
gradient is `q`. Starting on the unit circle in phase space, the exact
flow is a rotation:

``` r
library(SIMP)

grad_U <- function(q) q
out <- simp_integrate(q0 = c(1, 0), p0 = c(0, 1), grad_U = grad_U,
                      eps = 0.1, n_steps = 60)
out$q
#> [1]  0.9608666 -0.2773584
```

The result is a list with the endpoint position `q` and momentum `p`.
The default integrator is leapfrog, and `integrator` takes any name that
`simp_schemes()` lists:

``` r
simp_schemes()
#>   integrator order gradient_evals
#> 1   leapfrog     2              1
#> 2  minerror2     2              2
#> 3   yoshida4     4              3
#> 4   yoshida6     6              9
#> 5   yoshida8     8             27
```

`trajectory = TRUE` returns every step, and a potential `U` passed to
the call records the total energy along it.

## Sampling with Hamiltonian Monte Carlo

`simp_hmc()` takes the potential and its gradient. For a bivariate
Gaussian with correlation 0.8, the potential uses the precision matrix:

``` r
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
P <- chol2inv(chol(Sigma))
U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
grad_U <- function(q) as.numeric(P %*% q)

fit <- simp_hmc(U, grad_U, init = c(0, 0), n_iter = 2000, eps = 0.25,
                n_steps = 16, warmup = 500, integrator = "minerror2")
colMeans(fit$draws)
#> [1] -0.0005825939 -0.0005670790
fit$accept_rate
#> [1] 1
```

The sample means are close to the target mean of zero. [Symplectic
integration with SIMP](https://gcol33.github.io/SIMP/articles/SIMP.html)
compares the energy error of leapfrog, `minerror2` and `yoshida4` at the
same step size on this kind of target.

## Non-separable and split Hamiltonians

When the kinetic and potential parts do not separate, as with a
position-dependent mass, `simp_integrate_nonsep()` solves the two
implicit stages of the generalized leapfrog by fixed-point iteration.
For a separable Hamiltonian it reduces to the standard leapfrog.
`simp_integrate_mts()` handles a potential split into a stiff, cheap
part and a smooth, expensive part, with `m` inner substeps of the fast
part per outer step. [Symplectic integration with
SIMP](https://gcol33.github.io/SIMP/articles/SIMP.html) runs the
generalized leapfrog on a mass of `1 + 0.5 q^2`.

## Documentation

Start with [Symplectic integration with
SIMP](https://gcol33.github.io/SIMP/articles/SIMP.html): the schemes,
one trajectory, energy conservation, Hamiltonian Monte Carlo and a
non-separable Hamiltonian.

- [Full reference](https://gcol33.github.io/SIMP/reference/): every
  function, grouped by job.

## Getting help

Questions and bug reports go to the [GitHub
issues](https://github.com/gcol33/SIMP/issues). A bug report is easiest
to act on with a small reproducible example, for instance made with
[reprex](https://reprex.tidyverse.org/).

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

SIMP is licensed under MIT.

## Citation

To cite SIMP in publications use:

``` text
Colling G (2026). _SIMP: Symplectic Integrators for Monte-Carlo
Proposals_. R package version 0.3.0, <https://gcol33.github.io/SIMP/>.
```

BibTeX:

``` bibtex
@Manual{colling2026simp,
  title = {SIMP: Symplectic Integrators for Monte-Carlo Proposals},
  author = {Gilles Colling},
  year = {2026},
  note = {R package version 0.3.0},
  url = {https://gcol33.github.io/SIMP/},
}
```

Download the citation as
[BibTeX](https://gcol33.github.io/SIMP/SIMP.bib) or
[RIS](https://gcol33.github.io/SIMP/SIMP.ris).
