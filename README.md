# SIMP

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**S**ymplectic **I**ntegrators for **M**onte-carlo **P**roposals.

SIMP integrates separable Hamiltonian systems and builds the Hamiltonian Monte
Carlo proposals that run on them. Symplectic integrators conserve a nearby
energy, so a trajectory holds its total energy over long horizons and the
proposals it produces are accepted at high rates.

The integrators are splitting schemes: an alternating sequence of position and
momentum updates whose coefficients define the method. Leapfrog and the Yoshida
order 4, 6, and 8 members come from a single triple-jump composition, so the
whole family shares one stepper and one coefficient source. Adding a scheme is
one registry entry.

The stepper is a header-only, scalar-templated C++ core. Consumers link against
it with `LinkingTo: SIMP` and drive it with a native gradient under automatic
differentiation, so the same code path serves plain-double integration and
gradient-based inference engines.

## Install

```r
# install.packages("remotes")
remotes::install_github("gcol33/SIMP")
```

## Integrate a trajectory

```r
library(SIMP)

# Unit Gaussian target: U(q) = 0.5 q'q, so grad_U(q) = q.
grad_U <- function(q) q

out <- simp_integrate(
  q0 = c(1, 0), p0 = c(0, 1),
  grad_U = grad_U, eps = 0.1, n_steps = 60,
  integrator = "yoshida4"
)
out$q
```

`simp_schemes()` lists the built-in integrators with their order and gradient
cost per step.

## Sample with HMC

```r
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
P <- solve(Sigma)
U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
grad_U <- function(q) as.numeric(P %*% q)

fit <- simp_hmc(U, grad_U, init = c(0, 0), n_iter = 2000, eps = 0.2,
                n_steps = 20, warmup = 500, integrator = "leapfrog")

colMeans(fit$draws)
fit$accept_rate
```

Pass any name from `simp_schemes()` as `integrator` to compare schemes on the
same target.

## License

MIT (c) Gilles Colling.
