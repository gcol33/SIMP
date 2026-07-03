# Hamiltonian Monte Carlo with a chosen symplectic integrator

A minimal Hamiltonian Monte Carlo sampler that draws each proposal by
integrating the Hamiltonian flow of the target with
[`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md)
and accepting it by the Metropolis rule on the total energy. The
integrator is a free choice, so this is the smallest complete setting in
which the schemes can be compared on a real sampling problem.

## Usage

``` r
simp_hmc(
  U,
  grad_U,
  init,
  n_iter,
  eps,
  n_steps,
  integrator = "leapfrog",
  mass = NULL,
  warmup = 0L,
  divergence_threshold = 1000
)
```

## Arguments

- U:

  Function mapping a position to the potential energy, U = -log density
  of the target up to an additive constant.

- grad_U:

  Function mapping a position to the gradient of U.

- init:

  Numeric vector, the starting position.

- n_iter:

  Number of retained iterations.

- eps:

  Positive leapfrog step size.

- n_steps:

  Number of integrator steps per proposal.

- integrator:

  Name of a built-in scheme; see
  [`simp_schemes()`](https://gcol33.github.io/SIMP/reference/simp_schemes.md).

- mass:

  Mass matrix M, as in
  [`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md).
  Momenta are drawn from N(0, M).

- warmup:

  Number of initial iterations to discard.

- divergence_threshold:

  Energy-error magnitude above which a proposal is counted as divergent
  and rejected.

## Value

A list with the `draws` matrix (one row per retained iteration, one
column per coordinate), the `accept_rate`, the number of divergent
proposals `n_divergent`, and the `integrator` used.

## Examples

``` r
# Sample a correlated 2-D Gaussian.
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
P <- solve(Sigma)
U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
grad_U <- function(q) as.numeric(P %*% q)
fit <- simp_hmc(U, grad_U, init = c(0, 0), n_iter = 200, eps = 0.2,
                n_steps = 20, warmup = 100)
colMeans(fit$draws)
#> [1]  0.004591804 -0.010837943
```
