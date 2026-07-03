# Integrate a Hamiltonian trajectory

Advance a separable Hamiltonian system H(q, p) = U(q) + K(p), with
kinetic energy K(p) = 0.5 p' Minv p where Minv is the inverse mass
matrix, using a symplectic splitting integrator. The integrator
conserves a nearby (shadow) energy, so the total energy stays bounded
over long trajectories instead of drifting.

## Usage

``` r
simp_integrate(
  q0,
  p0,
  grad_U,
  eps,
  n_steps,
  integrator = "leapfrog",
  mass = NULL,
  trajectory = FALSE,
  U = NULL
)
```

## Arguments

- q0:

  Numeric vector, the starting position.

- p0:

  Numeric vector, the starting momentum, same length as `q0`.

- grad_U:

  Function mapping a position vector to the gradient of the potential
  energy U, that is dU/dq. For a target density with log-density `lp`,
  the potential is U = -lp and its gradient is -d(lp)/dq.

- eps:

  Positive step size.

- n_steps:

  Number of integrator steps.

- integrator:

  Name of a built-in scheme; see
  [`simp_schemes()`](https://gcol33.github.io/SIMP/reference/simp_schemes.md).
  Defaults to `"leapfrog"`.

- mass:

  Mass matrix M. `NULL` (the default) uses the identity; a numeric
  vector uses a diagonal mass; a matrix uses a dense mass.

- trajectory:

  If `TRUE`, return the position and momentum at every step. If `FALSE`
  (the default), return only the endpoint.

- U:

  Optional function mapping a position to the potential energy U. When
  supplied, the total energy is recorded along the returned points.

## Value

A list with the endpoint `q` and `p`. When `trajectory = TRUE` it also
holds `q_path` and `p_path`, each a matrix with one row per retained
point and one column per coordinate. When `U` is supplied it holds
`energy`, the total energy at the retained points.

## Examples

``` r
# A unit isotropic Gaussian target: U(q) = 0.5 q'q, grad_U(q) = q.
grad_U <- function(q) q
out <- simp_integrate(q0 = c(1, 0), p0 = c(0, 1), grad_U = grad_U,
                      eps = 0.1, n_steps = 60, integrator = "yoshida4")
out$q
#> [1]  0.9601592 -0.2794526
```
