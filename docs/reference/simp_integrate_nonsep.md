# Integrate a non-separable Hamiltonian (generalized leapfrog)

Advance a Hamiltonian H(q, p) whose kinetic and potential parts do not
separate, as in Riemannian dynamics with a position-dependent mass,
using the generalized (implicit) leapfrog. The step is symmetric, time
reversible, and symplectic; its two implicit stages are solved by
fixed-point iteration. For a separable H = U(q) + K(p) it reduces to
[`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md)
with the leapfrog scheme.

## Usage

``` r
simp_integrate_nonsep(
  q0,
  p0,
  dH_dq,
  dH_dp,
  eps,
  n_steps,
  tol = 1e-10,
  max_iter = 100L,
  trajectory = FALSE,
  H = NULL
)
```

## Arguments

- q0:

  Numeric vector, the starting position.

- p0:

  Numeric vector, the starting momentum, same length as `q0`.

- dH_dq:

  Function of `(q, p)` returning the partial gradient dH/dq.

- dH_dp:

  Function of `(q, p)` returning the partial gradient dH/dp.

- eps:

  Positive step size.

- n_steps:

  Number of integrator steps.

- tol:

  Fixed-point tolerance for the implicit stages.

- max_iter:

  Maximum fixed-point iterations per stage.

- trajectory:

  If `TRUE`, return the position and momentum at every step.

- H:

  Optional function of `(q, p)` returning the total energy, recorded
  along the returned points.

## Value

A list with the endpoint `q` and `p`, a `converged` flag reporting
whether every implicit stage met `tol`, and, when requested, `q_path`,
`p_path`, and `energy`.

## Examples

``` r
# Separable check: H = 0.5 p'p + 0.5 q'q reduces to leapfrog.
out <- simp_integrate_nonsep(
  q0 = c(1, 0), p0 = c(0, 1),
  dH_dq = function(q, p) q, dH_dp = function(q, p) p,
  eps = 0.1, n_steps = 60)
out$q
#> [1]  0.9608666 -0.2773584
```
