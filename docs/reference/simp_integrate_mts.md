# Integrate a split-potential Hamiltonian (multiple time-stepping)

Advance a separable Hamiltonian whose potential splits as U = U_slow +
U_fast, using the RESPA / Verlet-I multiple-time-stepping integrator.
The slow force is evaluated once per outer step while the fast force
drives `m` inner leapfrog substeps, so a stiff-but-cheap part (a
Gaussian latent prior, say) can take the small steps it needs without
paying for the expensive-but-smooth part (a likelihood) at the same
rate. The step is symmetric, time reversible, and symplectic.

## Usage

``` r
simp_integrate_mts(
  q0,
  p0,
  grad_slow,
  grad_fast,
  eps,
  n_steps,
  m,
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

- grad_slow, grad_fast:

  Functions mapping a position to the gradient of the slow and fast
  potential parts. Their sum is dU/dq.

- eps:

  Positive outer step size.

- n_steps:

  Number of outer steps.

- m:

  Number of inner fast substeps per outer step.

- mass:

  Mass matrix M, as in
  [`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md).

- trajectory:

  If `TRUE`, return the position and momentum at every outer step. If
  `FALSE` (the default), return only the endpoint.

- U:

  Optional function mapping a position to the total potential energy
  U_slow + U_fast. When supplied, the total energy is recorded along the
  returned points.

## Value

A list with the endpoint `q` and `p`, and, when requested, `q_path`,
`p_path`, and `energy`, as in
[`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md).

## Examples

``` r
# Stiff Gaussian prior (fast) plus a soft quadratic likelihood (slow).
grad_fast <- function(q) 100 * q
grad_slow <- function(q) 0.1 * q
out <- simp_integrate_mts(q0 = c(1, 0), p0 = c(0, 1),
                          grad_slow = grad_slow, grad_fast = grad_fast,
                          eps = 0.2, n_steps = 30, m = 10)
out$q
```
