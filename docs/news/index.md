# Changelog

## SIMP 0.2.0

- New integrator `"minerror2"`: a two-stage order-two scheme whose free
  kick weight is fixed at `(3 - sqrt(5)) / 4`, the value that cancels
  the leading energy-error term on a harmonic (Gaussian) target. It
  takes two gradient evaluations per step, conserves energy far better
  than leapfrog near a Gaussian, and admits larger steps without the
  stability cliff of the high-order composition schemes.
- [`simp_integrate_nonsep()`](https://gcol33.github.io/SIMP/reference/simp_integrate_nonsep.md)
  integrates non-separable Hamiltonians (for example Riemannian dynamics
  with a position-dependent mass) with the generalized (implicit)
  leapfrog, exposed in C++ as `simp::step_implicit` in
  `<simp/implicit.h>`. For a separable Hamiltonian it reduces to
  leapfrog.
- Recovery tests cover the new members: the `minerror2` energy
  cancellation and its order on a nonlinear problem, and the generalized
  leapfrog’s separable-reduction, energy conservation, and
  reversibility.

## SIMP 0.1.0

First release of SIMP as a symplectic-integration library.

- [`simp_integrate()`](https://gcol33.github.io/SIMP/reference/simp_integrate.md)
  advances a separable Hamiltonian system with a chosen symplectic
  splitting scheme, returning the endpoint or the full trajectory and,
  optionally, the total energy along it.
- [`simp_schemes()`](https://gcol33.github.io/SIMP/reference/simp_schemes.md)
  lists the built-in integrators: leapfrog and the Yoshida order 4, 6,
  and 8 members, generated from leapfrog by triple-jump composition.
- [`simp_hmc()`](https://gcol33.github.io/SIMP/reference/simp_hmc.md)
  runs Hamiltonian Monte Carlo with any of the schemes.
- Header-only, scalar-templated C++ core in `inst/include/simp/` for
  `LinkingTo: SIMP` consumers.
