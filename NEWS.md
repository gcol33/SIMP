# SIMP 0.1.0

First release of SIMP as a symplectic-integration library.

* `simp_integrate()` advances a separable Hamiltonian system with a chosen
  symplectic splitting scheme, returning the endpoint or the full trajectory
  and, optionally, the total energy along it.
* `simp_schemes()` lists the built-in integrators: leapfrog and the Yoshida
  order 4, 6, and 8 members, generated from leapfrog by triple-jump
  composition.
* `simp_hmc()` runs Hamiltonian Monte Carlo with any of the schemes.
* Header-only, scalar-templated C++ core in `inst/include/simp/` for
  `LinkingTo: SIMP` consumers.
