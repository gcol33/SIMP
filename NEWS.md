# SIMP 0.3.0

* Step-adapted minimum-error integrators. `simp_b_optimal()` returns the
  multistage coefficients that minimize the worst-case energy error a splitting
  integrator induces on a harmonic (Gaussian) target across a band of
  dimensionless steps `(0, nu_max]`, where `nu_max = omega_max * eps` follows
  from the highest target frequency (the largest mass-scaled Hessian
  eigenvalue). `simp_integrate()` and `simp_hmc()` take `integrator =
  "adaptive2"` / `"adaptive3"` with `nu_max` to use the tuned two- or
  three-stage scheme. Exposed in C++ as `simp::two_stage_adaptive` /
  `simp::three_stage_adaptive` in `<simp/adapt.h>`.
* Exact harmonic energy-error analysis backing the adaptation, in
  `<simp/harmonic.h>`: `simp::harmonic_energy_error(scheme, nu)` is the relative
  energy-oscillation amplitude the scheme induces on the unit oscillator at
  dimensionless step `nu`, and `simp::harmonic_stability_limit(scheme)` is the
  largest stable step. Both come from the scheme's exact 2x2 propagation matrix
  and its invariant quadratic form.
* Parametric multistage constructors `simp::two_stage(b)` and
  `simp::three_stage(b, a)` in `<simp/scheme.h>`; `minerror2()` is now the
  `b = (3 - sqrt(5)) / 4` member of `two_stage`.
* `simp_integrate_mts()` advances a split potential `U = U_slow + U_fast` with
  the multiple-time-stepping (RESPA / Verlet-I) integrator: the slow force is
  evaluated once per outer step while the fast force drives `m` inner substeps,
  so a stiff-but-cheap part takes small steps without the expensive part paying
  the same rate. Exposed in C++ as `simp::step_mts` / `simp::integrate_mts` in
  `<simp/mts.h>`.

# SIMP 0.2.0

* New integrator `"minerror2"`: a two-stage order-two scheme whose free kick
  weight is fixed at `(3 - sqrt(5)) / 4`, the value that cancels the leading
  energy-error term on a harmonic (Gaussian) target. It takes two gradient
  evaluations per step, conserves energy far better than leapfrog near a
  Gaussian, and admits larger steps without the stability cliff of the
  high-order composition schemes.
* `simp_integrate_nonsep()` integrates non-separable Hamiltonians (for example
  Riemannian dynamics with a position-dependent mass) with the generalized
  (implicit) leapfrog, exposed in C++ as `simp::step_implicit` in
  `<simp/implicit.h>`. For a separable Hamiltonian it reduces to leapfrog.
* Recovery tests cover the new members: the `minerror2` energy cancellation and
  its order on a nonlinear problem, and the generalized leapfrog's
  separable-reduction, energy conservation, and reversibility.

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
