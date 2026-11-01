# Story: SIMP

## Claim
SIMP integrates Hamiltonian systems with symplectic splitting schemes and builds Hamiltonian Monte
Carlo proposals on them, from a header-only, scalar-templated C++ core that other packages can link.

## Audience
Authors of samplers and inference engines (in R or C++) and statisticians who run Hamiltonian
Monte Carlo and want to choose the integrator. They know HMC and gradients; they may not know the
splitting-scheme family beyond leapfrog.

## Relation
- to: none as a parent
- phrase: none
- home: none
- context: none

## Arcs

### A1 One family of schemes, one stepper
- question: "which integrator should drive my HMC, and what does each cost?"
- answer: leapfrog, the minimum-error two-stage member and Yoshida 4, 6 and 8 from one triple-jump
  composition share one stepper; `simp_schemes()` lists order and gradient evaluations per step
- payoff: the reader picks a scheme by name and compares energy errors on the same target
- evidence: vignettes/SIMP.Rmd "Energy conservation"
- home: vignettes/SIMP.Rmd
- also: README, DESCRIPTION

### A2 Coefficients adapted to the step band
- question: "can the scheme be tuned to my target's frequencies?"
- answer: `"adaptive2"` / `"adaptive3"` and `simp_b_optimal()` minimise the worst-case harmonic
  energy error over `(0, nu_max]`
- payoff: the reader passes `nu_max` and gets the minimum-error member for the band
- evidence: ?simp_b_optimal
- home: ?simp_b_optimal
- also: README, DESCRIPTION

### A3 Proposals for HMC
- question: "how do I sample with it?"
- answer: `simp_hmc()` integrates the flow and accepts by the Metropolis rule on the energy, with
  any scheme as `integrator`
- payoff: the reader samples a target
- evidence: vignettes/SIMP.Rmd "Hamiltonian Monte Carlo"
- home: vignettes/SIMP.Rmd
- also: README

### A4 Beyond the separable case
- question: "what if the mass depends on position, or the potential has a stiff part?"
- answer: `simp_integrate_nonsep()` (generalized leapfrog) and `simp_integrate_mts()` (RESPA)
- payoff: the reader integrates Riemannian or split-potential dynamics
- evidence: vignettes/SIMP.Rmd "Non-separable Hamiltonians"; ?simp_integrate_mts
- home: vignettes/SIMP.Rmd
- also: README, DESCRIPTION

### A5 A C++ core for other packages
- question: "can my package use the steppers under automatic differentiation?"
- answer: header-only, scalar-templated core in `inst/include`, reached with `LinkingTo: SIMP`
- payoff: the reader drives the stepper from their own C++ with a native gradient
- evidence: inst/include/simp/
- home: none yet (no article)
- also: README, DESCRIPTION

## Non-goals
- No adaptation of step size or trajectory length (no NUTS).

## Documents
- DESCRIPTION: listing text; claim, A1, A2, A4, A5
- README / landing: entry point; claim, first example, one section per arc, routing
- Get started (SIMP.Rmd): tutorial; schemes, a trajectory, energy, HMC, non-separable

## Open
- A5 has no article showing a consumer package's C++.
- No R-CMD-check workflow in the repository.
