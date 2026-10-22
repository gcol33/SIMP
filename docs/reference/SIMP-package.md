# SIMP: Symplectic Integrators for Monte-Carlo Proposals

Symplectic integration of Hamiltonian systems and the Hamiltonian Monte
Carlo proposals built on them. Explicit splitting schemes for separable
Hamiltonians cover leapfrog, minimum-error two- and three-stage members,
and the Yoshida order 4, 6, and 8 members generated from a single
triple-jump composition, so the whole family shares one stepper and one
coefficient source. The multistage coefficients can be adapted to a
target's step band from an exact harmonic energy-error analysis, and a
multiple-time-stepping (RESPA) integrator advances a split potential at
two rates. A generalized (implicit) leapfrog handles non-separable
Hamiltonians, as in Riemannian dynamics with a position-dependent mass.
The steppers are a header-only, scalar-templated C++ core that other
packages can link against, so they run under automatic differentiation
as well as plain double precision.

## See also

Useful links:

- <https://gcol33.github.io/SIMP/>

- <https://github.com/gcol33/SIMP>

- Report bugs at <https://github.com/gcol33/SIMP/issues>

## Author

**Maintainer**: Gilles Colling <gilles.colling051@gmail.com>

Authors:

- Gilles Colling <gilles.colling051@gmail.com>
