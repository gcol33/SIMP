# SIMP: Symplectic Integrators for Monte-Carlo Proposals

Symplectic integration of Hamiltonian systems and the Hamiltonian Monte
Carlo proposals built on them. Explicit splitting schemes for separable
Hamiltonians cover leapfrog, a two-stage minimum-error member tuned to
cancel the leading energy error on a Gaussian target, and the Yoshida
order 4, 6, and 8 members generated from a single triple-jump
composition, so the whole family shares one stepper and one coefficient
source. A generalized (implicit) leapfrog handles non-separable
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
