#ifndef SIMP_INTEGRATE_H
#define SIMP_INTEGRATE_H

#include "scheme.h"
#include <Eigen/Dense>

namespace simp {

// One ABA step for a separable Hamiltonian. gradU(q) returns dU/dq and
// applyMinv(p) returns M^{-1} p; both are caller-supplied, so a diagonal,
// dense, or matrix-free mass all use this one stepper, as does any scalar T
// (double or an automatic-differentiation type) through the Eigen vectors.
// The mass is constant in q, which keeps H separable and the step symplectic.
template <typename T, typename GradU, typename ApplyMinv>
void step(Eigen::Matrix<T, Eigen::Dynamic, 1>& q,
          Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
          const T& eps, const Scheme& s,
          GradU&& gradU, ApplyMinv&& applyMinv) {
  const int m = s.stages();
  for (int j = 0; j < m; ++j) {
    q += (T(s.a[j]) * eps) * applyMinv(p);
    p -= (T(s.b[j]) * eps) * gradU(q);
  }
  q += (T(s.a[m]) * eps) * applyMinv(p);
}

// Advance n_steps steps in place.
template <typename T, typename GradU, typename ApplyMinv>
void integrate(Eigen::Matrix<T, Eigen::Dynamic, 1>& q,
               Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
               const T& eps, int n_steps, const Scheme& s,
               GradU&& gradU, ApplyMinv&& applyMinv) {
  for (int i = 0; i < n_steps; ++i) {
    step(q, p, eps, s, gradU, applyMinv);
  }
}

}  // namespace simp

#endif  // SIMP_INTEGRATE_H
