#ifndef SIMP_SCHEME_H
#define SIMP_SCHEME_H

#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace simp {

// An explicit symplectic integrator for a separable Hamiltonian
// H(q, p) = U(q) + K(p), stored as an ABA (drift-first, drift-last) splitting.
// One step of size eps advances (q, p) by
//
//   for j in 0 .. s-1:
//     q += a[j] * eps * Minv * p
//     p -= b[j] * eps * gradU(q)
//   q += a[s] * eps * Minv * p
//
// a holds s+1 drift coefficients, b holds s kick coefficients. Both are
// palindromic and normalised to sum to one, which makes the map time
// reversible and consistent. The kick count s equals the number of gradient
// evaluations per step, the cost unit for comparing schemes.
struct Scheme {
  std::string name;
  int order;
  std::vector<double> a;
  std::vector<double> b;

  int stages() const { return static_cast<int>(b.size()); }
};

// Leapfrog (Stormer-Verlet): order two, one gradient per step.
inline Scheme leapfrog() {
  return Scheme{"leapfrog", 2, {0.5, 0.5}, {1.0}};
}

// Scale a scheme so one step covers a fraction w of the step size. A drift or
// kick of size w * eps carries weight w * coefficient.
inline Scheme scale(const Scheme& s, double w) {
  Scheme out = s;
  for (double& x : out.a) x *= w;
  for (double& x : out.b) x *= w;
  return out;
}

// Run p then q, merging the trailing drift of p with the leading drift of q
// (two consecutive drifts add). The result is again in ABA form.
inline Scheme concat(const Scheme& p, const Scheme& q) {
  const int sp = p.stages();
  Scheme out;
  out.a = std::vector<double>(p.a.begin(), p.a.begin() + sp);   // p.a[0 .. sp-1]
  out.a.push_back(p.a[sp] + q.a[0]);                            // merged seam
  out.a.insert(out.a.end(), q.a.begin() + 1, q.a.end());        // q.a[1 .. end]
  out.b = p.b;
  out.b.insert(out.b.end(), q.b.begin(), q.b.end());
  return out;
}

// Yoshida triple jump: three scaled copies of a symmetric scheme of even order
// n raise the order to n+2. The outer copies use w1 = 1 / (2 - 2^(1/(n+1)))
// and the inner copy uses w0 = 1 - 2 * w1 (negative). The construction is
// exact; the coefficients are functions of n, not tabulated constants.
inline Scheme triple_jump(const Scheme& base, const std::string& name) {
  const double n = static_cast<double>(base.order);
  const double w1 = 1.0 / (2.0 - std::pow(2.0, 1.0 / (n + 1.0)));
  const double w0 = 1.0 - 2.0 * w1;
  Scheme out = concat(concat(scale(base, w1), scale(base, w0)), scale(base, w1));
  out.name = name;
  out.order = base.order + 2;
  return out;
}

// The built-in family: leapfrog plus the Yoshida members it generates.
inline std::vector<Scheme> registry() {
  Scheme lf = leapfrog();
  Scheme y4 = triple_jump(lf, "yoshida4");
  Scheme y6 = triple_jump(y4, "yoshida6");
  Scheme y8 = triple_jump(y6, "yoshida8");
  return {lf, y4, y6, y8};
}

inline Scheme scheme_by_name(const std::string& name) {
  for (const Scheme& s : registry()) {
    if (s.name == name) return s;
  }
  throw std::invalid_argument("unknown integrator: " + name);
}

}  // namespace simp

#endif  // SIMP_SCHEME_H
