#ifndef SIMP_AUTODIFF_H
#define SIMP_AUTODIFF_H

#include <RcppEigen.h>
#include <cmath>
#include <vector>

// Forward-mode automatic differentiation for gradient computation
// Simplified dual number implementation for SIMP

namespace simp {
namespace ad {

// Dual number for automatic differentiation
class Dual {
public:
  double val;  // Function value
  double grad; // Derivative

  Dual() : val(0.0), grad(0.0) {}
  Dual(double v) : val(v), grad(0.0) {}
  Dual(double v, double g) : val(v), grad(g) {}

  // Assignment
  Dual& operator=(double v) {
    val = v;
    grad = 0.0;
    return *this;
  }

  // Arithmetic operations
  Dual operator+(const Dual& other) const {
    return Dual(val + other.val, grad + other.grad);
  }

  Dual operator-(const Dual& other) const {
    return Dual(val - other.val, grad - other.grad);
  }

  Dual operator*(const Dual& other) const {
    return Dual(val * other.val, val * other.grad + grad * other.val);
  }

  Dual operator/(const Dual& other) const {
    double inv = 1.0 / other.val;
    return Dual(val * inv, (grad - val * inv * other.grad) * inv);
  }

  Dual operator-() const {
    return Dual(-val, -grad);
  }

  // Compound assignment
  Dual& operator+=(const Dual& other) {
    val += other.val;
    grad += other.grad;
    return *this;
  }

  Dual& operator-=(const Dual& other) {
    val -= other.val;
    grad -= other.grad;
    return *this;
  }

  Dual& operator*=(const Dual& other) {
    grad = val * other.grad + grad * other.val;
    val *= other.val;
    return *this;
  }

  // Scalar operations
  Dual operator+(double s) const { return Dual(val + s, grad); }
  Dual operator-(double s) const { return Dual(val - s, grad); }
  Dual operator*(double s) const { return Dual(val * s, grad * s); }
  Dual operator/(double s) const { return Dual(val / s, grad / s); }

  // Comparison
  bool operator<(const Dual& other) const { return val < other.val; }
  bool operator>(const Dual& other) const { return val > other.val; }
  bool operator<=(const Dual& other) const { return val <= other.val; }
  bool operator>=(const Dual& other) const { return val >= other.val; }
};

// Scalar-Dual operations
inline Dual operator+(double s, const Dual& d) { return Dual(s + d.val, d.grad); }
inline Dual operator-(double s, const Dual& d) { return Dual(s - d.val, -d.grad); }
inline Dual operator*(double s, const Dual& d) { return Dual(s * d.val, s * d.grad); }
inline Dual operator/(double s, const Dual& d) {
  double inv = 1.0 / d.val;
  return Dual(s * inv, -s * inv * inv * d.grad);
}

// Mathematical functions
inline Dual exp(const Dual& x) {
  double e = std::exp(x.val);
  return Dual(e, e * x.grad);
}

inline Dual log(const Dual& x) {
  return Dual(std::log(x.val), x.grad / x.val);
}

inline Dual sqrt(const Dual& x) {
  double s = std::sqrt(x.val);
  return Dual(s, 0.5 * x.grad / s);
}

inline Dual pow(const Dual& x, double p) {
  double xp = std::pow(x.val, p);
  return Dual(xp, p * std::pow(x.val, p - 1) * x.grad);
}

inline Dual pow(const Dual& x, const Dual& p) {
  double xp = std::pow(x.val, p.val);
  return Dual(xp, xp * (p.grad * std::log(x.val) + p.val * x.grad / x.val));
}

inline Dual abs(const Dual& x) {
  return x.val >= 0 ? x : -x;
}

inline Dual sin(const Dual& x) {
  return Dual(std::sin(x.val), std::cos(x.val) * x.grad);
}

inline Dual cos(const Dual& x) {
  return Dual(std::cos(x.val), -std::sin(x.val) * x.grad);
}

inline Dual tanh(const Dual& x) {
  double t = std::tanh(x.val);
  return Dual(t, (1 - t * t) * x.grad);
}

// Log-gamma function (lgamma)
inline Dual lgamma(const Dual& x) {
  return Dual(std::lgamma(x.val), R::digamma(x.val) * x.grad);
}

// Digamma function wrapper
inline Dual digamma(const Dual& x) {
  return Dual(R::digamma(x.val), R::trigamma(x.val) * x.grad);
}

// Logit and inverse logit
inline Dual logit(const Dual& x) {
  return log(x / (1.0 - x));
}

inline Dual inv_logit(const Dual& x) {
  double e = std::exp(-x.val);
  double p = 1.0 / (1.0 + e);
  return Dual(p, p * (1 - p) * x.grad);
}

// Softmax (for vectors)
inline std::vector<Dual> softmax(const std::vector<Dual>& x) {
  int n = x.size();
  std::vector<Dual> result(n);

  // Find max for numerical stability
  Dual max_x = x[0];
  for (int i = 1; i < n; i++) {
    if (x[i].val > max_x.val) max_x = x[i];
  }

  // Compute exp(x - max) and sum
  Dual sum(0.0);
  std::vector<Dual> exp_x(n);
  for (int i = 0; i < n; i++) {
    exp_x[i] = exp(x[i] - max_x);
    sum += exp_x[i];
  }

  // Normalize
  for (int i = 0; i < n; i++) {
    result[i] = exp_x[i] / sum;
  }

  return result;
}

// Compute gradient using forward-mode AD
// Evaluates function once per parameter
inline Eigen::VectorXd compute_gradient(
    const Eigen::VectorXd& params,
    std::function<double(const std::vector<Dual>&)> f) {

  int n = params.size();
  Eigen::VectorXd gradient(n);

  for (int i = 0; i < n; i++) {
    std::vector<Dual> dual_params(n);
    for (int j = 0; j < n; j++) {
      dual_params[j] = Dual(params(j), (i == j) ? 1.0 : 0.0);
    }
    Dual result = f(dual_params);
    gradient(i) = result.grad;
  }

  return gradient;
}

} // namespace ad
} // namespace simp

#endif // SIMP_AUTODIFF_H
