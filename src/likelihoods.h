#ifndef SIMP_LIKELIHOODS_H
#define SIMP_LIKELIHOODS_H

#include <RcppEigen.h>
#include "autodiff.h"
#include <cmath>

namespace simp {
namespace likelihood {

// -----------------------------------------------------------------------------
// Dirichlet distribution
// -----------------------------------------------------------------------------

// Log probability of Dirichlet distribution
// y: observed composition (sums to 1, all positive)
// alpha: concentration parameters (all positive)
inline double log_dirichlet(const Eigen::VectorXd& y,
                            const Eigen::VectorXd& alpha) {
  int K = y.size();
  double log_prob = std::lgamma(alpha.sum());

  for (int k = 0; k < K; k++) {
    log_prob -= std::lgamma(alpha(k));
    if (y(k) > 0) {
      log_prob += (alpha(k) - 1.0) * std::log(y(k));
    } else if (alpha(k) > 1.0) {
      // y(k) = 0 with alpha > 1 gives -infinity
      return -std::numeric_limits<double>::infinity();
    }
    // y(k) = 0 with alpha <= 1 is undefined, but we handle zeros separately
  }

  return log_prob;
}

// Gradient of log Dirichlet w.r.t. alpha
inline Eigen::VectorXd grad_log_dirichlet_alpha(const Eigen::VectorXd& y,
                                                  const Eigen::VectorXd& alpha) {
  int K = y.size();
  Eigen::VectorXd grad(K);

  double digamma_sum = R::digamma(alpha.sum());

  for (int k = 0; k < K; k++) {
    grad(k) = digamma_sum - R::digamma(alpha(k));
    if (y(k) > 0) {
      grad(k) += std::log(y(k));
    }
  }

  return grad;
}

// Dirichlet with softmax parameterization
// y ~ Dirichlet(softmax(eta) * phi)
// eta: linear predictor (unconstrained)
// phi: concentration parameter (positive)
inline double log_dirichlet_softmax(const Eigen::VectorXd& y,
                                     const Eigen::VectorXd& eta,
                                     double phi) {
  int K = y.size();

  // Compute alpha = softmax(eta) * phi
  Eigen::VectorXd alpha(K);
  double max_eta = eta.maxCoeff();
  double sum_exp = 0.0;

  for (int k = 0; k < K; k++) {
    alpha(k) = std::exp(eta(k) - max_eta);
    sum_exp += alpha(k);
  }

  for (int k = 0; k < K; k++) {
    alpha(k) = (alpha(k) / sum_exp) * phi;
  }

  return log_dirichlet(y, alpha);
}

// -----------------------------------------------------------------------------
// Logistic-normal distribution
// -----------------------------------------------------------------------------

// Additive log-ratio transform: alr(y)[k] = log(y[k] / y[ref])
inline Eigen::VectorXd alr_transform(const Eigen::VectorXd& y, int ref = -1) {
  int K = y.size();
  if (ref < 0) ref = K - 1;  // Default: last component as reference

  Eigen::VectorXd alr(K - 1);
  double log_ref = std::log(y(ref));

  int idx = 0;
  for (int k = 0; k < K; k++) {
    if (k != ref) {
      alr(idx++) = std::log(y(k)) - log_ref;
    }
  }

  return alr;
}

// Inverse ALR: softmax of extended vector
inline Eigen::VectorXd inv_alr_transform(const Eigen::VectorXd& alr, int ref = -1) {
  int Km1 = alr.size();
  int K = Km1 + 1;
  if (ref < 0) ref = K - 1;

  // Build extended vector with 0 at reference position
  Eigen::VectorXd eta(K);
  int idx = 0;
  for (int k = 0; k < K; k++) {
    if (k == ref) {
      eta(k) = 0.0;
    } else {
      eta(k) = alr(idx++);
    }
  }

  // Softmax
  double max_eta = eta.maxCoeff();
  Eigen::VectorXd y(K);
  double sum = 0.0;
  for (int k = 0; k < K; k++) {
    y(k) = std::exp(eta(k) - max_eta);
    sum += y(k);
  }
  y /= sum;

  return y;
}

// Log Jacobian of softmax transformation (for ALR -> composition)
// |dy/d(alr)| for change of variables
inline double log_jacobian_softmax(const Eigen::VectorXd& y) {
  double log_jac = 0.0;
  for (int k = 0; k < y.size(); k++) {
    log_jac += std::log(y(k));
  }
  return log_jac;
}

// Log probability of logistic-normal distribution
// y: observed composition
// mu: mean of ALR-transformed variable (K-1 vector)
// L: Cholesky factor of covariance (lower triangular, (K-1) x (K-1))
inline double log_logistic_normal(const Eigen::VectorXd& y,
                                   const Eigen::VectorXd& mu,
                                   const Eigen::MatrixXd& L) {
  int Km1 = mu.size();

  // Transform to ALR space
  Eigen::VectorXd alr_y = alr_transform(y);

  // Standardize: z = L^{-1} (alr_y - mu)
  Eigen::VectorXd centered = alr_y - mu;
  Eigen::VectorXd z = L.triangularView<Eigen::Lower>().solve(centered);

  // Log probability of standard normal
  double log_prob = -0.5 * z.dot(z);

  // Normalizing constant: -0.5 * K * log(2*pi) - log|L|
  log_prob -= 0.5 * Km1 * std::log(2.0 * M_PI);
  for (int i = 0; i < Km1; i++) {
    log_prob -= std::log(L(i, i));
  }

  // Jacobian of ALR -> composition transformation
  log_prob += log_jacobian_softmax(y);

  return log_prob;
}

// Simpler version with diagonal covariance
inline double log_logistic_normal_diag(const Eigen::VectorXd& y,
                                        const Eigen::VectorXd& mu,
                                        const Eigen::VectorXd& sigma) {
  int Km1 = mu.size();

  Eigen::VectorXd alr_y = alr_transform(y);
  Eigen::VectorXd z = (alr_y - mu).cwiseQuotient(sigma);

  double log_prob = -0.5 * z.dot(z);
  log_prob -= 0.5 * Km1 * std::log(2.0 * M_PI);
  log_prob -= sigma.array().log().sum();
  log_prob += log_jacobian_softmax(y);

  return log_prob;
}

// -----------------------------------------------------------------------------
// Hurdle-Dirichlet distribution
// -----------------------------------------------------------------------------

// Two-part model:
// 1. Binary model for P(y_k = 0) for each component
// 2. Dirichlet for non-zero components (renormalized)

inline double log_hurdle_dirichlet(const Eigen::VectorXd& y,
                                    const Eigen::VectorXd& alpha,
                                    const Eigen::VectorXd& zero_prob) {
  int K = y.size();
  double log_prob = 0.0;

  // Count non-zero components
  std::vector<int> nonzero_idx;
  for (int k = 0; k < K; k++) {
    if (y(k) > 0) {
      nonzero_idx.push_back(k);
    }
  }

  if (nonzero_idx.empty()) {
    // All zeros - this shouldn't happen for valid compositions
    return -std::numeric_limits<double>::infinity();
  }

  // Binary part
  for (int k = 0; k < K; k++) {
    if (y(k) == 0) {
      log_prob += std::log(zero_prob(k));
    } else {
      log_prob += std::log(1.0 - zero_prob(k));
    }
  }

  // Dirichlet part for non-zero components
  int n_nonzero = nonzero_idx.size();
  Eigen::VectorXd y_pos(n_nonzero);
  Eigen::VectorXd alpha_pos(n_nonzero);

  double y_sum = 0.0;
  for (int i = 0; i < n_nonzero; i++) {
    int k = nonzero_idx[i];
    y_pos(i) = y(k);
    alpha_pos(i) = alpha(k);
    y_sum += y(k);
  }

  // Renormalize to sum to 1
  y_pos /= y_sum;

  // Dirichlet log probability
  log_prob += log_dirichlet(y_pos, alpha_pos);

  // Jacobian of renormalization
  log_prob -= (n_nonzero - 1) * std::log(y_sum);

  return log_prob;
}

// -----------------------------------------------------------------------------
// Prior distributions
// -----------------------------------------------------------------------------

// Half-normal prior: p(x) \propto exp(-x^2 / (2*sigma^2)) for x > 0
inline double log_half_normal(double x, double sigma) {
  if (x < 0) return -std::numeric_limits<double>::infinity();
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) +
         std::log(2.0) - 0.5 * x * x / (sigma * sigma);
}

// Half-Cauchy prior
inline double log_half_cauchy(double x, double scale) {
  if (x < 0) return -std::numeric_limits<double>::infinity();
  return std::log(2.0) - std::log(M_PI) - std::log(scale) -
         std::log(1.0 + x * x / (scale * scale));
}

// Gamma prior
inline double log_gamma(double x, double shape, double rate) {
  if (x <= 0) return -std::numeric_limits<double>::infinity();
  return shape * std::log(rate) - std::lgamma(shape) +
         (shape - 1.0) * std::log(x) - rate * x;
}

// Normal prior
inline double log_normal(double x, double mean, double sd) {
  double z = (x - mean) / sd;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * z * z;
}

// LKJ prior for correlation matrix (via Cholesky factor)
// eta = 1 is uniform over correlation matrices
inline double log_lkj_cholesky(const Eigen::MatrixXd& L, double eta) {
  int K = L.rows();
  double log_prob = 0.0;

  for (int i = 1; i < K; i++) {
    log_prob += (K - i - 1 + 2.0 * eta - 2.0) * std::log(L(i, i));
  }

  return log_prob;
}

} // namespace likelihood
} // namespace simp

#endif // SIMP_LIKELIHOODS_H
