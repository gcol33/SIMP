#ifndef SIMP_HMC_CORE_H
#define SIMP_HMC_CORE_H

#include <RcppEigen.h>
#include <random>
#include <functional>

namespace simp {
namespace hmc {

// Function types for log probability and gradient
using LogProbFn = std::function<double(const Eigen::VectorXd&)>;
using GradientFn = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

// Sampler diagnostics
struct StepDiagnostics {
  bool accepted;
  bool divergent;
  int tree_depth;
  double step_size;
  double log_prob;
  double accept_prob;
};

// Sampler output for a single chain
struct ChainOutput {
  Eigen::MatrixXd samples;        // iterations x parameters
  Eigen::VectorXd log_probs;      // log probability at each iteration
  std::vector<StepDiagnostics> diagnostics;
  double mean_accept_rate;
  int n_divergent;
  double final_step_size;
  Eigen::VectorXd final_mass_matrix;
};

// Compute kinetic energy: 0.5 * p' * M^{-1} * p
inline double kinetic_energy(const Eigen::VectorXd& p,
                              const Eigen::VectorXd& inv_mass_diag) {
  return 0.5 * p.cwiseProduct(inv_mass_diag).dot(p);
}

// Compute Hamiltonian: H = -log_prob + kinetic_energy
inline double hamiltonian(double log_prob, const Eigen::VectorXd& p,
                          const Eigen::VectorXd& inv_mass_diag) {
  return -log_prob + kinetic_energy(p, inv_mass_diag);
}

// Single leapfrog step
inline void leapfrog_step(Eigen::VectorXd& q, Eigen::VectorXd& p,
                          double epsilon,
                          const Eigen::VectorXd& inv_mass_diag,
                          const GradientFn& gradient) {
  // Half step for momentum
  Eigen::VectorXd grad = gradient(q);
  p += 0.5 * epsilon * grad;

  // Full step for position
  q += epsilon * p.cwiseProduct(inv_mass_diag);

  // Half step for momentum
  grad = gradient(q);
  p += 0.5 * epsilon * grad;
}

// Multiple leapfrog steps
inline void leapfrog(Eigen::VectorXd& q, Eigen::VectorXd& p,
                     double epsilon, int L,
                     const Eigen::VectorXd& inv_mass_diag,
                     const GradientFn& gradient) {
  for (int i = 0; i < L; i++) {
    leapfrog_step(q, p, epsilon, inv_mass_diag, gradient);
  }
}

// Check for divergence (numerical issues)
inline bool is_divergent(double H_old, double H_new, double max_delta = 1000.0) {
  return std::isnan(H_new) || std::isinf(H_new) ||
         (H_new - H_old) > max_delta;
}

// Metropolis acceptance
inline bool metropolis_accept(double H_old, double H_new, std::mt19937& rng) {
  double log_accept = H_old - H_new;
  if (log_accept >= 0.0) return true;

  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  return std::log(uniform(rng)) < log_accept;
}

// Acceptance probability (clamped to [0, 1])
inline double acceptance_probability(double H_old, double H_new) {
  double log_accept = H_old - H_new;
  return std::min(1.0, std::exp(log_accept));
}

// Dual averaging for step size adaptation
class StepSizeAdapter {
public:
  double mu;           // log(10 * epsilon_0)
  double log_eps_bar;  // log(smoothed step size)
  double H_bar;        // smoothed accept stat
  double gamma;        // shrinkage toward mu
  double t0;           // iteration offset
  double kappa;        // decay rate
  double target_accept;

  StepSizeAdapter(double eps0 = 1.0, double target = 0.8)
    : mu(std::log(10.0 * eps0)),
      log_eps_bar(0.0),
      H_bar(0.0),
      gamma(0.05),
      t0(10.0),
      kappa(0.75),
      target_accept(target) {}

  double update(double accept_prob, int iter) {
    double w = 1.0 / (iter + t0);
    H_bar = (1.0 - w) * H_bar + w * (target_accept - accept_prob);

    double log_eps = mu - std::sqrt(iter) / gamma * H_bar;
    double eps = std::exp(log_eps);

    double iter_kappa = std::pow(iter, -kappa);
    log_eps_bar = iter_kappa * log_eps + (1.0 - iter_kappa) * log_eps_bar;

    return eps;
  }

  double final_step_size() const {
    return std::exp(log_eps_bar);
  }
};

// Mass matrix adaptation (diagonal)
class MassMatrixAdapter {
public:
  int n_params;
  Eigen::VectorXd mean;
  Eigen::VectorXd M2;  // Sum of squared deviations
  int n_samples;

  MassMatrixAdapter(int n) : n_params(n), n_samples(0) {
    mean = Eigen::VectorXd::Zero(n);
    M2 = Eigen::VectorXd::Zero(n);
  }

  void add_sample(const Eigen::VectorXd& sample) {
    n_samples++;
    Eigen::VectorXd delta = sample - mean;
    mean += delta / n_samples;
    Eigen::VectorXd delta2 = sample - mean;
    M2 += delta.cwiseProduct(delta2);
  }

  Eigen::VectorXd get_variance() const {
    if (n_samples < 2) {
      return Eigen::VectorXd::Ones(n_params);
    }
    Eigen::VectorXd var = M2 / (n_samples - 1);
    // Regularize small variances
    for (int i = 0; i < n_params; i++) {
      var(i) = std::max(var(i), 1e-8);
    }
    return var;
  }

  Eigen::VectorXd get_inv_mass_diag() const {
    return get_variance().cwiseInverse();
  }
};

// Sample momentum from standard normal * sqrt(mass)
inline Eigen::VectorXd sample_momentum(int n, const Eigen::VectorXd& mass_diag,
                                        std::mt19937& rng) {
  std::normal_distribution<double> normal(0.0, 1.0);
  Eigen::VectorXd p(n);
  for (int i = 0; i < n; i++) {
    p(i) = normal(rng) * std::sqrt(mass_diag(i));
  }
  return p;
}

} // namespace hmc
} // namespace simp

#endif // SIMP_HMC_CORE_H
