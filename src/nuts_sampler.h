#ifndef SIMP_NUTS_SAMPLER_H
#define SIMP_NUTS_SAMPLER_H

#include "hmc_core.h"
#include <random>
#include <cmath>

namespace simp {
namespace hmc {

// NUTS (No-U-Turn Sampler) implementation
class NUTSSampler {
public:
  LogProbFn log_prob_fn;
  GradientFn gradient_fn;
  int n_params;

  // Sampler settings
  int max_treedepth;
  double adapt_delta;

  // State
  Eigen::VectorXd current_q;
  double current_log_prob;
  Eigen::VectorXd mass_diag;
  Eigen::VectorXd inv_mass_diag;

  // Adaptation
  StepSizeAdapter step_adapter;
  MassMatrixAdapter mass_adapter;

  // RNG
  std::mt19937 rng;

  NUTSSampler(LogProbFn log_prob, GradientFn gradient, int n_params,
              int max_depth = 10, double target_accept = 0.8)
    : log_prob_fn(log_prob),
      gradient_fn(gradient),
      n_params(n_params),
      max_treedepth(max_depth),
      adapt_delta(target_accept),
      step_adapter(1.0, target_accept),
      mass_adapter(n_params) {
    mass_diag = Eigen::VectorXd::Ones(n_params);
    inv_mass_diag = Eigen::VectorXd::Ones(n_params);
  }

  void set_seed(unsigned int seed) {
    rng.seed(seed);
  }

  void set_initial_position(const Eigen::VectorXd& q0) {
    current_q = q0;
    current_log_prob = log_prob_fn(q0);
  }

  // Build tree recursively (simplified NUTS)
  struct TreeNode {
    Eigen::VectorXd q_minus;
    Eigen::VectorXd q_plus;
    Eigen::VectorXd p_minus;
    Eigen::VectorXd p_plus;
    Eigen::VectorXd q_prime;  // Proposed sample
    double log_prob_prime;
    int n_valid;              // Number of valid leaves
    bool divergent;
    bool stop;                // U-turn detected
  };

  // Check U-turn criterion
  bool check_uturn(const Eigen::VectorXd& p_minus,
                   const Eigen::VectorXd& p_plus,
                   const Eigen::VectorXd& q_diff) {
    return (q_diff.dot(p_minus) < 0) || (q_diff.dot(p_plus) < 0);
  }

  // Build tree of depth j
  TreeNode build_tree(const Eigen::VectorXd& q, const Eigen::VectorXd& p,
                      int direction, int depth, double epsilon,
                      double log_slice, double H0, double log_prob_current) {
    TreeNode result;
    result.divergent = false;
    result.stop = false;
    result.n_valid = 0;
    // Initialize endpoint vectors to input values (prevents uninitialized vectors on early return)
    result.q_minus = q;
    result.q_plus = q;
    result.p_minus = p;
    result.p_plus = p;
    result.q_prime = q;
    result.log_prob_prime = log_prob_current;

    if (depth == 0) {
      // Base case: single leapfrog step
      Eigen::VectorXd q_new = q;
      Eigen::VectorXd p_new = p;
      leapfrog_step(q_new, p_new, direction * epsilon, inv_mass_diag, gradient_fn);

      double log_prob_new = log_prob_fn(q_new);
      double H_new = hamiltonian(log_prob_new, p_new, inv_mass_diag);

      // Check for divergence - return with input vectors as endpoints
      if (is_divergent(H0, H_new)) {
        result.divergent = true;
        result.stop = true;
        result.n_valid = 0;
        return result;
      }

      // Check slice condition: H_new <= log_slice (i.e., -log_prob - KE <= slice)
      result.n_valid = (H_new <= log_slice) ? 1 : 0;

      result.q_minus = q_new;
      result.q_plus = q_new;
      result.p_minus = p_new;
      result.p_plus = p_new;
      result.q_prime = q_new;
      result.log_prob_prime = log_prob_new;

      return result;
    }

    // Recursive case: build left subtree
    TreeNode left = build_tree(q, p, direction, depth - 1, epsilon, log_slice, H0, log_prob_current);

    if (left.stop || left.divergent) {
      return left;
    }

    // Build right subtree from appropriate endpoint
    Eigen::VectorXd q_start, p_start;
    if (direction == 1) {
      q_start = left.q_plus;
      p_start = left.p_plus;
    } else {
      q_start = left.q_minus;
      p_start = left.p_minus;
    }

    TreeNode right = build_tree(q_start, p_start, direction, depth - 1,
                                 epsilon, log_slice, H0, left.log_prob_prime);

    if (right.stop || right.divergent) {
      result.divergent = right.divergent;
      result.stop = true;
      result.q_prime = left.q_prime;
      result.log_prob_prime = left.log_prob_prime;
      result.n_valid = left.n_valid;
      return result;
    }

    // Combine trees
    int total_valid = left.n_valid + right.n_valid;

    // Select sample
    if (total_valid > 0) {
      std::uniform_real_distribution<double> uniform(0.0, 1.0);
      if (uniform(rng) < (double)right.n_valid / total_valid) {
        result.q_prime = right.q_prime;
        result.log_prob_prime = right.log_prob_prime;
      } else {
        result.q_prime = left.q_prime;
        result.log_prob_prime = left.log_prob_prime;
      }
    } else {
      result.q_prime = left.q_prime;
      result.log_prob_prime = left.log_prob_prime;
    }

    result.n_valid = total_valid;

    // Update endpoints
    if (direction == 1) {
      result.q_minus = left.q_minus;
      result.p_minus = left.p_minus;
      result.q_plus = right.q_plus;
      result.p_plus = right.p_plus;
    } else {
      result.q_minus = right.q_minus;
      result.p_minus = right.p_minus;
      result.q_plus = left.q_plus;
      result.p_plus = left.p_plus;
    }

    // Check U-turn for combined tree
    Eigen::VectorXd q_diff = result.q_plus - result.q_minus;
    result.stop = check_uturn(result.p_minus, result.p_plus, q_diff);

    return result;
  }

  // Single NUTS step
  StepDiagnostics step(double epsilon, bool adapt, int iter) {
    StepDiagnostics diag;
    diag.step_size = epsilon;

    // Sample momentum
    Eigen::VectorXd p0 = sample_momentum(n_params, mass_diag, rng);

    // Compute initial Hamiltonian
    double H0 = hamiltonian(current_log_prob, p0, inv_mass_diag);

    // Sample slice variable
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double log_slice = H0 - std::log(uniform(rng));

    // Initialize tree
    Eigen::VectorXd q_minus = current_q;
    Eigen::VectorXd q_plus = current_q;
    Eigen::VectorXd p_minus = p0;
    Eigen::VectorXd p_plus = p0;

    Eigen::VectorXd q_prop = current_q;
    double log_prob_prop = current_log_prob;
    int n_valid = 1;

    int depth = 0;
    bool stop = false;
    diag.divergent = false;

    // Build tree until U-turn or max depth
    while (!stop && depth < max_treedepth) {
      // Choose direction
      int direction = (uniform(rng) < 0.5) ? -1 : 1;

      // Build tree in chosen direction
      TreeNode subtree;
      if (direction == 1) {
        subtree = build_tree(q_plus, p_plus, direction, depth, epsilon, log_slice, H0, log_prob_prop);
        // Only update endpoints if subtree is valid
        if (!subtree.divergent && !subtree.stop) {
          q_plus = subtree.q_plus;
          p_plus = subtree.p_plus;
        }
      } else {
        subtree = build_tree(q_minus, p_minus, direction, depth, epsilon, log_slice, H0, log_prob_prop);
        // Only update endpoints if subtree is valid
        if (!subtree.divergent && !subtree.stop) {
          q_minus = subtree.q_minus;
          p_minus = subtree.p_minus;
        }
      }

      if (subtree.divergent) {
        diag.divergent = true;
        stop = true;
      } else if (!subtree.stop && subtree.n_valid > 0) {
        // Accept new sample with probability n_valid_new / n_valid_total
        if (uniform(rng) < (double)subtree.n_valid / (n_valid + subtree.n_valid)) {
          q_prop = subtree.q_prime;
          log_prob_prop = subtree.log_prob_prime;
        }
        n_valid += subtree.n_valid;
      }

      // Check U-turn for full tree
      Eigen::VectorXd q_diff = q_plus - q_minus;
      stop = stop || subtree.stop || check_uturn(p_minus, p_plus, q_diff);

      depth++;
    }

    diag.tree_depth = depth;

    // Accept proposal
    current_q = q_prop;
    current_log_prob = log_prob_prop;
    diag.log_prob = current_log_prob;
    diag.accepted = true;  // NUTS always accepts (built into algorithm)

    // Compute acceptance probability for adaptation
    // If we diverged, acceptance probability should be 0 (to decrease step size)
    if (diag.divergent) {
      diag.accept_prob = 0.0;
    } else {
      double H_prop = hamiltonian(log_prob_prop, p0, inv_mass_diag);
      diag.accept_prob = acceptance_probability(H0, H_prop);
    }

    return diag;
  }

  // Run sampler
  ChainOutput run(int n_iter, int n_warmup, double init_epsilon = 0.01,
                  bool verbose = true) {
    ChainOutput output;
    output.samples.resize(n_iter, n_params);
    output.log_probs.resize(n_iter);
    output.diagnostics.reserve(n_iter);
    output.n_divergent = 0;

    double epsilon = init_epsilon;
    double accept_sum = 0.0;

    // Warmup phase
    int adapt_window_start = 0;
    int adapt_window_end = n_warmup;

    for (int i = 0; i < n_iter + n_warmup; i++) {
      bool in_warmup = (i < n_warmup);

      // Take NUTS step
      StepDiagnostics diag = step(epsilon, in_warmup, i + 1);

      if (in_warmup) {
        // Adapt step size
        epsilon = step_adapter.update(diag.accept_prob, i + 1);

        // Adapt mass matrix (during window)
        if (i >= adapt_window_start && i < adapt_window_end) {
          mass_adapter.add_sample(current_q);
        }

        // Update mass matrix at end of warmup
        if (i == n_warmup - 1) {
          mass_diag = mass_adapter.get_variance();
          inv_mass_diag = mass_diag.cwiseInverse();
          epsilon = step_adapter.final_step_size();
        }
      } else {
        // Sampling phase
        int sample_idx = i - n_warmup;
        output.samples.row(sample_idx) = current_q;
        output.log_probs(sample_idx) = current_log_prob;
        output.diagnostics.push_back(diag);

        if (diag.divergent) output.n_divergent++;
        accept_sum += diag.accept_prob;
      }

      // Progress
      if (verbose && ((i + 1) % 100 == 0 || i == n_iter + n_warmup - 1)) {
        if (in_warmup) {
          Rcpp::Rcout << "Warmup: " << (i + 1) << "/" << n_warmup
                      << " (eps=" << epsilon << ")" << std::endl;
        } else {
          Rcpp::Rcout << "Sampling: " << (i - n_warmup + 1) << "/" << n_iter
                      << std::endl;
        }
      }
    }

    output.mean_accept_rate = accept_sum / n_iter;
    output.final_step_size = epsilon;
    output.final_mass_matrix = mass_diag;

    return output;
  }
};

} // namespace hmc
} // namespace simp

#endif // SIMP_NUTS_SAMPLER_H
