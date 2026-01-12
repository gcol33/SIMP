// SIMP sampler interface for R
// Exposes NUTS sampler via Rcpp

#include <Rcpp.h>
#include <RcppEigen.h>
#include "nuts_sampler.h"
#include "likelihoods.h"

using namespace Rcpp;
using namespace simp;

// [[Rcpp::depends(RcppEigen)]]

// -----------------------------------------------------------------------------
// Model class that combines likelihood and priors
// -----------------------------------------------------------------------------

class DirichletModel {
public:
  Eigen::MatrixXd Y;      // Response matrix (n_obs x K)
  Eigen::MatrixXd X;      // Design matrix (n_obs x p)
  int n_obs;
  int K;                  // Number of components
  int p;                  // Number of fixed effects

  // Priors
  double beta_prior_sd;
  double phi_prior_shape;
  double phi_prior_rate;

  DirichletModel(const Eigen::MatrixXd& Y_, const Eigen::MatrixXd& X_,
                 double beta_sd = 2.5, double phi_shape = 2.0, double phi_rate = 0.1)
    : Y(Y_), X(X_),
      n_obs(Y_.rows()), K(Y_.cols()), p(X_.cols()),
      beta_prior_sd(beta_sd),
      phi_prior_shape(phi_shape), phi_prior_rate(phi_rate) {}

  // Parameter layout:
  // [0, p*K): beta coefficients (p predictors x K components)
  // [p*K]: log(phi) - concentration parameter

  int n_params() const {
    return p * K + 1;
  }

  // Extract parameters
  Eigen::MatrixXd get_beta(const Eigen::VectorXd& params) const {
    Eigen::MatrixXd beta(p, K);
    for (int j = 0; j < K; j++) {
      for (int i = 0; i < p; i++) {
        beta(i, j) = params(i * K + j);
      }
    }
    return beta;
  }

  double get_phi(const Eigen::VectorXd& params) const {
    return std::exp(params(p * K));
  }

  // Log posterior
  double log_posterior(const Eigen::VectorXd& params) const {
    Eigen::MatrixXd beta = get_beta(params);
    double log_phi = params(p * K);
    double phi = std::exp(log_phi);

    double log_post = 0.0;

    // Likelihood for each observation
    for (int i = 0; i < n_obs; i++) {
      // Linear predictor: eta = X * beta (result is 1 x K row vector)
      Eigen::RowVectorXd eta_row = X.row(i) * beta;
      Eigen::VectorXd eta = eta_row.transpose();

      // Get response as column vector
      Eigen::VectorXd y_i = Y.row(i).transpose();

      // Dirichlet likelihood with softmax parameterization
      log_post += likelihood::log_dirichlet_softmax(y_i, eta, phi);
    }

    // Priors on beta (normal)
    for (int j = 0; j < p * K; j++) {
      log_post += likelihood::log_normal(params(j), 0.0, beta_prior_sd);
    }

    // Prior on phi (gamma, with Jacobian for log transform)
    log_post += likelihood::log_gamma(phi, phi_prior_shape, phi_prior_rate);
    log_post += log_phi;  // Jacobian of log transform

    return log_post;
  }

  // Gradient of log posterior (numerical for now)
  Eigen::VectorXd gradient(const Eigen::VectorXd& params) const {
    int n = params.size();
    Eigen::VectorXd grad(n);
    double eps = 1e-6;

    double f0 = log_posterior(params);

    for (int i = 0; i < n; i++) {
      Eigen::VectorXd params_plus = params;
      params_plus(i) += eps;
      double f_plus = log_posterior(params_plus);
      grad(i) = (f_plus - f0) / eps;
    }

    return grad;
  }
};

// -----------------------------------------------------------------------------
// Rcpp interface functions
// -----------------------------------------------------------------------------

//' Fit Dirichlet model using NUTS
//'
//' Internal function to run NUTS sampler for Dirichlet compositional model.
//'
//' @param Y Response matrix (n_obs x K)
//' @param X Design matrix (n_obs x p)
//' @param n_iter Number of sampling iterations
//' @param n_warmup Number of warmup iterations
//' @param init Initial parameter values
//' @param seed Random seed
//' @param max_treedepth Maximum tree depth for NUTS
//' @param adapt_delta Target acceptance rate
//' @param beta_prior_sd Prior SD for beta coefficients
//' @param phi_prior_shape Gamma prior shape for phi
//' @param phi_prior_rate Gamma prior rate for phi
//' @param verbose Print progress
//' @return List with samples, log_probs, and diagnostics
//' @noRd
// [[Rcpp::export]]
List cpp_fit_dirichlet(
    Eigen::MatrixXd Y,
    Eigen::MatrixXd X,
    int n_iter,
    int n_warmup,
    Eigen::VectorXd init,
    int seed,
    int max_treedepth,
    double adapt_delta,
    double beta_prior_sd,
    double phi_prior_shape,
    double phi_prior_rate,
    bool verbose) {

  // Create model
  DirichletModel model(Y, X, beta_prior_sd, phi_prior_shape, phi_prior_rate);

  // Check dimensions
  if (init.size() != model.n_params()) {
    stop("Initial values must have length %d", model.n_params());
  }

  // Create log_prob and gradient functions
  auto log_prob_fn = [&model](const Eigen::VectorXd& params) -> double {
    return model.log_posterior(params);
  };

  auto gradient_fn = [&model](const Eigen::VectorXd& params) -> Eigen::VectorXd {
    return model.gradient(params);
  };

  // Create sampler
  hmc::NUTSSampler sampler(log_prob_fn, gradient_fn, model.n_params(),
                            max_treedepth, adapt_delta);

  sampler.set_seed(seed);
  sampler.set_initial_position(init);

  // Run sampler
  hmc::ChainOutput output = sampler.run(n_iter, n_warmup, 1.0, verbose);

  // Extract diagnostics
  IntegerVector tree_depths(n_iter);
  LogicalVector divergent(n_iter);
  NumericVector accept_probs(n_iter);

  for (int i = 0; i < n_iter; i++) {
    tree_depths[i] = output.diagnostics[i].tree_depth;
    divergent[i] = output.diagnostics[i].divergent;
    accept_probs[i] = output.diagnostics[i].accept_prob;
  }

  return List::create(
    Named("samples") = output.samples,
    Named("log_probs") = output.log_probs,
    Named("tree_depths") = tree_depths,
    Named("divergent") = divergent,
    Named("accept_probs") = accept_probs,
    Named("n_divergent") = output.n_divergent,
    Named("mean_accept_rate") = output.mean_accept_rate,
    Named("final_step_size") = output.final_step_size,
    Named("final_mass_matrix") = output.final_mass_matrix
  );
}

//' Get number of parameters for Dirichlet model
//' @noRd
// [[Rcpp::export]]
int cpp_dirichlet_n_params(int n_predictors, int n_components) {
  return n_predictors * n_components + 1;
}

//' Initialize parameters for Dirichlet model
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd cpp_dirichlet_init(int n_predictors, int n_components,
                                    double beta_init = 0.0,
                                    double log_phi_init = 2.0) {
  int n_params = n_predictors * n_components + 1;
  Eigen::VectorXd init(n_params);

  // Beta coefficients
  for (int i = 0; i < n_predictors * n_components; i++) {
    init(i) = beta_init;
  }

  // log(phi)
  init(n_predictors * n_components) = log_phi_init;

  return init;
}

//' Test log Dirichlet probability
//' @noRd
// [[Rcpp::export]]
double cpp_test_log_dirichlet(Eigen::VectorXd y, Eigen::VectorXd alpha) {
  return likelihood::log_dirichlet(y, alpha);
}

//' Test log Dirichlet with softmax
//' @noRd
// [[Rcpp::export]]
double cpp_test_log_dirichlet_softmax(Eigen::VectorXd y, Eigen::VectorXd eta,
                                       double phi) {
  return likelihood::log_dirichlet_softmax(y, eta, phi);
}
