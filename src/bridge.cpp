#include <RcppEigen.h>
#include <simp/simp.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::_;
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

namespace {

Vec as_vec(const Rcpp::NumericVector& x) {
  return Eigen::Map<const Vec>(x.begin(), x.size());
}

}  // namespace

// The built-in scheme family, so R and C++ read the same registry.
// [[Rcpp::export]]
Rcpp::DataFrame simp_schemes_cpp() {
  std::vector<simp::Scheme> reg = simp::registry();
  Rcpp::CharacterVector name(reg.size());
  Rcpp::IntegerVector order(reg.size());
  Rcpp::IntegerVector gevals(reg.size());
  for (std::size_t i = 0; i < reg.size(); ++i) {
    name[i] = reg[i].name;
    order[i] = reg[i].order;
    gevals[i] = reg[i].gradient_evals();
  }
  return Rcpp::DataFrame::create(_["integrator"] = name, _["order"] = order,
                                 _["gradient_evals"] = gevals,
                                 _["stringsAsFactors"] = false);
}

// Integrate a Hamiltonian trajectory. grad_U(q) returns dU/dq; Minv is the
// inverse mass matrix M^{-1}. When trajectory is true the full path is
// returned, otherwise only the endpoint. When U_fn is supplied the total
// energy H = U(q) + 0.5 p' M^{-1} p is recorded at each retained point.
// [[Rcpp::export]]
Rcpp::List simp_integrate_cpp(Rcpp::NumericVector q0, Rcpp::NumericVector p0,
                              Rcpp::Function grad_U, double eps, int n_steps,
                              std::string scheme_name, Rcpp::NumericMatrix Minv,
                              bool trajectory, Rcpp::Nullable<Rcpp::Function> U_fn) {
  const simp::Scheme s = simp::scheme_by_name(scheme_name);
  const Mat MinvE = Rcpp::as<Mat>(Minv);
  const bool want_energy = U_fn.isNotNull();
  Rcpp::Function U = want_energy ? Rcpp::Function(U_fn.get()) : grad_U;

  Vec q = as_vec(q0);
  Vec p = as_vec(p0);
  const int d = q.size();

  // The stepper works in the force convention (force = -dU/dq), so negate the
  // user's potential gradient.
  auto force = [&](const Vec& qq) -> Vec {
    Rcpp::NumericVector g = grad_U(Rcpp::wrap(qq));
    return -as_vec(g);
  };
  auto applyMinv = [&](const Vec& pp) -> Vec { return MinvE * pp; };

  auto energy_at = [&](const Vec& qq, const Vec& pp) -> double {
    double U_val = Rcpp::as<double>(U(Rcpp::wrap(qq)));
    return U_val + 0.5 * pp.dot(MinvE * pp);
  };

  const int n_keep = trajectory ? (n_steps + 1) : 1;
  Mat q_path(d, n_keep);
  Mat p_path(d, n_keep);
  Rcpp::NumericVector energy(want_energy ? n_keep : 0);

  if (trajectory) {
    q_path.col(0) = q;
    p_path.col(0) = p;
    if (want_energy) energy[0] = energy_at(q, p);
    for (int i = 0; i < n_steps; ++i) {
      simp::step(q, p, eps, s, force, applyMinv);
      q_path.col(i + 1) = q;
      p_path.col(i + 1) = p;
      if (want_energy) energy[i + 1] = energy_at(q, p);
    }
  } else {
    simp::integrate(q, p, eps, n_steps, s, force, applyMinv);
    q_path.col(0) = q;
    p_path.col(0) = p;
    if (want_energy) energy[0] = energy_at(q, p);
  }

  return Rcpp::List::create(_["q"] = Rcpp::wrap(q), _["p"] = Rcpp::wrap(p),
                            _["q_path"] = Rcpp::wrap(q_path),
                            _["p_path"] = Rcpp::wrap(p_path),
                            _["energy"] = energy);
}

// Integrate a non-separable Hamiltonian with the generalized (implicit)
// leapfrog. dH_dq(q, p) and dH_dp(q, p) return the partial gradients; H(q, p)
// is optional and records the total energy along the retained points.
// [[Rcpp::export]]
Rcpp::List simp_integrate_nonsep_cpp(Rcpp::NumericVector q0, Rcpp::NumericVector p0,
                                     Rcpp::Function dH_dq, Rcpp::Function dH_dp,
                                     double eps, int n_steps, double tol,
                                     int max_iter, bool trajectory,
                                     Rcpp::Nullable<Rcpp::Function> H_fn) {
  Vec q = as_vec(q0);
  Vec p = as_vec(p0);
  const int d = q.size();
  const bool want_energy = H_fn.isNotNull();
  Rcpp::Function H = want_energy ? Rcpp::Function(H_fn.get()) : dH_dq;

  auto dHdq = [&](const Vec& qq, const Vec& pp) -> Vec {
    Rcpp::NumericVector g = dH_dq(Rcpp::wrap(qq), Rcpp::wrap(pp));
    return as_vec(g);
  };
  auto dHdp = [&](const Vec& qq, const Vec& pp) -> Vec {
    Rcpp::NumericVector g = dH_dp(Rcpp::wrap(qq), Rcpp::wrap(pp));
    return as_vec(g);
  };
  auto energy_at = [&](const Vec& qq, const Vec& pp) -> double {
    return Rcpp::as<double>(H(Rcpp::wrap(qq), Rcpp::wrap(pp)));
  };

  const int n_keep = trajectory ? (n_steps + 1) : 1;
  Mat q_path(d, n_keep);
  Mat p_path(d, n_keep);
  Rcpp::NumericVector energy(want_energy ? n_keep : 0);
  bool all_converged = true;

  if (trajectory) {
    q_path.col(0) = q;
    p_path.col(0) = p;
    if (want_energy) energy[0] = energy_at(q, p);
  }
  for (int i = 0; i < n_steps; ++i) {
    simp::ImplicitResult res =
        simp::step_implicit(q, p, eps, dHdq, dHdp, tol, max_iter);
    all_converged = all_converged && res.converged;
    if (trajectory) {
      q_path.col(i + 1) = q;
      p_path.col(i + 1) = p;
      if (want_energy) energy[i + 1] = energy_at(q, p);
    }
  }
  if (!trajectory) {
    q_path.col(0) = q;
    p_path.col(0) = p;
    if (want_energy) energy[0] = energy_at(q, p);
  }

  return Rcpp::List::create(_["q"] = Rcpp::wrap(q), _["p"] = Rcpp::wrap(p),
                            _["q_path"] = Rcpp::wrap(q_path),
                            _["p_path"] = Rcpp::wrap(p_path),
                            _["energy"] = energy,
                            _["converged"] = all_converged);
}
