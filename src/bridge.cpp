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

// Reconstruct a scheme from an op sequence (kind: 1 = kick, 0 = drift; coef the
// weight of each op). Name and order are irrelevant to stepping, so any scheme
// -- built-in, parametric, or step-adapted -- reaches the same stepper this way.
simp::Scheme scheme_from_ops(const Rcpp::IntegerVector& kind,
                             const Rcpp::NumericVector& coef) {
  simp::Scheme s;
  s.name = "custom";
  s.order = 0;
  const int n = kind.size();
  s.ops.reserve(n);
  for (int i = 0; i < n; ++i) {
    s.ops.emplace_back(kind[i] == 1 ? simp::Op::Kick : simp::Op::Drift, coef[i]);
  }
  return s;
}

Rcpp::List ops_of(const simp::Scheme& s) {
  const int n = static_cast<int>(s.ops.size());
  Rcpp::IntegerVector kind(n);
  Rcpp::NumericVector coef(n);
  for (int i = 0; i < n; ++i) {
    kind[i] = (s.ops[i].first == simp::Op::Kick) ? 1 : 0;
    coef[i] = s.ops[i].second;
  }
  return Rcpp::List::create(_["kind"] = kind, _["coef"] = coef,
                            _["order"] = s.order,
                            _["gradient_evals"] = s.gradient_evals());
}

// Shared separable-Hamiltonian trajectory runner: the single body behind every
// named, parametric, and step-adapted scheme.
Rcpp::List run_separable(Vec q, Vec p, const Rcpp::Function& grad_U, double eps,
                         int n_steps, const simp::Scheme& s, const Mat& MinvE,
                         bool trajectory, Rcpp::Nullable<Rcpp::Function> U_fn) {
  const int d = q.size();
  const bool want_energy = U_fn.isNotNull();
  Rcpp::Function U = want_energy ? Rcpp::Function(U_fn.get()) : grad_U;

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

// Op sequence of a named built-in scheme.
// [[Rcpp::export]]
Rcpp::List simp_scheme_ops_cpp(std::string name) {
  return ops_of(simp::scheme_by_name(name));
}

// Op sequence of the step-adapted minimum-error scheme for a band of
// dimensionless steps up to nu_max. stages is 2 or 3.
// [[Rcpp::export]]
Rcpp::List simp_adaptive_ops_cpp(double nu_max, int stages) {
  simp::Scheme s = (stages == 3) ? simp::three_stage_adaptive(nu_max)
                                 : simp::two_stage_adaptive(nu_max);
  return ops_of(s);
}

// Minimum-error multistage coefficients for a band of dimensionless steps up to
// nu_max: a length-1 vector (b) for two stages, a length-2 vector (b, a) for
// three.
// [[Rcpp::export]]
Rcpp::NumericVector simp_b_optimal_cpp(double nu_max, int stages) {
  if (stages == 3) {
    double b = 0.25, a = 0.25;
    simp::three_stage_optimal(nu_max, b, a);
    return Rcpp::NumericVector::create(b, a);
  }
  return Rcpp::NumericVector::create(simp::two_stage_optimal_b(nu_max));
}

// Relative harmonic energy-error amplitude of an op sequence at dimensionless
// step nu (see simp::harmonic_energy_error).
// [[Rcpp::export]]
double simp_harmonic_error_cpp(Rcpp::IntegerVector kind, Rcpp::NumericVector coef,
                               double nu) {
  return simp::harmonic_energy_error(scheme_from_ops(kind, coef), nu);
}

// Largest stable dimensionless step of an op sequence (see
// simp::harmonic_stability_limit).
// [[Rcpp::export]]
double simp_stability_limit_cpp(Rcpp::IntegerVector kind,
                                Rcpp::NumericVector coef, double nu_hi) {
  return simp::harmonic_stability_limit(scheme_from_ops(kind, coef), nu_hi);
}

// Integrate a Hamiltonian trajectory with a scheme given as an op sequence.
// grad_U(q) returns dU/dq; Minv is the inverse mass M^{-1}. When trajectory is
// true the full path is returned, otherwise only the endpoint. When U_fn is
// supplied the total energy H = U(q) + 0.5 p' M^{-1} p is recorded at each
// retained point.
// [[Rcpp::export]]
Rcpp::List simp_integrate_ops_cpp(Rcpp::NumericVector q0, Rcpp::NumericVector p0,
                                  Rcpp::Function grad_U, double eps, int n_steps,
                                  Rcpp::IntegerVector kind,
                                  Rcpp::NumericVector coef,
                                  Rcpp::NumericMatrix Minv, bool trajectory,
                                  Rcpp::Nullable<Rcpp::Function> U_fn) {
  const simp::Scheme s = scheme_from_ops(kind, coef);
  const Mat MinvE = Rcpp::as<Mat>(Minv);
  return run_separable(as_vec(q0), as_vec(p0), grad_U, eps, n_steps, s, MinvE,
                       trajectory, U_fn);
}

// Integrate a split-potential Hamiltonian with the multiple-time-stepping
// (RESPA / Verlet-I) integrator. grad_slow and grad_fast each return the
// gradient dU/dq of their potential part; the fast part is inner-looped m times
// per outer step. U_fn, when supplied, returns the total potential U_slow +
// U_fast for the recorded energy.
// [[Rcpp::export]]
Rcpp::List simp_integrate_mts_cpp(Rcpp::NumericVector q0, Rcpp::NumericVector p0,
                                  Rcpp::Function grad_slow,
                                  Rcpp::Function grad_fast, double eps,
                                  int n_steps, int m, Rcpp::NumericMatrix Minv,
                                  bool trajectory,
                                  Rcpp::Nullable<Rcpp::Function> U_fn) {
  const Mat MinvE = Rcpp::as<Mat>(Minv);
  const bool want_energy = U_fn.isNotNull();
  Rcpp::Function U = want_energy ? Rcpp::Function(U_fn.get()) : grad_slow;

  Vec q = as_vec(q0);
  Vec p = as_vec(p0);
  const int d = q.size();

  auto force_slow = [&](const Vec& qq) -> Vec {
    Rcpp::NumericVector g = grad_slow(Rcpp::wrap(qq));
    return -as_vec(g);
  };
  auto force_fast = [&](const Vec& qq) -> Vec {
    Rcpp::NumericVector g = grad_fast(Rcpp::wrap(qq));
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
      simp::step_mts(q, p, eps, m, force_slow, force_fast, applyMinv);
      q_path.col(i + 1) = q;
      p_path.col(i + 1) = p;
      if (want_energy) energy[i + 1] = energy_at(q, p);
    }
  } else {
    simp::integrate_mts(q, p, eps, n_steps, m, force_slow, force_fast, applyMinv);
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
