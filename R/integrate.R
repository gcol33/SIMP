# Inverse mass matrix M^{-1} from the user's `mass` argument.
# NULL -> identity, a vector -> diagonal mass, a matrix -> dense mass.
.simp_minv <- function(mass, d) {
  if (is.null(mass)) {
    return(diag(d))
  }
  if (is.matrix(mass)) {
    if (nrow(mass) != d || ncol(mass) != d) {
      stop("mass must be a ", d, " by ", d, " matrix", call. = FALSE)
    }
    return(solve(mass))
  }
  mass <- as.numeric(mass)
  if (length(mass) != d) {
    stop("mass must have length ", d, call. = FALSE)
  }
  if (any(mass <= 0)) {
    stop("mass entries must be positive", call. = FALSE)
  }
  diag(1 / mass, nrow = d)
}

# Op sequence (kick/drift kinds and weights) of the requested integrator, built
# once in C++ so every scheme -- built-in, or step-adapted -- reaches the shared
# stepper the same way. `"adaptive2"` / `"adaptive3"` require `nu_max`.
.simp_ops <- function(integrator, nu_max) {
  if (integrator %in% c("adaptive2", "adaptive3")) {
    if (is.null(nu_max)) {
      stop("nu_max must be supplied for the adaptive integrators", call. = FALSE)
    }
    if (!is.numeric(nu_max) || length(nu_max) != 1 || nu_max <= 0) {
      stop("nu_max must be a single positive number", call. = FALSE)
    }
    stages <- if (integrator == "adaptive3") 3L else 2L
    return(simp_adaptive_ops_cpp(as.numeric(nu_max), stages))
  }
  simp_scheme_ops_cpp(integrator)  # errors "unknown integrator" for bad names
}

#' Integrate a Hamiltonian trajectory
#'
#' Advance a separable Hamiltonian system H(q, p) = U(q) + K(p), with kinetic
#' energy K(p) = 0.5 p' Minv p where Minv is the inverse mass matrix, using a
#' symplectic splitting integrator. The integrator conserves a nearby (shadow)
#' energy, so the total energy stays bounded over long trajectories instead of
#' drifting.
#'
#' @param q0 Numeric vector, the starting position.
#' @param p0 Numeric vector, the starting momentum, same length as `q0`.
#' @param grad_U Function mapping a position vector to the gradient of the
#'   potential energy U, that is dU/dq. For a target density with log-density
#'   `lp`, the potential is U = -lp and its gradient is -d(lp)/dq.
#' @param eps Positive step size.
#' @param n_steps Number of integrator steps.
#' @param integrator Name of a built-in scheme (see [simp_schemes()]), or
#'   `"adaptive2"` / `"adaptive3"` for the step-adapted minimum-error two- or
#'   three-stage scheme, which additionally needs `nu_max`. Defaults to
#'   `"leapfrog"`.
#' @param mass Mass matrix M. `NULL` (the default) uses the identity; a numeric
#'   vector uses a diagonal mass; a matrix uses a dense mass.
#' @param nu_max Largest dimensionless step `omega_max * eps` (with `omega_max`
#'   the highest target frequency) the adapted integrators tune for. Required
#'   when `integrator` is `"adaptive2"` or `"adaptive3"`, ignored otherwise.
#' @param trajectory If `TRUE`, return the position and momentum at every step.
#'   If `FALSE` (the default), return only the endpoint.
#' @param U Optional function mapping a position to the potential energy U. When
#'   supplied, the total energy is recorded along the returned points.
#'
#' @return A list with the endpoint `q` and `p`. When `trajectory = TRUE` it
#'   also holds `q_path` and `p_path`, each a matrix with one row per retained
#'   point and one column per coordinate. When `U` is supplied it holds
#'   `energy`, the total energy at the retained points.
#'
#' @examples
#' # A unit isotropic Gaussian target: U(q) = 0.5 q'q, grad_U(q) = q.
#' grad_U <- function(q) q
#' out <- simp_integrate(q0 = c(1, 0), p0 = c(0, 1), grad_U = grad_U,
#'                       eps = 0.1, n_steps = 60, integrator = "yoshida4")
#' out$q
#'
#' @export
simp_integrate <- function(q0, p0, grad_U, eps, n_steps,
                           integrator = "leapfrog", mass = NULL,
                           nu_max = NULL, trajectory = FALSE, U = NULL) {
  q0 <- as.numeric(q0)
  p0 <- as.numeric(p0)
  d <- length(q0)
  if (length(p0) != d) {
    stop("q0 and p0 must have the same length", call. = FALSE)
  }
  if (!is.function(grad_U)) {
    stop("grad_U must be a function", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("eps must be a single positive number", call. = FALSE)
  }
  n_steps <- as.integer(n_steps)
  if (is.na(n_steps) || n_steps < 1) {
    stop("n_steps must be a positive integer", call. = FALSE)
  }
  Minv <- .simp_minv(mass, d)
  ops <- .simp_ops(integrator, nu_max)

  res <- simp_integrate_ops_cpp(q0, p0, grad_U, eps, n_steps, ops$kind, ops$coef,
                                Minv, isTRUE(trajectory), U)

  out <- list(q = as.numeric(res$q), p = as.numeric(res$p))
  if (isTRUE(trajectory)) {
    out$q_path <- t(res$q_path)
    out$p_path <- t(res$p_path)
  }
  if (!is.null(U)) {
    out$energy <- as.numeric(res$energy)
  }
  out
}

#' Integrate a split-potential Hamiltonian (multiple time-stepping)
#'
#' Advance a separable Hamiltonian whose potential splits as
#' U = U_slow + U_fast, using the RESPA / Verlet-I multiple-time-stepping
#' integrator. The slow force is evaluated once per outer step while the fast
#' force drives `m` inner leapfrog substeps, so a stiff-but-cheap part (a
#' Gaussian latent prior, say) can take the small steps it needs without paying
#' for the expensive-but-smooth part (a likelihood) at the same rate. The step
#' is symmetric, time reversible, and symplectic.
#'
#' @param q0 Numeric vector, the starting position.
#' @param p0 Numeric vector, the starting momentum, same length as `q0`.
#' @param grad_slow,grad_fast Functions mapping a position to the gradient of
#'   the slow and fast potential parts. Their sum is dU/dq.
#' @param eps Positive outer step size.
#' @param n_steps Number of outer steps.
#' @param m Number of inner fast substeps per outer step.
#' @param mass Mass matrix M, as in [simp_integrate()].
#' @param trajectory If `TRUE`, return the position and momentum at every outer
#'   step. If `FALSE` (the default), return only the endpoint.
#' @param U Optional function mapping a position to the total potential energy
#'   U_slow + U_fast. When supplied, the total energy is recorded along the
#'   returned points.
#'
#' @return A list with the endpoint `q` and `p`, and, when requested, `q_path`,
#'   `p_path`, and `energy`, as in [simp_integrate()].
#'
#' @examples
#' # Stiff Gaussian prior (fast) plus a soft quadratic likelihood (slow).
#' grad_fast <- function(q) 100 * q
#' grad_slow <- function(q) 0.1 * q
#' out <- simp_integrate_mts(q0 = c(1, 0), p0 = c(0, 1),
#'                           grad_slow = grad_slow, grad_fast = grad_fast,
#'                           eps = 0.2, n_steps = 30, m = 10)
#' out$q
#'
#' @export
simp_integrate_mts <- function(q0, p0, grad_slow, grad_fast, eps, n_steps, m,
                               mass = NULL, trajectory = FALSE, U = NULL) {
  q0 <- as.numeric(q0)
  p0 <- as.numeric(p0)
  d <- length(q0)
  if (length(p0) != d) {
    stop("q0 and p0 must have the same length", call. = FALSE)
  }
  if (!is.function(grad_slow) || !is.function(grad_fast)) {
    stop("grad_slow and grad_fast must be functions", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("eps must be a single positive number", call. = FALSE)
  }
  n_steps <- as.integer(n_steps)
  if (is.na(n_steps) || n_steps < 1) {
    stop("n_steps must be a positive integer", call. = FALSE)
  }
  m <- as.integer(m)
  if (is.na(m) || m < 1) {
    stop("m must be a positive integer", call. = FALSE)
  }
  Minv <- .simp_minv(mass, d)

  res <- simp_integrate_mts_cpp(q0, p0, grad_slow, grad_fast, eps, n_steps, m,
                                Minv, isTRUE(trajectory), U)

  out <- list(q = as.numeric(res$q), p = as.numeric(res$p))
  if (isTRUE(trajectory)) {
    out$q_path <- t(res$q_path)
    out$p_path <- t(res$p_path)
  }
  if (!is.null(U)) {
    out$energy <- as.numeric(res$energy)
  }
  out
}

#' Integrate a non-separable Hamiltonian (generalized leapfrog)
#'
#' Advance a Hamiltonian H(q, p) whose kinetic and potential parts do not
#' separate, as in Riemannian dynamics with a position-dependent mass, using
#' the generalized (implicit) leapfrog. The step is symmetric, time reversible,
#' and symplectic; its two implicit stages are solved by fixed-point iteration.
#' For a separable H = U(q) + K(p) it reduces to [simp_integrate()] with the
#' leapfrog scheme.
#'
#' @param q0 Numeric vector, the starting position.
#' @param p0 Numeric vector, the starting momentum, same length as `q0`.
#' @param dH_dq Function of `(q, p)` returning the partial gradient dH/dq.
#' @param dH_dp Function of `(q, p)` returning the partial gradient dH/dp.
#' @param eps Positive step size.
#' @param n_steps Number of integrator steps.
#' @param tol Fixed-point tolerance for the implicit stages.
#' @param max_iter Maximum fixed-point iterations per stage.
#' @param trajectory If `TRUE`, return the position and momentum at every step.
#' @param H Optional function of `(q, p)` returning the total energy, recorded
#'   along the returned points.
#'
#' @return A list with the endpoint `q` and `p`, a `converged` flag reporting
#'   whether every implicit stage met `tol`, and, when requested, `q_path`,
#'   `p_path`, and `energy`.
#'
#' @examples
#' # Separable check: H = 0.5 p'p + 0.5 q'q reduces to leapfrog.
#' out <- simp_integrate_nonsep(
#'   q0 = c(1, 0), p0 = c(0, 1),
#'   dH_dq = function(q, p) q, dH_dp = function(q, p) p,
#'   eps = 0.1, n_steps = 60)
#' out$q
#'
#' @export
simp_integrate_nonsep <- function(q0, p0, dH_dq, dH_dp, eps, n_steps,
                                  tol = 1e-10, max_iter = 100L,
                                  trajectory = FALSE, H = NULL) {
  q0 <- as.numeric(q0)
  p0 <- as.numeric(p0)
  if (length(p0) != length(q0)) {
    stop("q0 and p0 must have the same length", call. = FALSE)
  }
  if (!is.function(dH_dq) || !is.function(dH_dp)) {
    stop("dH_dq and dH_dp must be functions of (q, p)", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("eps must be a single positive number", call. = FALSE)
  }
  n_steps <- as.integer(n_steps)
  if (is.na(n_steps) || n_steps < 1) {
    stop("n_steps must be a positive integer", call. = FALSE)
  }

  res <- simp_integrate_nonsep_cpp(q0, p0, dH_dq, dH_dp, eps, n_steps,
                                   as.numeric(tol), as.integer(max_iter),
                                   isTRUE(trajectory), H)

  out <- list(q = as.numeric(res$q), p = as.numeric(res$p),
              converged = isTRUE(res$converged))
  if (isTRUE(trajectory)) {
    out$q_path <- t(res$q_path)
    out$p_path <- t(res$p_path)
  }
  if (!is.null(H)) {
    out$energy <- as.numeric(res$energy)
  }
  out
}

#' Available symplectic integrators
#'
#' The built-in scheme family: leapfrog, the two-stage minimum-error member, and
#' the Yoshida members generated from leapfrog by triple-jump composition. The
#' step-adapted `"adaptive2"` / `"adaptive3"` schemes are not listed here because
#' their coefficients depend on the step band; see [simp_b_optimal()].
#'
#' @return A data frame with one row per integrator giving its name, order of
#'   accuracy, and number of gradient evaluations per step.
#'
#' @examples
#' simp_schemes()
#'
#' @export
simp_schemes <- function() {
  simp_schemes_cpp()
}

#' Minimum-error multistage coefficients for a step band
#'
#' The multistage kick and drift weights that minimize the worst-case energy
#' error a splitting integrator induces on a harmonic (Gaussian) target across
#' the band of dimensionless steps `(0, nu_max]`. A separable quadratic target
#' excites oscillator modes at `nu_k = omega_k * eps`, so with `nu_max =
#' omega_max * eps` (the largest such step, from the highest target frequency)
#' these coefficients give the minimum-error integrator for the problem at hand
#' rather than a fixed compromise. Feed the coefficients back through
#' [simp_integrate()] with `integrator = "adaptive2"` or `"adaptive3"`, which
#' recompute them internally from the same `nu_max`.
#'
#' @param nu_max Largest dimensionless step `omega_max * eps`, a single positive
#'   number.
#' @param stages `2` for the two-stage scheme (returns the kick weight `b`) or
#'   `3` for the three-stage scheme (returns `c(b, a)`, the kick and drift free
#'   weights).
#'
#' @return A numeric vector of the free coefficients: length 1 for two stages,
#'   length 2 for three.
#'
#' @examples
#' simp_b_optimal(nu_max = 1.5, stages = 2)
#'
#' @export
simp_b_optimal <- function(nu_max, stages = 2) {
  if (!is.numeric(nu_max) || length(nu_max) != 1 || nu_max <= 0) {
    stop("nu_max must be a single positive number", call. = FALSE)
  }
  stages <- as.integer(stages)
  if (is.na(stages) || !(stages %in% c(2L, 3L))) {
    stop("stages must be 2 or 3", call. = FALSE)
  }
  as.numeric(simp_b_optimal_cpp(as.numeric(nu_max), stages))
}
