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
#' @param integrator Name of a built-in scheme; see [simp_schemes()]. Defaults
#'   to `"leapfrog"`.
#' @param mass Mass matrix M. `NULL` (the default) uses the identity; a numeric
#'   vector uses a diagonal mass; a matrix uses a dense mass.
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
                           trajectory = FALSE, U = NULL) {
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

  res <- simp_integrate_cpp(q0, p0, grad_U, eps, n_steps, integrator, Minv,
                            isTRUE(trajectory), U)

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

#' Available symplectic integrators
#'
#' The built-in scheme family: leapfrog and the Yoshida members generated from
#' it by triple-jump composition.
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
