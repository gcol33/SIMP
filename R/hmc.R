#' Hamiltonian Monte Carlo with a chosen symplectic integrator
#'
#' A minimal Hamiltonian Monte Carlo sampler that draws each proposal by
#' integrating the Hamiltonian flow of the target with [simp_integrate()] and
#' accepting it by the Metropolis rule on the total energy. The integrator is a
#' free choice, so this is the smallest complete setting in which the schemes
#' can be compared on a real sampling problem.
#'
#' @param U Function mapping a position to the potential energy, U = -log
#'   density of the target up to an additive constant.
#' @param grad_U Function mapping a position to the gradient of U.
#' @param init Numeric vector, the starting position.
#' @param n_iter Number of retained iterations.
#' @param eps Positive leapfrog step size.
#' @param n_steps Number of integrator steps per proposal.
#' @param integrator Name of a built-in scheme; see [simp_schemes()].
#' @param mass Mass matrix M, as in [simp_integrate()]. Momenta are drawn from
#'   N(0, M).
#' @param warmup Number of initial iterations to discard.
#' @param divergence_threshold Energy-error magnitude above which a proposal is
#'   counted as divergent and rejected.
#'
#' @return A list with the `draws` matrix (one row per retained iteration, one
#'   column per coordinate), the `accept_rate`, the number of divergent
#'   proposals `n_divergent`, and the `integrator` used.
#'
#' @examples
#' # Sample a correlated 2-D Gaussian.
#' Sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
#' P <- solve(Sigma)
#' U <- function(q) 0.5 * as.numeric(t(q) %*% P %*% q)
#' grad_U <- function(q) as.numeric(P %*% q)
#' fit <- simp_hmc(U, grad_U, init = c(0, 0), n_iter = 200, eps = 0.2,
#'                 n_steps = 20, warmup = 100)
#' colMeans(fit$draws)
#'
#' @export
simp_hmc <- function(U, grad_U, init, n_iter, eps, n_steps,
                     integrator = "leapfrog", mass = NULL,
                     warmup = 0L, divergence_threshold = 1000) {
  if (!is.function(U) || !is.function(grad_U)) {
    stop("U and grad_U must be functions", call. = FALSE)
  }
  q <- as.numeric(init)
  d <- length(q)
  n_iter <- as.integer(n_iter)
  warmup <- as.integer(warmup)
  if (is.na(n_iter) || n_iter < 1) {
    stop("n_iter must be a positive integer", call. = FALSE)
  }
  Minv <- .simp_minv(mass, d)
  M <- solve(Minv)
  M_chol <- chol(M)  # upper triangular, M = t(M_chol) %*% M_chol

  kinetic <- function(p) 0.5 * sum(p * (Minv %*% p))

  draws <- matrix(NA_real_, nrow = n_iter, ncol = d)
  n_accept <- 0L
  n_divergent <- 0L
  total <- warmup + n_iter

  for (it in seq_len(total)) {
    p0 <- as.numeric(crossprod(M_chol, stats::rnorm(d)))  # ~ N(0, M)
    H0 <- U(q) + kinetic(p0)
    prop <- simp_integrate_cpp(q, p0, grad_U, eps, n_steps, integrator,
                               Minv, FALSE, NULL)
    q1 <- as.numeric(prop$q)
    p1 <- as.numeric(prop$p)
    dH <- U(q1) + kinetic(p1) - H0

    if (!is.finite(dH) || abs(dH) > divergence_threshold) {
      n_divergent <- n_divergent + 1L
    } else if (log(stats::runif(1)) < -dH) {
      q <- q1
      n_accept <- n_accept + 1L
    }
    if (it > warmup) {
      draws[it - warmup, ] <- q
    }
  }

  list(draws = draws, accept_rate = n_accept / total,
       n_divergent = n_divergent, integrator = integrator)
}
