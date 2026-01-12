# S3 class for SIMP fit objects
# Contains posterior samples and diagnostics

#' Create a new simp_fit object
#' @noRd
new_simp_fit <- function(draws, log_probs, diagnostics, model, timing) {
  structure(
    list(
      draws = draws,
      log_probs = log_probs,
      diagnostics = diagnostics,
      model = model,
      timing = timing
    ),
    class = "simp_fit"
  )
}

#' @export
print.simp_fit <- function(x, ...) {
  cli::cli_h3("SIMP Model Fit")
  cli::cli_rule()

  cat("Family:", x$model$family$family, "\n")
  cat("Observations:", x$model$response$n_obs, "\n")
  cat("Components:", x$model$response$n_components, "\n")

  cat("\nSampling:\n")
  cat(sprintf("  Chains: %d\n", x$diagnostics$n_chains))
  cat(sprintf("  Iterations: %d (warmup: %d)\n",
              x$diagnostics$n_iter, x$diagnostics$n_warmup))
  cat(sprintf("  Time: %.1f seconds\n", x$timing))

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Mean accept rate: %.2f\n", x$diagnostics$mean_accept_rate))

  if (x$diagnostics$n_divergent > 0) {
    cli::cli_alert_warning("Model has divergent transitions - interpret with caution")
  }

  cat("\nUse summary() for parameter estimates\n")

  invisible(x)
}

#' Summarize SIMP fit
#'
#' @param object A simp_fit object
#' @param ... Additional arguments (ignored)
#' @return A summary object
#' @export
summary.simp_fit <- function(object, ...) {
  # Compute summary statistics for each parameter
  n_params <- dim(object$draws)[3]
  param_names <- dimnames(object$draws)[[3]]

  # Combine chains
  draws_combined <- matrix(
    object$draws,
    nrow = prod(dim(object$draws)[1:2]),
    ncol = n_params
  )

  # Compute statistics
  stats <- data.frame(
    parameter = param_names,
    mean = apply(draws_combined, 2, mean),
    sd = apply(draws_combined, 2, sd),
    q2.5 = apply(draws_combined, 2, quantile, 0.025),
    q50 = apply(draws_combined, 2, quantile, 0.5),
    q97.5 = apply(draws_combined, 2, quantile, 0.975),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Compute Rhat (simplified split-R-hat)
  stats$rhat <- sapply(seq_len(n_params), function(i) {
    draws_matrix <- object$draws[, , i, drop = FALSE]
    dim(draws_matrix) <- dim(object$draws)[1:2]
    compute_rhat(draws_matrix)
  })

  # Compute ESS (simplified)
  stats$ess_bulk <- sapply(seq_len(n_params), function(i) {
    compute_ess(draws_combined[, i])
  })

  structure(
    list(
      stats = stats,
      diagnostics = object$diagnostics,
      model = object$model
    ),
    class = "summary.simp_fit"
  )
}

#' @export
print.summary.simp_fit <- function(x, digits = 2, ...) {
  cli::cli_h3("SIMP Model Summary")
  cli::cli_rule()

  cat("Family:", x$model$family$family, "\n")
  cat("Observations:", x$model$response$n_obs, "\n\n")

  # Print statistics
  stats <- x$stats
  stats$mean <- round(stats$mean, digits)
  stats$sd <- round(stats$sd, digits)
  stats$q2.5 <- round(stats$q2.5, digits)
  stats$q50 <- round(stats$q50, digits)
  stats$q97.5 <- round(stats$q97.5, digits)
  stats$rhat <- round(stats$rhat, 3)
  stats$ess_bulk <- round(stats$ess_bulk, 0)

  print(stats, row.names = FALSE)

  # Warnings
  if (any(stats$rhat > 1.01, na.rm = TRUE)) {
    cli::cli_alert_warning("Some Rhat values > 1.01 - chains may not have converged")
  }

  if (any(stats$ess_bulk < 100, na.rm = TRUE)) {
    cli::cli_alert_warning("Some ESS values < 100 - consider running longer chains")
  }

  if (x$diagnostics$n_divergent > 0) {
    cli::cli_alert_warning("{x$diagnostics$n_divergent} divergent transitions")
  }

  invisible(x)
}

#' Extract posterior draws
#'
#' @param x A simp_fit object
#' @param ... Additional arguments (ignored)
#' @return A data frame of posterior draws
#' @export
as.data.frame.simp_fit <- function(x, ...) {
  n_iter <- dim(x$draws)[1]
  n_chains <- dim(x$draws)[2]
  n_params <- dim(x$draws)[3]
  param_names <- dimnames(x$draws)[[3]]

  # Reshape to long format
  df <- data.frame(
    .chain = rep(rep(seq_len(n_chains), each = n_iter), n_params),
    .iteration = rep(rep(seq_len(n_iter), n_chains), n_params),
    .draw = rep(seq_len(n_iter * n_chains), n_params)
  )

  # Add parameter columns
  draws_matrix <- matrix(x$draws, nrow = n_iter * n_chains, ncol = n_params)
  colnames(draws_matrix) <- param_names

  cbind(df, as.data.frame(draws_matrix))
}

#' Extract draws as array
#'
#' @param x A simp_fit object
#' @param ... Additional arguments (ignored)
#' @return A 3D array: iterations x chains x parameters
#' @export
as.array.simp_fit <- function(x, ...) {
  x$draws
}

# -----------------------------------------------------------------------------
# Diagnostic functions
# -----------------------------------------------------------------------------

#' Compute split R-hat
#' @noRd
compute_rhat <- function(draws_matrix) {
  # draws_matrix: iterations x chains
  n_iter <- nrow(draws_matrix)
  n_chains <- ncol(draws_matrix)

  if (n_chains < 2) return(NA_real_)

  # Split each chain in half
  n_half <- floor(n_iter / 2)
  if (n_half < 2) return(NA_real_)

  split_draws <- matrix(nrow = n_half, ncol = n_chains * 2)
  for (i in seq_len(n_chains)) {
    split_draws[, 2 * i - 1] <- draws_matrix[1:n_half, i]
    split_draws[, 2 * i] <- draws_matrix[(n_half + 1):(2 * n_half), i]
  }

  # Between-chain and within-chain variance
  chain_means <- colMeans(split_draws)
  overall_mean <- mean(chain_means)

  B <- n_half * sum((chain_means - overall_mean)^2) / (2 * n_chains - 1)
  W <- mean(apply(split_draws, 2, var))

  # R-hat
  var_hat <- (1 - 1/n_half) * W + B / n_half
  sqrt(var_hat / W)
}

#' Compute effective sample size (simplified)
#' @noRd
compute_ess <- function(x) {
  n <- length(x)
  if (n < 4) return(NA_real_)

  # Check for constant values
  if (var(x) < .Machine$double.eps) return(NA_real_)

  # Simple autocorrelation-based ESS
  acf_vals <- acf(x, lag.max = min(n - 1, 500), plot = FALSE)$acf[-1]

  # Sum positive autocorrelations (handle NA values)
  rho_sum <- 0
  for (k in seq_along(acf_vals)) {
    if (is.na(acf_vals[k]) || acf_vals[k] < 0) break
    rho_sum <- rho_sum + acf_vals[k]
  }

  n / (1 + 2 * rho_sum)
}

# -----------------------------------------------------------------------------
# Accessor functions
# -----------------------------------------------------------------------------

#' Check if object is a simp_fit
#'
#' @param x Object to check
#' @return Logical
#' @export
is_simp_fit <- function(x) {
  inherits(x, "simp_fit")
}

#' Extract log probability from fit
#'
#' @param fit A simp_fit object
#' @return Matrix of log probabilities (iterations x chains)
#' @export
log_prob <- function(fit) {
  if (!inherits(fit, "simp_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls simp_fit} object")
  }
  fit$log_probs
}

#' Get number of divergent transitions
#'
#' @param fit A simp_fit object
#' @return Integer count
#' @export
n_divergent <- function(fit) {
  if (!inherits(fit, "simp_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls simp_fit} object")
  }
  fit$diagnostics$n_divergent
}
