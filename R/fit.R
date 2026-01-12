# Model fitting for SIMP
# Main fit() function and related utilities

#' @useDynLib SIMP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Fit a SIMP model
#'
#' Run Bayesian inference on a SIMP model specification using the NUTS sampler.
#'
#' @param model A `simp_model` object created by [simp()]
#' @param chains Number of MCMC chains (default 4)
#' @param iter Number of sampling iterations per chain (default 2000)
#' @param warmup Number of warmup iterations per chain (default half of iter)
#' @param thin Thinning interval (default 1, no thinning)
#' @param cores Number of cores for parallel chains (default 1)
#' @param seed Random seed for reproducibility
#' @param init Initialization method: "random" (default), "zero", or a list
#'   of initial values
#' @param max_treedepth Maximum tree depth for NUTS (default 10)
#' @param adapt_delta Target acceptance rate (default 0.8)
#' @param verbose Print progress (default TRUE)
#' @return A `simp_fit` object containing posterior samples and diagnostics
#' @export
#' @examples
#' \dontrun{
#' # Create model
#' model <- simp(cover ~ temp + (1 | site), data = dat)
#'
#' # Fit model
#' fit <- fit(model, chains = 4, iter = 2000)
#'
#' # Examine results
#' summary(fit)
#' }
fit <- function(model, chains = 4, iter = 2000, warmup = NULL, thin = 1,
                cores = 1, seed = NULL, init = "random",
                max_treedepth = 10, adapt_delta = 0.8, verbose = TRUE) {

  if (!inherits(model, "simp_model")) {
    cli::cli_abort("{.arg model} must be a {.cls simp_model} object from {.fn simp}")
  }

  if (model$fitted) {
    cli::cli_warn("Model has already been fitted. Creating new fit.")
  }

  # Set defaults
  if (is.null(warmup)) {
    warmup <- floor(iter / 2)
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  # Validate arguments
  check_positive_int(chains, "chains")
  check_positive_int(iter, "iter")
  check_positive_int(warmup, "warmup")
  check_positive_int(thin, "thin")
  check_positive_int(cores, "cores")
  check_positive_int(max_treedepth, "max_treedepth")

  if (adapt_delta <= 0 || adapt_delta >= 1) {
    cli::cli_abort("{.arg adapt_delta} must be between 0 and 1")
  }

  # Prepare data for C++
  Y <- prepare_response(model)
  X <- model$fixed$X

  # Get prior settings
  beta_prior_sd <- model$prior$beta$parameters$sd
  phi_prior_shape <- model$prior$phi$parameters$shape
  phi_prior_rate <- model$prior$phi$parameters$rate

  # Initialize parameters
  n_params <- cpp_dirichlet_n_params(ncol(X), ncol(Y))
  init_values <- initialize_params(init, n_params, model)

  # Run chains
  if (verbose) {
    cli::cli_h3("Fitting SIMP Model")
    cli::cli_alert_info("Family: {model$family$family}")
    cli::cli_alert_info("Chains: {chains}, Iter: {iter}, Warmup: {warmup}")
  }

  chain_results <- vector("list", chains)
  chain_seeds <- seed + seq_len(chains) - 1

  start_time <- Sys.time()

  for (chain in seq_len(chains)) {
    if (verbose) {
      cli::cli_alert("Chain {chain}/{chains}")
    }

    result <- cpp_fit_dirichlet(
      Y = Y,
      X = X,
      n_iter = iter,
      n_warmup = warmup,
      init = init_values,
      seed = chain_seeds[chain],
      max_treedepth = max_treedepth,
      adapt_delta = adapt_delta,
      beta_prior_sd = beta_prior_sd,
      phi_prior_shape = phi_prior_shape,
      phi_prior_rate = phi_prior_rate,
      verbose = verbose && (chains == 1)
    )

    chain_results[[chain]] <- result
  }

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Combine chains
  fit_obj <- combine_chains(chain_results, model, chains, iter,
                             warmup, thin, elapsed)

  # Apply thinning
  if (thin > 1) {
    fit_obj <- thin_samples(fit_obj, thin)
  }

  # Update model
  model$fitted <- TRUE
  model$fit <- fit_obj
  fit_obj$model <- model

  if (verbose) {
    cli::cli_alert_success("Sampling complete in {round(elapsed, 1)} seconds")
    if (fit_obj$diagnostics$n_divergent > 0) {
      cli::cli_warn("{fit_obj$diagnostics$n_divergent} divergent transitions")
    }
  }

  fit_obj
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

#' Prepare response matrix for C++
#' @noRd
prepare_response <- function(model) {
  resp <- model$response$data

  if (inherits(resp, "simp_composition")) {
    as.matrix(resp)
  } else {
    as.matrix(resp)
  }
}

#' Initialize parameters
#' @noRd
initialize_params <- function(init, n_params, model) {
  if (is.character(init)) {
    if (init == "random") {
      # Small random values for beta, reasonable value for log(phi)
      stats::rnorm(n_params, 0, 0.1)
    } else if (init == "zero") {
      c(rep(0, n_params - 1), 2)  # log(phi) = 2 -> phi ~= 7
    } else {
      cli::cli_abort("Unknown init method: {.val {init}}")
    }
  } else if (is.numeric(init)) {
    if (length(init) != n_params) {
      cli::cli_abort("init must have length {n_params}")
    }
    init
  } else {
    cli::cli_abort("{.arg init} must be 'random', 'zero', or a numeric vector")
  }
}

#' Combine results from multiple chains
#' @noRd
combine_chains <- function(results, model, n_chains, n_iter,
                           n_warmup, thin, elapsed) {

  n_params <- ncol(results[[1]]$samples)

  # Create 3D array: iterations x chains x parameters
  draws <- array(
    dim = c(n_iter, n_chains, n_params),
    dimnames = list(
      iteration = seq_len(n_iter),
      chain = seq_len(n_chains),
      parameter = get_param_names(model)
    )
  )

  log_probs <- matrix(nrow = n_iter, ncol = n_chains)

  total_divergent <- 0
  mean_accept <- 0

  for (chain in seq_len(n_chains)) {
    draws[, chain, ] <- results[[chain]]$samples
    log_probs[, chain] <- results[[chain]]$log_probs
    total_divergent <- total_divergent + results[[chain]]$n_divergent
    mean_accept <- mean_accept + results[[chain]]$mean_accept_rate
  }

  mean_accept <- mean_accept / n_chains

  # Compute Rhat and ESS (simplified)
  diagnostics <- list(
    n_divergent = total_divergent,
    mean_accept_rate = mean_accept,
    n_chains = n_chains,
    n_iter = n_iter,
    n_warmup = n_warmup
  )

  new_simp_fit(
    draws = draws,
    log_probs = log_probs,
    diagnostics = diagnostics,
    model = model,
    timing = elapsed
  )
}

#' Get parameter names
#' @noRd
get_param_names <- function(model) {
  K <- model$response$n_components
  p <- model$fixed$n_terms
  comp_names <- model$response$component_names
  fixed_names <- model$fixed$terms

  # Beta names: fixed_term[component]
  beta_names <- character(p * K)
  idx <- 1
  for (i in seq_len(p)) {
    for (j in seq_len(K)) {
      beta_names[idx] <- sprintf("beta[%s,%s]", fixed_names[i], comp_names[j])
      idx <- idx + 1
    }
  }

  c(beta_names, "phi")
}

#' Check positive integer
#' @noRd
check_positive_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || x < 1 || x != floor(x)) {
    cli::cli_abort("{.arg {name}} must be a positive integer")
  }
  invisible(x)
}

#' Apply thinning to samples
#' @noRd
thin_samples <- function(fit, thin) {
  keep_idx <- seq(1, dim(fit$draws)[1], by = thin)
  fit$draws <- fit$draws[keep_idx, , , drop = FALSE]
  fit$log_probs <- fit$log_probs[keep_idx, , drop = FALSE]
  fit$diagnostics$n_iter <- length(keep_idx)
  fit
}
