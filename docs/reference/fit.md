# Fit a SIMP model

Run Bayesian inference on a SIMP model specification using the NUTS
sampler.

## Usage

``` r
fit(
  model,
  chains = 4,
  iter = 2000,
  warmup = NULL,
  thin = 1,
  cores = 1,
  seed = NULL,
  init = "random",
  max_treedepth = 10,
  adapt_delta = 0.8,
  verbose = TRUE
)
```

## Arguments

- model:

  A `simp_model` object created by
  [`simp()`](https://gcol33.github.io/SIMP/reference/simp.md)

- chains:

  Number of MCMC chains (default 4)

- iter:

  Number of sampling iterations per chain (default 2000)

- warmup:

  Number of warmup iterations per chain (default half of iter)

- thin:

  Thinning interval (default 1, no thinning)

- cores:

  Number of cores for parallel chains (default 1)

- seed:

  Random seed for reproducibility

- init:

  Initialization method: "random" (default), "zero", or a list of
  initial values

- max_treedepth:

  Maximum tree depth for NUTS (default 10)

- adapt_delta:

  Target acceptance rate (default 0.8)

- verbose:

  Print progress (default TRUE)

## Value

A `simp_fit` object containing posterior samples and diagnostics

## Examples

``` r
if (FALSE) { # \dontrun{
# Create model
model <- simp(cover ~ temp + (1 | site), data = dat)

# Fit model
fit <- fit(model, chains = 4, iter = 2000)

# Examine results
summary(fit)
} # }
```
