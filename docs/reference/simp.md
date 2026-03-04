# Specify a SIMP model

Create a model specification for compositional data. This function
validates the formula, data, and priors, and constructs the model
object. The model is not fitted until
[`fit()`](https://gcol33.github.io/SIMP/reference/fit.md) is called
(v0.3.0+).

## Usage

``` r
simp(formula, data, family = simp_dirichlet(), prior = NULL, validate = TRUE)
```

## Arguments

- formula:

  A model formula. The response should be a `simp_composition` or
  `simp_cover` object, or a matrix/data frame of compositional data.
  Random effects use lme4-style syntax: `(1 | group)`,
  `(1 + x | group)`, etc.

- data:

  A data frame containing the predictor variables

- family:

  A family object specifying the likelihood. See
  [`simp_dirichlet()`](https://gcol33.github.io/SIMP/reference/simp_dirichlet.md),
  [`simp_logistic_normal()`](https://gcol33.github.io/SIMP/reference/simp_logistic_normal.md),
  [`simp_hurdle_dirichlet()`](https://gcol33.github.io/SIMP/reference/simp_hurdle_dirichlet.md),
  [`simp_stick_breaking()`](https://gcol33.github.io/SIMP/reference/simp_stick_breaking.md).

- prior:

  A `simp_priors` object, or NULL to use family defaults

- validate:

  Logical; validate the response data (default TRUE)

## Value

A `simp_model` object

## Examples

``` r
# Create some example data
set.seed(42)
n <- 50
site <- factor(rep(1:10, each = 5))
temp <- rnorm(n)

# Compositional response (3 components)
props <- matrix(runif(n * 3), ncol = 3)
props <- props / rowSums(props)
cover <- simp_composition(props, components = c("grass", "shrub", "bare"))

# Create data frame
dat <- data.frame(temp = temp, site = site)

# Basic model specification
model <- simp(cover ~ temp, data = dat)
print(model)

# With random effects
model <- simp(cover ~ temp + (1 | site), data = dat)

# With custom family and priors
model <- simp(
  cover ~ temp + (1 | site),
  data = dat,
  family = simp_dirichlet(),
  prior = simp_priors(beta = prior_normal(0, 5))
)
```
