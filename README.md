# SIMP

*models for compositions full of zeros*

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Bayesian framework for bounded, compositional, zero-heavy data, sampled with a native C++ NUTS backend (no Stan).**

Cover, market shares, budget shares, relative abundances: parts of a finite whole that sum to one, never go negative, and hit exact zero whenever a part is genuinely absent. SIMP treats those bounds and zeros as the data-generating process and builds the likelihood around them, so nothing has to be nudged off the boundary to fit a transform.

```r
library(SIMP)

# compositional response: each row sums to 1
props <- matrix(runif(50 * 3), ncol = 3)
props <- props / rowSums(props)
cover <- simp_composition(props, components = c("grass", "shrub", "bare"))

dat <- data.frame(temp = rnorm(50), site = factor(rep(1:10, each = 5)))

# specify a hierarchical Dirichlet model, lme4-style random effects
model <- simp(cover ~ temp + (1 | site), data = dat, family = simp_dirichlet())

# sample the posterior with the built-in NUTS sampler
fit <- fit(model, chains = 4, iter = 2000)
summary(fit)
```

## Zeros are kept exact

The usual route for compositional data is a log-ratio or logit transform onto an unbounded scale, a Gaussian model, and a back-transform. The transform is undefined at zero, so the zeros get a small constant added, replaced by half a detection limit, or dropped. Each fix moves the data off the boundary it actually sits on.

SIMP models the simplex directly. Boundedness, the unit-sum closure, and structural zeros are constraints the likelihood respects rather than artefacts to smooth away. `simp_hurdle_dirichlet()` gives zeros their own two-part model: a probability that a component is absent, and a Dirichlet over what remains.

## Families

- **`simp_dirichlet()`**: the natural distribution on the simplex; concentration estimated or fixed.
- **`simp_logistic_normal()`**: multivariate normal on additive log-ratios, for flexible correlation between components.
- **`simp_hurdle_dirichlet()`**: two-part model for structural zeros, with a covariate-dependent zero submodel.
- **`simp_stick_breaking()`**: sequential binary splits, for components with a natural order (e.g. canopy strata).

Pick a family, then `simp()` validates the response, parses the formula, and checks the family against the data (it warns if you hand zeros to a family that cannot hold them).

## Compositional and cover data structures

`simp_composition()` is the core class: a matrix whose rows sum to one, carrying its component names and zero pattern. `simp_cover()` extends it for vegetation surveys with layer (stratum) and residual support.

```r
df <- data.frame(grass = c(0.3, 0.4), shrub = c(0.5, 0.4), bare = c(0.2, 0.2))
comp <- simp_composition(df)

# get cover data onto the simplex before modelling
normalize_cover(raw)        # rescale rows that sum above or below one
collapse_layers(x)          # combine tree / shrub / herb strata
add_residual(x)             # close the composition with a residual category
identify_zeros(comp)        # classify zeros as structural or sampling
```

## Priors and formulas

Priors are family-aware closures (`prior_normal()`, `prior_gamma()`, `prior_lkj()`, `prior_pc()`, and more), assembled with `simp_priors()` or left to family defaults. Formulas use lme4-style random effects, including correlated slopes, nested (`region/site`), and crossed groupings; `parse_formula()` exposes the parsed structure.

```r
simp(
  cover ~ temp + precip + (1 + temp | site) + (1 | year),
  data = dat,
  family = simp_dirichlet(),
  prior = simp_priors(beta = prior_normal(0, 5))
)
```

## Status

SIMP is at v0.1.0 and under active development. The data structures, families, prior and formula machinery are in place, and `fit()` runs the native NUTS sampler for the Dirichlet family. Logistic-normal, hurdle, and stick-breaking likelihoods, posterior prediction, and spatial/temporal effects are on the [roadmap](roadmap/ROADMAP.md). The nine [theory articles](https://gcol33.github.io/SIMP/) develop the framework the implementation follows.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/SIMP")
```

## Documentation

- [Theory articles](https://gcol33.github.io/SIMP/) — the mathematical framework, in nine chapters
- [Reference](https://gcol33.github.io/SIMP/reference/) — function documentation

## License

MIT (see the LICENSE file)
