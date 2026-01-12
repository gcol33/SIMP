# SIMP Development Roadmap

## Overview

This roadmap outlines the development path from v0.0.1 (current: theory documentation only) to v1.0.0 (stable release with complete Bayesian inference for compositional cover data).

**Key Design Decision**: SIMP uses a native C++ HMC/NUTS backend ported from ratiod, with **no Stan dependency**. This provides:
- Zero external compilation requirements
- Faster startup (no Stan compilation)
- Full control over likelihood implementations
- Proven infrastructure (ratiod v1.2.0)

## Relationship to ratiod

SIMP builds on infrastructure from the ratiod package (v1.2.0), adapting its:

| Component | ratiod Source | SIMP Adaptation |
|-----------|---------------|-----------------|
| HMC/NUTS sampler | `src/hmc_core.cpp`, `hmc_sampler.cpp` | Adapt for Dirichlet/simplex likelihoods |
| Prior framework | `R/priors.R` | Direct port |
| Formula parsing | `R/formula.R` | Extend for multi-component responses |
| Random effects | `src/hmc_*.h` | Port intercepts, slopes, nested, crossed |
| Latent factors | `R/latent.R`, `src/hmc_latent.h` | Port for unmeasured confounders |
| Spatial effects | `R/spatial.R`, `src/hmc_car.h`, `hmc_gp.h` | Port CAR/BYM2/GP |
| Temporal effects | `R/temporal.R` | Port RW1/RW2/AR1 |
| Diagnostics | `R/plot_diagnostics.R` | Direct port |

## Versioning Strategy

| Version | Theme | Key Deliverable |
|---------|-------|-----------------|
| 0.1.0 | Foundation | Data structures, priors, formula parsing from ratiod |
| 0.2.0 | Model Specification | Families, random effects structure |
| 0.3.0 | Inference Engine | Native C++ HMC ported from ratiod |
| 0.4.0 | Posterior Extraction | Prediction, fitted values, contrasts |
| 0.5.0 | Diagnostics | Convergence, pp_check, model comparison |
| 0.6.0 | Real Data | Example datasets, case study vignettes |
| 0.7.0 | Performance | Threading, Laplace approximation |
| 0.8.0 | Advanced Features | Spatial (CAR/GP), temporal (RW/AR) |
| 0.9.0 | API Stabilization | Deprecation cycle, polish |
| 1.0.0 | Stable Release | Production-ready, full documentation |

## Current State (v0.0.1.1)

- 9 theory vignettes (1,039 lines) establishing mathematical framework
- pkgdown documentation site
- No exported functions
- No implementation code
- No tests or example data

## Design Principles

1. **No external compilation**: Native C++ via Rcpp, no Stan/JAGS/INLA dependency
2. **Constraints as features**: Boundedness, compositional closure, and structural zeros define the model
3. **Latent-to-observed separation**: Unconstrained latent processes mapped to constrained observations
4. **Shared structure default**: Components share random effects unless explicitly independent
5. **ratiod-proven infrastructure**: Build on tested code rather than reimplementing

## Random Effects Architecture (from ratiod)

SIMP inherits ratiod's flexible random effects system:

### Supported Structures

```r
# Random intercepts
cover ~ temp + (1 | site)

# Random slopes (uncorrelated with intercept)
cover ~ temp + (1 | site) + (0 + temp | site)

# Correlated random intercepts and slopes
cover ~ temp + (1 + temp | site)

# Nested random effects
cover ~ temp + (1 | region/site)
# Expands to: (1 | region) + (1 | region:site)

# Crossed random effects
cover ~ temp + (1 | site) + (1 | year) + (1 | observer)

# Multiple grouping factors with different terms
cover ~ temp + precip + (1 + temp | site) + (1 | year)
```

### Implementation (from ratiod)

Random effects use sparse Cholesky parameterization:
- `src/hmc_random_effects.h`: Core RE computations
- Automatic detection of nested vs crossed structures
- Efficient sparse matrix operations for large grouping factors
- Sum-to-zero constraints for identifiability

## File Index

- [v0.1.0 - Foundation](v0.1.0.md)
- [v0.2.0 - Model Specification](v0.2.0.md)
- [v0.3.0 - Inference Engine](v0.3.0.md)
- [v0.4.0 - Posterior Extraction](v0.4.0.md)
- [v0.5.0 - Diagnostics & Visualization](v0.5.0.md)
- [v0.6.0 - Datasets & Case Studies](v0.6.0.md)
- [v0.7.0 - Performance](v0.7.0.md)
- [v0.8.0 - Advanced Features](v0.8.0.md)
- [v0.9.0 - API Stabilization](v0.9.0.md)
- [v1.0.0 - Stable Release](v1.0.0.md)
