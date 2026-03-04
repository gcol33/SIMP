# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

SIMP (Structural Intensity Modeling Package) is an R package providing a
Bayesian modeling framework for bounded, compositional, and zero-heavy
intensity data. It treats constraints as defining features of the
data-generating process rather than transforming them away.

## Build and Test Commands

``` bash
# Run all tests
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::test()"

# Run a single test file
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "testthat::test_file('tests/testthat/test-simp_composition.R')"

# Check package (includes tests, examples, and documentation)
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::check()"

# Load package for interactive development
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::load_all()"

# Rebuild documentation
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "devtools::document()"

# Build pkgdown site
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "source('~/.R/build_pkgdown.R'); build_pkgdown_site()"
```

## Architecture

### Core Data Structures (S3 Classes)

- **`simp_composition`** (`R/simp_composition.R`): Base class for
  compositional data where rows sum to 1 and values are in \[0,1\].
  Stores component names, zero patterns, and observation counts as
  attributes.

- **`simp_cover`** (`R/simp_cover.R`): Subclass of `simp_composition`
  for vegetation cover data with layer (stratum) support. Adds species
  names, layer assignments, and residual component tracking.

- **`simp_prior`** / **`simp_priors`** (`R/priors.R`): Prior
  specification system with constructors for common distributions
  (normal, half-normal, half-Cauchy, gamma, beta, PC, LKJ).
  [`simp_priors()`](https://gcol33.github.io/SIMP/reference/simp_priors.md)
  bundles priors for model parameters.

- **`simp_parsed_formula`** (`R/formula.R`): Parsed model formula with
  extracted fixed effects, random effects (including
  nested/crossed/uncorrelated structures), and response information.

### Key Modules

| File | Purpose |
|----|----|
| `R/validate.R` | Diagnostic validation for compositional constraints, generates recommendations |
| `R/zeros.R` | Zero identification (structural vs sampling), multiplicative/additive imputation |
| `R/normalize.R` | Row normalization, residual component handling, log-ratio transforms (CLR/ALR/ILR) |
| `R/collapse.R` | Component aggregation by mapping or layer, rare component collapsing |
| `R/formula.R` | lme4-style formula parsing with random effect extraction |
| `R/utils.R` | Internal validation helpers (not exported) |

### Class Hierarchy

    simp_composition (matrix subclass)
        └── simp_cover (adds layer/species semantics)

Both classes preserve their type through subsetting operations via
custom `[` methods.

### Design Patterns

- Constructors validate by default but accept `validate = FALSE` for
  internal use
- S3 methods provide [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), `[`,
  [`rbind()`](https://rdrr.io/r/base/cbind.html),
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
- Error messages use `cli` package formatting throughout
- Zero patterns tracked as matrix attributes for efficiency

## Testing

Tests are organized by module in `tests/testthat/`. Test files mirror R
source files (e.g., `test-simp_composition.R` tests
`R/simp_composition.R`).

## Vignettes

Theory vignettes (`vignettes/theory-*.Rmd`) cover the mathematical
foundations. Focus is on compositional data constraints, dependence
geometry, and model admissibility.
