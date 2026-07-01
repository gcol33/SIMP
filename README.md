# SIMP (superseded)

[![Lifecycle: superseded](https://img.shields.io/badge/lifecycle-superseded-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#superseded)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

SIMP's compositional models now live in **[tulpaComp](https://github.com/gcol33/tulpaComp)**, where each family is a likelihood plugin on the [tulpa](https://github.com/gcol33/tulpa) inference engine. Running on the engine gives them shared exact-MCMC (NUTS) and Laplace tiers with autodiff gradients. This repository is archived for reference and no longer developed.

The four families SIMP set out to build are all implemented in tulpaComp against the engine, each with multi-seed parameter-recovery and credible-interval coverage tests:

- **Dirichlet**: the natural distribution on the simplex; concentration estimated or fixed.
- **logistic-normal**: multivariate normal on additive log-ratios, for flexible correlation between components.
- **hurdle-Dirichlet**: two-part model for structural zeros, with a covariate-dependent zero submodel.
- **stick-breaking**: sequential binary splits, for components with a natural order.

The zero-handling data structures (`composition()`, `identify_zeros()`, structural vs sampling zeros) carried over as well.

## Theory articles

The nine [theory articles](https://gcol33.github.io/SIMP/) that develop the mathematical framework remain available here. They document the simplex likelihoods, the zero taxonomy, and the prior and formula design that the tulpaComp implementation follows.

## License

MIT (see the LICENSE file)
