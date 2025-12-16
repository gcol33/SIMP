# Admissible Models for Compositional Data

## The constraint on models

The previous chapters established that cover data live on a constrained
space: the simplex. A composition of $`K`$ components is a vector
$`\mathbf{y} = (y_1, \ldots, y_K)`$ satisfying $`y_k \geq 0`$ for all
$`k`$ and $`\sum_{k=1}^{K} y_k = 1`$. This space is denoted
$`\mathcal{S}^{K-1}`$, the $`(K-1)`$-dimensional simplex embedded in
$`\mathbb{R}^K`$.

Any probability model for cover data must assign probability only to
points in this space. A distribution that places positive mass outside
$`\mathcal{S}^{K-1}`$ (on negative values, on values exceeding one, or
on vectors that do not sum to one) is an invalid model. It describes a
quantity other than the one measured.

This constraint must be satisfied before any other modelling
consideration.

## What admissibility requires

A probability distribution is admissible for compositional data if and
only if:

1.  Its support is exactly $`\mathcal{S}^{K-1}`$ (or a well-defined
    subset of it).
2.  It defines a proper probability measure on this space.
3.  It does not require the data to be modified before the distribution
    applies.

The third condition excludes models that are nominally defined on the
simplex but require transformations that fail at the boundary. A model
that requires adding constants to zeros, or that is undefined when any
component equals zero or one, is not admissible for data that include
such values.

These conditions are minimal. They say nothing about whether a model is
useful, flexible, or estimable. They say only whether the model is
coherent with the data it purports to describe.

## The Dirichlet family

The Dirichlet distribution is the canonical probability distribution on
the simplex. For a composition $`\mathbf{y} \in \mathcal{S}^{K-1}`$, the
Dirichlet density with concentration parameters
$`\boldsymbol{\alpha} = (\alpha_1, \ldots, \alpha_K)`$, where
$`\alpha_k > 0`$, is:

``` math
p(\mathbf{y} \mid \boldsymbol{\alpha}) = \frac{\Gamma\left(\sum_{k=1}^{K} \alpha_k\right)}{\prod_{k=1}^{K} \Gamma(\alpha_k)} \prod_{k=1}^{K} y_k^{\alpha_k - 1}
```

The Dirichlet is admissible. It is defined natively on the constrained
space, and no transformation is needed. Its behaviour near the simplex
boundary depends on the concentration parameters: when $`\alpha_k < 1`$,
the density is unbounded as $`y_k \to 0`$, concentrating mass near the
vertices and edges of the simplex; when $`\alpha_k = 1`$, the
distribution is uniform; when $`\alpha_k > 1`$, the density vanishes at
the boundary, concentrating mass in the interior. This flexibility
allows the Dirichlet to represent compositions that tend toward sparse
configurations (few dominant components) or toward balanced
configurations (mass spread across components).

The Dirichlet encodes specific assumptions about dependence. Components
are negatively correlated, but the correlation structure is constrained:
the correlation between any two components $`y_i`$ and $`y_j`$ is
determined entirely by the concentration parameters and is always
negative. There is no freedom to specify arbitrary pairwise
correlations. The covariance matrix of a Dirichlet-distributed vector
has rank $`K-1`$, reflecting the sum constraint, and its off-diagonal
elements are all negative.

This is a strong structural assumption. The Dirichlet assumes that all
dependence among components arises from the constraint and from a single
shared precision parameter (the sum of concentrations). Components that
share a simplex but have more complex relationships—pairs that are more
or less strongly coupled than the Dirichlet permits—cannot be
represented.

## The logistic-normal family

The logistic-normal distribution offers a different parameterisation. It
models compositions by assuming that log-ratios of components follow a
multivariate normal distribution. There are several equivalent
formulations. The additive log-ratio (alr) representation selects one
component as reference and models the log-ratios relative to it; the
centred log-ratio (clr) uses the geometric mean as reference; the
isometric log-ratio (ilr) constructs an orthonormal basis in log-ratio
space. These are different coordinate systems for the same family of
distributions.

In the alr formulation, if $`\mathbf{y} \in \mathcal{S}^{K-1}`$ and we
define $`\mathbf{z} = (\log(y_1/y_K), \ldots, \log(y_{K-1}/y_K))`$ using
the $`K`$-th component as reference, then the logistic-normal assumes
$`\mathbf{z} \sim \mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})`$
for some mean vector $`\boldsymbol{\mu}`$ and covariance matrix
$`\boldsymbol{\Sigma}`$. The choice of reference component is arbitrary
and does not affect the distribution on the simplex, only the
parameterisation.

The logistic-normal is admissible in a restricted sense. It defines a
proper distribution on the interior of the simplex,
$`\mathcal{S}^{K-1}_\circ = \{\mathbf{y} : y_k > 0 \text{ for all } k\}`$.
It assigns zero probability to the boundary—to any composition where at
least one component is exactly zero.

This is a serious limitation for cover data. Structural zeros are often
the majority of observations. A model that assigns zero probability to
exact zeros cannot represent the data-generating process. It treats a
central feature of the data as impossible.

Within its domain, the logistic-normal permits more flexible covariance
structures than the Dirichlet. The matrix $`\boldsymbol{\Sigma}`$ can
encode arbitrary correlations among log-ratios, allowing components to
be more or less tightly coupled than the Dirichlet would permit. This
flexibility is valuable when ecological or observational processes
induce non-exchangeable dependence.

But this flexibility operates only on compositions without zeros. For
data with structural zeros, the logistic-normal is inadmissible as a
complete model.

## Transformation-based approaches

Many analyses of compositional data proceed by transforming the data to
an unconstrained space, fitting a standard model (typically Gaussian),
and interpreting results on the original scale. Common transformations
include the additive log-ratio (alr), centred log-ratio (clr), and
isometric log-ratio (ilr).

These approaches are not admissible models for cover data. They fail on
two grounds.

First, all log-ratio transformations are undefined at zero. A
composition with any $`y_k = 0`$ cannot be transformed. The standard
remedy (adding a small constant before transformation) changes the data.
The resulting analysis is not an analysis of the observed composition
but of a perturbed version of it. The perturbation is arbitrary, and
results can depend on the constant chosen.

Second, fitting a Gaussian model to transformed compositions and
back-transforming does not yield a valid distribution on the simplex.
The Gaussian places positive probability on the entire real line.
Back-transformation maps some of this probability to impossible
compositions: negative proportions, proportions exceeding one, or
vectors that do not sum to one. The resulting “distribution” on the
original scale is not a probability distribution at all.

These approaches may produce numbers. They do not produce valid
probabilistic statements about compositions.

## What geometry forces

The simplex geometry constrains not just which distributions are
admissible, but what their covariance structures can look like.

Any random vector $`\mathbf{y}`$ on $`\mathcal{S}^{K-1}`$ must satisfy
$`\sum_k y_k = 1`$ with probability one. This implies that the
covariance matrix $`\text{Cov}(\mathbf{y})`$ is singular: it has rank at
most $`K-1`$. The rows and columns sum to zero. The constraint forces
this property; it does not arise from specific distributional choices.

Furthermore, not all covariance matrices with this structure are
achievable. The matrix must be consistent with the non-negativity
constraints $`y_k \geq 0`$. Consider a simple example: if $`K = 3`$ and
one component has very high variance while the others have low variance,
the high-variance component would frequently take values outside
$`[0, 1]`$ unless the covariances are carefully calibrated. More
generally, variances near the boundary of the simplex are constrained by
how close the mean is to that boundary. A component with mean 0.05
cannot have the same variance as a component with mean 0.50 without
violating the non-negativity constraint.

The space of valid covariance matrices for compositions is therefore
smaller than it might appear. Specifying an arbitrary
$`(K-1) \times (K-1)`$ covariance matrix for log-ratios, or an arbitrary
$`K \times K`$ covariance matrix for proportions, does not guarantee
that a corresponding distribution exists on the simplex. Model
specification must respect these geometric limits.

## What remains optional

Within the class of admissible models, genuine choices remain. These
include:

**The form of marginal distributions.** The Dirichlet implies Beta
marginals for individual components. Other admissible distributions on
the simplex can produce different marginal shapes—more or less peaked,
more or less skewed near the boundaries.

**The dependence structure beyond the constraint.** The sum constraint
forces some dependence, but does not determine all of it. Different
admissible models encode different assumptions about how tightly
components are coupled, whether some pairs are more dependent than
others, and how dependence changes across the simplex.

**Behaviour near boundaries.** Models differ in how they treat
compositions close to the simplex boundary, where one or more components
are near zero or near one. Some concentrate mass in the interior; others
allow mass to accumulate near edges or vertices.

These are modelling choices, not forced by geometry. But they are
choices within a constrained menu. The geometry rules out most of the
models one might naively consider; what remains is a small family with
specific, identifiable assumptions.

## The status of zeros

This chapter has focused on the simplex as a continuous space and on
distributions defined on its interior. The treatment of exact
zeros—structural absences—requires additional machinery and will be
addressed in the next chapter.

For now, it is sufficient to note that exact zeros create a boundary
problem for any continuous distribution on the simplex. The Dirichlet
can place mass near zero but not exactly at zero when
$`\alpha_k \geq 1`$; it can place mass at zero only in a limiting sense
when $`\alpha_k < 1`$. The logistic-normal excludes zeros entirely.
Neither family provides a complete model for compositions with
structural zeros.

This is not a minor gap. It is a fundamental limitation that any
complete framework must address.

## Summary

The requirement that probability models respect the simplex geometry is
a hard constraint that eliminates most familiar modelling options.

Gaussian models on raw or transformed compositions are inadmissible.
They place probability outside the feasible space.

Log-ratio transformations require modifying data with zeros. The
modifications are arbitrary and change what is being modelled.

The Dirichlet and logistic-normal families are admissible on the simplex
interior, but encode specific assumptions about dependence. The
Dirichlet permits only the dependence induced by shared concentration.
The logistic-normal permits flexible log-ratio covariance but assigns
zero probability to compositions with zeros.

Any valid model for cover data must be defined on the simplex, assign
probability only to feasible compositions, and handle structural zeros
without requiring data modification. The space of such models is small.
This clarifies what honest modelling requires.
