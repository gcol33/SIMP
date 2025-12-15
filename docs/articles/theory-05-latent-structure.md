# Latent Structure and Hierarchy

## The gap between process and observation

Vegetation cover as recorded in a plot is not the same thing as the
ecological process that generates it. The recorded value is a
constrained, discretised, and bounded summary of something more
fundamental: the capacity of a species to occupy space under the
conditions present at that location.

This distinction matters. The ecological quantity of interest—call it
suitability, intensity, or potential—is not inherently bounded between
zero and one. It is not inherently compositional. It does not sum to one
across species because ecology demands it. These constraints arise from
the act of measurement: we observe cover within a finite plot, so the
observed values must partition that plot.

The previous chapters established that constraints, dependence, and
zeros are unavoidable features of observed cover data. This chapter
argues that these features belong to the observation process, not to the
underlying ecological quantities. Modelling observed values directly
conflates process and measurement. Separating them requires latent
structure.

## Why direct modelling fails

Consider attempting to model observed cover proportions as the primary
process. The model must:

1.  Respect the simplex constraint (values non-negative, summing to
    one).
2.  Allow exact zeros with positive probability.
3.  Represent dependence among components that is partly geometric and
    partly ecological.
4.  Accommodate variation that differs in character across species.

No single distribution on the simplex accomplishes all of this. Chapter
3 established that admissible distributions on the simplex have
restricted covariance structures. Chapter 4 established that continuous
distributions cannot represent exact zeros. Attempting to build a model
that does everything at the observation level leads to either: (a) a
model that violates the constraints, (b) a model that cannot represent
zeros, or (c) a model with an impoverished dependence structure.

The problem is not that we lack sufficiently flexible distributions. The
problem is that we are asking one layer to do too much. The simplex
constraint, the zero structure, and the ecological variation are
distinct phenomena. Collapsing them into a single distributional choice
forces compromises that distort what the model represents.

## What latent variables accomplish

A latent variable is an unobserved quantity that mediates between the
process of interest and the data. In the context of cover data, latent
variables serve a specific structural role: they carry ecological
variation in an unconstrained space, before the observation process
imposes constraints.

Consider a latent quantity for each species $`k`$ at each plot—call it
an intensity, a suitability, or simply a latent value. The exact
parameterisation and interpretation will be specified in later chapters.
What matters here is the structural role: this latent quantity
represents the species’ tendency to occupy space under the conditions
present, before the observation process imposes constraints.

The latent values are not bounded. They are not required to sum to
anything across species. They can be modelled with familiar tools:
Gaussian distributions, linear predictors, covariance structures that
encode ecological similarity or spatial proximity.

The observed cover proportions are then generated from these latent
intensities through a mapping that imposes the simplex constraint and
determines presence or absence. The constraints that seemed burdensome
at the observation level become natural consequences of the mapping.
They are not forced onto the model; they emerge from structure.

This separation clarifies what each part of the model represents:

- **Latent layer**: Ecological variation, species-specific effects,
  environmental responses. Modelled in unconstrained space with flexible
  dependence.
- **Observation layer**: The constraint that cover must sum to one, the
  discretisation into presence/absence, the measurement of proportions
  within a bounded plot.

The observed data are projections of latent structure through an
observation process. The projection imposes constraints. The latent
structure is where ecological inference happens.

## Hierarchy as logical necessity

The term “hierarchical model” sometimes suggests an optional
elaboration—a way to add complexity when simpler models fail. For
constrained intensity data, this framing is backwards. Hierarchy is not
added to a simpler model. It is what remains when we refuse to conflate
distinct sources of variation.

The structure can be written schematically as a factorisation of the
joint distribution:

``` math
p(\text{data}, \text{latent}) = p(\text{data} \mid \text{latent}) \times p(\text{latent})
```

The first term, $`p(\text{data} \mid \text{latent})`$, is the
observation model. It describes how latent intensities become observed
proportions. This is where the simplex constraint lives. This is where
zeros are generated. This is where the compositional structure is
imposed.

The second term, $`p(\text{latent})`$, is the process model. It
describes how latent intensities vary across species, plots, and
conditions. This is where ecological covariance belongs. This is where
environmental predictors enter. This is where species-level parameters
are specified.

The factorisation is not a modelling convenience. It reflects the
structure of the problem. Observed cover data are conditionally
independent given the latent intensities: once we know the underlying
ecological state, the observation process is determined. The dependence
among observations arises from shared latent structure, not from direct
connections among observed values.

## Zeros revisited

Chapter 4 established that zeros require discrete structure to represent
exactly. In a hierarchical view, this structure has a natural home.

The latent intensity $`\lambda_k`$ for a species can be positive,
indicating some capacity to occupy space, or it can be such that the
species has zero probability of presence. The observation model then
includes a presence indicator: given the latent intensity, is the
species present in this plot? If yes, what proportion of space does it
occupy?

This two-stage structure—presence given latent intensity, then
allocation given presence—respects the compositional constraint while
placing probability mass exactly at zero. The zeros are not
approximations. They are outcomes of a discrete process that operates on
latent quantities.

Importantly, presence and allocation are not modelled independently
across species. The latent intensities are correlated. The presence
indicators depend on these correlated intensities. The allocation among
present species is constrained to sum to one minus the absent species’
share. The hierarchy couples these components without forcing a single
distribution to represent everything.

## What dependence means in this view

Chapter 2 established that compositional dependence is geometric and
unavoidable. Chapter 3 showed that different distributions on the
simplex encode different dependence assumptions. In a hierarchical view,
these facts find their proper place.

The observed dependence among cover proportions arises from two sources:

1.  **Latent covariance**: Species with similar ecological requirements
    will have correlated latent intensities. This is the ecological part
    of dependence—the part we want to learn about.

2.  **Constraint-induced dependence**: Even if latent intensities were
    independent, observed proportions would be dependent because they
    must sum to one. This is the geometric part—the part that is
    structurally guaranteed.

Direct modelling of observed proportions conflates these sources. The
Dirichlet, for instance, imposes a specific mixture of geometric and
ecological dependence that cannot be disentangled. The logistic-normal
allows more flexible dependence but only among present species.

In a hierarchical model, the sources are separated. Latent covariance
encodes ecological dependence. The observation model encodes geometric
dependence. Inference can target the latent covariance, which represents
the scientifically interesting structure, rather than the observed
covariance, which mixes signal with constraint.

## Conditional independence and its implications

A key feature of hierarchical structure is conditional independence.
Given the latent intensities, the observed cover values are generated
through a process that does not require further connections among plots
or among species beyond what the latent layer provides.

This does not mean observations are independent. It means their
dependence is mediated. Two plots with similar observed cover have
similar observations because they have similar latent structure, not
because one plot’s cover directly influences another’s.

Conditional independence has practical implications:

- The likelihood factorises into terms that depend on individual
  observations given their latent values.
- Complex dependence structures are built into the latent layer, not the
  observation layer.
- The observation model can be relatively simple—a mapping from latent
  intensities to proportions—while the process model carries the
  ecological complexity.

This is not a simplification for computational convenience. It is a
statement about how the world works: observations are conditionally
independent because the data-generating process goes through latent
ecological states.

## What hierarchy does not do

Introducing latent structure does not solve all problems. It clarifies
structure but does not guarantee identifiability. If two different
latent configurations produce the same distribution of observed data,
the data cannot distinguish them.

For compositional data with zeros, identifiability is a genuine concern.
When most species are absent from most plots, the latent values for
absent species are weakly informed by the data. We observe that the
species is absent, but many latent configurations are consistent with
absence. The data constrain the latent values to be “low enough” to
produce absence, but do not pinpoint them precisely.

Similarly, the allocation of cover among present species provides
information about relative latent values, but the absolute scale may be
poorly determined. The mapping from latent values to observed
proportions involves normalisation, which can absorb shifts in the
overall level of latent values.

A hierarchical view does not eliminate these limitations. It makes them
explicit. The latent layer is partly identified by the data (through
observed presences and cover values) and partly regularised by structure
(through the process model and any priors). What the data cannot tell
us, the model cannot recover—regardless of how it is parameterised. The
role of structural assumptions and priors is to provide coherent
behaviour in weakly identified regions, not to fabricate information.

This is a feature, not a bug. A model that obscures identifiability
problems is harder to diagnose than one that exposes them. Hierarchical
structure exposes them by separating what is observed from what is
inferred.

## Summary

The constraints, dependence, and zeros of cover data are features of
observation, not of the underlying ecological process. Modelling
observations directly requires a single distributional family to
represent all of these features simultaneously. No such family exists
that is both admissible and flexible.

The alternative is to separate process from observation through latent
structure. Latent intensities carry ecological variation in an
unconstrained space. The observation model maps these intensities to
constrained proportions, imposing the simplex structure and generating
zeros through a discrete mechanism.

This hierarchy is not a modelling convenience. It is the structure that
remains when we take the constraints seriously. The observed data are
projections of latent ecological states. The projections are
constrained; the states are not. Inference about ecology operates on the
latent layer.

This chapter establishes only that the structure is necessary. A
framework that operates directly on observed proportions, without
acknowledging latent mediation, has not separated process from
measurement. It has conflated them.
