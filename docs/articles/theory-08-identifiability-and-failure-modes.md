# Identifiability and Failure Modes

## What identifiability means

A parameter is identifiable if different values of that parameter imply
different distributions over observable data. If two distinct parameter
configurations produce exactly the same likelihood for all possible
datasets, the parameters are not identifiable. The data cannot
distinguish between them, regardless of sample size, estimation method,
or computational resources.

This definition concerns information, not estimation. An identifiable
parameter may still be difficult to estimate (requiring large samples,
careful priors, or specialised algorithms). An unidentifiable parameter
cannot be estimated at all, in any principled sense. No amount of data
resolves the ambiguity because there is no ambiguity to resolve: the
likelihood surface is flat in the relevant direction.

In hierarchical models with latent structure, identifiability becomes
more complex. The parameters of interest include not just fixed effects
or variance components but the entire latent configuration (the values
of unobserved quantities at each unit). Some of these quantities may be
well-identified; others may be weakly constrained; still others may be
fundamentally underdetermined given the model structure.

The previous chapters established that SIMP operates through latent
structure: ecological variation lives in an unconstrained latent space,
and observations are constrained projections of that space. This
architecture creates identifiability properties that must be understood,
not assumed away.

## The scale ambiguity

The observation model maps latent values to constrained proportions.
This mapping involves normalisation: the cover allocated to each present
species depends on relative latent values.

Consider a simple case. Suppose two species are present in a plot, with
latent values \\\lambda_1\\ and \\\lambda_2\\. The cover allocation
depends on how these values compare (perhaps through a softmax or
similar transformation). If we shift both values by the same constant,
the relative comparison is unchanged. If we scale both values by the
same factor, the relative comparison may also be unchanged, depending on
the functional form.

This creates a scale ambiguity. The data constrain relative latent
values but not their absolute level. A configuration where all latent
values are high produces the same observable proportions as one where
all are low, provided the relative structure is preserved.

Not a defect of any particular parameterisation. Structural. The simplex
constraint means that only ratios matter for allocation. The absolute
scale of latent values is not identified by cover data alone.

In practice, this ambiguity is resolved through constraints: fixing one
latent value to a reference, using sum-to-zero constraints, or imposing
priors that anchor the scale. These constraints are assumptions. They
make the model formally identified but do not add information. The data
remain silent about absolute scale; the constraints substitute a
modelling decision for missing information.

## Weak identification from absence

Chapter 4 established that zeros are structural features of cover data,
not nuisances. Chapter 5 noted that latent values for absent species are
weakly informed. This section makes the identifiability implications
explicit.

When a species is absent from a plot, the data tell us that its latent
value was “low enough” to produce absence. Depending on the observation
model, this may mean the value fell below some threshold, or that the
probability of presence given the latent value was low. Either way, the
observation is censored: we learn that the latent value is in some
region, but not where within that region.

A fundamental asymmetry. For present species, cover values constrain
latent values relatively tightly (higher cover implies higher relative
latent value). For absent species, the constraint is one-sided and weak.
The latent value could be slightly below the presence threshold or far
below it; the data do not distinguish.

When zeros dominate (when most species are absent from most plots) most
latent values are weakly identified. The full latent configuration has
high dimension, but the data inform only a small fraction of it. The
remainder is determined by assumptions: the prior distribution over
latent values, the structure imposed by covariates, the dependence
structure that borrows information across units.

Not a problem to be solved. A property of the data structure. Absence is
less informative than presence. A framework that represents zeros
honestly must also represent the limits of what zeros can tell us.

## Confounding between presence and allocation

The SIMP framework separates two processes: which species are present,
and how cover is allocated among present species. Both processes depend
on latent values. This creates potential for confounding.

Consider a covariate (say, soil moisture) that affects both presence and
allocation. Wetter sites may have higher probability of certain species
being present, and among present species, wetter sites may allocate more
cover to moisture-loving species. The covariate acts on both stages.

The two effects are not independently identified. If we observe that a
species has high cover in wet plots, this could reflect: (a) higher
presence probability in wet plots, so the species appears in more wet
plots; (b) higher allocation in wet plots, so the species gets more
cover when present; or (c) both. The observed pattern is consistent with
multiple combinations of presence and allocation effects.

This confounding is intrinsic to the two-stage structure. The stages are
not observed separately. We observe the final cover, which reflects the
composition of both processes. Separating their contributions requires
either: strong assumptions about how they differ in form, additional
data that informs one stage independently, or acceptance that only the
combined effect is identified.

The confounding is particularly severe when zeros are common. With few
observations of presence, the presence model is estimated from limited
data. The allocation model is estimated only from plots where the
species occurs. The two sources of information may be nearly
non-overlapping for rare species.

## Flat regions in the likelihood

The likelihood function summarises what the data say about parameters.
In well-identified models, the likelihood has a clear peak: moving away
from the maximum reduces the likelihood substantially. In weakly
identified models, the likelihood surface has flat regions, directions
along which the likelihood changes little or not at all.

For compositional data with latent structure, flat regions arise from
several sources.

**Scale invariance.** As noted above, normalisation removes information
about absolute scale. The likelihood is constant along directions that
shift all latent values uniformly.

**Absent species.** Latent values for absent species are constrained
only to be “low.” The likelihood varies weakly across configurations
that satisfy this constraint.

**Exchangeable structure.** If latent values share a common distribution
with estimated variance, the data constrain the variance but not
individual deviations from the mean. Many configurations of individual
values are consistent with the same variance estimate.

**Dependence parameters.** Structured dependence (spatial, group-level)
introduces parameters controlling how latent values covary. When most
values are weakly identified individually, the covariance structure is
also weakly identified. The data may constrain the overall strength of
dependence but not its detailed form.

These flat regions are not numerical artefacts. Properties of the
information structure. The likelihood is genuinely uninformative about
certain aspects of the parameter space.

## Ridges and near-ridges

Distinct from flat regions are ridges: directions along which different
parameter combinations yield nearly identical likelihoods. Ridges
indicate trade-offs. If one parameter increases, another can compensate
to maintain the same fit.

In hierarchical models, common ridges include:

**Mean-variance trade-offs.** A model with higher latent mean and higher
variance may fit comparably to one with lower mean and lower variance.
The data constrain a combination of these quantities but not each
separately.

**Intercept-covariate trade-offs.** A model with a high intercept and
negative covariate effect may fit comparably to one with a lower
intercept and less negative effect. The data constrain fitted values,
not the decomposition into components.

**Individual-group trade-offs.** In hierarchical models, variation can
be attributed to individual-level deviations or group-level effects. The
total variation is constrained; its partition is not.

Ridges create practical difficulties for estimation: algorithms may
wander along them, and posterior distributions may be elongated or
multimodal. But the underlying issue is informational. The data do not
contain enough information to resolve the trade-off. Any resolution
comes from assumptions: priors, constraints, or regularisation.

## How structured dependence can mislead

Chapter 6 established that structured dependence belongs at the latent
level. Chapter 7 established that any finite representation of
dependence is a modelling commitment. This section addresses how
dependence structure interacts with identifiability.

Structured dependence borrows information across units. If latent values
at nearby plots are assumed similar, observations at one plot inform
estimates at its neighbours. This borrowing can improve estimates when
the assumption is correct and the data are sparse.

But borrowing is not learning. It propagates assumptions as much as it
propagates data. If the dependence structure is misspecified (if the
true process has different range, different smoothness, or different
geometry than assumed) the borrowed information is systematically wrong.

The identifiability concern: misspecification may not be detectable. A
model with the wrong dependence structure can still fit the data well if
it adjusts other parameters to compensate. The data may not distinguish
between: (a) strong local dependence with weak covariate effects, and
(b) weak dependence with strong covariate effects. Both can produce
similar patterns of spatial similarity in the fitted values.

This confounding between dependence and fixed effects is well-known in
spatial statistics. Not unique to compositional data. But compositional
data exacerbate it because the weak identification of absent species
means there is less information to discipline the dependence structure.
The dependence parameters are estimated primarily from present species,
which may be a biased or unrepresentative subset.

## Why more data does not always help

A natural response to identifiability concerns is to collect more data.
More plots, more species, more covariates. For some limitations, this
helps. For others, it does not.

More data helps when the limitation is statistical uncertainty, when the
parameter is identified but estimated imprecisely due to small samples.
Additional observations tighten confidence intervals and resolve
ambiguity.

More data does not help when the limitation is structural. If two
parameter configurations produce identical likelihoods for all possible
datasets, more observations from the same data-generating process
provide no new information. The configurations remain indistinguishable.

For SIMP, several limitations are structural:

- The scale ambiguity is resolved by constraints, not by data.

- The weak identification of absent species’ latent values is a property
  of censoring, not sample size. Adding more plots where the species is
  absent does not inform its latent value; only plots where it is
  present do.

- Confounding between presence and allocation is a property of the
  two-stage structure, not the number of observations.

- Confounding between dependence and fixed effects depends on the design
  (how covariates vary in space), not just the number of units.

This does not mean data collection is futile. Data collection should be
strategic, targeting the specific identifiability limitations that
affect the scientific question. More of the same data may not resolve
ambiguities that are structural in nature.

## The role of assumptions

In regions where the likelihood is flat or nearly flat, inference is
determined by assumptions. In Bayesian terms, the posterior reflects the
prior. In frequentist terms, estimates are determined by constraints or
regularisation rather than data.

This is how inference works in partially identified models. The honest
response: acknowledge which components are informed by data and which by
assumptions, and assess sensitivity to those assumptions.

For SIMP, key assumptions that fill identifiability gaps include:

**Prior distributions on latent values.** These determine the baseline
level and spread of latent values when the data are uninformative. For
absent species, the prior is doing most of the work.

**Prior distributions on variance components.** These determine how
variation is partitioned when the data do not strongly constrain the
partition.

**The form of the dependence structure.** Assuming spatial smoothness,
exchangeable groups, or independence is a substantive choice that
affects estimates even in data-rich settings.

**The observation model linking latent values to cover.** The functional
form of this mapping affects which latent configurations are consistent
with observed data.

These assumptions are structural components of the model that determine
behaviour in weakly identified regions. A user who does not examine
sensitivity to these assumptions does not understand what the model is
doing.

## Failure modes

The identifiability properties described above produce characteristic
failure modes. Not bugs or rare pathologies. Predictable consequences of
the framework’s structure.

**Overconfident inference about rare species.** For species absent from
most plots, the apparent precision of estimates reflects the prior more
than the data. Narrow credible intervals may indicate strong prior
assumptions, not strong data support.

**Unstable covariate effects.** When presence-allocation confounding is
severe, estimated covariate effects may be sensitive to modelling
choices (which covariates are included, how the observation model is
specified, what priors are used). Different defensible choices yield
different conclusions.

**Spatial smoothing that masks local heterogeneity.** If the dependence
structure assumes smooth spatial variation, local discontinuities or
heterogeneity will be smoothed away. The model will produce smooth maps
even when the underlying process is patchy.

**Apparent covariate effects that reflect spatial confounding.** If a
covariate varies spatially and the dependence structure is misspecified,
the covariate may absorb variation that should be attributed to spatial
structure, or vice versa. The estimated effect conflates multiple
sources.

**Posterior multimodality or poor mixing.** When ridges or flat regions
are present, sampling algorithms may fail to explore the posterior.
Apparent convergence may mask unvisited regions of comparable posterior
density.

These failure modes are not always diagnosable from standard outputs. A
model may fit well, produce smooth predictions, and yield narrow
uncertainty intervals while suffering from one or more of these issues.
The failures are failures of interpretation and generalisation, not of
fit.

## What this means for use

This chapter has catalogued limitations, not solutions. The purpose:
establish what the framework can and cannot deliver.

SIMP is informative about:

- Relative latent values among present species within plots.

- Overall patterns of presence across species and plots, given
  assumptions about presence thresholds.

- Combined effects of covariates on cover, integrating presence and
  allocation.

- Broad patterns of spatial or group-level similarity, given assumptions
  about dependence structure.

SIMP is weakly informative about:

- Absolute latent values (resolved only by constraints).

- Latent values for absent species (resolved primarily by priors).

- Separate contributions of presence and allocation to covariate
  effects.

- Fine-scale spatial structure when data are sparse.

SIMP is silent about:

- Whether the dependence structure is correctly specified.

- Whether apparent covariate effects reflect confounding with unmeasured
  spatial structure.

- Whether apparent precision reflects data or prior assumptions.

A user who understands these boundaries can use the framework
appropriately, interpret outputs cautiously, and recognise when
conclusions are robust versus fragile. A user who does not understand
them risks mistaking assumption-driven results for data-driven
conclusions.

## Summary

Identifiability in SIMP is partial and uneven. Some components of the
latent structure are well-constrained by data; others are weakly
informed; still others are determined entirely by assumptions.

The scale of latent values is not identified; fixed by convention. The
latent values of absent species are weakly identified; determined
primarily by priors. The separation of presence and allocation effects
is not identified when both processes share covariates. The dependence
structure confounds with fixed effects when covariates vary spatially.

These limitations are structural. They arise from the geometry of
constrained observations, the censoring induced by zeros, and the
architecture of hierarchical models. They cannot be eliminated by better
algorithms, larger samples, or more careful implementation.

The appropriate response is transparency. Users should understand which
aspects of their inference are data-driven and which are
assumption-driven. Sensitivity analyses should probe the assumptions
that matter. Conclusions should be calibrated to the strength of the
evidence, not the precision of the output.

This chapter does not propose solutions to these limitations. It
establishes that they exist, that they are intrinsic to the framework,
and that ignoring them leads to misuse. A framework that is honest about
its limits is more trustworthy than one that obscures them.
