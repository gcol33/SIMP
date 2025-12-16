# The Nature of Cover Data

## What SIMP means by cover

Vegetation cover can mean different things depending on the measurement
protocol. Some protocols measure projected canopy area, allowing overlap
within a single layer and summing well above 100%. Others record
multiple vertical strata separately, each with its own cover values.
Still others estimate cover using ordinal scales (such as Braun-Blanquet
classes) whose midpoint conversions may not sum coherently.

SIMP adopts a specific definition: **cover as partition**. Under this
definition, a plot is exhaustively divided among mutually exclusive
categories. Each point in the plot belongs to exactly one category: a
species, bare ground, litter, or another residual class. The covers of
all categories sum to exactly one (or 100%). No overlap, no gaps.

This is a statement about how the data are interpreted, not about how
they were originally collected. Many datasets can be converted to this
format; some cannot. The framework requires the conversion to happen
before modelling, with the choices made explicit.

## Why this definition?

The partition definition is mathematically tractable. It places data on
the simplex, a well-understood geometric object with known properties.
The constraints that follow (non-negativity, unit sum) are structural
features that can be modelled coherently.

Other definitions of cover create different mathematical problems. If
covers can exceed 100% within a layer, the data live in the positive
orthant rather than the simplex. The geometry changes. The dependence
structure changes. The observation model changes. These are not small
adjustments; they require a different framework.

SIMP does one thing well: compositional data on the simplex with
structural zeros. Users whose data do not fit this structure should
either convert them (with appropriate functions) or use a different
approach.

## Scope beyond ecology

The mathematical structure of cover-as-partition appears in other
domains: market shares, budget allocations, time-use surveys, species
relative abundances in microbiome data, geochemical compositions. In
each case, the quantity represents a proportion of a finite whole, the
parts must sum to one, and zeros may indicate genuine absence rather
than small values.

The theory developed here applies to any such data. The examples and
language are ecological, but the framework is not ecology-specific. A
reader from another field can substitute their own domain; the
mathematics is identical.

## When data do not fit

Real vegetation data often arrive in forms that violate the partition
assumption:

**Covers summing above 100%.** This occurs when canopy overlap is
recorded within a layer, or when ordinal class midpoints are summed. The
data must be rescaled to sum to one before modelling. SIMP provides
`normalize_cover()` for this purpose.

**Multiple vertical layers.** Tree, shrub, and herb layers may each have
separate cover estimates. These must be collapsed into a single
composition. SIMP provides `collapse_layers()` with options for how to
combine them.

**Missing residual category.** If the recorded species do not exhaust
the plot (no bare ground or litter category), a residual must be added
to close the composition. SIMP provides `add_residual()` for this.

**Covers summing below 100% for unclear reasons.** This may indicate
missing data, recording errors, or a protocol that does not aim for
exhaustive coverage. The user must decide how to handle it; SIMP
provides `validate_composition()` to diagnose the issue.

Each of these functions makes explicit a choice that would otherwise be
hidden. The user sees what transformation was applied and can assess
whether it is appropriate for their data and question.

## The three defining constraints

Under the partition definition, cover data have three structural
properties.

**Boundedness.** Cover values live on a closed interval. For a single
category, this interval is \[0, 1\] when expressed as a proportion. A
value of 0.4 means something; a value of 1.3 does not. The bound is not
a statistical convenience. It is the edge of the quantity’s existence.

**Compositionality.** When multiple categories share a plot, their
covers are not independent. Total cover must sum to exactly one. If
species A occupies 60% of a plot, all other categories together occupy
exactly 40%. This follows from the definition of partition.

**Prevalence of structural zeros.** In most vegetation surveys, many
species are absent from most plots. These zeros are exact: the species
is not there. A sharp distinction between absence (cover = 0) and
presence (cover \> 0) that cannot be smoothed away.

These three properties interact. The bounds constrain the sum. The sum
constrains individual values. The zeros are structurally different from
small positive values. Any approach that addresses one property while
ignoring the others has misunderstood the data.

## Why transformation does not help

A common response to bounded data is to transform them onto an unbounded
scale, fit a standard model, and back-transform the results. For
compositional cover data, this strategy fails at a structural level.

Consider the logit transformation, which maps (0, 1) to the real line.
The transformation is undefined at 0 and 1. Since cover data routinely
include exact zeros (and sometimes exact ones, for dominant species),
the transformation cannot be applied without first modifying the data.
Common modifications include adding a small constant, replacing zeros
with half the detection limit, or treating zeros as a separate category.
Each modification changes what the data say. The zeros are no longer
exact. The bounds are no longer hard. The data have been altered to fit
the method, rather than the method chosen to fit the data.

The problem runs deeper than undefined values. The logit transformation
assumes that the quantity being transformed is a probability or
proportion that could, in principle, take any value in (0, 1). Cover
data do not behave this way. A species with 0.1% cover is present. Its
cover could increase up to 100%. The geometry near zero differs from the
geometry near the middle of the interval.

Compositional constraints create additional problems. If species covers
must sum to exactly one, any model that treats their transformed values
as independent has imposed an incorrect independence structure. The
constraint will be violated in expectation, in variance, or both.
Transforming each species separately and hoping the results combine
sensibly is not a modelling strategy. It is evasion.

## The structural nature of the failure

These are category errors, not approximation errors that shrink with
large samples.

When a model places positive probability on impossible values, it
describes a different quantity than the one measured. Predictions of
negative cover or cover exceeding 100% are meaningless. A model that
produces them, even occasionally, has failed to represent the
data-generating process.

When a model treats exact zeros as inconvenient, it refuses to engage
with a central feature of the data. In most vegetation datasets, zeros
carry as much information as positive values. A species being absent
from a plot is not a gap in the data. It is the data.

When a model ignores compositional constraints, it assumes that species
occupy space independently of one another. False. Space is finite. What
one species occupies, another cannot. Not a subtle biological
interaction that might or might not hold. Geometry.

## What this means for modelling

The constraints of cover data define the space in which the data exist.
Any model must respect this space as a foundational commitment, not as a
post-hoc correction.

This does not mean modelling cover data is impossible, or that no useful
work has been done. The starting point matters. A model built on the
assumption that cover is an unbounded continuous variable, later
constrained through truncation or transformation, begins from a false
premise. The constraints will reassert themselves as inconsistencies,
boundary effects, and predictions that make no sense.

The alternative is to begin with the constraints. To treat the bounded,
compositional, zero-structured nature of cover as the defining feature
of the problem.

This vignette makes no claim about how to do this. Only that it must be
done.
