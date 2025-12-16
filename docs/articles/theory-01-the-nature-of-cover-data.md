# The Nature of Cover Data

## What is vegetation cover?

Vegetation cover is a statement about space. When we say that a species
has 40% cover in a plot, we mean that 40% of the plot area is occupied
by that species. It is a direct claim about physical presence in a
bounded region.

This definition has consequences. Cover cannot be negative, because a
species cannot occupy less than no space. Cover cannot exceed 100%,
because a species cannot occupy more space than exists. These are
logical necessities that follow from what cover means.

The same structure appears in other domains: market shares, budget
allocations, time-use data, chemical concentrations in mixtures. In each
case, the quantity represents a proportion of a finite whole. The
finiteness defines the quantity.

## The three defining constraints

Cover data have three structural properties that any honest treatment
must respect.

**Boundedness.** Cover values live on a closed interval. For a single
species, this interval is \[0, 1\] when expressed as a proportion. A
value of 0.4 means something; a value of 1.3 does not. The bound is a
boundary beyond which the quantity ceases to exist.

**Compositionality.** When multiple species share a plot, their covers
are not independent. The total cover across all species must sum to
exactly one: the plot is partitioned, and every point in the plot
belongs to exactly one category (which may include bare ground or litter
as a residual class). If species A occupies 60% of a plot, then all
other categories together occupy exactly 40%. This is a consequence of
measuring area within a finite boundary.

**Prevalence of structural zeros.** In most vegetation surveys, many
species are absent from most plots. These zeros are exact: the species
is not there. This creates a sharp distinction between absence (cover =
0) and presence (cover \> 0) that cannot be smoothed away.

These three properties interact. The bounds constrain the sum. The sum
constrains individual values. The zeros are structurally different from
small positive values. Any approach that addresses one property while
ignoring the others has misunderstood the data.

## Why transformation does not help

A common response to bounded data is to transform them onto an unbounded
scale, fit a standard model, and back-transform the results. For cover
data, this strategy fails at a structural level.

Consider the logit transformation, which maps (0, 1) to the real line.
This transformation is undefined at 0 and 1. Since cover data routinely
include exact zeros (and sometimes exact ones, for dominant species),
the transformation cannot be applied without first modifying the data.
Common modifications include adding a small constant, replacing zeros
with half the detection limit, or treating zeros as a separate category.
Each of these modifications changes what the data say. The zeros are no
longer exact. The bounds are no longer hard. The data have been altered
to fit the method, rather than the method being chosen to fit the data.

The problem is deeper than undefined values. The logit transformation
assumes that the quantity being transformed is a probability or a
proportion that could, in principle, take any value in (0, 1). Cover
data do not behave this way. A species with 0.1% cover is present. Its
cover could increase without bound up to 100%. The geometry near zero is
fundamentally different from the geometry near the middle of the
interval.

Compositional constraints create additional problems. If species covers
must sum to exactly one, then any model that treats their transformed
values as independent has imposed an incorrect independence structure.
The constraint will be violated in expectation, in variance, or both.
Transforming each species separately and hoping the results combine
sensibly is not a modelling strategy. It is an evasion.

## The structural nature of the failure

These are category errors, not approximation errors that become
negligible with large samples.

When a model places positive probability on impossible values, it has
described a different quantity than the one measured. Predictions of
negative cover or cover exceeding 100% are meaningless. A model that
produces them, even occasionally, has failed to represent the
data-generating process.

When a model treats exact zeros as inconvenient, it has refused to
engage with a central feature of the data. In most vegetation datasets,
the zeros carry as much information as the positive values. A species
being absent from a plot is not a gap in the data. It is the data.

When a model ignores compositional constraints, it has assumed that
species occupy space independently of one another. This is false. Space
is finite. What one species occupies, another cannot. This is not a
subtle biological interaction that might or might not hold. It is
geometry.

## What this means for modelling

The constraints of cover data define the space in which the data exist.
Any model for cover data must respect this space as a foundational
commitment, not as a post-hoc correction.

This does not mean that modelling cover data is impossible, or that no
useful work has been done. It means that the starting point matters. A
model built on the assumption that cover is an unbounded continuous
variable, later constrained through truncation or transformation, begins
from a false premise. The constraints will reassert themselves as
inconsistencies, boundary effects, and predictions that do not make
sense.

The alternative is to begin with the constraints. To treat the bounded,
compositional, zero-structured nature of cover as the defining feature
of the problem.

This vignette makes no claim about how to do this. It claims only that
it must be done.
