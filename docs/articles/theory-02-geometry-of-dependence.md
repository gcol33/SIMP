# The Geometry of Dependence

## The space where cover data live

Cover data do not live in ordinary space. When we record the proportions
of five species in a plot, we do not have five independent numbers that
could each take any value. We have five numbers that must be
non-negative and must sum to one. This constraint defines a specific
region: a four-dimensional surface embedded in five-dimensional space.

This is not a technicality. It changes what the data are.

In unconstrained space, each dimension is independent of the others in a
structural sense. A point can move freely along any axis without
affecting its position on the other axes. The axes are orthogonal not
just geometrically but logically: knowing the value on one axis tells
you nothing about possible values on the others.

On a constrained surface, this independence disappears. A point cannot
move freely in all directions. Movement in one direction forces
compensating movement in others. The dimensions are coupled by the
constraint itself, before any data are observed and before any
biological process is invoked.

## Why constrained spaces induce dependence

Consider a plot with five cover categories: grasses, forbs, shrubs, bare
ground, and litter. The five proportions must sum to one. Suppose
grasses currently occupy 40% of the plot. What can we say about the
other categories?

We know that the remaining four categories together occupy exactly 60%.
If grasses increase to 50%, the total available to all other categories
drops to 50%. If grasses occupy 80%, then forbs, shrubs, bare ground,
and litter together can occupy at most 20%—regardless of their
individual tendencies or how they were measured.

This is not a correlation we estimated from data. It is not a biological
interaction between plant functional types. It is a logical consequence
of the constraint. The cover of grasses restricts the possible covers of
every other category, and vice versa, because all draw from the same
finite total.

The same logic applies to any pair of components in any compositional
dataset. Knowing the value of one component constrains the possible
values of all others. This induced dependence is present regardless of
how the data were generated, how many observations we have, or whether
the components interact in any meaningful biological sense.

## The impossibility of independent components

In unconstrained settings, it is meaningful to ask whether two variables
are independent. We can imagine a data-generating process where each
variable is determined by its own separate mechanism, with no connection
between them. Independence is a coherent hypothesis that could be true
or false.

For components of a composition, this question is malformed.
Independence would require that knowing the value of one component tells
us nothing about the possible values of another. But the sum constraint
guarantees that it does. If species A occupies 80% of a plot, then all
other species together occupy exactly 20%. This is not information we
infer; it is information the constraint imposes.

One might object: perhaps independence holds in some underlying latent
space, before the constraint is applied. This is a reasonable line of
inquiry, and later chapters will examine what can and cannot be said
about latent structure. But it does not rescue component-wise thinking
at the level of observation. The data we observe are the constrained
proportions, not the latent quantities. Whatever independence might hold
in an underlying unconstrained space, the observed components remain
dependent. The constraint transforms independent latent variation into
dependent observed variation.

This is why “testing for dependence” among cover components is
conceptually confused. The dependence is guaranteed. The question is not
whether it exists, but what form it takes beyond the minimum induced by
the constraint.

## Negative correlations without interaction

A striking consequence of compositional structure is that components
will tend to be negatively correlated, even when there is no process
driving them apart.

Consider the simplest case: a plot with only two categories, species A
and bare ground. These two covers must sum to one, so they are perfectly
negatively correlated. If species A has high cover, bare ground must be
low. This correlation is not estimated; it is determined.

With more categories, the pattern is less rigid but persists. If we
generate compositional data by any reasonable process and compute
correlations among components, we will find negative correlations
dominating. This is not because the components are competing. It is
because they are drawing from the same finite pool.

A researcher who observes negative correlations among species covers
might interpret this as evidence of competitive exclusion. In some
cases, it may be. But the compositional structure alone guarantees some
negative correlation. Separating the geometric effect from the
biological effect requires understanding what dependence the constraint
induces—and this is not possible if the constraint is ignored.

## Why marginal thinking fails

Much of applied statistics operates marginally: we consider one variable
at a time, or pairs of variables, and build up an understanding from
these components. For compositional data, this approach is unreliable.

The marginal distribution of a single cover component—say, species A
considered in isolation—does not tell us how species A behaves in the
composition. Species A might have high variance marginally but low
variance conditionally on the other components. It might appear
symmetric marginally but be highly skewed in context. The marginal view
is incomplete.

Pairwise relationships are equally misleading. The correlation between
species A and species B depends on what else is in the composition. Add
a third species with high and variable cover, and the correlation
between A and B will change—not because their relationship changed, but
because the geometric context changed. Pairwise correlations are not
stable properties of the components; they are artifacts of the full
compositional structure.

This is why standard multivariate methods, which often assume that
marginal or pairwise information is sufficient, cannot be trusted for
compositional data. The constraint couples all components
simultaneously. Understanding one or two at a time misses the structure.

## The constraint is not a nuisance

It is tempting to treat the sum constraint as an inconvenience—something
that complicates analysis but can be managed with appropriate
adjustments. This view is mistaken.

The constraint is not separate from the data. It defines what the data
are. A vegetation plot with 30% grass, 25% forbs, 20% shrubs, and 25%
bare ground is a single object: a partition of space. The four numbers
are not four independent measurements that happen to sum to one. They
are four aspects of one thing.

Attempting to remove the constraint—by dropping one component, by using
ratios, by transforming to an unconstrained space—does not make the
dependence disappear. It relocates the dependence into a form that may
be harder to see but is no less present. The structure of the
constrained space follows the data into whatever representation we
choose.

## What this means

Dependence in cover data is not a hypothesis to be tested or an effect
to be estimated. It is a structural fact imposed by the nature of the
measurement. Any two components of a composition are dependent,
necessarily, because they share a common bound.

This dependence is geometric, not biological. It arises from the
constraint, not from species interactions, environmental filtering, or
ecological dynamics. These biological processes may add further
structure, but they operate on top of a dependence that is already
guaranteed.

Ignoring this geometry does not make it go away. It makes the resulting
analysis incoherent. A model that treats cover components as independent
has not made a simplifying assumption. It has made a false statement
about what the data are.

This vignette does not propose a solution to this problem. It
establishes that the problem exists and cannot be avoided. Any serious
treatment of cover data must begin by accepting that the components are
coupled, and that this coupling is not optional.
