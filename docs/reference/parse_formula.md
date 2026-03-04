# Parse a SIMP model formula

Extracts fixed effects, random effects, and response information from a
model formula. Supports all standard random effect structures.

## Usage

``` r
parse_formula(formula, data = NULL)
```

## Arguments

- formula:

  A model formula

- data:

  Optional data frame for validation

## Value

A list with parsed formula components

## Examples

``` r
# Fixed effects only
parse_formula(y ~ x1 + x2)

# Random intercepts
parse_formula(y ~ x1 + (1 | group))

# Random slopes
parse_formula(y ~ x1 + (1 + x1 | group))

# Nested random effects
parse_formula(y ~ x1 + (1 | site/plot))

# Crossed random effects
parse_formula(y ~ x1 + (1 | site) + (1 | year))

# Uncorrelated random effects
parse_formula(y ~ x1 + (1 + x1 || group))
```
