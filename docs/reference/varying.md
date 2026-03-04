# Specify a varying effect

Mark a predictor as having component-specific (varying) coefficients.
Each component of the compositional response will have its own
coefficient for this predictor.

## Usage

``` r
varying(x)
```

## Arguments

- x:

  A variable name or expression

## Value

A term marker for formula parsing

## Examples

``` r
# Each component gets its own temperature coefficient
# cover ~ varying(temp) + (1 | site)

# Equivalent to the default behavior for fixed effects
# cover ~ temp + (1 | site)
```
