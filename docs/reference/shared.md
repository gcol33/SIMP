# Specify a shared effect

Mark a predictor as having a shared coefficient across all components.
All components of the compositional response will use the same
coefficient for this predictor.

## Usage

``` r
shared(x)
```

## Arguments

- x:

  A variable name or expression

## Value

A term marker for formula parsing

## Examples

``` r
# All components share the same temperature coefficient
# cover ~ shared(temp) + (1 | site)

# Mix shared and varying effects
# cover ~ shared(temp) + varying(precip) + (1 | site)
```
