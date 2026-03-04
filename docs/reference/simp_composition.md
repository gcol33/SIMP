# Create a composition object

Constructor for the `simp_composition` class, representing data where
each row sums to 1 and all values are in \[0, 1\].

## Usage

``` r
simp_composition(x, components = NULL, validate = TRUE, normalize = FALSE)
```

## Arguments

- x:

  A matrix or data frame of compositional data

- components:

  Character vector of component names. If NULL, uses column names from x
  or generates default names.

- validate:

  Logical; validate the composition (default TRUE)

- normalize:

  Logical; normalize rows to sum to 1 if needed (default FALSE)

## Value

A `simp_composition` object

## Examples

``` r
# From matrix
x <- matrix(c(0.3, 0.5, 0.2,
              0.4, 0.4, 0.2,
              0.5, 0.3, 0.2), nrow = 3, byrow = TRUE)
comp <- simp_composition(x, components = c("sp1", "sp2", "sp3"))

# From data frame
df <- data.frame(grass = c(0.3, 0.4), shrub = c(0.5, 0.4), tree = c(0.2, 0.2))
comp <- simp_composition(df)
```
