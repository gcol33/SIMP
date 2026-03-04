# Impute zeros in compositional data

Replace zeros with small positive values using various imputation
methods. This is necessary for log-ratio transformations which require
non-zero values.

## Usage

``` r
impute_zeros(
  x,
  method = c("multiplicative", "additive", "dl"),
  delta = NULL,
  preserve_structural = FALSE,
  structural_threshold = 0.8
)
```

## Arguments

- x:

  A simp_composition or simp_cover object

- method:

  Imputation method:

  - `"multiplicative"`: Multiplicative replacement (default)

  - `"additive"`: Simple additive replacement

  - `"dl"`: Detection limit replacement

- delta:

  Replacement value. If NULL, automatically determined based on method.
  For multiplicative, uses 0.65 times minimum non-zero value.

- preserve_structural:

  Logical; if TRUE, don't impute zeros that appear to be structural
  (default FALSE)

- structural_threshold:

  Threshold for structural zeros (default 0.8)

## Value

Object with imputed values

## References

Martin-Fernandez, J.A., Barcelo-Vidal, C., Pawlowsky-Glahn, V. (2003)
Dealing with zeros and missing values in compositional data sets using
nonparametric imputation. Mathematical Geology, 35(3), 253-278.

## Examples

``` r
x <- simp_composition(
  matrix(c(0.4, 0.4, 0.2, 0,
           0.5, 0.3, 0.2, 0,
           0.3, 0.5, 0, 0.2), nrow = 3, byrow = TRUE),
  components = c("sp1", "sp2", "sp3", "sp4")
)

# Multiplicative replacement
impute_zeros(x)

# Additive with custom delta
impute_zeros(x, method = "additive", delta = 0.001)
```
