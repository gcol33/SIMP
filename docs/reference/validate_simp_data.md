# Validate that data is suitable for SIMP modeling

Performs comprehensive checks beyond basic composition validity.

## Usage

``` r
validate_simp_data(x, min_obs = 10, min_components = 2, max_zero_pct = 95)
```

## Arguments

- x:

  Data to validate

- min_obs:

  Minimum observations required (default 10)

- min_components:

  Minimum components required (default 2)

- max_zero_pct:

  Maximum allowed zero percentage (default 95)

## Value

TRUE invisibly if valid, error otherwise
