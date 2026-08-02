# Convert between survival analysis measures

Converts a survival probability (and optionally a hazard) to other
survival analysis metrics.

## Usage

``` r
convert_survival_measures(survival, hazard = NULL, measure)
```

## Arguments

- survival:

  Numeric survival probability or vector of probabilities.

- hazard:

  Numeric hazard value or vector. Required for `"hazard"` and
  `"density"` measures.

- measure:

  Character string specifying the desired measure. One of: `"survival"`,
  `"risk"`, `"cumulative_hazard"`, `"hazard"`, or `"density"`.

## Value

Numeric value or vector of the requested measure.

## Examples

``` r
# Risk is the complement of survival
convert_survival_measures(0.8, measure = "risk")
#> [1] 0.2

# The density is the product of the hazard and the survival, so the
# "density" and "hazard" measures need the hazard as well.
convert_survival_measures(
  c(0.9, 0.8, 0.6),
  hazard = c(0.05, 0.06, 0.08),
  measure = "density"
)
#> [1] 0.045 0.048 0.048
```
