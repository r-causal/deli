# Digamma function

Evaluates the digamma function, the first derivative of the log gamma
function.

This is the equivalent of
[`base::digamma()`](https://rdrr.io/r/base/Special.html), and unlike
most of the deli utilities it is not required anywhere:
[`digamma()`](https://rdrr.io/r/base/Special.html) belongs to the `Math`
group generic and has a tangent rule, so it differentiates under
`deriv_method = "exact"` just as `deli_digamma()` does. The two agree
value for value, including at the poles and on missing input. The only
difference is that `deli_digamma()` returns `NaN` at the non-positive
integers, where the digamma function has poles, without base R's
`NaNs produced` warning. `deli_digamma()` exists so that code translated
from Python delicatessen, where the function is called
[`digamma()`](https://rdrr.io/r/base/Special.html), can keep its shape.

## Usage

``` r
deli_digamma(z)
```

## Arguments

- z:

  Numeric value or vector.

## Value

Numeric digamma values.

## See also

[`base::digamma()`](https://rdrr.io/r/base/Special.html), which is
usable in the same positions, and
[`deli_polygamma()`](https://r-causal.github.io/deli/reference/deli_polygamma.md)
for the higher-order derivatives.
[`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md)
has the full list of deli utilities and their base R counterparts,
including the five that base R cannot stand in for under exact
differentiation.

## Examples

``` r
deli_digamma(1)
#> [1] -0.5772157
deli_digamma(c(0.5, 1, 2))
#> [1] -1.9635100 -0.5772157  0.4227843

# The same values as the base R counterpart
all.equal(deli_digamma(c(0.5, 1, 2)), digamma(c(0.5, 1, 2)))
#> [1] TRUE

# NaN at the poles, where base R warns, and missing values pass through
deli_digamma(c(-1, 0, 2.5, NA))
#> [1]       NaN       NaN 0.7031566        NA
```
