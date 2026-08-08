# Display methods for deli estimators

Methods for [`base::print()`](https://rdrr.io/r/base/print.html) and
[`base::summary()`](https://rdrr.io/r/base/summary.html) so that a
fitted estimator shows its coefficients at the console and reports its
parameter-level inference as a table.

## Usage

``` r
# S3 method for class '`deli::deli_estimator`'
print(x, ..., subset = NULL)

# S3 method for class '`deli::deli_estimator`'
summary(object, alpha = 0.05, subset = NULL, ...)
```

## Arguments

- x, object:

  A fitted `MEstimator` or `GMMEstimator` object. Named `x` for
  [`print()`](https://rdrr.io/r/base/print.html) and `object` for
  [`summary()`](https://rdrr.io/r/base/summary.html), because
  [`base::print()`](https://rdrr.io/r/base/print.html) and
  [`base::summary()`](https://rdrr.io/r/base/summary.html) name their
  first argument that.

- ...:

  Not used. Must be empty, so that a name neither method recognizes is
  an error rather than silently ignored: a misspelled `alpha` would
  report limits at the default width and a misspelled `subset` would
  display every parameter, each while the call still read as the one
  that was meant.

- subset:

  Integer vector of parameter indices to display, or `NULL` (default) to
  display all of them.

- alpha:

  Numeric significance level for the confidence limits reported by
  [`summary()`](https://rdrr.io/r/base/summary.html), between 0 and 1.
  Default `0.05` for 95% limits.

## Value

- [`print()`](https://rdrr.io/r/base/print.html): its input, invisibly.

- [`summary()`](https://rdrr.io/r/base/summary.html): an object holding
  the estimates, standard errors, Z-scores, confidence limits, P-values,
  and S-values, which prints as a table.

## Details

[`print()`](https://rdrr.io/r/base/print.html) reports the parameter and
observation counts followed by the estimates, rounded to four decimal
places. An estimator that has not been through
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
has no estimates to show, so it reports its parameter count and says so.

[`summary()`](https://rdrr.io/r/base/summary.html) collects the
estimates, their standard errors, Z-scores, confidence limits, P-values,
and S-values into one object, which prints as a table with one row per
parameter. The values are those
[`confint()`](https://r-causal.github.io/deli/reference/deli-generics.md),
[`z_scores()`](https://r-causal.github.io/deli/reference/z_scores.md),
[`p_values()`](https://r-causal.github.io/deli/reference/p_values.md),
and
[`s_values()`](https://r-causal.github.io/deli/reference/s_values.md)
return individually, so `alpha` here means what it means there: `0.05`
gives 95% limits.

An over-identified `GMMEstimator` fit reports Hansen's J-statistic above
the table, with its degrees of freedom and its P-value, because it
judges the fit as a whole rather than any one parameter. See
[`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
for what it means and where its reference distribution holds. No other
fit has one, so no other output carries the line.

The `S` column reads `Inf` when a P-value underflows to exactly zero,
for the reason
[`s_values()`](https://r-causal.github.io/deli/reference/s_values.md)
gives. The `P` column reports that same underflow as `<2e-16`:
[`base::format.pval()`](https://rdrr.io/r/base/format.pval.html) stops
printing digits below the `eps` it is given, and the table gives it
`2.2e-16`.
[`tidy()`](https://r-causal.github.io/deli/reference/deli-tidiers.md)
returns the literal `0` instead.

`subset` restricts which parameters are displayed and nothing else. The
reported parameter count and every reported value are computed from the
whole fit, so displaying a subset of a stacked estimator is a way to
read the parameters of interest without the nuisance parameters, not a
way to refit without them. Row labels keep the names the full fit gave
them.

## See also

[deli-generics](https://r-causal.github.io/deli/reference/deli-generics.md)
for the accessors that return these quantities as plain vectors and
matrices, and
[deli-tidiers](https://r-causal.github.io/deli/reference/deli-tidiers.md)
for the same results as a data frame.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

fit
#> <MEstimator>
#>   Parameters: 3
#>   Observations: 32
#> Coefficients:
#> (Intercept): 37.2273
#> wt: -3.8778
#> hp: -0.0318

# `subset` is an argument of the method rather than of the fit, so showing it
# here means calling the generic by name.
print(fit, subset = 2:3)
#> <MEstimator>
#>   Parameters: 3
#>   Observations: 32
#> Coefficients:
#> wt: -3.8778
#> hp: -0.0318

summary(fit)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   270.5102
#> wt             -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    31.2310
#> hp             -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    19.1270

# `summary()` takes `subset` as well, alongside `alpha` for the width of the
# reported limits.
summary(fit, subset = 2:3, alpha = 0.1)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>         Estimate    Std.Err    Z-score    90% LCL    90% UCL    P-value    S-value
#> wt       -3.8778     0.6199    -6.2553    -4.8975    -2.8581   3.97e-10    31.2310
#> hp       -0.0318     0.0066    -4.7807    -0.0427    -0.0208   1.75e-06    19.1270
```
