# Augment data with predictions from a fitted deli estimator

A
[`generics::augment()`](https://generics.r-lib.org/reference/augment.html)
method for `MEstimator` and `GMMEstimator` objects fitted through the
formula interface. It returns the model frame the fit was built from, or
`newdata` when supplied, with the fitted values, their standard errors,
a Wald confidence interval, and the residuals as columns beside it.

## Arguments

- x:

  A fitted `MEstimator` or `GMMEstimator` object made with the formula
  interface (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- newdata:

  A data frame of covariate values to predict at, or `NULL` (default) to
  report the model frame the fit was built from. That frame holds the
  variables the formula named, so a transformed term appears as the
  column the transformation produced and a column of the fitting data
  the formula did not name is not there.

- type.predict:

  Character string. `"link"` (default) puts `.fitted` and its interval
  on the scale of the linear predictor; `"response"` puts them on the
  scale of the response.

- conf.level:

  Numeric confidence level for `.lower` and `.upper`. Default `0.95`.

- ...:

  Not used. Must be empty, so that a name that is not one of the
  documented arguments is an error rather than silently ignored. A
  misspelled `newdata` would otherwise augment the fitted rows while the
  caller believed they had asked for rows of their own.

## Value

A data frame: the model frame, or `newdata`, followed by the columns
`.fitted`, `.se.fit`, `.lower`, `.upper`, and `.resid`. The last is
absent when `newdata` is supplied.

## Details

The columns added are `.fitted`, `.se.fit`, `.lower`, `.upper`, and,
when `newdata` is not supplied, `.resid`. They are exactly what
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md)
and
[`residuals()`](https://r-causal.github.io/deli/reference/deli-generics.md)
return for the same fit, so `.fitted` is
[`predict()`](https://rdrr.io/r/stats/predict.html), `.lower` and
`.upper` are `predict(interval = "confidence")`, and `.resid` is
[`residuals()`](https://rdrr.io/r/stats/residuals.html).

`.fitted` is on the link scale by default, matching both
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md)
and `broom::augment()` on a `glm`, and `type.predict = "response"` puts
it and its interval on the scale of the response. `.resid` is the
response residual either way, since a residual measured against a linear
predictor would be a response minus a quantity the response is not
measured in.

Rows the fit dropped for missing data are not reported, so the result
has one row per
[`nobs()`](https://r-causal.github.io/deli/reference/deli-generics.md)
and its row names are those of the retained rows. A `newdata` row with a
missing value is kept, with `NA` in the added columns, so that the
result lines up with the rows handed in.

`augment()` covers the estimating equations whose linear predictor
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md)
forms, and refuses the same fits with the same reasons; see
[deli-predict](https://r-causal.github.io/deli/reference/deli-predict.md).
It has no counterpart to that method's `times` argument: a survival
measure is one value per row of the data and time rather than one per
row, so the predictions do not go beside the data as columns.

## See also

[deli-predict](https://r-causal.github.io/deli/reference/deli-predict.md)
for the predictions themselves and the equations they are available for,
and
[deli-tidiers](https://r-causal.github.io/deli/reference/deli-tidiers.md)
for the parameter-level and model-level summaries.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

# Attaching generics or broom makes the shorter augment(fit) call work
# as well
head(generics::augment(fit))
#>                    mpg    wt  hp  .fitted   .se.fit   .lower   .upper
#> Mazda RX4         21.0 2.620 110 23.57233 0.6045332 22.38747 24.75719
#> Mazda RX4 Wag     21.0 2.875 110 22.58348 0.5531278 21.49937 23.66759
#> Datsun 710        22.8 2.320  93 25.27582 0.7364506 23.83240 26.71924
#> Hornet 4 Drive    21.4 3.215 110 21.26502 0.5516788 20.18375 22.34629
#> Hornet Sportabout 18.7 3.440 175 18.32727 0.4282785 17.48786 19.16668
#> Valiant           18.1 3.460 105 20.47382 0.6221296 19.25446 21.69317
#>                       .resid
#> Mazda RX4         -2.5723294
#> Mazda RX4 Wag     -1.5834826
#> Datsun 710        -2.4758187
#> Hornet 4 Drive     0.1349799
#> Hornet Sportabout  0.3727334
#> Valiant           -2.3738163

# New covariate patterns come back with no residual column, since they carry
# no response to residualize against.
generics::augment(fit, newdata = data.frame(wt = c(2, 3, 4), hp = 110))
#>   wt  hp  .fitted   .se.fit   .lower   .upper
#> 1  2 110 25.97658 0.8477261 24.31507 27.63810
#> 2  3 110 22.09875 0.5431699 21.03416 23.16335
#> 3  4 110 18.22092 0.8000295 16.65289 19.78895
```
