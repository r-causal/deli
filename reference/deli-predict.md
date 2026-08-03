# Predictions from a fitted deli estimator

A method for [`stats::predict()`](https://rdrr.io/r/stats/predict.html)
that predicts from a fit made through the formula interface, with
sandwich standard errors and Wald confidence intervals. It forms the
linear predictor of a regression fit on either the link or the response
scale, and evaluates a survival measure at a set of `times` for a fit of
[`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) or
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md).

## Usage

``` r
# S3 method for class '`deli::deli_estimator`'
predict(object, newdata = NULL,
  type = c("link", "response"), se.fit = FALSE,
  interval = c("none", "confidence"), level = 0.95, times = NULL,
  measure = "survival", deriv_method = "capprox", dx = 1e-9, ...)
```

## Arguments

- object:

  A fitted `MEstimator` or `GMMEstimator` object made with the formula
  interface (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- newdata:

  A data frame of covariate values to predict at, or `NULL` (default) to
  predict the rows the fit was made on. An offset written into the
  formula with [`offset()`](https://rdrr.io/r/stats/offset.html) is
  evaluated on `newdata`; an offset supplied through `...` at fit time
  is one value per fitted observation, so predicting on `newdata` is an
  error rather than a prediction at an offset of zero.

- type:

  Character string. `"link"` (default) returns the linear predictor;
  `"response"` applies the inverse link, giving the conditional mean of
  the response. Cannot be supplied beside `times`.

- se.fit:

  Logical. Return standard errors beside the predictions? Default
  `FALSE`.

- interval:

  Character string. `"none"` (default) or `"confidence"` for Wald
  intervals at `level`.

- level:

  The confidence level for `interval`. Default `0.95`. The critical
  value comes from the standard normal distribution, or from the
  t-distribution with \\n - p\\ degrees of freedom when a
  `finite_correction` is set on the fit, matching
  [`confint()`](https://r-causal.github.io/deli/reference/deli-generics.md).

- times:

  Numeric vector of times to predict a survival measure at, or `NULL`
  (default) to predict the linear predictor instead. Supported for a fit
  of [`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md)
  or
  [`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md).

- measure:

  Character string naming the survival measure, one of `"survival"`
  (default), `"risk"`, `"cumulative_hazard"`, `"hazard"`, or
  `"density"`; see
  [`convert_survival_measures()`](https://r-causal.github.io/deli/reference/convert_survival_measures.md).
  Only meaningful beside `times`.

- deriv_method:

  Character string for the derivative method used to build the Jacobian
  of the survival measure. One of `"capprox"` (central difference, the
  default), `"fapprox"`, `"bapprox"`, or `"exact"` (forward-mode
  automatic differentiation, available for an AFT fit and an error for a
  pooled logistic one). Only meaningful beside `times`.

- dx:

  Numeric step size for the finite-difference methods, ignored when
  `deriv_method = "exact"`. Default `1e-9`. Only meaningful beside
  `times`. Must be a single positive finite number, which is checked
  wherever it is supplied beside `times`, including a prediction asking
  for no standard error, which takes no step at all.

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

Without `times`, a named numeric vector of predictions, one per row of
`newdata` or of the fitted design, and with `interval = "confidence"` a
matrix with columns `"fit"`, `"lwr"`, and `"upr"`.

With `times`, a data frame with one row per row of the design and time,
every time for the first row before any time for the second, and columns
`.row`, the row label of the design, `time`, and `fit`. With
`interval = "confidence"` it also has `"lwr"` and `"upr"`.

With `se.fit = TRUE`, either shape becomes a list whose `fit` element is
whichever of the two the other arguments call for and whose `se.fit`
element is a numeric vector of standard errors.

## Two surfaces

`times` chooses what is predicted. Without it,
[`predict()`](https://rdrr.io/r/stats/predict.html) returns the linear
predictor on the scale `type` names. With it,
[`predict()`](https://rdrr.io/r/stats/predict.html) returns the survival
measure `measure` names at each of the times, one prediction per row of
the design at each time.

The two cover disjoint sets of fits. An equation reaches the first
surface by having a linear predictor that is a conditional mean of the
response, and the two equations on the second surface are there because
theirs is not:
[`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) puts
its linear predictor on the log-time scale, and
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
has one linear predictor per person and time interval rather than one
per person. Supplying `type` beside `times`, or `measure`,
`deriv_method`, or `dx` without `times`, is an error rather than an
argument silently ignored, since no fit takes both sets.

## The linear predictor

The design for `newdata` is rebuilt through the terms, factor levels,
and contrasts the fit recorded, so a factor whose new values cover only
some of the fitted levels still produces the fitted set of columns, and
a data-dependent term such as `poly(x, 2)` or `scale(x)` is evaluated
with the coefficients it was fitted with rather than refitted to
`newdata`. A factor level the fit never saw, or a predictor `newdata`
does not carry, is an error.

The standard error is the delta-method standard error of the linear
predictor, \\\sqrt{\mathrm{diag}(X \hat{V} X^{T})}\\, formed from the
coefficient block of the sandwich variance. On the response scale it is
scaled by the derivative of the inverse link, which is exact because the
inverse link is applied elementwise. The intervals are Wald intervals on
the scale asked for, so they are symmetric about `fit` on that scale
rather than transformed from the link scale.

[`predict()`](https://rdrr.io/r/stats/predict.html) needs to know which
parameters are coefficients on the design and what takes the linear
predictor to the mean of the response. It supports
[`ee_regression()`](https://r-causal.github.io/deli/reference/ee_regression.md),
[`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md),
[`ee_robust_regression()`](https://r-causal.github.io/deli/reference/ee_robust_regression.md),
[`ee_beta_regression()`](https://r-causal.github.io/deli/reference/ee_beta_regression.md),
and the five penalized regressions
[`ee_bridge_regression()`](https://r-causal.github.io/deli/reference/ee_bridge_regression.md),
[`ee_ridge_regression()`](https://r-causal.github.io/deli/reference/ee_ridge_regression.md),
[`ee_lasso_regression()`](https://r-causal.github.io/deli/reference/ee_lasso_regression.md),
[`ee_dlasso_regression()`](https://r-causal.github.io/deli/reference/ee_dlasso_regression.md),
and
[`ee_elasticnet_regression()`](https://r-causal.github.io/deli/reference/ee_elasticnet_regression.md),
whose parameters are one coefficient per design column followed by at
most one parameter of the outcome distribution. Any other estimating
equation, and any fit built from a `stacked_equations` function, is an
error naming the reason;
[`regression_predictions()`](https://r-causal.github.io/deli/reference/regression_predictions.md)
takes a design, estimates, and a covariance matrix directly and can be
used wherever this method declines.

## Survival measures at a set of times

`times` is supported for a fit of
[`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) or of
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md),
and `measure` names any of the measures
[`convert_survival_measures()`](https://r-causal.github.io/deli/reference/convert_survival_measures.md)
defines. The point estimates are those of
[`aft_predictions_individual()`](https://r-causal.github.io/deli/reference/aft_predictions_individual.md)
and
[`plogit_predict()`](https://r-causal.github.io/deli/reference/plogit_predict.md)
for the same fit.

[`ee_survival_model()`](https://r-causal.github.io/deli/reference/ee_survival_model.md)
has no surface here because the formula interface cannot drive it: it
takes no design matrix, and the interface always passes the one it
built. Predict from such a fit with
[`survival_predictions()`](https://r-causal.github.io/deli/reference/survival_predictions.md).

Each covariate pattern has its own variance, so each row gets its own
interval. The variance of the measure for a row at a time is the
delta-method variance \\G \hat{V} G^{T}\\ of that one prediction, where
\\G\\ is its row of the Jacobian of the whole grid with respect to the
parameters, built by `deriv_method`. A pooled logistic fit is predicted
through
[`plogit_predict()`](https://r-causal.github.io/deli/reference/plogit_predict.md),
whose matrix products cannot carry tangents, so `deriv_method = "exact"`
is available for an AFT fit alone. Asking for it on a pooled logistic
fit is an error whether or not a standard error is wanted, since the
Jacobian that could not be built is only built when one is.

## See also

[deli-augment](https://r-causal.github.io/deli/reference/deli-augment.md),
which returns the linear-predictor predictions as columns beside the
data;
[`regression_predictions()`](https://r-causal.github.io/deli/reference/regression_predictions.md),
which computes the same quantities from a design matrix, a vector of
estimates, and a covariance matrix, for fits this method does not cover;
and
[`aft_predictions_individual()`](https://r-causal.github.io/deli/reference/aft_predictions_individual.md),
[`aft_predictions_function()`](https://r-causal.github.io/deli/reference/aft_predictions_function.md),
[`plogit_predict()`](https://r-causal.github.io/deli/reference/plogit_predict.md),
and
[`survival_predictions()`](https://r-causal.github.io/deli/reference/survival_predictions.md),
which compute the survival measures from estimates directly.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

head(predict(fit))
#>         Mazda RX4     Mazda RX4 Wag        Datsun 710    Hornet 4 Drive 
#>          23.57233          22.58348          25.27582          21.26502 
#> Hornet Sportabout           Valiant 
#>          18.32727          20.47382 

# New covariate patterns, with sandwich standard errors.
at <- data.frame(wt = c(2, 3, 4), hp = 110)

predict(fit, newdata = at, se.fit = TRUE)
#> $fit
#>        1        2        3 
#> 25.97658 22.09875 18.22092 
#> 
#> $se.fit
#>         1         2         3 
#> 0.8477261 0.5431699 0.8000295 
#> 

predict(fit, newdata = at, interval = "confidence")
#>        fit      lwr      upr
#> 1 25.97658 24.31507 27.63810
#> 2 22.09875 21.03416 23.16335
#> 3 18.22092 16.65289 19.78895

# A logistic fit predicts log-odds by default and probabilities on the
# response scale.
vs_fit <- m_estimate(vs ~ wt, data = mtcars, .ee = ee_regression,
                     model = "logistic")

predict(vs_fit, newdata = data.frame(wt = c(2, 3)), type = "response")
#>         1         2 
#> 0.8691639 0.4957597 

# A pooled logistic fit has no linear predictor to put on a scale, and
# predicts a survival measure at a set of times instead. Every row of the
# design gets its own interval.
bladder <- collett_bladder
bladder$novel <- bladder$treat - 1
k <- length(unique(bladder$time[bladder$delta == 1]))

plogit_fit <- m_estimate(
  time ~ novel + init + size - 1, data = bladder, .ee = ee_plogit,
  event = delta, init = c(rep(0, 3), -4, rep(0, k - 1))
)

head(predict(plogit_fit, times = c(12, 24), interval = "confidence"))
#>   .row time       fit       lwr       upr
#> 1    1   12 0.6310973 0.4870471 0.7751474
#> 2    1   24 0.5302241 0.3640916 0.6963566
#> 3    2   12 0.5882252 0.4497605 0.7266899
#> 4    2   24 0.4810942 0.3313945 0.6307940
#> 5    3   12 0.5534843 0.4042987 0.7026700
#> 6    3   24 0.4422455 0.2794547 0.6050364
```
