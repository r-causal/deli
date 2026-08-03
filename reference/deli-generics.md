# Standard S3 generics for deli estimators

Methods for base R generics
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
[`stats::confint()`](https://rdrr.io/r/stats/confint.html),
[`stats::nobs()`](https://rdrr.io/r/stats/nobs.html),
[`stats::df.residual()`](https://rdrr.io/r/stats/df.residual.html),
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`stats::residuals()`](https://rdrr.io/r/stats/residuals.html),
[`stats::weights()`](https://rdrr.io/r/stats/weights.html),
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html),
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
[`stats::formula()`](https://rdrr.io/r/stats/formula.html),
[`stats::terms()`](https://rdrr.io/r/stats/terms.html),
[`stats::sigma()`](https://rdrr.io/r/stats/sigma.html),
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html), and
[`stats::deviance()`](https://rdrr.io/r/stats/deviance.html) so that
deli estimator objects interoperate with the broader R modeling
ecosystem.

## Arguments

- object:

  A fitted `MEstimator` or `GMMEstimator` object (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  The confidence level required. Default `0.95`.

- type:

  Character string naming the residual
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) returns. Only
  `"response"` is available, since it is the one residual defined for
  every estimating equation
  [`predict()`](https://rdrr.io/r/stats/predict.html) supports. Any
  other value is an error rather than a response residual under another
  name.

- x:

  A fitted `MEstimator` or `GMMEstimator` object. Named `x` because
  [`stats::formula()`](https://rdrr.io/r/stats/formula.html) and
  [`stats::terms()`](https://rdrr.io/r/stats/terms.html) name their
  first argument that.

- formula:

  A fitted `MEstimator` or `GMMEstimator` object. Named `formula`
  because
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  names its first argument that; the method takes a fit, not a formula.

- data:

  A data frame to build the model frame or the design matrix of, or
  `NULL` (default) to report the one the fit was solved on.
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) needs it
  to carry the response as well as the predictors;
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) names no
  response, so the predictors alone are enough there.

- subset, na.action, drop.unused.levels, xlev:

  Passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html),
  and meaning there what they mean for any other model frame. `xlev`
  defaults to the factor levels the fit recorded. All four describe how
  to build a frame from `data`, so supplying one without `data` is an
  error rather than a silent no-op.

- ...:

  Not used. [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html), and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) require
  them to be empty, so that a name none of them recognizes is an error
  rather than silently ignored: a misspelled `level` would report the
  default limits, a misspelled `type` would report response residuals as
  though they had been asked for, and a misspelled `data` would report
  on the fitted rows while the caller believed they had asked for rows
  of their own. The rest ignore them, which is the convention for a base
  generic with no optional argument for a wrong name to displace.

## Value

- [`coef()`](https://rdrr.io/r/stats/coef.html): Named numeric vector of
  parameter estimates.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html): Named
  variance-covariance matrix.

- [`confint()`](https://rdrr.io/r/stats/confint.html): Matrix with
  columns `"lower"` and `"upper"`.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html): Integer number of
  observations.

- [`df.residual()`](https://rdrr.io/r/stats/df.residual.html): Integer
  residual degrees of freedom, the number of observations less the
  number of parameters.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): Named
  numeric vector of fitted values on the response scale, one per
  observation the fit was solved on.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): Named numeric
  vector of response residuals.

- [`weights()`](https://rdrr.io/r/stats/weights.html): The observation
  weights the fit was solved with, as they were recorded, or `NULL` for
  a formula fit specified without any. Weights that vary over time are
  reported as the matrix they were supplied as rather than as a vector:
  [`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
  takes an n-by-K matrix carrying one column per time interval, and an
  observation weighted differently in each interval has no one weight
  for a vector to hold.

- [`model.frame()`](https://rdrr.io/r/stats/model.frame.html): The model
  frame the fit was built from, with the rows dropped for missing data
  already removed, or the model frame of `data` when one is supplied.

- [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html): The
  design matrix the fit was solved on, or the design matrix of `data`
  when one is supplied, coded with the contrasts and factor levels the
  fit recorded.

- [`formula()`](https://rdrr.io/r/stats/formula.html): The model
  formula.

- [`terms()`](https://rdrr.io/r/stats/terms.html): The `terms` object of
  the model frame, carrying the response index, any offset, and the
  `predvars` of a data-dependent term.

- [`sigma()`](https://rdrr.io/r/stats/sigma.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html): Nothing. Each
  raises an error saying that an M-estimator states no likelihood and
  records no residual scale.

## Details

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`weights()`](https://rdrr.io/r/stats/weights.html),
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
[`formula()`](https://rdrr.io/r/stats/formula.html), and
[`terms()`](https://rdrr.io/r/stats/terms.html) report on the model a
fit was specified as, which only the formula interface records, so each
of them is an error for a fit built from a `stacked_equations` function.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are also an
error for a formula fit of an estimating equation
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md)
does not support, since a fitted value is a prediction; see
[deli-predict](https://r-causal.github.io/deli/reference/deli-predict.md)
for the equations it covers.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is the
conditional mean of the response, so it is on the response scale rather
than the link scale, matching
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a
`glm` object. [`residuals()`](https://rdrr.io/r/stats/residuals.html) is
the response minus that mean, the residual a GLM calls a response
residual, and the only `type` it takes. Deviance and Pearson residuals
are not offered: the first needs a likelihood and the second a variance
function, and an M-estimator need not state either, while the response
residual is defined the same way for every equation
[`predict()`](https://rdrr.io/r/stats/predict.html) supports. Under a
non-identity link it is heteroscedastic by construction.

[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) returns the
frame the fit was solved on, or builds the frame of `data` when one is
supplied, through the terms and factor levels the fit recorded. A
data-dependent term such as `poly(x, 2)` is therefore evaluated at its
fitted basis and a factor gets its fitted levels, as they are under
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md).
`data` has to carry the response, since a model frame holds every
variable the formula names; covariate values that carry no response are
what [`predict()`](https://rdrr.io/r/stats/predict.html) and
[deli-augment](https://r-causal.github.io/deli/reference/deli-augment.md)
take a `newdata` for.

[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) returns
the design the fit was solved on, or the design of `data` when one is
supplied, coded with the contrasts and factor levels the fit recorded
rather than with whatever `getOption("contrasts")` says when the call is
made. A design matrix is a property of the fit, so a fit made under one
setting of that option answers with the coding it was solved on under
any other, which is what
[`predict()`](https://r-causal.github.io/deli/reference/deli-predict.md)
does as well. A design names no response, so `data` needs only the
predictors.

[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) is the
number of observations less the number of parameters. Both counts belong
to the fit rather than to a specification, so it answers for a fit built
from a `stacked_equations` function as readily as for a formula one.
[`weights()`](https://rdrr.io/r/stats/weights.html) returns the
observation weights the fit was solved with, which reach a fit only
through the formula interface. A formula fit specified without weights
returns `NULL`, as
[`stats::weights()`](https://rdrr.io/r/stats/weights.html) does for an
unweighted [`stats::lm()`](https://rdrr.io/r/stats/lm.html) fit: a
vector of ones would be indistinguishable from a fit weighted by ones.

[`sigma()`](https://rdrr.io/r/stats/sigma.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), and
[`deviance()`](https://rdrr.io/r/stats/deviance.html) are errors for
every fit, and say why. An M-estimator is defined by its estimating
equations, which need come from no likelihood at all, so there is no
log-likelihood to return and no deviance to take from one;
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) reach a fit through
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and are refused with
it. A residual standard deviation belongs to the model an equation
states rather than to the solve: a linear equation has one, a logistic
equation has none, and nothing in the estimates says which was solved.
Reporting an error is the only answer that does not invent a quantity
the fit never had.

## See also

[deli-predict](https://r-causal.github.io/deli/reference/deli-predict.md)
for predictions at new covariate values, and
[deli-augment](https://r-causal.github.io/deli/reference/deli-augment.md)
for the fitted values, intervals, and residuals as columns beside the
data.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

coef(fit)
#> (Intercept)          wt          hp 
#> 37.22727012 -3.87783074 -0.03177295 

vcov(fit)
#>              (Intercept)           wt            hp
#> (Intercept)  3.759414267 -0.991168177 -1.918899e-03
#> wt          -0.991168177  0.384310278 -1.649189e-03
#> hp          -0.001918899 -0.001649189  4.417009e-05

confint(fit)
#>                   lower       upper
#> (Intercept) 33.42705498 41.02748525
#> wt          -5.09286659 -2.66279490
#> hp          -0.04479898 -0.01874691

nobs(fit)
#> [1] 32

df.residual(fit)
#> [1] 29

head(fitted(fit))
#>         Mazda RX4     Mazda RX4 Wag        Datsun 710    Hornet 4 Drive 
#>          23.57233          22.58348          25.27582          21.26502 
#> Hornet Sportabout           Valiant 
#>          18.32727          20.47382 

head(residuals(fit))
#>         Mazda RX4     Mazda RX4 Wag        Datsun 710    Hornet 4 Drive 
#>        -2.5723294        -1.5834826        -2.4758187         0.1349799 
#> Hornet Sportabout           Valiant 
#>         0.3727334        -2.3738163 

formula(fit)
#> mpg ~ wt + hp
#> <environment: 0x55ed2f344b88>

# Weights reach a fit through the formula interface, and `weights()` reports
# the vector the fit was solved with.
weighted <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                       model = "linear", weights = 1 / hp)

head(weights(weighted))
#> [1] 0.009090909 0.009090909 0.010752688 0.009090909 0.005714286 0.009523810
```
