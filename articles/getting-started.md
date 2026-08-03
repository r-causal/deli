# Getting Started with deli

## Overview

deli is an R package for M-estimation and empirical sandwich variance
estimation. It provides:

1.  [`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md),
    a one-step fitting function with a formula interface for the
    pre-built estimating equations and a function interface for your own
2.  Pre-built estimating equations for common models (means, regression,
    causal inference, survival analysis, and more)
3.  Automatic sandwich variance estimation
4.  The `MEstimator` and `GMMEstimator` classes, which solve arbitrary
    user-specified estimating equations and hold the results

deli is an R port of the Python
[Delicatessen](https://github.com/pzivich/Delicatessen) library.

## Fitting your first model

The quickest way to fit a model is
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md),
which takes a formula, a data frame, and a pre-built estimating
equation. It constructs the estimator and solves it in a single call,
much as [`lm()`](https://rdrr.io/r/stats/lm.html) fits a linear model in
one step. Here is a linear regression on the `mtcars` data:

``` r

library(deli)

fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_regression,
  model = "linear"
)

coef(fit)
#> (Intercept)          wt          hp 
#> 37.22727012 -3.87783074 -0.03177295
```

The `.ee` argument names the estimating equation, and any further
arguments (such as `model = "linear"`) are passed on to it. The formula
and data supply the response and design matrix, so you do not need to
build them by hand.

## Inspecting the fit

Fitted estimators support the standard R accessors.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the point
estimates, [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the
sandwich variance-covariance matrix, and
[`confint()`](https://rdrr.io/r/stats/confint.html) returns Wald
confidence intervals:

``` r

vcov(fit)
#>              (Intercept)           wt            hp
#> (Intercept)  3.759412749 -0.991167246 -1.918909e-03
#> wt          -0.991167246  0.384309827 -1.649185e-03
#> hp          -0.001918909 -0.001649185  4.417008e-05

confint(fit)
#>                   lower       upper
#> (Intercept) 33.42705575 41.02748449
#> wt          -5.09286587 -2.66279561
#> hp          -0.04479898 -0.01874691
```

[`summary()`](https://rdrr.io/r/base/summary.html) collects the
estimates, standard errors, test statistics, and confidence intervals
into a single table:

``` r

summary(fit)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   270.5102
#> wt             -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    31.2310
#> hp             -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    19.1270
```

deli also provides [broom](https://broom.tidymodels.org/) tidiers, so
results flow into tidyverse pipelines.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
row per parameter, and
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
one-row model summary:

``` r

tidy(fit, conf.int = TRUE)
#>          term    estimate   std.error statistic      p.value   s.value
#> 1 (Intercept) 37.22727012 1.938920511 19.199998 3.701056e-82 270.51017
#> 2          wt -3.87783074 0.619927275 -6.255299 3.967539e-10  31.23104
#> 3          hp -0.03177295 0.006646057 -4.780721 1.746674e-06  19.12696
#>      conf.low   conf.high
#> 1 33.42705575 41.02748449
#> 2 -5.09286587 -2.66279561
#> 3 -0.04479898 -0.01874691

glance(fit)
#>   nobs npar  estimator finite_correction j_statistic j_df j_p_value
#> 1   32    3 MEstimator              <NA>          NA   NA        NA
```

If you prefer deli’s own inference helpers,
[`confidence_intervals()`](https://r-causal.github.io/deli/reference/confidence_intervals.md),
[`z_scores()`](https://r-causal.github.io/deli/reference/z_scores.md),
[`p_values()`](https://r-causal.github.io/deli/reference/p_values.md),
and
[`s_values()`](https://r-causal.github.io/deli/reference/s_values.md)
each return the corresponding quantities directly:

``` r

z_scores(fit)
#> (Intercept)          wt          hp 
#>   19.199998   -6.255299   -4.780721
p_values(fit)
#>  (Intercept)           wt           hp 
#> 3.701056e-82 3.967539e-10 1.746674e-06
```

## Other regression models

The pre-built regression equation covers several model families. Switch
from linear to logistic regression by changing the `model` argument.
Here `vs` is a binary indicator of engine type:

``` r

fit_logistic <- m_estimate(
  vs ~ mpg + wt,
  data = mtcars,
  .ee = ee_regression,
  model = "logistic"
)

summary(fit_logistic)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)   -12.5412     5.6235    -2.2301   -23.5631    -1.5193     0.0257     5.2799
#> mpg             0.5241     0.2041     2.5676     0.1240     0.9241     0.0102     6.6095
#> wt              0.5829     0.7276     0.8011    -0.8432     2.0089      0.423     1.2410
```

The sandwich variance estimator supplies robust standard errors
automatically, without any additional arguments.
[`ee_regression()`](https://r-causal.github.io/deli/reference/ee_regression.md)
covers linear, logistic, and Poisson models. See
[`?ee_glm`](https://r-causal.github.io/deli/reference/ee_glm.md) for the
wider family of generalized linear models, and the dedicated functions
[`?ee_robust_regression`](https://r-causal.github.io/deli/reference/ee_robust_regression.md),
[`?ee_ridge_regression`](https://r-causal.github.io/deli/reference/ee_ridge_regression.md),
[`?ee_lasso_regression`](https://r-causal.github.io/deli/reference/ee_lasso_regression.md),
[`?ee_elasticnet_regression`](https://r-causal.github.io/deli/reference/ee_elasticnet_regression.md),
[`?ee_beta_regression`](https://r-causal.github.io/deli/reference/ee_beta_regression.md),
and [`?ee_tobit`](https://r-causal.github.io/deli/reference/ee_tobit.md)
for robust, ridge, LASSO, elastic net, beta, and tobit regression.

## The estimator objects

[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
builds and solves an `MEstimator` for you. When a pre-built estimating
equation and the formula interface do not fit your problem, for example
when you supply a custom estimating equation, give
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
a function in place of the formula. The workflow has two steps:

1.  Define a `psi` function that returns the estimating equation
    contributions
2.  Pass it to
    [`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
    with starting values for the parameters

Here is a simple example that estimates a mean:

``` r

y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)

# Define the estimating equation
psi <- function(theta) {
  ee_mean(theta, y = y)
}

m <- m_estimate(psi, init = 0)

coef(m)
#> theta_1 
#>     3.4
```

The estimated mean matches `mean(y)`:

``` r

mean(y)
#> [1] 3.4
```

The result is an `MEstimator` object, the same class the formula
interface returns, so the same accessors work here. The point estimate
and variance are also available on the object as `m@theta` and
`m@variance`:

``` r

m@theta       # Point estimate (the mean)
#> theta_1 
#>     3.4
m@variance    # Sandwich variance estimate
#>           theta_1
#> theta_1 0.3839999
```

## Stacking estimating equations

A key strength of M-estimation is the ability to stack estimating
equations. When you stack equations, the sandwich variance correctly
accounts for all sources of uncertainty. Here is an example that
estimates a mean and variance simultaneously:

``` r

y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)

psi <- function(theta) {
  ee_mean_variance(theta, y = y)
}

m <- m_estimate(stacked_equations = psi, init = c(0, 0))
m@theta  # c(mean, variance)
#>     mean variance 
#>     3.40     3.84
```

For more complex stacking with custom equations, see
[`vignette("custom-estimating-equations")`](https://r-causal.github.io/deli/articles/custom-estimating-equations.md).

## Causal inference: inverse probability weighting

deli includes estimating equations for causal inference. Because they
combine several components, they are a natural fit for the function
interface. Here is an example that uses inverse probability weighting
(IPW) to estimate an average treatment effect:

``` r

set.seed(42)
n <- 500
w <- rbinom(n, 1, 0.5)               # Binary confounder
A <- rbinom(n, 1, plogis(-0.5 + w))  # Treatment depends on w
Y <- 1 + 2 * A + w + rnorm(n)        # Outcome
W <- cbind(1, w)                     # Propensity score design matrix

psi <- function(theta) {
  ee_ipw(theta, y = Y, A = A, W = W)
}

# theta: ACE, E[Y(1)], E[Y(0)], beta0, beta1
m <- m_estimate(stacked_equations = psi, init = c(0, 0, 0, 0, 0))

# ACE (average causal effect) ~ 2
m
#> <MEstimator>
#>   Parameters: 5
#>   Observations: 500
#> Coefficients:
#> ACE: 1.9909
#> E[Y^1]: 3.4609
#> E[Y^0]: 1.4701
#> W_1: -0.4688
#> W_2: 1.1492
```

## Predictions

After fitting a model, use
[`augment()`](https://generics.r-lib.org/reference/augment.html) for
predicted values with confidence intervals. Give it a data frame of new
covariate values as `newdata`, with one column for each covariate the
formula names, and it returns that frame with `.fitted`, `.se.fit`,
`.lower`, and `.upper` beside it:

``` r

set.seed(42)
n <- 200
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n)
d <- data.frame(x, y)

m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")

# Predict at new values
augment(m, newdata = data.frame(x = seq(-2, 2, by = 1)))
#>    x   .fitted    .se.fit     .lower     .upper
#> 1 -2 -2.834668 0.14108320 -3.1111856 -2.5581497
#> 2 -1 -0.912765 0.08435576 -1.0780992 -0.7474307
#> 3  0  1.009138 0.06706354  0.8776956  1.1405798
#> 4  1  2.931040 0.10976950  2.7158961  3.1461846
#> 5  2  4.852943 0.17254381  4.5147634  5.1911227
```

Called without `newdata`,
[`augment()`](https://generics.r-lib.org/reference/augment.html) reports
the rows the model was fitted to and adds a `.resid` column as well.

## Clustered data

For clustered or grouped data, use
[`aggregate_efuncs()`](https://r-causal.github.io/deli/reference/aggregate_efuncs.md)
inside your `psi` function to get cluster-robust variance estimates:

``` r

set.seed(42)
n <- 200
n_groups <- 50
group <- rep(1:n_groups, each = n / n_groups)
group_effect <- rnorm(n_groups, sd = 2)
y <- group_effect[group] + rnorm(n)

# Cluster-robust variance
psi <- function(theta) {
  ef <- ee_mean(theta, y = y)
  aggregate_efuncs(ef, group = group)
}

m <- m_estimate(stacked_equations = psi, init = mean(y))

m@theta
#>     theta_1 
#> -0.08792509
m@variance  # Accounts for within-cluster correlation
#>           theta_1
#> theta_1 0.1058511
```

## Exact differentiation

The delta method needs the derivative (Jacobian) of a transform, and the
sandwich variance needs the derivative of the estimating equations. By
default deli computes these with central finite differences
(`deriv_method = "capprox"`), which introduces a small step-size
approximation. Passing `deriv_method = "exact"` switches to forward-mode
automatic differentiation, which evaluates the derivative in closed
form. It returns exact derivatives with no step size to tune. In Python
Delicatessen exact differentiation is the default only for
`delta_method`; the estimators and the sandwich variance default to a
forward finite difference.

The
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md)
function accepts `deriv_method` directly, and the survival prediction
helpers
[`survival_predictions()`](https://r-causal.github.io/deli/reference/survival_predictions.md)
and
[`aft_predictions_function()`](https://r-causal.github.io/deli/reference/aft_predictions_function.md)
forward it to the delta method. The following fits a logistic model and
uses exact differentiation to obtain the variance of the fitted
probability at a covariate pattern:

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
pr <- plogis(0.5 + x)
y <- rbinom(n, 1, pr)
d <- data.frame(x, y)

m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "logistic")

# Variance of the predicted probability at x = 1, via the delta method.
# The log-odds at the pattern (intercept = 1, x = 1) is theta[1] + theta[2].
transform <- function(theta) inverse_logit(theta[1] + theta[2])
delta_method(m, transform = transform, deriv_method = "exact")
#>              [,1]
#> [1,] 0.0005451013
```

That chunk computes the inverse logit twice with two different
functions, and the difference is deliberate. The simulation uses
[`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html), because it
is ordinary numeric work on plain doubles. The transform uses deli’s
[`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md),
because [`plogis()`](https://rdrr.io/r/stats/Logistic.html) cannot be
differentiated exactly and would stop the delta method with an error.
The two functions return identical numbers; only one of them survives
`deriv_method = "exact"`.

## Next steps

- [`vignette("custom-estimating-equations")`](https://r-causal.github.io/deli/articles/custom-estimating-equations.md):
  Write your own estimating equations and stack them
- [`?m_estimate`](https://r-causal.github.io/deli/reference/m_estimate.md):
  the one-step formula and function interface
- [`?ee_regression`](https://r-causal.github.io/deli/reference/ee_regression.md):
  linear, logistic, and Poisson regression
- [`?ee_glm`](https://r-causal.github.io/deli/reference/ee_glm.md):
  generalized linear models, plus
  [`?ee_ridge_regression`](https://r-causal.github.io/deli/reference/ee_ridge_regression.md),
  [`?ee_lasso_regression`](https://r-causal.github.io/deli/reference/ee_lasso_regression.md),
  and the other penalized and robust regression functions
- [`?ee_gformula`](https://r-causal.github.io/deli/reference/ee_gformula.md):
  G-computation / standardization
- [`?ee_aipw`](https://r-causal.github.io/deli/reference/ee_aipw.md):
  Augmented inverse probability weighting
- [`?ee_aft`](https://r-causal.github.io/deli/reference/ee_aft.md):
  Accelerated failure time models
- [`?ee_survival_model`](https://r-causal.github.io/deli/reference/ee_survival_model.md):
  Parametric survival models
