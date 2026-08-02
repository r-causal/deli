# Regression Models

``` r

library(deli)
```

## Overview

deli provides estimating equations for a wide range of regression
models. Every model in this vignette is fitted the same way: pass a
formula, a data frame, and the estimating equation to
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md),
which builds the design matrix, constructs the estimator, and solves it
in one call. Arguments belonging to the estimating equation, such as
`model` or `penalty`, go in the same call and are forwarded to it, while
`init`, `subset`, `solver`, and the other solver controls are arguments
of
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
itself. The sandwich variance estimator automatically provides robust
standard errors.

The parameters are labeled from the first of three sources that names
every one of them: the names on `init`, the columns of the design
matrix, and the row names the estimating equation writes on its own
return. A default `init` carries no names, so the design columns label
the fit, together with any parameter the equation estimates on top of
the coefficients; an explicit `init` is labeled the same way when its
length accounts for one of those two shapes. Where it accounts for
neither, the equation’s row names label the fit instead, and where
nothing labels every parameter they are numbered `theta_1` through
`theta_p`. Name the elements of `init` yourself to label the parameters
any other way.

`init` defaults to a zero vector with one element per design matrix
column. That is the right length for most estimating equations, but not
for the few that estimate a parameter of their own on top of the
coefficients; the gamma and negative binomial GLMs below are the cases
you are likely to meet.

Every model shown here can also be fitted through the function
interface, where you write a `psi` function that supplies the design
matrix and response yourself. That form is what you need for a design
the formula notation cannot express, for a response that is not a vector
([`ee_mlogit()`](https://r-causal.github.io/deli/reference/ee_mlogit.md)),
and for custom or stacked equations; see
[`vignette("custom-estimating-equations")`](https://r-causal.github.io/deli/articles/custom-estimating-equations.md).

## Linear regression

The most basic regression model uses
[`ee_regression()`](https://r-causal.github.io/deli/reference/ee_regression.md)
with `model = "linear"`:

``` r

set.seed(42)
n <- 300
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(n)
d <- data.frame(x1, x2, y)

m <- m_estimate(y ~ x1 + x2, data = d, .ee = ee_regression, model = "linear")
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 300
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)     0.8984     0.0838    10.7262     0.7342     1.0625     <2e-16    86.7533
#> x1              2.0730     0.0609    34.0581     1.9537     2.1923     <2e-16   842.1482
#> x2             -0.4307     0.1146    -3.7576    -0.6554    -0.2061   0.000172    12.5091
```

## Logistic regression

For binary outcomes, use `model = "logistic"`:

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
y <- rbinom(n, 1, plogis(0.5 + x))
d <- data.frame(x, y)

m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "logistic")
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 500
#> Parameters: 2
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)     0.5943     0.1030     5.7688     0.3924     0.7962   7.99e-09    26.9000
#> x               1.0916     0.1270     8.5919     0.8426     1.3406     <2e-16    56.6986
```

## Poisson regression

For count data, use `model = "poisson"`:

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
y <- rpois(n, lambda = exp(0.5 + 0.3 * x))
d <- data.frame(x, y)

m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "poisson")
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 500
#> Parameters: 2
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)     0.4729     0.0357    13.2496     0.4030     0.5429     <2e-16   130.6947
#> x               0.2918     0.0358     8.1392     0.2215     0.3620   3.98e-16    51.1581
```

## GLM: generalized linear models

[`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md)
provides a more flexible interface where you specify the distribution
and link function separately:

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
y <- rpois(n, lambda = exp(0.5 + 0.3 * x))
d <- data.frame(x, y)

m <- m_estimate(
  y ~ x,
  data = d,
  .ee = ee_glm,
  distribution = "poisson",
  link = "log"
)
m@theta
#> (Intercept)           x 
#>   0.4729439   0.2917693
```

Available distributions: `"normal"`, `"binomial"`, `"poisson"`,
`"gamma"`, `"negative_binomial"`, `"inverse_gaussian"`, and `"tweedie"`.
`"tweedie"` takes a variance-power `hyperparameter`; see
[`?ee_glm`](https://r-causal.github.io/deli/reference/ee_glm.md) for the
details. Available links: `"identity"`, `"log"`, `"logit"`, `"probit"`,
`"cauchy"` (alias `"cauchit"`), `"loglog"`, `"cloglog"`, `"inverse"`,
and `"sqrt"`.

`"gamma"` and `"negative_binomial"` estimate one parameter beyond the
regression coefficients: the log of the gamma shape, or the log of the
negative binomial dispersion. The automatic `init` has one element per
design matrix column and so is one element short for these two, which
makes an explicit `init` necessary. It needs no names of its own: the
formula interface labels an unnamed `init` of that length from the model
matrix columns and the extra parameter.

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
mu <- exp(0.5 + 0.3 * x)
y <- rgamma(n, shape = 2, scale = mu / 2)
d <- data.frame(x, y)

m <- m_estimate(
  y ~ x,
  data = d,
  .ee = ee_glm,
  distribution = "gamma",
  link = "log",
  init = c(0, 0, 0)
)
m@theta
#> (Intercept)           x   log_shape 
#>   0.4945950   0.2546939   0.7053180
```

The extra parameter is on the log scale, so
`exp(m@theta[["log_shape"]])` recovers the shape, here close to the
value of 2 used in the simulation.

## Penalized regression

deli supports several penalized regression methods. These add a penalty
term to the estimating equations.

### Ridge regression

L2 penalty shrinks coefficients toward zero:

``` r

set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 0.5 * x1 + 0.3 * x2 + rnorm(n)
d <- data.frame(x1, x2, y)

m <- m_estimate(
  y ~ x1 + x2,
  data = d,
  .ee = ee_ridge_regression,
  model = "linear",
  penalty = 0.5
)

summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 200
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)     0.9377     0.0729    12.8704     0.7949     1.0805     <2e-16   123.5096
#> x1              0.4468     0.0852     5.2462     0.2799     0.6137   1.55e-07    22.6184
#> x2              0.3661     0.0802     4.5648     0.2089     0.5232      5e-06    17.6096
```

### LASSO regression

L1 penalty produces sparse solutions:

``` r

# Not differentiable, so the sandwich variance should not be trusted here
m <- m_estimate(
  y ~ x1 + x2,
  data = d,
  .ee = ee_lasso_regression,
  model = "linear",
  penalty = 0.1
)
#> Warning: The estimating equation for the chosen penalized regression model is not always
#> differentiable. Therefore, the bread matrix is not always defined for finite
#> samples, and the sandwich should not be used to estimate the variance.
m@theta
#> (Intercept)          x1          x2 
#>   0.9395739   0.4475518   0.3665738
```

### Elastic net

Combines L1 and L2 penalties:

``` r

# The L1 half is not differentiable, so again distrust the sandwich variance
m <- m_estimate(
  y ~ x1 + x2,
  data = d,
  .ee = ee_elasticnet_regression,
  model = "linear",
  penalty = 0.1,
  ratio = 0.5
)
#> Warning: The estimating equation for the chosen penalized regression model is not always
#> differentiable. Therefore, the bread matrix is not always defined for finite
#> samples, and the sandwich should not be used to estimate the variance.
m@theta
#> (Intercept)          x1          x2 
#>   0.9395912   0.4477125   0.3667634
```

## Robust regression

[`ee_robust_regression()`](https://r-causal.github.io/deli/reference/ee_robust_regression.md)
replaces the squared loss with a robust loss function, providing
resistance to outliers:

``` r

set.seed(42)
n <- 200
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n)
# Add some outliers
y[1:5] <- y[1:5] + 20
d <- data.frame(x, y)

# The Huber loss is convex, so its estimating function has a single root. What
# makes the seed necessary is that the Huber psi is bounded: far from the
# solution every residual is past the tuning constant k, every contribution
# saturates at k, and the estimating function is constant with a Jacobian of
# exactly zero. Starting from zero lands in that flat region, where the solver
# has no slope to follow, so start from a least-squares fit.
start <- coef(lm(y ~ x, data = d))

# Huber loss with k = 1.345
m <- m_estimate(
  y ~ x,
  data = d,
  .ee = ee_robust_regression,
  model = "linear",
  loss = "huber",
  k = 1.345,
  init = start
)

# Compare with OLS (affected by outliers)
m_ols <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")

rbind(robust = m@theta, ols = m_ols@theta)
#>        (Intercept)        x
#> robust    1.033468 1.949058
#> ols       1.515954 2.169921
```

Available loss functions: `"huber"`, `"tukey"`, `"andrew"`, `"hampel"`.

## Weighted regression

All regression EEs support observation weights via the `weights`
argument. Arguments passed through
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
are evaluated in the data frame, so the column name is enough:

``` r

set.seed(42)
n <- 200
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n)
w <- runif(n, 0.5, 1.5)
d <- data.frame(x, y, w)

m <- m_estimate(
  y ~ x,
  data = d,
  .ee = ee_regression,
  model = "linear",
  weights = w
)
m@theta
#> (Intercept)           x 
#>   0.9894503   1.9304272
```

## Predictions after regression

After fitting any regression model, use
[`generics::augment()`](https://generics.r-lib.org/reference/augment.html)
for predicted values with confidence intervals. Give it a data frame of
new covariate values as `newdata`, with one column for each covariate
the formula names, and it returns that frame with `.fitted`, `.se.fit`,
`.lower`, and `.upper` beside it:

``` r

set.seed(42)
n <- 300
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n)
d <- data.frame(x, y)

m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")

# Predict at new covariate values
generics::augment(m, newdata = data.frame(x = seq(-2, 2, length.out = 5)))
#>    x   .fitted    .se.fit     .lower     .upper
#> 1 -2 -3.077760 0.11004179 -3.2934376 -2.8620817
#> 2 -1 -1.052358 0.07100359 -1.1915228 -0.9131938
#> 3  0  0.973043 0.05711591  0.8610979  1.0849882
#> 4  1  2.998444 0.08228820  2.8371625  3.1597263
#> 5  2  5.023846 0.12477529  4.7792907  5.2684008
```

Called without `newdata`, `augment()` reports the rows the model was
fitted to and adds a `.resid` column as well. For a model with a
non-identity link, `type.predict = "response"` moves `.fitted` and its
interval from the linear predictor to the scale of the response.
