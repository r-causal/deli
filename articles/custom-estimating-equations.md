# Custom Estimating Equations

``` r

library(deli)
```

## Overview

The deli package lets you define custom estimating equations, stack them
with built-in ones, and get valid sandwich variance estimates for the
full system. This vignette covers three topics:

1.  Writing a custom estimating equation from scratch
2.  Stacking custom and built-in estimating equations
3.  Using the delta method for post-hoc transformations

## Writing a custom estimating equation

The core requirement is simple: your `psi` function must take a
parameter vector `theta` and return a **p-by-n matrix**, where p is the
number of parameters and n is the number of observations. Each row
corresponds to one parameter’s estimating equation, and each column
corresponds to one observation’s contribution.

### Example: estimating a ratio of means

Suppose we observe two variables `y1` and `y2` and want to estimate the
ratio of their means, `mu1 / mu2`. We need three parameters:

- `theta[1]`: the mean of `y1` (i.e., `mu1`)
- `theta[2]`: the mean of `y2` (i.e., `mu2`)
- `theta[3]`: the ratio `mu1 / mu2`

The estimating equations are:

1.  `y1_i - theta[1]` (solves for the mean of `y1`)
2.  `y2_i - theta[2]` (solves for the mean of `y2`)
3.  `theta[1] / theta[2] - theta[3]` (solves for the ratio)

The third equation has no variation across observations, so we repeat
the same value n times.

``` r

set.seed(42)
n <- 200
y1 <- rnorm(n, mean = 4, sd = 1)
y2 <- rnorm(n, mean = 2, sd = 1)

psi <- function(theta) {
  # Row 1: estimating equation for the mean of y1
  mu1 <- y1 - theta[1]
  # Row 2: estimating equation for the mean of y2
  mu2 <- y2 - theta[2]
  # Row 3: estimating equation for the ratio (repeated n times)
  ratio <- rep(theta[1] / theta[2] - theta[3], n)
  # Stack into a 3-by-n matrix. `rbind()` labels each row with the name of the
  # variable it came from, and those labels become the parameter names.
  rbind(mu1, mu2, ratio)
}

m <- m_estimate(stacked_equations = psi, init = c(1, 1, 1))
m@theta
#>      mu1      mu2    ratio 
#> 3.972516 2.011284 1.975114
```

The third element of `theta` is the estimated ratio. Because the ratio
is estimated jointly with the means, the sandwich variance accounts for
the uncertainty in all three parameters:

``` r

summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 200
#> Parameters: 3
#> 
#>         Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> mu1       3.9725     0.0687    57.7895     3.8378     4.1072     <2e-16        Inf
#> mu2       2.0113     0.0668    30.1086     1.8804     2.1422     <2e-16   659.1592
#> ratio     1.9751     0.0764    25.8637     1.8254     2.1248     <2e-16   487.5522
```

Compare the point estimate to the naive ratio:

``` r

mean(y1) / mean(y2)
#> [1] 1.975114
```

### Key points for custom equations

- **Always return a matrix.** Even for a single parameter, return a
  1-by-n matrix (e.g., `matrix(..., nrow = 1)`).
- **Rows = parameters, columns = observations.** The number of rows must
  match the length of `theta`, and the number of columns must equal n.
- **Deterministic equations** (like the ratio above) should repeat the
  same value across all n columns. This is necessary so the matrix
  dimensions are consistent with the other rows.
- **Row names name the parameters.** The row names of the returned
  matrix become the labels on
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html). Names on `init`
  take precedence when it has any. Row names are read only when every
  parameter is labeled and no two labels are alike, so a stack that
  names some rows and not others, or that repeats a name, is numbered
  `theta_1` through `theta_p` instead.
  [`rbind()`](https://rdrr.io/r/base/cbind.html) supplies a label for
  each plain vector it is given, taken from the variable name, as in the
  example above; a matrix argument contributes whatever row names it
  already carries. Many built-in `ee_*()` functions name their rows and
  many do not. The ones that do say so under **Value** on their help
  pages, along with what the names are, so a stack that mixes in one of
  the others is numbered unless you name it yourself. Most of the
  regression estimating equations leave their rows unnamed, because
  there is nothing to name: a regression coefficient has no name apart
  from the design column it multiplies, and deli drops a design’s column
  headings, so the labels left would say no more than `theta_1` already
  does. Four returns are different.
  [`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md)
  under `distribution = "gamma"` and
  `distribution = "negative_binomial"`, along with
  [`ee_tobit()`](https://r-causal.github.io/deli/reference/ee_tobit.md)
  and
  [`ee_beta_regression()`](https://r-causal.github.io/deli/reference/ee_beta_regression.md),
  each hold one row more than the design has columns, and that last row
  is a parameter of the outcome distribution: a log shape, a log
  dispersion, a log scale, or a log precision. It is the row you are
  most likely to read as a coefficient, so each of the four names it
  `log_shape`, `log_dispersion`, `log_sigma`, or `log_phi`, with the
  design rows labeled `X_1` through `X_p` beside it. Set the names
  yourself with [`rownames()`](https://rdrr.io/r/base/colnames.html)
  where the ones you want are not the ones you get. The assignment works
  under `deriv_method = "exact"` as well, unlike the reshaping helpers
  described in
  [`vignette("getting-started")`](https://r-causal.github.io/deli/articles/getting-started.md).

## Stacking custom EEs with built-in ones

A major strength of M-estimation is stacking: you can combine built-in
estimating equations with custom ones using
[`rbind()`](https://rdrr.io/r/base/cbind.html). The sandwich variance
estimator then correctly propagates uncertainty through the entire
system.

### Example: log odds-ratio from logistic regression

Suppose we fit a logistic regression and want to estimate the odds ratio
for a coefficient, along with a proper confidence interval. We can stack
the regression estimating equations with a custom equation that
exponentiates the log-odds coefficient.

``` r

set.seed(42)
n <- 500
x <- rnorm(n)
pr <- plogis(-0.5 + 0.8 * x)
y <- rbinom(n, 1, pr)
X <- cbind(1, x)

psi <- function(theta) {
  # theta[1:2]: logistic regression coefficients (intercept, slope)
  # theta[3]: odds ratio = exp(theta[2])
  beta <- theta[1:2]
  or <- theta[3]

  # Built-in logistic regression EE (returns a 2-by-n matrix)
  ee_reg <- ee_regression(beta, X = X, y = y, model = "logistic")

  # Custom EE for the odds ratio (deterministic, repeated n times)
  ee_or <- matrix(rep(exp(theta[2]) - or, n), nrow = 1)

  # Stack: 3-by-n matrix
  rbind(ee_reg, ee_or)
}

m <- m_estimate(stacked_equations = psi, init = c(0, 0, 1))
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 500
#> Parameters: 3
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1    -0.4364     0.0981    -4.4492    -0.6286    -0.2442   8.62e-06    16.8241
#> theta_2     0.8223     0.1097     7.4930     0.6072     1.0374   6.73e-14    43.7558
#> theta_3     2.2758     0.2498     9.1117     1.7863     2.7654     <2e-16    63.4195
```

The third row of the summary gives the odds ratio with a sandwich-based
confidence interval that correctly accounts for the estimation
uncertainty in the regression coefficients. Compare with a manual
calculation:

``` r

exp(m@theta[2])
#>  theta_2 
#> 2.275827
```

### How stacking works

When you [`rbind()`](https://rdrr.io/r/base/cbind.html) a 2-by-n matrix
from
[`ee_regression()`](https://r-causal.github.io/deli/reference/ee_regression.md)
with a 1-by-n matrix from a custom equation, you get a 3-by-n matrix.
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
solves the full 3-parameter system simultaneously and computes the
sandwich variance for all parameters at once. This is what makes the
variance estimates valid: the covariance between the regression
coefficients and the odds ratio is captured automatically.

## Using the delta method

The delta method is an alternative to stacking for obtaining variance
estimates of transformed parameters. Instead of adding extra equations
to the system, you apply a transformation after estimation and use a
first-order approximation to compute the variance.

### Example: odds ratio via the delta method

Using the same logistic regression as above, we can get the odds ratio
variance without stacking. Nothing custom is stacked onto the regression
here, so the formula interface applies: give
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
a formula and a data frame and it builds the design matrix and response
for you.

``` r

# Fit the logistic regression only
d <- data.frame(x, y)

m_reg <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "logistic")
m_reg@theta
#> (Intercept)           x 
#>  -0.4363908   0.8223434
```

Now apply the delta method. The `transform` function takes the full
`theta` vector and returns the transformed quantity:

``` r

# Transform: exponentiate the second coefficient to get the odds ratio
dm_var <- delta_method(m_reg, transform = function(theta) exp(theta[2]))
dm_var
#>            [,1]
#> [1,] 0.06238429
```

The result is the variance of the odds ratio. We can compute a
confidence interval:

``` r

or_est <- exp(m_reg@theta[2])
or_se <- sqrt(dm_var[1, 1])
ci <- or_est + c(-1, 1) * qnorm(0.975) * or_se
round(c(or = unname(or_est), lower = ci[1], upper = ci[2]), 3)
#>    or lower upper 
#> 2.276 1.786 2.765
```

### When to use the delta method vs. stacking

- **Stacking** is preferred when the transformation involves additional
  data or when you want the transformed parameter to be part of the
  root-finding system. Stacking also makes it easy to add further
  equations that depend on the transformed parameter.
- **The delta method** is convenient for simple post-hoc transformations
  of already-estimated parameters. It avoids modifying the estimating
  equation system and does not require specifying initial values for the
  extra parameters.

Both approaches give the same asymptotic variance. Use whichever is more
natural for your problem.

### Using the delta method with raw estimates

You can also call
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md)
directly on a numeric vector of estimates and a covariance matrix,
without a fitted estimator object:

``` r

delta_method(
  m_reg@theta,
  transform = function(theta) exp(theta[2]),
  covariance = m_reg@variance
)
#>            [,1]
#> [1,] 0.06238429
```

This is useful when you have estimates and covariances from another
source (e.g., from [`glm()`](https://rdrr.io/r/stats/glm.html) or
another package) and want to apply the delta method.
