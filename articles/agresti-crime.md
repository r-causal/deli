# Agresti & Finlay (2009): Robust Regression

> **Note**
>
> This article is translated from the [Agresti & Finlay (2009): Robust
> Regression
> example](https://deli.readthedocs.io/en/latest/Examples/Agresti-Crime.html)
> in the documentation of [delicatessen](https://deli.readthedocs.io/),
> deli’s Python counterpart.

``` r

library(deli)
```

This example demonstrates robust regression using data from *Statistical
Methods for the Social Sciences*, 4th Edition by Alan Agresti and
Barbara Finlay. Violent crime is modeled as a function of poverty and
the percentage of single parents, at the US state level. The applied
example follows the [UCLA
tutorial](https://stats.oarc.ucla.edu/r/dae/robust-regression/) with
data from the 2005 *Statistical Abstract of the United States*,
available [here](https://users.stat.ufl.edu/~aa/social/data.html).

## Data

The data contain 51 observations (50 US states plus Washington DC). We
load the data and standardize the crime outcome.

``` r

d <- crime
d$intercept <- 1
d$crime_std <- (d$crime - mean(d$crime)) / sd(d$crime)

X <- as.matrix(d[, c("intercept", "poverty", "single")])
y <- d$crime_std
```

## Linear Regression

We start by modeling the standardized crime rate as a function of
poverty and single parent percentage using ordinary least squares. The
estimating equation for linear regression is

\sum\_{i=1}^{n} \begin{bmatrix} (Y_i - \hat{Y}\_i) \\ (Y_i - \hat{Y}\_i)
P_i \\ (Y_i - \hat{Y}\_i) S_i \end{bmatrix} = 0

where Y is the crime rate, P is poverty, S the proportion of single
parents, and \hat{Y}\_i = \theta_0 + \theta_1 P_i + \theta_2 S_i.

``` r

estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  },
  init = c(0, 0, 0)
)

# Save OLS estimates for use as initial values in robust regression
ols_theta <- unname(coef(estr))

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.49 -5.72 -3.26
#> theta_2   Poverty  0.02 -0.03  0.06
#> theta_3    Single  0.38  0.26  0.49
```

## Robust Regression

Robust regression limits the influence of outliers by applying a
function f to the residuals. The resulting estimating equation is

\sum\_{i=1}^{n} f_k(Y_i - X_i^T \theta)\\ X_i = 0

where f_k is the chosen robust loss function with tuning parameter k.
Different choices of f_k impose different constraints on how much any
single observation can shift the solution.

Huber and Fair are monotone: a large residual contributes at most a
fixed amount, but it never stops contributing. The other six losses are
redescending, driving the contribution of a large residual back toward
zero. Andrew’s Sine, Tukey’s biweight, and Hampel reach zero at a hard
cut-off and hold it there, while Cauchy, Ullah, and Welsch approach zero
smoothly and never quite arrive.

A redescending estimating equation is non-convex and can have several
roots, so the starting values decide which root a solver reaches.
Cauchy, Andrew’s Sine, and Hampel do not solve from a zero start on this
data, and Ullah solves from one to a root no other fit here is near, so
those four fits begin at the least-squares estimates above instead.
Andrew’s Sine, Ullah, and Hampel use the Levenberg-Marquardt solver,
matching the reference analysis, and take the bread derivative by exact
automatic differentiation rather than by finite differencing.

### Huber

The Huber function caps residual contributions at \pm k. While k = 1.345
is commonly recommended, we use k = 1 here.

``` r

estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "huber", k = 1)
  },
  init = c(0, 0, 0)
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.67 -5.71 -3.62
#> theta_2   Poverty  0.02 -0.02  0.06
#> theta_3    Single  0.39  0.29  0.49
```

### Fair

Like Huber, the Fair function is monotonic but tapers continuously
rather than having a sharp hinge point. We set k = 1.

``` r

estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "fair", k = 1)
  },
  init = c(0, 0, 0)
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.56 -6.08 -3.04
#> theta_2   Poverty  0.02 -0.02  0.06
#> theta_3    Single  0.38  0.24  0.52
```

### Cauchy

The Cauchy function begins to reduce outlier influence beyond a certain
point but does not shrink contributions fully to zero. We set k = 1.

``` r

# Use OLS estimates as starting values. Redescending loss functions have
# multiple roots, and from c(0,0,0) this one reaches none of them: the solver
# runs the intercept out past five million and warns that the estimating
# functions are not solved where it stopped.
estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "cauchy", k = 1)
  },
  init = ols_theta
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.67 -6.54 -2.81
#> theta_2   Poverty  0.02 -0.02  0.06
#> theta_3    Single  0.39  0.21  0.56
```

### Andrew’s Sine

Andrew’s Sine allows the same maximum influence as Huber at k but varies
smoothly and shrinks outlier contributions to zero at high levels. We
set k = 1.

``` r

# The redescending losses are non-convex, and a zero start does not solve this
# one, so the fit begins at the OLS estimates. The Levenberg-Marquardt solver
# ("lm") follows the reference analysis. Andrew's score is zero beyond a hard
# cut-off, so its bread derivative is discontinuous there; we take the
# derivative by exact automatic differentiation to get a deterministic variance
# at that discontinuity.
estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "andrew", k = 1)
  },
  init = ols_theta,
  solver = "lm",
  deriv_method = "exact"
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.63 -5.92 -3.35
#> theta_2   Poverty  0.02 -0.02  0.06
#> theta_3    Single  0.38  0.27  0.50
```

### Welsch

The Welsch function has a similar shape to Andrew’s Sine but differs in
steepness at the same value of k. We set k = 1.

``` r

estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "welsch", k = 1)
  },
  init = c(0, 0, 0)
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.78 -6.11 -3.46
#> theta_2   Poverty  0.02 -0.02  0.07
#> theta_3    Single  0.39  0.27  0.52
```

### Ullah

The Ullah function shares a similar shape to the previous redescending
functions. We set k = 1.

``` r

# A zero start does solve this one, but to a root that sits well away from every
# other fit here, so the fit begins at the OLS estimates for the same reason the
# three above do. The Levenberg-Marquardt solver matches the reference analysis,
# and the derivative is taken by exact automatic differentiation.
estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "ullah", k = 1)
  },
  init = ols_theta,
  solver = "lm",
  deriv_method = "exact"
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.49 -6.81 -2.17
#> theta_2   Poverty  0.02 -0.02  0.07
#> theta_3    Single  0.37  0.15  0.58
```

### Tukey’s Biweight

Tukey’s biweight substantially reduces the influence of outlying
observations. Small values of k produce large standard errors, so we use
k = 5.

``` r

estr <- m_estimate(
  stacked_equations = function(theta) {
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "tukey", k = 5)
  },
  init = c(0, 0, 0)
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.55 -5.82 -3.28
#> theta_2   Poverty  0.02 -0.03  0.06
#> theta_3    Single  0.38  0.26  0.50
```

### Hampel

The Hampel function mimics Huber until a set point, then linearly
decreases contributions to zero. Unlike the other functions, it takes
three tuning parameters. In R these are passed as `k = c(a, b, c)` where
`a` is the first hinge, `b` the plateau upper limit, and `c` the zero
point.

``` r

# A zero start does not solve this one either, so the fit begins at the OLS
# estimates, with the Levenberg-Marquardt solver matching the reference analysis
# and the derivative taken by exact automatic differentiation. The tuning
# constants k = c(1/3, 4/3, 2) correspond to the reference parameters
# k = 2, lower = 1/3, upper = 4/3.
estr <- m_estimate(
  stacked_equations = function(theta) {
    # k = c(a, b, c) where a < b < c
    ee_robust_regression(theta, X = X, y = y, model = "linear",
                         loss = "hampel", k = c(1/3, 4/3, 2))
  },
  init = ols_theta,
  solver = "lm",
  deriv_method = "exact"
)

data.frame(
  Param = c("Intercept", "Poverty", "Single"),
  Coef = round(coef(estr), 2),
  LCL = round(confint(estr)[, 1], 2),
  UCL = round(confint(estr)[, 2], 2)
)
#>             Param  Coef   LCL   UCL
#> theta_1 Intercept -4.76 -7.58 -1.94
#> theta_2   Poverty  0.03 -0.02  0.07
#> theta_3    Single  0.39  0.12  0.66
```

## References

Agresti A & Finlay B. (2009). *Statistical Methods for the Social
Sciences*, 4th Edition. Pearson.

UCLA Advanced Research Computing. Robust regression.
<https://stats.oarc.ucla.edu/r/dae/robust-regression/>
