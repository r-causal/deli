# Delta method for variance of transformed parameters

Computes the variance-covariance matrix for a transformation of
parameters using the Delta Method: \$\$Var\[g(\theta)\] \approx G \Sigma
G^T\$\$ where \\G\\ is the Jacobian of \\g\\ and \\\Sigma\\ is the
covariance matrix of \\\theta\\.

## Usage

``` r
delta_method(
  object,
  transform,
  covariance = NULL,
  deriv_method = "capprox",
  dx = 1e-09,
  ...
)
```

## Arguments

- object:

  A fitted `MEstimator` object, or a numeric vector of parameter
  estimates.

- transform:

  Function that takes `theta` and returns a numeric vector.

- covariance:

  Numeric covariance matrix (only used when `object` is a numeric
  vector).

- deriv_method:

  Character string for the derivative method used to build the Jacobian
  of `transform`. One of `"capprox"` (central difference), `"fapprox"`
  (forward difference), `"bapprox"` (backward difference), or `"exact"`
  (forward-mode automatic differentiation). Default `"capprox"`.

- dx:

  Numeric step size for the finite-difference methods; ignored when
  `deriv_method = "exact"`. Default `1e-9`. Must be a single positive
  finite number, which is checked whichever `deriv_method` is in force.
  The step is absolute and is floored at the floating-point resolution
  of each estimate, so a large parameter magnitude cannot silently
  reduce it to nothing; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored. Exact names matter
  here because `deriv_method` selects how the Jacobian is built, and a
  dropped misspelling would leave the default in place and return a
  different variance with nothing to signal the substitution.

## Value

A covariance matrix for `transform(theta)`.

## Examples

``` r
fit <- m_estimate(vs ~ mpg, data = mtcars, .ee = ee_regression,
                  model = "logistic")

# Variance of the odds ratio for mpg, exponentiating the log-odds coefficient
delta_method(fit, transform = function(theta) exp(theta[2]))
#>            [,1]
#> [1,] 0.06420381

# The same variance from the estimates and covariance alone
delta_method(coef(fit), transform = function(theta) exp(theta[2]),
             covariance = vcov(fit))
#>            [,1]
#> [1,] 0.06420381
```
