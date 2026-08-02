# Estimating equation for beta regression

Beta regression for outcomes in (0, 1) using mean-precision
parameterization. The last element of theta is log(phi), the log
precision parameter.

## Usage

``` r
ee_beta_regression(theta, X, y, weights = NULL, offset = NULL)
```

## Arguments

- theta:

  Numeric vector of length `b + 1`.

- X:

  Numeric n-by-b design matrix.

- y:

  Numeric vector of n outcomes in (0, 1).

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A `(b+1)`-by-n matrix. The rows are named `X_1` through `X_b` for the
columns of `X`, and the final row is named `log_phi`.

## Examples

``` r
set.seed(42)
n <- 50
W <- rnorm(n)
mu <- 1 / (1 + exp(-(0.5 + 0.3 * W)))
y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
X <- cbind(1, W)

psi <- function(theta) ee_beta_regression(theta, X = X, y = y)

# The last parameter is log(phi), started here at a precision of 10. The
# default solver diverges on this equation, so use nleqslv.
m <- m_estimate(
  stacked_equations = psi,
  init = c(0, 0, log(10)),
  solver = "nleqslv"
)
coef(m)
#>       X_1       X_2   log_phi 
#> 0.4620901 0.2783017 2.5566274 
```
