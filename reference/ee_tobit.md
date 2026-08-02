# Estimating equation for Tobit regression (Type I)

Handles left and/or right censored outcomes using standard normal
PDF/CDF. Theta is `(beta, log(sigma))`, a vector of length `b + 1`.

## Usage

``` r
ee_tobit(
  theta,
  X,
  y,
  lower = NULL,
  upper = NULL,
  weights = NULL,
  offset = NULL
)
```

## Arguments

- theta:

  Numeric vector of length `b + 1`.

- X:

  Numeric n-by-b design matrix.

- y:

  Numeric vector of n observed (possibly censored) outcome values.

- lower:

  Numeric lower censoring limit, or `NULL` (default, no left censoring).

- upper:

  Numeric upper censoring limit, or `NULL` (default, no right
  censoring).

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A `(b+1)`-by-n matrix. The rows are named `X_1` through `X_b` for the
columns of `X`, and the final row is named `log_sigma`.

## Examples

``` r
# A latent outcome observed only down to zero, so the negative values are
# left censored at the limit.
set.seed(123)
n <- 200
X <- cbind(1, rnorm(n))
y <- pmax(1 + 0.5 * X[, 2] + rnorm(n), 0)

psi <- function(theta) ee_tobit(theta, X = X, y = y, lower = 0)

# The last parameter is log(sigma), started at the log of the observed
# standard deviation.
m <- m_estimate(
  stacked_equations = psi,
  init = c(mean(y), 0, log(sd(y)))
)
coef(m)
#>          X_1          X_2    log_sigma 
#>  1.040942234  0.477997107 -0.002076965 
```
