# Estimating equation for additive regression (GAM)

Generalized Additive Model via L2-penalized splines. Internally expands
`X` using
[`additive_design_matrix()`](https://r-causal.github.io/deli/reference/additive_design_matrix.md)
and delegates to
[`ee_bridge_regression()`](https://r-causal.github.io/deli/reference/ee_bridge_regression.md)
with `gamma = 2` (ridge penalty). The penalty only applies to the spline
basis terms, not to the original linear terms.

## Usage

``` r
ee_additive_regression(
  theta,
  X,
  y,
  specifications,
  model,
  weights = NULL,
  offset = NULL
)
```

## Arguments

- theta:

  Numeric vector of length equal to the number of columns in the
  expanded additive design matrix.

- X:

  Numeric n-by-b design matrix (before spline expansion).

- y:

  Numeric vector of n observed outcome values.

- specifications:

  A list of length `b` controlling spline generation. Each element is
  either `NULL` (no spline) or a list with keys `knots`, and optionally
  `natural`, `power`, `penalty`, `normalized`. See
  [`additive_design_matrix()`](https://r-causal.github.io/deli/reference/additive_design_matrix.md)
  for details.

- model:

  Character string: `"linear"`, `"logistic"`, or `"poisson"`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix, where p is the number of columns in the expanded
additive design matrix.

## Examples

``` r
set.seed(42)
n <- 200
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = 0.3)
X <- cbind(1, x)

# No spline on the intercept column, a penalized spline on x.
specs <- list(NULL, list(knots = c(-2, -1, 0, 1, 2), penalty = 5))

psi <- function(theta) {
  ee_additive_regression(
    theta,
    X = X,
    y = y,
    specifications = specs,
    model = "linear"
  )
}

# One parameter per column of the expanded design matrix.
m <- m_estimate(
  stacked_equations = psi,
  init = rep(0, ncol(additive_design_matrix(X, specs)))
)
coef(m)
#>     theta_1     theta_2     theta_3     theta_4     theta_5     theta_6 
#> -1.34755291 -0.28958141  0.20753205 -0.39959325  0.04045066  0.07544165 
```
