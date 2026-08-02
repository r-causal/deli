# Estimating equations for IPW marginal structural model

Estimates the parameters of a marginal structural model using inverse
probability weighting with a logistic propensity score model.

## Usage

``` r
ee_ipw_msm(
  theta,
  y,
  A,
  W,
  V,
  distribution,
  link,
  hyperparameter = NULL,
  truncate = NULL,
  weights = NULL
)
```

## Arguments

- theta:

  Numeric vector of length `c + b`, where `c` is the number of MSM
  parameters and `b` is the number of propensity score model parameters.

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n binary treatment indicators (0/1).

- W:

  Numeric n-by-b design matrix for the propensity score model.

- V:

  Numeric n-by-c design matrix for the marginal structural model.

- distribution:

  Character string for the GLM distribution of the outcome model. Every
  name [`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md)
  takes except the three that estimate a nuisance parameter, which the
  `theta` partition here reserves no slot for: `"normal"` (or
  `"gaussian"`), `"binomial"` (or `"bernoulli"`, or `"bin"`),
  `"poisson"`, `"inverse_normal"` (or `"inverse_gaussian"`), and
  `"tweedie"`. See `hyperparameter` for what naming one of the three
  does.

- link:

  Character string for the GLM link function. See
  [`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md) for
  options.

- hyperparameter:

  Optional numeric hyperparameter passed straight through to the
  marginal structural model's
  [`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md)
  call. Used by the tweedie outcome distribution, where it fixes the
  variance power `v(mu) = mu^hyperparameter`. Default `NULL`. Note that
  the theta partition reserves exactly `ncol(V)` slots for the MSM and
  no slot for an estimated nuisance parameter, so outcome families that
  carry one (`"gamma"`, `"negative_binomial"`, `"nb"`) cannot be used as
  the MSM outcome model. Naming one is refused by name, before the
  outcome model is formed and whatever `hyperparameter` is set to.
  Python Delicatessen cannot fit those families either, where the
  attempt fails as a shape error.

- truncate:

  Optional length-2 numeric vector `c(lower, upper)` to clip estimated
  propensity scores. Bounds must be in ascending order
  (`lower <= upper`). Default `NULL`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A `(c+b)`-by-n matrix of estimating equation contributions. The marginal
structural model rows are named `MSM alpha_0` through `MSM alpha_(c-1)`,
matching the zero-based subscripts the literature gives those
parameters. The propensity score rows are named `W_1` through `W_b` for
the columns of `W`.

## Examples

``` r
# A confounded binary treatment whose true effect on the outcome is -2.
set.seed(42)
n <- 500
W <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, 0.25 + 0.5 * W)
Y <- 5 + 2 * W - 2 * A + rnorm(n)

W_ps <- cbind(1, W) # Propensity score design matrix
V <- cbind(1, A) # Marginal structural model design matrix

psi <- function(theta) {
  ee_ipw_msm(
    theta,
    y = Y,
    A = A,
    W = W_ps,
    V = V,
    distribution = "normal",
    link = "identity"
  )
}

# theta holds the two marginal structural model coefficients, whose slope is
# the causal effect, followed by the two propensity score coefficients.
m <- m_estimate(stacked_equations = psi, init = rep(0, 4))
coef(m)
#> MSM alpha_0 MSM alpha_1         W_1         W_2 
#>    5.953823   -2.066997   -1.024504    2.234017 
```
