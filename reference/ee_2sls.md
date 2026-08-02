# Estimating equations for Two-Stage Least Squares (2SLS)

Estimates the causal effect using two-stage least squares for
instrumental variable analysis. The first stage regresses the treatment
on the instruments (and exogenous variables), and the second stage
regresses the outcome on predicted treatment (and exogenous variables).

## Usage

``` r
ee_2sls(theta, y, A, Z, W = NULL, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length `1 + b + 2c`, where `b` is the number of
  instruments in `Z` and `c` is the number of exogenous variables in
  `W`. The first `1 + c` parameters are for the second-stage model; the
  remainder are for the first-stage model.

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n observed treatment values.

- Z:

  Numeric n-by-b matrix of instrumental variable(s).

- W:

  Optional n-by-c matrix of exogenous variables included in both stages.
  Default `NULL`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A `(1+b+2c)`-by-n matrix of estimating equation contributions. The
second-stage rows are named `stage2_A` for the fitted treatment and
`stage2_W_1` through `stage2_W_c`; the first-stage rows are named
`stage1_Z_1` through `stage1_Z_b` and `stage1_W_1` through `stage1_W_c`.
The columns of `W` appear in both stages, so each name carries the stage
it belongs to.

## Examples

``` r
# A continuous treatment confounded by an unmeasured U, a strong instrument
# Z, and one measured exogenous covariate. The true causal effect is 2.
set.seed(42)
n <- 500
W1 <- rnorm(n)
Z <- cbind(rbinom(n, 1, 0.5))
U <- rnorm(n)
A <- 1.5 * Z[, 1] + 0.3 * W1 + U + rnorm(n, sd = 0.5)
Y <- 2 * A - U + 0.5 * W1 + rnorm(n)
W <- cbind(1, W1)

# At a first stage of all zeros the fitted treatment is identically zero, so
# the second stage degenerates: its leading design column vanishes and the
# coefficient on it, the causal effect, has nothing left to be estimated
# from. Seed the starting values with the ordinary least squares fit of each
# stage instead: the first stage regresses A on the instrument and the
# exogenous covariates, the second regresses Y on the fitted A and the same
# covariates.
alpha_init <- as.numeric(coef(lm(A ~ cbind(Z, W) - 1)))
a_hat <- as.numeric(cbind(Z, W) %*% alpha_init)
beta_init <- as.numeric(coef(lm(Y ~ cbind(a_hat, W) - 1)))

psi <- function(theta) ee_2sls(theta, y = Y, A = A, Z = Z, W = W)

# theta holds the three second-stage coefficients followed by the three
# first-stage coefficients.
m <- m_estimate(
  stacked_equations = psi,
  init = c(beta_init, alpha_init)
)

# theta_1 is the coefficient on the fitted treatment, the causal effect.
coef(m)
#>    stage2_A  stage2_W_1  stage2_W_2  stage1_Z_1  stage1_W_1  stage1_W_2 
#>  1.91139044  0.03769661  0.54088789  1.63627118 -0.06199778  0.30488707 
```
