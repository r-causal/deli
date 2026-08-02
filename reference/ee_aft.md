# Estimating equation for accelerated failure time models

Returns a p-by-n matrix of estimating equation contributions for
accelerated failure time (AFT) models. Supports Weibull, exponential,
log-logistic, and log-normal distributions.

## Usage

``` r
ee_aft(theta, X, time, event, distribution, weights = NULL, offset = NULL)
```

## Arguments

- theta:

  Numeric vector of length `b + 1` (or `b` for exponential). The first
  `b` elements are regression coefficients and the last element is
  `log(1/sigma)` (the log-inverse scale).

- X:

  Numeric n-by-b design matrix (should include intercept column).

- time:

  Numeric vector of n observed (possibly censored) times.

- event:

  Numeric vector of n event indicators (1 = event, 0 = censored).

- distribution:

  Character string: `"weibull"`, `"exponential"`, `"log-logistic"`, or
  `"log-normal"`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix where p is the number of parameters. For the exponential
distribution, whose parameters are all regression coefficients, the rows
are unnamed. For every other distribution the regression rows are named
`X_1` through `X_b` for the columns of `X` and the final row is named
`log_inv_scale`.

## Details

The AFT model uses the parameterization where \\Z_i = (\log(t_i) - X_i
\beta) / \sigma\\ and the last element of theta is \\\log(1/\sigma)\\.

## Examples

``` r
# A Weibull AFT model for times generated from one binary covariate, with
# some of the times right censored by an independent exponential.
set.seed(1)
n <- 200
x <- rbinom(n, 1, 0.5)
Xd <- cbind(1, x)
eps <- log(-log(runif(n)))
t_event <- exp(2 + 0.5 * x + 0.8 * eps)
t_censor <- rexp(n, rate = 0.02)
t_obs <- pmin(t_event, t_censor)
delta <- as.numeric(t_event <= t_censor)

psi <- function(theta) {
  ee_aft(theta, X = Xd, time = t_obs, event = delta, distribution = "weibull")
}

# The default rootSolve solver does not converge here, so nleqslv is used.
# The last parameter is log(1/sigma), not a regression coefficient.
m <- m_estimate(
  stacked_equations = psi,
  init = c(mean(log(t_obs)), 0, 0),
  solver = "nleqslv"
)
coef(m)
#>           X_1           X_2 log_inv_scale 
#>     2.0952700     0.3677968     0.3182687 
```
