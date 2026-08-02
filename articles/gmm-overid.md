# GMM for Over-Identified Parameters

> **Note**
>
> This article is translated from the [GMM for Over-Identified
> Parameters
> example](https://deli.readthedocs.io/en/latest/Examples/GMM-OverID.html)
> in the documentation of [delicatessen](https://deli.readthedocs.io/),
> deli’s Python counterpart.

``` r

library(deli)
```

When there are more estimating equations than parameters, the system is
*over-identified*. The standard `MEstimator` requires a square system
(number of equations equals number of parameters). `GMMEstimator`
handles over-identified systems by minimizing a quadratic form of the
estimating equations, producing valid point estimates and sandwich
variance estimates.

This document walks through instrumental variable (IV) examples where
combining multiple instruments creates over-identification.

## IV Example 1: Separate instruments (just-identified)

Consider estimating the causal effect of treatment A on outcome Y using
two instrumental variables Z_1 and Z_2. When each instrument gets its
own parameter, the system is just-identified: 2 equations, 2 parameters.

### Data

``` r

# Generate data
set.seed(777)
n <- 500

# Covariate
W <- rbinom(n, size = 1, prob = 0.25)

# Two instruments
Z1 <- rnorm(n, mean = 0, sd = 0.5)
Z2 <- rnorm(n, mean = 0, sd = 0.5)

# Treatment depends on both instruments
A <- Z1 + Z2 + rnorm(n)

# Outcome with effect modification by W
Y <- 2 * A - W * A + rnorm(n)
```

### Estimation with separate parameters

Each instrument defines its own estimating equation with its own
parameter. This is a square system (2 equations, 2 parameters), so
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
solves it.

``` r

# Two instruments, two parameters (just-identified)
estr <- m_estimate(
  stacked_equations = function(theta) {
    # Z1 moment condition with beta_1
    ee1 <- matrix(Z1 * (Y - theta[1] * A), nrow = 1)
    # Z2 moment condition with beta_2
    ee2 <- matrix(Z2 * (Y - theta[2] * A), nrow = 1)
    rbind(ee1, ee2)
  },
  init = c(0, 0)
)

data.frame(
  Param = c("beta_1 (Z1)", "beta_2 (Z2)"),
  Coef = round(estr@theta, 4),
  LCL = round(confint(estr)[, 1], 4),
  UCL = round(confint(estr)[, 2], 4)
)
#>               Param   Coef    LCL    UCL
#> theta_1 beta_1 (Z1) 1.6420 1.4662 1.8178
#> theta_2 beta_2 (Z2) 1.6787 1.4704 1.8871
```

Each instrument gives a separate estimate of the effect. Both are valid,
but they use the data inefficiently since each equation only leverages
one instrument.

## IV Example 1: Combined instruments (over-identified)

Now we combine both instruments to estimate a single parameter \beta.
This creates an over-identified system: 2 equations, 1 parameter. Both
moment conditions E\[Z_1 (Y - \beta A)\] = 0 and E\[Z_2 (Y - \beta A)\]
= 0 must hold simultaneously for the same \beta.

### MEstimator fails

`MEstimator` requires a square system (number of equations = number of
parameters). With 2 equations and 1 parameter, it cannot solve the
system.

``` r

# Over-identified: 2 equations, 1 parameter
estr_fail <- MEstimator(
  stacked_equations = function(theta) {
    # Both instruments, single beta
    ee1 <- matrix(Z1 * (Y - theta[1] * A), nrow = 1)
    ee2 <- matrix(Z2 * (Y - theta[1] * A), nrow = 1)
    rbind(ee1, ee2)
  },
  init = 0
) |>
  estimate()
#> Error in `method(estimate, deli::MEstimator)`:
#> ! `stacked_equations` returned 2 estimating equations at the initial
#>   values, but `init` has 1 parameter.
#> ℹ M-estimation requires one estimating equation per parameter (a 1-by-n
#>   matrix).
```

### GMMEstimator works

`GMMEstimator` handles the over-identified case by minimizing a
quadratic form Q(\theta) = \bar{g}(\theta)^T W \bar{g}(\theta), where
\bar{g}(\theta) is the vector of average estimating equations and W is a
weight matrix. It iterates to find the optimal weight matrix.
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
builds one and solves it in a single call.

``` r

# Over-identified: 2 equations, 1 parameter
estr_gmm <- gmm_estimate(
  stacked_equations = function(theta) {
    # Both instruments, single beta
    ee1 <- matrix(Z1 * (Y - theta[1] * A), nrow = 1)
    ee2 <- matrix(Z2 * (Y - theta[1] * A), nrow = 1)
    rbind(ee1, ee2)
  },
  init = 0
)

data.frame(
  Param = "beta (combined)",
  Coef = round(estr_gmm@theta, 4),
  LCL = round(confint(estr_gmm)[, 1], 4),
  UCL = round(confint(estr_gmm)[, 2], 4)
)
#>                   Param   Coef    LCL UCL
#> theta_1 beta (combined) 1.6567 1.5135 1.8
```

By combining both instruments into a single estimating equation system,
`GMMEstimator` efficiently pools information from both Z_1 and Z_2. This
typically produces narrower confidence intervals compared to using
either instrument alone.

## IV Example 2: With transportability

Now consider a setting with two populations: an external study
population (S = 0) where instrumental variable analyses can be
conducted, and a target population (S = 1) where we want to transport
the results. We use inverse odds of sampling weights to transport the IV
estimates, combined with over-identified moment conditions.

### Data

``` r

set.seed(777)

# External population (S=0)
n0 <- 400
W0 <- rbinom(n0, size = 1, prob = 0.25)
Z1_0 <- rnorm(n0, mean = 0, sd = 0.5)
Z2_0 <- rnorm(n0, mean = 0, sd = 0.5)
A0 <- Z1_0 + Z2_0 + rnorm(n0)
Y0 <- 2 * A0 - W0 * A0 + rnorm(n0)
S0 <- rep(0, n0)

# Target population (S=1)
n1 <- 100
W1 <- rbinom(n1, size = 1, prob = 0.5)
Z1_1 <- rep(0, n1)
Z2_1 <- rep(0, n1)
A1 <- rnorm(n1)
Y1 <- rep(0, n1)
S1 <- rep(1, n1)

# Combine into single data frame
d <- data.frame(
  W = c(W0, W1),
  Z1 = c(Z1_0, Z1_1),
  Z2 = c(Z2_0, Z2_1),
  A = c(A0, A1),
  Y = c(Y0, Y1),
  S = c(S0, S1)
)
```

### Estimation with transportability

The estimating equations have three parameters:

- \beta: the IV effect estimate, transported to the target population
- \alpha_0, \alpha_1: logistic regression coefficients for the sampling
  model P(S = 1 \mid W)

The IV moment conditions are weighted by inverse odds of sampling
weights \frac{P(S=1 \mid W)}{P(S=0 \mid W)} to transport the estimates
from the external study to the target population.

``` r

# Extract vectors for use in estimating equations
S <- d$S
W_vec <- d$W
Z1_all <- d$Z1
Z2_all <- d$Z2
A_all <- d$A
Y_all <- d$Y
X_s <- cbind(1, W_vec)  # Design matrix for sampling model

estr_transport <- gmm_estimate(
  stacked_equations = function(theta) {
    beta <- theta[1]
    alpha <- theta[2:3]
    n <- length(Y_all)

    # Sampling model: P(S=1 | W)
    ee_s <- ee_regression(alpha, X = X_s, y = S, model = "logistic")

    # Inverse odds of sampling weights
    pr_s <- inverse_logit(as.numeric(X_s %*% alpha))
    iosw <- pr_s / (1 - pr_s)

    # Weighted IV moment conditions (only external population contributes)
    ee1 <- matrix((1 - S) * iosw * Z1_all * (Y_all - beta * A_all), nrow = 1)
    ee2 <- matrix((1 - S) * iosw * Z2_all * (Y_all - beta * A_all), nrow = 1)

    rbind(ee1, ee2, ee_s)
  },
  init = c(0, 0, 0)
)

data.frame(
  Param = c("beta (transported)", "alpha_0", "alpha_1"),
  Coef = round(estr_transport@theta, 4),
  LCL = round(confint(estr_transport)[, 1], 4),
  UCL = round(confint(estr_transport)[, 2], 4)
)
#>                      Param    Coef     LCL     UCL
#> theta_1 beta (transported)  1.4655  1.2613  1.6696
#> theta_2            alpha_0 -1.9085 -2.2169 -1.6001
#> theta_3            alpha_1  1.4005  0.9407  1.8603
```

This example has 4 estimating equations (2 IV moment conditions + 2
sampling model equations) but only 3 parameters (\beta, \alpha_0,
\alpha_1), making it over-identified. The `GMMEstimator` efficiently
combines both instruments while correctly propagating uncertainty from
the sampling model through the sandwich variance estimator.

## References

Hansen LP. (1982). Large sample properties of generalized method of
moments estimators. *Econometrica*, 50(4), 1029-1054.
