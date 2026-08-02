# Causal Inference

``` r

library(deli)
```

## Overview

deli provides estimating equations for several causal inference methods.
A key advantage of the M-estimation framework is that all nuisance model
parameters (e.g., propensity scores, outcome models) are estimated
jointly, so the sandwich variance correctly propagates uncertainty from
all stages.

This vignette demonstrates:

- G-formula (g-computation / standardization)
- Inverse probability weighting (IPW)
- Augmented IPW (AIPW / doubly robust)

## Simulated data

We simulate data with a binary treatment `A`, a confounder `W`, and a
continuous outcome `Y`:

``` r

set.seed(42)
n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, 1, 0.4)
A <- rbinom(n, 1, plogis(-0.5 + 0.5 * W1 + 0.3 * W2))
# True ATE = 1.5
Y <- 2 + 1.5 * A + W1 - 0.5 * W2 + rnorm(n)
```

## G-formula (g-computation)

The g-formula estimates the average causal effect by fitting an outcome
model and standardizing (averaging predictions over the covariate
distribution):

``` r

X <- cbind(1, A, W1, W2)      # Observed design matrix
X1 <- cbind(1, 1, W1, W2)     # Counterfactual: all treated
X0 <- cbind(1, 0, W1, W2)     # Counterfactual: all untreated

psi_gformula <- function(theta) {
  ee_gformula(theta, y = Y, X = X, X1 = X1, X0 = X0)
}

# theta: ACE, E[Y(1)], E[Y(0)], beta0, beta1, beta2, beta3
m_gformula <- m_estimate(
  stacked_equations = psi_gformula,
  init = c(0, 0, 0, 0, 0, 0, 0)
)

# ACE estimate (true = 1.5)
m_gformula@theta[1]
#>      ACE 
#> 1.578994

# With confidence interval
summary(m_gformula)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 1000
#> Parameters: 7
#> 
#>          Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> ACE        1.5790     0.0668    23.6319     1.4480     1.7100     <2e-16   407.7400
#> E[Y^1]     3.3099     0.0591    55.9685     3.1940     3.4259     <2e-16        Inf
#> E[Y^0]     1.7309     0.0539    32.0888     1.6252     1.8367     <2e-16   748.0967
#> X_1        1.9658     0.0503    39.0490     1.8671     2.0644     <2e-16        Inf
#> X_2        1.5790     0.0668    23.6319     1.4480     1.7100     <2e-16   407.7400
#> X_3        0.9694     0.0328    29.5917     0.9052     1.0337     <2e-16   636.8753
#> X_4       -0.5067     0.0660    -7.6826    -0.6360    -0.3775   1.56e-14    45.8669
```

The first parameter is the average causal effect (ACE).

## Inverse probability weighting (IPW)

IPW reweights the observed data by the inverse of the treatment
probability (propensity score):

``` r

W_ps <- cbind(1, W1, W2)  # Propensity score design matrix

psi_ipw <- function(theta) {
  ee_ipw(theta, y = Y, A = A, W = W_ps)
}

# theta: ACE, E[Y(1)], E[Y(0)], alpha0, alpha1, alpha2
m_ipw <- m_estimate(stacked_equations = psi_ipw, init = c(0, 0, 0, 0, 0, 0))

# ACE
m_ipw@theta[1]
#>      ACE 
#> 1.575945
summary(m_ipw)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 1000
#> Parameters: 6
#> 
#>          Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> ACE        1.5759     0.0671    23.4958     1.4445     1.7074     <2e-16   403.1047
#> E[Y^1]     3.3135     0.0591    56.0548     3.1977     3.4294     <2e-16        Inf
#> E[Y^0]     1.7376     0.0551    31.5497     1.6296     1.8455     <2e-16   723.3260
#> W_1       -0.3774     0.0865    -4.3625    -0.5470    -0.2079   1.29e-05    16.2468
#> W_2        0.4685     0.0687     6.8215     0.3339     0.6031   9.01e-12    36.6917
#> W_3        0.0997     0.1333     0.7476    -0.1616     0.3609      0.455     1.1370
```

### Truncating propensity scores

Propensity scores near zero or one produce very large weights, and a
handful of large weights can dominate an IPW estimate. `truncate` clips
the fitted scores to a range before they are inverted. Whether that
changes anything depends on the scores themselves, so look at them
first:

``` r

# Propensity scores implied by the alphas in the fitted stack
ps <- plogis(drop(W_ps %*% m_ipw@theta[4:6]))
range(ps)
#> [1] 0.1350166 0.7790256

# How many scores a tight c(0.3, 0.7) range would clip, by tail
c(lower = sum(ps < 0.3), upper = sum(ps > 0.7))
#> lower upper 
#>   157     9
```

Nothing here is close to zero or one, so a conventional range such as
`c(0.1, 0.9)` would leave every score untouched and the fit unchanged. A
deliberately tight range shows what truncation does. `c(0.3, 0.7)` clips
166 of the 1000 scores, all but nine of them in the lower tail:

``` r

psi_ipw_trunc <- function(theta) {
  ee_ipw(theta, y = Y, A = A, W = W_ps,
         truncate = c(0.3, 0.7))
}

m_ipw_trunc <- m_estimate(
  stacked_equations = psi_ipw_trunc,
  init = c(0, 0, 0, 0, 0, 0)
)
summary(m_ipw_trunc)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 1000
#> Parameters: 6
#> 
#>          Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> ACE        1.5457     0.0700    22.0790     1.4085     1.6829     <2e-16   356.4365
#> E[Y^1]     3.2804     0.0641    51.1496     3.1547     3.4060     <2e-16        Inf
#> E[Y^0]     1.7346     0.0554    31.2919     1.6260     1.8433     <2e-16   711.6262
#> W_1       -0.3774     0.0865    -4.3625    -0.5470    -0.2079   1.29e-05    16.2468
#> W_2        0.4685     0.0687     6.8215     0.3339     0.6031   9.01e-12    36.6917
#> W_3        0.0997     0.1333     0.7476    -0.1616     0.3609      0.455     1.1370
```

Both the ACE and its standard error move relative to the untruncated
fit:

``` r

data.frame(
  fit = c("Untruncated", "Truncated"),
  estimate = round(c(m_ipw@theta[[1]], m_ipw_trunc@theta[[1]]), 4),
  std_err = round(
    sqrt(c(vcov(m_ipw)[1, 1], vcov(m_ipw_trunc)[1, 1])),
    4
  )
)
#>           fit estimate std_err
#> 1 Untruncated   1.5759  0.0671
#> 2   Truncated   1.5457  0.0700
```

Truncation is usually motivated as a variance reduction, so the rise in
the standard error is worth dwelling on. Clipping pays for itself only
when the weights it caps are genuinely extreme, which is not the case
here. It also costs precision that joint estimation would otherwise
supply: a clipped score no longer moves with the propensity score
parameters, so the sandwich variance credits less of the gain from
estimating those parameters rather than treating them as known. Compare
a truncated fit against the untruncated one instead of reaching for a
range by default.

## Augmented IPW (doubly robust)

AIPW combines outcome modeling and IPW. It is “doubly robust”: the
estimate is consistent if either the outcome model or the propensity
score model is correctly specified. It reuses both sets of design
matrices built above:

``` r

psi_aipw <- function(theta) {
  ee_aipw(theta, y = Y, A = A,
          W = W_ps,             # Propensity score model
          X = X, X1 = X1, X0 = X0)  # Outcome model
}

# theta: ACE, E[Y(1)], E[Y(0)], alpha (3), beta (4)
m_aipw <- m_estimate(
  stacked_equations = psi_aipw,
  init = c(0, 0, 0, rep(0, 3), rep(0, 4))
)

summary(m_aipw)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 1000
#> Parameters: 10
#> 
#>          Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> ACE        1.5750     0.0666    23.6451     1.4445     1.7056     <2e-16   408.1897
#> E[Y^1]     3.3069     0.0595    55.6246     3.1904     3.4235     <2e-16        Inf
#> E[Y^0]     1.7319     0.0542    31.9673     1.6257     1.8381     <2e-16   742.4752
#> W_1       -0.3774     0.0865    -4.3625    -0.5470    -0.2079   1.29e-05    16.2468
#> W_2        0.4685     0.0687     6.8215     0.3339     0.6031   9.01e-12    36.6917
#> W_3        0.0997     0.1333     0.7476    -0.1616     0.3609      0.455     1.1370
#> X_1        1.9658     0.0503    39.0490     1.8671     2.0644     <2e-16        Inf
#> X_2        1.5790     0.0668    23.6319     1.4480     1.7100     <2e-16   407.7400
#> X_3        0.9694     0.0328    29.5917     0.9052     1.0337     <2e-16   636.8753
#> X_4       -0.5067     0.0660    -7.6826    -0.6360    -0.3775   1.56e-14    45.8669
```

## Comparing methods

The g-formula, IPW, and AIPW fits above all target the same average
causal effect. Both the outcome model and the propensity score model are
correctly specified in this simulation, so the estimates agree closely
with one another and with the true value:

``` r

data.frame(
  method = c("G-formula", "IPW", "AIPW"),
  estimate = round(
    c(m_gformula@theta[[1]], m_ipw@theta[[1]], m_aipw@theta[[1]]),
    3
  ),
  true_ate = 1.5
)
#>      method estimate true_ate
#> 1 G-formula    1.579      1.5
#> 2       IPW    1.576      1.5
#> 3      AIPW    1.575      1.5
```

## Further reading

- [`?ee_ipw_msm`](https://r-causal.github.io/deli/reference/ee_ipw_msm.md):
  Marginal structural models with IPW
- [`?ee_gestimation_snmm`](https://r-causal.github.io/deli/reference/ee_gestimation_snmm.md):
  G-estimation for structural nested mean models
- [`?ee_iv_causal`](https://r-causal.github.io/deli/reference/ee_iv_causal.md):
  Instrumental variable estimation
- [`?ee_2sls`](https://r-causal.github.io/deli/reference/ee_2sls.md):
  Two-stage least squares
- [`?ee_mean_sensitivity_analysis`](https://r-causal.github.io/deli/reference/ee_mean_sensitivity_analysis.md):
  Sensitivity analysis for unmeasured confounding
