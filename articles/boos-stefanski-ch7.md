# Boos & Stefanski (2013): M-Estimation (Estimating Equations)

> **Note**
>
> This article is translated from the [Boos & Stefanski (2013):
> M-Estimation (Estimating Equations)
> example](https://deli.readthedocs.io/en/latest/Examples/Boos-Stefanski-Ch7.html)
> in the documentation of [delicatessen](https://deli.readthedocs.io/),
> deli’s Python counterpart.

``` r

library(deli)
```

Selected examples from Chapter 7 of Boos & Stefanski (2013), *Essential
Statistical Inference*, demonstrating M-estimation with the empirical
sandwich variance estimator.

## 7.2.2 Mean and Variance

Simultaneous estimation of the mean and variance using stacked
estimating equations.

### Data

``` r

set.seed(80950841)
n <- 200
y_mv <- rnorm(n, mean = 10, sd = sqrt(2))
x_mv <- rnorm(n, mean = 5, sd = 1)
```

### By-hand estimating equations

``` r

# By-hand mean and variance EE
psi <- function(theta) {
  rbind(
    matrix(y_mv - theta[1], nrow = 1),            # EE for the mean
    matrix((y_mv - theta[1])^2 - theta[2], nrow = 1) # EE for the variance
  )
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0))
estr@theta
#>  theta_1  theta_2 
#> 9.906689 2.405768
confint(estr)
#>            lower    upper
#> theta_1 9.691728 10.12165
#> theta_2 2.039566  2.77197
```

### Built-in `ee_mean_variance`

``` r

psi <- function(theta) {
  ee_mean_variance(theta, y = y_mv)
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0))
estr@theta
#>     mean variance 
#> 9.906689 2.405768
confint(estr)
#>             lower    upper
#> mean     9.691728 10.12165
#> variance 2.039566  2.77197
```

## 7.2.3 Ratio Estimator

Estimating the ratio \theta = E\[Y\] / E\[X\] via the single estimating
equation \psi_i(\theta) = Y_i - \theta X_i.

### Single-equation version

``` r

psi <- function(theta) {
  matrix(y_mv - theta[1] * x_mv, nrow = 1)
}

estr <- m_estimate(stacked_equations = psi, init = c(1))
estr@theta
#>  theta_1 
#> 2.002463
confint(estr)
#>            lower    upper
#> theta_1 1.929498 2.075427
```

### Three-equation version

Stacking mean of Y, mean of X, and the ratio:

``` r

psi <- function(theta) {
  rbind(
    matrix(y_mv - theta[1], nrow = 1),             # EE for E[Y]
    matrix(x_mv - theta[2], nrow = 1),             # EE for E[X]
    matrix(rep(theta[1] - theta[3] * theta[2], n), nrow = 1) # ratio constraint
  )
}

estr <- m_estimate(stacked_equations = psi, init = c(10, 5, 2))
estr@theta
#>  theta_1  theta_2  theta_3 
#> 9.906689 4.947253 2.002463
confint(estr)
#>            lower     upper
#> theta_1 9.691728 10.121650
#> theta_2 4.808608  5.085899
#> theta_3 1.929498  2.075427
```

## 7.2.4 Delta Method via Stacking

Rather than applying the delta method after estimation, additional rows
can be stacked to propagate uncertainty through transformations. Here we
add \sqrt{\sigma^2} and \log(\sigma^2) as additional parameters.

``` r

psi <- function(theta) {
  rbind(
    matrix(y_mv - theta[1], nrow = 1),                      # mean
    matrix((y_mv - theta[1])^2 - theta[2], nrow = 1),       # variance
    matrix(rep(theta[3] - sqrt(theta[2]), n), nrow = 1),     # sqrt(var)
    matrix(rep(theta[4] - log(theta[2]), n), nrow = 1)       # log(var)
  )
}

estr <- m_estimate(stacked_equations = psi, init = c(10, 2, 1, 1))
estr@theta
#>   theta_1   theta_2   theta_3   theta_4 
#> 9.9066890 2.4057679 1.5510538 0.8778691
confint(estr)
#>             lower     upper
#> theta_1 9.6917279 10.121650
#> theta_2 2.0395657  2.771970
#> theta_3 1.4330043  1.669103
#> theta_4 0.7256507  1.030088
```

## 7.2.6 Instrumental Variable

Correcting for measurement error using an instrumental variable.
Consider the model Y = \alpha + \beta X^\* + \epsilon where X^\* is
measured with error as X = X^\* + U. An instrument W satisfies
\text{Cov}(W, \epsilon) = 0.

### Data

``` r

set.seed(12345)
n_iv <- 500
xstar <- rnorm(n_iv, mean = 0, sd = 1)           # true X
e <- rnorm(n_iv, mean = 0, sd = 1)                # outcome error
u <- rnorm(n_iv, mean = 0, sd = 0.5)              # measurement error
w <- xstar + rnorm(n_iv, mean = 0, sd = 0.5)      # instrument
x_iv <- xstar + u                                  # mismeasured X
y_iv <- 1 + 0.5 * xstar + e                        # outcome
```

### IV estimator: two-equation formulation

Using (1, W) as instruments:

``` r

psi <- function(theta) {
  resid <- y_iv - theta[1] - theta[2] * x_iv
  rbind(
    matrix(resid, nrow = 1),           # instrument: constant
    matrix(resid * w, nrow = 1)        # instrument: W
  )
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0))
estr@theta
#>   theta_1   theta_2 
#> 1.0219527 0.4812642
confint(estr)
#>             lower    upper
#> theta_1 0.9313164 1.112589
#> theta_2 0.3815893 0.580939
```

### IV estimator: four-equation formulation

Including the mean of X and W to allow further inference:

``` r

psi <- function(theta) {
  resid <- y_iv - theta[1] - theta[2] * x_iv
  rbind(
    matrix(resid, nrow = 1),              # instrument: constant
    matrix(resid * w, nrow = 1),          # instrument: W
    matrix(x_iv - theta[3], nrow = 1),    # mean of X
    matrix(w - theta[4], nrow = 1)        # mean of W
  )
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0, 0, 0))
estr@theta
#>    theta_1    theta_2    theta_3    theta_4 
#> 1.02195273 0.48126418 0.06070097 0.07382810
confint(estr)
#>               lower     upper
#> theta_1  0.93131640 1.1125891
#> theta_2  0.38158935 0.5809390
#> theta_3 -0.03928175 0.1606837
#> theta_4 -0.02223276 0.1698890
```

## 7.4.1 Robust Location (Huber)

The Huber loss function down-weights large residuals using a tuning
constant k. The estimating equation is \psi_i(\theta) = h_k(Y_i -
\theta) where h_k is the Huber function.

### Data

``` r

set.seed(42)
y_rob <- rnorm(250)
```

### By-hand Huber robust mean

``` r

k <- 3

psi <- function(theta) {
  resid <- y_rob - theta[1]
  # Huber loss: clip residuals to [-k, k]
  clipped <- pmax(pmin(resid, k), -k)
  matrix(clipped, nrow = 1)
}

estr <- m_estimate(stacked_equations = psi, init = 0)
estr@theta
#>     theta_1 
#> -0.02039958
confint(estr)
#>              lower     upper
#> theta_1 -0.1410716 0.1002724
```

### Built-in `ee_mean_robust`

``` r

psi <- function(theta) {
  ee_mean_robust(theta, y = y_rob, k = 3, loss = "huber")
}

estr <- m_estimate(stacked_equations = psi, init = 0)
estr@theta
#>     theta_1 
#> -0.02039958
confint(estr)
#>              lower     upper
#> theta_1 -0.1410716 0.1002724
```

## 7.5.1 Linear Regression

Linear regression Y_i = \beta_0 + \beta_1 X_i + \beta_2 Z_i + \epsilon_i
with the estimating equation \psi_i(\beta) = (Y_i - X_i^T \beta) X_i.

### Data

``` r

set.seed(99)
n_lm <- 500
x_lm <- rnorm(n_lm)
z_lm <- rbinom(n_lm, 1, 0.5)
y_lm <- 0.5 + 2 * x_lm - z_lm + rnorm(n_lm)
X_lm <- cbind(1, x_lm, z_lm)
```

### By-hand linear regression

``` r

psi <- function(theta) {
  yhat <- X_lm %*% theta                 # predicted values (n-by-1)
  resid <- y_lm - yhat                   # residuals (n-by-1)
  t(X_lm * as.numeric(resid))            # score: p-by-n
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0, 0))
estr@theta
#>    theta_1    theta_2    theta_3 
#>  0.5342416  1.9374201 -0.9359291
confint(estr)
#>              lower      upper
#> theta_1  0.4121997  0.6562835
#> theta_2  1.8425613  2.0322790
#> theta_3 -1.1134217 -0.7584366
```

### Built-in `ee_regression`

``` r

psi <- function(theta) {
  ee_regression(theta, X = X_lm, y = y_lm, model = "linear")
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0, 0))
estr@theta
#>    theta_1    theta_2    theta_3 
#>  0.5342416  1.9374201 -0.9359291
confint(estr)
#>              lower      upper
#> theta_1  0.4121997  0.6562835
#> theta_2  1.8425613  2.0322790
#> theta_3 -1.1134217 -0.7584366
```

## 7.5.4 Robust Regression

Huber robust regression with tuning constant k = 1.345. The estimating
equation is \psi_i(\beta) = h_k(Y_i - X_i^T \beta) X_i.

### By-hand robust regression

``` r

k_rr <- 1.345

psi <- function(theta) {
  yhat <- X_lm %*% theta
  resid <- y_lm - as.numeric(yhat)
  # Huber loss applied to residuals
  clipped <- pmax(pmin(resid, k_rr), -k_rr)
  t(X_lm * clipped)
}

estr <- m_estimate(
  stacked_equations = psi,
  init = c(0, 0, 0),
  solver = "nleqslv"
)
estr@theta
#>    theta_1    theta_2    theta_3 
#>  0.5245726  1.9337397 -0.9112539
confint(estr)
#>             lower      upper
#> theta_1  0.402750  0.6463953
#> theta_2  1.834690  2.0327898
#> theta_3 -1.095599 -0.7269085
```

### Built-in `ee_robust_regression`

``` r

psi <- function(theta) {
  ee_robust_regression(theta, X = X_lm, y = y_lm,
                       model = "linear", k = 1.345, loss = "huber")
}

estr <- m_estimate(
  stacked_equations = psi,
  init = c(0, 0, 0),
  solver = "nleqslv"
)
estr@theta
#>    theta_1    theta_2    theta_3 
#>  0.5245726  1.9337397 -0.9112539
confint(estr)
#>             lower      upper
#> theta_1  0.402750  0.6463953
#> theta_2  1.834690  2.0327898
#> theta_3 -1.095599 -0.7269085
```

## 7.5.5 Generalized Linear Models

Using `ee_glm` with different distribution/link combinations to estimate
common epidemiologic effect measures from binary data.

### Data

``` r

set.seed(321)
n_glm <- 500
x_glm <- rnorm(n_glm)
pr <- plogis(-0.5 + x_glm)
y_glm <- rbinom(n_glm, 1, pr)
X_glm <- cbind(1, x_glm)
```

### Logistic regression (odds ratio)

``` r

psi <- function(theta) {
  ee_glm(theta, X = X_glm, y = y_glm,
         distribution = "binomial", link = "logit")
}

estr <- m_estimate(stacked_equations = psi, init = c(0, 0))
estr@theta
#>    theta_1    theta_2 
#> -0.6230037  1.0244966
confint(estr)
#>              lower     upper
#> theta_1 -0.8234004 -0.422607
#> theta_2  0.7965750  1.252418
```

### Log-binomial regression (risk ratio)

``` r

psi <- function(theta) {
  ee_glm(theta, X = X_glm, y = y_glm,
         distribution = "binomial", link = "log")
}

estr <- m_estimate(stacked_equations = psi, init = c(-1, 0))
estr@theta
#>   theta_1   theta_2 
#> -1.404725  1.185587
confint(estr)
#>              lower     upper
#> theta_1 -1.6444513 -1.164999
#> theta_2  0.9745417  1.396632
```

### Identity-binomial regression (risk difference)

``` r

psi <- function(theta) {
  ee_glm(theta, X = X_glm, y = y_glm,
         distribution = "binomial", link = "identity")
}

estr <- m_estimate(stacked_equations = psi, init = c(0.3, 0))
estr@theta
#>   theta_1   theta_2 
#> 0.3681685 0.2098780
confint(estr)
#>             lower    upper
#> theta_1 0.3333530 0.402984
#> theta_2 0.1828439 0.236912
```
