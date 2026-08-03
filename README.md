
<!-- README.md is generated from README.Rmd. Please edit that file -->

# deli

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/r-causal/deli/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-causal/deli/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/r-causal/deli/graph/badge.svg)](https://app.codecov.io/gh/r-causal/deli)
<!-- badges: end -->

deli provides M-estimation and empirical sandwich variance estimation in
R.

M-estimators express a wide range of statistical procedures as the
solution to a set of estimating equations, and the empirical sandwich
estimator supplies their variance without further derivation. deli
offers a general interface for both custom and built-in estimating
equations, covering basic statistics, regression, causal inference,
survival analysis, measurement error, and pharmacokinetics.

deli is an R port of the Python
[delicatessen](https://github.com/pzivich/Delicatessen) library. The
[Translating from
Python](https://r-causal.github.io/deli/articles/translating-from-python.html)
article maps the Python interface onto deli and documents where the R
surface differs by convention.

## Installation

Install deli from CRAN with:

``` r
install.packages("deli")
```

Install the development version of deli from
[GitHub](https://github.com/r-causal/deli) with:

``` r
# install.packages("pak")
pak::pak("r-causal/deli")
```

## Example

The quickest way to fit a model is `m_estimate()`, which takes a
formula, a data frame, and a built-in estimating equation, then
constructs and solves the estimator in a single call. Here is a linear
regression on the `mtcars` data:

``` r
library(deli)

fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_regression,
  model = "linear"
)

fit
#> <MEstimator>
#>   Parameters: 3
#>   Observations: 32
#> Coefficients:
#> (Intercept): 37.2273
#> wt: -3.8778
#> hp: -0.0318
```

Standard errors come from the sandwich variance estimator, so the usual
accessors report robust inference without additional arguments:

``` r
summary(fit)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   270.5101
#> wt             -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    31.2310
#> hp             -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    19.1269
```

`vcov()` returns the sandwich variance-covariance matrix and `confint()`
returns Wald confidence intervals. deli also supplies
[broom](https://broom.tidymodels.org) tidiers, `tidy()`, `augment()`,
and `glance()`.

## Learn more

The [package website](https://r-causal.github.io/deli/) collects the
reference documentation and articles. Start with [Getting
Started](https://r-causal.github.io/deli/articles/getting-started.html),
then see the articles on [custom estimating
equations](https://r-causal.github.io/deli/articles/custom-estimating-equations.html),
[regression
models](https://r-causal.github.io/deli/articles/regression-models.html),
and [causal
inference](https://r-causal.github.io/deli/articles/causal-inference.html)
for worked examples.
