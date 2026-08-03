# One-step M-estimation

Creates an `MEstimator` and estimates it in one call, analogous to how
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) creates and fits a
model in a single step. Supports both a formula interface (for
regression-family estimating equations) and a function interface (for
custom estimating equations).

## Usage

``` r
m_estimate(stacked_equations, ...)

# S3 method for class 'formula'
m_estimate(
  stacked_equations,
  data,
  .ee,
  ...,
  init = NULL,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-09,
  deriv_method = "capprox",
  dx = 1e-09,
  allow_pinv = TRUE
)

# Default S3 method
m_estimate(
  stacked_equations,
  ...,
  init,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-09,
  deriv_method = "capprox",
  dx = 1e-09,
  allow_pinv = TRUE
)
```

## Arguments

- stacked_equations:

  A formula or a function. When a formula, `data` and `.ee` must also be
  provided. When a function, it should take a numeric vector `theta` and
  return a p-by-n matrix, whose row names name the parameters when
  `init` has none; see
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md).
  A two-level factor or character response is converted to a 0/1
  indicator against its first level; an
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the formula
  is passed to `.ee` through its `offset` argument.

- ...:

  For the formula interface, additional arguments passed to `.ee`. These
  are evaluated with tidy evaluation in the context of `data`, so column
  names can be used directly (e.g., `event = status`). If the model
  frame drops rows for missing data, any such argument that spans the
  full data is subset to the same rows so it stays aligned with the
  design matrix and response. The function interface forwards nothing,
  so it requires `...` to be empty.

- data:

  A data frame (required when `stacked_equations` is a formula).

- .ee:

  An estimating equation function that accepts `theta`, `X`, and the
  response as its third argument, plus optionally additional arguments
  (required when `stacked_equations` is a formula). The formula response
  is passed positionally, so it reaches whatever the function calls that
  argument (`y` for
  [ee_regression](https://r-causal.github.io/deli/reference/ee_regression.md)
  or [ee_glm](https://r-causal.github.io/deli/reference/ee_glm.md),
  `time` for
  [ee_aft](https://r-causal.github.io/deli/reference/ee_aft.md)). An
  equation whose arguments leave any of those nowhere to go cannot be
  driven by a formula and is refused before anything is estimated:
  [ee_survival_model](https://r-causal.github.io/deli/reference/ee_survival_model.md)
  takes no design matrix, so it is fitted through the function interface
  instead.

- init:

  Numeric vector of initial parameter values. When `NULL` (default) and
  using the formula interface, a zero vector with names from the model
  matrix columns is generated automatically. Names on it label the
  parameters and take precedence over the row names of
  `stacked_equations`. An explicit `init` with no names of its own takes
  the same model matrix names on the formula interface, together with
  the parameter the estimating equation estimates beyond the design
  coefficients (`log_shape` for
  [ee_glm](https://r-causal.github.io/deli/reference/ee_glm.md) with
  `"gamma"`, `log_dispersion` with `"negative_binomial"`,
  `log_inv_scale` for a non-exponential
  [ee_aft](https://r-causal.github.io/deli/reference/ee_aft.md),
  `log_sigma` for
  [ee_tobit](https://r-causal.github.io/deli/reference/ee_tobit.md), and
  `log_phi` for
  [ee_beta_regression](https://r-causal.github.io/deli/reference/ee_beta_regression.md)).
  An `init` of any other length is left unnamed, which passes the
  labeling to the row names of the estimating functions; see
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  for that channel and for when the parameters are numbered instead.

- subset:

  Integer vector of parameter indices to solve for, or `NULL` (default)
  to solve for all parameters. Indices are 1-based; parameters not
  listed are held fixed at their `init` values while the rest are
  solved. The equations outside the subset are set aside along with the
  parameters they estimate, so the subset parameters are the root of the
  subset equations alone and the rest of the stack has no say in where
  they land: give a three-equation linear regression stack `subset = 1L`
  and the intercept comes back as the mean of the response less what the
  slopes held at their `init` values account for, because the first
  equation on its own is the estimating equation for a mean. Held at
  zero, which is what an unset `init` usually means, they account for
  nothing and the intercept is the mean of the response itself.
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  and
  [`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
  read the argument differently, since the GMM objective sums every
  equation whether the subset lists it or not, so the same stack and the
  same `subset` give the two different values. The variance estimator
  ignores `subset`.

- finite_correction:

  Character string for finite-sample correction (e.g., `"HC1"`), or
  `NULL` (default) for no correction. When set, the meat matrix is
  rescaled and inference switches to the t-distribution with
  `df = n_obs - n_params`. Passed straight through to the estimator
  constructor.

- solver:

  Character string or function for the solver. Default `NULL` uses
  `"rootSolve"` for M-estimation and `"BFGS"` for GMM. See
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  for the full list of solvers and for how the returned point is judged
  against the estimating equations.

- maxiter:

  Integer maximum iterations (default 5000). Must be a single positive
  whole number.

- tolerance:

  Numeric convergence tolerance (default 1e-9).

- deriv_method:

  Character string for numerical differentiation method (default
  `"capprox"`).

- dx:

  Numeric step size for differentiation (default 1e-9). Must be a single
  positive finite number, which is checked whichever `deriv_method` is
  in force. The step is absolute and is floored at the floating-point
  resolution of each estimate, so a large parameter magnitude cannot
  silently reduce it to nothing; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

- allow_pinv:

  Logical. Use pseudo-inverse if bread is singular? Default `TRUE`.

## Value

A fitted `MEstimator` object with populated `theta`, `variance`, etc.
Use
[`coef()`](https://r-causal.github.io/deli/reference/deli-generics.md),
[`vcov()`](https://r-causal.github.io/deli/reference/deli-generics.md),
[`confint()`](https://r-causal.github.io/deli/reference/deli-generics.md),
[`summary()`](https://r-causal.github.io/deli/reference/deli-display.md),
or [`tidy()`](https://r-causal.github.io/deli/reference/deli-tidiers.md)
to extract results.

## Details

Both interfaces place `...` ahead of `init`, `subset`, `tolerance`, and
the rest of the settings, so each of those must be named in full: R does
not partially match a supplied name against an argument that follows
`...`.
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
and the inference generics take `...` last, so they still accept R's
usual abbreviations.

What becomes of a name that matches no argument depends on the
interface. The function interface has no estimating equation to forward
`...` to, so it requires `...` to be empty and reports an unrecognized
name as an error rather than silently ignoring it. The formula interface
forwards `...` to `.ee` and matches each name against that function's
arguments exactly, naming the argument a refused name was probably meant
for. A name that merely abbreviates one is refused too, rather than
partially matched to it, since a fit that quietly took a misspelling for
`weights` reports different numbers and says nothing about it. An `.ee`
that takes `...` of its own accepts any name.

## Examples

``` r
# Formula interface
m <- m_estimate(mpg ~ wt + hp, data = mtcars,
                .ee = ee_regression, model = "linear")
coef(m)
#> (Intercept)          wt          hp 
#> 37.22727012 -3.87783074 -0.03177295 
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   270.5102
#> wt             -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    31.2310
#> hp             -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    19.1270

# Function interface
y <- c(1, 2, 3, 4, 5)
m2 <- m_estimate(
  stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
  init = c(mean = 0)
)
coef(m2)
#> mean 
#>    3 
```
