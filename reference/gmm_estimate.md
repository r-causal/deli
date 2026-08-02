# One-step GMM estimation

Creates a `GMMEstimator` and estimates it in one call. Supports both a
formula interface and a function interface, parallel to
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md).

## Usage

``` r
gmm_estimate(stacked_equations, ...)

# S3 method for class 'formula'
gmm_estimate(
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
  allow_pinv = TRUE,
  overid_maxiter = 200L,
  overid_tolerance = 1e-09
)

# Default S3 method
gmm_estimate(
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
  allow_pinv = TRUE,
  overid_maxiter = 200L,
  overid_tolerance = 1e-09
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
  solved. The objective is a quadratic form in every moment condition
  and `subset` changes only which parameters are free to move within it,
  so the conditions outside the subset are still summed in and still
  pull on the free parameters. A subset fit is therefore not the fit of
  the subset equations on their own, which is what
  [`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
  and
  [`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
  return, and the same stack and the same `subset` give the two
  different values. The variance estimator ignores `subset`.

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

- overid_maxiter:

  Integer maximum iterations for the two-step iterative procedure for
  over-identified problems. Default `200L`. The update converges
  linearly rather than quadratically, so a well-identified system
  commonly needs tens of passes to reach `overid_tolerance` and a weakly
  identified one can need hundreds.

- overid_tolerance:

  Numeric tolerance for convergence of the two-step iterative procedure.
  Default `1e-9`.

## Value

A fitted `GMMEstimator` object.

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

## Moment quality of an over-identified fit

A just-identified system has as many moment conditions as parameters, so
the moments vanish at a solution and the size of what is left over says
whether the fit succeeded. An over-identified system has no such
reading: no value of the parameters drives every condition to zero, and
a residual moment is expected rather than diagnostic. Hansen's
J-statistic is the reading that is available there. It is n times the
GMM objective at the minimum, \\J = n \bar{g}(\hat{\theta})' W
\bar{g}(\hat{\theta})\\, where \\\bar{g}\\ averages the moment
conditions over the observations and \\W\\ is the weight matrix the fit
finished with. Under correct specification it is asymptotically
chi-squared on as many degrees of freedom as the system has moment
conditions beyond parameters, so its size can be judged against a
reference distribution rather than against the scale of the data.

[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
records it in the `j_statistic` property of an over-identified fit, and
[`summary()`](https://r-causal.github.io/deli/reference/deli-display.md)
reports it with its degrees of freedom and its P-value. A
just-identified fit has no degrees of freedom left over and leaves the
property `NULL`; its moments are judged directly instead, as
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
describes. A `subset` fit holds the parameters outside the subset at
their initial values rather than estimating them, which the reference
distribution does not allow for, so it is left `NULL` too.

A P-value the reference distribution all but rules out warns with the
class `deli_gmm_moments_rejected`, which usually means the moment
conditions cannot all hold at one value of the parameters. The weight
matrix is what makes J comparable across problems, so the warning is
raised only where the two-step update settled: a fit that exhausted
`overid_maxiter` has already warned about that, and its J has no
reference distribution to be judged against. The property still records
the statistic in that case, as it does for `overid_maxiter = 0`, which
leaves the identity weight matrix in place and so leaves J an
unstandardized sum of squared moments.

The reading J cannot make is the opposite failure. Moment conditions
that are linearly dependent, one of them repeating what the others
already say, leave the covariance the weight matrix inverts singular,
and the update falls through to the pseudo-inverse; the fit that comes
back is the fit of the independent conditions alone. J is silent about
it, because a condition the others account for agrees with them wherever
the parameters sit and so adds nothing for J to measure, which drives J
toward zero rather than away from it. That case warns with the class
`deli_gmm_moments_dependent` instead, naming the conditions the
factorization found redundant.

## Examples

``` r
# Two instruments for a single treatment effect, confounded by an
# unmeasured U. Two moment conditions for one parameter leave the system
# over-identified, which is the case GMM is for: `m_estimate()` requires one
# estimating equation per parameter.
set.seed(42)
n <- 200
d <- data.frame(Z1 = rbinom(n, 1, 0.5), Z2 = rnorm(n))
U <- rnorm(n)
d$A <- 0.5 * d$Z1 + 0.3 * d$Z2 + U + rnorm(n)
d$Y <- 2 * d$A - U + rnorm(n)

# One moment condition per instrument: an instrument should be uncorrelated
# with the residual of the outcome on the treatment.
ee_iv_moments <- function(theta, X, y, Z) {
  t(Z * (y - as.numeric(X %*% theta)))
}

# The formula interface reads the design and the starting values off the
# model and passes anything else, here the instruments, on to `.ee`.
g <- gmm_estimate(
  Y ~ A - 1,
  data = d,
  .ee = ee_iv_moments,
  Z = cbind(Z1, Z2)
)

# The instruments move the estimate toward the treatment effect of 2 that
# generated the data. The least-squares fit of Y on A ignores the confounding
# and stays further from it.
coef(g)
#>        A 
#> 1.837705 
coef(lm(Y ~ A - 1, data = d))
#>        A 
#> 1.585141 

# The function interface takes a `stacked_equations` closure instead, for a
# system no formula describes. It has no design to read parameter names from,
# so names on `init` are what label the results.
psi_iv <- function(theta) {
  residual <- d$Y - theta[1] * d$A
  rbind(d$Z1 * residual, d$Z2 * residual)
}

g2 <- gmm_estimate(
  stacked_equations = psi_iv,
  init = c(effect = 0)
)
coef(g2)
#>   effect 
#> 1.837705 
```
