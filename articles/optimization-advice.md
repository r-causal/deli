# Optimization advice

> **Note**
>
> This article is translated from the [Optimization Advice
> page](https://deli.readthedocs.io/en/latest/Optimization%20Advice.html)
> in the documentation of [delicatessen](https://deli.readthedocs.io/),
> deli’s Python counterpart.

``` r

library(deli)
```

This article offers guidance for when root-finding or minimization is
difficult for a given set of estimating functions. The flexibility of
the generic estimating-equation interface comes at a cost: deli does not
provide the fastest or most specialized routines for any one model,
because it aims to solve them all. When a fit stalls, diverges, or
returns a variance full of `NaN`, the remedies below usually help. Each
is grounded in how the solve actually behaves.

The controls referenced throughout are `init`, `solver`, `maxiter`,
`tolerance`, `deriv_method`, and `dx`.
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
takes all of them in a single call. In the two-step form, `init` belongs
to
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
and travels with the object, while the rest are arguments of
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md).

## Choose good starting values

Optimization improves when the starting values in `init` are close to
the true parameter vector. We rarely know those values in advance, so
the practical rule is to start in a reasonable part of the parameter
space. Three properties matter:

1.  The starting values are numerically sensible.
2.  They lie inside the range of possible parameter values. A proportion
    is bound to \[0, 1\], so an initial value of -10 can break the
    search. A variance is non-negative, so it should not start below
    zero.
3.  They are as close to the estimate as available information allows.

For most regression coefficients, a starting value of zero is a
reasonable default in the absence of other information. Some estimating
equations are far more sensitive to the starting point than a linear
model is. Robust regression, the tobit model, and accelerated failure
time models can stall or diverge when started far from the solution,
because their estimating functions are flat or steep in regions away
from the optimum.

`init` is supplied when the estimator is constructed and travels with
the object. Constructing one without solving it and reading `init` back
off the object shows this directly:

``` r

data(robust_regress)
robust_regress$height_c <- robust_regress$height - mean(robust_regress$height)

X <- cbind(1, robust_regress$height_c)
y <- robust_regress$weight

psi_robust <- function(theta) {
  ee_robust_regression(theta, X = X, y = y, model = "linear", k = 1.345)
}

# Seed the intercept at the mean of the outcome and the slope at zero
m <- MEstimator(psi_robust, init = c(mean(y), 0))

m@init
#> [1] 65.61107  0.00000
```

## Seed from a simpler fit

When an estimator stacks several estimating functions, some parameters
can be estimated on their own, outside the joint solve. Fitting those
pieces first and passing the results as `init` lets the joint
optimization concentrate on the parameters that cannot be solved
separately. This is especially useful for regression models, where
dedicated and numerically stable fitting routines exist: seeding `init`
with the coefficients from
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) lets deli borrow the
strength of those routines.

A lighter version of the same idea does not require solving any
coefficients. For a regression model, set the intercept to the mean of
the outcome (or its link-transformed value for a generalized linear
model) and leave the remaining coefficients at zero. The seed in the
previous section does exactly this.

[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
and
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
also accept a `subset` argument, which holds some parameters fixed at
their initial values and solves only for the rest. When an estimator has
hundreds of parameters and many can be pinned down outside the
procedure, restricting the search to the remaining `subset` reduces the
problem the solver has to face. The indices are one-based in deli.

## Increase the iteration budget

If better starting values are not enough, raise the iteration budget.
`maxiter` defaults to `5000`; increasing it lets the solver take more
steps before giving up.

``` r

m_estimate(psi_robust, init = c(mean(y), 0), maxiter = 20000)
```

## Try a different solver

deli offers three built-in solvers for `MEstimator`, selected through
the `solver` argument. Each has different operating characteristics, and
a problem that defeats one may yield to another.

| `solver` | Routine | Method | Notes |
|----|----|----|----|
| `"rootSolve"` | \[rootSolve::multiroot()\] | Newton | Default. No extra package required. |
| `"nleqslv"` | \[nleqslv::nleqslv()\] | Newton with line search and trust region | Alternative global strategies. |
| `"lm"` | \[minpack.lm::nls.lm()\] | Levenberg-Marquardt | Damped least squares; robust to a singular Jacobian at the start. |

`"rootSolve"` is the default. It applies a Newton method and converges
on the great majority of well-conditioned square systems without any
additional package. `"nleqslv"` is an alternative Newton solver with
line-search and trust-region strategies that can succeed where a plain
Newton step overshoots. `"lm"` is the damped least-squares
Levenberg-Marquardt algorithm from MINPACK, the same routine Python
delicatessen calls by default; the damping makes it tolerant of an
ill-conditioned or singular Jacobian near the starting point, where the
Newton solvers stall.

The robust regression above is a concrete case. Started from the
mean-of-outcome seed, the Newton default makes almost no progress,
because near that point the Huber score clips most residuals and the
Jacobian is nearly singular. The Levenberg-Marquardt solver converges
from the same seed.

``` r

m_newton <- MEstimator(psi_robust, init = c(mean(y), 0)) |>
  estimate(solver = "rootSolve")
#> Warning: ! rootSolve did not converge within 5000 iterations.
#> ℹ The estimating functions are not solved at the returned values (achieved
#>   precision 15.5).
#> ℹ Results may be unreliable. Consider increasing `maxiter` or using the "lm"
#>   solver.
m_lm <- MEstimator(psi_robust, init = c(mean(y), 0)) |>
  estimate(solver = "lm")

coef(m_newton)
#>     theta_1     theta_2 
#> 64.82739739 -0.05702882
coef(m_lm)
#>   theta_1   theta_2 
#> 65.543970  0.619315
```

The Newton fit has barely moved from the seed, while the
Levenberg-Marquardt fit has reached the robust solution. The next
section shows how to tell the two outcomes apart without knowing the
answer in advance.

For `GMMEstimator`, the solve is a minimization rather than a
root-finding problem, and `solver` accepts any method name understood by
\[stats::optim()\]. The default is `"BFGS"`.

A custom solver is also accepted. Pass a function that takes
`stacked_equations` and `init` arguments and returns the solved
parameter vector; some problems have a bespoke algorithm that
outperforms the general-purpose routines.

## Scale the predictors

When covariates sit on very different numeric scales, the bread matrix
(the derivative of the estimating functions) becomes ill-conditioned,
which slows convergence and degrades the numerical stability of the
bread inversion that produces the variance. Centering and scaling the
predictors before fitting keeps the bread well-conditioned.

The effect is easy to see in the condition number of the bread, where a
larger value means a matrix that is closer to singular and harder to
invert reliably:

``` r

set.seed(1)
n <- 200
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1e4)
y_sim <- 1 + 0.5 * x1 + 1e-4 * x2 + rnorm(n)

# Raw predictors, on wildly different scales
X_raw <- cbind(1, x1, x2)
fit_raw <- m_estimate(
  function(theta) ee_regression(theta, X = X_raw, y = y_sim, model = "linear"),
  init = c(0, 0, 0)
)
kappa(fit_raw@bread, exact = TRUE)
#> [1] 119644040

# Standardized predictors
X_scaled <- cbind(1, scale(x1), scale(x2))
fit_scaled <- m_estimate(
  function(theta) ee_regression(theta, X = X_scaled, y = y_sim, model = "linear"),
  init = c(0, 0, 0)
)
kappa(fit_scaled@bread, exact = TRUE)
#> [1] 1.043349
```

Standardizing the predictors reduces the condition number by many orders
of magnitude. The point estimates and their standard errors are on the
scaled coefficients, so translate them back to the original units
afterward if the raw scale is the one you want to report.

## Choose a derivative method

The bread is a derivative, and deli computes it by the method named in
`deriv_method`.

| `deriv_method` | Meaning                                 |
|----------------|-----------------------------------------|
| `"capprox"`    | Central finite difference (default).    |
| `"fapprox"`    | Forward finite difference.              |
| `"bapprox"`    | Backward finite difference.             |
| `"exact"`      | Forward-mode automatic differentiation. |

`"capprox"`, the default, is a central finite difference. It is accurate
and robust, and it is the recommended starting point. `"fapprox"` and
`"bapprox"` are one-sided differences, useful when the estimating
function cannot be evaluated on one side of the current parameter value.
The `dx` argument sets the finite-difference step for all three; it
defaults to `1e-9` and is the one tuning knob these methods expose.

`dx` is an absolute step, so it interacts with the magnitude of the
estimates. The spacing of the doubles around an estimate grows with the
estimate, so `1e-9` spans fewer and fewer of them as a parameter grows,
and past roughly `1.7e7` adding it lands back on the same double and
leaves nothing to difference. deli floors the step at each parameter’s
floating-point resolution so that the derivative cannot collapse. That
floor engages as soon as `dx` spans fewer than ten to twenty thousand
representable values, the count depending on where the estimate sits
between consecutive powers of two, which at the default `dx` is any
estimate past about `450`: from there up the step actually applied is
the floor rather than `dx`, and the bread is accurate to a few parts in
ten thousand rather than to roughly `1e-7`. Rescaling the data, or
moving to `deriv_method = "exact"`, which takes no step at all, restores
full accuracy.

`"exact"` is forward-mode automatic differentiation. It carries
derivatives through the computation exactly, so it is the most accurate
option and has no step size to choose. Exact differentiation now works
for every built-in estimating equation except the two that are not
differentiable, `ee_percentile` and `ee_positive_mean_deviation`. It
requires that the estimating function be built from operations the
automatic differentiation supports; the pre-built equations satisfy
this, and custom equations do as long as they use supported operations.

``` r

m_exact <- m_estimate(
  psi_robust,
  init = c(mean(y), 0),
  solver = "lm",
  deriv_method = "exact"
)
sqrt(diag(vcov(m_exact)))
#>   theta_1   theta_2 
#> 0.5713395 0.3348466
```

The derivative method affects only the bread, and therefore the
variance; it does not change the point estimate.

## Diagnose convergence

deli does not print a running log while it solves, so convergence is
confirmed after the fit. Three checks cover most cases.

**The estimating functions should be near zero at the solution.** A root
of the estimating equations makes the estimating functions sum to
approximately zero across observations. The estimating-function matrix
is parameters-by-observations, so this sum is the row sums,
`rowSums(psi(...))`. Evaluating the stacked equations at the fitted
parameters and summing across observations reveals a solver that stopped
short:

``` r

# Stalled Newton fit from the solver comparison above
rowSums(psi_robust(coef(m_newton)))
#> [1]  5.253443 22.089524

# Converged Levenberg-Marquardt fit
rowSums(psi_robust(coef(m_lm)))
#> [1] 1.421085e-14 1.344671e-14
```

The Newton values are far from zero, which flags the stalled fit; the
Levenberg-Marquardt values are zero to numerical precision.

**A fit warns when the solver stops without solving the equations.**
Each built-in routine reports a status code, but the code on its own is
not the test. Every convergence test on offer is relative to something
of the solver’s own: the residuals it started from, the step it last
took, or an absolute bound on the function values that the scale of the
estimating functions can put out of reach. A solver can therefore report
success while sitting at the starting values. deli judges the returned
point as well, and warns when the estimating functions are not solved
there. A warning of that kind is a signal to revisit the starting
values, the scaling, or the solver.

A custom solver passed to `solver` is judged the same way. It reports no
status of its own, so the returned point is all there is to go on, and
it receives exactly the reading a built-in solver’s point receives.

[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
is judged the same way. It minimizes a quadratic in the summed moments
rather than root-finding, and
[`optim()`](https://rdrr.io/r/stats/optim.html) reports convergence on
the flat tail of that objective, well before the moments reach zero.

The judgement has two parts. The first asks whether any single equation
is in a state no root can be in. Where the per-observation contributions
of an equation carry both signs, that means they barely cancel. Where
they all carry one sign they cannot cancel at all, which is how a
relation among the parameters behaves, so what must vanish instead is
the mean contribution, measured against the scale of the problem. This
is the check most fits that fail meet, and on a
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
fit the message names the equation at fault.

The second part measures the Newton step the returned point would still
take towards a root and warns when that step is larger than the
estimates themselves. It is what catches a fit built from a design
matrix, where the contributions are mixed-sign and so cancel well
however wrong the parameters are. It is measured with the bread, so an
estimating equation whose Jacobian does not exist, such as a median, is
left unjudged by it, as is a fit whose finite-difference bread has
collapsed to round-off.

An over-identified
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
fit cannot drive every moment to zero, and a `subset` one holds some
parameters fixed while still summing every equation into the objective,
so neither is judged either way; for those, inspect
[`rowSums()`](https://rdrr.io/r/base/colSums.html) of the estimating
functions as above.

**The variance should be finite and of a reasonable magnitude.** A
singular or near-singular bread produces a variance full of `NaN`,
`Inf`, or implausibly large entries. Inspecting the standard errors
catches this:

``` r

sqrt(diag(vcov(m_lm)))
#>   theta_1   theta_2 
#> 0.5713376 0.3348466
```

If these are `NaN` or enormous, the bread could not be inverted
reliably. Scaling the predictors, improving the starting values, or
checking the model specification for a redundant parameter usually
resolves it.

## A word on tolerance

The `tolerance` argument controls how close to a root the solver must
get before it stops; it defaults to `1e-9`. It is tempting to loosen
`tolerance` to force a fit to succeed, but this only permits a larger
optimization error: the solver stops further from the true root of the
estimating equation, and the point estimate and variance inherit that
error. Do not loosen `tolerance` to get a fit to converge unless you are
certain the looser value is within acceptable error for your problem.
Better starting values, scaling, or a different solver address the cause
rather than the symptom.
