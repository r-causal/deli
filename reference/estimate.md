# Estimate parameters and sandwich variance

Solves the estimating equations for the parameter vector `theta` and
computes the empirical sandwich variance estimator.

## Usage

``` r
estimate(
  object,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-09,
  deriv_method = "capprox",
  dx = 1e-09,
  allow_pinv = TRUE,
  ...
)

## S7 method for class <deli::GMMEstimator>
estimate(
  object,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-09,
  deriv_method = "capprox",
  dx = 1e-09,
  allow_pinv = TRUE,
  ...
)
```

## Arguments

- object:

  An `MEstimator` or `GMMEstimator` object.

- solver:

  Character string specifying the solver algorithm, or a custom
  function. When `NULL` (default), uses `"rootSolve"` for `MEstimator`
  ([`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html))
  and `"BFGS"` for `GMMEstimator`
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html)). Other
  options for `MEstimator`: `"lm"`, the Levenberg-Marquardt algorithm
  ([`minpack.lm::nls.lm()`](https://rdrr.io/pkg/minpack.lm/man/nls.lm.html)),
  which mirrors the default solver of Python `delicatessen`
  (`scipy.optimize.root(method = "lm")`); and `"nleqslv"` (uses
  [`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)).
  A custom function must accept `stacked_equations` and `init` arguments
  and return the solved theta vector. Whichever solver is used, the
  point it returns is judged against the estimating equations themselves
  and a warning is raised when they are not solved there. A custom
  function reports no status of its own, so its point is judged exactly
  as a built-in solver's is. Two `GMMEstimator` fits are exceptions. An
  over-identified one cannot drive every moment to zero, so its moments
  are read against the J-statistic instead, which
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  describes. A `subset` one is judged neither way and does not warn,
  because what it minimizes is neither driven to zero nor read against
  that statistic; see `subset` in
  [`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
  for what such a fit estimates, and inspect
  [`rowSums()`](https://rdrr.io/r/base/colSums.html) of the estimating
  functions at the returned values.
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  cannot run inside itself, so a fit whose estimating function fits a
  second M-estimator cannot leave both on the default solver; that is
  refused rather than attempted, and naming `"nleqslv"` for either of
  the two fits resolves it.

  That refusal covers the solves deli makes and cannot cover a solve
  made inside a custom function, whose body deli does not read. **A
  custom solver that calls
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  itself must not be used while another `rootSolve` solve is running.**
  The inner call overwrites the environment the outer one is still
  using, and neither reports it: the outer solve was measured returning
  the inner fit's estimates under nothing louder than a generic
  non-convergence warning. Give either the outer fit or the custom
  solver a different algorithm, `"nleqslv"` or `"lm"`, so that no two
  `rootSolve` solves are open at once.

  What every solver is given is the estimating equations summed across
  the observations, since a root of the system is a root of those sums,
  and the bread is differentiated from the same reduction. A
  `GMMEstimator` minimizes a quadratic form in them instead and reduces
  them just as often. Both reductions are derived from the full p-by-n
  return unless the estimator carries a `summed_equations` property,
  which replaces them with the closed form the caller supplied and
  leaves the full evaluation to the meat and to the validation at the
  starting values; see
  [`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
  and
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  for what each class reduces and what it cannot.

  A `"rootSolve"` solve discards anything printed to standard output
  while it runs, because
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  prints Fortran diagnostics there that report on its internals rather
  than on the fit. The sink covers the whole solve rather than the
  solver alone, so a [`print()`](https://rdrr.io/r/base/print.html) or a
  [`cat()`](https://rdrr.io/r/base/cat.html) in the estimating function
  is discarded along with them and its output does not reach the
  console. Warnings and messages are untouched, since neither goes to
  standard output, so [`message()`](https://rdrr.io/r/base/message.html)
  traces an estimating function under every solver.

- maxiter:

  Integer maximum iterations for the solver (default 5000). Must be a
  single positive whole number, which is checked whichever solver is in
  force.

- tolerance:

  Numeric tolerance for the solver (default 1e-9).

- deriv_method:

  Character string for the derivative method used to compute the bread
  matrix. One of `"capprox"` (central difference), `"fapprox"` (forward
  difference), `"bapprox"` (backward difference), or `"exact"`
  (forward-mode automatic differentiation). Default `"capprox"`. Exact
  differentiation removes the finite-difference step size but requires
  that the estimating equation is composed of operations the autodiff
  supports (see
  [`auto_differentiation()`](https://r-causal.github.io/deli/reference/auto_differentiation.md)).
  For post-estimation transforms and the survival prediction helpers,
  exact differentiation is also available through
  [`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md).

- dx:

  Numeric step size for numerical differentiation; ignored when
  `deriv_method = "exact"` (default 1e-9). Must be a single positive
  finite number, which is checked whichever `deriv_method` is in force.
  The step is absolute and is floored at the floating-point resolution
  of each estimate, so a large parameter magnitude cannot silently
  reduce it to nothing; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

- allow_pinv:

  Logical. Use pseudo-inverse if bread is singular? Default `TRUE`.

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

A modified `MEstimator` object with populated `theta`, `bread`, `meat`,
`variance`, and `asymptotic_variance` properties.

## Details

The estimates and every matrix built from them carry parameter names,
taken from the first of two channels that supplies them. Names on `init`
come first: a parameter the caller named keeps that name, and one left
unnamed is numbered by position, so `init = c(mu = 0, 1)` gives
`c("mu", "theta_2")`. Where `init` carries no names at all, the row
names of the estimating functions are read instead. That is how a
`stacked_equations` function names the parameters it defines, since the
estimator otherwise sees only an opaque closure returning a matrix. The
formula interface names `init` from the model matrix columns where it
can account for its length, one parameter per design column or one more
for the parameter the equation appends, so the row names carry the
function interface and every formula fit of another length.
[`ee_gformula()`](https://r-causal.github.io/deli/reference/ee_gformula.md)
is one of those: its extra parameter leads rather than trails, so the
formula interface leaves its `init` unnamed and the fit takes the labels
the equation writes on its own rows.

Row names are read only when they label every parameter distinctly: one
name per parameter, none empty, none missing, and no two alike. An
incomplete or repetitive set is discarded rather than patched up,
because it is usually an accident of how the stack was built.
[`rbind()`](https://rdrr.io/r/base/cbind.html) pads the rows of an
unlabeled block with empty strings, and `t(X * resid)` on a design whose
intercept column has no name produces the same shape;
[`rbind()`](https://rdrr.io/r/base/cbind.html) also names each row after
the variable that supplied it, so two blocks that each begin with a
variable of the same name repeat a label. The count has to match as
well, which is what keeps an over-identified `GMMEstimator` out: its
rows are moment conditions and outnumber the parameters, so their labels
describe the equations. Where neither channel applies, the parameters
are numbered `theta_1` through `theta_p`.

Row names survive exact differentiation. Assigning them is ignored while
a value carries derivatives, since the labels are read from the plain
evaluation at the solved values. `rownames<-` hands its work to
`dimnames<-`, which is a generic, so deli registers the setter there and
an estimating function can label its rows from anywhere.

Many of the built-in estimating equations name their own rows, so a fit
that passes one of them a wholly unnamed `init` comes back labeled
rather than numbered. Each documents its labels under **Value**. The
rule above still governs a stack built from them: two blocks that name
the same parameter, as two
[`ee_ipw()`](https://r-causal.github.io/deli/reference/ee_ipw.md) blocks
do, repeat a label and the fit is numbered instead, and one named block
stacked with an unnamed one is incomplete and is numbered as well. Name
`init` where a stack needs labels the blocks cannot agree on.

## Examples

``` r
psi <- function(theta) {
  y <- c(1, 2, 3, 4, 5)
  matrix(y - theta[1], nrow = 1)
}
m <- MEstimator(stacked_equations = psi, init = 0) |>
  estimate()
coef(m)
#> theta_1 
#>       3 
```
