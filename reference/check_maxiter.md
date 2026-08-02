# Check that an iteration budget is a single positive whole number

Validates `maxiter`, the number of iterations a solver is allowed.
Nothing downstream reports a budget that cannot be counted with, and
each solver mishandles a bad one in its own way. A vector budget reached
[`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html),
whose own guard tests it with an `if` and failed with
`the condition has length > 1` against an expression no caller wrote,
and where the budget survived to be reported the non-convergence message
pluralized on the length of the vector rather than on the budget. A
budget that is not a number failed as `'maxiter' must be numeric` from
the same place. A budget below one asks for a solve with no iterations
in it.

## Usage

``` r
check_maxiter(maxiter, call = rlang::caller_env())
```

## Arguments

- maxiter:

  The iteration budget supplied by the caller.

- call:

  The frame to report the refusal against, on the same terms as
  [`check_dx()`](https://r-causal.github.io/deli/reference/check_dx.md)'s.

## Value

Invisible `NULL`. Raises an error if the budget is not a single positive
whole number.

## Details

Judged where it is supplied, on the same terms as
[`check_dx()`](https://r-causal.github.io/deli/reference/check_dx.md),
so a budget that cannot be used is reported whether or not the solver
this call reaches would have looked at it.
