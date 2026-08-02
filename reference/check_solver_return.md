# Check the return of a custom solver

Validates that a user-supplied solver returned a numeric vector of the
expected length, so a malformed return produces an informative error
rather than an opaque failure while assembling the sandwich components.

## Usage

``` r
check_solver_return(theta, n_params)
```

## Arguments

- theta:

  The value returned by the custom solver.

- n_params:

  The expected number of solved parameters.

## Value

Invisible `NULL`. Raises an error if the return is invalid.
