# Check that a finite-difference step is a single positive number

Validates `dx`, the absolute perturbation the finite-difference methods
apply. Nothing downstream reports a step that cannot be taken, because
[`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md)
floors each parameter's step at that magnitude's floating-point
resolution and the floor absorbs the two values a caller is most likely
to mean something by. A `dx` of zero or a negative `dx` becomes the
floor at every parameter away from zero, silently substituting a step
the caller did not ask for, and leaves a division by zero at a parameter
of exactly zero, where the floor is zero as well. A longer vector is
recycled one element per parameter, so each parameter is differentiated
with a different step and no two rows of the Jacobian are comparable.

## Usage

``` r
check_dx(dx, call = rlang::caller_env())
```

## Arguments

- dx:

  The finite-difference step supplied by the caller.

- call:

  The frame to report the refusal against, on the same terms as
  [`check_deriv_method()`](https://r-causal.github.io/deli/reference/check_deriv_method.md)'s.

## Value

Invisible `NULL`. Raises an error if the step is not a single positive
finite number.

## Details

The step is validated wherever it is supplied, including under
`deriv_method = "exact"` and on a prediction that asks for no standard
error, neither of which takes a step at all. A value that cannot be a
step is worth reporting whether or not this particular call would have
used it, since the alternative is accepting it in silence and rejecting
it on the next call.
