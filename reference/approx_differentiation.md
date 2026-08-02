# Numerical differentiation via finite differences

Computes the Jacobian matrix of a vector-valued function using forward,
backward, or central difference approximation.

## Usage

``` r
approx_differentiation(func, theta, method = "capprox", dx = 1e-09)
```

## Arguments

- func:

  A function that takes a numeric vector and returns a numeric vector.

- theta:

  Numeric vector of parameter values at which to evaluate the
  derivative.

- method:

  Character string specifying the approximation method. One of
  `"capprox"` (central difference, default), `"fapprox"` (forward
  difference), or `"bapprox"` (backward difference).

- dx:

  Numeric step size for the finite difference (default `1e-9`). The step
  is absolute, floored at the floating-point resolution of each
  parameter.

## Value

A matrix where element `[i, j]` is the partial derivative of the `i`-th
output with respect to the `j`-th parameter. A Jacobian holding an entry
whose difference was lost against the magnitude of the values of `func`
is returned all the same, with a warning carrying the class
`deli_finite_difference_lost`. One warning covers the call however many
entries were lost, and its wording is the same at every call, so a
caller wrapping several calls in `without_repeated_warnings()` reports
it once for the operation.

## Details

`dx` is an absolute perturbation, matching the `epsilon` argument of the
Python `delicatessen` library, so a `dx` carried over from there means
the same thing here. An absolute step has a limit, though: the spacing
of the doubles surrounding `theta` grows with `theta`, so a fixed `dx`
spans fewer and fewer of them as a parameter grows. At the default `dx`
it spans about 17,600 of them at `|theta| = 450`, 69 at
`|theta| = 1.3e5`, and barely one by
`|theta| = dx / .Machine$double.eps`, about `4.5e6`. The perturbation is
rounded away with them, at first losing significant digits and past
roughly `1.7e7` leaving `theta + dx` equal to `theta`, at which point
every difference is zero and the Jacobian collapses.

Each parameter's step is therefore floored at that magnitude's
floating-point resolution: it always spans at least ten thousand
representable values, so the perturbation actually applied reproduces
the one intended to about four significant digits. The floor is
`1e4 * .Machine$double.eps * |theta|`, and
`.Machine$double.eps * |theta|` runs from one spacing to two as
`|theta|` climbs from one power of two to the next, so the floor engages
wherever `dx` would span fewer than ten to twenty thousand of them:
every `|theta|` above `dx / (1e4 * .Machine$double.eps)`, about `450` at
the default `dx`. Below that magnitude the step is `dx` exactly and
nothing about the result changes. Above it the step applied is the floor
rather than `dx`, and the derivative moves with it: at `|theta| = 500`,
where `dx` still spans 17,592 values and is resolved to five significant
digits, the step is `1.11e-9`.

Where the floor does engage, the derivative is accurate to a few parts
in ten thousand rather than to the roughly `1e-7` a representable step
reaches. `deriv_method = "exact"` takes no step at all and is unaffected
by parameter magnitude, so it is the better choice for a badly scaled
problem.

The floor rescues a step lost against a large `theta`. A difference is
lost a second way that no reading of the step can see, because nothing
about the step is wrong: where the step is applied exactly but the
values of `func` are large, the change it produces falls below the
spacing of the doubles holding them, both evaluations round to the same
double, and the quotient is exactly zero. That is indistinguishable, in
the returned value, from a function that is genuinely flat. What
separates the two is the significance of the difference against the
magnitude of the values it was taken between, which is what
`deli_finite_difference_lost` reports; see **Value**.
