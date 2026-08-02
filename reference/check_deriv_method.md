# Normalize and validate a derivative-method argument

Lowercases `deriv_method` and validates it against the supported
options, mirroring Python Delicatessen, which lowercases every method
comparison and accepts any case. Returns the normalized value so callers
can branch on it directly.

## Usage

``` r
check_deriv_method(deriv_method, call = rlang::caller_env())
```

## Arguments

- deriv_method:

  The derivative method supplied by the caller. One of `"capprox"`,
  `"fapprox"`, `"bapprox"`, or `"exact"`, in any case.

- call:

  The frame to report the refusal against. This function judges an
  argument the caller wrote and appears in no man page, so the default
  is the frame one up rather than this one. Every caller runs the check
  in its own body, so that frame is the entry point the caller reached,
  except in `delta_method_impl()`, which is itself a worker and passes
  the frame it was given.

## Value

The lower-cased method string. Raises an error if the value is not a
single supported string.
