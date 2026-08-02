# Contributing to deli

Thanks for your interest in deli. Bug reports and pull requests are
welcome at <https://github.com/r-causal/deli>. For bugs, please include
a minimal reproducible example.

## Development workflow

deli uses roxygen2 for documentation, testthat (edition 3) for tests,
and [air](https://posit-dev.github.io/air/) for formatting. After
changing any roxygen block, regenerate the documentation:

``` r

devtools::document()
```

Before opening a pull request, format the R sources and run the test
suite:

``` sh
air format .
```

``` r

devtools::test()
```

## Running R CMD check locally

``` r

rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))
```

`rcmdcheck()` defaults `libpath` to
[`.libPaths()`](https://rdrr.io/r/base/libPaths.html), so the packages
listed in `Suggests` are already visible to the check process and the
tests guarded by `skip_if_not_installed()` run. Passing
`libpath = .libPaths()` explicitly is the same thing and changes
nothing. If those guards do start skipping, the cause is a missing
installation rather than the check’s library path, so install the
`Suggests` and run again.

The check still reports fewer tests than `devtools::test()`, and the
difference is `skip_on_cran()`. `devtools::test()` sets `NOT_CRAN=true`,
while `R CMD check` leaves it unset, so anything gated on CRAN skips.
That covers the snapshot tests, which testthat skips on CRAN by default,
along with the two documentation tests that call `skip_on_cran()`
directly. The gap is expected, not a coverage hole. To exercise those
tests, run the suite with `NOT_CRAN` set, which is what
`devtools::test()` does.

The documentation tests carry a second guard. They read `_pkgdown.yml`,
which is listed in `.Rbuildignore` and so is absent from a built
tarball. They skip when the file is missing, by design, and run against
the source tree.

Continuous integration installs the full dependency set, including
`Suggests`, through `r-lib/actions/setup-r-dependencies`, so every job
in the matrix resolves these guards the same way a local check does.
