# Tests for the internal helpers in R/utils.R

# `%||%` -----------------------------------------------------------------------

test_that("`%||%` returns its left side unless that side is NULL", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
  expect_null(NULL %||% NULL)
  # An empty vector is not NULL, so it passes through.
  expect_equal(character(0) %||% "fallback", character(0))
})

test_that("`%||%` is defined by this package rather than inherited from base", {
  # `base::%||%` was added in R 4.4.0, and DESCRIPTION declares R (>= 4.3), so
  # relying on the base operator would leave the package broken on R 4.3.x. The
  # binding must therefore live in the package namespace itself.
  expect_true(exists("%||%", envir = asNamespace("deli"), inherits = FALSE))
})
