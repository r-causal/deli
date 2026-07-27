# ---- logit() ----------------------------------------------------------------

test_that("logit() matches qlogis() for various probabilities", {
  probs <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  for (p in probs) {
    expect_equal(logit(p), qlogis(p))
  }
})

test_that("logit() works with vector input", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(logit(probs), qlogis(probs))
})

test_that("logit() handles edge cases", {
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  expect_equal(logit(0.5), 0)
  # The documented correspondence with qlogis() has to hold at the boundaries
  # too, where the log-odds are infinite.
  expect_equal(logit(c(0, 0.5, 1)), qlogis(c(0, 0.5, 1)))
})

# ---- logit() autodiff compatibility -----------------------------
#
# logit() was qlogis(prob), which cannot accept a tangent-carrying argument, so
# a transform using logit() under deriv_method = "exact" errored with
# "Non-numeric argument to mathematical function" (inverse_logit already worked).
# Writing logit() as log(prob / (1 - prob)) lets the autodiff Math and Ops rules
# carry tangents, with identical numerics for numeric input. The derivative of
# log(p / (1 - p)) is 1 / (p (1 - p)).

test_that("logit() differentiates under exact autodiff", {
  x0 <- 0.3
  f <- function(x) logit(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  analytic <- 1 / (x0 * (1 - x0))
  expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
})

# ---- inverse_logit() --------------------------------------------------------

test_that("inverse_logit() matches plogis() for various log-odds", {
  logodds <- c(-5, -2, -1, 0, 1, 2, 5)
  for (lo in logodds) {
    expect_equal(inverse_logit(lo), plogis(lo))
  }
})

test_that("inverse_logit() works with vector input", {
  logodds <- c(-5, -2, -1, 0, 1, 2, 5)
  expect_equal(inverse_logit(logodds), plogis(logodds))
})

test_that("inverse_logit(0) returns 0.5", {
  expect_equal(inverse_logit(0), 0.5)
})

# ---- roundtrip logit / inverse_logit ----------------------------------------

test_that("inverse_logit(logit(p)) returns p (roundtrip)", {
  probs <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  expect_equal(inverse_logit(logit(probs)), probs)
})

test_that("logit(inverse_logit(x)) returns x (roundtrip)", {
  logodds <- c(-5, -2, -1, 0, 1, 2, 5)
  expect_equal(logit(inverse_logit(logodds)), logodds)
})

# ---- identity_transform() ---------------------------------------------------

test_that("identity_transform() returns scalar unchanged", {
  expect_equal(identity_transform(42), 42)
  expect_equal(identity_transform(-3.14), -3.14)
  expect_equal(identity_transform(0), 0)
})

test_that("identity_transform() returns vector unchanged", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(identity_transform(x), x)
})

test_that("identity_transform() returns matrix unchanged", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  expect_equal(identity_transform(m), m)
})

test_that("identity_transform() preserves character input", {
  expect_equal(identity_transform("hello"), "hello")
})

test_that("identity_transform() matches base::identity()", {
  # The documented counterpart of identity_transform() is base::identity(), which
  # unlike most of the utility functions is usable everywhere, including under
  # exact-mode autodiff. The two must be interchangeable on every input type.
  expect_identical(identity_transform(42), identity(42))
  expect_identical(identity_transform(c(1, 2, 3)), identity(c(1, 2, 3)))
  expect_identical(
    identity_transform(matrix(1:6, nrow = 2, ncol = 3)),
    identity(matrix(1:6, nrow = 2, ncol = 3))
  )
  expect_identical(identity_transform("hello"), identity("hello"))
})
