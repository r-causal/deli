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
