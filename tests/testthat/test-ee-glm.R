# Tests for ee_glm() estimating equations using Python fixtures
# ee_glm() is a general GLM estimating equation that takes theta, X, y,
# distribution, and link arguments. For canonical links (identity/normal,
# logit/binomial, log/poisson), the deriv/variance ratio simplifies to 1,
# so the results should match ee_regression.

test_that("ee_glm normal/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with normal distribution and identity link is equivalent to linear
  # regression. The estimating equation simplifies to:
  # psi(theta) = (y - X %*% theta) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "normal", link = "identity")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm gaussian/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # "gaussian" is an alias for "normal" in the Python implementation
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gaussian", link = "identity")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm binomial/logit matches ee_regression logistic", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with binomial distribution and logit link is equivalent to logistic

  # regression. For this canonical link, deriv/variance = 1, so the
  # estimating equation simplifies to: psi(theta) = (y - expit(X %*% theta)) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "logit")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm binomial/logistic link alias matches ee_regression logistic", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # "logistic" is an alias for "logit" link in the Python implementation
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "logistic")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm poisson/log matches ee_regression poisson", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with Poisson distribution and log link is the canonical Poisson
  # regression. For this canonical link, deriv/variance = 1, so the
  # estimating equation simplifies to: psi(theta) = (y - exp(X %*% theta)) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "poisson", link = "log")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_glm rejects a y whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_glm(
      init,
      X = X,
      y = ref$y[seq_len(nrow(X) / 2)],
      distribution = "normal",
      link = "identity"
    ),
    "same length as the data"
  )
})

test_that("ee_glm rejects an offset whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_glm(
      init,
      X = X,
      y = ref$y,
      distribution = "normal",
      link = "identity",
      offset = rep(0, nrow(X) / 2)
    ),
    "same length as the data"
  )
})

# The family names and the frame their refusals report -------------------------
#
# `distribution` and `link` are read in helpers that are not exported and appear
# in no man page, so an unsupported name was reported against `inverse_link()`
# or `distribution_variance()` rather than against the estimating equation the
# caller wrote. And neither helper is reached at all by a `distribution` that is
# `NULL` or longer than one: the branch in ee_glm() that partitions `theta` tests
# it first, so base R reported `argument is of length zero` or `the condition has
# length > 1` against an `if` the caller never wrote.

reported_call <- function(err) {
  paste(deparse(conditionCall(err)), collapse = " ")
}

glm_case <- function() {
  set.seed(3)
  list(X = cbind(1, stats::rnorm(20)), y = stats::rbinom(20, 1, 0.5))
}

test_that("an unsupported distribution and link name ee_glm()", {
  case <- glm_case()

  err <- expect_error(
    ee_glm(
      c(0, 0),
      X = case$X,
      y = case$y,
      distribution = "nope",
      link = "identity"
    ),
    "not supported"
  )
  expect_match(reported_call(err), "ee_glm(", fixed = TRUE)
  expect_false(grepl("distribution_variance", reported_call(err), fixed = TRUE))

  err <- expect_error(
    ee_glm(
      c(0, 0),
      X = case$X,
      y = case$y,
      distribution = "normal",
      link = "nope"
    ),
    "not supported"
  )
  expect_match(reported_call(err), "ee_glm(", fixed = TRUE)
  expect_false(grepl("inverse_link", reported_call(err), fixed = TRUE))
})

test_that("ee_glm refuses a distribution that is not a single string", {
  case <- glm_case()

  for (value in list(NULL, c("normal", "poisson"), 1, NA_character_)) {
    err <- expect_error(
      ee_glm(
        c(0, 0),
        X = case$X,
        y = case$y,
        distribution = value,
        link = "identity"
      ),
      "distribution"
    )
    expect_match(conditionMessage(err), "single string")
    expect_match(reported_call(err), "ee_glm(", fixed = TRUE)
  }
})

test_that("ee_glm refuses a link that is not a single string", {
  case <- glm_case()

  for (value in list(NULL, c("identity", "log"), 1, NA_character_)) {
    err <- expect_error(
      ee_glm(
        c(0, 0),
        X = case$X,
        y = case$y,
        distribution = "normal",
        link = value
      ),
      "link"
    )
    expect_match(conditionMessage(err), "single string")
    expect_match(reported_call(err), "ee_glm(", fixed = TRUE)
  }
})

test_that("ee_ipw_msm refuses a distribution that is not a single string", {
  set.seed(4)
  n <- 30
  W <- cbind(1, stats::rnorm(n))
  V <- cbind(rep(1, n))
  A <- stats::rbinom(n, 1, 0.5)
  y <- stats::rnorm(n)

  for (value in list(NULL, c("normal", "poisson"))) {
    err <- expect_error(
      ee_ipw_msm(
        c(0, 0, 0),
        y = y,
        A = A,
        W = W,
        V = V,
        distribution = value,
        link = "identity"
      ),
      "distribution"
    )
    expect_match(conditionMessage(err), "single string")
    expect_match(reported_call(err), "ee_ipw_msm(", fixed = TRUE)
  }
})
