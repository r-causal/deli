# Tests for default_param_names() and the internal `%||%` operator

# ---- Helpers ----

# Two parameters, the mean and the variance, so a partially named `init` has
# one named and one unnamed entry to distinguish.
mean_variance_psi <- function(y) {
  function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }
}

mean_variance_y <- function() {
  c(1.2, 2.4, 3.1, 4.8, 5.5, 6.9, 7.2, 8.1)
}

# A summary object built directly, so the printing path can be reached with a
# parameter vector carrying whatever names the test needs. A fitted estimator
# always names its parameters, so this is the only way to exercise the fallback
# in the printing code.
summary_with_theta <- function(theta) {
  p <- length(theta)
  EstimatorSummary(
    theta = theta,
    se = rep(1, p),
    ci = cbind(lower = theta - 1, upper = theta + 1),
    z = theta,
    p = rep(0.5, p),
    s = rep(1, p),
    n_obs = 10L,
    n_params = as.integer(p),
    alpha = 0.05,
    estimator = "MEstimator"
  )
}

# default_param_names ----------------------------------------------------------

test_that("default_param_names() numbers the parameters when no names are given", {
  expect_equal(default_param_names(NULL, 3), c("theta_1", "theta_2", "theta_3"))
})

test_that("default_param_names() returns a complete name vector unchanged", {
  expect_equal(default_param_names(c("a", "b", "c"), 3), c("a", "b", "c"))
})

test_that("default_param_names() fills the empty entries of a partial vector", {
  # `c(a = 1, 2, 3)` gives names `c("a", "", "")`, and an empty string is not a
  # usable label, so each empty entry takes the default for its position.
  expect_equal(
    default_param_names(c("a", "", "c"), 3),
    c("a", "theta_2", "c")
  )
})

test_that("default_param_names() fills missing entries of a partial vector", {
  expect_equal(
    default_param_names(c("a", NA, "c"), 3),
    c("a", "theta_2", "c")
  )
})

test_that("default_param_names() returns character(0) for zero parameters", {
  expect_equal(default_param_names(NULL, 0), character(0))
  expect_equal(default_param_names(character(0), 0), character(0))
})

test_that("default_param_names() always returns one name per parameter", {
  # A name vector shorter than the parameter count leaves the tail to the
  # defaults rather than returning a vector too short to label the estimates.
  expect_equal(
    default_param_names(c("a", "b"), 4),
    c("a", "b", "theta_3", "theta_4")
  )
  expect_length(default_param_names(c("a", "b", "c", "d"), 2), 2)
})

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

# Parameter names on fitted estimators -----------------------------------------

test_that("estimate() fills the unnamed entries of a partially named init", {
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = mean_variance_psi(y),
    init = c(mu = 0, 1)
  )

  expect_equal(names(coef(m)), c("mu", "theta_2"))
  expect_equal(dimnames(vcov(m)), list(c("mu", "theta_2"), c("mu", "theta_2")))
  expect_equal(rownames(confint(m)), c("mu", "theta_2"))
  expect_equal(tidy(m)$term, c("mu", "theta_2"))
})

test_that("gmm_estimate() fills the unnamed entries of a partially named init", {
  y <- mean_variance_y()
  g <- gmm_estimate(
    stacked_equations = mean_variance_psi(y),
    init = c(mu = 0, 1)
  )

  expect_equal(names(coef(g)), c("mu", "theta_2"))
  expect_equal(dimnames(vcov(g)), list(c("mu", "theta_2"), c("mu", "theta_2")))
})

test_that("estimate() numbers the parameters of a wholly unnamed init", {
  y <- mean_variance_y()
  m <- m_estimate(stacked_equations = mean_variance_psi(y), init = c(0, 1))

  expect_equal(names(coef(m)), c("theta_1", "theta_2"))
})

test_that("summary() labels a subset when the estimates have lost their names", {
  # A fitted estimator always names its parameters, so the fallback in
  # summarize_estimator() is defensive. Dropping the names by hand reaches it.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = mean_variance_psi(y),
    init = c(mu = 0, sigma2 = 1)
  )
  m@theta <- unname(m@theta)

  expect_equal(names(summary(m, subset = 2L)@theta), "theta_2")
})

test_that("summary() labels a subset when only some estimates are named", {
  # As above, defensive: nothing reaches summarize_estimator() with a partially
  # named theta on its own.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = mean_variance_psi(y),
    init = c(mu = 0, sigma2 = 1)
  )
  names(m@theta) <- c("mu", "")

  expect_equal(names(summary(m, subset = 2L)@theta), "theta_2")
  expect_equal(names(summary(m, subset = 1L)@theta), "mu")
})

test_that("printing a summary numbers rows when the estimates carry no names", {
  # cli writes to the message connection, so capture.output() sees nothing and
  # cli_fmt() is what collects the printed lines.
  output <- cli::cli_fmt(print(summary_with_theta(c(1, 2))))

  expect_true(any(grepl("^theta_1 ", output)))
  expect_true(any(grepl("^theta_2 ", output)))
})

test_that("printing a summary fills the unnamed rows of a partial name vector", {
  theta <- c(1, 2)
  names(theta) <- c("mu", "")

  output <- cli::cli_fmt(print(summary_with_theta(theta)))

  expect_true(any(grepl("^mu ", output)))
  expect_true(any(grepl("^theta_2 ", output)))
})
