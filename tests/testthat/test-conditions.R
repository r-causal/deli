# Tests for the warning de-duplication scope in R/conditions.R and for the
# operations that install it: both estimate() methods and compute_sandwich().

# ---- Helpers ----

# Collects every warning signaled while `expr` runs and muffles each one, so
# the count is of what the scope delivered rather than of what reached the
# console. Muffling short-circuits R's default warning handling, so the count is
# the same under `options(warn = 2)`, which the whole suite runs under.
collect_warnings <- function(expr) {
  caught <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  caught
}

warning_messages <- function(caught) {
  vapply(caught, conditionMessage, character(1))
}

warning_classes <- function(caught) {
  vapply(caught, function(w) class(w)[[1L]], character(1))
}

percentile_data <- function() {
  set.seed(1)
  stats::rnorm(50)
}

percentile_psi <- function(y) {
  function(theta) ee_percentile(theta, y = y, q = 0.5)
}

# A fitted mean, whose estimating function is quiet, so the only warnings in the
# delta_method() tests are the ones the transform raises.
mean_fit <- function() {
  y <- percentile_data()
  m_estimate(function(theta) matrix(y - theta[1], nrow = 1), init = 0)
}

warning_transform <- function(theta) {
  cli::cli_warn("the transform is not differentiable")
  exp(theta[1])
}

# A LASSO fit that raises two different warnings: the non-differentiability
# warning from the penalty, once per evaluation of the estimating function, and
# a single report that rootSolve exhausted its iteration budget. The budget is
# capped so the fit reaches that report quickly.
lasso_fit <- function() {
  set.seed(2)
  n <- 60
  x <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  y <- as.numeric(x %*% c(1, 2, 0) + stats::rnorm(n))
  m_estimate(
    function(theta) {
      ee_lasso_regression(theta, X = x, y = y, model = "linear", penalty = 100)
    },
    init = c(0, 0, 0),
    maxiter = 20
  )
}

# without_repeated_warnings -----------------------------------------------

test_that("without_repeated_warnings delivers a repeated warning once", {
  caught <- collect_warnings(without_repeated_warnings({
    cli::cli_warn("the same thing")
    cli::cli_warn("the same thing")
    cli::cli_warn("the same thing")
  }))

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "the same thing")
})

test_that("without_repeated_warnings delivers each distinct message", {
  caught <- collect_warnings(without_repeated_warnings({
    cli::cli_warn("the first thing")
    cli::cli_warn("the second thing")
    cli::cli_warn("the first thing")
  }))

  expect_length(caught, 2)
  expect_match(warning_messages(caught), "first|second")
})

test_that("without_repeated_warnings returns the value of the expression", {
  expect_equal(without_repeated_warnings(1 + 1), 2)
})

test_that("without_repeated_warnings starts a fresh seen-set on every call", {
  once <- function() without_repeated_warnings(cli::cli_warn("the same thing"))

  caught <- collect_warnings({
    once()
    once()
  })

  expect_length(caught, 2)
})

# The key is the condition class alongside the message, so two warnings that
# read alike but carry different classes are both delivered. A class is part of
# what a condition is, and a caller matching on one class must be able to see it
# even when another class reached the same words first. Keying on the message
# alone would silently drop the second.
test_that("without_repeated_warnings delivers one warning per class", {
  caught <- collect_warnings(without_repeated_warnings({
    cli::cli_warn("shared wording", class = "deli_test_warning_a")
    cli::cli_warn("shared wording", class = "deli_test_warning_b")
    cli::cli_warn("shared wording", class = "deli_test_warning_a")
  }))

  expect_length(caught, 2)
  expect_equal(
    warning_classes(caught),
    c("deli_test_warning_a", "deli_test_warning_b")
  )
})

# Nested scopes are governed by the outermost one: the inner scope sees a
# warning first, records it and lets it through, and the outer scope then muffles
# whatever it has already delivered. So one operation nested inside another
# still reports each distinct warning once in total, which is the same promise
# the outer operation makes on its own.
test_that("without_repeated_warnings lets the outermost scope govern", {
  caught <- collect_warnings(without_repeated_warnings({
    without_repeated_warnings(cli::cli_warn("the same thing"))
    cli::cli_warn("the same thing")
  }))

  expect_length(caught, 1)
})

# A condition can inherit from "warning" and still be signaled without the
# muffleWarning restart, so the handler muffles with rlang::cnd_muffle(), which
# returns `FALSE` when there is no restart and leaves the condition to
# propagate. invokeRestart("muffleWarning") would turn the second signal into an
# error instead. The absorbing handler below is established outside the scope so
# that it runs after it, calling handlers running most recently established
# first, and it keeps the signaled condition from reaching testthat's own
# warning handler.
test_that("without_repeated_warnings muffles a warning that has no restart", {
  cnd <- structure(
    class = c("deli_test_restartless", "warning", "condition"),
    list(message = "no restart here", call = NULL)
  )
  signal_twice <- function() {
    withRestarts(signalCondition(cnd), deli_test_absorb = function() NULL)
    withRestarts(signalCondition(cnd), deli_test_absorb = function() NULL)
    "value"
  }

  result <- withCallingHandlers(
    without_repeated_warnings(signal_twice()),
    deli_test_restartless = function(w) invokeRestart("deli_test_absorb")
  )

  expect_equal(result, "value")
})

# The key is the class and the message, not the call, although base R's own
# warnings() keys on the call. Two warnings that read alike but were raised from
# different frames therefore collapse to one, and the second is dropped rather
# than merged. That loss is deliberate; R/conditions.R records what including
# the call would cost instead, and the next test measures the part of it that
# runs through a fit.
test_that("without_repeated_warnings keys on the message and not the call", {
  f1 <- function() base::warning("same text")
  f2 <- function() base::warning("same text")

  caught <- collect_warnings(without_repeated_warnings({
    f1()
    f2()
  }))

  expect_length(caught, 1)
  expect_equal(deparse(conditionCall(caught[[1L]])), "f1()")
})

# The other half of that decision is about this package's own warnings: the call
# would separate none of them, because none of them carries one.
# cli::cli_warn() records a call only when one is passed through to
# rlang::warn(), which no site here does. The warnings collected below cover an
# estimating equation's, a penalty's, the solver's, and the bread's.
test_that("no warning this package raises carries a call", {
  y <- percentile_data()
  na_bread_psi <- function(theta) matrix(rep(NA_real_, 3) * theta[1], nrow = 1)

  caught <- collect_warnings({
    ee_percentile(0, y = y, q = 0.5)
    lasso_fit()
    compute_bread(na_bread_psi, theta = 1)
  })

  expect_length(caught, 4)
  expect_true(any(warning_classes(caught) == "deli_solver_not_converged"))
  expect_match(
    warning_messages(caught),
    "bread matrix contains NA",
    all = FALSE
  )
  expect_true(all(vapply(caught, function(w) is.null(conditionCall(w)), NA)))
})

# estimate() --------------------------------------------------------------

test_that("estimate() delivers a repeated estimating-equation warning once", {
  y <- percentile_data()

  caught <- collect_warnings(
    m_estimate(percentile_psi(y), init = stats::median(y))
  )

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "not differentiable")
})

test_that("estimate() de-duplicates under exact differentiation too", {
  y <- percentile_data()

  caught <- collect_warnings(
    m_estimate(
      percentile_psi(y),
      init = stats::median(y),
      deriv_method = "exact"
    )
  )

  expect_length(caught, 1)
})

test_that("estimate() delivers one warning per distinct message", {
  caught <- collect_warnings(lasso_fit())

  expect_length(caught, 2)
  expect_length(unique(warning_messages(caught)), 2)
  expect_match(
    warning_messages(caught),
    "not always\\s+differentiable|did not converge"
  )
})

# The two warnings the LASSO fit above raises differ in class as well as in
# wording, so that test cannot tell whether the message is part of the key. These
# two differ in wording alone, which is the common case: cli attaches the same
# classes to every warning raised without one of its own.
test_that("estimate() delivers two warnings that differ only in wording", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    cli::cli_warn("the first thing")
    cli::cli_warn("the second thing")
    matrix(y - theta[1], nrow = 1)
  }

  caught <- collect_warnings(m_estimate(psi, init = 0))

  expect_length(caught, 2)
  expect_length(unique(warning_classes(caught)), 1)
  expect_match(warning_messages(caught), "first|second")
})

# What ignoring the call in the key buys. A warning() raised in the body of an
# estimating function reports the call to the frame that raised it, and a fit
# evaluates its estimating function from a different frame at each stage:
# `psi(init)` at the initial values, `stacked_equations(full_theta)` after the
# solve, and `stacked_equations(input_theta)` inside the bread. Keying on the
# call as well would deliver this warning three times rather than once.
test_that("estimate() delivers a base warning raised in a psi body once", {
  y <- percentile_data()
  psi <- function(theta) {
    base::warning("same text")
    matrix(y - theta[1], nrow = 1)
  }

  caught <- collect_warnings(m_estimate(psi, init = 0))

  expect_length(caught, 1)
})

test_that("estimate() starts a fresh seen-set on every fit", {
  y <- percentile_data()
  psi <- percentile_psi(y)

  caught <- collect_warnings({
    m_estimate(psi, init = stats::median(y))
    m_estimate(psi, init = stats::median(y))
  })

  expect_length(caught, 2)
})

# The outer fit uses a custom solver because rootSolve::multiroot() cannot be
# called from inside another call to itself: its C code rejects the environment
# the inner call would evaluate in. That is a rootSolve limitation and has
# nothing to do with the scope under test, so the outer fit stays out of
# rootSolve and the inner one keeps the default solver.
test_that("a nested fit reports its estimating equation's warning once", {
  y <- c(1, 2, 3, 4, 5, 6, 7)
  inner <- percentile_psi(y)
  outer <- function(theta) {
    m_estimate(inner, init = stats::median(y))
    ee_percentile(theta, y = y, q = 0.5)
  }

  caught <- collect_warnings(
    m_estimate(
      outer,
      init = stats::median(y),
      solver = function(stacked_equations, init) init
    )
  )

  expect_length(caught, 1)
})

test_that("a GMM fit delivers a repeated warning once", {
  y <- percentile_data()

  caught <- collect_warnings(
    gmm_estimate(percentile_psi(y), init = stats::median(y))
  )

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "not differentiable")
})

# Calling handlers run before R converts a warning into an error, so the scope
# sees the first warning and lets it through to be converted exactly as it would
# be without the scope. The fit therefore stops at the first evaluation of the
# estimating function rather than running to completion, and the error carries
# the estimating equation's own wording.
test_that("under options(warn = 2) a fit errors on the first warning", {
  y <- percentile_data()
  calls <- 0L
  psi <- function(theta) {
    calls <<- calls + 1L
    ee_percentile(theta, y = y, q = 0.5)
  }

  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)

  expect_error(m_estimate(psi, init = stats::median(y)), "not differentiable")
  expect_equal(calls, 1L)
})

# compute_sandwich() ------------------------------------------------------

test_that("compute_sandwich() delivers a repeated warning once", {
  y <- percentile_data()

  caught <- collect_warnings(
    compute_sandwich(percentile_psi(y), theta = stats::median(y))
  )

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "not differentiable")
})

# The bread's own warning travels through the scope like any other. One call to
# compute_sandwich() builds one bread, so there is nothing here for the scope to
# collapse, and what this pins is that the warning is neither repeated nor
# swallowed by the scope that surrounds it.
test_that("compute_sandwich() delivers the NA bread warning once", {
  # Finite where the meat reads it and not finite where the difference quotient
  # does, so the bread is the only part of the sandwich holding an NA. An
  # estimating function that is already not finite at the point it is evaluated
  # at is a different failure, judged before any of this is built.
  psi <- function(theta) {
    matrix(rep(if (theta[1] == 1) 1 else NA_real_, 3), nrow = 1)
  }

  caught <- collect_warnings(compute_sandwich(psi, theta = 1))

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "bread matrix contains NA")
})

test_that("compute_sandwich() starts a fresh seen-set on every call", {
  y <- percentile_data()
  psi <- percentile_psi(y)

  caught <- collect_warnings({
    compute_sandwich(psi, theta = stats::median(y))
    compute_sandwich(psi, theta = stats::median(y))
  })

  expect_length(caught, 2)
})

# delta_method() ----------------------------------------------------------

test_that("delta_method() delivers a repeated transform warning once", {
  caught <- collect_warnings(
    delta_method(mean_fit(), transform = warning_transform)
  )

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "not differentiable")
})

test_that("delta_method() de-duplicates on a numeric vector too", {
  fit <- mean_fit()

  caught <- collect_warnings(
    delta_method(
      coef(fit),
      transform = warning_transform,
      covariance = vcov(fit)
    )
  )

  expect_length(caught, 1)
  expect_match(warning_messages(caught), "not differentiable")
})

test_that("delta_method() starts a fresh seen-set on every call", {
  fit <- mean_fit()

  caught <- collect_warnings({
    delta_method(fit, transform = warning_transform)
    delta_method(fit, transform = warning_transform)
  })

  expect_length(caught, 2)
})

test_that("delta_method() de-duplicates under exact differentiation too", {
  caught <- collect_warnings(
    delta_method(
      mean_fit(),
      transform = warning_transform,
      deriv_method = "exact"
    )
  )

  expect_length(caught, 1)
})

# Direct calls ------------------------------------------------------------

# Nothing de-duplicates outside an operation that installs the scope. A direct
# call to an estimating equation warns every time, which is the behavior Python
# delicatessen has and which tests/testthat/test-ee-percentile.R pins.
test_that("a direct ee_percentile() call warns on every call", {
  y <- percentile_data()

  caught <- collect_warnings({
    ee_percentile(0, y = y, q = 0.5)
    ee_percentile(0, y = y, q = 0.5)
    ee_percentile(0, y = y, q = 0.5)
  })

  expect_length(caught, 3)
})

test_that("a direct ee_positive_mean_deviation() call warns on every call", {
  y <- percentile_data()

  caught <- collect_warnings({
    ee_positive_mean_deviation(c(0, 0), y = y)
    ee_positive_mean_deviation(c(0, 0), y = y)
  })

  expect_length(caught, 2)
})
