# Tests for regression estimating equations pipeline using Python fixtures

test_that("ee_regression linear matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_linear")

  # fromJSON with simplifyVector produces a matrix for X directly
  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Linear regression estimating equation: psi(theta) = X_i * (y_i - X_i^T theta)
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - Xb
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_regression logistic matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Inverse logit (expit) link function
  inverse_logit <- function(x) {
    1 / (1 + exp(-x))
  }

  # Logistic regression estimating equation
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - inverse_logit(Xb)
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_regression poisson matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Poisson regression estimating equation (exp link)
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - exp(Xb)
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_ridge_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_ridge_regression")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty
  n <- length(y)

  # Ridge regression: linear regression + L2 penalty
  # Python: penalty_terms = (penalty / n) * theta (for gamma=2, center=0)
  # Return: ((y - pred_y) * X).T - penalty_terms[:, None]
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - Xb

    # Regression score: p-by-n matrix
    score <- t(X * as.numeric(residuals))

    # L2 penalty term: (penalty / n) * theta, broadcast across all n columns
    penalty_terms <- (penalty / n) * theta
    score - penalty_terms # recycled across columns
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

# Behaviour-preservation pins for the performance refactors -------------------
#
# These reconstruct the estimating-function matrix element by element from the
# mathematical definition, independent of the internal transpose and coercion
# structure, so a refactor that changes how the p-by-n matrix is assembled or
# how design and outcome arguments are coerced is caught to machine precision.

manual_ee_regression <- function(
  theta,
  X,
  y,
  model,
  weights = NULL,
  offset = NULL
) {
  n <- nrow(X)
  p <- ncol(X)
  eta <- as.numeric(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + offset
  }
  pred <- switch(
    model,
    linear = eta,
    logistic = 1 / (1 + exp(-eta)),
    poisson = exp(eta)
  )
  w <- if (is.null(weights)) rep(1, n) else weights
  resid <- y - pred
  out <- matrix(0, nrow = p, ncol = n)
  for (i in seq_len(n)) {
    out[, i] <- X[i, ] * w[i] * resid[i]
  }
  out
}

test_that("ee_regression matrix matches the element-wise definition with weights and offset", {
  set.seed(11)
  n <- 60
  X <- cbind(1, rnorm(n), rnorm(n))
  theta <- c(0.3, -0.5, 0.8)
  weights <- runif(n, 0.25, 3)
  offset <- rnorm(n, sd = 0.2)

  for (model in c("linear", "logistic", "poisson")) {
    y <- switch(
      model,
      linear = as.numeric(X %*% theta) + rnorm(n),
      logistic = rbinom(n, 1, 0.5),
      poisson = rpois(n, 2)
    )
    expect_equal(
      ee_regression(
        theta,
        X = X,
        y = y,
        model = model,
        weights = weights,
        offset = offset
      ),
      manual_ee_regression(
        theta,
        X = X,
        y = y,
        model = model,
        weights = weights,
        offset = offset
      ),
      tolerance = 1e-12
    )
  }
})

test_that("ee_regression coerces a data frame design and an integer outcome without changing results", {
  set.seed(12)
  n <- 40
  Xm <- cbind(1, rnorm(n), rnorm(n))
  theta <- c(0.2, 0.4, -0.1)
  y_int <- as.integer(rpois(n, 2))

  ref <- ee_regression(theta, X = Xm, y = as.numeric(y_int), model = "poisson")

  # A data-frame design carries column names through as.matrix that a bare
  # matrix does not, so compare the numeric contributions after unname().
  Xdf <- as.data.frame(Xm)
  expect_equal(
    unname(ee_regression(theta, X = Xdf, y = y_int, model = "poisson")),
    unname(ref),
    tolerance = 1e-12
  )
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_regression rejects a y whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_regression(
      init,
      X = X,
      y = ref$y[seq_len(nrow(X) / 2)],
      model = "linear"
    ),
    "same length as the data"
  )
})

test_that("ee_regression rejects an offset whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_regression(
      init,
      X = X,
      y = ref$y,
      model = "linear",
      offset = rep(0, nrow(X) / 2)
    ),
    "same length as the data"
  )
})

test_that("ee_regression returns an unnamed score for both matrix and data-frame X", {
  set.seed(1)
  n <- 20
  df <- data.frame(int = 1, x1 = rnorm(n), x2 = rnorm(n))
  y <- rnorm(n)
  theta <- c(0.1, 0.2, 0.3)

  out_df <- ee_regression(theta, X = df, y = y, model = "linear")
  out_mat <- ee_regression(
    theta,
    X = unname(as.matrix(df)),
    y = y,
    model = "linear"
  )

  expect_null(dimnames(out_df))
  expect_identical(out_df, out_mat)
})
