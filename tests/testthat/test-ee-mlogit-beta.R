# Tests for ee_mlogit and ee_beta_regression

# ee_mlogit tests -------------------------------------------------------------

test_that("ee_mlogit returns correct shape for 3-category outcome", {
  n <- 20
  set.seed(99)
  X <- cbind(1, rnorm(n))
  y_cat <- sample(1:3, n, replace = TRUE)
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  b <- ncol(X)
  k <- ncol(y)
  theta <- rep(0, b * (k - 1))

  result <- ee_mlogit(theta, X = X, y = y)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b * (k - 1))
  expect_equal(ncol(result), n)
})

test_that("ee_mlogit solves via MEstimator", {
  set.seed(123)
  n <- 50
  W <- rbinom(n, 1, 0.5)
  probs <- cbind(
    0.5 - 0.2 * W,
    0.3 + 0.1 * W,
    0.2 + 0.1 * W
  )
  y_cat <- sapply(1:n, function(i) sample(1:3, 1, prob = probs[i, ]))
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  X <- cbind(1, W)

  psi <- function(theta) {
    ee_mlogit(theta, X = X, y = y)
  }
  m <- MEstimator(stacked_equations = psi, init = rep(0, 4))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_mlogit errors on parameter count mismatch", {
  X <- cbind(1, rnorm(10))
  y <- cbind(
    c(1, 0, 0, 1, 0, 0, 1, 0, 1, 0),
    c(0, 1, 0, 0, 1, 0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  )
  expect_error(ee_mlogit(rep(0, 3), X = X, y = y), "mismatch")
})

# The outcome of a multinomial fit is a matrix rather than a vector, so its
# coercion at the head of ee_mlogit() is its own case. These pin what every
# accepted form of that matrix produces.

mlogit_indicator_data <- function() {
  set.seed(4321)
  n <- 60
  w <- rbinom(n, 1, 0.5)
  probs <- cbind(0.5 - 0.2 * w, 0.3 + 0.1 * w, 0.2 + 0.1 * w)
  category <- vapply(
    seq_len(n),
    function(i) sample(1:3, 1, prob = probs[i, ]),
    integer(1)
  )
  list(
    X = cbind(1, w),
    y = cbind(
      as.integer(category == 1),
      as.integer(category == 2),
      as.integer(category == 3)
    ),
    category = category
  )
}

fit_mlogit <- function(X, y, ...) {
  m_estimate(
    stacked_equations = function(theta) ee_mlogit(theta, X = X, y = y),
    init = rep(0, 4),
    ...
  )
}

test_that("ee_mlogit fits an indicator matrix of counted outcomes", {
  d <- mlogit_indicator_data()
  m <- fit_mlogit(d$X, d$y)

  expect_equal(
    unname(coef(m)),
    c(0.944461608841, -0.839101093183, 0.251314428281, -0.502628856562),
    tolerance = 1e-8
  )
  expect_equal(unname(vcov(m)[1, 1]), 0.198412487438, tolerance = 1e-8)
  expect_equal(names(coef(m)), c("theta_1", "theta_2", "theta_3", "theta_4"))
})

test_that("ee_mlogit reads a data frame of indicators as the same matrix", {
  d <- mlogit_indicator_data()
  y_df <- as.data.frame(d$y)
  names(y_df) <- c("reference", "second", "third")

  expect_identical(
    ee_mlogit(rep(0, 4), X = d$X, y = y_df),
    ee_mlogit(rep(0, 4), X = d$X, y = d$y)
  )
  expect_equal(coef(fit_mlogit(d$X, y_df)), coef(fit_mlogit(d$X, d$y)))
})

test_that("ee_mlogit reads a logical indicator matrix as the same matrix", {
  d <- mlogit_indicator_data()
  y_logical <- cbind(d$category == 1, d$category == 2, d$category == 3)

  expect_identical(
    ee_mlogit(rep(0, 4), X = d$X, y = y_logical),
    ee_mlogit(rep(0, 4), X = d$X, y = d$y)
  )
  expect_equal(coef(fit_mlogit(d$X, y_logical)), coef(fit_mlogit(d$X, d$y)))
})

test_that("ee_mlogit does not read parameter names off the outcome labels", {
  d <- mlogit_indicator_data()
  y_labeled <- d$y
  dimnames(y_labeled) <- list(
    paste0("obs", seq_len(nrow(d$y))),
    c("reference", "second", "third")
  )

  score <- ee_mlogit(rep(0, 4), X = d$X, y = y_labeled)
  expect_null(dimnames(score))
  expect_identical(score, ee_mlogit(rep(0, 4), X = d$X, y = d$y))
  expect_equal(
    names(coef(fit_mlogit(d$X, y_labeled))),
    c("theta_1", "theta_2", "theta_3", "theta_4")
  )
})

test_that("ee_mlogit counts a vector outcome as a single category", {
  d <- mlogit_indicator_data()

  # A vector becomes a one-column matrix, which leaves no non-reference
  # category to estimate.
  expect_error(
    ee_mlogit(rep(0, 4), X = d$X, y = as.numeric(d$category)),
    "For 1 categories"
  )
})

test_that("ee_mlogit reaches the same fit under exact differentiation", {
  d <- mlogit_indicator_data()
  y_df <- as.data.frame(d$y)

  expect_equal(
    coef(fit_mlogit(d$X, d$y, deriv_method = "exact")),
    coef(fit_mlogit(d$X, d$y))
  )
  expect_equal(
    coef(fit_mlogit(d$X, y_df, deriv_method = "exact")),
    coef(fit_mlogit(d$X, d$y, deriv_method = "exact"))
  )
})

# The softmax exponentiates every linear predictor, and exp() overflows to Inf
# a little past 709. These pin the probabilities the equation reads on both
# sides of that point, and the values it reads well below it.

test_that("ee_mlogit stays finite where exp() of the linear predictor overflows", {
  X <- cbind(1, c(0, 1))
  # The first observation is in the reference category and the second is in the
  # last category.
  y <- cbind(c(1, 0), c(0, 0), c(0, 1))

  # The first observation holds every linear predictor at zero, so its three
  # categories are equally probable. The second gives the first non-reference
  # category all of the probability, leaving the reference and last categories
  # at zero.
  expected <- rbind(
    c(1 / 3, -1),
    c(0, -1),
    c(1 / 3, 1),
    c(0, 1)
  )

  eta_max <- c(709, 710, 1000)
  scores <- lapply(
    eta_max,
    function(e) ee_mlogit(c(0, e, 0, 0), X = X, y = y)
  )

  expect_true(all(is.finite(unlist(scores))))
  expect_equal(scores, rep(list(expected), length(eta_max)))
})

test_that("ee_mlogit reads the same scores where nothing overflows", {
  d <- mlogit_indicator_data()
  score <- ee_mlogit(c(0.8, -0.5, 0.3, 0.2), X = d$X, y = d$y)

  expect_equal(
    rowSums(score),
    c(
      4.75135497492361,
      3.72049382506904,
      -2.68470787602859,
      -1.22280201723337
    ),
    tolerance = 1e-12
  )
  expect_equal(
    score[, 1],
    c(0.295025327936899, 0, -0.513585466435153, 0),
    tolerance = 1e-12
  )
})

# Below the shift threshold the equation exponentiates the linear predictors as
# they stand, so the softmax denominator itself grows to whatever exp() reads.
# A linear predictor a little past 354 puts that denominator past 1.3e154,
# where the probabilities are ordinary values but the square of the denominator
# is not. Differentiating exactly through the divisions that read those
# probabilities has to reach the same derivative the finite differences do.

test_that("ee_mlogit differentiates exactly where the softmax denominator is huge", {
  X <- cbind(1, c(0, 1))
  # The first observation is in the reference category and the second is in the
  # last category.
  y <- cbind(c(1, 0), c(0, 0), c(0, 1))
  # The second observation carries a linear predictor of 500 in the first
  # non-reference category, which is under the threshold the softmax shift
  # applies at, log(.Machine$double.xmax) - log(k + 1).
  theta <- c(0, 500, 0, 0)
  psi <- function(th) ee_mlogit(th, X = X, y = y)

  expect_true(all(is.finite(psi(theta))))

  bread <- expect_no_warning(compute_bread(psi, theta, deriv_method = "exact"))
  expect_true(all(is.finite(bread)))
  expect_equal(
    bread,
    compute_bread(psi, theta, deriv_method = "capprox"),
    tolerance = 1e-6
  )

  sandwich <- compute_sandwich(psi, theta, deriv_method = "exact")
  expect_true(all(is.finite(sandwich)))
})

# ee_beta_regression tests ----------------------------------------------------

test_that("ee_beta_regression returns correct shape", {
  set.seed(42)
  n <- 30
  X <- cbind(1, rnorm(n))
  y <- runif(n, 0.01, 0.99)
  b <- ncol(X)

  theta <- rep(0, b + 1)
  result <- ee_beta_regression(theta, X = X, y = y)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

test_that("ee_beta_regression solves via MEstimator", {
  set.seed(42)
  n <- 50
  W <- rnorm(n)
  mu <- 1 / (1 + exp(-(0.5 + 0.3 * W)))
  y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
  y <- pmin(pmax(y, 0.001), 0.999)
  X <- cbind(1, W)

  psi <- function(theta) {
    ee_beta_regression(theta, X = X, y = y)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0, log(10)))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
  expect_true(exp(m@theta[3]) > 0)
})

test_that("ee_beta_regression intercept-only model recovers mean", {
  set.seed(123)
  n <- 100
  mu_true <- 0.6
  phi_true <- 20
  y <- rbeta(n, shape1 = mu_true * phi_true, shape2 = (1 - mu_true) * phi_true)
  y <- pmin(pmax(y, 0.001), 0.999)
  X <- cbind(rep(1, n))

  psi <- function(theta) {
    ee_beta_regression(theta, X = X, y = y)
  }
  # Use logit(mean(y)) as init for better convergence
  init_beta <- log(mean(y) / (1 - mean(y)))
  m <- MEstimator(stacked_equations = psi, init = c(init_beta, log(phi_true)))
  m <- estimate(m, solver = "nleqslv")

  estimated_mean <- 1 / (1 + exp(-m@theta[1]))
  expect_equal(unname(estimated_mean), mu_true, tolerance = 0.1)
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_mlogit rejects a y whose rows differ from the rows of X", {
  n <- 20
  set.seed(99)
  X <- cbind(1, rnorm(n))
  y_cat <- sample(1:3, n, replace = TRUE)
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  theta <- rep(0, ncol(X) * (ncol(y) - 1))

  expect_error(
    ee_mlogit(theta, X = X, y = y[seq_len(n / 2), , drop = FALSE]),
    "same number of rows|same length as the data"
  )
})

test_that("ee_beta_regression rejects a y whose length differs from the rows of X", {
  set.seed(11)
  n <- 40
  X <- cbind(1, rnorm(n))
  y <- runif(n, 0.1, 0.9)
  theta <- c(0, 0, 0)

  expect_error(
    ee_beta_regression(theta, X = X, y = y[seq_len(n / 2)]),
    "same length as the data"
  )
})
