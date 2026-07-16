# Tests for survival estimating equations (ee_aft, ee_plogit)

# ---- ee_aft shape tests ----

test_that("ee_aft returns correct shape for weibull", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  b <- ncol(X)

  # theta: beta (b elements) + log(1/sigma) (1 element)
  theta <- c(rep(0, b), 0)

  result <- ee_aft(theta, X = X, time = y, event = delta,
                   distribution = "weibull")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

test_that("ee_aft returns correct shape for exponential", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  b <- ncol(X)

  # Exponential has no sigma parameter

  theta <- rep(0, b)

  result <- ee_aft(theta, X = X, time = y, event = delta,
                   distribution = "exponential")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b)
  expect_equal(ncol(result), n)
})

test_that("ee_aft returns correct shape for log-logistic", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  b <- ncol(X)

  theta <- c(rep(0, b), 0)

  result <- ee_aft(theta, X = X, time = y, event = delta,
                   distribution = "log-logistic")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

test_that("ee_aft returns correct shape for log-normal", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  b <- ncol(X)

  theta <- c(rep(0, b), 0)

  result <- ee_aft(theta, X = X, time = y, event = delta,
                   distribution = "log-normal")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

# ---- ee_aft solve tests ----

test_that("ee_aft weibull solves via MEstimator", {
  set.seed(101)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  # Generate Weibull survival times
  shape <- 0.75
  scale_factor <- exp(2 + 0.5 * x)
  y_true <- scale_factor * rweibull(n, shape = shape)
  # Censoring
  censor <- rexp(n, rate = 0.1)
  y <- pmin(y_true, censor)
  delta <- as.numeric(y_true <= censor)
  X <- cbind(1, x)

  psi <- function(theta) {
    ee_aft(theta, X = X, time = y, event = delta, distribution = "weibull")
  }

  m <- MEstimator(stacked_equations = psi, init = c(2, 0, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_aft exponential solves via MEstimator", {
  set.seed(202)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  # Exponential survival times (Weibull with shape = 1)
  y_true <- exp(2 + 0.5 * x) * rexp(n)
  censor <- rexp(n, rate = 0.05)
  y <- pmin(y_true, censor)
  delta <- as.numeric(y_true <= censor)
  X <- cbind(1, x)

  psi <- function(theta) {
    ee_aft(theta, X = X, time = y, event = delta, distribution = "exponential")
  }

  m <- MEstimator(stacked_equations = psi, init = c(2, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_aft log-logistic solves via MEstimator", {
  set.seed(313)
  n <- 300
  x <- rbinom(n, 1, 0.5)
  # Log-logistic times: generate via inverse CDF
  u <- runif(n)
  scale_factor <- exp(2 + 0.5 * x)
  shape_param <- 2
  y_true <- scale_factor * (u / (1 - u))^(1 / shape_param)
  censor <- rexp(n, rate = 0.02)
  y <- pmin(y_true, censor)
  delta <- as.numeric(y_true <= censor)
  X <- cbind(1, x)

  psi <- function(theta) {
    ee_aft(theta, X = X, time = y, event = delta, distribution = "log-logistic")
  }

  m <- MEstimator(stacked_equations = psi, init = c(2, 0.5, 0.5))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_aft log-normal solves via MEstimator", {
  set.seed(404)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  # Log-normal survival times
  y_true <- exp(2 + 0.5 * x + rnorm(n))
  censor <- rexp(n, rate = 0.02)
  y <- pmin(y_true, censor)
  delta <- as.numeric(y_true <= censor)
  X <- cbind(1, x)

  psi <- function(theta) {
    ee_aft(theta, X = X, time = y, event = delta, distribution = "log-normal")
  }

  m <- MEstimator(stacked_equations = psi, init = c(2, 0, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
})

# ---- ee_aft error tests ----

test_that("ee_aft errors on unsupported distribution", {
  X <- cbind(1, rnorm(10))
  y <- rexp(10)
  delta <- rbinom(10, 1, 0.7)

  expect_error(
    ee_aft(c(0, 0, 0), X = X, time = y, event = delta,
           distribution = "gamma"),
    "not supported"
  )
})

test_that("ee_aft errors on invalid delta values", {
  X <- cbind(1, rnorm(10))
  y <- rexp(10)
  delta <- c(0, 1, 2, 0, 1, 0, 1, 0, 1, 0)

  expect_error(
    ee_aft(c(0, 0, 0), X = X, time = y, event = delta,
           distribution = "weibull"),
    "zero or one"
  )
})

# ---- ee_aft weights and offset tests ----

test_that("ee_aft supports weights", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  wts <- runif(n, 0.5, 1.5)

  result <- ee_aft(c(0, 0, 0), X = X, time = y, event = delta,
                   distribution = "weibull", weights = wts)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

test_that("ee_aft supports offset", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)
  off <- rep(0.1, n)

  result <- ee_aft(c(0, 0, 0), X = X, time = y, event = delta,
                   distribution = "weibull", offset = off)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

# ---- ee_plogit shape tests ----

test_that("ee_plogit returns correct shape with default S", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n))
  # Discrete survival times
  y <- sample(1:5, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  # Number of unique event times
  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  b <- ncol(X)
  theta <- rep(0, b + k)

  result <- ee_plogit(theta, X = X, time = y, event = delta)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + k)
  expect_equal(ncol(result), n)
})

test_that("ee_plogit returns correct shape with custom S", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n))
  y <- sample(1:5, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  # Custom S: intercept + linear time
  max_t <- max(y)
  t_steps <- seq_len(max_t)
  S_mat <- cbind(1, t_steps)

  b <- ncol(X)
  # theta: b covariate params + ncol(S_mat) time params
  theta <- rep(0, b + ncol(S_mat))

  result <- ee_plogit(theta, X = X, time = y, event = delta, S = S_mat)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + ncol(S_mat))
  expect_equal(ncol(result), n)
})

# ---- ee_plogit solve tests ----

test_that("ee_plogit solves with disjoint indicators via MEstimator", {
  set.seed(505)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  # Generate discrete survival times
  hazard <- inverse_logit(-3 + 0.5 * x)
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  max_time <- 10
  for (i in seq_len(n)) {
    for (k in seq_len(max_time)) {
      if (runif(1) < hazard[i]) {
        y[i] <- k
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  X <- cbind(x)
  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }

  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
})

test_that("ee_plogit solves with custom S matrix (Gompertz-like)", {
  set.seed(606)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  # Generate discrete survival times with time-varying hazard
  max_time <- 8
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  for (i in seq_len(n)) {
    for (k in seq_len(max_time)) {
      h <- inverse_logit(-4 + 0.3 * x[i] + 0.1 * k)
      if (runif(1) < h) {
        y[i] <- k
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  X <- cbind(x)
  # Gompertz-like S: intercept + linear time
  t_steps <- seq_len(max_time)
  S_mat <- cbind(1, t_steps)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta, S = S_mat)
  }

  inits <- c(0, -4, 0.1)
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
})

# ---- ee_plogit error tests ----

test_that("ee_plogit errors on S dimension mismatch", {
  set.seed(42)
  n <- 20
  X <- cbind(rnorm(n))
  y <- sample(1:5, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  # S with wrong number of rows
  S_bad <- cbind(1, 1:3)  # 3 rows but max(y) = 5

  expect_error(
    ee_plogit(c(0, 0, 0), X = X, time = y, event = delta, S = S_bad),
    "mismatch"
  )
})

# ---- ee_plogit weights test ----

test_that("ee_plogit supports weights", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n))
  y <- sample(1:5, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)
  wts <- runif(n, 0.5, 1.5)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  theta <- rep(0, ncol(X) + k)
  result <- ee_plogit(theta, X = X, time = y, event = delta, weights = wts)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), n)
})

# ---- ee_survival_model shape tests ----

test_that("ee_survival_model returns correct shape for weibull", {
  set.seed(42)
  n <- 50
  t_obs <- rweibull(n, shape = 1.5, scale = 2)
  delta <- rbinom(n, 1, 0.7)

  theta <- c(0.5, 1.5)
  result <- ee_survival_model(theta, time = t_obs, event = delta,
                              distribution = "weibull")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), n)
})

test_that("ee_survival_model returns correct shape for exponential", {
  set.seed(42)
  n <- 50
  t_obs <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)

  theta <- c(0.5)
  result <- ee_survival_model(theta, time = t_obs, event = delta,
                              distribution = "exponential")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), n)
})

test_that("ee_survival_model returns correct shape for gompertz", {
  set.seed(42)
  n <- 50
  t_obs <- rexp(n, rate = 0.5)
  delta <- rbinom(n, 1, 0.7)

  theta <- c(0.5, 0.1)
  result <- ee_survival_model(theta, time = t_obs, event = delta,
                              distribution = "gompertz")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), n)
})

# ---- ee_survival_model solve tests ----

test_that("ee_survival_model weibull solves via MEstimator", {
  set.seed(123)
  n <- 300
  true_lambda <- 0.5
  true_gamma <- 1.2
  # Generate Weibull times using the hazard parameterization:
  # h(t) = lambda * gamma * t^(gamma-1)
  # S(t) = exp(-lambda * t^gamma)
  # Use inverse CDF: t = (-log(U)/lambda)^(1/gamma)
  u <- runif(n)
  t_true <- (-log(u) / true_lambda)^(1 / true_gamma)
  # Censoring
  censor <- rexp(n, rate = 0.2)
  t_obs <- pmin(t_true, censor)
  delta <- as.numeric(t_true <= censor)

  psi <- function(theta) {
    ee_survival_model(theta, time = t_obs, event = delta,
                      distribution = "weibull")
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.3, 1.0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
  # Parameters should be in reasonable range of truth
  expect_true(abs(m@theta[1] - true_lambda) < 0.3)
  expect_true(abs(m@theta[2] - true_gamma) < 0.3)
})

test_that("ee_survival_model exponential solves via MEstimator", {
  set.seed(234)
  n <- 300
  true_lambda <- 0.8
  # Exponential: S(t) = exp(-lambda * t)
  t_true <- rexp(n, rate = true_lambda)
  censor <- rexp(n, rate = 0.3)
  t_obs <- pmin(t_true, censor)
  delta <- as.numeric(t_true <= censor)

  psi <- function(theta) {
    ee_survival_model(theta, time = t_obs, event = delta,
                      distribution = "exponential")
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.5))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
  expect_true(abs(m@theta[1] - true_lambda) < 0.4)
})

test_that("ee_survival_model gompertz solves via MEstimator", {
  set.seed(345)
  n <- 300
  true_lambda <- 0.3
  true_gamma <- 0.5
  # Gompertz: S(t) = exp(-(lambda/gamma)*(exp(gamma*t) - 1))
  # Inverse CDF: t = (1/gamma) * log(1 - (gamma/lambda)*log(U))
  u <- runif(n)
  t_true <- (1 / true_gamma) * log(1 - (true_gamma / true_lambda) * log(u))
  # Remove any invalid times (should be rare with these params)
  valid <- is.finite(t_true) & t_true > 0
  t_true <- t_true[valid]
  n_valid <- length(t_true)

  censor <- rexp(n_valid, rate = 0.3)
  t_obs <- pmin(t_true, censor)
  delta <- as.numeric(t_true <= censor)

  psi <- function(theta) {
    ee_survival_model(theta, time = t_obs, event = delta,
                      distribution = "gompertz")
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.2, 0.3))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) > 0))
})

# ---- ee_survival_model error tests ----

test_that("ee_survival_model errors on unsupported distribution", {
  t_obs <- rexp(10)
  delta <- rbinom(10, 1, 0.7)

  expect_error(
    ee_survival_model(c(0.5, 1), time = t_obs, event = delta,
                      distribution = "gamma"),
    "not supported"
  )
})

test_that("ee_survival_model errors on invalid delta values", {
  t_obs <- rexp(10)
  delta <- c(0, 1, 2, 0, 1, 0, 1, 0, 1, 0)

  expect_error(
    ee_survival_model(c(0.5, 1), time = t_obs, event = delta,
                      distribution = "weibull"),
    "zero or one"
  )
})

# ---- consistent time/event argument names (bd-9ejm.4 contract) --------------
#
# These tests enforce the naming contract that the survival functions share a
# single idiomatic pair of argument names, with NO deprecated aliases:
#
#   ee_aft(theta, X, time, event, distribution, weights, offset)
#   ee_survival_model(theta, time, event, distribution)
#   ee_plogit(theta, X, time, event, S, unique_times, weights, offset)
#   plogit_predict(theta, time, event, X, S, times_to_predict, measure,
#                  unique_times)
#
# `time` is always the observed (possibly censored) time and `event` is
# always the 1/0 event indicator. The prediction-grid arguments `times`
# (aft_predictions_*) and `times_to_predict` (plogit_predict) are a distinct
# concept and are left unchanged.
#
# The previous signatures used inconsistent names for the same two vectors,
# so the rejection tests below assert that those names no longer exist:
#
#   ee_aft(theta, X, y, delta, ...)          y = time,  delta = event
#   ee_survival_model(theta, t, delta, ...)  t = time,  delta = event
#   ee_plogit(theta, X, y, t, ...)           y = time,  t     = event
#   plogit_predict(theta, t, delta, X, ...)  t = time,  delta = event
#
# ee_plogit had reversed the meaning of `t` relative to plogit_predict and to
# Python Delicatessen, and no two functions agreed on the names.
#
# Each acceptance test below builds a dataset where the time and event
# vectors differ, so an implementation that renames the arguments but keeps
# the values wired to the wrong internal slot produces a different result and
# fails. The expected values are computed independently from the estimating
# equation math, not by re-calling the function under test.

test_that("ee_survival_model accepts time/event with correct semantics", {
  # Exponential: contribution is event / lambda - time (gamma fixed at 1)
  lambda <- 2
  time <- c(3, 5, 7)
  event <- c(1, 0, 1)
  expected <- matrix(event / lambda - time, nrow = 1)

  result <- ee_survival_model(lambda, time = time, event = event,
                              distribution = "exponential")
  expect_equal(result, expected)
})

test_that("ee_survival_model rejects the old t/delta argument names", {
  time <- c(3, 5, 7)
  event <- c(1, 0, 1)

  expect_error(
    ee_survival_model(2, t = time, delta = event, distribution = "exponential")
  )
})

test_that("ee_aft accepts time/event with correct semantics", {
  # Exponential AFT: sigma = 1, so Z_i = log(time) - X beta.
  # lambda_eps = event (1 - exp(Z)) - (1 - event) exp(Z)
  # score      = -(1 / sigma) lambda_eps X, weight 1
  time <- c(2, 4, 6, 8)
  event <- c(1, 0, 1, 0)
  X <- cbind(1, c(0.5, -0.5, 1, -1))
  beta <- c(0.3, -0.2)

  z <- as.numeric(log(time) - X %*% beta)
  lambda_eps <- event * (1 - exp(z)) - (1 - event) * exp(z)
  expected <- t(X * (-1 * lambda_eps))

  result <- ee_aft(beta, X = X, time = time, event = event,
                   distribution = "exponential")
  expect_equal(result, expected)
})

test_that("ee_aft rejects the old y/delta argument names", {
  time <- c(2, 4, 6, 8)
  event <- c(1, 0, 1, 0)
  X <- cbind(1, c(0.5, -0.5, 1, -1))

  expect_error(
    ee_aft(c(0.3, -0.2), X = X, y = time, delta = event,
           distribution = "exponential"),
    "unused argument"
  )
})

test_that("ee_plogit accepts time/event with correct semantics", {
  # Default disjoint indicators. The output is reproduced independently from
  # the risk-set construction so that swapping time and event fails.
  X <- matrix(c(0, 1, 0, 1), ncol = 1)
  time <- c(1, 2, 1, 3)
  event <- c(1, 1, 0, 1)
  theta <- c(0.2, -0.5, 0.1, 0.3)

  unique_times <- sort(unique(time[event == 1]))   # c(1, 2, 3)
  k <- length(unique_times)
  time_design <- diag(k)
  time_design[, 1] <- 1

  log_odds_w <- as.numeric(X %*% theta[1])
  log_odds_t <- as.numeric(time_design %*% theta[2:4])
  log_odds_w_matrix <- matrix(rep(log_odds_w, each = k), nrow = k)

  y_obs <- rep(event, each = k) * outer(unique_times, time, "==")
  y_pred <- inverse_logit(log_odds_w_matrix + log_odds_t)
  in_risk_set <- outer(unique_times, time, "<=")
  storage.mode(in_risk_set) <- "double"
  residual <- (y_obs - y_pred) * in_risk_set
  expected <- rbind(t(X * colSums(residual)), residual)

  result <- ee_plogit(theta, X = X, time = time, event = event)
  expect_equal(result, expected)
})

test_that("ee_plogit rejects the old y/t argument names", {
  X <- matrix(c(0, 1, 0, 1), ncol = 1)
  time <- c(1, 2, 1, 3)
  event <- c(1, 1, 0, 1)
  theta <- c(0.2, -0.5, 0.1, 0.3)

  # In the current signature `y` is the time and `t` is the event indicator;
  # after the rename neither name exists.
  expect_error(
    ee_plogit(theta, X = X, y = time, t = event)
  )
})

test_that("ee_plogit + MEstimator with time/event matches Python fixture", {
  # Cross-checks that the renamed arguments carry the correct semantics all
  # the way through a fit: swapped meanings would not converge to the Python
  # theta and variance.
  ref <- load_fixture("ee_plogit")
  X <- ref$X
  S <- ref$S

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = ref$t, event = ref$delta, S = S)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_plogit", tolerance = 1e-4)
})

# ---- ee_plogit time-varying (2D) weight matrix (bd-9ejm.5) -------------------
# ee_plogit accepts an n-by-K weight matrix, one column per time interval of S.
# generate_weights passes the matrix through (validating its row count) and
# ee_plogit validates the column count (K must equal the number of intervals)
# before transposing and multiplying it into the residual matrix. These tests
# verify that end-to-end path against Python and confirm the weights genuinely
# change the fit. Fixture provenance: generate_plogit_weight_matrix_fixture.py,
# breast cancer data with a Gompertz-like S (K = max(t) = 225) and a fixed-seed
# n-by-K weight matrix that varies across both rows and columns.
test_that("ee_plogit + MEstimator with a 2D weight matrix matches Python", {
  ref <- load_fixture("ee_plogit_weight_matrix")
  X <- ref$X
  S <- ref$S
  W <- ref$weights          # loads as an n-by-K matrix

  # Guard the fixture geometry so a silently reshaped matrix cannot pass.
  expect_equal(dim(W), c(ref$n, ref$k))

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = ref$t, event = ref$delta, S = S, weights = W)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  # Tolerance 1e-4 matches the other plogit fixture: the pooled-logistic solve
  # over 225 intervals does not reproduce Python to 1e-6.
  expect_python_match(m, "ee_plogit_weight_matrix", tolerance = 1e-4)
})

test_that("ee_plogit 2D weight matrix genuinely changes the fit", {
  # Discrimination guard: the matrix-weighted theta must differ from the
  # unit-weighted theta by far more than the 1e-4 fixture tolerance, otherwise
  # a fit that dropped or mishandled the weights could still "match". The
  # recorded per-parameter gaps are about 0.00296 (statin), 0.219 (time
  # intercept), and 0.00321 (time slope); the smallest, 0.00296, is ~30x the
  # tolerance. We require every component to move by at least 2e-3.
  ref <- load_fixture("ee_plogit_weight_matrix")
  X <- ref$X
  S <- ref$S
  W <- ref$weights

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = ref$t, event = ref$delta, S = S, weights = W)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  gap <- abs(unname(m@theta) - ref$theta_unweighted)
  expect_true(all(gap > 2e-3))
})
