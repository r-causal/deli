# Extended tests for ee_glm() estimating equations
# These tests cover non-canonical link/distribution combinations and verify
# equivalence with ee_regression for canonical links using fixture data.

test_that("ee_glm probit/binomial converges", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Probit link with binomial distribution is a non-canonical GLM.
  # We verify that estimation converges and produces reasonable results.
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "probit")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Should converge: theta should be finite

  expect_true(all(is.finite(m@theta)))
  # Variance should be positive on diagonal
  expect_true(all(diag(m@variance) > 0))
  # Theta length should match init
  expect_length(m@theta, length(init))
})

test_that("ee_glm cloglog/binomial converges", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Complementary log-log link with binomial distribution is a non-canonical
  # GLM. We verify that estimation converges and produces reasonable results.
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "cloglog")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Should converge: theta should be finite
  expect_true(all(is.finite(m@theta)))
  # Variance should be positive on diagonal
  expect_true(all(diag(m@variance) > 0))
  # Theta length should match init
  expect_length(m@theta, length(init))
})

test_that("ee_glm poisson/log matches ee_regression poisson", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Poisson/log is a canonical link, so ee_glm should give the same results
  # as ee_regression with model = "poisson"
  psi_glm <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "poisson", link = "log")
  }

  psi_reg <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "poisson")
  }

  m_glm <- MEstimator(stacked_equations = psi_glm, init = init)
  m_glm <- estimate(m_glm)

  m_reg <- MEstimator(stacked_equations = psi_reg, init = init)
  m_reg <- estimate(m_reg)

  expect_equal(m_glm@theta, m_reg@theta, tolerance = 1e-6)
  expect_equal(diag(m_glm@variance), diag(m_reg@variance), tolerance = 1e-6)
})

test_that("ee_glm normal/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Normal/identity is a canonical link, so ee_glm should give the same
  # results as ee_regression with model = "linear"
  psi_glm <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "normal", link = "identity")
  }

  psi_reg <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }

  m_glm <- MEstimator(stacked_equations = psi_glm, init = init)
  m_glm <- estimate(m_glm)

  m_reg <- MEstimator(stacked_equations = psi_reg, init = init)
  m_reg <- estimate(m_reg)

  expect_equal(m_glm@theta, m_reg@theta, tolerance = 1e-6)
  expect_equal(diag(m_glm@variance), diag(m_reg@variance), tolerance = 1e-6)
})

# ---- the loglog link ---------------------------------------------------------
#
# The inverse loglog link is exp(-exp(-eta)), whose derivative with respect to
# eta is exp(-eta - exp(-eta)), a positive quantity: the inverse link is
# increasing. Python delicatessen returns the negative of it and deli ported the
# sign, so these tests separate the two claims that sign touches.
#
# The derivative is reported as a derivative, so its sign is part of the
# contract. Nothing deli returns depends on it, though, which is why fixing it
# changes no result. M-estimation is equivariant to scaling every row of psi by
# one nonzero constant: the root does not move, the bread changes sign, the meat
# is a cross-product and does not, and the two sign changes cancel in
# B^-1 M B^-T. predict() squares the derivative on its way to the response-scale
# standard error, so it does not see the sign either. The tests below the first
# hold whichever sign the link carries, and they are what says so.

# A binomial outcome with a loglog mean, small enough to fit in well under a
# second and seeded so the estimates are fixed.
loglog_data <- function() {
  set.seed(4)
  n <- 300
  d <- data.frame(
    x1 = round(stats::rnorm(n), 4),
    x2 = stats::rbinom(n, 1, 0.5)
  )
  eta <- 0.3 + 0.6 * d$x1 - 0.4 * d$x2
  d$y <- stats::rbinom(n, 1, exp(-exp(-eta)))
  d
}

# The loglog link as a glm() family component, which base R does not supply. It
# is the reference for the estimates rather than a restatement of deli's own
# arithmetic: glm() reaches them by iteratively reweighted least squares where
# deli solves the estimating equations.
loglog_link <- function() {
  structure(
    list(
      linkfun = function(mu) -log(-log(mu)),
      linkinv = function(eta) exp(-exp(-eta)),
      mu.eta = function(eta) exp(-eta - exp(-eta)),
      valideta = function(eta) TRUE,
      name = "loglog"
    ),
    class = "link-glm"
  )
}

test_that("inverse_link() returns the signed derivative of its own inverse link", {
  # Every link is checked here rather than the loglog alone, so the contract is
  # stated for all of them. test-predict.R states the same contract from the side
  # of the caller that reads the derivative, over the same roster. Zero is left
  # out of the grid because the "inverse" link is undefined there.
  eta <- c(-1.7, -0.4, 0.6, 1.3)
  for (link in inverse_link_names()) {
    expect_equal(
      inverse_link(eta, link)$dmu,
      numeric_inverse_link_deriv(eta, link),
      tolerance = 1e-6,
      info = link
    )
  }
})

test_that("ee_glm loglog/binomial matches glm() with a loglog link", {
  d <- loglog_data()
  fit <- m_estimate(
    y ~ x1 + x2,
    data = d,
    .ee = ee_glm,
    distribution = "binomial",
    link = "loglog"
  )
  ref <- stats::glm(
    y ~ x1 + x2,
    data = d,
    family = stats::binomial(link = loglog_link())
  )
  expect_equal(coef(fit), coef(ref), tolerance = 1e-5)
  expect_true(all(diag(vcov(fit)) > 0))
})

test_that("the loglog fit is unchanged by the sign of the link derivative", {
  # The same score written out with the derivative negated, which is the sign the
  # Python source this package was ported from carries. Every row of the two
  # estimating functions then differs by the one constant factor -1, so the
  # equivariance of M-estimation makes the estimates and the sandwich variance
  # identical rather than merely close. Writing the score with the sign
  # ee_glm() now uses would compare the fit against itself and demonstrate
  # nothing.
  d <- loglog_data()
  X <- cbind(1, d$x1, d$x2)
  y <- d$y
  negated_score <- function(theta) {
    eta <- as.numeric(X %*% theta)
    mu <- exp(-exp(-eta))
    dmu <- -exp(-eta - exp(-eta))
    t(X * ((y - mu) * dmu / (mu * (1 - mu))))
  }
  by_hand <- m_estimate(stacked_equations = negated_score, init = c(0, 0, 0))
  fit <- m_estimate(
    y ~ x1 + x2,
    data = d,
    .ee = ee_glm,
    distribution = "binomial",
    link = "loglog"
  )

  expect_equal(unname(coef(fit)), unname(coef(by_hand)), tolerance = 1e-8)
  expect_equal(unname(vcov(fit)), unname(vcov(by_hand)), tolerance = 1e-8)
})

test_that("the loglog fit differentiates exactly as it does by finite differences", {
  d <- loglog_data()
  args <- list(
    y ~ x1 + x2,
    data = d,
    .ee = ee_glm,
    distribution = "binomial",
    link = "loglog"
  )
  approx <- do.call(m_estimate, args)
  exact <- do.call(m_estimate, c(args, list(deriv_method = "exact")))
  expect_equal(coef(exact), coef(approx), tolerance = 1e-8)
  expect_equal(vcov(exact), vcov(approx), tolerance = 1e-6)
})

test_that("predict() on a loglog fit reports the mean and its scaled error", {
  d <- loglog_data()
  fit <- m_estimate(
    y ~ x1 + x2,
    data = d,
    .ee = ee_glm,
    distribution = "binomial",
    link = "loglog"
  )
  eta <- predict(fit)
  response <- predict(fit, type = "response", se.fit = TRUE)

  # The response scale is the inverse link of the linear predictor, and the
  # standard error on it is the magnitude of the derivative times the standard
  # error of the linear predictor, which is where the sign stops mattering.
  expect_equal(response$fit, exp(-exp(-eta)))
  expect_equal(
    response$se.fit,
    exp(-eta - exp(-eta)) * predict(fit, se.fit = TRUE)$se.fit
  )

  ref <- stats::glm(
    y ~ x1 + x2,
    data = d,
    family = stats::binomial(link = loglog_link())
  )
  expect_equal(
    unname(response$fit),
    unname(stats::fitted(ref)),
    tolerance = 1e-5
  )
})
