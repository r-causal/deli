# Exact-mode (forward-mode autodiff) sandwich variance for the regression
# estimating equations. Under deriv_method = "exact", theta enters an
# estimating function as a tangent-carrying object, and the linear predictor
# X %*% theta must flatten to a vector without dropping the tangent. These tests
# pin two invariants: the exact-mode bread, variance, and intervals match the
# Python Delicatessen reference computed with deriv_method='exact', and they
# agree with the finite-difference derivative on these smooth problems.
#
# Tolerance note: the internal exact-versus-capprox checks use 1e-6 because a
# central difference with dx = 1e-9 carries a rounding floor near 1e-7, so exact
# and capprox derivatives cannot be expected to agree more tightly than that.
# The fixture comparisons also use 1e-6: the reference theta comes from Python's
# Levenberg-Marquardt solver while the R fit uses rootSolve, and the small
# difference in the solved root propagates into the bread at that scale.

# ---- fixture parity: exact-mode results match Python ------------------------

test_that("ee_regression exact mode matches the Python reference", {
  for (model in c("linear", "logistic", "poisson")) {
    ref <- load_fixture(paste0("ee_regression_exact_", model))
    psi <- function(theta) {
      ee_regression(theta, X = ref$X, y = ref$y, model = model)
    }
    m <- MEstimator(stacked_equations = psi, init = ref$init)
    m <- estimate(m, deriv_method = "exact")
    expect_python_match(m, paste0("ee_regression_exact_", model), tolerance = 1e-6)
  }
})

test_that("ee_glm gamma exact mode matches the Python reference", {
  ref <- load_fixture("ee_glm_exact_gamma_log")
  psi <- function(theta) {
    ee_glm(theta, X = ref$X, y = ref$y, distribution = "gamma", link = "log")
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_glm_exact_gamma_log", tolerance = 1e-6)
})

test_that("ee_ridge_regression exact mode matches the Python reference", {
  ref <- load_fixture("ee_ridge_exact")
  psi <- function(theta) {
    ee_ridge_regression(theta, X = ref$X, y = ref$y, model = "linear",
                        penalty = ref$penalty)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_ridge_exact", tolerance = 1e-6)
})

test_that("ee_dlasso_regression exact mode matches the Python reference", {
  ref <- load_fixture("ee_dlasso_exact")
  psi <- function(theta) {
    ee_dlasso_regression(theta, X = ref$X, y = ref$y, model = "linear",
                         penalty = ref$penalty)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_dlasso_exact", tolerance = 1e-6)
})

test_that("ee_robust_regression cauchy loss exact mode matches the Python reference", {
  ref <- load_fixture("ee_robust_exact_cauchy")
  psi <- function(theta) {
    ee_robust_regression(theta, X = ref$X, y = ref$y, model = "linear",
                         k = ref$k, loss = "cauchy")
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_robust_exact_cauchy", tolerance = 1e-6)
})

# ---- internal consistency: exact agrees with the finite difference ----------

# Fit each equation twice from a shared setup and compare the exact bread and
# standard errors against the central-difference bread and standard errors. The
# point estimate is identical between the two fits (the derivative method only
# enters the bread), so this isolates the derivative computation.
# The bread is the direct target of the exact-mode fix, so it is always
# compared. The standard errors invert the bread, which amplifies any residual
# difference; for the L1 approximations (lasso, elastic net) the bread is nearly
# singular at the penalty kink, so a sub-1e-6 bread difference blows up on the
# standard-error scale. Those equations therefore compare the bread only, in
# keeping with the caution that the L1 penalties are not differentiable at the
# kink and the sandwich is not defined there.
expect_exact_matches_capprox <- function(psi, init, check_se = TRUE) {
  m_exact <- estimate(MEstimator(stacked_equations = psi, init = init),
                      deriv_method = "exact")
  m_cap <- estimate(MEstimator(stacked_equations = psi, init = init),
                    deriv_method = "capprox")
  expect_equal(unname(m_exact@bread), unname(m_cap@bread), tolerance = 1e-6)
  if (check_se) {
    expect_equal(sqrt(diag(m_exact@variance)), sqrt(diag(m_cap@variance)),
                 tolerance = 1e-6)
  }
}

test_that("exact and capprox derivatives agree for the regression families", {
  set.seed(2024)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  eta <- 0.6 + 1.1 * x1 - 0.5 * x2
  y_lin <- eta + rnorm(n)
  y_bin <- rbinom(n, 1, 1 / (1 + exp(-eta)))
  y_pois <- rpois(n, exp(eta - 0.5))
  y_gam <- rgamma(n, shape = 2, scale = exp(eta) / 2)

  expect_exact_matches_capprox(
    function(t) ee_regression(t, X = X, y = y_lin, model = "linear"),
    c(0, 0, 0)
  )
  expect_exact_matches_capprox(
    function(t) ee_regression(t, X = X, y = y_bin, model = "logistic"),
    c(0, 0, 0)
  )
  expect_exact_matches_capprox(
    function(t) ee_regression(t, X = X, y = y_pois, model = "poisson"),
    c(0, 0, 0)
  )
  expect_exact_matches_capprox(
    function(t) ee_glm(t, X = X, y = y_gam, distribution = "gamma", link = "log"),
    c(0, 0, 0, 0)
  )
  expect_exact_matches_capprox(
    function(t) ee_ridge_regression(t, X = X, y = y_lin, model = "linear",
                                    penalty = c(0, 5, 5)),
    c(0, 0, 0)
  )
  expect_exact_matches_capprox(
    function(t) suppressWarnings(
      ee_lasso_regression(t, X = X, y = y_lin, model = "linear",
                          penalty = c(0, 5, 5))
    ),
    c(0.1, 0.1, 0.1),
    check_se = FALSE
  )
  expect_exact_matches_capprox(
    function(t) suppressWarnings(
      ee_elasticnet_regression(t, X = X, y = y_lin, model = "linear",
                               penalty = c(0, 5, 5), ratio = 0.5)
    ),
    c(0.1, 0.1, 0.1),
    check_se = FALSE
  )
  # The cauchy loss redescends, so it needs a starting value near the solution
  # to avoid a degenerate fit that downweights every observation.
  expect_exact_matches_capprox(
    function(t) ee_robust_regression(t, X = X, y = y_lin, model = "linear",
                                     k = 1.345, loss = "cauchy"),
    c(0.6, 1.1, -0.5)
  )
})

# ---- multinomial and beta regression ----------------------------------------

# Multinomial logistic regression stacks a score per non-reference category, and
# the shared denominator 1 + sum(exp(X beta_j)) couples the categories, so every
# parameter carries tangents through the exponentials under exact mode. Beta
# regression runs the predicted mean through the digamma function twice and
# reshapes the tangent-carrying precision score into a single-row matrix. These
# fixtures pin the exact-mode bread, variance, and intervals against Python
# Delicatessen computed with deriv_method='exact'.
test_that("ee_mlogit exact mode matches the Python reference", {
  ref <- load_fixture("ee_mlogit_exact")
  psi <- function(theta) {
    ee_mlogit(theta, X = ref$X, y = ref$y)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_mlogit_exact", tolerance = 1e-6)
})

test_that("ee_beta_regression exact mode matches the Python reference", {
  ref <- load_fixture("ee_beta_regression_exact")
  psi <- function(theta) {
    ee_beta_regression(theta, X = ref$X, y = ref$y)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_beta_regression_exact", tolerance = 1e-6)
})

# The exact bread must also agree with the central-difference bread on the same
# fit, which isolates the tangent flow through the multinomial denominator and
# the beta digamma residual from the solver.
test_that("exact and capprox derivatives agree for mlogit and beta", {
  mlo <- load_fixture("ee_mlogit_exact")
  expect_exact_matches_capprox(
    function(t) ee_mlogit(t, X = mlo$X, y = mlo$y),
    mlo$init
  )
  bet <- load_fixture("ee_beta_regression_exact")
  expect_exact_matches_capprox(
    function(t) ee_beta_regression(t, X = bet$X, y = bet$y),
    bet$init
  )
})

# ---- clipping and branching robust losses -----------------------------------

# The Huber score clips residuals with pmin/pmax, and the Tukey, Andrew, and
# Hampel scores branch on the residual magnitude. Under exact mode the residual
# carries tangents, so the clip and the branch selection must be tangent-aware:
# the tangent follows the clipped or selected arm, matching numpy's np.clip and
# np.where on primal-tangent pairs. These fixtures pin the exact-mode bread,
# variance, and intervals against Python Delicatessen computed with
# deriv_method='exact' for each loss.
test_that("ee_robust_regression clipping and branching losses match Python", {
  for (loss in c("huber", "tukey", "andrew", "hampel")) {
    ref <- load_fixture(paste0("ee_robust_exact_", loss))
    psi <- function(theta) {
      ee_robust_regression(theta, X = ref$X, y = ref$y, model = "linear",
                           k = ref$k, loss = loss)
    }
    m <- MEstimator(stacked_equations = psi, init = ref$init)
    m <- estimate(m, deriv_method = "exact")
    expect_python_match(m, paste0("ee_robust_exact_", loss), tolerance = 1e-6)
  }
})

# The exact bread must also agree with the central-difference bread on the same
# fit, which isolates the tangent-aware clip and branch selection from the
# solver. Each fit starts near the robust solution so the redescending losses
# settle on the intended root.
test_that("exact and capprox derivatives agree for clipping and branching losses", {
  set.seed(2025)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  y <- 0.7 + 1.2 * x1 - 0.6 * x2 + rnorm(n, sd = 2)

  expect_exact_matches_capprox(
    function(t) ee_robust_regression(t, X = X, y = y, model = "linear",
                                     k = 1.345, loss = "huber"),
    c(0.7, 1.2, -0.6)
  )
  expect_exact_matches_capprox(
    function(t) ee_robust_regression(t, X = X, y = y, model = "linear",
                                     k = 4.685, loss = "tukey"),
    c(0.7, 1.2, -0.6)
  )
  expect_exact_matches_capprox(
    function(t) ee_robust_regression(t, X = X, y = y, model = "linear",
                                     k = 1.339, loss = "andrew"),
    c(0.7, 1.2, -0.6)
  )
  expect_exact_matches_capprox(
    function(t) ee_robust_regression(t, X = X, y = y, model = "linear",
                                     k = c(1.5, 3, 6), loss = "hampel"),
    c(0.7, 1.2, -0.6)
  )
})

# ---- normal CDF and PDF on theta-derived quantities -------------------------

# The Tobit score, the log-normal AFT survival term, and the probit inverse link
# evaluate the standard normal CDF and PDF on quantities derived from theta.
# Under exact mode those quantities carry tangents, so the CDF and PDF must be
# tangent-aware (the CDF tangent is the PDF, the PDF tangent is -x times the
# PDF). The Tobit score also clips the CDF away from zero with pmax, which must
# carry the tangent of the selected side. These fixtures pin the exact-mode
# bread, variance, and intervals against Python Delicatessen computed with
# deriv_method='exact'.
test_that("ee_tobit exact mode matches the Python reference", {
  ref <- load_fixture("ee_tobit_exact")
  psi <- function(theta) {
    ee_tobit(theta, X = ref$X, y = ref$y, lower = ref$lower)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_tobit_exact", tolerance = 1e-6)
})

test_that("ee_aft log-normal exact mode matches the Python reference", {
  ref <- load_fixture("ee_aft_exact_lognormal")
  psi <- function(theta) {
    ee_aft(theta, X = ref$X, time = ref$t, event = ref$delta,
           distribution = "log-normal")
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_aft_exact_lognormal", tolerance = 1e-6)
})

test_that("ee_glm probit exact mode matches the Python reference", {
  ref <- load_fixture("ee_glm_exact_probit")
  psi <- function(theta) {
    ee_glm(theta, X = ref$X, y = ref$y, distribution = "binomial",
           link = "probit")
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_glm_exact_probit", tolerance = 1e-6)
})

# The exact bread must also agree with the central-difference bread on the same
# fit, which isolates the tangent-aware normal CDF and PDF from the solver. Each
# fit starts from the fixture's initial values.
test_that("exact and capprox derivatives agree for the normal CDF and PDF families", {
  tob <- load_fixture("ee_tobit_exact")
  expect_exact_matches_capprox(
    function(t) ee_tobit(t, X = tob$X, y = tob$y, lower = tob$lower),
    tob$init
  )
  aft <- load_fixture("ee_aft_exact_lognormal")
  expect_exact_matches_capprox(
    function(t) ee_aft(t, X = aft$X, time = aft$t, event = aft$delta,
                       distribution = "log-normal"),
    aft$init
  )
  pro <- load_fixture("ee_glm_exact_probit")
  expect_exact_matches_capprox(
    function(t) ee_glm(t, X = pro$X, y = pro$y, distribution = "binomial",
                       link = "probit"),
    pro$init
  )
})
