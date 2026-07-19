# Tests for advanced causal inference estimating equations
# ee_ipw_msm, ee_gestimation_snmm, ee_iv_causal

# Helper: generate simple confounded data for MSM / g-estimation
make_causal_data_advanced <- function(n = 200, seed = 42) {
  set.seed(seed)
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.25 + 0.5 * W)
  Ya0 <- 5 + 2 * W + rnorm(n)
  Ya1 <- Ya0 - 2
  Y <- (1 - A) * Ya0 + A * Ya1

  list(
    Y = Y,
    A = A,
    W = W,
    # Propensity score design matrix: intercept, W
    Wmat = cbind(1, W),
    # MSM design matrix: intercept, A
    Vmsm = cbind(1, A),
    # SMM design matrix (V for g-estimation): intercept only
    Vsmm = cbind(Intercept = rep(1, n)),
    n = n
  )
}

# ee_ipw_msm tests ------------------------------------------------------------

test_that("ee_ipw_msm returns (c+b)-by-n matrix", {
  d <- make_causal_data_advanced()
  c_params <- ncol(d$Vmsm)
  b_params <- ncol(d$Wmat)
  theta <- rep(0, c_params + b_params)

  result <- ee_ipw_msm(
    theta,
    y = d$Y,
    A = d$A,
    W = d$Wmat,
    V = d$Vmsm,
    distribution = "normal",
    link = "identity"
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), c_params + b_params)
  expect_equal(ncol(result), d$n)
})

test_that("ee_ipw_msm solves via MEstimator with normal/identity", {
  d <- make_causal_data_advanced()
  c_params <- ncol(d$Vmsm)
  b_params <- ncol(d$Wmat)

  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vmsm,
      distribution = "normal",
      link = "identity"
    )
  }

  m <- MEstimator(stacked_equations = psi, init = rep(0, c_params + b_params))
  m <- estimate(m)

  # Should converge to finite values

  expect_true(all(is.finite(m@theta)))
  # MSM coefficient for A (theta[2]) should be negative (true effect = -2)
  expect_true(m@theta[2] < 0)
  # Should have positive variance
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_ipw_msm with binomial/logit distribution works", {
  set.seed(42)
  n <- 300
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.25 + 0.5 * W)
  Y <- rbinom(n, 1, plogis(0.5 - 0.5 * A - 0.3 * W))
  Wmat <- cbind(1, W)
  Vmsm <- cbind(1, A)

  c_params <- ncol(Vmsm)
  b_params <- ncol(Wmat)

  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = Y,
      A = A,
      W = Wmat,
      V = Vmsm,
      distribution = "binomial",
      link = "logit"
    )
  }

  m <- MEstimator(stacked_equations = psi, init = rep(0, c_params + b_params))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_ipw_msm with truncation works", {
  d <- make_causal_data_advanced()
  c_params <- ncol(d$Vmsm)
  b_params <- ncol(d$Wmat)

  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vmsm,
      distribution = "normal",
      link = "identity",
      truncate = c(0.1, 0.9)
    )
  }

  m <- MEstimator(stacked_equations = psi, init = rep(0, c_params + b_params))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
})

# ee_ipw_msm hyperparameter (tweedie MSM) -------------------------------------
#
# ee_ipw_msm carries a `hyperparameter` argument that flows straight through to
# the marginal structural model's weighted GLM: distribution = "tweedie" sets
# the variance function v(mu) = mu^p, where p is the fixed (not estimated)
# hyperparameter. The stacked theta is partitioned as ncol(V) MSM coefficients
# followed by ncol(W) propensity score coefficients, with no slot for an
# estimated nuisance parameter. Tweedie estimates no nuisance parameter, so this
# partition is exact; the gamma and negative binomial families, which read the
# last theta element as a log-nuisance parameter, cannot be used as an MSM
# outcome family here (see the non-support test below), matching Python
# Delicatessen, which raises a shape-alignment error for those families.
#
# These tests are backed by fixtures from generate_ipw_msm_hyperparameter_fixtures.py.
# The MSM design includes a continuous effect-modifier column so that the
# tweedie power genuinely reweights observations: with a design saturated in the
# binary exposure the power cancels and every hyperparameter gives the same
# inverse-probability-weighted group means. Fixture comparisons use tolerance
# 1e-5, matching the established causal-fixture precedent (the sibling
# ee_ipw_msm fixture test in test-python-fixtures.R uses 1e-5): the ee_ipw_msm
# sandwich stacks a numerically differentiated bread across the propensity model
# and the IPW-weighted GLM, so cross-platform agreement of the full variance is
# not guaranteed to the tighter 1e-6.

test_that("ee_ipw_msm tweedie stacked EE has ncol(V)+ncol(W) rows at init", {
  ref <- load_fixture("ee_ipw_msm_tweedie_p15")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  V <- ref$V

  # The raw stacked estimating function must return exactly ncol(V)+ncol(W) rows:
  # the tweedie MSM contributes ncol(V) score rows (no nuisance row appended) and
  # the propensity model contributes ncol(W). This requires the hyperparameter to
  # reach ee_glm; without it tweedie has no variance power and cannot evaluate.
  ee <- ee_ipw_msm(
    ref$init,
    y = y,
    A = A,
    W = W,
    V = V,
    distribution = "tweedie",
    link = "log",
    hyperparameter = ref$hyperparameter
  )

  expect_equal(nrow(ee), ref$ee_nrow)
  expect_equal(nrow(ee), ncol(V) + ncol(W))
  expect_equal(ncol(ee), length(y))

  # Row sums at the starting values pin the exact form of every stacked equation.
  # This is a direct evaluation with no sandwich, so the tighter 1e-6 applies.
  expect_equal(rowSums(ee), ref$ee_sum_at_init, tolerance = 1e-6)
})

test_that("ee_ipw_msm tweedie hyperparameter=1.5 matches Python Delicatessen", {
  ref <- load_fixture("ee_ipw_msm_tweedie_p15")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  V <- ref$V

  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = y,
      A = A,
      W = W,
      V = V,
      distribution = "tweedie",
      link = "log",
      hyperparameter = ref$hyperparameter
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  # theta has length ncol(V) + ncol(W); tweedie estimates no extra parameter.
  expect_length(m@theta, ncol(V) + ncol(W))
  expect_python_match(m, "ee_ipw_msm_tweedie_p15", tolerance = 1e-5)
})

test_that("ee_ipw_msm tweedie hyperparameter changes the MSM fit", {
  # Same data, two tweedie powers. If the hyperparameter were ignored the two
  # fits would be identical; instead the MSM coefficients must move while the
  # propensity score coefficients stay put (the propensity model never sees the
  # outcome distribution).
  ref15 <- load_fixture("ee_ipw_msm_tweedie_p15")
  ref18 <- load_fixture("ee_ipw_msm_tweedie_p18")

  y <- ref15$y
  A <- ref15$a
  W <- ref15$W
  V <- ref15$V

  fit <- function(p) {
    estimate(MEstimator(
      stacked_equations = function(theta) {
        ee_ipw_msm(
          theta,
          y = y,
          A = A,
          W = W,
          V = V,
          distribution = "tweedie",
          link = "log",
          hyperparameter = p
        )
      },
      init = ref15$init
    ))
  }

  m15 <- fit(ref15$hyperparameter)
  m18 <- fit(ref18$hyperparameter)

  # The two solutions genuinely differ. The largest coefficient shift is about
  # 5e-3, three orders of magnitude above the lm solver's ~1e-9 root-finding
  # noise, so this is not a tolerance-floor artifact.
  n_msm <- ncol(V)
  expect_gt(max(abs(unname(m15@theta) - unname(m18@theta))), 1e-3)

  # The change is confined to the MSM coefficients; the propensity score
  # coefficients (the trailing ncol(W)) are invariant to the outcome power to
  # full precision, since the tweedie power enters only the MSM's variance
  # function.
  ps_idx <- seq.int(n_msm + 1L, length(m15@theta))
  expect_equal(
    unname(m15@theta)[ps_idx],
    unname(m18@theta)[ps_idx],
    tolerance = 1e-8
  )

  # And each fit lands on its own Python reference.
  expect_python_match(m15, "ee_ipw_msm_tweedie_p15", tolerance = 1e-5)
  expect_python_match(m18, "ee_ipw_msm_tweedie_p18", tolerance = 1e-5)
})

test_that("ee_ipw_msm rejects nuisance-parameter outcome families", {
  # Python's ee_ipw_msm cannot fit a gamma or negative binomial MSM: those GLM
  # families read the last theta element as a log-nuisance parameter, but the
  # theta partition reserves exactly ncol(V) slots for the MSM with no room for
  # it, so ee_glm's design matrix multiplies a coefficient vector one element
  # short and the multiplication is non-conformable (Python raises the same
  # shape-alignment failure). The R implementation must mirror this non-support
  # rather than silently invent a capability Python lacks.
  ref <- load_fixture("ee_ipw_msm_tweedie_p15")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  V <- ref$V

  # init has length ncol(V) + ncol(W), exactly as for a supported family.
  expect_error(
    ee_ipw_msm(
      ref$init,
      y = y,
      A = A,
      W = W,
      V = V,
      distribution = "negative_binomial",
      link = "log"
    ),
    "conformable"
  )
  expect_error(
    ee_ipw_msm(
      ref$init,
      y = y,
      A = A,
      W = W,
      V = V,
      distribution = "gamma",
      link = "log"
    ),
    "conformable"
  )
})

# ee_gestimation_snmm tests ---------------------------------------------------

test_that("ee_gestimation_snmm inefficient returns (b+c)-by-n matrix", {
  d <- make_causal_data_advanced()
  b_params <- ncol(d$Vsmm) # SMM parameters
  c_params <- ncol(d$Wmat) # PS model parameters
  theta <- rep(0, b_params + c_params)

  result <- ee_gestimation_snmm(
    theta,
    y = d$Y,
    A = d$A,
    W = d$Wmat,
    V = d$Vsmm
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), b_params + c_params)
  expect_equal(ncol(result), d$n)
})

test_that("ee_gestimation_snmm efficient returns (b+c+d)-by-n matrix", {
  d <- make_causal_data_advanced()
  b_params <- ncol(d$Vsmm) # SMM parameters
  c_params <- ncol(d$Wmat) # PS model parameters
  Xmat <- d$Wmat # Outcome model design matrix
  d_params <- ncol(Xmat)
  theta <- rep(0, b_params + c_params + d_params)

  result <- ee_gestimation_snmm(
    theta,
    y = d$Y,
    A = d$A,
    W = d$Wmat,
    V = d$Vsmm,
    X = Xmat
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), b_params + c_params + d_params)
  expect_equal(ncol(result), d$n)
})

test_that("ee_gestimation_snmm inefficient linear solves via MEstimator", {
  d <- make_causal_data_advanced()
  b_params <- ncol(d$Vsmm)
  c_params <- ncol(d$Wmat)

  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vsmm,
      model = "linear"
    )
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = rep(0, b_params + c_params)
  )
  m <- estimate(m)

  # Should converge
  expect_true(all(is.finite(m@theta)))
  # SMM parameter (theta[1]) should be near -2 (true causal effect)
  expect_true(m@theta[1] < 0)
  # Variance should be positive
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_gestimation_snmm efficient linear solves via MEstimator", {
  d <- make_causal_data_advanced()
  b_params <- ncol(d$Vsmm)
  c_params <- ncol(d$Wmat)
  Xmat <- d$Wmat
  d_params <- ncol(Xmat)

  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vsmm,
      X = Xmat,
      model = "linear"
    )
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = rep(0, b_params + c_params + d_params)
  )
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  # SMM parameter should be near -2
  expect_true(m@theta[1] < 0)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_gestimation_snmm with V having multiple columns works", {
  set.seed(42)
  n <- 200
  W1 <- rnorm(n)
  V1 <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.25 + 0.5 * V1 + W1))
  Y <- 5 - 2 * A - 0.5 * A * V1 + 3 * V1 + W1 + rnorm(n)
  Wmat <- cbind(1, V1, W1)
  Vsmm <- cbind(1, V1)

  b_params <- ncol(Vsmm)
  c_params <- ncol(Wmat)

  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = Y,
      A = A,
      W = Wmat,
      V = Vsmm,
      model = "linear"
    )
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = rep(0, b_params + c_params)
  )
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  # theta[1] should be near -2, theta[2] near -0.5
  expect_true(m@theta[1] < 0)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_gestimation_snmm poisson model works", {
  set.seed(42)
  n <- 300
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * W))
  # Positive outcomes for Poisson SMM
  Y <- exp(1 + 0.5 * W - 0.3 * A + rnorm(n, sd = 0.3))
  Wmat <- cbind(1, W)
  Vsmm <- cbind(Intercept = rep(1, n))

  b_params <- ncol(Vsmm)
  c_params <- ncol(Wmat)

  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = Y,
      A = A,
      W = Wmat,
      V = Vsmm,
      model = "poisson"
    )
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = rep(0, b_params + c_params)
  )
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  # Should estimate a negative effect
  expect_true(m@theta[1] < 0)
})

test_that("ee_gestimation_snmm errors on unsupported model", {
  d <- make_causal_data_advanced()
  b_params <- ncol(d$Vsmm)
  c_params <- ncol(d$Wmat)
  theta <- rep(0, b_params + c_params)

  expect_error(
    ee_gestimation_snmm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vsmm,
      model = "logistic"
    ),
    "not supported"
  )
})

# ee_iv_causal tests -----------------------------------------------------------

test_that("ee_iv_causal returns 2-by-n matrix", {
  set.seed(42)
  n <- 200
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(Z))
  Y <- 2 * A + rnorm(n)

  theta <- c(0, 0.5)
  result <- ee_iv_causal(theta, y = Y, A = A, Z = Z)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), n)
})

test_that("ee_iv_causal solves via MEstimator", {
  set.seed(42)
  n <- 500
  Z <- rbinom(n, 1, 0.5)
  U <- rnorm(n)
  A <- rbinom(n, 1, plogis(U + Z))
  Y <- 2 * A - U + rnorm(n)

  psi <- function(theta) {
    ee_iv_causal(theta, y = Y, A = A, Z = Z)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5))
  m <- estimate(m)

  # Should converge to finite values
  expect_true(all(is.finite(m@theta)))
  # theta[1] is the causal effect (true = 2), should be positive
  expect_true(m@theta[1] > 0)
  # theta[2] is mean of Z, should be near 0.5
  expect_true(abs(m@theta[2] - 0.5) < 0.15)
  # Should have positive variance
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_iv_causal with weights works", {
  set.seed(42)
  n <- 300
  Z <- rbinom(n, 1, 0.5)
  U <- rnorm(n)
  A <- rbinom(n, 1, plogis(U + Z))
  Y <- 2 * A - U + rnorm(n)
  wts <- runif(n, 0.5, 1.5)

  psi <- function(theta) {
    ee_iv_causal(theta, y = Y, A = A, Z = Z, weights = wts)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(m@theta[1] > 0)
})

test_that("ee_iv_causal recovers known causal effect", {
  # Use a larger sample for more precise estimation
  set.seed(123)
  n <- 2000
  Z <- rbinom(n, 1, 0.5)
  U <- rnorm(n)
  # Strong instrument
  A <- rbinom(n, 1, plogis(-1 + 3 * Z + U))
  # True causal effect = 3
  Y <- 3 * A - U + rnorm(n, sd = 0.5)

  psi <- function(theta) {
    ee_iv_causal(theta, y = Y, A = A, Z = Z)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5))
  m <- estimate(m)

  # Should be reasonably close to 3
  expect_true(abs(m@theta[1] - 3) < 1.5)
  # Mean of Z should be near 0.5
  expect_true(abs(m@theta[2] - 0.5) < 0.05)
})

# ee_2sls tests ----------------------------------------------------------------

test_that("ee_2sls returns correct dimension matrix without W", {
  set.seed(42)
  n <- 200
  Z <- cbind(1, rbinom(n, 1, 0.5))
  U <- rnorm(n)
  A <- as.numeric(plogis(U + Z[, 2]) > runif(n))
  Y <- 2 * A - U + rnorm(n)

  # Without W: 1 second-stage param + 2 first-stage params = 3
  theta <- rep(0, 3)
  result <- ee_2sls(theta, y = Y, A = A, Z = Z)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

test_that("ee_2sls returns correct dimension matrix with W", {
  set.seed(42)
  n <- 200
  Z <- cbind(rbinom(n, 1, 0.5))
  W <- cbind(1, rnorm(n))
  U <- rnorm(n)
  A <- as.numeric(plogis(U + Z[, 1]) > runif(n))
  Y <- 2 * A - U + rnorm(n)

  # With W (2 cols): second-stage = 1+2=3, first-stage = 1+2=3, total = 6
  theta <- rep(0, 6)
  result <- ee_2sls(theta, y = Y, A = A, Z = Z, W = W)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 6)
  expect_equal(ncol(result), n)
})

test_that("ee_2sls solves via MEstimator with W (intercept + covariate)", {
  set.seed(42)
  n <- 500
  X_cov <- rnorm(n)
  Z_raw <- rbinom(n, 1, 0.5)
  Z <- cbind(Z_raw)
  U <- rnorm(n)
  # Continuous treatment with strong instrument
  A <- 0.5 * Z_raw + 0.3 * X_cov + U + rnorm(n, sd = 0.5)
  Y <- 2 * A - U + 0.5 * X_cov + rnorm(n)
  W <- cbind(1, X_cov)

  # Get OLS-based initial values for reliable convergence
  dmatrix1 <- cbind(Z, W)
  alpha_init <- as.numeric(coef(lm(A ~ dmatrix1 - 1)))
  a_hat_init <- as.numeric(dmatrix1 %*% alpha_init)
  dmatrix2 <- cbind(a_hat_init, W)
  beta_init <- as.numeric(coef(lm(Y ~ dmatrix2 - 1)))

  psi <- function(theta) {
    ee_2sls(theta, y = Y, A = A, Z = Z, W = W)
  }

  m <- MEstimator(stacked_equations = psi, init = c(beta_init, alpha_init))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Second-stage coefficient for A_hat (theta[1]) should be positive
  expect_true(m@theta[1] > 0)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_2sls with multiple instruments and W works", {
  set.seed(42)
  n <- 500
  Z1 <- rbinom(n, 1, 0.5)
  Z2 <- rnorm(n)
  U <- rnorm(n)
  # Continuous treatment with strong instruments
  A <- 0.5 * Z1 + 0.3 * Z2 + U + rnorm(n, sd = 0.5)
  Y <- 2 * A - U + rnorm(n)
  Z <- cbind(Z1, Z2)
  W <- cbind(rep(1, n))

  # Get OLS-based initial values
  dmatrix1 <- cbind(Z, W)
  alpha_init <- as.numeric(coef(lm(A ~ dmatrix1 - 1)))
  a_hat_init <- as.numeric(dmatrix1 %*% alpha_init)
  dmatrix2 <- cbind(a_hat_init, W)
  beta_init <- as.numeric(coef(lm(Y ~ dmatrix2 - 1)))

  psi <- function(theta) {
    ee_2sls(theta, y = Y, A = A, Z = Z, W = W)
  }

  m <- MEstimator(stacked_equations = psi, init = c(beta_init, alpha_init))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(m@theta[1] > 0)
  expect_true(all(diag(m@variance) > 0))
})

# ee_mean_sensitivity_analysis tests -------------------------------------------

test_that("ee_mean_sensitivity_analysis returns (1+b)-by-n matrix", {
  set.seed(42)
  n <- 200
  W <- rbinom(n, 1, 0.5)
  Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
  M <- rbinom(n, 1, plogis(2 + W - 0.01 * Y_full))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)

  q_eval <- ifelse(is.na(Y), 0, -0.01 * Y)

  theta <- c(180, rep(0, ncol(X)))
  result <- ee_mean_sensitivity_analysis(
    theta,
    y = ifelse(is.na(Y), 0, Y),
    delta = M,
    X = X,
    q_eval = q_eval,
    H_function = inverse_logit
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1 + ncol(X))
  expect_equal(ncol(result), n)
})

test_that("ee_mean_sensitivity_analysis solves via MEstimator (alpha=0, MAR)", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
  M <- rbinom(n, 1, plogis(2 + W))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)

  # alpha = 0 means no sensitivity (MAR)
  psi <- function(theta) {
    y_filled <- ifelse(is.na(Y), 0, Y)
    q_eval <- rep(0, n)
    ee_mean_sensitivity_analysis(
      theta,
      y = y_filled,
      delta = M,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(180, 0, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Mean (theta[1]) should be in a reasonable range
  expect_true(m@theta[1] > 150 && m@theta[1] < 220)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_mean_sensitivity_analysis with nonzero alpha works", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
  M <- rbinom(n, 1, plogis(2 + W - 0.01 * Y_full))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)

  psi <- function(theta) {
    y_filled <- ifelse(is.na(Y), 0, Y)
    q_eval <- -0.01 * y_filled
    ee_mean_sensitivity_analysis(
      theta,
      y = y_filled,
      delta = M,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(180, 0, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(m@theta[1] > 150 && m@theta[1] < 220)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_mean_sensitivity_analysis sensitivity changes mean estimate", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
  M <- rbinom(n, 1, plogis(2 + W - 0.01 * Y_full))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)
  y_filled <- ifelse(is.na(Y), 0, Y)

  # alpha = 0
  psi0 <- function(theta) {
    q_eval <- rep(0, n)
    ee_mean_sensitivity_analysis(
      theta,
      y = y_filled,
      delta = M,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }
  m0 <- MEstimator(stacked_equations = psi0, init = c(180, 0, 0))
  m0 <- estimate(m0, solver = "nleqslv")

  # alpha = -0.05 (stronger sensitivity)
  psi1 <- function(theta) {
    q_eval <- -0.05 * y_filled
    ee_mean_sensitivity_analysis(
      theta,
      y = y_filled,
      delta = M,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }
  m1 <- MEstimator(stacked_equations = psi1, init = c(180, 0, 0))
  expect_warning(
    m1 <- estimate(m1, solver = "nleqslv"),
    "nleqslv did not converge"
  )

  # Different alphas should produce different mean estimates
  expect_false(abs(m0@theta[1] - m1@theta[1]) < 0.001)
})

test_that("ee_mean_sensitivity_analysis with a scalar intercept design emits no array-recycling warning", {
  set.seed(42)
  n <- 100
  y_full <- 200 + rnorm(n, sd = 5)
  r <- rbinom(n, 1, 0.7)
  y_filled <- ifelse(r == 0, 0, y_full)
  qy <- rep(0.01, n)

  recycling_warnings <- character()
  result <- withCallingHandlers(
    ee_mean_sensitivity_analysis(
      c(180, 0),
      y = y_filled,
      delta = r,
      X = 1,
      q_eval = qy,
      H_function = inverse_logit
    ),
    warning = function(w) {
      if (grepl("Recycling array", conditionMessage(w))) {
        recycling_warnings <<- c(recycling_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )

  expect_length(recycling_warnings, 0)

  # A scalar intercept design must agree with an explicit column of ones.
  explicit <- ee_mean_sensitivity_analysis(
    c(180, 0),
    y = y_filled,
    delta = r,
    X = matrix(1, nrow = n, ncol = 1),
    q_eval = qy,
    H_function = inverse_logit
  )
  expect_equal(result, explicit)
})

test_that("ee_mean_sensitivity_analysis with binary outcome works", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  Y_full <- rbinom(n, 1, plogis(0.5 - 0.5 * W))
  M <- rbinom(n, 1, plogis(1 + 0.5 * W))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)
  y_filled <- ifelse(is.na(Y), 0, Y)

  psi <- function(theta) {
    q_eval <- rep(0, n)
    ee_mean_sensitivity_analysis(
      theta,
      y = y_filled,
      delta = M,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.5, 0, 0))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Binary mean should be between 0 and 1
  expect_true(m@theta[1] > 0 && m@theta[1] < 1)
  expect_true(all(diag(m@variance) > 0))
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_iv_causal rejects weights of the wrong length", {
  set.seed(3)
  n <- 40
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.3 + 0.4 * Z)
  Y <- rnorm(n, mean = 1 + 0.5 * A)

  expect_error(
    ee_iv_causal(
      c(0.5, 0.5),
      y = Y,
      A = A,
      Z = Z,
      weights = rep(1, n / 2)
    ),
    "same length as the data"
  )
})

test_that("ee_ipw_msm rejects weights of the wrong length", {
  d <- make_causal_data_advanced(n = 40)
  c_params <- ncol(d$Vmsm)
  b_params <- ncol(d$Wmat)
  theta <- rep(0, c_params + b_params)

  expect_error(
    ee_ipw_msm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = d$Vmsm,
      distribution = "normal",
      link = "identity",
      weights = rep(1, d$n / 2)
    ),
    "same length as the data"
  )
})

test_that("ee_mean_sensitivity_analysis rejects a short q_eval", {
  ref <- load_fixture("ee_mean_sensitivity_analysis")
  y <- ref$y
  delta <- ref$delta
  X <- ref$X
  half <- seq_len(length(y) / 2)

  expect_error(
    ee_mean_sensitivity_analysis(
      ref$init,
      y = y,
      delta = delta,
      X = X,
      q_eval = ref$q_eval[half],
      H_function = inverse_logit
    ),
    "same length as the data"
  )
})

test_that("ee_mean_sensitivity_analysis rejects a short delta", {
  ref <- load_fixture("ee_mean_sensitivity_analysis")
  y <- ref$y
  X <- ref$X
  half <- seq_len(length(y) / 2)

  expect_error(
    ee_mean_sensitivity_analysis(
      ref$init,
      y = y,
      delta = ref$delta[half],
      X = X,
      q_eval = ref$q_eval,
      H_function = inverse_logit
    ),
    "same length as the data"
  )
})
