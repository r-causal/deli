# The parameter names the built-in estimating equations write on their own
# returns, and the labels those row names produce on a fit.
#
# Each equation covered here labels every row it returns, the nuisance
# regression blocks included. That is not decoration: psi_param_names() discards
# a row-name vector that leaves any parameter unlabeled, so a blank nuisance
# block would cost the causal estimates stacked above it their names as well.
# The nuisance labels are positional and say only what is knowable, since
# coerce_design() drops the caller's column headings, so `W_2` is the
# coefficient on the second column of `W` and claims nothing further.
#
# test-param-names.R holds the rules this file exercises, and
# test-exact-mode-acceptance.R holds the exhaustive exact-versus-capprox
# comparison for every family. What is pinned here is the names.

# ---- Helpers ----

# Right-censored times from a Weibull with shape 1.2, large enough for all three
# parametric survival distributions to solve from a plain start.
survival_times <- function() {
  set.seed(123)
  n <- 300
  t_true <- (-log(runif(n)) / 0.5)^(1 / 1.2)
  censor <- rexp(n, rate = 0.2)
  list(
    time = pmin(t_true, censor),
    event = as.numeric(t_true <= censor)
  )
}

# A hyperbolic E-max dose-response, and a four-parameter log-logistic one on a
# log-spaced dose grid.
emax_doses <- function() {
  set.seed(42)
  n <- 50
  dose <- runif(n, 0, 100)
  list(dose = dose, response = 1 + (25 * dose) / (10 + dose) + rnorm(n, sd = 1))
}

loglogistic_doses <- function() {
  set.seed(42)
  n <- 40
  dose <- exp(runif(n, log(0.1), log(100)))
  rho <- (dose / 10)^2
  list(dose = dose, response = 5 + (25 - 5) / (1 + rho) + rnorm(n, sd = 0.5))
}

# ee-basic.R --------------------------------------------------------------

test_that("ee_mean_variance names the mean and the variance", {
  y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
  psi <- function(theta) ee_mean_variance(theta, y = y)

  expect_equal(rownames(psi(c(0, 1))), c("mean", "variance"))

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("mean", "variance"))
})

test_that("ee_positive_mean_deviation names the deviation and the median", {
  skip_if_not_installed("minpack.lm")

  y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
  # The equation warns on every call that its median row is not
  # differentiable, so both the direct call and the fit suppress it.
  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }

  expect_equal(rownames(psi(c(1, 3))), c("positive_mean_deviation", "median"))

  # The documented starting values and solver, since no solver can search for
  # the median.
  init <- c(mean(2 * (y - median(y)) * (y > median(y))), median(y))
  m <- m_estimate(stacked_equations = psi, init = init, solver = "lm")

  expect_equal(names(coef(m)), c("positive_mean_deviation", "median"))
})

# ee-regression.R ---------------------------------------------------------

test_that("ee_tobit names the design rows and the log scale", {
  ref <- load_fixture("ee_tobit")
  psi <- function(theta) {
    ee_tobit(theta, X = ref$X, y = ref$y, lower = ref$lower)
  }
  expected <- c("X_1", "X_2", "log_sigma")

  expect_equal(ncol(ref$X), 2L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_beta_regression names the design rows and the log precision", {
  ref <- load_fixture("ee_beta_regression")
  psi <- function(theta) ee_beta_regression(theta, X = ref$X, y = ref$y)
  expected <- c("X_1", "X_2", "log_phi")

  expect_equal(ncol(ref$X), 2L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_regression leaves its rows unnamed", {
  # Every parameter is a coefficient on a design column, and coerce_design()
  # has already dropped the caller's column headings, so positional labels
  # would say no more than the numbering they would replace.
  ref <- load_fixture("ee_tobit")
  ee <- ee_regression(c(0, 0), X = ref$X, y = ref$y, model = "linear")

  expect_null(rownames(ee))
})

# ee-glm.R ----------------------------------------------------------------

test_that("ee_glm names the design rows and the log shape under the gamma", {
  ref <- load_fixture("ee_glm_gamma_log")
  psi <- function(theta) {
    ee_glm(theta, X = ref$X, y = ref$y, distribution = "gamma", link = "log")
  }
  expected <- c("X_1", "X_2", "X_3", "log_shape")

  expect_equal(ncol(ref$X), 3L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_glm names the log dispersion under the negative binomial", {
  ref <- load_fixture("ee_glm_negative_binomial_log")
  psi <- function(theta) {
    ee_glm(
      theta,
      X = ref$X,
      y = ref$y,
      distribution = "negative_binomial",
      link = "log"
    )
  }
  expected <- c("X_1", "X_2", "X_3", "log_dispersion")

  expect_equal(ncol(ref$X), 3L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("the nb alias names its rows as negative_binomial does", {
  ref <- load_fixture("ee_glm_nb_alias_log")
  psi <- function(theta) {
    ee_glm(theta, X = ref$X, y = ref$y, distribution = "nb", link = "log")
  }
  expected <- c("X_1", "X_2", "X_3", "log_dispersion")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_glm leaves the tweedie rows unnamed", {
  # The tweedie power is a fixed hyperparameter rather than an estimated
  # parameter, so nothing is stacked beneath the beta scores and every
  # parameter is a coefficient on a design column.
  ref <- load_fixture("ee_glm_tweedie_p15_log")
  psi <- function(theta) {
    ee_glm(
      theta,
      X = ref$X,
      y = ref$y,
      distribution = "tweedie",
      link = "log",
      hyperparameter = ref$hyperparameter
    )
  }

  expect_equal(length(ref$init), ncol(ref$X))
  expect_null(rownames(psi(ref$init)))

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), paste0("theta_", 1:3))
})

test_that("ee_glm leaves the distributions with no extra parameter unnamed", {
  for (fixture in c("ee_glm_normal_identity", "ee_glm_poisson_log")) {
    ref <- load_fixture(fixture)
    ee <- ee_glm(
      ref$init,
      X = ref$X,
      y = ref$y,
      distribution = ref$distribution,
      link = ref$link
    )

    expect_null(rownames(ee), label = fixture)
  }
})

# ee-pharma.R -------------------------------------------------------------

test_that("ee_emax names the zero-dose response, maximum, and ED50", {
  d <- emax_doses()
  psi <- function(theta) ee_emax(theta, dose = d$dose, response = d$response)

  expect_equal(rownames(psi(c(1, 25, 10))), c("e0", "emax", "ed50"))

  m <- m_estimate(stacked_equations = psi, init = c(1, 25, 10))

  expect_equal(names(coef(m)), c("e0", "emax", "ed50"))
})

test_that("ee_loglogistic names its four dose-response parameters", {
  d <- loglogistic_doses()
  psi <- function(theta) {
    ee_loglogistic(theta, dose = d$dose, response = d$response)
  }

  expect_equal(
    rownames(psi(c(5, 25, 10, 2))),
    c("lower", "upper", "ed50", "steepness")
  )

  m <- m_estimate(stacked_equations = psi, init = c(5, 25, 10, 2))

  expect_equal(names(coef(m)), c("lower", "upper", "ed50", "steepness"))
})

# ee-survival.R -----------------------------------------------------------

test_that("ee_survival_model names lambda alone under the exponential", {
  d <- survival_times()
  psi <- function(theta) {
    ee_survival_model(
      theta,
      time = d$time,
      event = d$event,
      distribution = "exponential"
    )
  }

  expect_equal(rownames(psi(0.5)), "lambda")

  m <- m_estimate(stacked_equations = psi, init = 0.5)

  expect_equal(names(coef(m)), "lambda")
})

test_that("ee_survival_model names lambda and gamma under the Weibull", {
  d <- survival_times()
  psi <- function(theta) {
    ee_survival_model(
      theta,
      time = d$time,
      event = d$event,
      distribution = "weibull"
    )
  }

  expect_equal(rownames(psi(c(0.3, 1))), c("lambda", "gamma"))

  m <- m_estimate(stacked_equations = psi, init = c(0.3, 1))

  expect_equal(names(coef(m)), c("lambda", "gamma"))
})

test_that("ee_survival_model names lambda and gamma under the Gompertz", {
  set.seed(345)
  n <- 300
  t_true <- (1 / 0.5) * log(1 - (0.5 / 0.3) * log(runif(n)))
  t_true <- t_true[is.finite(t_true) & t_true > 0]
  censor <- rexp(length(t_true), rate = 0.3)
  time <- pmin(t_true, censor)
  event <- as.numeric(t_true <= censor)

  psi <- function(theta) {
    ee_survival_model(
      theta,
      time = time,
      event = event,
      distribution = "gompertz"
    )
  }

  expect_equal(rownames(psi(c(0.2, 0.3))), c("lambda", "gamma"))

  m <- m_estimate(stacked_equations = psi, init = c(0.2, 0.3))

  expect_equal(names(coef(m)), c("lambda", "gamma"))
})

test_that("ee_aft names the design rows and the log-inverse scale", {
  ref <- load_fixture("ee_aft_weibull")
  psi <- function(theta) {
    ee_aft(
      theta,
      X = ref$X,
      time = ref$t,
      event = ref$delta,
      distribution = "weibull"
    )
  }
  expected <- c("X_1", "X_2", "log_inv_scale")

  expect_equal(ncol(ref$X), 2L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_aft leaves the exponential distribution's rows unnamed", {
  # Sigma is fixed at 1 for the exponential, so there is no log(1/sigma) row and
  # every parameter is a coefficient on a column of X. The counterpart in
  # test-param-names.R pins what a fit made from it is called.
  ref <- load_fixture("ee_aft_weibull")
  psi <- function(theta) {
    ee_aft(
      theta,
      X = ref$X,
      time = ref$t,
      event = ref$delta,
      distribution = "exponential"
    )
  }

  expect_null(rownames(psi(ref$init[1:2])))
})

# ee-measurement.R --------------------------------------------------------

test_that("ee_rogan_gladen names the two proportions and the two accuracies", {
  ref <- load_fixture("ee_rogan_gladen")
  psi <- function(theta) {
    ee_rogan_gladen(theta, y = ref$y, y_star = ref$y_star, r = ref$r)
  }
  expected <- c(
    "corrected_proportion",
    "naive_proportion",
    "sensitivity",
    "specificity"
  )

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_rogan_gladen_extended names its two accuracy models apart", {
  # The one nuisance block in the package that earns real names: the two models
  # fit the same design, so a shared prefix would repeat every label and cost
  # the return its names entirely.
  ref <- load_fixture("ee_rogan_gladen_extended_exact")
  psi <- function(theta) {
    ee_rogan_gladen_extended(
      theta,
      y = ref$y,
      y_star = ref$y_star,
      r = ref$r,
      X = ref$X
    )
  }
  expected <- c("corrected_proportion", "sens_1", "sens_2", "spec_1", "spec_2")

  expect_equal(ncol(ref$X), 2L)
  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_regression_calibration names the default calibration design", {
  # With no X the calibration design is the mismeasured action and an
  # intercept, both supplied by the function, so both are named for what they
  # are rather than positionally.
  ref <- load_fixture("ee_regression_calibration_exact")
  psi <- function(theta) {
    ee_regression_calibration(
      theta,
      beta = ref$beta,
      a = ref$a,
      a_star = ref$a_star,
      r = ref$r
    )
  }
  expected <- c("corrected_beta", "a_star", "intercept")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_regression_calibration names a supplied calibration design", {
  ref <- load_fixture("ee_regression_calibration_exact")
  set.seed(551)
  X <- base::cbind(1, rnorm(length(ref$a)))
  psi <- function(theta) {
    ee_regression_calibration(
      theta,
      beta = ref$beta,
      a = ref$a,
      a_star = ref$a_star,
      r = ref$r,
      X = X
    )
  }
  expected <- c("corrected_beta", "a_star", "X_1", "X_2")

  expect_equal(rownames(psi(c(ref$init, 0))), expected)

  m <- m_estimate(stacked_equations = psi, init = c(ref$init, 0))

  expect_equal(names(coef(m)), expected)
})

# ee-causal.R -------------------------------------------------------------

test_that("ee_gformula names a single plan's estimate the causal mean", {
  # The superscripts the two-plan branch uses index the pair of plans, and one
  # plan has no second element to index, so the row takes a bare name.
  ref <- load_fixture("ee_gformula_exact")
  psi <- function(theta) {
    ee_gformula(theta, y = ref$y, X = ref$X, X1 = ref$X1)
  }
  expected <- c("causal_mean", "X_1", "X_2", "X_3", "X_4")

  # Drop two of the three causal scalars from the two-plan starting values.
  init <- ref$init[-(1:2)]

  expect_equal(rownames(psi(init)), expected)

  m <- m_estimate(stacked_equations = psi, init = init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_gformula names the two-plan contrast and both plan means", {
  ref <- load_fixture("ee_gformula_exact")
  psi <- function(theta) {
    ee_gformula(theta, y = ref$y, X = ref$X, X1 = ref$X1, X0 = ref$X0)
  }
  expected <- c("ACE", "E[Y^1]", "E[Y^0]", "X_1", "X_2", "X_3", "X_4")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_ipw names the causal scalars ahead of the propensity model", {
  ref <- load_fixture("ee_ipw_exact")
  psi <- function(theta) ee_ipw(theta, y = ref$y, A = ref$a, W = ref$W)
  expected <- c("ACE", "E[Y^1]", "E[Y^0]", "W_1", "W_2", "W_3")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_aipw names the propensity and outcome blocks separately", {
  ref <- load_fixture("ee_aipw_exact")
  psi <- function(theta) {
    ee_aipw(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      X = ref$X,
      X1 = ref$X1,
      X0 = ref$X0
    )
  }
  expected <- c(
    "ACE",
    "E[Y^1]",
    "E[Y^0]",
    "W_1",
    "W_2",
    "W_3",
    "X_1",
    "X_2",
    "X_3",
    "X_4"
  )

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_ipw_msm subscripts the MSM parameters from zero", {
  # The same shape as ee_gestimation_snmm: a causal-model block over V stacked
  # on a propensity score block over W, with the block's own symbol taken from
  # the literature and the nuisance block named positionally.
  ref <- load_fixture("ee_ipw_msm")
  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      distribution = ref$distribution,
      link = ref$link
    )
  }
  expected <- c("MSM alpha_0", "MSM alpha_1", "W_1", "W_2", "W_3")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_gestimation_snmm subscripts the SNM parameters from zero", {
  ref <- load_fixture("ee_gestimation_exact")
  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      model = ref$model
    )
  }
  expected <- c("SNM phi_0", "SNM phi_1", "W_1", "W_2", "W_3")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_gestimation_snmm names the efficient estimator's outcome model", {
  ref <- load_fixture("ee_gestimation_efficient_exact")
  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      X = ref$X,
      model = ref$model
    )
  }
  expected <- c(
    "SNM phi_0",
    "SNM phi_1",
    "W_1",
    "W_2",
    "W_3",
    "X_1",
    "X_2",
    "X_3"
  )

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_iv_causal names the causal effect and the instrument mean", {
  # An unmeasured confounder biases A on Y, and Z reaches Y only through A.
  set.seed(123)
  n <- 500
  Z <- rbinom(n, 1, 0.5)
  U <- rnorm(n)
  A <- rbinom(n, 1, plogis(-1 + 3 * Z + U))
  y <- 3 * A - U + rnorm(n, sd = 0.5)
  psi <- function(theta) ee_iv_causal(theta, y = y, A = A, Z = Z)

  expect_equal(rownames(psi(c(0, 0.5))), c("causal_effect", "mean_Z"))

  m <- m_estimate(stacked_equations = psi, init = c(0, 0.5))

  expect_equal(names(coef(m)), c("causal_effect", "mean_Z"))
})

test_that("ee_2sls carries the stage in every name", {
  # The columns of W are fitted once per stage, so an unqualified name would
  # repeat and the whole return would lose its names.
  ref <- load_fixture("ee_2sls_exact_w")
  psi <- function(theta) {
    ee_2sls(theta, y = ref$y, A = ref$a, Z = ref$Z, W = ref$W)
  }
  expected <- c(
    "stage2_A",
    "stage2_W_1",
    "stage2_W_2",
    "stage1_Z_1",
    "stage1_W_1",
    "stage1_W_2"
  )

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_2sls keeps the stage prefixes without exogenous covariates", {
  ref <- load_fixture("ee_2sls_exact_now")
  psi <- function(theta) ee_2sls(theta, y = ref$y, A = ref$a, Z = ref$Z)
  expected <- c("stage2_A", "stage1_Z_1", "stage1_Z_2")

  expect_equal(rownames(psi(ref$init)), expected)

  m <- m_estimate(stacked_equations = psi, init = ref$init)

  expect_equal(names(coef(m)), expected)
})

test_that("ee_mean_sensitivity_analysis names the corrected mean first", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  delta <- rbinom(n, 1, plogis(2 + W))
  y <- base::ifelse(delta == 1, 200 - 35 * W + rnorm(n, sd = 5), 0)
  X <- base::cbind(1, W)
  psi <- function(theta) {
    ee_mean_sensitivity_analysis(
      theta,
      y = y,
      delta = delta,
      X = X,
      q_eval = rep(0, n),
      H_function = inverse_logit
    )
  }
  expected <- c("corrected_mean", "X_1", "X_2")

  expect_equal(rownames(psi(c(180, 0, 0))), expected)

  m <- m_estimate(stacked_equations = psi, init = c(180, 0, 0))

  expect_equal(names(coef(m)), expected)
})

test_that("a scalar missingness design still names its one coefficient", {
  # A scalar design arrives as a 1-by-1 matrix and reaches the branch that
  # strips the dim off X, after which ncol() answers NULL. The column count is
  # read before that branch, and this is the case that would abort without the
  # hoist. A one-column matrix does not reach it: the branch tests length(X),
  # which is n rows for a matrix and 1 only for the scalar.
  set.seed(88)
  n <- 200
  delta <- rbinom(n, 1, 0.85)
  y <- base::ifelse(delta == 1, 12 + rnorm(n), 0)
  psi <- function(theta) {
    ee_mean_sensitivity_analysis(
      theta,
      y = y,
      delta = delta,
      X = 1,
      q_eval = rep(0, n),
      H_function = inverse_logit
    )
  }

  expect_equal(rownames(psi(c(12, 0))), c("corrected_mean", "X_1"))

  m <- m_estimate(stacked_equations = psi, init = c(12, 0))

  expect_equal(names(coef(m)), c("corrected_mean", "X_1"))
})

# The naming rules applied to the built-in equations ----------------------

test_that("names on init override the ones a built-in equation writes", {
  ref <- load_fixture("ee_ipw_exact")
  psi <- function(theta) ee_ipw(theta, y = ref$y, A = ref$a, W = ref$W)
  init <- ref$init
  names(init) <- c("ate", "y1", "y0", "b0", "b1", "b2")

  m <- m_estimate(stacked_equations = psi, init = init)

  expect_equal(names(coef(m)), c("ate", "y1", "y0", "b0", "b1", "b2"))
})

test_that("stacking two ee_ipw blocks repeats the labels and numbers the fit", {
  # Two blocks of the same equation name the same six parameters, and a
  # repeated label is worse than none: coef() would show it twice and
  # confint(m)["ACE", ] would silently return only the first. The stack is
  # numbered instead. Name `init` to label a stack like this.
  ref <- load_fixture("ee_ipw_exact")
  # base:: because the test environment is a clone of the package namespace and
  # would otherwise reach deli's masked binder, which a user stacking two blocks
  # in the global environment does not.
  psi <- function(theta) {
    base::rbind(
      ee_ipw(theta[1:6], y = ref$y, A = ref$a, W = ref$W),
      ee_ipw(theta[7:12], y = 2 * ref$y, A = ref$a, W = ref$W)
    )
  }
  init <- c(ref$init, ref$init)
  one_block <- c("ACE", "E[Y^1]", "E[Y^0]", "W_1", "W_2", "W_3")

  expect_equal(rownames(psi(init)), rep(one_block, 2))

  m <- m_estimate(stacked_equations = psi, init = init)

  expect_equal(names(coef(m)), paste0("theta_", 1:12))
})

test_that("a named built-in block stacked with an unnamed one is numbered", {
  # ee_emax_ed() writes no row names, so the stack is incomplete and every
  # parameter is numbered, exactly as it was before ee_emax() named its rows.
  d <- emax_doses()
  psi <- function(theta) {
    base::rbind(
      ee_emax(theta[1:3], dose = d$dose, response = d$response),
      ee_emax_ed(theta[4], dose = d$dose, delta = 0.8, ed50 = theta[3])
    )
  }

  expect_equal(rownames(psi(c(1, 25, 10, 40))), c("e0", "emax", "ed50", ""))

  m <- m_estimate(stacked_equations = psi, init = c(1, 25, 10, 40))

  expect_equal(names(coef(m)), paste0("theta_", 1:4))
})

# Exact differentiation ---------------------------------------------------

# Assigning row names is a no-op while a value carries derivatives, and the
# labels are read from the plain evaluation at the solved values instead, so
# naming a return cannot change what exact mode computes. One family per
# ee-*.R file is checked here, together with all four of the returns whose
# extra row is a parameter of the outcome distribution rather than a
# coefficient on a design column, since those are the rows a wrong label would
# shift. test-exact-mode-acceptance.R compares exact against capprox for every
# family.

test_that("naming the rows leaves the exact-mode fits unchanged", {
  compare_methods <- function(psi, init, expected, label, solver = NULL) {
    fit <- function(method) {
      m <- MEstimator(stacked_equations = psi, init = init)
      if (is.null(solver)) {
        estimate(m, deriv_method = method)
      } else {
        estimate(m, deriv_method = method, solver = solver)
      }
    }
    exact <- fit("exact")
    approx <- fit("capprox")

    expect_equal(names(coef(exact)), expected, label = label)
    expect_equal(names(coef(approx)), expected, label = label)
    expect_equal(coef(exact), coef(approx), tolerance = 1e-8, label = label)
    expect_equal(vcov(exact), vcov(approx), tolerance = 1e-5, label = label)
  }

  y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
  compare_methods(
    function(theta) ee_mean_variance(theta, y = y),
    c(0, 1),
    c("mean", "variance"),
    "ee_mean_variance"
  )

  gf <- load_fixture("ee_gformula_exact")
  compare_methods(
    function(theta) {
      ee_gformula(theta, y = gf$y, X = gf$X, X1 = gf$X1, X0 = gf$X0)
    },
    gf$init,
    c("ACE", "E[Y^1]", "E[Y^0]", "X_1", "X_2", "X_3", "X_4"),
    "ee_gformula"
  )

  rg <- load_fixture("ee_rogan_gladen_extended_exact")
  compare_methods(
    function(theta) {
      ee_rogan_gladen_extended(
        theta,
        y = rg$y,
        y_star = rg$y_star,
        r = rg$r,
        X = rg$X
      )
    },
    rg$init,
    c("corrected_proportion", "sens_1", "sens_2", "spec_1", "spec_2"),
    "ee_rogan_gladen_extended"
  )

  surv <- survival_times()
  compare_methods(
    function(theta) {
      ee_survival_model(
        theta,
        time = surv$time,
        event = surv$event,
        distribution = "weibull"
      )
    },
    c(0.3, 1),
    c("lambda", "gamma"),
    "ee_survival_model"
  )

  d <- emax_doses()
  compare_methods(
    function(theta) ee_emax(theta, dose = d$dose, response = d$response),
    c(1, 25, 10),
    c("e0", "emax", "ed50"),
    "ee_emax"
  )

  gam <- load_fixture("ee_glm_exact_gamma_log")
  compare_methods(
    function(theta) {
      ee_glm(theta, X = gam$X, y = gam$y, distribution = "gamma", link = "log")
    },
    gam$init,
    c("X_1", "X_2", "X_3", "log_shape"),
    "ee_glm"
  )

  nb <- load_fixture("ee_glm_negative_binomial_log")
  compare_methods(
    function(theta) {
      ee_glm(
        theta,
        X = nb$X,
        y = nb$y,
        distribution = "negative_binomial",
        link = "log"
      )
    },
    nb$init,
    c("X_1", "X_2", "X_3", "log_dispersion"),
    "ee_glm (negative binomial)"
  )

  tob <- load_fixture("ee_tobit_exact")
  compare_methods(
    function(theta) {
      ee_tobit(theta, X = tob$X, y = tob$y, lower = tob$lower)
    },
    tob$init,
    c("X_1", "X_2", "X_3", "log_sigma"),
    "ee_tobit"
  )

  bet <- load_fixture("ee_beta_regression_exact")
  compare_methods(
    function(theta) ee_beta_regression(theta, X = bet$X, y = bet$y),
    bet$init,
    c("X_1", "X_2", "X_3", "log_phi"),
    "ee_beta_regression"
  )
})
