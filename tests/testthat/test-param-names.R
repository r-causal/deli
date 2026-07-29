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

# The same two equations with the parameters labeled on the rows of the return.
# Row names are the naming channel a `stacked_equations` function has, since it
# is an opaque closure to the estimator and nothing else it returns is per-row.
labelled_mean_variance_psi <- function(y, labels = c("mu", "sigma2")) {
  function(theta) {
    out <- rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
    rownames(out) <- labels
    out
  }
}

# A two-parameter linear regression and a psi that labels its rows, written the
# way a user writes one. Every operation in the psi is an S3 method, so the
# return carries tangents under exact differentiation no matter which
# environment the function was written in; a `rbind()` stack would reach base
# R's binder from outside the package and lose them.
regression_fixture <- function() {
  list(
    X = base::cbind(1, c(-1.1, 0.3, 1.7, -0.4, 0.9, 2.2, -1.8, 0.5)),
    y = c(0.8, 2.1, 4.6, 1.2, 3.0, 5.4, -0.6, 2.4)
  )
}

labelled_regression_psi <- function(d, labels = c("intercept", "slope")) {
  function(theta) {
    resid <- d$y - c(d$X %*% theta)
    out <- t(d$X * resid)
    rownames(out) <- labels
    out
  }
}

# An accelerated failure time fit whose design carries column headings, so an
# `ee_*()` function that inherited them would name the parameters after them.
# The exponential distribution is the case that would supply a complete set:
# it has exactly one parameter per design column.
aft_fixture <- function() {
  list(
    X = base::cbind(intercept = 1, treated = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1)),
    time = c(2.1, 5.4, 1.2, 6.8, 3.3, 4.9, 0.8, 7.2, 2.7, 5.1),
    event = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 0)
  )
}

# A pooled logistic fit whose covariate design and time design both carry
# column headings, so between them they would supply a complete set. The
# observations come from the exact-mode reference fixture, which is large enough
# for the fit to converge from a zero start; the intercept lives in the time
# design alone, since one in each design would make the two collinear.
plogit_fixture <- function() {
  ref <- load_fixture("ee_plogit_exact")
  n <- length(ref$time)
  list(
    X = base::cbind(
      exposure = as.numeric(ref$X),
      treated = rep(c(0, 1), length.out = n)
    ),
    time = ref$time,
    event = ref$event,
    S = base::cbind(t_intercept = 1, t_linear = 1:max(ref$time))
  )
}

# A small frame with one covariate, so its model matrix has the two columns
# `(Intercept)` and `x`. prepare_formula_psi() builds the closure without
# calling it, so the same frame serves every estimating equation below whatever
# outcome that equation expects.
formula_frame <- function() {
  data.frame(
    x = c(-1.1, 0.3, 1.7, -0.4, 0.9),
    y = c(0.2, 0.5, 0.8, 0.3, 0.7)
  )
}

# The names the formula interface resolves for an `init`, read off the vector it
# hands to the estimator.
resolved_init_names <- function(.ee, init, ...) {
  prep <- prepare_formula_psi(
    y ~ x,
    data = formula_frame(),
    .ee = .ee,
    dots = rlang::quos(...),
    init = init
  )
  names(prep$init)
}

# A gamma outcome with one covariate, the smallest fit that appends a parameter
# beyond the design columns and still solves from a zero start.
formula_gamma_data <- function() {
  set.seed(6)
  d <- data.frame(x = stats::rnorm(50))
  d$y <- stats::rgamma(50, shape = 2, rate = exp(-(0.2 + 0.3 * d$x)))
  d
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

# Parameter names from estimating-function row names ---------------------------

test_that("estimate() reads parameter names from the row names of psi", {
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(y),
    init = c(0, 1)
  )

  expect_equal(names(coef(m)), c("mu", "sigma2"))
  expect_equal(
    dimnames(vcov(m)),
    list(c("mu", "sigma2"), c("mu", "sigma2"))
  )
  expect_equal(rownames(confint(m)), c("mu", "sigma2"))
  expect_equal(tidy(m)$term, c("mu", "sigma2"))
})

test_that("gmm_estimate() reads parameter names from the row names of psi", {
  y <- mean_variance_y()
  g <- gmm_estimate(
    stacked_equations = labelled_mean_variance_psi(y),
    init = c(0, 1)
  )

  expect_equal(names(coef(g)), c("mu", "sigma2"))
  expect_equal(
    dimnames(vcov(g)),
    list(c("mu", "sigma2"), c("mu", "sigma2"))
  )
})

test_that("names on init override the row names of psi", {
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(y),
    init = c(a = 0, b = 1)
  )

  expect_equal(names(coef(m)), c("a", "b"))
})

test_that("a partially named init closes the row-name channel", {
  # The two channels are never merged. Once `init` carries any name at all it
  # is the only one read, and its unnamed entries are numbered by position, so
  # a reader can tell which channel produced a label from `init` alone.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(y),
    init = c(a = 0, 1)
  )

  expect_equal(names(coef(m)), c("a", "theta_2"))
})

test_that("a partially labelled psi stack is discarded rather than filled", {
  # `rbind()` pads the rows of an unlabeled block with empty strings, so a
  # stack that labels only some of its blocks returns an incomplete vector.
  # Unlike a partially named `init`, that is rarely deliberate, so the whole
  # vector is dropped and every parameter is numbered.
  y <- mean_variance_y()
  psi <- function(theta) {
    labelled <- matrix(y - theta[1], nrow = 1, dimnames = list("mu", NULL))
    unlabelled <- matrix((y - theta[1])^2 - theta[2], nrow = 1)
    rbind(labelled, unlabelled)
  }

  expect_equal(rownames(psi(c(0, 1))), c("mu", ""))

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("theta_1", "theta_2"))
})

test_that("base rbind's variable-name labels reach the parameter names", {
  # base::rbind() labels a row with the name of the variable supplying it
  # whenever that variable is a plain vector, so a stack of vectors written
  # outside the package arrives already labeled and those labels are read.
  # base:: is named explicitly because deli's own rbind() is masked inside the
  # namespace and forwards deparse.level = 0, which suppresses exactly these
  # labels; a test written against the masked binder would assert nothing about
  # what a user of the installed package gets.
  y <- mean_variance_y()
  psi <- function(theta) {
    mu <- y - theta[1]
    sigma2 <- (y - theta[1])^2 - theta[2]
    base::rbind(mu, sigma2)
  }

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("mu", "sigma2"))
})

test_that("a matrix argument to base rbind carries no variable-name label", {
  # The labels come only from vector arguments. A stack mixing a labeled vector
  # with a matrix is therefore incomplete and is discarded, unless the matrix
  # brought row names of its own. That is what keeps the unlabeled built-in
  # equations, the regression family among them, out of the naming channel: each
  # returns a matrix, and a matrix contributes only the row names it was given.
  # The equations that do write their own are covered in
  # test-param-names-builtin.R.
  y <- mean_variance_y()
  psi <- function(theta) {
    mu <- y - theta[1]
    sigma2 <- base::matrix((y - theta[1])^2 - theta[2], nrow = 1)
    base::rbind(mu, sigma2)
  }

  expect_equal(rownames(psi(c(0, 1))), c("mu", ""))

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("theta_1", "theta_2"))
})

test_that("an over-identified GMM fit numbers its parameters", {
  # Two moment conditions label the rows against a single parameter, so the row
  # names describe the equations rather than the parameters and are not read.
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    out <- rbind(
      y - theta[1],
      (y - theta[1])^2 - 2
    )
    rownames(out) <- c("mean_eq", "var_eq")
    out
  }

  g <- gmm_estimate(stacked_equations = psi, init = 2)

  expect_equal(names(coef(g)), "theta_1")
  # Row names on psi reached @meat and @weight_matrix before this naming
  # channel existed, carried there by tcrossprod() and then by solve(). Neither
  # matrix is indexed by the parameters in an over-identified fit, so both are
  # now unlabeled, as @bread already was.
  expect_null(dimnames(g@meat))
  expect_null(dimnames(g@weight_matrix))
  expect_null(dimnames(g@bread))
})

test_that("a row-labelled psi gives the same fit under exact differentiation", {
  # `rownames<-` is an ordinary closure in base R and cannot take a method, so
  # deli registers its no-op on `dimnames<-`, which `rownames<-` delegates to
  # and which is a generic. That dispatch is what lets an estimating function
  # written outside the package label its rows under exact mode, and the
  # assertion below is what keeps this test honest about it: testthat's test
  # environment is a clone of the package namespace, so a mask on the setter
  # would stand in for the dispatch and the test would pass while asserting
  # nothing about what a user of the installed package gets.
  expect_identical(`rownames<-`, base::`rownames<-`)

  d <- regression_fixture()
  psi <- labelled_regression_psi(d)

  expect_equal(rownames(psi(c(0, 0))), c("intercept", "slope"))

  approx <- m_estimate(stacked_equations = psi, init = c(0, 0))
  exact <- m_estimate(
    stacked_equations = psi,
    init = c(0, 0),
    deriv_method = "exact"
  )

  expect_equal(names(coef(exact)), c("intercept", "slope"))
  expect_equal(coef(exact), coef(approx), tolerance = 1e-6)
  expect_equal(vcov(exact), vcov(approx), tolerance = 1e-6)
})

test_that("a psi that repeats a row name numbers the parameters", {
  # A repeated label is worse than no label: coef() would show one name twice
  # and confint(m)["mu", ] would silently return only the first row. The stack
  # that produces one need not look unusual, since rbind() names each row after
  # the variable that supplied it, so it is discarded whole as an incomplete
  # set is.
  y <- mean_variance_y()
  psi <- labelled_mean_variance_psi(y, labels = c("mu", "mu"))

  expect_equal(rownames(psi(c(0, 1))), c("mu", "mu"))

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("theta_1", "theta_2"))
})

test_that("names the caller typed on init may repeat", {
  # The rule above belongs to the row-name channel alone. A caller who writes
  # the repetition out on `init` has said what they meant, and that has always
  # been allowed, so it stays allowed.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = mean_variance_psi(y),
    init = c(mu = 0, mu = 1)
  )

  expect_equal(names(coef(m)), c("mu", "mu"))
})

# Parameter names on the formula interface -------------------------------------

test_that("an unnamed init takes the model matrix column names", {
  # The formula interface knows what the caller called the design columns, which
  # coerce_design() drops before any estimating equation sees them, so no row
  # name can carry them.
  expect_equal(
    resolved_init_names(ee_regression, c(0, 0), model = "linear"),
    c("(Intercept)", "x")
  )
})

test_that("an unnamed init takes the parameter its equation appends", {
  expect_equal(
    resolved_init_names(
      ee_glm,
      rep(0, 3),
      distribution = "gamma",
      link = "log"
    ),
    c("(Intercept)", "x", "log_shape")
  )
  expect_equal(
    resolved_init_names(
      ee_glm,
      rep(0, 3),
      distribution = "negative_binomial",
      link = "log"
    ),
    c("(Intercept)", "x", "log_dispersion")
  )
  expect_equal(
    resolved_init_names(ee_glm, rep(0, 3), distribution = "nb", link = "log"),
    c("(Intercept)", "x", "log_dispersion")
  )
  expect_equal(
    resolved_init_names(ee_tobit, rep(0, 3), lower = 0),
    c("(Intercept)", "x", "log_sigma")
  )
  expect_equal(
    resolved_init_names(ee_beta_regression, rep(0, 3)),
    c("(Intercept)", "x", "log_phi")
  )
  expect_equal(
    resolved_init_names(
      ee_aft,
      rep(0, 3),
      event = 1,
      distribution = "weibull"
    ),
    c("(Intercept)", "x", "log_inv_scale")
  )
})

test_that("an equation that appends nothing names only the design columns", {
  # The tweedie power is a fixed hyperparameter and the exponential AFT scale is
  # fixed at 1, so both estimate one parameter per design column and neither
  # accounts for a third.
  expect_equal(
    resolved_init_names(
      ee_glm,
      c(0, 0),
      distribution = "tweedie",
      link = "log",
      hyperparameter = 1.5
    ),
    c("(Intercept)", "x")
  )
  expect_null(
    resolved_init_names(
      ee_glm,
      rep(0, 3),
      distribution = "tweedie",
      link = "log",
      hyperparameter = 1.5
    )
  )
  expect_equal(
    resolved_init_names(
      ee_aft,
      c(0, 0),
      event = 1,
      distribution = "exponential"
    ),
    c("(Intercept)", "x")
  )
  expect_null(
    resolved_init_names(
      ee_aft,
      rep(0, 3),
      event = 1,
      distribution = "exponential"
    )
  )
})

test_that("an init of a length the design cannot account for is left unnamed", {
  # Naming here would attach the intercept's label to whatever the first
  # parameter happens to be, so the vector is passed on as the caller wrote it
  # and the row names of the estimating functions label the fit instead, or the
  # parameters are numbered where those row names do not label every one.
  expect_null(resolved_init_names(ee_regression, rep(0, 3), model = "linear"))
  expect_null(resolved_init_names(ee_regression, 0, model = "linear"))
  expect_null(
    resolved_init_names(ee_glm, rep(0, 4), distribution = "gamma", link = "log")
  )
})

test_that("an unrecognized equation names only the design columns", {
  # An equation this package does not know may append anything or nothing, so
  # only the length that needs no such knowledge is named.
  ee_custom <- function(theta, X, y, ...) t(X * as.vector(y - X %*% theta))

  expect_equal(
    resolved_init_names(ee_custom, c(0, 0)),
    c("(Intercept)", "x")
  )
  expect_null(resolved_init_names(ee_custom, rep(0, 3)))
})

test_that("an init the caller named reaches the estimator unaltered", {
  expect_equal(
    resolved_init_names(ee_regression, c(a = 0, b = 0), model = "linear"),
    c("a", "b")
  )
  # A partially named `init` is left alone as well, so its gaps are filled by
  # position in the one place that fills them.
  expect_equal(
    resolved_init_names(
      ee_glm,
      c(a = 0, 0, 0),
      distribution = "gamma",
      link = "log"
    ),
    c("a", "", "")
  )
})

test_that("a formula fit labels an unnamed init through to coef()", {
  d <- formula_gamma_data()
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )
  expected <- c("(Intercept)", "x", "log_shape")

  expect_equal(names(coef(m)), expected)
  expect_equal(dimnames(vcov(m)), list(expected, expected))
  expect_equal(rownames(confint(m)), expected)
  expect_equal(tidy(m)$term, expected)
})

test_that("labeling an unnamed init leaves the estimates unchanged", {
  # Names never enter the arithmetic, so the fit is the one the same call made
  # before the labels were resolved.
  d <- formula_gamma_data()
  fit <- function(init) {
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_glm,
      distribution = "gamma",
      link = "log",
      init = init
    )
  }
  unnamed <- fit(c(0, 0, 0))
  named <- fit(c(b0 = 0, b1 = 0, shape = 0))

  expect_identical(unname(coef(unnamed)), unname(coef(named)))
  expect_identical(unname(vcov(unnamed)), unname(vcov(named)))
  expect_equal(names(coef(named)), c("b0", "b1", "shape"))
})

test_that("an init the design cannot account for fails in the equation", {
  # Assigning the design column names to a vector of another length would raise
  # a names-length error before the fit began, reporting a problem with the
  # wrong thing. The failure has to stay the one the estimating equation
  # raises.
  d <- formula_gamma_data()

  expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(0, 0, 0)
    ),
    "non-conformable"
  )
})

test_that("an equation whose extra parameter leads keeps its own row names", {
  # ee_gformula() holds one parameter more than the design has columns, the same
  # count as ee_tobit(), but the extra one leads rather than trails. That is why
  # the appended parameter is written out for each equation rather than inferred
  # from the count: a fit labeled by count alone would call the causal mean
  # `(Intercept)` and shift every coefficient's label by one.
  set.seed(31)
  n <- 200
  d <- data.frame(x = stats::rbinom(n, 1, 0.5))
  d$y <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.8 * d$x))
  X1 <- cbind(1, rep(1, n))

  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_gformula,
    X1 = X1,
    init = c(0.5, 0, 0)
  )

  expect_equal(names(coef(m)), c("causal_mean", "X_1", "X_2"))
})

test_that("a GMM formula fit labels an unnamed init as well", {
  d <- formula_gamma_data()
  # The plain residual stacked beneath the two score rows leaves the two
  # coefficients over-identified, which is the case GMM is for.
  ee_over <- function(theta, X, y, ...) {
    r <- as.vector(y - X %*% theta)
    base::rbind(t(X * r), r)
  }

  g <- gmm_estimate(y ~ x, data = d, .ee = ee_over, init = c(0, 0))

  expect_equal(names(coef(g)), c("(Intercept)", "x"))
})

# ee_*() functions and the naming channel --------------------------------------

test_that("ee_aft() does not name the parameters from its design columns", {
  # Every built-in estimating equation routes its design through
  # coerce_design(), which drops the dimnames, so the caller's column headings
  # cannot reach the row names of the return. ee_aft() under the exponential
  # distribution is the case that would otherwise supply a complete set: it has
  # one parameter per design column and builds its return as t(X * ...). Sigma
  # is fixed at 1 there, so every parameter is a design coefficient and the
  # return is left unnamed. The other distributions carry log(1/sigma) as well
  # and do name their rows; test-param-names-builtin.R covers that.
  d <- aft_fixture()
  psi <- function(theta) {
    ee_aft(
      theta,
      X = d$X,
      time = d$time,
      event = d$event,
      distribution = "exponential"
    )
  }

  expect_equal(colnames(d$X), c("intercept", "treated"))
  expect_null(rownames(psi(c(0, 0))))

  m <- m_estimate(stacked_equations = psi, init = c(0, 0))

  expect_equal(names(coef(m)), c("theta_1", "theta_2"))
})

test_that("ee_plogit() does not name the parameters from either design", {
  # ee_plogit() stacks a covariate design and a time design, and both used to
  # contribute their column headings, so a named `S` completed the set that a
  # named `X` started.
  d <- plogit_fixture()
  psi <- function(theta) {
    ee_plogit(
      theta,
      X = d$X,
      time = d$time,
      event = d$event,
      S = d$S
    )
  }

  expect_equal(colnames(d$X), c("exposure", "treated"))
  expect_equal(colnames(d$S), c("t_intercept", "t_linear"))
  expect_null(rownames(psi(rep(0, 4))))

  m <- m_estimate(stacked_equations = psi, init = rep(0, 4))

  expect_equal(names(coef(m)), paste0("theta_", 1:4))
})

test_that("a caller still names an ee_*() fit through init", {
  # Closing the design-column route leaves the deliberate channels intact.
  d <- aft_fixture()
  psi <- function(theta) {
    ee_aft(
      theta,
      X = d$X,
      time = d$time,
      event = d$event,
      distribution = "exponential"
    )
  }

  m <- m_estimate(
    stacked_equations = psi,
    init = c(intercept = 0, treated = 0)
  )

  expect_equal(names(coef(m)), c("intercept", "treated"))
})

test_that("a clustered fit keeps the row names through aggregate_efuncs()", {
  set.seed(42)
  group <- rep(1:20, each = 4)
  y <- rnorm(80, mean = 3)
  psi <- function(theta) {
    out <- rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
    rownames(out) <- c("mu", "sigma2")
    aggregate_efuncs(out, group = group)
  }

  m <- m_estimate(stacked_equations = psi, init = c(0, 1))

  expect_equal(names(coef(m)), c("mu", "sigma2"))
  expect_equal(m@n_obs, 20L)
})

test_that("a brace in a parameter name is escaped rather than interpolated", {
  # cli reads a brace in a label as the start of an interpolated expression.
  # Some brace-bearing names error there, but several corrupt silently:
  # "theta{1}" prints as "theta1" and "E[Y^{a=1}]" as "E[Y^1: ]". The names
  # themselves are kept exactly as the estimating function wrote them.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(
      y,
      labels = c("theta{1}", "E[Y^{a=1}]")
    ),
    init = c(0, 1)
  )

  expect_equal(names(coef(m)), c("theta{1}", "E[Y^{a=1}]"))

  output <- cli::cli_fmt(print(m))
  expect_true(any(grepl("theta{1}", output, fixed = TRUE)))
  expect_true(any(grepl("E[Y^{a=1}]", output, fixed = TRUE)))
})

test_that("a parameter name cli would reject prints literally", {
  # "{foo}" and "beta{" are the brace-bearing names that do raise a cli error,
  # so escaping has to reach them as well as the ones that corrupt quietly.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(
      y,
      labels = c("{foo}", "beta{")
    ),
    init = c(0, 1)
  )

  output <- expect_no_error(cli::cli_fmt(print(m)))
  expect_true(any(grepl("{foo}", output, fixed = TRUE)))
  expect_true(any(grepl("beta{", output, fixed = TRUE)))
})

test_that("printing an estimator with no parameter names still fails loudly", {
  # A fitted estimator always names its parameters, so an unnamed @theta is
  # reachable only by assigning one. The escape step has to pass NULL straight
  # back: gsub() answers it with character(0), and stats::setNames() pads a
  # short name vector with NA, so an escaped NULL would relabel every
  # coefficient NA instead of leaving cli to refuse an unnamed list.
  y <- mean_variance_y()
  m <- m_estimate(stacked_equations = mean_variance_psi(y), init = c(0, 1))
  m@theta <- unname(m@theta)

  expect_null(escape_cli_braces(NULL))
  expect_error(cli::cli_fmt(print(m)), "must be a named character vector")
})

test_that("a summary prints a brace-bearing parameter name unaltered", {
  # The summary table is written with cli_verbatim(), which does not
  # interpolate, so its rows are safe as they stand. Pinned so a later switch
  # to an interpolating cli function has to account for the names.
  y <- mean_variance_y()
  m <- m_estimate(
    stacked_equations = labelled_mean_variance_psi(
      y,
      labels = c("theta{1}", "E[Y^{a=1}]")
    ),
    init = c(0, 1)
  )

  output <- cli::cli_fmt(print(summary(m)))

  expect_true(any(grepl("theta{1}", output, fixed = TRUE)))
  expect_true(any(grepl("E[Y^{a=1}]", output, fixed = TRUE)))
})
