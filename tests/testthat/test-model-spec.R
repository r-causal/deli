# ---- Helpers ----

# One data frame drives most of the tests in this file. It carries a factor with
# three levels, a continuous covariate to interact with it, and a positive
# exposure to offset by, which is the combination the round-trip test needs.
model_spec_data <- function() {
  set.seed(20240712)
  n <- 60
  data.frame(
    count = rpois(n, 4) + 1,
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b", "c"), length.out = n)),
    exposure = round(stats::runif(n, 1, 5), 3),
    w = rep(c(0.5, 1.5), length.out = n)
  )
}

# Subsetting a model matrix drops the assign and contrasts attributes that
# model.matrix() puts on it, so a rebuild compared against a subset of the
# fitted design has to shed them too.
unclass_design <- function(X) {
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  X
}

model_spec_fit <- function(data = model_spec_data()) {
  m_estimate(
    count ~ x * g + offset(log(exposure)),
    data = data,
    .ee = ee_glm,
    distribution = "poisson",
    link = "log"
  )
}

# formula_model_spec ---------------------------------------------------------

test_that("formula_model_spec() records every field of the model spec", {
  m <- model_spec_fit()
  spec <- m@model_spec

  expect_type(spec, "list")
  expect_setequal(
    names(spec),
    c(
      "terms",
      "xlevels",
      "contrasts",
      "ee",
      "ee_spec_args",
      "ee_obs_args",
      "n_coef",
      "model_frame",
      "X",
      "y",
      "offset",
      "response_levels"
    )
  )
})

test_that("formula_model_spec() records the terms with response and offset", {
  m <- model_spec_fit()
  model_terms <- m@model_spec$terms

  expect_s3_class(model_terms, "terms")
  expect_identical(attr(model_terms, "response"), 1L)
  # offset(log(exposure)) is the fourth variable in the terms, after the
  # response, x, and g.
  expect_identical(attr(model_terms, "offset"), 4L)
  expect_identical(attr(model_terms, "term.labels"), c("x", "g", "x:g"))
})

test_that("formula_model_spec() records xlevels and contrasts for the factor", {
  m <- model_spec_fit()
  spec <- m@model_spec

  expect_identical(spec$xlevels, list(g = c("a", "b", "c")))
  expect_identical(spec$contrasts, list(g = "contr.treatment"))
})

test_that("formula_model_spec() records the estimating function itself", {
  m <- model_spec_fit()
  expect_identical(m@model_spec$ee, ee_glm)
})

test_that("formula_model_spec() records the design, response, and offset", {
  data <- model_spec_data()
  m <- model_spec_fit(data)
  spec <- m@model_spec

  expect_identical(dim(spec$X), c(60L, 6L))
  expect_identical(
    colnames(spec$X),
    c("(Intercept)", "x", "gb", "gc", "x:gb", "x:gc")
  )
  # The response keeps the model frame's row names.
  expect_identical(unname(spec$y), data$count)
  expect_identical(spec$offset, log(data$exposure))
})

test_that("formula_model_spec() records the design the closure was built on", {
  m <- model_spec_fit()
  # The recorded design is the same object the estimating-function closure
  # holds, not a second copy of it.
  expect_identical(m@model_spec$X, environment(m@stacked_equations)$X)
  expect_identical(m@model_spec$y, environment(m@stacked_equations)$y)
  expect_identical(
    m@model_spec$model_frame,
    environment(m@stacked_equations)$mf
  )
})

test_that("formula_model_spec() records the model frame behind the design", {
  data <- model_spec_data()
  m <- model_spec_fit(data)
  mf <- m@model_spec$model_frame

  expect_s3_class(mf, "data.frame")
  expect_identical(
    names(mf),
    c("count", "x", "g", "offset(log(exposure))")
  )
  expect_identical(nrow(mf), nrow(data))
  # The frame keeps the variables the formula named, before the factor was
  # coded and the offset lifted out, which the design no longer shows.
  expect_identical(mf$g, data$g)
  expect_identical(attr(mf, "terms"), m@model_spec$terms)
})

test_that("formula_model_spec() records n_coef as the design column count", {
  m <- model_spec_fit()
  expect_identical(m@model_spec$n_coef, 6L)
  expect_identical(m@model_spec$n_coef, ncol(m@model_spec$X))
  expect_identical(m@model_spec$n_coef, m@n_params)
})

test_that("formula_model_spec() records n_coef short of an appended parameter", {
  m <- m_estimate(
    count ~ x,
    data = model_spec_data(),
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )
  # Gamma theta is c(beta, log_shape), so the design coefficients stop one
  # short of the parameter count and X %*% theta is non-conformable.
  expect_identical(m@model_spec$n_coef, 2L)
  expect_identical(m@n_params, 3L)
  expect_identical(names(m@theta), c("(Intercept)", "x", "log_shape"))
})

test_that("formula_model_spec() records levels for a factor response", {
  data <- model_spec_data()
  data$high <- factor(
    ifelse(data$count > 4, "yes", "no"),
    levels = c("no", "yes")
  )
  m <- m_estimate(
    high ~ x,
    data = data,
    .ee = ee_regression,
    model = "logistic"
  )

  expect_identical(m@model_spec$response_levels, c("no", "yes"))
  # coerce_formula_response() scores against the first level.
  expect_identical(unname(m@model_spec$y), as.numeric(data$high == "yes"))
})

test_that("formula_model_spec() records levels for a character response", {
  data <- model_spec_data()
  data$high <- ifelse(data$count > 4, "yes", "no")
  m <- m_estimate(
    high ~ x,
    data = data,
    .ee = ee_regression,
    model = "logistic"
  )

  # `coerce_formula_response()` and `formula_response_levels()` each turn a
  # character response into a factor before reading it, so the levels reported
  # are the alphabetical ones the coercion scores against. Their docstrings
  # assert that coupling and nothing else pins it.
  expect_identical(m@model_spec$response_levels, c("no", "yes"))
  expect_identical(unname(m@model_spec$y), as.numeric(data$high == "yes"))
})

test_that("formula_model_spec() records no levels for a numeric response", {
  m <- model_spec_fit()
  expect_null(m@model_spec$response_levels)
})

test_that("formula_model_spec() records no offset when the formula has none", {
  m <- m_estimate(
    count ~ x,
    data = model_spec_data(),
    .ee = ee_regression,
    model = "poisson"
  )
  expect_null(m@model_spec$offset)
})

# spec_ee_args ---------------------------------------------------------------

test_that("spec_ee_args() keeps the model specification arguments", {
  m <- model_spec_fit()
  expect_identical(
    m@model_spec$ee_spec_args,
    list(distribution = "poisson", link = "log")
  )
})

test_that("spec_ee_args() drops the per-observation arguments", {
  data <- model_spec_data()
  m <- m_estimate(
    count ~ x + g + offset(log(exposure)),
    data = data,
    .ee = ee_glm,
    distribution = "poisson",
    link = "log",
    weights = w
  )
  spec <- m@model_spec

  # weights and the formula offset are per-observation data, so they do not
  # survive into the reusable specification.
  expect_identical(
    spec$ee_spec_args,
    list(distribution = "poisson", link = "log")
  )
  expect_false("weights" %in% names(spec$ee_spec_args))
  expect_false("offset" %in% names(spec$ee_spec_args))
  # The offset the fit used is still recorded on its own.
  expect_identical(spec$offset, log(data$exposure))
})

test_that("spec_ee_args() drops the event indicator of a survival fit", {
  set.seed(414)
  n <- 200
  x <- round(stats::rnorm(n), 3)
  event_time <- stats::rexp(n, rate = exp(-(1 + 0.5 * x)))
  censor_time <- stats::rexp(n, rate = 0.15)
  data <- data.frame(
    time = round(pmin(event_time, censor_time), 4),
    status = as.numeric(event_time <= censor_time),
    x = x
  )
  m <- m_estimate(
    time ~ x,
    data = data,
    .ee = ee_aft,
    distribution = "exponential",
    event = status
  )
  spec <- m@model_spec

  expect_identical(spec$ee_spec_args, list(distribution = "exponential"))
  expect_false("event" %in% names(spec$ee_spec_args))

  # It is recorded as one of the fit's per-observation values instead, which is
  # where a later caller predicting on the fitted sample reads it.
  expect_identical(spec$ee_obs_args, list(event = data$status))
})

test_that("spec_obs_args() keeps exactly what spec_ee_args() drops", {
  mixed <- list(
    distribution = "weibull",
    event = c(1, 0, 1),
    weights = c(1, 2, 1),
    link = "log",
    offset = c(0, 0, 0)
  )
  expect_identical(
    spec_obs_args(mixed),
    list(event = c(1, 0, 1), weights = c(1, 2, 1), offset = c(0, 0, 0))
  )
  # The two halves partition the list, so nothing forwarded is lost and nothing
  # is recorded twice.
  expect_identical(
    sort(c(names(spec_ee_args(mixed)), names(spec_obs_args(mixed)))),
    sort(names(mixed))
  )
})

test_that("spec_ee_args() records nothing when nothing forwarded is a spec", {
  data <- model_spec_data()
  m <- m_estimate(count ~ x, data = data, .ee = ee_regression, model = "linear")
  expect_identical(m@model_spec$ee_spec_args, list(model = "linear"))

  # An estimating equation that needs no arguments beyond the design and the
  # response forwards nothing at all.
  linear <- function(theta, X, y) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m2 <- m_estimate(count ~ x, data = data, .ee = linear)
  expect_type(m2@model_spec$ee_spec_args, "list")
  expect_length(m2@model_spec$ee_spec_args, 0L)
  expect_identical(m2@model_spec$ee, linear)
})

test_that("spec_ee_args() splits a mixed argument list by name", {
  expect_identical(
    spec_ee_args(list(
      distribution = "weibull",
      event = c(1, 0, 1),
      weights = c(1, 2, 1),
      link = "log",
      offset = c(0, 0, 0)
    )),
    list(distribution = "weibull", link = "log")
  )
})

# Round trip -----------------------------------------------------------------

test_that("the recorded terms, xlevels, and contrasts reproduce X exactly", {
  data <- model_spec_data()
  m <- model_spec_fit(data)
  spec <- m@model_spec

  rebuilt <- stats::model.matrix(
    spec$terms,
    stats::model.frame(spec$terms, data, xlev = spec$xlevels),
    contrasts.arg = spec$contrasts
  )

  expect_identical(rebuilt, spec$X)
})

test_that("the recorded terms reproduce a data-dependent transformation", {
  set.seed(7)
  n <- 50
  data <- data.frame(
    y = round(stats::rnorm(n, 5), 3),
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b"), length.out = n))
  )
  m <- m_estimate(
    y ~ poly(x, 2) + g,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  spec <- m@model_spec
  newdata <- data[1:5, ]

  # The terms come off the model frame, so they carry the predvars attribute
  # holding the fitted orthogonal-polynomial coefficients. That is what makes
  # the basis reproducible for data the fit never saw.
  expect_false(is.null(attr(spec$terms, "predvars")))
  rebuilt <- stats::model.matrix(
    spec$terms,
    stats::model.frame(spec$terms, newdata, xlev = spec$xlevels),
    contrasts.arg = spec$contrasts
  )
  expect_equal(
    unname(unclass_design(rebuilt)),
    unname(spec$X[1:5, , drop = FALSE])
  )

  # Re-fitting the basis from the formula alone recomputes it on whatever data
  # it is handed, which is the failure the recorded terms avoid.
  naive <- stats::model.matrix(y ~ poly(x, 2) + g, data = newdata)
  expect_false(isTRUE(all.equal(
    unname(unclass_design(naive)),
    unname(spec$X[1:5, , drop = FALSE])
  )))
})

test_that("the recorded spec rebuilds a design for data missing a level", {
  data <- model_spec_data()
  m <- model_spec_fit(data)
  spec <- m@model_spec

  # New data holding only level "a" would otherwise produce a one-column
  # factor encoding. The recorded xlevels keep all three columns in place.
  newdata <- data[data$g == "a", , drop = FALSE]
  rebuilt <- stats::model.matrix(
    spec$terms,
    stats::model.frame(spec$terms, newdata, xlev = spec$xlevels),
    contrasts.arg = spec$contrasts
  )

  # Subsetting a design drops its assign and contrasts attributes, so the
  # comparison is on the columns the rebuild produced.
  expect_identical(colnames(rebuilt), colnames(spec$X))
  expect_identical(
    unclass_design(rebuilt),
    spec$X[rownames(rebuilt), , drop = FALSE]
  )
})

test_that("the recorded terms rebuild a design from data without a response", {
  data <- model_spec_data()
  m <- model_spec_fit(data)
  spec <- m@model_spec

  newdata <- data[1:5, ]
  newdata$count <- NULL
  bare <- stats::delete.response(spec$terms)
  rebuilt <- stats::model.matrix(
    bare,
    stats::model.frame(bare, newdata, xlev = spec$xlevels),
    contrasts.arg = spec$contrasts
  )

  expect_identical(colnames(rebuilt), colnames(spec$X))
  expect_identical(
    unname(unclass_design(rebuilt)),
    unname(spec$X[1:5, , drop = FALSE])
  )
})

# GMM ------------------------------------------------------------------------

test_that("a GMM formula fit records the model spec too", {
  data <- model_spec_data()
  g <- gmm_estimate(
    count ~ x * g,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  spec <- g@model_spec

  expect_type(spec, "list")
  expect_identical(spec$ee, ee_regression)
  expect_identical(spec$ee_spec_args, list(model = "linear"))
  expect_identical(spec$n_coef, 6L)
  expect_identical(spec$xlevels, list(g = c("a", "b", "c")))

  rebuilt <- stats::model.matrix(
    spec$terms,
    stats::model.frame(spec$terms, data, xlev = spec$xlevels),
    contrasts.arg = spec$contrasts
  )
  expect_identical(rebuilt, spec$X)
})

# Function interface ---------------------------------------------------------

test_that("a .default fit records no model spec", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)

  m <- m_estimate(stacked_equations = psi, init = c(mean = 0))
  expect_null(m@model_spec)

  g <- gmm_estimate(stacked_equations = psi, init = c(mean = 0))
  expect_null(g@model_spec)
})

test_that("a hand-built estimator records no model spec", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)

  expect_null(MEstimator(stacked_equations = psi, init = 0)@model_spec)
  expect_null(GMMEstimator(stacked_equations = psi, init = 0)@model_spec)
  expect_null(
    estimate(MEstimator(stacked_equations = psi, init = 0))@model_spec
  )
})

test_that("model_spec is not a constructor argument", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)

  expect_error(
    MEstimator(stacked_equations = psi, init = 0, model_spec = list()),
    "unused argument"
  )
  expect_error(
    GMMEstimator(stacked_equations = psi, init = 0, model_spec = list()),
    "unused argument"
  )
})

# Estimates unchanged --------------------------------------------------------

test_that("recording the model spec leaves the M-estimator results alone", {
  m <- model_spec_fit()

  expect_identical(
    names(coef(m)),
    c("(Intercept)", "x", "gb", "gc", "x:gb", "x:gc")
  )
  expect_equal(
    unname(coef(m)),
    c(
      0.5653513564,
      0.0553160519,
      0.0419520443,
      0.0187500025,
      -0.0227253976,
      0.0819771366
    ),
    tolerance = 1e-8
  )
  expect_equal(
    unname(diag(vcov(m))),
    c(
      0.016669649,
      0.0333824712,
      0.0313827868,
      0.0303101337,
      0.0421783596,
      0.0495552941
    ),
    tolerance = 1e-8
  )
  expect_equal(
    unname(vcov(m)[1, ]),
    c(
      0.016669649,
      -0.0027070672,
      -0.0166696489,
      -0.0166696498,
      0.0027070672,
      0.0027070675
    ),
    tolerance = 1e-8
  )
})

test_that("recording the model spec leaves the GMM results alone", {
  g <- gmm_estimate(
    count ~ x * g,
    data = model_spec_data(),
    .ee = ee_regression,
    model = "linear"
  )

  expect_equal(
    unname(coef(g)),
    c(
      5.0507739591,
      0.2035395084,
      0.5485521177,
      -0.0792514873,
      0.0532323549,
      -0.3966075224
    ),
    tolerance = 1e-8
  )
  expect_equal(
    unname(diag(vcov(g))),
    c(
      0.1253267943,
      0.1619423422,
      0.3448506965,
      0.2934704533,
      0.2814206426,
      0.3963427768
    ),
    tolerance = 1e-8
  )
})
