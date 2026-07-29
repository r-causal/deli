# ---- Helpers ----

# The same shape of data the other prediction tests use: a continuous
# covariate, a three-level factor, and a continuous, a binary, and a count
# response.
augment_data <- function() {
  set.seed(20240728)
  n <- 50
  d <- data.frame(
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b", "c"), length.out = n))
  )
  d$y <- round(1 + 0.7 * d$x + as.numeric(d$g) + stats::rnorm(n), 3)
  d$hit <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x))
  d
}

augment_fit <- function(data = augment_data()) {
  m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
}

flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

# The columns augment() appends -----------------------------------------------

test_that("augment() returns the model frame beside the fitted columns", {
  data <- augment_data()
  m <- augment_fit(data)
  a <- generics::augment(m)

  expect_s3_class(a, "data.frame")
  expect_identical(class(a), "data.frame")
  expect_identical(
    names(a),
    c("y", "x", "g", ".fitted", ".se.fit", ".lower", ".upper", ".resid")
  )
  expect_identical(a$y, data$y)
  expect_identical(a$g, data$g)
})

test_that("augment() returns one row per observation the fit used", {
  m <- augment_fit()
  expect_identical(nrow(generics::augment(m)), nobs(m))
})

test_that("augment() takes .fitted, its interval, and its error from predict()", {
  m <- augment_fit()
  a <- generics::augment(m)
  p <- predict(m, se.fit = TRUE, interval = "confidence")

  expect_identical(a$.fitted, unname(p$fit[, "fit"]))
  expect_identical(a$.se.fit, unname(p$se.fit))
  expect_identical(a$.lower, unname(p$fit[, "lwr"]))
  expect_identical(a$.upper, unname(p$fit[, "upr"]))
  expect_identical(a$.fitted, unname(predict(m)))
})

test_that("augment() takes .resid from residuals()", {
  m <- augment_fit()
  expect_identical(generics::augment(m)$.resid, unname(residuals(m)))
})

test_that("augment() returns plain columns rather than a model frame", {
  m <- augment_fit()
  a <- generics::augment(m)

  expect_null(attr(a, "terms"))
  expect_null(attr(a, "na.action"))
  expect_null(names(a$.fitted))
})

# The prediction scale ---------------------------------------------------------

test_that("augment() puts .fitted on the link scale by default", {
  data <- augment_data()
  m <- m_estimate(hit ~ x, data = data, .ee = ee_regression, model = "logistic")
  a <- generics::augment(m)

  expect_identical(a$.fitted, unname(predict(m, type = "link")))
  expect_true(any(a$.fitted < 0))
})

test_that("augment(type.predict = 'response') puts .fitted on the mean scale", {
  data <- augment_data()
  m <- m_estimate(hit ~ x, data = data, .ee = ee_regression, model = "logistic")
  a <- generics::augment(m, type.predict = "response")

  expect_identical(a$.fitted, unname(predict(m, type = "response")))
  expect_true(all(a$.fitted > 0 & a$.fitted < 1))
  expect_identical(
    a$.lower,
    unname(predict(m, type = "response", interval = "confidence")[, "lwr"])
  )
  expect_identical(
    a$.se.fit,
    unname(
      predict(m, type = "response", se.fit = TRUE)$se.fit
    )
  )
})

test_that("augment() leaves .resid on the response scale at either scale", {
  data <- augment_data()
  m <- m_estimate(hit ~ x, data = data, .ee = ee_regression, model = "logistic")

  link <- generics::augment(m)
  response <- generics::augment(m, type.predict = "response")

  # The residual is the response minus the conditional mean, which is what
  # residuals() answers, whichever scale .fitted is asked for on.
  expect_identical(link$.resid, response$.resid)
  expect_identical(link$.resid, unname(residuals(m)))
  # On the response scale, and only there, .resid is the response minus
  # .fitted.
  expect_equal(response$.resid, data$hit - response$.fitted)
  expect_false(isTRUE(all.equal(link$.resid, data$hit - link$.fitted)))
})

test_that("augment() rejects a prediction scale it does not offer", {
  m <- augment_fit()
  expect_error(generics::augment(m, type.predict = "terms"), "arg")
})

# The confidence level ---------------------------------------------------------

test_that("augment() takes its interval at conf.level", {
  m <- augment_fit()
  a <- generics::augment(m, conf.level = 0.8)
  p <- predict(m, interval = "confidence", level = 0.8)

  expect_identical(a$.lower, unname(p[, "lwr"]))
  expect_identical(a$.upper, unname(p[, "upr"]))
})

test_that("augment() widens its interval as conf.level rises", {
  m <- augment_fit()
  narrow <- generics::augment(m, conf.level = 0.8)
  wide <- generics::augment(m, conf.level = 0.99)

  expect_true(all(wide$.lower < narrow$.lower))
  expect_true(all(wide$.upper > narrow$.upper))
})

test_that("augment() rejects a conf.level outside the unit interval", {
  m <- augment_fit()
  # The argument the message names is the one the caller wrote, rather than
  # `level`, which is what predict() calls it.
  expect_error(generics::augment(m, conf.level = 0), "conf.level")
  expect_error(generics::augment(m, conf.level = c(0.9, 0.95)), "conf.level")
})

# newdata ----------------------------------------------------------------------

test_that("augment() on newdata returns newdata beside the fitted columns", {
  data <- augment_data()
  m <- augment_fit(data)
  newdata <- data.frame(
    x = c(-1, 0, 1),
    g = factor("b", levels = levels(data$g))
  )
  a <- generics::augment(m, newdata = newdata)

  expect_identical(
    names(a),
    c("x", "g", ".fitted", ".se.fit", ".lower", ".upper")
  )
  expect_identical(nrow(a), 3L)
  expect_identical(a$x, newdata$x)
  expect_identical(a$.fitted, unname(predict(m, newdata = newdata)))
})

test_that("augment() on newdata has no .resid to report", {
  data <- augment_data()
  m <- augment_fit(data)
  a <- generics::augment(m, newdata = data[1:5, , drop = FALSE])

  # There is a response in these rows, but newdata in general carries none, so
  # the column is absent rather than present for some data and not others.
  expect_false(".resid" %in% names(a))
  expect_true(".fitted" %in% names(a))
})

test_that("augment() on newdata keeps every column newdata carried", {
  data <- augment_data()
  m <- augment_fit(data)
  newdata <- data[1:5, , drop = FALSE]
  a <- generics::augment(m, newdata = newdata)

  expect_identical(a[names(newdata)], newdata)
  expect_identical(rownames(a), rownames(newdata))
})

test_that("augment() keeps a newdata row with a missing value", {
  data <- augment_data()
  m <- augment_fit(data)
  newdata <- data.frame(
    x = c(-1, NA, 1),
    g = factor(c("a", "b", "c"), levels = levels(data$g))
  )
  a <- generics::augment(m, newdata = newdata)

  expect_identical(nrow(a), 3L)
  expect_identical(is.na(a$.fitted), c(FALSE, TRUE, FALSE))
  expect_identical(is.na(a$.se.fit), c(FALSE, TRUE, FALSE))
})

test_that("augment() requires newdata to be a data frame", {
  m <- augment_fit()
  expect_error(
    generics::augment(m, newdata = list(x = 1, g = "a")),
    "data frame"
  )
})

# Missing values in the fitted data --------------------------------------------

test_that("augment() reports the complete rows the fit was solved on", {
  data <- augment_data()
  data$x[c(3, 9)] <- NA
  m <- m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
  a <- generics::augment(m)

  expect_identical(nrow(a), nobs(m))
  expect_identical(nrow(a), nrow(data) - 2L)
  expect_identical(rownames(a), rownames(stats::model.frame(m)))
  expect_identical(a$.fitted, unname(predict(m)))
})

# GMM --------------------------------------------------------------------------

test_that("augment() answers for a GMM formula fit", {
  data <- augment_data()
  g <- gmm_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  a <- generics::augment(g)

  expect_identical(
    names(a),
    c("y", "x", "g", ".fitted", ".se.fit", ".lower", ".upper", ".resid")
  )
  expect_identical(nrow(a), nobs(g))
  expect_identical(a$.fitted, unname(predict(g)))
})

# Fits augment() cannot answer for ---------------------------------------------

test_that("augment() aborts for a fit built from a stacked_equations function", {
  y <- mtcars$mpg
  X <- cbind(1, mtcars$wt)
  m <- m_estimate(
    stacked_equations = function(theta) {
      ee_regression(theta, X = X, y = y, model = "linear")
    },
    init = c(0, 0)
  )

  err <- expect_error(generics::augment(m), "formula interface")
  flat <- flatten_message(err)
  expect_match(flat, "augment()", fixed = TRUE)
  expect_match(flat, "regression_predictions()", fixed = TRUE)
})

test_that("augment() aborts for an equation predict() does not support", {
  set.seed(414)
  n <- 150
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

  # The abort names augment() rather than predict(), which it reaches the
  # supported-equation table through.
  err <- expect_error(generics::augment(m), "linear predictor")
  expect_match(flatten_message(err), "augment()", fixed = TRUE)
})

test_that("augment() aborts before an estimator has been estimated", {
  y <- c(1, 2, 3, 4, 5)
  m <- MEstimator(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(0)
  )
  expect_error(generics::augment(m), "before calling")
})

# Argument handling ------------------------------------------------------------

test_that("augment() rejects an argument it does not take", {
  data <- augment_data()
  m <- augment_fit(data)

  # A name augment() does not take would otherwise be discarded: `new_data`
  # would augment the fitted rows while the caller believed they had asked for
  # rows of their own.
  err <- expect_error(
    generics::augment(m, new_data = data),
    class = "rlib_error_dots_nonempty"
  )
  expect_match(flatten_message(err), "new_data", fixed = TRUE)
})
