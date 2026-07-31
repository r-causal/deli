# ---- Helpers ----

# One data frame drives most of the tests in this file: a continuous covariate,
# a three-level factor so the model frame carries something that is not numeric,
# and a continuous, a binary, and a count response.
accessor_data <- function() {
  set.seed(20240726)
  n <- 60
  d <- data.frame(
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b", "c"), length.out = n))
  )
  d$y <- round(1 + 0.7 * d$x + as.numeric(d$g) + stats::rnorm(n), 3)
  d$hit <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x))
  d$count <- stats::rpois(n, exp(0.4 + 0.2 * d$x))
  d
}

accessor_fit <- function(data = accessor_data()) {
  m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
}

# Three rows of covariate values, with a response, since a model frame holds
# every variable the formula names.
accessor_newdata <- function() {
  data.frame(
    y = c(0, 0, 0),
    x = c(-1, 0, 1),
    g = factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  )
}

# A fit built from a stacked_equations function, which records no model
# specification for any of these accessors to report.
function_fit <- function() {
  y <- mtcars$mpg
  X <- cbind(1, mtcars$wt)
  m_estimate(
    stacked_equations = function(theta) {
      ee_regression(theta, X = X, y = y, model = "linear")
    },
    init = c(0, 0)
  )
}

# A formula fit of an equation predict() does not support: the linear predictor
# of an accelerated failure time model is on the log-time scale, so it is not a
# link-scale mean of the response.
aft_fit <- function() {
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
  m_estimate(
    time ~ x,
    data = data,
    .ee = ee_aft,
    distribution = "exponential",
    event = status
  )
}

flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

# fitted -----------------------------------------------------------------------

test_that("fitted() returns the response-scale predictions of the fitted rows", {
  data <- accessor_data()
  m <- accessor_fit(data)

  expect_identical(fitted(m), predict(m, type = "response"))
  expect_length(fitted(m), nobs(m))
  expect_identical(names(fitted(m)), rownames(m@model_spec$X))
})

test_that("fitted() of a linear fit matches lm()", {
  data <- accessor_data()
  m <- accessor_fit(data)
  reference <- stats::lm(y ~ x + g, data = data)

  expect_equal(fitted(m), stats::fitted(reference), tolerance = 1e-6)
})

test_that("fitted() of a logistic fit is on the probability scale", {
  data <- accessor_data()
  m <- m_estimate(hit ~ x, data = data, .ee = ee_regression, model = "logistic")

  expect_equal(fitted(m), stats::plogis(predict(m)))
  expect_true(all(fitted(m) > 0 & fitted(m) < 1))
  # The response scale is the one that differs from the linear predictor, which
  # is what predict() answers by default.
  expect_false(isTRUE(all.equal(unname(fitted(m)), unname(predict(m)))))
})

test_that("fitted() of a poisson glm fit is on the count scale", {
  data <- accessor_data()
  m <- m_estimate(
    count ~ x,
    data = data,
    .ee = ee_glm,
    distribution = "poisson",
    link = "log"
  )
  expect_equal(fitted(m), exp(predict(m)))
})

test_that("fitted.values() reaches the same method", {
  m <- accessor_fit()
  expect_identical(stats::fitted.values(m), fitted(m))
})

# residuals --------------------------------------------------------------------

test_that("residuals() is the response minus the fitted values", {
  data <- accessor_data()
  m <- accessor_fit(data)

  expect_equal(unname(residuals(m)), data$y - unname(fitted(m)))
  expect_identical(names(residuals(m)), names(fitted(m)))
  expect_length(residuals(m), nobs(m))
})

test_that("residuals() of a linear fit matches lm()", {
  data <- accessor_data()
  m <- accessor_fit(data)
  reference <- stats::lm(y ~ x + g, data = data)

  expect_equal(residuals(m), stats::residuals(reference), tolerance = 1e-6)
})

test_that("residuals() of a logistic fit is on the response scale", {
  data <- accessor_data()
  m <- m_estimate(hit ~ x, data = data, .ee = ee_regression, model = "logistic")

  # The residual is the response minus the fitted probability, which is the
  # conditional mean of the response. Subtracting the linear predictor instead
  # would residualize against log-odds, which the response is not measured in.
  expect_equal(
    unname(residuals(m)),
    data$hit -
      unname(stats::plogis(
        predict(m)
      ))
  )
  expect_false(isTRUE(all.equal(
    unname(residuals(m)),
    data$hit - unname(predict(m))
  )))
})

test_that("residuals() of a poisson glm fit is on the count scale", {
  data <- accessor_data()
  m <- m_estimate(
    count ~ x,
    data = data,
    .ee = ee_glm,
    distribution = "poisson",
    link = "log"
  )
  expect_equal(unname(residuals(m)), data$count - unname(exp(predict(m))))
})

test_that("residuals() of a fit with a factor response scores the same way", {
  data <- accessor_data()
  data$high <- factor(ifelse(data$hit == 1, "yes", "no"), c("no", "yes"))
  m <- m_estimate(
    high ~ x,
    data = data,
    .ee = ee_regression,
    model = "logistic"
  )

  # coerce_formula_response() scores the response against its first level, so
  # the residual is measured against the same 0/1 coding the fit solved for.
  expect_equal(
    unname(residuals(m)),
    as.numeric(data$high == "yes") - unname(fitted(m))
  )
})

test_that("resid() reaches the same method", {
  m <- accessor_fit()
  expect_identical(stats::resid(m), residuals(m))
})

test_that("residuals() takes the one type it computes", {
  m <- accessor_fit()
  expect_identical(residuals(m, type = "response"), residuals(m))
})

test_that("residuals() aborts on a type it does not compute", {
  m <- accessor_fit()

  # A caller migrating from glm() writes type = reflexively, and would
  # otherwise be handed the response residual under the name they asked for.
  for (type in c("pearson", "deviance", "working", "partial")) {
    err <- expect_error(residuals(m, type = type), "type")
    expect_match(flatten_message(err), "response residual", fixed = TRUE)
  }
  expect_match(
    flatten_message(expect_error(residuals(m, type = "pearson"))),
    "variance function",
    fixed = TRUE
  )
  expect_match(
    flatten_message(expect_error(residuals(m, type = "deviance"))),
    "likelihood",
    fixed = TRUE
  )
})

test_that("the refused residual type names the call the caller wrote", {
  # The guard beside it, `check_dots_empty()`, already reports the
  # `residuals()` call, and this one reported no call at all, so one message
  # about the same call carried a frame and the other did not.
  m <- accessor_fit()
  err <- expect_error(residuals(m, type = "pearson"))
  expect_identical(conditionCall(err), quote(residuals(m, type = "pearson")))
})

# model.frame ------------------------------------------------------------------

test_that("model.frame() returns the frame the fit was built from", {
  data <- accessor_data()
  m <- accessor_fit(data)
  mf <- stats::model.frame(m)

  expect_s3_class(mf, "data.frame")
  expect_identical(names(mf), c("y", "x", "g"))
  expect_identical(nrow(mf), nobs(m))
  expect_identical(mf$y, data$y)
  expect_identical(mf$g, data$g)
})

test_that("model.frame() carries the terms the design rebuilds from", {
  m <- accessor_fit()
  rebuilt <- stats::model.matrix(
    stats::terms(m),
    stats::model.frame(m),
    contrasts.arg = m@model_spec$contrasts
  )
  expect_identical(rebuilt, m@model_spec$X)
})

test_that("model.matrix() answers with the design the fit was solved on", {
  m <- accessor_fit()

  # With no data to build one from, the recorded design is reported as recorded,
  # which under the default contrasts is also what the terms and the model frame
  # above rebuild to. The assertion is here so that the three stay the same
  # matrix.
  expect_identical(stats::model.matrix(m), m@model_spec$X)
})

test_that("model.frame() keeps the offset term of the formula", {
  data <- accessor_data()
  data$exposure <- round(stats::runif(nrow(data), 1, 5), 3)
  m <- m_estimate(
    count ~ x + offset(log(exposure)),
    data = data,
    .ee = ee_regression,
    model = "poisson"
  )
  mf <- stats::model.frame(m)

  expect_identical(names(mf), c("count", "x", "offset(log(exposure))"))
  expect_identical(stats::model.offset(mf), log(data$exposure))
})

# model.frame(data = ) ---------------------------------------------------------

test_that("model.frame() builds the frame of the data it is given", {
  m <- accessor_fit()
  nd <- accessor_newdata()
  mf <- stats::model.frame(m, data = nd)

  expect_identical(nrow(mf), nrow(nd))
  expect_identical(names(mf), c("y", "x", "g"))
  expect_identical(mf$x, nd$x)
  # The fitted frame is the answer only when no data is given, so the two
  # cannot be confused for one another.
  expect_false(identical(nrow(mf), nrow(stats::model.frame(m))))
})

test_that("model.matrix() honors the data passed through to model.frame()", {
  m <- accessor_fit()
  nd <- accessor_newdata()
  mm <- stats::model.matrix(m, data = nd)

  expect_identical(nrow(mm), nrow(nd))
  expect_identical(colnames(mm), colnames(m@model_spec$X))
  expect_false(identical(dim(mm), dim(m@model_spec$X)))
  # The design at the caller's covariate values is the one predict() forms, so
  # multiplying it by the coefficients gives the predictions.
  expect_equal(
    unname(drop(mm %*% coef(m))),
    unname(predict(m, newdata = nd))
  )
})

test_that("model.frame() keeps the fitted factor levels for new data", {
  m <- accessor_fit()
  nd <- accessor_newdata()[1, , drop = FALSE]
  nd$g <- factor("a")
  mf <- stats::model.frame(m, data = nd)

  expect_identical(levels(mf$g), levels(accessor_data()$g))
  expect_identical(
    colnames(stats::model.matrix(m, data = nd)),
    colnames(m@model_spec$X)
  )
})

test_that("model.frame() evaluates a data-dependent term at its fitted basis", {
  set.seed(7)
  n <- 50
  data <- data.frame(x = round(stats::rnorm(n), 3))
  data$y <- round(1 + 0.5 * data$x - 0.4 * data$x^2 + stats::rnorm(n), 3)
  m <- m_estimate(
    y ~ poly(x, 2),
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  nd <- data.frame(y = c(0, 0, 0), x = c(-1, 0, 1))

  # The recorded terms carry the predvars of poly(x, 2), so the basis is the
  # one the fit was made with rather than one refitted to three rows.
  expect_equal(
    unname(drop(stats::model.matrix(m, data = nd) %*% coef(m))),
    unname(predict(m, newdata = nd))
  )
})

# model.matrix() ---------------------------------------------------------------
#
# A design matrix is the fit's, not the session's. Factor coding comes from the
# contrasts the fit recorded, so a fit made under one setting of
# `options(contrasts = )` and read back under another answers with the coding it
# was solved on. This is the rule `predict()` already follows, which rebuilds
# its design with the recorded contrasts.
#
# A design matrix also names no response, so covariate values that carry none
# are enough to build one. That is where it parts company with `model.frame()`,
# which holds every variable the formula names and says so when the response is
# absent.

test_that("model.matrix() codes a factor the way the fit coded it", {
  data <- accessor_data()

  # Fitted under sum contrasts and read back under the default treatment
  # contrasts. The coding belongs to the fit, which recorded it, and not to the
  # option in force when the design is asked for. The restore is immediate
  # because the read has to happen under the default; the exit handler covers
  # the path where the fit fails before reaching it.
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  m <- accessor_fit(data)
  options(old)

  mm <- stats::model.matrix(m)

  expect_identical(attr(mm, "contrasts"), m@model_spec$contrasts)
  expect_identical(colnames(mm), c("(Intercept)", "x", "g1", "g2"))
  expect_identical(mm, m@model_spec$X)
})

test_that("model.matrix() builds a design from covariate values alone", {
  m <- accessor_fit()
  covariates <- accessor_newdata()[, c("x", "g")]

  mm <- stats::model.matrix(m, data = covariates)

  expect_identical(nrow(mm), nrow(covariates))
  expect_identical(colnames(mm), colnames(m@model_spec$X))
  # The design at the caller's covariate values is the one predict() forms, so
  # multiplying it by the coefficients gives the predictions.
  expect_equal(
    unname(drop(mm %*% coef(m))),
    unname(predict(m, newdata = covariates))
  )
})

test_that("model.matrix() codes new data from the recorded contrasts", {
  data <- accessor_data()

  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  m <- accessor_fit(data)
  options(old)

  covariates <- accessor_newdata()[, c("x", "g")]
  mm <- stats::model.matrix(m, data = covariates)

  expect_identical(attr(mm, "contrasts"), m@model_spec$contrasts)
  expect_identical(colnames(mm), colnames(m@model_spec$X))
  expect_equal(
    unname(drop(mm %*% coef(m))),
    unname(predict(m, newdata = covariates))
  )
})

test_that("model.matrix() aborts for a fit built from a function", {
  # A fit with no recorded specification has no design to describe, the same
  # answer the other specification accessors give.
  expect_error(stats::model.matrix(function_fit()), "formula interface")
})

test_that("model.frame() honors subset and na.action beside data", {
  m <- accessor_fit()
  nd <- accessor_newdata()

  expect_identical(nrow(stats::model.frame(m, data = nd, subset = 1:2)), 2L)
  # The expression form is evaluated in the data, as stats::model.frame() does.
  expect_identical(
    nrow(stats::model.frame(m, data = nd, subset = g == "a")),
    1L
  )

  nd$x[2] <- NA
  expect_identical(nrow(stats::model.frame(m, data = nd)), 2L)
  expect_identical(
    nrow(stats::model.frame(m, data = nd, na.action = stats::na.pass)),
    3L
  )
})

test_that("model.frame() aborts on data that carries no response", {
  m <- accessor_fit()
  nd <- accessor_newdata()[, c("x", "g")]

  err <- expect_error(stats::model.frame(m, data = nd), "response")
  expect_match(flatten_message(err), "y", fixed = TRUE)
  expect_match(flatten_message(err), "newdata", fixed = TRUE)
})

test_that("model.frame() aborts on an argument that needs data", {
  m <- accessor_fit()

  expect_error(stats::model.frame(m, subset = 1:2), "subset")
  expect_error(
    stats::model.frame(m, na.action = stats::na.pass),
    "na\\.action"
  )
  expect_error(
    stats::model.frame(m, drop.unused.levels = TRUE),
    "drop\\.unused\\.levels"
  )
  expect_error(
    stats::model.frame(m, xlev = list(g = c("a", "b", "c"))),
    "xlev"
  )
})

test_that("model.frame() aborts on an argument it does not recognize", {
  m <- accessor_fit()

  # A misspelled `data` would otherwise report the fitted frame while the
  # caller believed they had asked for a frame of their own.
  expect_error(
    stats::model.frame(m, newdata = accessor_newdata()),
    "newdata"
  )
})

# formula ----------------------------------------------------------------------

test_that("formula() returns the formula the fit was written with", {
  m <- accessor_fit()
  f <- stats::formula(m)

  expect_s3_class(f, "formula")
  expect_identical(class(f), "formula")
  expect_identical(deparse1(f), "y ~ x + g")
})

test_that("formula() keeps an offset term written into the formula", {
  data <- accessor_data()
  data$exposure <- round(stats::runif(nrow(data), 1, 5), 3)
  m <- m_estimate(
    count ~ x + offset(log(exposure)),
    data = data,
    .ee = ee_regression,
    model = "poisson"
  )
  expect_identical(
    deparse1(stats::formula(m)),
    "count ~ x + offset(log(exposure))"
  )
})

# terms ------------------------------------------------------------------------

test_that("terms() returns the recorded terms", {
  m <- accessor_fit()
  tt <- stats::terms(m)

  expect_s3_class(tt, "terms")
  expect_identical(tt, m@model_spec$terms)
  expect_identical(attr(tt, "term.labels"), c("x", "g"))
  expect_identical(attr(tt, "response"), 1L)
})

test_that("terms() carries the predvars of a data-dependent term", {
  set.seed(7)
  n <- 50
  data <- data.frame(x = round(stats::rnorm(n), 3))
  data$y <- round(1 + 0.5 * data$x - 0.4 * data$x^2 + stats::rnorm(n), 3)
  m <- m_estimate(
    y ~ poly(x, 2),
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  expect_false(is.null(attr(stats::terms(m), "predvars")))
})

# Missing values ---------------------------------------------------------------

test_that("the accessors report the complete rows a fit was solved on", {
  data <- accessor_data()
  data$x[c(5, 11)] <- NA
  m <- m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
  reference <- stats::lm(y ~ x + g, data = data)

  # The fit was solved on the complete rows, so every accessor answers for
  # those alone, with the row names lm() gives them under the default na.omit.
  expect_identical(nobs(m), nrow(data) - 2L)
  expect_length(fitted(m), nobs(m))
  expect_identical(names(fitted(m)), names(stats::fitted(reference)))
  expect_identical(nrow(stats::model.frame(m)), nobs(m))
  expect_identical(
    rownames(stats::model.frame(m)),
    names(stats::fitted(reference))
  )
  expect_equal(residuals(m), stats::residuals(reference), tolerance = 1e-6)
})

# GMM --------------------------------------------------------------------------

test_that("the accessors answer for a GMM formula fit", {
  data <- accessor_data()
  g <- gmm_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )

  expect_length(fitted(g), nobs(g))
  expect_equal(unname(residuals(g)), data$y - unname(fitted(g)))
  expect_identical(nrow(stats::model.frame(g)), nobs(g))
  expect_identical(deparse1(stats::formula(g)), "y ~ x + g")
  expect_identical(attr(stats::terms(g), "term.labels"), c("x", "g"))
})

# Fits the accessors cannot answer for -----------------------------------------

test_that("fitted() and residuals() abort for a fit built from a function", {
  m <- function_fit()

  for (f in list(fitted, residuals)) {
    err <- expect_error(f(m), "formula interface")
    expect_match(flatten_message(err), "regression_predictions()", fixed = TRUE)
  }
  expect_match(
    flatten_message(expect_error(fitted(m))),
    "fitted()",
    fixed = TRUE
  )
  expect_match(
    flatten_message(expect_error(residuals(m))),
    "residuals()",
    fixed = TRUE
  )
})

test_that("model.frame(), formula(), and terms() abort for a function fit", {
  m <- function_fit()

  expect_match(
    flatten_message(expect_error(stats::model.frame(m), "formula interface")),
    "model.frame()",
    fixed = TRUE
  )
  expect_match(
    flatten_message(expect_error(stats::formula(m), "formula interface")),
    "formula()",
    fixed = TRUE
  )
  expect_match(
    flatten_message(expect_error(stats::terms(m), "formula interface")),
    "terms()",
    fixed = TRUE
  )
})

test_that("fitted() and residuals() abort for an equation predict() declines", {
  m <- aft_fit()

  # The abort names the function the user called rather than predict(), which
  # every one of them reaches the supported-equation table through.
  expect_match(
    flatten_message(expect_error(fitted(m), "linear predictor")),
    "fitted()",
    fixed = TRUE
  )
  expect_match(
    flatten_message(expect_error(residuals(m), "linear predictor")),
    "residuals()",
    fixed = TRUE
  )
})

test_that("the specification accessors answer for an equation predict() declines", {
  m <- aft_fit()

  # model.frame(), formula(), and terms() report what the fit was specified as,
  # which every formula fit records whether or not it can be predicted from.
  expect_identical(deparse1(stats::formula(m)), "time ~ x")
  expect_identical(attr(stats::terms(m), "term.labels"), "x")
  expect_identical(names(stats::model.frame(m)), c("time", "x"))
})

test_that("fitted() and residuals() abort before an estimator is estimated", {
  y <- c(1, 2, 3, 4, 5)
  m <- MEstimator(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(0)
  )
  expect_error(fitted(m), "before calling")
  expect_error(residuals(m), "before calling")
})

# df.residual, weights, sigma, logLik, and deviance ----------------------------
#
# Five more base generics a fit is asked in the ordinary course of using a model
# object, and the two answers they can give divide the same way the accessors
# above do.
#
# `df.residual()` is the residual degrees of freedom, the observations less the
# parameters, and both counts are properties of every fit however it was built,
# so it answers for a `stacked_equations` fit as readily as for a formula one.
# `weights()` reports the weights the fit was solved with, which reach a fit
# only through the formula interface, so it is a specification accessor and
# aborts like the rest of them for a fit that records none. An unweighted
# formula fit returns `NULL`, which is what `weights()` answers for an
# unweighted [stats::lm()] fit; the alternative, a vector of ones, would be
# indistinguishable from a fit weighted by ones and would say a fit had weights
# when the caller supplied none.
#
# `sigma()`, `logLik()`, and `deviance()` are the three that have no answer at
# all, and each has to say so rather than report a wrong one. An M-estimator is
# defined by its estimating equations, which need come from no likelihood, so
# there is no log-likelihood to return, no deviance to take from one, and no
# residual standard deviation the object records: a linear equation has one and a
# logistic equation has none, and nothing in `theta` says which was solved. Every
# one of those defaults is worse than an error. `deviance.default()` reads
# `object$deviance`, so on an S7 object the reported failure is about `$` and
# names a property the class does not have; `sigma.default()` reaches the same
# failure through it; and `logLik()` has no default at all, so it reports that no
# method applies to a class the package wrote the rest of its methods for.

accessor_weighted_fit <- function(data = accessor_data()) {
  data$w <- rep(c(0.5, 1.5), length.out = nrow(data))
  m_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
}

test_that("df.residual() is the observations less the parameters", {
  m <- accessor_fit()

  expect_equal(stats::df.residual(m), 56)
  expect_equal(stats::df.residual(m), nobs(m) - m@n_params)
})

test_that("df.residual() answers for a GMMEstimator", {
  m <- gmm_estimate(
    y ~ x + g,
    data = accessor_data(),
    .ee = ee_regression,
    model = "linear"
  )

  expect_equal(stats::df.residual(m), 56)
})

test_that("df.residual() answers for a fit built from a function", {
  m <- function_fit()

  # Both counts are properties of the fit rather than of a specification, so
  # this is not one of the accessors a function fit has no answer for.
  expect_equal(stats::df.residual(m), nobs(m) - m@n_params)
  expect_equal(stats::df.residual(m), 30)
})

test_that("weights() returns the weights a formula fit was solved with", {
  data <- accessor_data()
  m <- accessor_weighted_fit(data)
  used <- rep(c(0.5, 1.5), length.out = nrow(data))

  expect_equal(unname(stats::weights(m)), used)
  expect_length(stats::weights(m), nobs(m))
})

test_that("weights() reports time-varying weights as the matrix they were", {
  # `ee_plogit()` takes an n-by-K matrix of weights, one column per time
  # interval, and an observation weighted differently in each interval has no
  # one weight a vector could hold. The matrix is reported as it was supplied,
  # which the documented return says and this pins: a reshape into a vector
  # would lose which interval each weight belonged to.
  #
  # Four unit-time intervals with events in three of them, and censoring in the
  # last, so the time parameters are all finite and the fit solves from the
  # start below without a solver warning.
  data <- data.frame(
    x = rep(c(0, 1), length.out = 120),
    time = rep(c(1, 2, 3, 4, 2, 4), length.out = 120),
    event = rep(c(1, 1, 0, 1, 1, 0), length.out = 120)
  )
  k <- length(unique(data$time[data$event == 1]))
  w <- matrix(
    seq_len(nrow(data) * k) / (nrow(data) * k),
    nrow = nrow(data),
    ncol = k
  )

  # No intercept: the time parameters carry one per interval, so an intercept
  # column beside them has nothing left to identify it.
  m <- m_estimate(
    time ~ x - 1,
    data = data,
    .ee = ee_plogit,
    event = data$event,
    weights = w,
    init = c(0, rep(-1, k))
  )

  reported <- stats::weights(m)
  expect_true(is.matrix(reported))
  expect_identical(dim(reported), c(nobs(m), k))
  expect_equal(unname(reported), unname(w))
})

test_that("weights() returns NULL for an unweighted formula fit", {
  m <- accessor_fit()

  # The precedent is weights() on an unweighted lm() fit, which is NULL rather
  # than a vector of ones.
  expect_null(stats::weights(m))
})

test_that("weights() aborts for a fit built from a function", {
  m <- function_fit()

  expect_match(
    flatten_message(expect_error(stats::weights(m), "formula interface")),
    "weights()",
    fixed = TRUE
  )
})

test_that("sigma() reports that no residual standard deviation is defined", {
  for (m in list(accessor_fit(), function_fit())) {
    err <- expect_error(stats::sigma(m), "standard deviation")
    expect_match(flatten_message(err), "sigma()", fixed = TRUE)
  }
})

test_that("logLik() reports that an M-estimator defines no likelihood", {
  for (m in list(accessor_fit(), function_fit())) {
    err <- expect_error(stats::logLik(m), "likelihood")
    expect_match(flatten_message(err), "logLik()", fixed = TRUE)
  }
})

test_that("deviance() reports that an M-estimator defines no likelihood", {
  # A deviance is a likelihood quantity, so it has no more of an answer here
  # than the log-likelihood it would be taken from, and it needs a method of its
  # own: without one the reported failure is `sigma.default()`'s, about `$` on a
  # property the class does not have.
  for (m in list(accessor_fit(), function_fit())) {
    err <- expect_error(stats::deviance(m), "likelihood")
    expect_match(flatten_message(err), "deviance()", fixed = TRUE)
  }
})

test_that("nobs() and the specification accessors are unchanged", {
  m <- accessor_fit()

  expect_equal(nobs(m), 60L)
  expect_length(fitted(m), 60L)
  expect_length(residuals(m), 60L)
  expect_identical(names(stats::model.frame(m)), c("y", "x", "g"))
  expect_identical(deparse1(stats::formula(m)), "y ~ x + g")
  expect_identical(attr(stats::terms(m), "term.labels"), c("x", "g"))
})
