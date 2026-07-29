# ---- Helpers ----

# One data frame drives most of the tests in this file. It carries a continuous
# covariate, a three-level factor so that contrasts and levels are exercised, a
# positive exposure to offset by, and a continuous, a binary, a count, and a
# positive response.
predict_data <- function() {
  set.seed(20240719)
  n <- 80
  d <- data.frame(
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b", "c"), length.out = n)),
    exposure = round(stats::runif(n, 1, 5), 3)
  )
  d$y <- round(1 + 0.7 * d$x + as.numeric(d$g) + stats::rnorm(n), 3)
  d$hit <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x))
  d$count <- stats::rpois(n, exp(0.4 + 0.2 * d$x))
  d$wear <- round(
    stats::rgamma(n, shape = 4, rate = 4 / exp(0.5 + 0.3 * d$x)),
    4
  )
  d
}

predict_fit <- function(data = predict_data()) {
  m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
}

flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

# predict() on the link scale --------------------------------------------------

test_that("predict() returns the design matrix times the coefficients", {
  m <- predict_fit()
  expect_identical(
    unname(predict(m)),
    unname(drop(m@model_spec$X %*% coef(m)))
  )
})

test_that("predict() with no newdata predicts the rows the fit was made on", {
  data <- predict_data()
  m <- predict_fit(data)
  p <- predict(m)

  expect_type(p, "double")
  expect_null(dim(p))
  expect_length(p, nrow(data))
  expect_identical(names(p), rownames(m@model_spec$X))
  # Handing back the data the fit was made on rebuilds the same design, so it
  # answers with the same predictions.
  expect_identical(p, predict(m, newdata = data))
})

test_that("predict() rebuilds the design for newdata the fit never saw", {
  data <- predict_data()
  m <- predict_fit(data)
  newdata <- data.frame(
    x = c(-1, 0, 1),
    g = factor("b", levels = levels(data$g))
  )

  expect_identical(
    unname(predict(m, newdata = newdata)),
    unname(drop(cbind(1, c(-1, 0, 1), 1, 0) %*% coef(m)))
  )
})

test_that("predict() rebuilds a data-dependent term through the recorded terms", {
  set.seed(88)
  n <- 60
  data <- data.frame(x = round(stats::rnorm(n), 3))
  data$y <- round(1 + 0.5 * data$x - 0.4 * data$x^2 + stats::rnorm(n), 3)
  m <- m_estimate(
    y ~ poly(x, 2),
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  newdata <- data[1:5, , drop = FALSE]

  # poly() carries its fitted coefficients in the terms' predvars attribute, so
  # predicting on a subset reproduces the corresponding fitted values. Rebuilding
  # the basis from the formula alone would recompute it on the five rows handed
  # over and disagree.
  expect_equal(
    unname(predict(m, newdata = newdata)),
    unname(predict(m)[1:5])
  )
  naive <- stats::model.matrix(y ~ poly(x, 2), data = newdata)
  expect_false(isTRUE(all.equal(
    unname(drop(naive %*% coef(m))),
    unname(predict(m)[1:5])
  )))
})

# predict() on the response scale ----------------------------------------------

test_that("predict() on the response scale of a logistic fit is a probability", {
  m <- m_estimate(
    hit ~ x,
    data = predict_data(),
    .ee = ee_regression,
    model = "logistic"
  )
  eta <- predict(m)

  expect_equal(predict(m, type = "response"), stats::plogis(eta))
  expect_true(all(predict(m, type = "response") > 0))
  expect_true(all(predict(m, type = "response") < 1))
})

test_that("predict() reads the response scale of a glm fit from its link", {
  m <- m_estimate(
    count ~ x,
    data = predict_data(),
    .ee = ee_glm,
    distribution = "poisson",
    link = "log"
  )
  expect_equal(predict(m, type = "response"), exp(predict(m)))
})

test_that("predict() scales the response variance by the link derivative", {
  m <- m_estimate(
    hit ~ x,
    data = predict_data(),
    .ee = ee_regression,
    model = "logistic"
  )
  eta <- predict(m)
  se_link <- predict(m, se.fit = TRUE)$se.fit
  mu <- stats::plogis(eta)

  # The delta-method standard error on the response scale is |dmu/deta| times
  # the standard error of the linear predictor.
  expect_equal(
    predict(m, type = "response", se.fit = TRUE)$se.fit,
    mu * (1 - mu) * se_link
  )
})

test_that("inverse_link() returns the derivative of its own inverse link", {
  # predict() forms the response-scale variance as dmu^2 times the variance of
  # the linear predictor, which is only correct if dmu is d(mu)/d(eta). The
  # comparison is on magnitudes because the "loglog" branch returns the negative
  # of the derivative, matching the Python source this package was ported from.
  # Squaring makes predict() insensitive to that, and nothing here relies on the
  # sign. Zero is left out of the grid because the "inverse" link is undefined
  # there.
  eta <- c(-1.7, -0.4, 0.6, 1.3)
  step <- 1e-6
  links <- c(
    "identity",
    "log",
    "logit",
    "probit",
    "cauchit",
    "loglog",
    "cloglog",
    "inverse",
    "sqrt"
  )
  for (link in links) {
    numeric_deriv <- (inverse_link(eta + step, link)$mu -
      inverse_link(eta - step, link)$mu) /
      (2 * step)
    expect_equal(
      abs(inverse_link(eta, link)$dmu),
      abs(numeric_deriv),
      tolerance = 1e-6,
      info = link
    )
  }
})

test_that("the three regression models are the three links predict() uses", {
  # predict() maps a `model` argument onto a link and evaluates every scale
  # through inverse_link(), so the two must agree wherever both are defined.
  eta <- c(-2, -0.5, 0, 0.5, 2)
  expect_identical(
    model_transform(eta, "linear"),
    inverse_link(eta, "identity")$mu
  )
  expect_identical(
    model_transform(eta, "logistic"),
    inverse_link(eta, "logit")$mu
  )
  expect_identical(model_transform(eta, "poisson"), inverse_link(eta, "log")$mu)
})

# Standard errors and intervals ------------------------------------------------

test_that("predict(se.fit = TRUE) returns the fit beside its standard error", {
  data <- predict_data()
  m <- predict_fit(data)
  p <- predict(m, se.fit = TRUE)

  expect_type(p, "list")
  expect_named(p, c("fit", "se.fit"))
  expect_identical(p$fit, predict(m))
  expect_length(p$se.fit, nrow(data))
  expect_true(all(p$se.fit > 0))
})

test_that("predict(interval = 'confidence') returns a fit, lwr, upr matrix", {
  data <- predict_data()
  m <- predict_fit(data)
  p <- predict(m, interval = "confidence")

  expect_true(is.matrix(p))
  expect_identical(dim(p), c(nrow(data), 3L))
  expect_identical(colnames(p), c("fit", "lwr", "upr"))
  expect_identical(rownames(p), rownames(m@model_spec$X))
  expect_identical(p[, "fit"], predict(m))
  expect_true(all(p[, "lwr"] < p[, "upr"]))
})

test_that("predict() with both se.fit and interval returns the matrix in fit", {
  m <- predict_fit()
  p <- predict(m, se.fit = TRUE, interval = "confidence")

  expect_named(p, c("fit", "se.fit"))
  expect_identical(p$fit, predict(m, interval = "confidence"))
  expect_identical(p$se.fit, predict(m, se.fit = TRUE)$se.fit)
})

test_that("predict() widens its interval as the level rises", {
  m <- predict_fit()
  narrow <- predict(m, interval = "confidence", level = 0.8)
  wide <- predict(m, interval = "confidence", level = 0.99)

  expect_true(all(wide[, "lwr"] < narrow[, "lwr"]))
  expect_true(all(wide[, "upr"] > narrow[, "upr"]))
})

test_that("predict() on the response scale is Wald on the response scale", {
  m <- m_estimate(
    hit ~ x,
    data = predict_data(),
    .ee = ee_regression,
    model = "logistic"
  )
  p <- predict(m, type = "response", interval = "confidence", se.fit = TRUE)
  cv <- stats::qnorm(0.975)

  # The bounds are the response-scale fit plus and minus the response-scale
  # standard error, so they are symmetric about the fitted probability.
  expect_identical(
    unname(p$fit[, "lwr"]),
    unname(p$fit[, "fit"] - cv * p$se.fit)
  )
  expect_identical(
    unname(p$fit[, "upr"]),
    unname(p$fit[, "fit"] + cv * p$se.fit)
  )

  # That is a different interval from the link-scale one put through the
  # inverse link. The two are told apart here by a bound a transformed interval
  # could not produce: this one leaves the unit interval at both ends, while
  # the inverse logit of anything is a probability.
  link <- predict(m, interval = "confidence")
  expect_false(isTRUE(all.equal(
    unname(p$fit[, "lwr"]),
    unname(stats::plogis(link[, "lwr"]))
  )))
  expect_true(any(p$fit[, "lwr"] < 0))
  expect_true(any(p$fit[, "upr"] > 1))
})

test_that("predict() takes a t critical value under a finite correction", {
  data <- predict_data()
  m <- m_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_regression,
    model = "linear",
    finite_correction = "HC1"
  )
  p <- predict(m, interval = "confidence", se.fit = TRUE)
  cv <- stats::qt(0.975, df = m@n_obs - m@n_params)

  expect_identical(
    unname(p$fit[, "lwr"]),
    unname(p$fit[, "fit"] - cv * p$se.fit)
  )
  expect_identical(
    unname(p$fit[, "upr"]),
    unname(p$fit[, "fit"] + cv * p$se.fit)
  )

  # It is the critical value confint() uses, read back off the width of its
  # bounds in standard errors.
  ci <- confint(m)
  expect_equal(
    unname((ci[, "upper"] - coef(m)) / sqrt(diag(vcov(m)))),
    rep(cv, m@n_params)
  )

  # Without a correction the critical value is the normal quantile, so the
  # correction is doing the work rather than the two quantiles agreeing.
  expect_gt(cv, stats::qnorm(0.975))
  plain <- predict(predict_fit(data), interval = "confidence", se.fit = TRUE)
  expect_identical(
    unname(plain$fit[, "upr"]),
    unname(plain$fit[, "fit"] + stats::qnorm(0.975) * plain$se.fit)
  )
})

# Agreement with regression_predictions() --------------------------------------

test_that("predict() matches regression_predictions() exactly", {
  data <- predict_data()
  m <- predict_fit(data)
  reference <- regression_predictions(
    m@model_spec$X,
    theta = coef(m),
    covariance = vcov(m)
  )
  p <- predict(m, se.fit = TRUE, interval = "confidence")

  expect_identical(unname(p$fit[, "fit"]), reference$predicted)
  expect_identical(unname(p$fit[, "lwr"]), reference$lower)
  expect_identical(unname(p$fit[, "upr"]), reference$upper)
  expect_identical(unname(p$se.fit), sqrt(reference$variance))
})

test_that("predict() on newdata matches regression_predictions() exactly", {
  data <- predict_data()
  m <- predict_fit(data)
  newdata <- data[c(3, 17, 42), , drop = FALSE]
  X_new <- stats::model.matrix(
    stats::delete.response(m@model_spec$terms),
    stats::model.frame(
      stats::delete.response(m@model_spec$terms),
      newdata,
      xlev = m@model_spec$xlevels
    ),
    contrasts.arg = m@model_spec$contrasts
  )
  reference <- regression_predictions(
    X_new,
    theta = coef(m),
    covariance = vcov(m)
  )
  p <- predict(m, newdata = newdata, se.fit = TRUE, interval = "confidence")

  expect_identical(unname(p$fit[, "fit"]), reference$predicted)
  expect_identical(unname(p$fit[, "lwr"]), reference$lower)
  expect_identical(unname(p$fit[, "upr"]), reference$upper)
  expect_identical(unname(p$se.fit), sqrt(reference$variance))
})

test_that("predict() matches regression_predictions() under an offset", {
  data <- predict_data()
  m <- m_estimate(
    count ~ x + offset(log(exposure)),
    data = data,
    .ee = ee_regression,
    model = "poisson"
  )
  reference <- regression_predictions(
    m@model_spec$X,
    theta = coef(m),
    covariance = vcov(m),
    offset = log(data$exposure)
  )
  p <- predict(m, se.fit = TRUE, interval = "confidence")

  expect_identical(unname(p$fit[, "fit"]), reference$predicted)
  expect_identical(unname(p$fit[, "lwr"]), reference$lower)
  expect_identical(unname(p$fit[, "upr"]), reference$upper)
  expect_identical(unname(p$se.fit), sqrt(reference$variance))
})

# The supported equations ------------------------------------------------------

test_that("predict_link_name() names the link of every equation it supports", {
  expect_identical(
    predict_link_name(ee_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_regression, list(model = "logistic")),
    "logit"
  )
  expect_identical(
    predict_link_name(ee_regression, list(model = "poisson")),
    "log"
  )
  expect_identical(
    predict_link_name(ee_ridge_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_bridge_regression, list(model = "logistic")),
    "logit"
  )
  expect_identical(
    predict_link_name(ee_lasso_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_dlasso_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_elasticnet_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_robust_regression, list(model = "linear")),
    "identity"
  )
  expect_identical(
    predict_link_name(ee_glm, list(distribution = "gamma", link = "log")),
    "log"
  )
  expect_identical(
    predict_link_name(ee_glm, list(distribution = "binomial", link = "probit")),
    "probit"
  )
  expect_identical(predict_link_name(ee_beta_regression, list()), "logit")
})

test_that("predict_link_name() declines every equation it does not support", {
  expect_null(predict_link_name(ee_gformula, list()))
  expect_null(predict_link_name(ee_additive_regression, list(model = "linear")))
  expect_null(predict_link_name(ee_tobit, list()))
  expect_null(predict_link_name(ee_mlogit, list()))
  expect_null(predict_link_name(ee_aft, list(distribution = "weibull")))
  expect_null(predict_link_name(ee_plogit, list()))
  expect_null(predict_link_name(function(theta, X, y) NULL, list()))
  # No fit reaches this, since model_transform() rejects an unsupported model
  # while the fit is being solved, but the table answers for it all the same.
  expect_null(predict_link_name(ee_regression, list(model = "probit")))
})

test_that("predict() predicts from a penalized regression fit", {
  data <- predict_data()
  m <- m_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_ridge_regression,
    model = "linear",
    penalty = c(0, 5, 5, 5)
  )
  expect_identical(
    unname(predict(m)),
    unname(drop(m@model_spec$X %*% coef(m)))
  )
})

test_that("predict() predicts from a robust regression fit", {
  start <- coef(stats::lm(weight ~ height, data = robust_regress))
  m <- m_estimate(
    weight ~ height,
    data = robust_regress,
    .ee = ee_robust_regression,
    model = "linear",
    k = 1.345,
    loss = "huber",
    init = start
  )
  expect_identical(
    unname(predict(m)),
    unname(drop(m@model_spec$X %*% coef(m)))
  )
})

test_that("predict() predicts from a GMM fit made with a formula", {
  data <- predict_data()
  m <- gmm_estimate(
    y ~ x + g,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
  expect_identical(
    unname(predict(m)),
    unname(drop(m@model_spec$X %*% coef(m)))
  )
})

# The coefficient block --------------------------------------------------------

test_that("predict() on a gamma glm fit uses only the coefficient block", {
  data <- predict_data()
  m <- m_estimate(
    wear ~ x,
    data = data,
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )
  # theta is c(beta, log_shape), so slicing to the parameter count would put the
  # log shape into the linear predictor.
  expect_identical(names(coef(m)), c("(Intercept)", "x", "log_shape"))
  beta <- coef(m)[1:2]
  X <- m@model_spec$X

  expect_identical(unname(predict(m)), unname(drop(X %*% beta)))
  expect_identical(
    unname(predict(m, se.fit = TRUE)$se.fit),
    unname(sqrt(rowSums((X %*% vcov(m)[1:2, 1:2]) * X)))
  )
})

test_that("predict() on a beta regression fit uses only the mean model", {
  skip_if_not_installed("nleqslv")
  set.seed(517)
  n <- 120
  data <- data.frame(x = round(stats::rnorm(n), 3))
  mu <- stats::plogis(0.5 + 0.8 * data$x)
  data$prop <- round(
    stats::rbeta(n, shape1 = mu * 12, shape2 = (1 - mu) * 12),
    5
  )
  # The default solver diverges on this equation, as ee_beta_regression()'s own
  # example notes.
  m <- m_estimate(
    prop ~ x,
    data = data,
    .ee = ee_beta_regression,
    init = c(0, 0, log(10)),
    solver = "nleqslv"
  )

  # theta is c(beta, log_phi), one parameter more than the design has columns,
  # so slicing to the parameter count would put the log precision into the
  # linear predictor.
  expect_identical(names(coef(m)), c("(Intercept)", "x", "log_phi"))
  beta <- coef(m)[1:2]
  X <- m@model_spec$X

  expect_identical(unname(predict(m)), unname(drop(X %*% beta)))
  expect_identical(
    unname(predict(m, se.fit = TRUE)$se.fit),
    unname(sqrt(rowSums((X %*% vcov(m)[1:2, 1:2]) * X)))
  )
  # The mean model is logit by construction, so the response scale is the
  # fitted proportion.
  expect_identical(predict(m, type = "response"), stats::plogis(predict(m)))
})

test_that("predict() rejects a parameter count the design cannot account for", {
  m <- predict_fit()
  spec <- m@model_spec

  # The design coefficient count and the parameter count are recorded
  # separately, so a fit whose parameters the design cannot account for is
  # caught rather than sliced.
  expect_error(
    predict_coef_index(spec, m@n_params + 2L),
    "does not account for"
  )
  expect_identical(predict_coef_index(spec, m@n_params), seq_len(spec$n_coef))
})

# Offsets ----------------------------------------------------------------------

test_that("predict() applies the offset a formula term supplies for newdata", {
  data <- predict_data()
  m <- m_estimate(
    count ~ x + offset(log(exposure)),
    data = data,
    .ee = ee_regression,
    model = "poisson"
  )
  newdata <- data.frame(x = c(-0.5, 0.5), exposure = c(2, 4))

  expect_identical(
    unname(predict(m, newdata = newdata)),
    unname(drop(cbind(1, c(-0.5, 0.5)) %*% coef(m)) + log(c(2, 4)))
  )
  # A different exposure moves the prediction by the difference in log exposure.
  shifted <- data.frame(x = c(-0.5, 0.5), exposure = c(4, 8))
  expect_equal(
    predict(m, newdata = shifted) - predict(m, newdata = newdata),
    stats::setNames(rep(log(2), 2), c("1", "2"))
  )
})

test_that("predict() aborts when newdata cannot supply the formula offset", {
  data <- predict_data()
  m <- m_estimate(
    count ~ x + offset(log(exposure)),
    data = data,
    .ee = ee_regression,
    model = "poisson"
  )
  expect_error(
    predict(m, newdata = data.frame(x = c(-0.5, 0.5))),
    "exposure"
  )
})

test_that("predict() aborts on newdata when the offset came through the dots", {
  data <- predict_data()
  m <- m_estimate(
    count ~ x,
    data = data,
    .ee = ee_regression,
    model = "poisson",
    offset = log(exposure)
  )
  # The fitted rows keep the offset the fit was given.
  expect_identical(
    unname(predict(m)),
    unname(drop(m@model_spec$X %*% coef(m)) + log(data$exposure))
  )

  err <- expect_error(
    predict(m, newdata = data.frame(x = c(-0.5, 0.5))),
    "offset"
  )
  flat <- flatten_message(err)
  expect_match(flat, "one value per fitted observation", fixed = TRUE)
  expect_match(flat, "offset()", fixed = TRUE)
})

# newdata that cannot rebuild the design ---------------------------------------

test_that("predict() aborts for newdata carrying an unseen factor level", {
  data <- predict_data()
  m <- predict_fit(data)
  newdata <- data.frame(
    x = c(0, 1),
    g = factor(c("d", "d"), levels = c(levels(data$g), "d"))
  )
  expect_error(predict(m, newdata = newdata), "new level")
})

test_that("predict() aborts for newdata missing a predictor", {
  data <- predict_data()
  m <- predict_fit(data)
  expect_error(
    predict(m, newdata = data.frame(x = c(0, 1))),
    "'g' not found"
  )
})

# Missing values ---------------------------------------------------------------

test_that("predict() answers NA for a newdata row with a missing value", {
  data <- predict_data()
  m <- predict_fit(data)
  newdata <- data.frame(
    x = c(-1, NA, 1),
    g = factor(c("a", "b", "c"), levels = levels(data$g))
  )
  p <- predict(m, newdata = newdata, se.fit = TRUE, interval = "confidence")

  # The row passes through rather than being dropped, so the predictions line
  # up with the rows of newdata and the complete rows answer as they would on
  # their own.
  expect_identical(nrow(p$fit), 3L)
  expect_identical(unname(is.na(p$fit[, "fit"])), c(FALSE, TRUE, FALSE))
  expect_identical(unname(is.na(p$fit[, "lwr"])), c(FALSE, TRUE, FALSE))
  expect_identical(unname(is.na(p$se.fit)), c(FALSE, TRUE, FALSE))
  expect_identical(
    unname(p$fit[c(1, 3), "fit"]),
    unname(predict(m, newdata = newdata[c(1, 3), , drop = FALSE]))
  )
})

test_that("predict() on a fit made with missing data answers its fitted rows", {
  data <- predict_data()
  data$x[c(5, 11)] <- NA
  m <- m_estimate(y ~ x + g, data = data, .ee = ee_regression, model = "linear")
  reference <- stats::lm(y ~ x + g, data = data)

  # The fit was solved on the complete rows, so predicting the rows it was made
  # on answers for those alone, with the row names predict.lm() gives them
  # under the default na.omit.
  expect_length(predict(m), nrow(data) - 2L)
  expect_identical(names(predict(m)), names(stats::predict(reference)))
  expect_equal(unname(predict(m)), unname(stats::predict(reference)))
})

# Fits predict() cannot make predictions from -----------------------------------

test_that("predict() aborts for a fit built from a stacked_equations function", {
  y <- mtcars$mpg
  X <- cbind(1, mtcars$wt)
  m <- m_estimate(
    stacked_equations = function(theta) {
      ee_regression(theta, X = X, y = y, model = "linear")
    },
    init = c(0, 0)
  )

  err <- expect_error(predict(m), "formula interface")
  flat <- flatten_message(err)
  expect_match(flat, "regression_predictions()", fixed = TRUE)
})

test_that("predict() aborts for a g-formula fit and names the reason", {
  data <- predict_data()
  data$a <- stats::rbinom(nrow(data), 1, 0.5)
  data$out <- round(2 + 1.5 * data$a + data$x + stats::rnorm(nrow(data)), 3)
  X1 <- cbind(1, 1, data$x)
  m <- m_estimate(
    out ~ a + x,
    data = data,
    .ee = ee_gformula,
    X1 = X1,
    init = rep(0, 4)
  )

  err <- expect_error(predict(m), "ee_gformula")
  flat <- flatten_message(err)
  expect_match(flat, "causal mean", fixed = TRUE)
  expect_match(flat, "regression_predictions()", fixed = TRUE)
})

test_that("predict() aborts for an additive fit and names the reason", {
  data <- predict_data()
  specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))
  n_expanded <- ncol(additive_design_matrix(cbind(1, data$x), specs))
  m <- m_estimate(
    y ~ x,
    data = data,
    .ee = ee_additive_regression,
    specifications = specs,
    model = "linear",
    init = rep(0, n_expanded)
  )

  err <- expect_error(predict(m), "ee_additive_regression")
  flat <- flatten_message(err)
  expect_match(flat, "additive_design_matrix()", fixed = TRUE)
  expect_match(flat, "regression_predictions()", fixed = TRUE)
})

test_that("predict() aborts for a formula fit of an equation it does not know", {
  data <- predict_data()
  linear <- function(theta, X, y) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m <- m_estimate(y ~ x, data = data, .ee = linear)

  err <- expect_error(predict(m), "estimating equation")
  flat <- flatten_message(err)
  expect_match(flat, "ee_regression()", fixed = TRUE)
  expect_match(flat, "regression_predictions()", fixed = TRUE)
})

test_that("predict() aborts for an aft fit, whose response is not a mean", {
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
  # It declines the linear predictor while still offering the `times` surface,
  # so the message says which of the two is unavailable and why.
  err <- expect_error(predict(m), "linear predictor")
  expect_match(flatten_message(err), "log-time")
  expect_match(flatten_message(err), "times", fixed = TRUE)
})

test_that("predict() aborts before an estimator has been estimated", {
  y <- c(1, 2, 3, 4, 5)
  m <- MEstimator(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(0)
  )
  expect_error(predict(m), "before calling")
})

# Argument handling ------------------------------------------------------------

test_that("predict() rejects an unrecognized type or interval", {
  m <- predict_fit()
  expect_error(predict(m, type = "terms"), "arg")
  expect_error(predict(m, interval = "prediction"), "arg")
})

test_that("predict() rejects a level outside the unit interval", {
  m <- predict_fit()
  expect_error(predict(m, interval = "confidence", level = 0), "level")
  expect_error(predict(m, interval = "confidence", level = 1), "level")
  expect_error(
    predict(m, interval = "confidence", level = c(0.9, 0.95)),
    "level"
  )
})

test_that("predict() rejects an argument it does not take", {
  m <- predict_fit()
  err <- expect_error(
    predict(m, scale = 2),
    class = "rlib_error_dots_nonempty"
  )
  expect_match(flatten_message(err), "scale = 2", fixed = TRUE)
})

test_that("predict() attributes an unrecognized argument to the predict call", {
  m <- predict_fit()

  # The call the error names is the one the caller wrote, rather than this
  # method's own frame or the function the call sits inside.
  err <- expect_error(
    predict(m, scale = 2),
    class = "rlib_error_dots_nonempty"
  )
  expect_identical(conditionCall(err), quote(predict(m, scale = 2)))

  caller <- function(fit) predict(fit, scale = 2)
  from_caller <- expect_error(caller(m), class = "rlib_error_dots_nonempty")
  expect_identical(conditionCall(from_caller), quote(predict(fit, scale = 2)))
})
