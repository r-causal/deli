# ---- Helpers ----

# One data frame drives the AFT tests. The errors are logistic, which is the
# log-logistic AFT's own error distribution, and every one of the four
# distributions converges on it under the default solver.
aft_data <- function() {
  set.seed(414)
  n <- 200
  d <- data.frame(
    x = round(stats::rnorm(n), 3),
    g = factor(rep(c("a", "b"), length.out = n))
  )
  eps <- stats::rlogis(n)
  event_time <- exp(1 + 0.5 * d$x + 0.3 * (d$g == "b") + 0.7 * eps)
  censor_time <- stats::rexp(n, rate = 0.03)
  d$time <- round(pmin(event_time, censor_time), 4)
  d$status <- as.numeric(event_time <= censor_time)
  d
}

aft_distributions <- function() {
  c("exponential", "weibull", "log-normal", "log-logistic")
}

aft_measures <- function() {
  c("survival", "risk", "cumulative_hazard", "hazard", "density")
}

aft_fit <- function(distribution = "weibull", data = aft_data()) {
  n_beta <- 3L
  init <- c(mean(log(data$time)), rep(0, n_beta - 1L))
  if (distribution != "exponential") {
    init <- c(init, 0)
  }
  m_estimate(
    time ~ x + g,
    data = data,
    .ee = ee_aft,
    distribution = distribution,
    event = status,
    init = init
  )
}

plogit_data <- function() {
  d <- collett_bladder
  d$novel <- d$treat - 1
  d
}

plogit_fit <- function(data = plogit_data()) {
  k <- length(unique(data$time[data$delta == 1]))
  m_estimate(
    time ~ novel + init + size - 1,
    data = data,
    .ee = ee_plogit,
    event = delta,
    init = c(rep(0, 3), -4, rep(0, k - 1))
  )
}

regression_fit <- function() {
  m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
}

# The predictions `aft_predictions_individual()` makes for the same design,
# flattened one observation at a time so they line up with the rows of the frame
# `predict()` returns.
#
# The two agree to about twenty-five units in the last place rather than to the
# bit, and the difference is the linear predictor. `aft_predictions_individual()`
# forms it with a matrix product, and `predict()` accumulates it one design
# column at a time so that the same transform can be differentiated exactly. The
# worst relative difference measured over all four distributions and all five
# measures is 5.3e-15, so these comparisons are made at 1e-13, which is still two
# orders of magnitude tighter than the Python parity tests elsewhere in the
# suite. `aft_predictions_function()` accumulates the linear predictor the same
# way `predict()` does and is matched more tightly below.
aft_value_tolerance <- 1e-13

aft_helper_long <- function(m, times, measure, X = m@model_spec$X) {
  wide <- aft_predictions_individual(
    X = X,
    times = times,
    theta = m@theta,
    distribution = m@model_spec$ee_spec_args$distribution,
    measure = measure
  )
  as.numeric(t(as.matrix(wide)))
}

flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

# predict() point values against the helpers -----------------------------------

test_that("predict() matches aft_predictions_individual() per distribution", {
  data <- aft_data()
  times <- c(2, 8, 25)

  for (distribution in aft_distributions()) {
    m <- aft_fit(distribution, data)
    p <- predict(m, times = times, measure = "survival")
    expect_equal(
      p$fit,
      aft_helper_long(m, times, "survival"),
      tolerance = aft_value_tolerance,
      label = distribution
    )
  }
})

test_that("predict() matches aft_predictions_individual() per measure", {
  data <- aft_data()
  times <- c(2, 8, 25)
  m <- aft_fit("weibull", data)

  for (measure in aft_measures()) {
    p <- predict(m, times = times, measure = measure)
    expect_equal(
      p$fit,
      aft_helper_long(m, times, measure),
      tolerance = aft_value_tolerance,
      label = measure
    )
  }
})

test_that("predict() defaults to the survival measure", {
  m <- aft_fit()
  expect_identical(
    predict(m, times = c(3, 9)),
    predict(m, times = c(3, 9), measure = "survival")
  )
})

test_that("predict() matches plogit_predict() per measure", {
  data <- plogit_data()
  m <- plogit_fit(data)
  times <- c(12, 24, 59)
  W <- m@model_spec$X

  for (measure in aft_measures()) {
    p <- predict(m, times = times, measure = measure)
    expected <- plogit_predict(
      m@theta,
      time = data$time,
      event = data$delta,
      X = W,
      times_to_predict = times,
      measure = measure
    )
    # plogit_predict() is times-by-individuals, so its column-major flattening
    # is already one observation at a time.
    expect_equal(
      p$fit,
      as.numeric(expected),
      tolerance = 1e-14,
      label = measure
    )
  }
})

# The shape of the survival surface --------------------------------------------

test_that("predict() returns one row per observation and time", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  times <- c(2, 8, 25)
  p <- predict(m, times = times)

  expect_s3_class(p, "data.frame")
  expect_identical(names(p), c(".row", "time", "fit"))
  expect_identical(nrow(p), nrow(data) * length(times))
  # Observation-major: every time for the first row, then every time for the
  # second, so a one-row design comes back as one row per time.
  expect_identical(p$.row, rep(rownames(m@model_spec$X), each = length(times)))
  expect_identical(p$time, rep(times, times = nrow(data)))
})

test_that("predict() predicts newdata the fit never saw", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(x = c(-1, 0.5), g = factor("b", levels = levels(data$g)))
  times <- c(4, 16)

  p <- predict(m, newdata = at, times = times)
  expect_identical(nrow(p), 4L)
  expect_equal(
    p$fit,
    aft_helper_long(m, times, "survival", X = cbind(1, c(-1, 0.5), 1)),
    tolerance = aft_value_tolerance
  )
})

test_that("predict() reports NA for the newdata rows that carry one", {
  # `predict_design()` builds the design with `na.action = na.pass` so the
  # predictions line up with the rows of `newdata` rather than with whichever
  # of them survived. A missing covariate therefore travels through the
  # transform and the Jacobian into the row it belongs to, and only into that
  # row: an interval is reported for the complete patterns either side of it.
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(
    x = c(-1, NA, 0.5),
    g = factor("b", levels = levels(data$g))
  )
  times <- c(4, 16)
  p <- predict(
    m,
    newdata = at,
    times = times,
    se.fit = TRUE,
    interval = "confidence"
  )

  missing_rows <- p$fit$.row == "2"
  expect_identical(sum(missing_rows), length(times))
  expect_true(all(is.na(p$fit[missing_rows, c("fit", "lwr", "upr")])))
  expect_true(all(is.na(p$se.fit[missing_rows])))
  expect_false(anyNA(p$fit[!missing_rows, c("fit", "lwr", "upr")]))
  expect_false(anyNA(p$se.fit[!missing_rows]))

  # The complete rows are the predictions they would have got on their own, so
  # the missing row costs its neighbours nothing.
  complete <- predict(
    m,
    newdata = at[c(1, 3), ],
    times = times,
    se.fit = TRUE,
    interval = "confidence"
  )
  expect_equal(
    unname(as.matrix(p$fit[!missing_rows, c("fit", "lwr", "upr")])),
    unname(as.matrix(complete$fit[, c("fit", "lwr", "upr")])),
    tolerance = 1e-14
  )
  expect_equal(
    unname(p$se.fit[!missing_rows]),
    unname(complete$se.fit),
    tolerance = 1e-14
  )
})

test_that("predict() carries an offset written into the formula", {
  data <- aft_data()
  data$shift <- round(seq(-0.2, 0.2, length.out = nrow(data)), 4)
  m <- m_estimate(
    time ~ x + offset(shift),
    data = data,
    .ee = ee_aft,
    distribution = "weibull",
    event = status,
    init = c(mean(log(data$time)), 0, 0)
  )
  p <- predict(m, times = 10)

  # The offset shifts the linear predictor, so the same design without it
  # predicts something else.
  plain <- aft_predictions_individual(
    X = m@model_spec$X,
    times = 10,
    theta = m@theta,
    distribution = "weibull",
    measure = "survival"
  )
  expect_false(isTRUE(all.equal(p$fit, as.numeric(plain[[1]]))))

  shifted <- aft_measure_at_time(
    10,
    as.numeric(m@model_spec$X %*% m@theta[1:2]) + data$shift,
    exp(-m@theta[3]),
    "weibull",
    "survival"
  )
  expect_equal(p$fit, shifted, tolerance = 1e-14)
})

# Standard errors and intervals ------------------------------------------------

test_that("predict() adds interval columns without changing the fit", {
  m <- aft_fit()
  times <- c(3, 12)
  plain <- predict(m, times = times)
  banded <- predict(m, times = times, interval = "confidence")

  expect_identical(names(banded), c(".row", "time", "fit", "lwr", "upr"))
  expect_identical(banded$fit, plain$fit)
  expect_identical(banded[[".row"]], plain[[".row"]])
  expect_true(all(banded$lwr < banded$upr))
})

test_that("predict() returns standard errors beside the frame", {
  m <- aft_fit()
  times <- c(3, 12)
  p <- predict(m, times = times, se.fit = TRUE, interval = "confidence")

  expect_identical(names(p), c("fit", "se.fit"))
  expect_s3_class(p$fit, "data.frame")
  expect_length(p$se.fit, nrow(p$fit))
  # The band is the Wald band on those standard errors.
  cv <- stats::qnorm(0.975)
  expect_equal(p$fit$lwr, p$fit$fit - cv * p$se.fit, tolerance = 1e-14)
  expect_equal(p$fit$upr, p$fit$fit + cv * p$se.fit, tolerance = 1e-14)
})

test_that("predict() widens the interval at a lower level", {
  m <- aft_fit()
  narrow <- predict(m, times = 10, interval = "confidence", level = 0.8)
  wide <- predict(m, times = 10, interval = "confidence", level = 0.99)
  expect_true(all(wide$upr - wide$lwr > narrow$upr - narrow$lwr))
})

test_that("predict() gives each covariate pattern its own standard error", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(
    x = c(-1.5, 0, 1.5),
    g = factor(c("a", "b", "a"), levels = levels(data$g))
  )
  times <- c(4, 16)
  se <- predict(m, newdata = at, times = times, se.fit = TRUE)$se.fit

  # Three patterns at two times give six standard errors, and no pattern
  # borrows another's.
  expect_length(se, 6L)
  by_pattern <- split(se, rep(1:3, each = 2))
  expect_false(isTRUE(all.equal(by_pattern[[1]], by_pattern[[2]])))
  expect_false(isTRUE(all.equal(by_pattern[[2]], by_pattern[[3]])))
})

test_that("predict() gives each row what a one-row call gives it", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(
    x = c(-1.5, 0, 1.5),
    g = factor(c("a", "b", "a"), levels = levels(data$g))
  )
  times <- c(4, 16)
  joint <- predict(m, newdata = at, times = times, interval = "confidence")

  for (i in seq_len(nrow(at))) {
    one <- predict(
      m,
      newdata = at[i, , drop = FALSE],
      times = times,
      interval = "confidence"
    )
    rows <- joint[joint$.row == rownames(at)[i], c("fit", "lwr", "upr")]
    expect_equal(
      unname(as.matrix(rows)),
      unname(as.matrix(one[, c("fit", "lwr", "upr")])),
      tolerance = 1e-14,
      label = paste("row", i)
    )
  }
})

test_that("predict() matches aft_predictions_function() for one pattern", {
  data <- aft_data()
  times <- c(2, 8, 25)

  for (distribution in aft_distributions()) {
    m <- aft_fit(distribution, data)
    pattern <- matrix(c(1, 0.75, 1), nrow = 1)
    at <- data.frame(x = 0.75, g = factor("b", levels = levels(data$g)))
    p <- predict(m, newdata = at, times = times, interval = "confidence")
    expected <- aft_predictions_function(
      X = pattern,
      times = times,
      theta = m@theta,
      covariance = m@variance,
      distribution = distribution,
      measure = "survival"
    )
    expect_equal(
      p$fit,
      expected$predicted,
      tolerance = 1e-14,
      label = distribution
    )
    expect_equal(p$lwr, expected$lower, tolerance = 1e-12, label = distribution)
    expect_equal(p$upr, expected$upper, tolerance = 1e-12, label = distribution)
  }
})

test_that("predict() differentiates exactly when asked", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(x = 0.75, g = factor("b", levels = levels(data$g)))
  times <- c(2, 8, 25)

  p <- predict(
    m,
    newdata = at,
    times = times,
    interval = "confidence",
    deriv_method = "exact"
  )
  expected <- aft_predictions_function(
    X = matrix(c(1, 0.75, 1), nrow = 1),
    times = times,
    theta = m@theta,
    covariance = m@variance,
    distribution = "weibull",
    measure = "survival",
    deriv_method = "exact"
  )
  expect_equal(p$fit, expected$predicted, tolerance = 1e-14)
  expect_equal(p$lwr, expected$lower, tolerance = 1e-12)
  expect_equal(p$upr, expected$upper, tolerance = 1e-12)

  # The finite-difference default agrees with it to the precision of the step.
  approx <- predict(m, newdata = at, times = times, se.fit = TRUE)
  exact <- predict(
    m,
    newdata = at,
    times = times,
    se.fit = TRUE,
    deriv_method = "exact"
  )
  expect_equal(approx$se.fit, exact$se.fit, tolerance = 1e-6)
})

test_that("predict() differentiates exactly across a multi-row design", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  at <- data.frame(
    x = c(-1.5, 0, 1.5),
    g = factor(c("a", "b", "a"), levels = levels(data$g))
  )
  times <- c(4, 16)
  joint <- predict(
    m,
    newdata = at,
    times = times,
    se.fit = TRUE,
    interval = "confidence",
    deriv_method = "exact"
  )

  # A one-row design degenerates every mechanism the exact path uses on a real
  # grid: `X[, j] * th[beta_index[j]]` is a scalar tangent pair rather than a
  # vector of them, `c()` concatenates scalars rather than columns, and the
  # observation-major reordering is the identity permutation. Each row of the
  # joint call is therefore pinned against the one-row call that covers it
  # alone, which is the case those degenerate shapes do cover.
  for (i in seq_len(nrow(at))) {
    one <- predict(
      m,
      newdata = at[i, , drop = FALSE],
      times = times,
      se.fit = TRUE,
      interval = "confidence",
      deriv_method = "exact"
    )
    rows <- joint$fit$.row == rownames(at)[i]
    expect_equal(
      unname(as.matrix(joint$fit[rows, c("fit", "lwr", "upr")])),
      unname(as.matrix(one$fit[, c("fit", "lwr", "upr")])),
      tolerance = 1e-14,
      label = paste("row", i)
    )
    expect_equal(
      unname(joint$se.fit[rows]),
      unname(one$se.fit),
      tolerance = 1e-14,
      label = paste("row", i)
    )
  }

  # And over the whole fitted design, where the grid is 200 rows by two times,
  # against the finite-difference method that has no such shapes to get right.
  exact <- predict(m, times = times, se.fit = TRUE, deriv_method = "exact")
  approx <- predict(m, times = times, se.fit = TRUE)
  expect_length(exact$se.fit, nrow(data) * length(times))
  expect_equal(exact$fit$fit, approx$fit$fit, tolerance = 1e-14)
  expect_equal(exact$se.fit, approx$se.fit, tolerance = 1e-6)
})

test_that("predict() varies the step of the finite-difference methods", {
  m <- aft_fit()
  coarse <- predict(m, times = 10, se.fit = TRUE, dx = 1e-4)$se.fit
  fine <- predict(m, times = 10, se.fit = TRUE, dx = 1e-9)$se.fit
  expect_equal(coarse, fine, tolerance = 1e-5)
  expect_false(identical(coarse, fine))
})

test_that("predict() gives a pooled logistic fit per-row standard errors", {
  data <- plogit_data()
  m <- plogit_fit(data)
  times <- c(12, 24)
  p <- predict(m, times = times, interval = "confidence")

  expect_identical(nrow(p), nrow(data) * 2L)
  expect_true(all(p$lwr < p$upr))
  # Two people whose covariates differ do not share a standard error.
  se <- predict(m, times = times, se.fit = TRUE)$se.fit
  first <- se[1:2]
  other <- se[(which(data$size != data$size[1])[1] - 1) * 2 + 1:2]
  expect_false(isTRUE(all.equal(first, other)))
})

test_that("predict() predicts a pooled logistic fit at new covariates", {
  data <- plogit_data()
  m <- plogit_fit(data)
  at <- data.frame(novel = c(0, 1), init = 2, size = 3)
  times <- c(12, 24)
  p <- predict(m, newdata = at, times = times)

  expect_identical(nrow(p), 4L)
  # The time grid is the fitted sample's, so a new covariate pattern is
  # predicted on it rather than on one of its own.
  expected <- plogit_predict(
    m@theta,
    time = data$time,
    event = data$delta,
    X = cbind(c(0, 1), 2, 3),
    times_to_predict = times,
    measure = "survival"
  )
  expect_equal(p$fit, as.numeric(expected), tolerance = 1e-14)
})

test_that("predict() will not carry an offset into a pooled logistic fit", {
  data <- plogit_data()
  data$shift <- round(seq(-0.1, 0.1, length.out = nrow(data)), 4)
  k <- length(unique(data$time[data$delta == 1]))
  m <- m_estimate(
    time ~ novel + init + size + offset(shift) - 1,
    data = data,
    .ee = ee_plogit,
    event = delta,
    init = c(rep(0, 3), -4, rep(0, k - 1))
  )
  expect_error(predict(m, times = 12), "offset")
})

test_that("predict() refuses an offset a pooled logistic fit took through ...", {
  data <- plogit_data()
  data$shift <- round(seq(-0.1, 0.1, length.out = nrow(data)), 4)
  k <- length(unique(data$time[data$delta == 1]))
  m <- m_estimate(
    time ~ novel + init + size - 1,
    data = data,
    .ee = ee_plogit,
    event = delta,
    offset = shift,
    init = c(rep(0, 3), -4, rep(0, k - 1))
  )
  at <- data.frame(novel = c(0, 1), init = 2, size = 3)

  # An offset supplied through `...` cannot be evaluated on newdata, and the
  # message about that says to write it into the formula instead. For this fit
  # that is a dead end, since the formula offset is refused too, so the refusal
  # that applies is the one reported.
  err <- expect_error(predict(m, newdata = at, times = 12), "offset")
  expect_match(flatten_message(err), "plogit_predict()", fixed = TRUE)
  expect_match(flatten_message(err), "offset of zero", fixed = TRUE)
  expect_false(grepl("evaluated on", flatten_message(err), fixed = TRUE))
})

test_that("predict() cannot differentiate a pooled logistic fit exactly", {
  m <- plogit_fit()
  # Refused whether or not a standard error is asked for: the Jacobian that
  # cannot be built is only built when one is, so the answer would otherwise
  # depend on which other arguments were written.
  expect_error(
    predict(m, times = 12, se.fit = TRUE, deriv_method = "exact"),
    "deriv_method"
  )
  err <- expect_error(
    predict(m, times = 12, deriv_method = "exact"),
    "deriv_method"
  )
  expect_match(flatten_message(err), "ee_plogit()", fixed = TRUE)
  # The finite-difference methods still answer for the same fit.
  expect_no_error(predict(m, times = 12, deriv_method = "fapprox"))
})

# Which equations route and which do not ---------------------------------------

test_that("predict_survival_kind() names only the equations with a surface", {
  expect_identical(predict_survival_kind(ee_aft), "aft")
  expect_identical(predict_survival_kind(ee_plogit), "plogit")
  expect_null(predict_survival_kind(ee_mlogit))
  expect_null(predict_survival_kind(ee_survival_model))
  expect_null(predict_survival_kind(ee_regression))
  expect_null(predict_survival_kind(ee_gformula))
  expect_null(predict_survival_kind("not a function"))
})

test_that("ee_survival_model() cannot be driven by the formula interface", {
  # It takes no design matrix, so the formula interface has nowhere to put the
  # one it builds. There is therefore no fit of it for predict() to route, and
  # survival_predictions() is the way to predict from one.
  data <- aft_data()
  expect_error(
    m_estimate(
      time ~ 1,
      data = data,
      .ee = ee_survival_model,
      event = status,
      distribution = "weibull",
      init = c(0.1, 0.1)
    ),
    class = "deli_formula_ee_signature_error"
  )
})

test_that("an ee_mlogit() fit has no survival surface but does have fits", {
  # Its indicator-matrix response might suggest the formula interface cannot
  # drive it, as it cannot drive `ee_survival_model()`.
  # `coerce_formula_response()` rejects only a factor or character response with
  # other than two levels, and a matrix is neither, so the response passes
  # through and the fit succeeds. The equation is refused for having no survival
  # function rather than for having no fits to refuse.
  set.seed(414)
  n <- 120L
  d <- data.frame(x = round(stats::rnorm(n), 3))
  category <- sample(3L, n, replace = TRUE)
  d$y1 <- as.numeric(category == 1L)
  d$y2 <- as.numeric(category == 2L)
  d$y3 <- as.numeric(category == 3L)

  m <- m_estimate(
    cbind(y1, y2, y3) ~ x,
    data = d,
    .ee = ee_mlogit,
    init = rep(0, 4)
  )
  expect_identical(dim(m@model_spec$y), c(n, 3L))
  expect_error(predict(m, times = 5), "survival measure")
})

test_that("predict() aborts on times for an equation with no survival measure", {
  m <- regression_fit()
  err <- expect_error(predict(m, times = c(1, 2)), "times")
  expect_match(flatten_message(err), "ee_aft()", fixed = TRUE)
  expect_match(flatten_message(err), "ee_plogit()", fixed = TRUE)
})

test_that("predict() aborts on times for a fit with no model specification", {
  y <- c(1, 2, 3, 4, 5)
  m <- m_estimate(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(0)
  )
  expect_error(predict(m, times = 1), "formula interface")
})

test_that("predict() names the times surface when it declines an aft fit", {
  m <- aft_fit()
  err <- expect_error(predict(m), "log-time")
  expect_match(flatten_message(err), "times", fixed = TRUE)
})

test_that("predict() names the times surface when it declines a plogit fit", {
  m <- plogit_fit()
  err <- expect_error(predict(m), "ee_plogit()", fixed = TRUE)
  expect_match(flatten_message(err), "times", fixed = TRUE)
})

# How type and times coexist ---------------------------------------------------

test_that("predict() rejects type beside times", {
  m <- aft_fit()
  err <- expect_error(predict(m, times = 5, type = "response"), "type")
  expect_match(flatten_message(err), "times", fixed = TRUE)
  # The link scale is refused too: `type` names a scale of the linear
  # predictor, and the survival surface has none.
  expect_error(predict(m, times = 5, type = "link"), "type")
})

test_that("predict() rejects the survival arguments without times", {
  m <- regression_fit()
  expect_error(predict(m, measure = "risk"), "measure")
  expect_error(predict(m, deriv_method = "exact"), "deriv_method")
  expect_error(predict(m, dx = 1e-6), "dx")
  err <- expect_error(predict(m, measure = "risk", dx = 1e-6), "measure")
  expect_match(flatten_message(err), "times", fixed = TRUE)
})

test_that("predict() keeps the linear-predictor surface it had", {
  m <- regression_fit()
  expect_identical(
    predict(m, type = "response"),
    predict(m, newdata = mtcars, type = "response")
  )
  expect_length(predict(m), nrow(mtcars))
})

# Argument validation ----------------------------------------------------------

test_that("predict() rejects times that are not usable numbers", {
  m <- aft_fit()
  expect_error(predict(m, times = "10"), "times")
  expect_error(predict(m, times = numeric(0)), "times")
  expect_error(predict(m, times = c(5, NA)), "times")
  expect_error(predict(m, times = c(5, Inf)), "times")
})

test_that("predict() rejects a measure it has no formula for", {
  m <- aft_fit()
  expect_error(predict(m, times = 5, measure = "odds"), "measure")
  expect_error(
    predict(m, times = 5, measure = c("risk", "survival")),
    "measure"
  )
})

test_that("a pooled logistic fit rejects a measure on the same terms", {
  # `check_predict_measure()` only asks for a single string; which measures
  # exist is `convert_survival_measures()`'s to state, and both paths reach it,
  # the AFT one through `aft_measure_at_time()` and this one through
  # `plogit_predict()`. The shape check is common to the two surfaces and the
  # name check is not, so the pooled logistic path is exercised on its own.
  m <- plogit_fit()
  err <- expect_error(predict(m, times = 12, measure = "odds"), "measure")
  expect_match(flatten_message(err), "not supported", fixed = TRUE)
  expect_error(
    predict(m, times = 12, measure = c("risk", "survival")),
    "measure"
  )
})

test_that("predict() rejects a derivative method it has no rule for", {
  m <- aft_fit()
  expect_error(
    predict(m, times = 5, deriv_method = "quadrature"),
    "deriv_method"
  )
})

test_that("predict() rejects a level outside the unit interval on times", {
  m <- aft_fit()
  expect_error(
    predict(m, times = 5, interval = "confidence", level = 0),
    "level"
  )
})

test_that("predict() rejects an argument it does not take beside times", {
  m <- aft_fit()
  expect_error(predict(m, times = 5, mesure = "risk"), "mesure")
})

test_that("predict() does not route a user function that wraps ee_aft", {
  data <- aft_data()
  # Routing is by identity, so a wrapper is not the equation it calls. The
  # wrapper may reorder or reinterpret the parameters, and nothing about it says
  # whether it did.
  weibull_aft <- function(theta, X, time, event) {
    ee_aft(theta, X = X, time = time, event = event, distribution = "weibull")
  }
  m <- m_estimate(
    time ~ x + g,
    data = data,
    .ee = weibull_aft,
    event = status,
    init = c(mean(log(data$time)), 0, 0, 0)
  )
  # The same fit by every number in it, and still not routed.
  expect_equal(
    unname(coef(m)),
    unname(coef(aft_fit("weibull", data))),
    tolerance = 1e-12
  )
  expect_error(predict(m, times = 5), "estimating equation")
})

# The helpers keep their shapes ------------------------------------------------

test_that("aft_predictions_individual() keeps its wide shape and no intervals", {
  data <- aft_data()
  m <- aft_fit("weibull", data)
  times <- c(2, 8, 25)
  wide <- aft_predictions_individual(
    X = m@model_spec$X,
    times = times,
    theta = m@theta,
    distribution = "weibull",
    measure = "survival"
  )

  # Twelve call sites in vignettes/articles/collett-survival.qmd read a one-row
  # call through as.numeric(), so an extra column would lengthen those vectors
  # rather than error. The shape is pinned here directly.
  expect_s3_class(wide, "data.frame")
  expect_identical(dim(wide), c(nrow(data), length(times)))
  expect_identical(names(wide), c("t_2", "t_8", "t_25"))
  expect_true(all(vapply(wide, is.numeric, logical(1))))

  one <- aft_predictions_individual(
    X = m@model_spec$X[1, , drop = FALSE],
    times = times,
    theta = m@theta,
    distribution = "weibull",
    measure = "survival"
  )
  expect_identical(dim(one), c(1L, 3L))
  expect_length(as.numeric(one), 3L)
})

test_that("aft_predictions_function() still takes one covariate pattern", {
  m <- aft_fit()
  expect_error(
    aft_predictions_function(
      X = m@model_spec$X[1:2, ],
      times = c(5, 10),
      theta = m@theta,
      covariance = m@variance,
      distribution = "weibull"
    ),
    "one row"
  )
})
