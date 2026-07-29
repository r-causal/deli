# Formula-interface behavior for m_estimate() / gmm_estimate(): offset() terms,
# factor and character responses, NA-filtered alignment of dots-supplied
# vectors, and the diagnostic raised when the automatically generated init does
# not fit the estimating equation.

# ---- offset() terms ----------------------------------------------

make_offset_data <- function() {
  set.seed(1)
  n <- 40
  d <- data.frame(x = stats::rnorm(n), w = stats::runif(n))
  lp <- -0.5 + 0.6 * d$x + d$w
  d$y <- stats::rbinom(n, 1, 1 / (1 + exp(-lp)))
  d
}

test_that("m_estimate() honors an offset() term in the formula", {
  d <- make_offset_data()
  oracle <- stats::glm(y ~ x, family = stats::binomial, offset = w, data = d)

  m <- m_estimate(
    y ~ x + offset(w),
    data = d,
    .ee = ee_regression,
    model = "logistic"
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-5)
})

test_that("gmm_estimate() honors an offset() term in the formula", {
  d <- make_offset_data()
  oracle <- stats::glm(y ~ x, family = stats::binomial, offset = w, data = d)

  g <- gmm_estimate(
    y ~ x + offset(w),
    data = d,
    .ee = ee_regression,
    model = "logistic"
  )
  expect_equal(unname(coef(g)), unname(coef(oracle)), tolerance = 1e-5)
})

test_that("m_estimate() rejects an offset supplied twice", {
  d <- make_offset_data()
  expect_error(
    m_estimate(
      y ~ x + offset(w),
      data = d,
      .ee = ee_regression,
      model = "logistic",
      offset = w
    ),
    regexp = "both"
  )
})

# ---- factor and character responses -----------------------------

make_binary_data <- function() {
  set.seed(2)
  n <- 60
  d <- data.frame(x = stats::rnorm(n))
  lp <- 0.3 + 0.8 * d$x
  d$y01 <- stats::rbinom(n, 1, 1 / (1 + exp(-lp)))
  d$yf <- factor(ifelse(d$y01 == 1, "yes", "no"))
  d$yc <- ifelse(d$y01 == 1, "yes", "no")
  d
}

test_that("m_estimate() converts a two-level factor response", {
  d <- make_binary_data()
  oracle <- stats::glm(y01 ~ x, family = stats::binomial, data = d)

  m <- m_estimate(yf ~ x, data = d, .ee = ee_regression, model = "logistic")
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-5)
})

test_that("m_estimate() converts a two-level character response", {
  d <- make_binary_data()
  oracle <- stats::glm(y01 ~ x, family = stats::binomial, data = d)

  m <- m_estimate(yc ~ x, data = d, .ee = ee_regression, model = "logistic")
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-5)
})

test_that("m_estimate() rejects a factor response with more than two levels", {
  set.seed(3)
  n <- 60
  d <- data.frame(
    x = stats::rnorm(n),
    y = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_regression, model = "logistic"),
    regexp = "two-level"
  )
})

# ---- NA-filtered alignment of dots vectors ----------------------

test_that("m_estimate() aligns dots vectors with the NA-filtered model frame", {
  set.seed(4)
  n <- 30
  d <- data.frame(x = stats::rnorm(n), w = stats::runif(n))
  d$y <- 1 + 2 * d$x + d$w + stats::rnorm(n, sd = 0.1)
  d$x[5] <- NA

  complete <- d[!is.na(d$x), ]
  oracle <- stats::lm(y ~ x, offset = w, data = complete)

  m <- expect_no_warning(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      offset = w
    )
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
})

test_that("gmm_estimate() aligns dots vectors with the NA-filtered frame", {
  set.seed(5)
  n <- 30
  d <- data.frame(x = stats::rnorm(n), w = stats::runif(n))
  d$y <- 1 + 2 * d$x + d$w + stats::rnorm(n, sd = 0.1)
  d$x[7] <- NA

  complete <- d[!is.na(d$x), ]
  oracle <- stats::lm(y ~ x, offset = w, data = complete)

  g <- expect_no_warning(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      offset = w
    )
  )
  expect_equal(unname(coef(g)), unname(coef(oracle)), tolerance = 1e-6)
})

# ---- auto-init length for extra-parameter EEs -------------------

make_gamma_data <- function() {
  set.seed(6)
  n <- 50
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rgamma(n, shape = 2, rate = exp(-(0.2 + 0.3 * d$x)))
  d
}

test_that("m_estimate() rejects a too-short auto-init informatively", {
  d <- make_gamma_data()
  expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_glm,
      distribution = "gamma",
      link = "log"
    ),
    regexp = "init"
  )
})

test_that("m_estimate() gamma fit succeeds with an explicit init", {
  d <- make_gamma_data()
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(`(Intercept)` = 0, x = 0, dispersion = 1)
  )
  expect_s3_class(m, "deli::MEstimator")
  expect_equal(length(coef(m)), 3L)
})

# ---- the automatic-init diagnostic ------------------------------

# cli wraps a message to the console width, so a phrase that reads as one line
# in the source can arrive with a newline in the middle of it. Collapsing runs
# of whitespace lets a test match the wording without pinning the wrapping.
flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

collect_warnings <- function(expr) {
  seen <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  seen
}

make_line_data <- function() {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(20))
  d$y <- 1 + 2 * d$x + stats::rnorm(20, sd = 0.5)
  d
}

# The score of a linear model, one row per coefficient.
ee_line <- function(theta, X, y, ...) {
  t(X * as.vector(y - X %*% theta))
}

# The same model with the plain residual stacked on as a third moment, which
# leaves the two coefficients over-identified and so needs GMM.
ee_line_over <- function(theta, X, y, ...) {
  r <- as.vector(y - X %*% theta)
  rbind(t(X * r), r)
}

test_that("m_estimate() reports the automatic init length for a gamma fit", {
  d <- make_gamma_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_glm,
      distribution = "gamma",
      link = "log"
    ),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_match(
    flat,
    "the automatic zero `init` of length 2 (the number of model-matrix columns)",
    fixed = TRUE
  )
  # The equation is one this package recognizes, so the parameter the automatic
  # length leaves out is named rather than described.
  expect_match(
    flat,
    "one parameter beyond the design coefficients, \"log_shape\"",
    fixed = TRUE
  )
  expect_match(flat, "Supply an explicit `init` of length 3", fixed = TRUE)
})

test_that("the automatic-init diagnostic names a wrong-shaped return's parameter", {
  # ee_tobit() reads its log scale off an element the automatic init does not
  # reach, which makes every value NA rather than raising, so the return is
  # judged on its shape. The parameter is named on that path too.
  set.seed(11)
  d <- data.frame(x = stats::rnorm(40))
  d$y <- pmax(1 + 0.5 * d$x + stats::rnorm(40), 0)

  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_tobit, lower = 0),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)

  expect_match(
    flat,
    "The automatic zero `init` has length 2, the number of model-matrix columns",
    fixed = TRUE
  )
  expect_match(
    flat,
    "one parameter beyond the design coefficients, \"log_sigma\"",
    fixed = TRUE
  )
})

test_that("the automatic-init diagnostic stays general for an unrecognized equation", {
  # Nothing is known about a custom equation's parameter layout, so the hint
  # describes the common cause instead of naming a parameter.
  d <- make_line_data()
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }

  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)

  expect_match(
    flat,
    "`ee_glm()` with \"gamma\" or \"negative_binomial\" append an extra parameter",
    fixed = TRUE
  )
  expect_match(flat, "Supply an explicit `init`", fixed = TRUE)
})

test_that("the automatic-init diagnostic keeps the original error as its cause", {
  d <- make_line_data()
  ee_broken <- function(theta, X, y, ...) {
    cli::cli_abort("the equation itself gave up")
  }
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_broken))
  flat <- flatten_message(err)
  expect_match(
    flat,
    "Evaluating the estimating function at the automatic zero `init` of length 2",
    fixed = TRUE
  )
  expect_match(flat, "failed", fixed = TRUE)
  expect_match(flat, "A length mismatch is the most common cause", fixed = TRUE)
  expect_match(
    conditionMessage(err$parent),
    "the equation itself gave up",
    fixed = TRUE
  )
})

test_that("gmm_estimate() allows more equations than the automatic init", {
  d <- make_line_data()
  g <- gmm_estimate(y ~ x, data = d, .ee = ee_line_over)
  expect_s3_class(g, "deli::GMMEstimator")
  expect_equal(length(coef(g)), 2L)
})

test_that("gmm_estimate() rejects fewer equations than the automatic init", {
  d <- make_line_data()
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  expect_match(
    flatten_message(err),
    "The automatic zero `init` has length 2, the number of model-matrix columns",
    fixed = TRUE
  )
})

test_that("m_estimate()'s formula interface evaluates the estimating function as often as the function interface", {
  d <- make_line_data()
  evaluations <- 0L
  ee_counted <- function(theta, X, y, ...) {
    evaluations <<- evaluations + 1L
    ee_line(theta, X, y)
  }

  m_formula <- m_estimate(y ~ x, data = d, .ee = ee_counted)
  formula_evaluations <- evaluations

  X <- stats::model.matrix(y ~ x, data = d)
  evaluations <- 0L
  m_direct <- estimate(MEstimator(
    stacked_equations = function(theta) ee_counted(theta, X = X, y = d$y),
    init = stats::setNames(rep(0, ncol(X)), colnames(X))
  ))

  expect_equal(formula_evaluations, evaluations)
  expect_equal(coef(m_formula), coef(m_direct))
})

test_that("gmm_estimate()'s formula interface evaluates the estimating function as often as the function interface", {
  d <- make_line_data()
  evaluations <- 0L
  ee_counted <- function(theta, X, y, ...) {
    evaluations <<- evaluations + 1L
    ee_line_over(theta, X, y)
  }

  g_formula <- gmm_estimate(y ~ x, data = d, .ee = ee_counted)
  formula_evaluations <- evaluations

  X <- stats::model.matrix(y ~ x, data = d)
  evaluations <- 0L
  g_direct <- estimate(GMMEstimator(
    stacked_equations = function(theta) ee_counted(theta, X = X, y = d$y),
    init = stats::setNames(rep(0, ncol(X)), colnames(X))
  ))

  expect_equal(formula_evaluations, evaluations)
  expect_equal(coef(g_formula), coef(g_direct))
})

test_that("a warning signalled at the automatic init is emitted once", {
  d <- make_line_data()
  ee_warns <- function(theta, X, y, ...) {
    if (all(theta == 0)) {
      cli::cli_warn("the starting values are on a boundary")
    }
    ee_line(theta, X, y)
  }
  seen <- collect_warnings(m_estimate(y ~ x, data = d, .ee = ee_warns))
  expect_length(seen, 1L)
  expect_match(conditionMessage(seen[[1]]), "on a boundary", fixed = TRUE)
})

test_that("a NULL return at the automatic init keeps its own message", {
  d <- make_line_data()
  ee_null <- function(theta, X, y, ...) NULL
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_null))
  flat <- flatten_message(err)
  expect_match(
    flat,
    "`stacked_equations` returned \"NULL\" at the initial values",
    fixed = TRUE
  )
  expect_false(grepl("automatic zero", flat, fixed = TRUE))
})

# The row count is what separates the next two tests, and it is the whole
# distinction the diagnostic turns on. ee_line() returns one row per
# coefficient, so this return fits the length-2 automatic init and the only
# thing wrong with it is the value. A return that also has the wrong number of
# rows is reframed instead; see the test below.
test_that("a non-finite return that fits the automatic init keeps its own message", {
  d <- make_line_data()
  ee_non_finite <- function(theta, X, y, ...) {
    vals <- ee_line(theta, X, y)
    vals[1, 1] <- NaN
    vals
  }
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_non_finite))
  flat <- flatten_message(err)
  expect_match(
    flat,
    "`stacked_equations` returned non-finite values at the initial values",
    fixed = TRUE
  )
  expect_false(grepl("automatic zero", flat, fixed = TRUE))
})

# The commonest length mismatch there is: the estimating function appends a
# scale parameter after the coefficients, so at the length-2 automatic init
# theta[3] is NA and the return is both three rows against two parameters and
# full of NA. Reporting the NAs would name a symptom of the short init rather
# than the init itself, so the row count is judged first.
test_that("a non-finite return that does not fit the automatic init reports the length", {
  d <- make_line_data()
  ee_appended <- function(theta, X, y, ...) {
    r <- as.vector(y - X %*% theta[1:2])
    rbind(t(X * r), r^2 - theta[3])
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_appended),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_match(
    flat,
    "The automatic zero `init` has length 2, the number of model-matrix columns",
    fixed = TRUE
  )
  expect_match(flat, "Supply an explicit `init`", fixed = TRUE)
  expect_false(grepl("non-finite", flat, fixed = TRUE))
})

# A short init cannot turn a numeric return into a character one, so the type
# is reported even when the row count is wrong too. These two pin the type
# description as well, and pin that the return's own values stay out of the
# message: the offending object is an estimating-function return, so pasting it
# in would print one entry per observation.
test_that("a non-numeric return that does not fit the automatic init keeps its own message", {
  d <- make_line_data()
  ee_character <- function(theta, X, y, ...) {
    matrix(as.character(y), nrow = 1)
  }
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_character))
  flat <- flatten_message(err)
  expect_match(
    flat,
    "`stacked_equations` must return a numeric vector or matrix at the initial values, not a character matrix.",
    fixed = TRUE
  )
  expect_false(grepl(as.character(d$y[[1]]), flat, fixed = TRUE))
  expect_false(grepl("automatic zero", flat, fixed = TRUE))
})

# is.finite() accepts a character vector and calls every element non-finite, so
# the type has to be judged before the values or this return is reported as
# non-finite. The row count fits here, which isolates that one boundary.
test_that("a non-numeric return that fits the automatic init is not called non-finite", {
  d <- make_line_data()
  ee_character <- function(theta, X, y, ...) {
    matrix(as.character(rep(y, 2)), nrow = 2)
  }
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_character))
  flat <- flatten_message(err)
  expect_match(
    flat,
    "`stacked_equations` must return a numeric vector or matrix at the initial values, not a character matrix.",
    fixed = TRUE
  )
  expect_false(grepl(as.character(d$y[[1]]), flat, fixed = TRUE))
  expect_false(grepl("non-finite", flat, fixed = TRUE))
})

test_that("a failure away from the initial values is not blamed on the automatic init", {
  d <- make_line_data()
  ee_late <- function(theta, X, y, ...) {
    if (!all(theta == 0)) {
      cli::cli_abort("this equation only evaluates at the origin")
    }
    ee_line(theta, X, y)
  }
  err <- expect_error(m_estimate(y ~ x, data = d, .ee = ee_late))
  flat <- flatten_message(err)
  expect_match(flat, "only evaluates at the origin", fixed = TRUE)
  expect_false(grepl("automatic zero", flat, fixed = TRUE))
})

test_that("an explicit init is not described as the automatic one", {
  d <- make_gamma_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_glm,
      distribution = "gamma",
      link = "log",
      init = c(`(Intercept)` = 0, x = 0, dispersion = 1, extra = 0)
    )
  )
  expect_false(grepl("automatic zero", flatten_message(err), fixed = TRUE))
})
