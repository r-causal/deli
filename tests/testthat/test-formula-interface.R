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

# The entry point an abort reports. `conditionCall()` is `NULL` for an error
# raised with `call = NULL`, and rlang maps a method dispatched through
# `UseMethod()` back to its generic, so the head of the reported call is the
# function the caller typed rather than the `.formula` method or any of the
# helpers beneath it.
reported_entry_point <- function(err) {
  call <- conditionCall(err)
  if (is.null(call)) NULL else call[[1]]
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

# The same model with the residual against the squared covariate stacked on as a
# third moment, which leaves the two coefficients over-identified and so needs
# GMM. The plain residual will not do here: the design carries an intercept, so
# its own score row is already the plain residual, and a third row repeating it
# leaves the moment covariance singular and the system over-identified in name
# only.
ee_line_over <- function(theta, X, y, ...) {
  r <- as.vector(y - X %*% theta)
  rbind(t(X * r), r * X[, 2]^2)
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
  expect_match(
    conditionMessage(err$parent),
    "the equation itself gave up",
    fixed = TRUE
  )
})

# An estimating function that fails at the automatic init fails for a reason,
# and the reason is not always the length. The three tests below are the parent
# errors a short automatic `init` really does produce, and the fourth is a
# failure that has nothing to do with the length; only the first three are told
# that a length mismatch is the likely cause.

test_that("the automatic-init diagnostic keeps the length hint for a non-conformable failure", {
  d <- make_line_data()
  # The shape `ee_glm()` with "gamma" fails in: the equation multiplies the
  # design matrix by a parameter vector one longer than the automatic `init`
  # reaches, so the product is non-conformable.
  ee_nonconformable <- function(theta, X, y, ...) {
    r <- as.vector(y - X %*% theta[1:3])
    rbind(t(X * r), r)
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_nonconformable),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "A length mismatch is the most common cause", fixed = TRUE)
  expect_match(
    conditionMessage(err$parent),
    "non-conformable arguments",
    fixed = TRUE
  )
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("the automatic-init diagnostic keeps the length hint for a subscript out of bounds", {
  d <- make_line_data()
  # `[[` past the end of the parameter vector is the other way an equation that
  # needs one more parameter fails outright rather than returning `NA`s.
  ee_out_of_bounds <- function(theta, X, y, ...) {
    r <- as.vector(y - X %*% theta) - theta[[3]]
    rbind(t(X * r), r)
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_out_of_bounds),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "A length mismatch is the most common cause", fixed = TRUE)
  expect_match(
    conditionMessage(err$parent),
    "subscript out of bounds",
    fixed = TRUE
  )
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("the automatic-init diagnostic keeps the length hint for an equation that names the length itself", {
  d <- make_line_data()
  # An equation that checks the length of `theta` before using it says so in
  # its own words, which the hint agrees with rather than contradicting.
  ee_checks_length <- function(theta, X, y, ...) {
    cli::cli_abort("{.arg theta} must have length 3, not {length(theta)}.")
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_checks_length),
    class = "deli_formula_auto_init_error"
  )
  expect_match(
    flatten_message(err),
    "A length mismatch is the most common cause",
    fixed = TRUE
  )
})

test_that("the automatic-init diagnostic claims no length mismatch for an unrelated failure", {
  d <- make_line_data()
  ee_unrelated <- function(theta, X, y, ...) {
    cli::cli_abort("the observations reached the equation out of order")
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_unrelated),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl("length mismatch", flat, fixed = TRUE))
  # Nor the example of an equation that appends a parameter, which is the same
  # claim by illustration.
  expect_false(grepl("negative_binomial", flat, fixed = TRUE))
  expect_match(
    flat,
    "The estimating function itself failed at those starting values",
    fixed = TRUE
  )
  expect_match(
    flat,
    "Supplying an explicit `init` rules the automatic one out",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(err$parent),
    "out of order",
    fixed = TRUE
  )
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("gmm_estimate() reports an unrelated automatic-init failure against itself", {
  d <- make_line_data()
  ee_unrelated <- function(theta, X, y, ...) {
    cli::cli_abort("the observations reached the equation out of order")
  }
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_unrelated),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl("length mismatch", flat, fixed = TRUE))
  expect_match(
    flat,
    "The estimating function itself failed at those starting values",
    fixed = TRUE
  )
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
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

# ---- an under-identified GMM system --------------------------------------
#
# The shortfall branch and the append-a-parameter branch are opposite failures
# that arrive at the same reframing. An estimating equation that estimates a
# parameter beyond the design coefficients returns too few rows because `init`
# is too short, and the fix is a longer `init`. A GMM system with fewer moment
# conditions than parameters returns too few rows because the system is
# under-identified, and a longer `init` makes it worse: GMM needs at least as
# many moment conditions as parameters, so the fix is more conditions or fewer
# parameters. Sending a user with the second failure after the first fix is the
# one direction the hint must not point.

test_that("an under-identified formula GMM fit is not told to lengthen init", {
  d <- make_line_data()
  # One moment condition against the two columns of the model matrix.
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl("one longer", flat, fixed = TRUE))
  expect_false(grepl("append an extra parameter", flat, fixed = TRUE))
  expect_false(grepl("ee_glm", flat, fixed = TRUE))
})

test_that("an under-identified formula GMM fit is told what it is short of", {
  d <- make_line_data()
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  flat <- flatten_message(err)
  # Both counts, each against the thing it counts, so the reader can see which
  # side is short and a bare digit from elsewhere in the message cannot stand in
  # for either.
  expect_match(flat, "1 moment condition", fixed = TRUE)
  expect_match(flat, "2 parameters", fixed = TRUE)
  expect_match(flat, "moment condition", ignore.case = TRUE)
  expect_match(flat, "under-identified", ignore.case = TRUE)
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

test_that("the append-a-parameter branch keeps its own hint", {
  # The opposite failure, unchanged: a short `init` really is the cause here and
  # a longer one really is the fix.
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
    "one parameter beyond the design coefficients",
    fixed = TRUE
  )
  expect_match(flat, "Supply an explicit `init` of length 3", fixed = TRUE)
  expect_false(grepl("under-identified", flat, fixed = TRUE))
  expect_false(grepl("moment condition", flat, fixed = TRUE))
})

test_that("an M-estimation row shortfall is not called under-identified", {
  # M-estimation needs one equation per parameter exactly, so a shortfall there
  # is a length mismatch rather than an identification failure, and it keeps the
  # length hint.
  d <- make_line_data()
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  expect_false(grepl("under-identified", flatten_message(err), fixed = TRUE))
})

test_that("the over-identified formula GMM path is unaffected", {
  d <- make_line_data()
  g <- gmm_estimate(y ~ x, data = d, .ee = ee_line_over)
  expect_s3_class(g, "deli::GMMEstimator")
  expect_equal(length(coef(g)), 2L)
})

test_that("the just-identified formula GMM path is unaffected", {
  d <- make_line_data()
  g <- gmm_estimate(y ~ x, data = d, .ee = ee_line)
  expect_equal(
    unname(coef(g)),
    unname(stats::coef(stats::lm(y ~ x, data = d))),
    tolerance = 1e-6
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

test_that("a warning signaled at the automatic init is emitted once", {
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

# ---- exact-name matching of the dots forwarded to .ee -----------
#
# R matches a supplied name that is a prefix of exactly one formal to that
# formal, and the estimating equations take no `...` of their own for a name to
# fall into. Forwarding `...` with `do.call()` therefore resolved an abbreviation
# or a prefix typo against the equation's arguments: `weight = w` reached
# `ee_regression()`'s `weights` and returned the weighted estimates, and
# `mod = "linear"` supplied its `model`. Both fitted silently, so a misspelling
# changed a reported number with nothing in the output to say so. The names in
# `...` are matched against the arguments of `.ee` exactly instead, and anything
# else is refused before the equation is evaluated.
#
# None of the built-in `ee_*` equations takes `...`, so the exact match applies
# to every one of them. A caller's own `.ee` may take one, and a name that its
# `...` would absorb stays acceptable.

make_weighted_data <- function() {
  set.seed(21)
  n <- 40
  d <- data.frame(x = stats::rnorm(n), w = stats::runif(n, 0.5, 2))
  d$y <- 1 + 2 * d$x + stats::rnorm(n, sd = 0.5)
  d
}

test_that("m_estimate() refuses a dots name that only prefixes an argument of .ee", {
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`weight`", fixed = TRUE)
  # The refusal is this package's, not the unused-argument error base R raises
  # for a name that matches nothing at all.
  expect_false(grepl("unused argument", flat, fixed = TRUE))
  # The name is in the caller's own call, so that is the call to report.
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("a prefix typo does not silently apply weights", {
  # What the refusal is worth: the weighted fit differs from the unweighted one
  # in every coefficient, so a prefix typo reaching `weights` changes every
  # reported estimate and signals nothing. The two fits are pinned as a pair so
  # that a fit which quietly ignored the argument instead would not satisfy
  # this either.
  d <- make_weighted_data()
  unweighted <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_regression,
    model = "linear"
  )
  weighted <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
  expect_false(isTRUE(all.equal(coef(unweighted), coef(weighted))))

  expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
})

test_that("m_estimate() refuses a dots name that prefixes no argument of .ee", {
  # A transposition is a prefix of nothing, so base R already refused it, but as
  # an unused-argument error carrying the whole vector. The refusal is the same
  # one a prefix gets.
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      wieghts = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`wieghts`", fixed = TRUE)
  expect_false(grepl("unused argument", flat, fixed = TRUE))
})

test_that("the refused-argument message names the argument that was meant", {
  d <- make_weighted_data()

  # A prefix says which argument was meant exactly, since R would have matched
  # it there.
  prefix <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  expect_match(flatten_message(prefix), "Did you mean", fixed = TRUE)
  expect_match(flatten_message(prefix), "`weights`", fixed = TRUE)

  # A near miss that is a prefix of nothing is answered from the spelling.
  near <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      wieghts = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  expect_match(flatten_message(near), "Did you mean", fixed = TRUE)
  expect_match(flatten_message(near), "`weights`", fixed = TRUE)

  # A required argument reached by an abbreviation is named the same way.
  abbreviated <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_regression, mod = "linear"),
    class = "deli_formula_ee_argument_error"
  )
  expect_match(flatten_message(abbreviated), "`model`", fixed = TRUE)
})

test_that("the suggestion never names an argument the interface fills itself", {
  # `theta`, `X`, and the response are supplied by the interface, so a caller
  # cannot pass them and a suggestion naming one would send them after an
  # argument they are not allowed to write. Both names below are one edit from
  # an interface-filled argument (`thta` from `theta`, `yy` from `y`) and four
  # or more from every argument a caller may pass, so the refusal carries no
  # suggestion at all.
  d <- make_weighted_data()

  from_theta <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      thta = 1
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat_theta <- flatten_message(from_theta)
  expect_match(flat_theta, "`thta`", fixed = TRUE)
  expect_false(grepl("Did you mean", flat_theta, fixed = TRUE))

  from_response <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      yy = 1
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat_response <- flatten_message(from_response)
  expect_match(flat_response, "`yy`", fixed = TRUE)
  expect_false(grepl("Did you mean", flat_response, fixed = TRUE))
})

test_that("a name resembling no argument is refused without a suggestion", {
  # Every argument of the equation is a long way from this name, and a
  # suggestion nothing supports would send the caller after the wrong argument.
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      favorite_color = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`favorite_color`", fixed = TRUE)
  expect_false(grepl("Did you mean", flat, fixed = TRUE))
})

test_that("the refused-argument message carries no data", {
  # The offending argument is a column of the data, so pasting its value in
  # would print one number per observation. Base R's unused-argument error does
  # exactly that.
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, x = 0),
      wieghts = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl(as.character(d$w[[1]]), flat, fixed = TRUE))
  expect_lt(nchar(flat), 400L)
})

test_that("exact dots names still reach the estimating equation", {
  d <- make_weighted_data()
  oracle <- stats::lm(y ~ x, data = d, weights = w)
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
})

test_that("an .ee taking dots keeps accepting a name of its own choosing", {
  # A caller's equation with a `...` has somewhere to put any name, so the exact
  # match against its arguments does not apply to it.
  d <- make_weighted_data()
  ee_dots <- function(theta, X, y, ...) {
    t(X * as.vector(y - X %*% theta))
  }
  oracle <- stats::lm(y ~ x, data = d)
  m <- expect_no_warning(
    m_estimate(y ~ x, data = d, .ee = ee_dots, favorite_color = "blue")
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
})

test_that("gmm_estimate() refuses a dots name the same way", {
  # Both formula methods build their estimating function through the same
  # helper, so the check reaches this one too.
  d <- make_weighted_data()
  err <- expect_error(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  expect_match(flatten_message(err), "`weights`", fixed = TRUE)
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

# ---- the dots the formula interface cannot forward ---------------
#
# The exact match above answers a name that belongs to no argument. Three
# further ways of writing `...` reach an argument the caller did not mean, and
# all of them turn on the same fact: the response is passed positionally, so
# every slot after `theta` and `X` is already spoken for.
#
# An argument with no name at all is matched by position, into whichever
# argument the response did not take. `.ee = ee_regression` with a bare
# `"linear"` fitted the linear model, because the string landed in `model` after
# the response filled `y`; a bare vector meant for `weights` would have landed
# there instead and fitted something else again, silently either way.
#
# A name the interface fills itself is worse than unmatched. `theta` and `X` are
# passed by name, so supplying either produced R's "matched by multiple actual
# arguments" from inside the estimating function, which the automatic-`init`
# diagnostic then reframed as a length problem. The argument the response is
# passed to is not passed by name, so supplying it displaced the response into
# the next free argument: `y = something` against `ee_regression()` sent the
# response into `weights` and fitted a weighted model whose every coefficient
# differed, with nothing in the output to say so.
#
# All three are refused where the names are read, before the equation is
# evaluated.

test_that("m_estimate() refuses an argument forwarded without a name", {
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_regression, "linear"),
    class = "deli_formula_ee_unnamed_argument"
  )
  flat <- flatten_message(err)
  # The arguments the caller may write are listed, so the report says where the
  # unnamed one was supposed to go.
  expect_match(flat, "`model`", fixed = TRUE)
  expect_match(flat, "`weights`", fixed = TRUE)
  # Neither the arguments the interface fills nor the value itself belongs in a
  # report about a missing name.
  expect_false(grepl("`theta`", flat, fixed = TRUE))
  expect_false(grepl("linear", flat, fixed = TRUE))
  expect_identical(reported_entry_point(err), quote(m_estimate))
  # The general class the other dots refusals carry is on it too, so a caller
  # matching that catches every fault in `...`.
  expect_s3_class(err, "deli_formula_ee_argument_error")
})

test_that("an unnamed argument does not silently fill the next argument", {
  # What the refusal is worth: positional matching made the bare string reach
  # `model` and fit, so a caller who meant it for any other argument got a
  # different model reported as though they had asked for it.
  d <- make_weighted_data()
  named <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")
  expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_regression, "linear"),
    class = "deli_formula_ee_unnamed_argument"
  )
  # The named spelling is the one that goes on working.
  expect_equal(unname(coef(named)), unname(coef(stats::lm(y ~ x, data = d))))
})

test_that("gmm_estimate() refuses an unnamed argument the same way", {
  d <- make_weighted_data()
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_regression, "linear"),
    class = "deli_formula_ee_unnamed_argument"
  )
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

test_that("m_estimate() refuses theta= and X= in the forwarded arguments", {
  d <- make_weighted_data()

  from_theta <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      theta = c(0, 0)
    ),
    class = "deli_formula_ee_reserved_argument"
  )
  flat_theta <- flatten_message(from_theta)
  expect_match(flat_theta, "`theta`", fixed = TRUE)
  # The report is about the argument, not about the starting values, which the
  # automatic-`init` diagnostic blamed instead.
  expect_false(grepl("automatic zero", flat_theta, fixed = TRUE))
  expect_false(grepl("multiple actual arguments", flat_theta, fixed = TRUE))
  expect_identical(reported_entry_point(from_theta), quote(m_estimate))
  expect_s3_class(from_theta, "deli_formula_ee_argument_error")

  from_design <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear", X = 1),
    class = "deli_formula_ee_reserved_argument"
  )
  flat_design <- flatten_message(from_design)
  expect_match(flat_design, "`X`", fixed = TRUE)
  expect_false(grepl("automatic zero", flat_design, fixed = TRUE))
})

test_that("m_estimate() refuses a name that shadows the response argument", {
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      y = d$y
    ),
    class = "deli_formula_ee_reserved_argument"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`y`", fixed = TRUE)
  # The response comes from the formula, which is what makes the name a
  # duplicate rather than a misspelling.
  expect_match(flat, "response", fixed = TRUE)
  expect_false(grepl("Did you mean", flat, fixed = TRUE))
  # The value is a column of the data and has no place in the report.
  expect_false(grepl(as.character(d$y[[1]]), flat, fixed = TRUE))
})

test_that("a response-shadowing name does not silently reach another argument", {
  # The displaced response became the weights, so every coefficient moved. The
  # two fits are pinned as a pair, so a fit that quietly ignored the name would
  # not satisfy this either.
  d <- make_weighted_data()
  plain <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")
  as_weights <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = d$y
  )
  expect_false(isTRUE(all.equal(coef(plain), coef(as_weights))))

  expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      y = d$y
    ),
    class = "deli_formula_ee_reserved_argument"
  )
})

test_that("gmm_estimate() refuses the arguments the interface fills too", {
  d <- make_weighted_data()

  from_theta <- expect_error(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      theta = c(0, 0)
    ),
    class = "deli_formula_ee_reserved_argument"
  )
  expect_identical(reported_entry_point(from_theta), quote(gmm_estimate))

  from_response <- expect_error(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = ee_regression,
      model = "linear",
      y = d$y
    ),
    class = "deli_formula_ee_reserved_argument"
  )
  expect_identical(reported_entry_point(from_response), quote(gmm_estimate))
})

test_that("an .ee taking dots is still refused the arguments the interface fills", {
  # A `...` of the equation's own gives a name somewhere to go, which is why the
  # exact match does not apply to such an equation. It does not give the
  # interface's own arguments anywhere else to go: `theta` still arrives twice,
  # and the response still lands past the argument it was meant for.
  d <- make_weighted_data()
  ee_dots <- function(theta, X, y, ...) {
    t(X * as.vector(y - X %*% theta))
  }
  expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_dots, theta = c(0, 0)),
    class = "deli_formula_ee_reserved_argument"
  )
  expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_dots, y = d$y),
    class = "deli_formula_ee_reserved_argument"
  )
  # The name its own `...` absorbs still reaches it.
  expect_no_error(
    m_estimate(y ~ x, data = d, .ee = ee_dots, favorite_color = "blue")
  )
})

test_that("the response slot stops at the equation's own dots", {
  # The response is passed positionally, and R matches a positional argument
  # only against formals that precede `...`. An argument written after the dots
  # can be filled by name and by nothing else, so reading past them named an
  # argument the interface cannot reach.
  expect_identical(formula_ee_response_slot(c("theta", "X", "y")), "y")
  expect_identical(formula_ee_response_slot(c("theta", "X", "y", "...")), "y")
  expect_null(formula_ee_response_slot(c("theta", "X", "...", "y")))
  expect_null(formula_ee_response_slot(c("theta", "...", "X", "y")))
  expect_null(formula_ee_response_slot(c("theta", "X", "...")))
})

test_that("an equation whose response argument follows its dots is refused", {
  # The response lands in the equation's `...` here and `y` is never filled, so
  # the fit failed with `argument "y" is missing, with no default` wrapped in
  # the automatic-`init` diagnostic, which named the starting values and said
  # nothing about the signature. It is refused before the equation is called.
  d <- make_weighted_data()
  ee_dots_first <- function(theta, X, ..., y) {
    t(X * as.vector(y - X %*% theta))
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_dots_first),
    class = "deli_formula_ee_signature_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "response", fixed = TRUE)
  expect_false(grepl("automatic zero", flat, fixed = TRUE))
})

test_that("an equation whose dots follow its response still fits", {
  # The truncation must not cost the ordinary dots-last signature, where the
  # response reaches the argument before the dots exactly as before.
  d <- make_weighted_data()
  ee_dots_last <- function(theta, X, y, ...) {
    t(X * as.vector(y - X %*% theta))
  }
  oracle <- stats::lm(y ~ x, data = d)
  m <- m_estimate(y ~ x, data = d, .ee = ee_dots_last)
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
})

test_that("the arguments the formula interface fills are still forwarded", {
  # The refusals above must not cost the ordinary route: the response comes from
  # the formula, the design from the formula and `data`, and a legitimate name
  # reaches the equation exactly as before.
  d <- make_weighted_data()
  oracle <- stats::lm(y ~ x, data = d, weights = w)
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
  expect_named(coef(m), c("(Intercept)", "x"))
})

# ---- the .ee signature the formula interface requires ------------
#
# The interface fills four arguments of its own: `theta`, the model matrix it
# passes as `X`, the response it passes positionally, and the offset it takes
# from an `offset()` term in the formula. An equation whose arguments leave any
# of them nowhere to go cannot be driven by a formula at all.
# `ee_survival_model(theta, time, event, distribution)` is the built-in example,
# having no `X`. Reaching the estimating equation anyway produced base R's
# unused-argument error, whose message pastes in the whole offending vector, so
# the report was both uninformative about the cause and one line per observation
# long. The signature is checked before the equation is called instead.
#
# This check and the exact match on the names in `...` are both up-front
# validations in the same helper, and they divide by whose argument is at fault:
# this one covers the arguments the interface itself fills, and that one the
# names the caller supplied.

make_survival_data <- function() {
  set.seed(31)
  n <- 50
  d <- data.frame(time = stats::rexp(n, 0.5))
  d$status <- rep(1, n)
  d
}

test_that("m_estimate() refuses an .ee that takes no design matrix", {
  d <- make_survival_data()
  err <- expect_error(
    m_estimate(
      time ~ 1,
      data = d,
      .ee = ee_survival_model,
      event = status,
      distribution = "weibull",
      init = c(0.1, 0.1)
    ),
    class = "deli_formula_ee_signature_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`X`", fixed = TRUE)
  expect_false(grepl("unused argument", flat, fixed = TRUE))
  # An equation passed as a name can be named, which says which of the
  # arguments in the call is the one to change.
  expect_match(flat, "ee_survival_model", fixed = TRUE)
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("the refused-signature message carries no data", {
  # The argument at fault is the design the interface built, so the message
  # must describe it rather than print it.
  d <- make_survival_data()
  err <- expect_error(
    m_estimate(
      time ~ 1,
      data = d,
      .ee = ee_survival_model,
      event = status,
      distribution = "weibull",
      init = c(0.1, 0.1)
    ),
    class = "deli_formula_ee_signature_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl("c(1, 1", flat, fixed = TRUE))
  expect_lt(nchar(flat), 500L)
})

test_that("the signature check applies to any .ee lacking an X argument", {
  # Nothing about the survival equation is special. Any equation whose
  # arguments leave the design nowhere to go is refused the same way.
  d <- make_weighted_data()
  ee_no_design <- function(theta, y) {
    matrix(y - theta[1], nrow = 1)
  }
  err <- expect_error(
    m_estimate(y ~ 1, data = d, .ee = ee_no_design, init = c(0)),
    class = "deli_formula_ee_signature_error"
  )
  expect_match(flatten_message(err), "`X`", fixed = TRUE)
})

test_that("the signature check covers the argument the response is passed to", {
  # The response is passed positionally, into the first argument that is
  # neither `theta` nor `X`. An equation with no such argument leaves it nowhere
  # to go, and base R reported that by pasting the whole response vector in
  # without even naming an argument.
  d <- make_weighted_data()
  ee_no_response <- function(theta, X) {
    t(X * as.vector(X %*% theta))
  }
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_no_response,
      init = c(`(Intercept)` = 0, x = 0)
    ),
    class = "deli_formula_ee_signature_error"
  )
  flat <- flatten_message(err)
  expect_false(grepl(as.character(d$y[[1]]), flat, fixed = TRUE))
  expect_lt(nchar(flat), 500L)
})

test_that("the signature check covers theta", {
  # `theta` is passed by name like the design is, so an equation without it is
  # refused for the same reason and by the same check.
  d <- make_weighted_data()
  ee_no_theta <- function(X, y) {
    t(X * as.vector(y))
  }
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_no_theta, init = c(0, 0)),
    class = "deli_formula_ee_signature_error"
  )
  expect_match(flatten_message(err), "`theta`", fixed = TRUE)
})

test_that("an offset() term is refused against an .ee that takes no offset", {
  # The offset is the fourth argument the interface fills, and the only one it
  # fills from the formula rather than from the model frame. The caller wrote no
  # `offset` name, so the report has to say the offset came from the formula's
  # `offset()` term and must not read as a misspelling to correct.
  d <- make_weighted_data()
  ee_no_offset <- function(theta, X, y) {
    t(X * as.vector(y - X %*% theta))
  }
  err <- expect_error(
    m_estimate(
      y ~ x + offset(w),
      data = d,
      .ee = ee_no_offset,
      init = c(`(Intercept)` = 0, x = 0)
    )
  )
  expect_true(
    inherits(err, "deli_formula_ee_signature_error") ||
      inherits(err, "deli_formula_ee_argument_error")
  )
  flat <- flatten_message(err)
  expect_match(flat, "`offset`", fixed = TRUE)
  # Named as the formula term it came from, which `{.arg offset}` alone does not
  # say.
  expect_match(flat, "offset()", fixed = TRUE)
  expect_match(flat, "formula", fixed = TRUE)
  expect_false(grepl("Did you mean", flat, fixed = TRUE))
  expect_false(grepl(as.character(d$w[[1]]), flat, fixed = TRUE))
  expect_lt(nchar(flat), 400L)
})

test_that("the signature check precedes the automatic-init diagnostic", {
  # The signature is wrong whatever the starting values are, so the report is
  # about the signature rather than about the automatic length that the failure
  # to evaluate would otherwise be blamed on.
  d <- make_weighted_data()
  ee_no_design <- function(theta, y) {
    matrix(y - theta[1], nrow = 1)
  }
  err <- expect_error(
    m_estimate(y ~ 1, data = d, .ee = ee_no_design),
    class = "deli_formula_ee_signature_error"
  )
  expect_false(grepl("automatic zero", flatten_message(err), fixed = TRUE))
})

test_that("gmm_estimate() refuses an .ee that takes no design matrix", {
  d <- make_survival_data()
  err <- expect_error(
    gmm_estimate(
      time ~ 1,
      data = d,
      .ee = ee_survival_model,
      event = status,
      distribution = "weibull",
      init = c(0.1, 0.1)
    ),
    class = "deli_formula_ee_signature_error"
  )
  expect_match(flatten_message(err), "`X`", fixed = TRUE)
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

test_that("an .ee with both faults is refused by one of the two checks", {
  # The two checks coexist, and neither is skipped because the other applies.
  # Which of them reports first is not the point; that base R's unused-argument
  # error is no longer what the caller sees is.
  d <- make_survival_data()
  ee_no_design <- function(theta, time, event) {
    matrix(time - theta[1], nrow = 1)
  }
  err <- expect_error(
    m_estimate(
      time ~ 1,
      data = d,
      .ee = ee_no_design,
      evnt = status,
      init = c(0.1)
    )
  )
  expect_true(
    inherits(err, "deli_formula_ee_signature_error") ||
      inherits(err, "deli_formula_ee_argument_error")
  )
  expect_false(grepl("unused argument", flatten_message(err), fixed = TRUE))
})

# ---- an .ee written as the name of a function --------------------
#
# `do.call()` takes the name of a function as readily as the function itself, so
# a character `.ee` reaches the estimating equation and has to pass the checks on
# the way. Asking whether `.ee` was a function let it past both of them:
# `.ee = "ee_regression"` with `weight = w` fitted weighted, the one outcome the
# exact match exists to remove. The name is resolved to the function it names
# before anything is checked, so one path covers either spelling, and the rest of
# the interface sees the function too.

test_that("a character .ee is held to the same exact match", {
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = "ee_regression",
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`weight`", fixed = TRUE)
  expect_match(flat, "`weights`", fixed = TRUE)
  # A string is the name of what it resolves to, so the report names it.
  expect_match(flat, "ee_regression", fixed = TRUE)
  expect_identical(reported_entry_point(err), quote(m_estimate))
})

test_that("a character .ee naming a prefix does not silently apply weights", {
  # The pair the exact match is worth, as for the equation passed as a function:
  # a fit that took `weight` for `weights` reports every coefficient of the
  # weighted fit and says nothing.
  d <- make_weighted_data()
  weighted <- m_estimate(
    y ~ x,
    data = d,
    .ee = "ee_regression",
    model = "linear",
    weights = w
  )
  unweighted <- m_estimate(
    y ~ x,
    data = d,
    .ee = "ee_regression",
    model = "linear"
  )
  expect_false(isTRUE(all.equal(coef(weighted), coef(unweighted))))
  expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = "ee_regression",
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
})

test_that("a character .ee taking the names it was given still fits", {
  d <- make_weighted_data()
  oracle <- stats::lm(y ~ x, data = d, weights = w)
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = "ee_regression",
    model = "linear",
    weights = w
  )
  expect_equal(unname(coef(m)), unname(coef(oracle)), tolerance = 1e-6)
})

test_that("gmm_estimate() holds a character .ee to the same match", {
  d <- make_weighted_data()
  err <- expect_error(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = "ee_regression",
      model = "linear",
      weight = w
    ),
    class = "deli_formula_ee_argument_error"
  )
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

test_that("a character .ee naming no function is refused by the interface", {
  # Resolving the name is the interface's own step, so a name that resolves to
  # nothing is the interface's own report, naming the string and the call the
  # caller typed rather than the lookup underneath.
  d <- make_weighted_data()
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = "ee_nope", model = "linear")
  )
  flat <- flatten_message(err)
  expect_match(flat, "ee_nope", fixed = TRUE)
  expect_match(flat, "`.ee`", fixed = TRUE)
  expect_identical(reported_entry_point(err), quote(m_estimate))
  expect_s3_class(err, "deli_formula_ee_lookup_error")
})

test_that("a character .ee must name exactly one function", {
  # `match.fun()` resolves one name, so a character vector of any other length
  # was left as it arrived and reached `do.call()`, which refused it as `'what'
  # must be a function or character string` from inside the estimating function.
  # The automatic-`init` diagnostic wraps whatever fails there, so the report
  # named the starting values and buried the cause. The length is checked where
  # the name is resolved instead.
  d <- make_weighted_data()

  several <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = c("ee_regression", "ee_glm"),
      model = "linear"
    ),
    class = "deli_formula_ee_lookup_error"
  )
  flat_several <- flatten_message(several)
  expect_match(flat_several, "`.ee`", fixed = TRUE)
  expect_match(flat_several, "2", fixed = TRUE)
  expect_false(grepl("automatic zero", flat_several, fixed = TRUE))
  expect_false(grepl("'what' must be", flat_several, fixed = TRUE))
  expect_identical(reported_entry_point(several), quote(m_estimate))

  none <- expect_error(
    m_estimate(y ~ x, data = d, .ee = character(0), model = "linear"),
    class = "deli_formula_ee_lookup_error"
  )
  flat_none <- flatten_message(none)
  expect_match(flat_none, "`.ee`", fixed = TRUE)
  expect_false(grepl("automatic zero", flat_none, fixed = TRUE))
})

test_that("gmm_estimate() holds a character .ee to the same length", {
  d <- make_weighted_data()
  err <- expect_error(
    gmm_estimate(
      y ~ x,
      data = d,
      .ee = c("ee_regression", "ee_glm"),
      model = "linear"
    ),
    class = "deli_formula_ee_lookup_error"
  )
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
})

test_that("an .ee that is neither a function nor a name is refused up front", {
  # The two accepted forms are a function and the name of one. Anything else
  # reached `do.call()` from inside the estimating function, where it failed as
  # `'what' must be a function or character string` and the automatic-`init`
  # diagnostic wrapped it, so the report named the starting values and buried
  # the cause. The type is judged before the interface builds anything.
  d <- make_weighted_data()

  for (entry in list(m_estimate, gmm_estimate)) {
    err <- expect_error(
      entry(y ~ x, data = d, .ee = 42),
      class = "deli_formula_ee_lookup_error"
    )
    flat <- flatten_message(err)
    expect_match(flat, "`.ee`", fixed = TRUE)
    expect_false(grepl("automatic zero", flat, fixed = TRUE))
    expect_false(grepl("'what' must be", flat, fixed = TRUE))
  }

  from_m <- expect_error(m_estimate(y ~ x, data = d, .ee = 42))
  expect_identical(reported_entry_point(from_m), quote(m_estimate))
  from_gmm <- expect_error(gmm_estimate(y ~ x, data = d, .ee = 42))
  expect_identical(reported_entry_point(from_gmm), quote(gmm_estimate))
})

test_that("a character .ee reaches the signature check as well", {
  d <- make_survival_data()
  err <- expect_error(
    m_estimate(
      time ~ 1,
      data = d,
      .ee = "ee_survival_model",
      event = status,
      distribution = "weibull",
      init = c(0.1, 0.1)
    ),
    class = "deli_formula_ee_signature_error"
  )
  expect_match(flatten_message(err), "ee_survival_model", fixed = TRUE)
})

test_that("a character .ee is the equation it names throughout", {
  # Resolving the name serves the rest of the interface too: the parameter
  # `ee_glm()` appends past the design coefficients is recognized from the
  # function, which the string naming it is not, so the labels are the ones a
  # gamma fit driven by the function itself gets.
  d <- make_gamma_data()
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = "ee_glm",
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )
  expect_named(coef(m), c("(Intercept)", "x", "log_shape"))
})

test_that("neither check inspects an .ee whose arguments cannot be read", {
  # `args()` has no argument list to give for a primitive such as `[`, and
  # `formals(NULL)` warns rather than reporting anything. Neither check speaks
  # about a function whose arguments it cannot read, so neither reads them and
  # neither warns: nothing about a primitive is an estimating equation, and
  # calling one fails on its own. Asserted on the checks themselves because the
  # warning is raised on the way into them, before any call they could carry it
  # out of.
  expect_no_warning(check_formula_ee_signature(`[`, ee_args = list()))
  expect_no_warning(check_formula_ee_dots(`[`, ee_args = list(weights = 1)))
})

# ---- an abbreviated distribution= on the formula path ------------
#
# `appended_param_name()` reads `distribution` out of the forwarded `...` by
# name, so a prefix that R matched to `ee_glm()`'s `distribution` fitted the
# gamma model but labeled the coefficients from the equation's own row names,
# `X_1` and `X_2`, instead of the model-matrix headings. Spelling the argument
# in full changed the labels and nothing else, which is the same silence the
# exact match on the forwarded names removes.
#
# One case of the labels degrading remains and is intended. A caller's wrapper
# around `ee_glm()` is a different function, so `appended_param_name()` does not
# recognize it and the labels come from the row names. Nothing the interface can
# read says which distribution such a wrapper fixed.

test_that("m_estimate() refuses an abbreviated distribution= on the formula path", {
  d <- make_gamma_data()
  err <- expect_error(
    m_estimate(
      y ~ x,
      data = d,
      .ee = ee_glm,
      dist = "gamma",
      link = "log",
      init = c(0, 0, 0)
    ),
    class = "deli_formula_ee_argument_error"
  )
  flat <- flatten_message(err)
  expect_match(flat, "`dist`", fixed = TRUE)
  expect_match(flat, "`distribution`", fixed = TRUE)
})

test_that("a fully spelled distribution= names the appended parameter", {
  d <- make_gamma_data()
  m <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )
  expect_named(coef(m), c("(Intercept)", "x", "log_shape"))
})

test_that("a wrapper around ee_glm() keeps the equation's own row names", {
  d <- make_gamma_data()
  wrapped_gamma <- function(...) ee_glm(..., distribution = "gamma")

  wrapped <- expect_no_warning(
    m_estimate(
      y ~ x,
      data = d,
      .ee = wrapped_gamma,
      link = "log",
      init = c(0, 0, 0)
    )
  )
  spelled <- m_estimate(
    y ~ x,
    data = d,
    .ee = ee_glm,
    distribution = "gamma",
    link = "log",
    init = c(0, 0, 0)
  )

  expect_named(coef(wrapped), c("X_1", "X_2", "log_shape"))
  # The labels are all that the wrapper costs.
  expect_equal(unname(coef(wrapped)), unname(coef(spelled)))
})

# ---- the call the formula aborts report --------------------------
#
# The automatic-`init` diagnostic is raised where the estimating function is
# first evaluated, several frames below the method the caller reached, and it
# named no call at all, so it printed a bare `Error:` with nothing in it that
# appears in the caller's code. Every abort the formula interface raises names
# the entry point instead: `m_estimate()` on one path and `gmm_estimate()` on
# the other, whichever internal frame the abort is raised in.

test_that("the automatic-init diagnostic names m_estimate()", {
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
  expect_identical(reported_entry_point(err), quote(m_estimate))
  # Neither the method the generic dispatched to nor the helper that evaluates
  # the estimating function is a call the caller can act on.
  reported <- reported_call(err)
  expect_false(grepl("m_estimate.formula", reported, fixed = TRUE))
  expect_false(grepl("eval_psi_at_init", reported, fixed = TRUE))
})

test_that("the automatic-init diagnostic names gmm_estimate()", {
  d <- make_line_data()
  ee_short <- function(theta, X, y, ...) {
    matrix(as.vector(y - X %*% theta), nrow = 1)
  }
  err <- expect_error(
    gmm_estimate(y ~ x, data = d, .ee = ee_short),
    class = "deli_formula_auto_init_error"
  )
  expect_identical(reported_entry_point(err), quote(gmm_estimate))
  reported <- reported_call(err)
  expect_false(grepl("gmm_estimate.formula", reported, fixed = TRUE))
})
