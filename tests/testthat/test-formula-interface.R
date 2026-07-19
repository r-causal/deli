# Formula-interface behavior for m_estimate() / gmm_estimate(): offset() terms,
# factor and character responses, NA-filtered alignment of dots-supplied
# vectors, and auto-init length checking for extra-parameter estimating
# equations.

# ---- offset() terms (bd-2x8.3) ----------------------------------------------

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

# ---- factor and character responses (bd-2x8.20) -----------------------------

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

# ---- NA-filtered alignment of dots vectors (bd-2x8.21) ----------------------

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

# ---- auto-init length for extra-parameter EEs (bd-2x8.32) -------------------

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
