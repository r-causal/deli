# ---- huber loss --------------------------------------------------------------

test_that("huber loss clips residuals at k", {
  r <- c(-3, -1.345, -1, 0, 0.5, 1.345, 2, 5)
  k <- 1.345
  expected <- ifelse(abs(r) <= k, r, k * sign(r))
  result <- robust_loss_functions(r, loss = "huber", k = k)
  expect_equal(result, expected)
})

# ---- tukey loss --------------------------------------------------------------

test_that("tukey loss applies biweight function", {
  r <- c(-6, -4.685, -2, 0, 1.5, 4.685, 5)
  k <- 4.685
  expected <- ifelse(abs(r) <= k, r * (1 - (r / k)^2)^2, 0)
  result <- robust_loss_functions(r, loss = "tukey", k = k)
  expect_equal(result, expected)
})

# ---- andrew loss -------------------------------------------------------------

test_that("andrew loss applies sine function", {
  k <- 1.339
  r <- c(-k * pi - 0.1, -k * pi, -2, 0, 1, k * pi, k * pi + 0.1)
  expected <- ifelse(abs(r) <= k * pi, sin(r / k), 0)
  result <- robust_loss_functions(r, loss = "andrew", k = k)
  expect_equal(result, expected)
})

# ---- hampel loss -------------------------------------------------------------

test_that("hampel loss applies three-part redescending function", {
  a <- 2
  b <- 4
  cc <- 8
  r <- c(-10, -8, -6, -3, -1, 0, 1, 3, 6, 8, 10)
  expected <- ifelse(
    abs(r) < a,
    r,
    ifelse(
      abs(r) < b,
      a * sign(r),
      ifelse(
        abs(r) < cc,
        a * sign(r) * (cc - abs(r)) / (cc - b),
        0
      )
    )
  )
  result <- robust_loss_functions(r, loss = "hampel", k = c(a, b, cc))
  expect_equal(result, expected)
})

test_that("hampel requires k of length 3", {
  expect_error(
    robust_loss_functions(c(1, 2), loss = "hampel", k = 2),
    "length-3"
  )
})

test_that("hampel requires a < b < c", {
  expect_error(
    robust_loss_functions(c(1, 2), loss = "hampel", k = c(4, 2, 8)),
    "a < b < c"
  )
})

# ---- fair loss ---------------------------------------------------------------

test_that("fair loss computes correctly", {
  r <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  k <- 1.4
  expected <- r / (1 + abs(r) / k)
  result <- robust_loss_functions(r, loss = "fair", k = k)
  expect_equal(result, expected)
})

# ---- cauchy loss -------------------------------------------------------------

test_that("cauchy loss computes correctly", {
  r <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  k <- 2.385
  expected <- r / (1 + (r / k)^2)
  result <- robust_loss_functions(r, loss = "cauchy", k = k)
  expect_equal(result, expected)
})

# ---- ullah loss --------------------------------------------------------------

test_that("ullah loss computes correctly", {
  r <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  k <- 1.5
  # Python: r * (1 + (r/k)^4)^(-2)
  expected <- r * (1 + (r / k)^4)^(-2)
  result <- robust_loss_functions(r, loss = "ullah", k = k)
  expect_equal(result, expected)
})

# ---- welsch loss -------------------------------------------------------------

test_that("welsch loss computes correctly", {
  r <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  k <- 2.11
  # Python: r * exp(-r^2 / (2*k^2))
  expected <- r * exp(-r^2 / (2 * k^2))
  result <- robust_loss_functions(r, loss = "welsch", k = k)
  expect_equal(result, expected)
})

# ---- zero residuals ----------------------------------------------------------

test_that("all loss types return 0 when residuals are 0", {
  r <- c(0, 0, 0)
  expect_equal(robust_loss_functions(r, "huber", k = 1.345), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "tukey", k = 4.685), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "andrew", k = 1.339), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "fair", k = 1.4), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "cauchy", k = 2.385), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "ullah", k = 1.5), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "welsch", k = 2.11), c(0, 0, 0))
  expect_equal(robust_loss_functions(r, "hampel", k = c(2, 4, 8)), c(0, 0, 0))
})

# ---- vector input ------------------------------------------------------------

test_that("all loss types handle vector input correctly", {
  r <- c(-10, -5, -1, 0, 1, 5, 10)
  expect_length(robust_loss_functions(r, "huber", k = 1.345), length(r))
  expect_length(robust_loss_functions(r, "tukey", k = 4.685), length(r))
  expect_length(robust_loss_functions(r, "andrew", k = 1.339), length(r))
  expect_length(robust_loss_functions(r, "fair", k = 1.4), length(r))
  expect_length(robust_loss_functions(r, "cauchy", k = 2.385), length(r))
  expect_length(robust_loss_functions(r, "ullah", k = 1.5), length(r))
  expect_length(robust_loss_functions(r, "welsch", k = 2.11), length(r))
  expect_length(robust_loss_functions(r, "hampel", k = c(2, 4, 8)), length(r))
})

# ---- invalid loss type -------------------------------------------------------

test_that("invalid loss type raises an error", {
  expect_error(robust_loss_functions(c(1, 2, 3), loss = "invalid", k = 1))
})

test_that("non-string loss raises an error", {
  expect_error(robust_loss_functions(c(1, 2, 3), loss = 42, k = 1))
})

# ---- symmetry ----------------------------------------------------------------

test_that("all loss types are odd functions (antisymmetric)", {
  r <- c(0.5, 1.5, 3, 7)
  losses <- list(
    list("huber", 1.345),
    list("tukey", 4.685),
    list("andrew", 1.339),
    list("fair", 1.4),
    list("cauchy", 2.385),
    list("ullah", 1.5),
    list("welsch", 2.11),
    list("hampel", c(2, 4, 8))
  )
  for (l in losses) {
    pos <- robust_loss_functions(r, loss = l[[1]], k = l[[2]])
    neg <- robust_loss_functions(-r, loss = l[[1]], k = l[[2]])
    expect_equal(pos, -neg, info = paste("loss =", l[[1]]))
  }
})

# ---- case insensitivity ------------------------------------------------------

test_that("loss name is case insensitive", {
  r <- c(-1, 0, 1)
  expect_equal(
    robust_loss_functions(r, "Huber", k = 1.345),
    robust_loss_functions(r, "huber", k = 1.345)
  )
  expect_equal(
    robust_loss_functions(r, "TUKEY", k = 4.685),
    robust_loss_functions(r, "tukey", k = 4.685)
  )
})
