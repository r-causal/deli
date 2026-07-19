# Tests for deli_spline() and additive_design_matrix()

# --- deli_spline() ---

test_that("unrestricted cubic spline returns correct dimensions", {
  x <- seq(-3, 3, length.out = 50)
  knots <- c(-1, 0, 1)
  result <- deli_spline(x, knots, power = 3, restricted = FALSE)

  # Unrestricted: same number of columns as knots

  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 3)
})

test_that("restricted cubic spline returns one fewer column than knots", {
  x <- seq(-3, 3, length.out = 50)
  knots <- c(-1, 0, 1)
  result <- deli_spline(x, knots, power = 3, restricted = TRUE)

  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 2)
})

test_that("unrestricted spline values are correct", {
  x <- c(-2, -0.5, 0, 0.5, 2)
  knots <- c(-1, 0, 1)

  result <- deli_spline(x, knots, power = 3, restricted = FALSE)

  # x = -2: all knots are greater than x, so all zero
  expect_equal(result[1, ], c(0, 0, 0))

  # x = -0.5: only greater than knot -1, so ((-0.5) - (-1))^3 = 0.5^3 = 0.125

  expect_equal(result[2, ], c(0.125, 0, 0))

  # x = 2: greater than all knots
  # (2 - (-1))^3 = 27, (2 - 0)^3 = 8, (2 - 1)^3 = 1
  expect_equal(result[5, ], c(27, 8, 1))
})

test_that("restricted spline values are correct", {
  x <- c(-2, -0.5, 0, 0.5, 2)
  knots <- c(-1, 0, 1)

  result <- deli_spline(x, knots, power = 3, restricted = TRUE)

  # x = -2: all zero (no knots exceeded)
  expect_equal(result[1, ], c(0, 0))

  # x = -0.5: s1 = 0.125, sK (knot=1) = 0 -> restricted = 0.125 - 0 = 0.125
  expect_equal(result[2, ], c(0.125, 0))

  # x = 2: s1 = 27, s2 = 8, sK = 1
  # r1 = s1 - sK = 27 - 1 = 26
  # r2 = s2 - sK = 8 - 1 = 7
  expect_equal(result[5, ], c(26, 7))
})

test_that("spline with power = 2 produces quadratic terms", {
  x <- c(0, 1.5, 3)
  knots <- c(1, 2)
  result <- deli_spline(x, knots, power = 2, restricted = FALSE)

  # x = 0: below both knots
  expect_equal(result[1, ], c(0, 0))

  # x = 1.5: (1.5 - 1)^2 = 0.25, below knot 2
  expect_equal(result[2, ], c(0.25, 0))

  # x = 3: (3 - 1)^2 = 4, (3 - 2)^2 = 1
  expect_equal(result[3, ], c(4, 1))
})

test_that("normalized spline divides by range^power", {
  x <- c(0, 3)
  knots <- c(1, 2)
  unnorm <- deli_spline(
    x,
    knots,
    power = 3,
    restricted = FALSE,
    normalized = FALSE
  )
  normed <- deli_spline(
    x,
    knots,
    power = 3,
    restricted = FALSE,
    normalized = TRUE
  )

  # Divisor = (2 - 1)^3 = 1 in this case
  divisor <- (2 - 1)^3
  expect_equal(normed, unnorm / divisor)

  # Now with wider knot range
  knots2 <- c(0, 3)
  unnorm2 <- deli_spline(
    x,
    knots2,
    power = 3,
    restricted = FALSE,
    normalized = FALSE
  )
  normed2 <- deli_spline(
    x,
    knots2,
    power = 3,
    restricted = FALSE,
    normalized = TRUE
  )
  divisor2 <- (3 - 0)^3
  expect_equal(normed2, unnorm2 / divisor2)
})

test_that("spline with single knot and normalized works", {
  x <- c(0, 2, 4)
  knots <- c(1)
  result <- deli_spline(
    x,
    knots,
    power = 3,
    restricted = FALSE,
    normalized = TRUE
  )
  # Divisor = 1^3 = 1
  expect_equal(ncol(result), 1)
  # x = 0: below knot
  expect_equal(result[1, 1], 0)
  # x = 2: (2-1)^3 / 1 = 1
  expect_equal(result[2, 1], 1)
})

test_that("single-knot normalization divides by knot^power", {
  # With a single knot the range of knots is zero, so the normalization
  # divisor is knots^power (here 2^2 = 4) rather than the range.
  result <- deli_spline(
    c(1, 3, 5),
    knots = 2,
    power = 2,
    restricted = FALSE,
    normalized = TRUE
  )
  expect_equal(result[, 1], c(0, 0.25, 2.25))
})

test_that("spline handles NA values", {
  x <- c(1, NA, 3)
  knots <- c(2)
  result <- deli_spline(x, knots, power = 3, restricted = FALSE)

  expect_equal(result[1, 1], 0)
  expect_true(is.na(result[2, 1]))
  expect_equal(result[3, 1], 1)
})

test_that("spline knots are sorted internally", {
  x <- c(0, 2, 4)
  knots_sorted <- c(1, 3)
  knots_unsorted <- c(3, 1)

  r1 <- deli_spline(x, knots_sorted, power = 3, restricted = FALSE)
  r2 <- deli_spline(x, knots_unsorted, power = 3, restricted = FALSE)
  expect_equal(r1, r2)
})


# --- additive_design_matrix() ---

test_that("additive_design_matrix with all NULL specs returns X unchanged", {
  X <- matrix(c(1, 1, 1, 2, 3, 4), ncol = 2)
  specs <- list(NULL, NULL)

  result <- additive_design_matrix(X, specs)
  expect_equal(result, X)
})

test_that("additive_design_matrix adds spline columns", {
  n <- 20
  X <- cbind(rep(1, n), seq(-2, 2, length.out = n))
  specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))

  result <- additive_design_matrix(X, specs)

  # Original 2 columns + 2 restricted spline columns (3 knots - 1)
  expect_equal(ncol(result), 4)
  expect_equal(nrow(result), n)

  # First two columns should match X
  expect_equal(result[, 1], X[, 1])
  expect_equal(result[, 2], X[, 2])
})

test_that("additive_design_matrix returns penalty vector", {
  n <- 20
  X <- cbind(rep(1, n), seq(-2, 2, length.out = n))
  specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))

  result <- additive_design_matrix(X, specs, return_penalty = TRUE)

  expect_true(is.list(result))
  expect_equal(ncol(result$X), 4)
  # Penalties: 0 (intercept), 0 (linear X), 5, 5 (two spline terms)
  expect_equal(result$penalty, c(0, 0, 5, 5))
})

test_that("additive_design_matrix with unrestricted spline", {
  n <- 10
  X <- cbind(rep(1, n), seq(0, 1, length.out = n))
  knots <- c(0.25, 0.5, 0.75)
  specs <- list(NULL, list(knots = knots, natural = FALSE, penalty = 1))

  result <- additive_design_matrix(X, specs, return_penalty = TRUE)

  # 2 original + 3 unrestricted spline columns
  expect_equal(ncol(result$X), 5)
  # Penalties: 0, 0, 1, 1, 1
  expect_equal(result$penalty, c(0, 0, 1, 1, 1))
})

test_that("additive_design_matrix errors on mismatched spec length", {
  X <- matrix(1:6, ncol = 2)
  specs <- list(NULL)

  expect_error(
    additive_design_matrix(X, specs),
    "specifications"
  )
})

test_that("additive_design_matrix with multiple spline columns", {
  n <- 30
  X <- cbind(rep(1, n), rnorm(n), rnorm(n))
  specs <- list(
    NULL,
    list(knots = c(-1, 0, 1), penalty = 10),
    list(knots = c(-0.5, 0.5), penalty = 5)
  )

  result <- additive_design_matrix(X, specs, return_penalty = TRUE)

  # 3 original + 2 (restricted, 3 knots) + 1 (restricted, 2 knots)
  expect_equal(ncol(result$X), 6)
  expect_equal(result$penalty, c(0, 0, 10, 10, 0, 5))
})

test_that("additive_design_matrix with NULL specifications returns X", {
  X <- matrix(1:9, ncol = 3)
  result <- additive_design_matrix(X, specifications = NULL)
  expect_equal(result, X)
})

test_that("process_spline_spec fills defaults", {
  spec <- list(knots = c(1, 2, 3))
  processed <- process_spline_spec(spec)

  expect_equal(processed$natural, TRUE)
  expect_equal(processed$power, 3)
  expect_equal(processed$penalty, 0)
  expect_equal(processed$normalized, FALSE)
})

test_that("process_spline_spec errors when knots missing", {
  expect_error(
    process_spline_spec(list(penalty = 5)),
    "knots"
  )
})

test_that("process_spline_spec warns on unexpected keys", {
  expect_warning(
    process_spline_spec(list(knots = c(1), foo = "bar")),
    "foo"
  )
})
