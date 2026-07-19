# ---- NULL weights returns vector of ones -------------------------------------
test_that("NULL weights returns rep(1, n) for various n", {
  expect_equal(generate_weights(5, NULL), rep(1, 5))
  expect_equal(generate_weights(10, NULL), rep(1, 10))
  expect_equal(generate_weights(100, NULL), rep(1, 100))
})

# ---- provided weights are returned unchanged --------------------------------
test_that("provided weights are returned unchanged", {
  w <- c(1, 2, 3, 4, 5)
  expect_equal(generate_weights(5, w), w)

  w2 <- c(10, 20, 30)
  expect_equal(generate_weights(3, w2), w2)
})

# ---- length mismatch raises error -------------------------------------------
test_that("length mismatch raises error", {
  expect_error(
    generate_weights(5, c(1, 2, 3)),
    "length.*weights.*must.*match|weights.*length"
  )
  expect_error(
    generate_weights(3, c(1, 2, 3, 4, 5)),
    "length.*weights.*must.*match|weights.*length"
  )
})

# ---- works with non-integer weights -----------------------------------------
test_that("works with non-integer weights (e.g., probability weights)", {
  w <- c(0.1, 0.5, 0.9, 0.3)
  expect_equal(generate_weights(4, w), w)

  w2 <- c(1.5, 2.7, 0.01)
  expect_equal(generate_weights(3, w2), w2)
})

# ---- works with n=1 ---------------------------------------------------------
test_that("works with n = 1", {
  expect_equal(generate_weights(1, NULL), 1)
  expect_equal(generate_weights(1, 5), 5)
  expect_equal(generate_weights(1, 0.5), 0.5)
})

# ---- all-zero weights are accepted -------------------------------------------
test_that("all-zero weights are accepted (no validation on values)", {
  w <- rep(0, 5)
  expect_equal(generate_weights(5, w), w)
})

# ---- 2D weight matrices pass through -----------------------------------------
# ee_plogit accepts an n-by-K time-varying weight matrix (one column per time
# interval). generate_weights passes such a matrix through with its dimensions
# preserved, validating only the row count (nrow == n) since it alone knows the
# observation count; the column-count check lives downstream in ee_plogit.
# Python's generate_weights performs no validation at all (a bare np.asarray),
# so the row check here is a deliberate strictness addition consistent with the
# parity policy: it catches a mis-shaped matrix at the earliest point that n is
# known. These tests verify the pass-through, the dimension preservation, and
# the row-count rejection.
test_that("an n-by-K weight matrix passes through with dimensions preserved", {
  n <- 6
  k <- 4
  w <- matrix(seq_len(n * k) / 10, nrow = n, ncol = k)

  result <- generate_weights(n, w)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n, k))
  expect_equal(unname(result), unname(w))
})

test_that("a single-column weight matrix (n-by-1) passes through", {
  n <- 5
  w <- matrix(c(0.5, 1, 1.5, 2, 2.5), nrow = n, ncol = 1)

  result <- generate_weights(n, w)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n, 1))
  expect_equal(unname(result), unname(w))
})

test_that("a weight matrix with the wrong number of rows is rejected", {
  # Matrix rows must equal the observation count. The message must speak to
  # rows, not overall length: a plain length check would wrongly accept an
  # (n-1)-by-(K+1) matrix whose element count happens to equal n.
  n <- 6
  k <- 4
  w_bad <- matrix(seq_len((n - 1) * k) / 10, nrow = n - 1, ncol = k)

  expect_error(generate_weights(n, w_bad), "row")
})

# ---- tangent-carrying weights pass through -----------------------------------
# Under deriv_method = "exact", ee_ipw_msm forms IPW weights from the propensity
# model (theta-derived, so tangent-carrying) and passes them as the weights
# argument to a weighted GLM. generate_weights must preserve the primal-tangent
# structure so the downstream weighted score keeps its derivative; the plain
# as.numeric() coercion strips the tangent and aborts. These tests verify the
# pass-through, that a subsequent w * score keeps the tangent, and that the
# length check still applies (validated on the primal).
test_that("a tangent-carrying weights vector passes through without stripping the class", {
  primal <- c(1.5, 2.0, 0.5, 3.0)
  tangent <- c(0.1, -0.2, 0.05, 0.4)
  w_pt <- primal_tangent_array(primal, tangent)

  result <- generate_weights(4, w_pt)
  expect_true(is_pt_array(result))
  expect_equal(result$primal, primal)
  expect_equal(result$tangent, tangent)

  # w * score is exactly the operation ee_glm applies to the weights; the
  # tangent must survive so the weighted score differentiates.
  scored <- result * c(2, 2, 2, 2)
  expect_true(is_pt_array(scored))
  expect_equal(scored$tangent, tangent * 2)
})

test_that("a tangent-carrying weights vector of the wrong length is rejected", {
  w_pt <- primal_tangent_array(c(1, 2, 3), c(0.1, 0.2, 0.3))
  expect_error(generate_weights(5, w_pt), "length")
})
