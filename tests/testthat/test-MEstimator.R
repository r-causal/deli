# ---- MEstimator construction -------------------------------------------------

test_that("MEstimator constructs with valid inputs", {
  psi <- function(theta) {
    y <- c(1, 2, 3, 4, 5)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_s3_class(m, "deli::MEstimator")
  expect_equal(m@init, 0)
  expect_null(m@theta)
  expect_null(m@variance)
})

test_that("MEstimator stores stacked_equations as function", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_type(m@stacked_equations, "closure")
})

test_that("MEstimator stores init vector", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  expect_equal(m@init, c(0, 1))
})

test_that("MEstimator accepts subset parameter", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 1), subset = c(1L))
  expect_equal(m@subset, 1L)
})

test_that("MEstimator subset is sorted", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2], y * theta[3])
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(0, 1, 1),
    subset = c(3L, 1L)
  )
  expect_equal(m@subset, c(1L, 3L))
})

test_that("MEstimator subset defaults to NULL", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@subset)
})

test_that("MEstimator accepts finite_correction parameter", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(0),
    finite_correction = "HC1"
  )
  expect_equal(m@finite_correction, "HC1")
})

test_that("MEstimator finite_correction defaults to NULL", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@finite_correction)
})

# ---- MEstimator input validation ---------------------------------------------

test_that("MEstimator requires stacked_equations to be a function", {
  expect_error(
    MEstimator(stacked_equations = "not a function", init = c(0)),
    "function"
  )
})

test_that("MEstimator requires init to be numeric", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    MEstimator(stacked_equations = psi, init = "not numeric"),
    "numeric|double"
  )
})

# ---- MEstimator pre-estimation state -----------------------------------------

test_that("MEstimator theta is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@theta)
})

test_that("MEstimator variance is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@variance)
})

test_that("MEstimator bread is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@bread)
})

test_that("MEstimator meat is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@meat)
})

test_that("MEstimator n_params is set from init length", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  expect_equal(m@n_params, 2L)
})

test_that("MEstimator n_obs is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_null(m@n_obs)
})

# ---- MEstimator constructor validation ---------------------------------------
#
# These tests pin construction-time validation of the finite_correction, init,
# and subset arguments. The supported finite_correction values are NULL and the
# HC1 string (see finite_sample_correction() in R/sandwich.R). init must contain
# at least one value. subset must be whole-number parameter indices between 1
# and the number of parameters (1-based indexing), which is knowable at
# construction from length(init). The validators these tests target are added in
# bd-286o.4.

test_that("MEstimator rejects an unsupported finite_correction string", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  # The message must name the supported set, of which HC1 is the only non-NULL
  # option.
  expect_error(
    MEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = "nonsense"
    ),
    "HC1"
  )
})

test_that("MEstimator rejects a non-character finite_correction", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    MEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = 5
    ),
    "HC1"
  )
})

test_that("MEstimator accepts the supported finite_correction values", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_no_error(
    MEstimator(stacked_equations = psi, init = c(0), finite_correction = "HC1")
  )
  expect_no_error(
    MEstimator(stacked_equations = psi, init = c(0), finite_correction = NULL)
  )
})

test_that("MEstimator rejects zero-length init", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  # Target message: init must contain at least one value.
  expect_error(
    MEstimator(stacked_equations = psi, init = numeric(0)),
    "at least one"
  )
})

test_that("MEstimator rejects a subset index outside the parameter range", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  # Only two parameters, so index 5 is out of range.
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = 5L),
    "subset"
  )
})

test_that("MEstimator rejects a zero subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  # Indices are 1-based, so 0 is invalid.
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = 0L),
    "subset"
  )
})

test_that("MEstimator rejects a negative subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = -1L),
    "subset"
  )
})

test_that("MEstimator rejects a non-numeric subset", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = "a"),
    "subset"
  )
})

test_that("MEstimator rejects a fractional subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  # Fractional indices are not valid parameter positions.
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = 1.5),
    "subset"
  )
})

test_that("MEstimator accepts a valid in-range subset", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_no_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = c(1L, 2L))
  )
})

test_that("MEstimator rejects a subset with duplicated indices", {
  # Duplicated indices make the subsetted system rank-deficient, so the solver
  # returns the starting values as if they were estimates. Reject them at
  # construction, naming the duplicates.
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    MEstimator(stacked_equations = psi, init = c(0, 1), subset = c(1L, 1L)),
    "duplicated"
  )
})
