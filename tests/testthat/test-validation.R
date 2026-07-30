# ---- Helpers ----

# cli wraps a message at the console width, so a type description can straddle
# a line break. Collapsing the whitespace keeps the assertions on it independent
# of where the break falls.
flatten_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

# check_alpha -------------------------------------------------------------

test_that("check_alpha accepts valid alpha values", {
  expect_no_error(check_alpha(0.05))
  expect_no_error(check_alpha(0.5))
  expect_no_error(check_alpha(0.99))
  expect_no_error(check_alpha(0.01))
  expect_no_error(check_alpha(0.001))
})

test_that("check_alpha rejects alpha outside (0, 1)", {
  expect_error(check_alpha(0), "alpha.*0.*1")
  expect_error(check_alpha(1), "alpha.*0.*1")
  expect_error(check_alpha(-0.5), "alpha.*0.*1")
  expect_error(check_alpha(1.5), "alpha.*0.*1")
  expect_error(check_alpha(-1), "alpha.*0.*1")
  expect_error(check_alpha(2), "alpha.*0.*1")
})

test_that("check_alpha rejects boundary values exactly at 0 and 1", {
  expect_error(check_alpha(0), "alpha")
  expect_error(check_alpha(1), "alpha")
})

# check_estimated ---------------------------------------------------------

test_that("check_estimated accepts a fitted estimator", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  m <- estimate(MEstimator(stacked_equations = psi, init = c(0)))
  expect_no_error(check_estimated(m))
})

test_that("check_estimated rejects an estimator that has not been fitted", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(check_estimated(m), "before calling")
})

# A fit that ran and came back without a variance is a different state from one
# that was never fitted, and the two are told apart by the estimates: a solve
# records them whatever became of the sandwich. Sending a caller who fitted the
# model to fit it again names the wrong problem and the remedy it implies does
# not exist.
na_bread_fit <- function() {
  y <- c(-1, 0, 1)
  # Finite at the root and not finite anywhere else, so the meat is built and
  # the finite differences behind the bread are not.
  psi <- function(theta) {
    if (theta[[1]] == 0) {
      matrix(y - theta[[1]], nrow = 1)
    } else {
      matrix(NA_real_, nrow = 1, ncol = length(y))
    }
  }
  at_the_root <- function(stacked_equations, init) 0
  suppressWarnings(estimate(
    MEstimator(stacked_equations = psi, init = 0),
    solver = at_the_root
  ))
}

test_that("check_estimated rejects a fit whose variance could not be built", {
  m <- na_bread_fit()
  expect_null(m@variance)
  expect_equal(unname(m@theta), 0)

  err <- expect_error(check_estimated(m))

  reported <- flatten_message(err)
  expect_no_match(reported, "before calling", fixed = TRUE)
  expect_match(reported, "variance", fixed = TRUE)
  expect_match(reported, "bread", fixed = TRUE)
})

test_that("the accessors report the missing variance rather than a missing fit", {
  m <- na_bread_fit()
  for (accessor in list(vcov, confint, summary)) {
    err <- expect_error(accessor(m))
    expect_no_match(flatten_message(err), "before calling", fixed = TRUE)
  }
})

# What the message above sends the reader to has to be there when they arrive.
# The estimates and the counts a solve records are properties of the solve
# rather than of the sandwich, so every accessor that reads nothing else answers
# on a fit whose variance could not be built, and only the ones that read the
# variance refuse.

test_that("the estimates and the counts survive a variance that could not be built", {
  m <- na_bread_fit()
  expect_equal(unname(coef(m)), 0)
  expect_identical(nobs(m), 3L)
  expect_identical(df.residual(m), 2L)
  expect_identical(generics::glance(m)$nobs, 3L)
})

# The same state on a fit made through the formula interface, which is what the
# accessors that route through `predict()` need. `na_bread_fit()` is what shows
# the state is reachable; an estimating equation the formula interface builds is
# finite wherever its data are, so the state is reached here by dropping the
# variance the solve did build.
no_variance_formula_fit <- function() {
  d <- data.frame(y = c(1, 3, 2, 5, 4), x = c(1, 2, 3, 4, 5))
  m <- m_estimate(y ~ x, data = d, .ee = ee_regression, model = "linear")
  m@variance <- NULL
  m@asymptotic_variance <- NULL
  m
}

test_that("a point prediction is made from a fit whose variance could not be built", {
  m <- no_variance_formula_fit()
  fitted_values <- unname(predict(m, type = "response"))

  expect_length(fitted_values, 5L)
  expect_equal(unname(fitted(m)), fitted_values)
  expect_equal(unname(residuals(m)), c(1, 3, 2, 5, 4) - fitted_values)
})

test_that("a standard error asked of the same fit names the variance", {
  m <- no_variance_formula_fit()
  for (asked in list(
    function(x) predict(x, se.fit = TRUE),
    function(x) predict(x, interval = "confidence")
  )) {
    err <- expect_error(asked(m))
    reported <- flatten_message(err)
    expect_no_match(reported, "before calling", fixed = TRUE)
    expect_match(reported, "variance", fixed = TRUE)
  }
})

# check_estimator_init ----------------------------------------------------

# Each of the four type-naming tests in this file asserts that the type
# description is present and that the offending object's own values are not.
# Only the second half separates a working message from one built with
# `{.obj_type_of}`, which cli does not recognize as an inline class: rather than
# refusing the name it falls back to pasting the values in, so a message meant
# to read "not a character matrix" reads out the whole object instead.
test_that("check_estimator_init names the type of a non-numeric init", {
  values <- c("alpha", "beta", "gamma", "delta")
  err <- expect_error(check_estimator_init(matrix(values, nrow = 2)))
  flat <- flatten_message(err)
  expect_match(flat, "not a character matrix", fixed = TRUE)
  expect_false(grepl("alpha", flat, fixed = TRUE))
})

# check_estimator_subset --------------------------------------------------

test_that("check_estimator_subset accepts NULL and valid indices", {
  expect_no_error(check_estimator_subset(NULL, 3))
  expect_no_error(check_estimator_subset(c(1L, 3L), 3))
  expect_no_error(check_estimator_subset(2L, 3))
})

test_that("check_estimator_subset rejects out-of-range and non-whole indices", {
  expect_error(check_estimator_subset(5L, 3), "between 1 and 3")
  expect_error(check_estimator_subset(0L, 3), "between 1 and 3")
  expect_error(check_estimator_subset(1.5, 3), "between 1 and 3")
})

test_that("check_estimator_subset rejects duplicated indices and names them", {
  expect_error(check_estimator_subset(c(1L, 1L), 3), "duplicated")
  expect_error(check_estimator_subset(c(2L, 3L, 2L), 3), "duplicated")
  # The offending index is reported.
  expect_error(check_estimator_subset(c(2L, 2L), 3), "2")
})

# check_psi_at_init -------------------------------------------------------

test_that("check_psi_at_init accepts a well-formed return", {
  expect_no_error(check_psi_at_init(matrix(1:5, nrow = 1), init = c(0)))
  expect_no_error(check_psi_at_init(
    rbind(1:3, 4:6),
    init = c(0, 0)
  ))
  # A plain vector is the single-parameter form.
  expect_no_error(check_psi_at_init(1:5, init = c(0)))
})

test_that("check_psi_at_init rejects a NULL return", {
  expect_error(check_psi_at_init(NULL, init = c(0)), "NULL")
})

test_that("check_psi_at_init reports a NULL return as NULL, not as a row count", {
  # `NULL` has no dimensions and so counts as one estimating equation, which
  # reads as a row-count mismatch against any `init` longer than one. The test
  # above supplies a length-one `init`, where the counts agree and the shape
  # branch stays quiet, so this is what pins the `NULL` check above it.
  err <- expect_error(check_psi_at_init(NULL, init = c(0, 0)))
  expect_match(conditionMessage(err), "\"NULL\"", fixed = TRUE)
  expect_false(inherits(err, "deli_psi_shape_error"))
})

test_that("check_psi_at_init rejects a non-numeric return", {
  expect_error(
    check_psi_at_init(matrix("a", nrow = 1), init = c(0)),
    "numeric"
  )
})

# The offending object here is an estimating-function return, so it holds one
# value per observation. Pasting it into the message is the worst case of the
# `{.obj_type_of}` fallback in the package, and the length of the message grows
# with the data.
test_that("check_psi_at_init names the type of a non-numeric return", {
  values <- c("alpha", "beta", "gamma")
  err <- expect_error(
    check_psi_at_init(matrix(values, nrow = 1), init = c(0, 0, 0))
  )
  flat <- flatten_message(err)
  expect_match(flat, "not a character matrix", fixed = TRUE)
  expect_false(grepl("alpha", flat, fixed = TRUE))
})

test_that("check_psi_at_init rejects non-finite values at init", {
  expect_error(
    check_psi_at_init(matrix(c(1, NaN, 3), nrow = 1), init = c(0)),
    "non-finite"
  )
  expect_error(
    check_psi_at_init(matrix(c(1, NA, 3), nrow = 1), init = c(0)),
    "non-finite"
  )
  expect_error(
    check_psi_at_init(matrix(c(1, Inf, 3), nrow = 1), init = c(0)),
    "non-finite"
  )
})

test_that("check_psi_at_init rejects a dimension mismatch for M-estimation", {
  expect_error(
    check_psi_at_init(matrix(1:5, nrow = 1), init = c(0, 0)),
    "1 estimating equation.*2 parameter"
  )
})

test_that("check_psi_at_init classes a dimension mismatch as a shape error", {
  # The formula interface catches this class alone, so that a NULL,
  # non-numeric, or non-finite return keeps its own more accurate message.
  expect_error(
    check_psi_at_init(matrix(1:5, nrow = 1), init = c(0, 0)),
    class = "deli_psi_shape_error"
  )
})

test_that("check_psi_at_init allows more equations than parameters when over-identified", {
  # GMM allows more estimating equations than parameters.
  expect_no_error(
    check_psi_at_init(
      rbind(1:5, 6:10),
      init = c(0),
      allow_over_identification = TRUE
    )
  )
})

test_that("check_psi_at_init rejects fewer equations than parameters when over-identified", {
  expect_error(
    check_psi_at_init(
      matrix(1:5, nrow = 1),
      init = c(0, 0),
      allow_over_identification = TRUE
    ),
    "less than or equal"
  )
  expect_error(
    check_psi_at_init(
      matrix(1:5, nrow = 1),
      init = c(0, 0),
      allow_over_identification = TRUE
    ),
    class = "deli_psi_shape_error"
  )
})

# The counts on their own say what the shape is and not what is wrong with it. A
# GMM system with fewer moment conditions than parameters is under-identified,
# which is a property of the system rather than of the starting values, so the
# fix is more moment conditions or fewer parameters and never a different `init`
# of the same length. The formula path reframes this failure against its
# automatic `init` and explains it there; the function path reports it as it
# stands and had nothing to say beyond the two counts.

test_that("the GMM shortfall explains that the system is under-identified", {
  err <- expect_error(
    check_psi_at_init(
      matrix(1:5, nrow = 1),
      init = c(0, 0),
      allow_over_identification = TRUE
    ),
    class = "deli_psi_shape_error"
  )
  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_match(flat, "under-identified", ignore.case = TRUE)
  expect_match(flat, "moment condition", ignore.case = TRUE)
})

test_that("the GMM shortfall does not advise a longer init", {
  # A longer `init` widens the shortfall, so it is the one direction the hint
  # must not point.
  err <- expect_error(
    check_psi_at_init(
      matrix(1:5, nrow = 1),
      init = c(0, 0, 0),
      allow_over_identification = TRUE
    ),
    class = "deli_psi_shape_error"
  )
  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_false(grepl("one longer", flat, fixed = TRUE))
  expect_false(grepl("longer", flat, fixed = TRUE))
})

test_that("an under-identified direct GMM fit is explained, not just counted", {
  # The path a user reaches through the constructor and estimate(), which
  # reported the two counts and stopped.
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) matrix(y - theta[1] - theta[2], nrow = 1)
  err <- expect_error(
    estimate(GMMEstimator(stacked_equations = psi, init = c(0, 0))),
    class = "deli_psi_shape_error"
  )
  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_match(flat, "under-identified", ignore.case = TRUE)
  expect_match(flat, "moment condition", ignore.case = TRUE)
  # The counts stay, because they are what the reader checks the diagnosis
  # against.
  expect_match(flat, "1", fixed = TRUE)
  expect_match(flat, "2", fixed = TRUE)
})

test_that("gmm_estimate()'s function interface explains an under-identified system", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) matrix(y - theta[1] - theta[2], nrow = 1)
  err <- expect_error(
    gmm_estimate(stacked_equations = psi, init = c(0, 0)),
    class = "deli_psi_shape_error"
  )
  expect_match(
    gsub("\\s+", " ", conditionMessage(err)),
    "under-identified",
    ignore.case = TRUE
  )
})

test_that("an M-estimation shape mismatch is not called under-identified", {
  # One equation per parameter is a requirement of M-estimation rather than an
  # identification condition, so this branch keeps its own wording.
  err <- expect_error(
    check_psi_at_init(matrix(1:5, nrow = 1), init = c(0, 0)),
    class = "deli_psi_shape_error"
  )
  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_false(grepl("under-identified", flat, fixed = TRUE))
  expect_match(flat, "one estimating equation per parameter", fixed = TRUE)
})

test_that("a GMM system with more equations than parameters is still accepted", {
  expect_no_error(
    check_psi_at_init(
      rbind(1:5, 6:10, 11:15),
      init = c(0, 0),
      allow_over_identification = TRUE
    )
  )
  expect_no_error(
    check_psi_at_init(
      rbind(1:5, 6:10),
      init = c(0, 0),
      allow_over_identification = TRUE
    )
  )
})

# The two shape branches are guarded separately and are separately movable, so
# each is pinned against the neighbors it sits above. The M-estimation one is
# pinned below the `NULL` check by a test further up; the finite check is what
# remains, and the GMM one is pinned against all three of its neighbors in turn.

test_that("check_psi_at_init reports an M-estimation mismatch ahead of a non-finite value", {
  # An estimating function reading a parameter beyond a short `init` returns
  # both the wrong number of rows and `NA`s, so the two failures arrive
  # together and the row count is the accurate diagnosis of the pair. This pins
  # the M-estimation branch above the finite check.
  err <- expect_error(
    check_psi_at_init(matrix(c(1, NaN, 3), nrow = 1), init = c(0, 0)),
    class = "deli_psi_shape_error"
  )
  expect_false(grepl("non-finite", conditionMessage(err), fixed = TRUE))
})

test_that("check_psi_at_init reports a GMM shortfall ahead of a non-finite value", {
  # An under-identified system has no solution at any starting values, so the
  # shortfall is the durable diagnosis and a non-finite value at the ones
  # supplied is not. This pins the GMM branch above the finite check.
  err <- expect_error(
    check_psi_at_init(
      matrix(c(1, NaN, 3), nrow = 1),
      init = c(0, 0),
      allow_over_identification = TRUE
    ),
    class = "deli_psi_shape_error"
  )
  expect_false(grepl("non-finite", conditionMessage(err), fixed = TRUE))
})

test_that("check_psi_at_init reports a NULL GMM return as NULL, not as a shortfall", {
  # `NULL` is dimensionless and so counts as one estimating equation, which
  # reads as a shortfall against any longer `init`. This pins the GMM branch
  # below the `NULL` check.
  err <- expect_error(
    check_psi_at_init(NULL, init = c(0, 0), allow_over_identification = TRUE)
  )
  expect_match(conditionMessage(err), "\"NULL\"", fixed = TRUE)
  expect_false(inherits(err, "deli_psi_shape_error"))
})

test_that("check_psi_at_init reports a non-numeric GMM return by its type, not as a shortfall", {
  # An `init` longer than the return cannot turn a numeric return into a
  # character one, so the type is never a symptom of the shortfall and counting
  # rows on a character matrix would misdiagnose it. This pins the GMM branch
  # below the numeric check.
  err <- expect_error(
    check_psi_at_init(
      matrix("a", nrow = 1),
      init = c(0, 0),
      allow_over_identification = TRUE
    ),
    "numeric"
  )
  expect_false(inherits(err, "deli_psi_shape_error"))
})

# eval_psi_at_init --------------------------------------------------------

test_that("eval_psi_at_init returns the estimating-function value", {
  psi <- function(theta) matrix(1:5 - theta[1], nrow = 1)
  expect_equal(eval_psi_at_init(psi, init = c(0)), matrix(1:5, nrow = 1))
})

test_that("eval_psi_at_init reframes only at the recorded automatic init", {
  psi <- function(theta) matrix(1:5 - theta[1], nrow = 1)
  attr(psi, "deli_auto_init") <- c(0, 0)
  expect_error(eval_psi_at_init(psi, init = c(0, 0)), "automatic zero")
  # Any other starting values were the caller's own choice, so the failure is
  # described in its own terms.
  err <- expect_error(eval_psi_at_init(psi, init = c(1, 1)))
  expect_false(grepl("automatic zero", conditionMessage(err), fixed = TRUE))
})

# The frame the psi-at-init aborts report ---------------------------------
#
# check_psi_at_init() judges a return the caller's own estimating function
# produced, from several frames below the method the caller reached, and every
# one of its aborts reported its own frame. That frame names the helper's
# parameters rather than anything in the call, and the helper appears in no man
# page, so nothing in the report is a call the caller can act on. Each abort
# names the frame that was reached instead, the way the delta method worker and
# the solver dispatcher already do.

reported_call <- function(err) {
  paste(deparse(conditionCall(err)), collapse = " ")
}

test_that("a NULL return reports the MEstimator method the caller reached", {
  m <- MEstimator(stacked_equations = function(theta) NULL, init = c(0))
  err <- expect_error(estimate(m))

  expect_match(conditionMessage(err), "NULL")
  expect_match(
    reported_call(err),
    "method(estimate, deli::MEstimator)",
    fixed = TRUE
  )
  expect_false(grepl("check_psi_at_init", reported_call(err), fixed = TRUE))
  expect_false(grepl("eval_psi_at_init", reported_call(err), fixed = TRUE))
})

test_that("a GMM shortfall reports the GMMEstimator method the caller reached", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  g <- GMMEstimator(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(0, 0)
  )
  err <- expect_error(estimate(g), class = "deli_psi_shape_error")

  expect_match(
    reported_call(err),
    "method(estimate, deli::GMMEstimator)",
    fixed = TRUE
  )
  expect_false(grepl("check_psi_at_init", reported_call(err), fixed = TRUE))
})

test_that("a formula fit reports the method its wrapper delegated to", {
  # `m_estimate()` builds the estimator and calls `estimate()` on it, so the
  # nearest frame a report raised below can name is the method the wrapper
  # reached rather than the wrapper itself. The failures a formula fit reframes
  # against an `init` it generated are different: those carry the entry point
  # the caller typed, and test-formula-interface.R pins them. This one is
  # raised at an `init` the caller supplied, so nothing reframes it.
  d <- data.frame(
    x = c(-1.1, 0.3, 1.7, -0.4, 0.9),
    y = c(0.2, 0.5, 0.8, 0.3, 0.7)
  )
  ee_null <- function(theta, X, y, ...) NULL
  err <- expect_error(
    m_estimate(y ~ x, data = d, .ee = ee_null, init = c(0, 0))
  )

  expect_match(conditionMessage(err), "NULL")
  expect_match(
    reported_call(err),
    "method(estimate, deli::MEstimator)",
    fixed = TRUE
  )
  expect_false(grepl("check_psi_at_init", reported_call(err), fixed = TRUE))
})

# error_looks_length_related ----------------------------------------------

# Most of what the predicate judges is the message, so the condition each
# message is carried in is the same throughout.
reads_as_length_problem <- function(message) {
  error_looks_length_related(simpleError(message))
}

test_that("error_looks_length_related recognizes the failures a short init causes", {
  expect_true(reads_as_length_problem("non-conformable arguments"))
  expect_true(reads_as_length_problem("non-conformable arrays"))
  expect_true(reads_as_length_problem("subscript out of bounds"))
  expect_true(reads_as_length_problem("length(theta) == 3 is not TRUE"))
  expect_true(reads_as_length_problem("`theta` must have length 3, not 2."))
})

test_that("error_looks_length_related reads an out-of-bounds subscript off the class", {
  # The class is what makes this family recognizable in a locale whose messages
  # are translated, so it is matched without reading the message at all.
  cnd <- structure(
    list(message = "not the English wording", call = NULL),
    class = c("subscriptOutOfBoundsError", "error", "condition")
  )
  expect_true(error_looks_length_related(cnd))
})

test_that("error_looks_length_related rejects failures that say nothing about a length", {
  expect_false(reads_as_length_problem("bad data"))
  expect_false(reads_as_length_problem("object 'z' not found"))
  # `%*%` on a non-numeric argument: a conformability error's neighbor, and not
  # a length problem.
  expect_false(
    reads_as_length_problem("requires numeric/complex matrix/vector arguments")
  )
})

# check_solver_return -----------------------------------------------------

test_that("check_solver_return accepts a numeric vector of the right length", {
  expect_no_error(check_solver_return(c(1, 2), 2))
})

test_that("check_solver_return rejects NULL, non-numeric, and wrong length", {
  expect_error(check_solver_return(NULL, 2), "numeric vector of length 2")
  expect_error(check_solver_return("a", 2), "numeric vector of length 2")
  expect_error(check_solver_return(c(1, 2, 3), 2), "numeric vector of length 2")
})

test_that("check_solver_return names the type of what the solver returned", {
  values <- c("alpha", "beta", "gamma")
  err <- expect_error(check_solver_return(values, 2))
  flat <- flatten_message(err)
  expect_match(flat, "It returned a character vector of length 3", fixed = TRUE)
  expect_false(grepl("alpha", flat, fixed = TRUE))
})

# check_deriv_method ------------------------------------------------------

test_that("check_deriv_method names the type of a non-string method", {
  err <- expect_error(check_deriv_method(c("capprox", "exact")))
  flat <- flatten_message(err)
  expect_match(flat, "not a character vector", fixed = TRUE)
  expect_false(grepl("capprox", flat, fixed = TRUE))
})

# check_penalty_shape -----------------------------------------------------

test_that("check_penalty_shape accepts scalar penalty and center", {
  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(theta, penalty = 0.5, center = 0))
})

test_that("check_penalty_shape accepts penalty and center matching theta length", {
  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(
    theta,
    penalty = c(0.1, 0.2, 0.3),
    center = c(0, 0, 0)
  ))
})

test_that("check_penalty_shape accepts mixed scalar and vector arguments", {
  theta <- c(1, 2, 3)
  # scalar penalty, vector center

  expect_no_error(check_penalty_shape(
    theta,
    penalty = 0.5,
    center = c(0, 0, 0)
  ))
  # vector penalty, scalar center
  expect_no_error(check_penalty_shape(
    theta,
    penalty = c(0.1, 0.2, 0.3),
    center = 0
  ))
})

test_that("check_penalty_shape rejects penalty with wrong length", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, 0.2), center = 0),
    "penalty.*single number.*same length"
  )
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, 0.2, 0.3, 0.4), center = 0),
    "penalty.*single number.*same length"
  )
})

test_that("check_penalty_shape rejects center with wrong length", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = 0.5, center = c(0, 0)),
    "center.*single number.*same length"
  )
  expect_error(
    check_penalty_shape(theta, penalty = 0.5, center = c(0, 0, 0, 0)),
    "center.*single number.*same length"
  )
})

test_that("check_penalty_shape rejects negative penalty values", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = -1, center = 0),
    "penalty.*non-negative"
  )
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, -0.2, 0.3), center = 0),
    "penalty.*non-negative"
  )
})

test_that("check_penalty_shape accepts zero penalty values", {
  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(theta, penalty = 0, center = 0))
  expect_no_error(check_penalty_shape(theta, penalty = c(0, 0, 0), center = 0))
})

# check_survival_data_valid -----------------------------------------------

test_that("check_survival_data_valid accepts valid survival data", {
  expect_no_error(check_survival_data_valid(
    delta = c(0, 1, 1, 0),
    time = c(1, 2, 3, 4)
  ))
  expect_no_error(check_survival_data_valid(
    delta = c(1, 1, 1),
    time = c(0.5, 1.5, 10)
  ))
  expect_no_error(check_survival_data_valid(
    delta = c(0, 0, 0),
    time = c(1, 2, 3)
  ))
})

test_that("check_survival_data_valid accepts data with NAs", {
  # NAs in delta should be ignored

  expect_no_error(check_survival_data_valid(
    delta = c(0, 1, NA, 0),
    time = c(1, 2, 3, 4)
  ))
  # NAs in time should be ignored
  expect_no_error(check_survival_data_valid(
    delta = c(0, 1, 1, 0),
    time = c(1, NA, 3, 4)
  ))
  # NAs in both

  expect_no_error(check_survival_data_valid(
    delta = c(0, NA, 1, 0),
    time = c(1, NA, 3, 4)
  ))
})

test_that("check_survival_data_valid rejects delta values other than 0 or 1", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 2), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
  expect_error(
    check_survival_data_valid(delta = c(-1, 0, 1), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
  expect_error(
    check_survival_data_valid(delta = c(0.5, 1, 0), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
})

test_that("check_survival_data_valid rejects delta with invalid values even when NAs present", {
  expect_error(
    check_survival_data_valid(delta = c(0, NA, 2), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
})

test_that("check_survival_data_valid rejects non-positive times", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(0, 2, 3)),
    "times.*positive"
  )
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(-1, 2, 3)),
    "times.*positive"
  )
})

test_that("check_survival_data_valid rejects non-positive times even when NAs present", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(NA, -1, 3)),
    "times.*positive"
  )
})

# check_dx ----------------------------------------------------------------

# Every value the finite-difference step may not take, in one place, so that the
# helper and each entry point taking `dx` are held to the same set. None of them
# is caught anywhere at present: a zero or negative step is replaced by the step
# floor of approx_differentiation() at any non-zero parameter and leaves a
# division by zero at a parameter of exactly zero, a longer vector is recycled
# one element per parameter, and a missing or infinite step reaches the quotient
# and returns an NA or NaN derivative.
invalid_dx <- list(
  zero = 0,
  negative = -1e-9,
  vector = c(1e-9, 1e-8),
  character = "1e-9",
  missing = NA_real_,
  infinite = Inf
)

# Each entry point is held to the whole set, and every rejection has to name the
# argument rather than surface a base arithmetic error from inside the
# difference quotient.
expect_dx_rejected <- function(call_with) {
  for (case in names(invalid_dx)) {
    expect_error(
      call_with(invalid_dx[[case]]),
      regexp = "dx",
      class = "rlang_error",
      info = case
    )
  }
  invisible(NULL)
}

test_that("check_dx accepts a single positive finite number", {
  expect_no_error(check_dx(1e-9))
  expect_no_error(check_dx(1))
  # `ee_percentile()` is differentiated with a step of 1 to smooth an indicator,
  # so a large step is as valid as a small one.
  expect_no_error(check_dx(1e-3))
})

test_that("check_dx rejects every step that is not one positive finite number", {
  expect_dx_rejected(check_dx)
})

test_that("check_dx says what a step has to be", {
  # Three separate reasons to reject, each of which has to reach the message,
  # however they are worded together.
  expect_error(check_dx(0), "positive")
  expect_error(check_dx(-1e-9), "positive")
  expect_error(check_dx(c(1e-9, 1e-8)), "single")
  expect_error(check_dx(numeric(0)), "single")
  expect_error(check_dx(NA_real_), "finite")
  expect_error(check_dx(Inf), "finite")
})

# `dx` at the user entry points -------------------------------------------

test_that("estimate() rejects an invalid dx", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  expect_dx_rejected(function(dx) {
    estimate(MEstimator(stacked_equations = psi, init = 0), dx = dx)
  })
  expect_dx_rejected(function(dx) {
    estimate(GMMEstimator(stacked_equations = psi, init = 0), dx = dx)
  })
})

test_that("m_estimate() rejects an invalid dx on both interfaces", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  expect_dx_rejected(function(dx) {
    m_estimate(stacked_equations = psi, init = 0, dx = dx)
  })
  expect_dx_rejected(function(dx) {
    m_estimate(
      mpg ~ wt,
      data = mtcars,
      .ee = ee_regression,
      model = "linear",
      dx = dx
    )
  })
})

test_that("gmm_estimate() rejects an invalid dx on both interfaces", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  expect_dx_rejected(function(dx) {
    gmm_estimate(stacked_equations = psi, init = 0, dx = dx)
  })
  expect_dx_rejected(function(dx) {
    gmm_estimate(
      mpg ~ wt,
      data = mtcars,
      .ee = ee_regression,
      model = "linear",
      dx = dx
    )
  })
})

test_that("compute_sandwich() rejects an invalid dx", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  expect_dx_rejected(function(dx) {
    compute_sandwich(psi, theta = mean(y), dx = dx)
  })
})

test_that("delta_method() rejects an invalid dx on both methods", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  fit <- m_estimate(stacked_equations = psi, init = 0)
  expect_dx_rejected(function(dx) {
    delta_method(fit, transform = function(theta) exp(theta[1]), dx = dx)
  })
  expect_dx_rejected(function(dx) {
    delta_method(
      coef(fit),
      transform = function(theta) exp(theta[1]),
      covariance = vcov(fit),
      dx = dx
    )
  })
})

test_that("the survival prediction helpers reject an invalid dx", {
  expect_dx_rejected(function(dx) {
    survival_predictions(
      times = c(1, 2),
      theta = c(0.1, 0.2),
      covariance = diag(2),
      distribution = "weibull",
      dx = dx
    )
  })
  expect_dx_rejected(function(dx) {
    aft_predictions_function(
      X = matrix(c(1, 1), nrow = 1),
      times = c(5, 10),
      theta = c(1, 0.5, 0),
      covariance = diag(3),
      distribution = "weibull",
      dx = dx
    )
  })
})

test_that("predict() rejects an invalid dx", {
  set.seed(414)
  n <- 200
  d <- data.frame(x = round(stats::rnorm(n), 3))
  event_time <- exp(1 + 0.5 * d$x + 0.7 * stats::rlogis(n))
  censor_time <- stats::rexp(n, rate = 0.03)
  d$time <- round(pmin(event_time, censor_time), 4)
  d$status <- as.numeric(event_time <= censor_time)
  fit <- m_estimate(
    time ~ x,
    data = d,
    .ee = ee_aft,
    distribution = "weibull",
    event = status,
    init = c(mean(log(d$time)), 0, 0)
  )

  expect_dx_rejected(function(dx) {
    predict(fit, times = c(5, 10), measure = "survival", se.fit = TRUE, dx = dx)
  })
  # The step is validated where it is supplied rather than where it is used, so
  # a prediction asking for no standard error reports it too rather than
  # accepting a step it would never have taken.
  expect_dx_rejected(function(dx) {
    predict(fit, times = c(5, 10), measure = "survival", dx = dx)
  })
})
