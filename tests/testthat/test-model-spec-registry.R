# The registry of observation-indexed argument names, pinned at both ends.
#
# `per_observation_ee_args()` is a written-out list of names, and nothing in the
# package ties it to the arguments the estimating equations actually take. It is
# written out for the reason `appended_param_name()` is: the property being
# asked about is not visible in what an equation returns, and two arguments of
# the same length can be a vector of censoring indicators and a length-2 vector
# of truncation bounds. So the list cannot be derived, and a new `ee_*()`
# argument would join the package with nothing to say which half of a model
# specification it belongs in.
#
# What replaces the derivation is the same both-ends pin that test-param-names.R
# and test-param-names-builtin.R give the appended-parameter table: one end is a
# second written-out table, here rather than in `R/`, classifying every argument
# of every exported `ee_*()` function; the other end is the registry and the two
# functions that read it. A new argument fails the first test in this file until
# it is classified, and a classification the registry disagrees with fails the
# second, so neither table can be updated without the other.
#
# The three kinds divide by what the argument is and by whether the split ever
# sees it.
#
# `theta` and `X` are filled by the formula interface itself.
# `prepare_formula_psi()` passes both by name, so a caller who names either in
# `...` makes `do.call()` match one formal twice and the fit fails there. Neither
# can reach the argument list the split is handed, and the registry says nothing
# about them. `X` is the one observation-indexed design argument the registry
# omits, where `W`, `V`, `X1`, and `X0` are all on it; the last test in this file
# is why that omission cannot be reached.
#
# `y` is on the registry although the interface fills it too, and it has to be:
# the response is passed positionally into the first argument that is neither
# `theta` nor `X` nor named in `...`, so a caller who does name `y` displaces the
# response to the next free argument and their `y` is forwarded like any other.
#
# Everything else is either one value per row of the sample the fit was made on,
# which is what the registry names, or a description of the model, which is what
# is left once the registry's names are removed.

# ---- Helpers ----

exported_ee_functions <- function() {
  sort(grep("^ee_", getNamespaceExports("deli"), value = TRUE))
}

ee_formal_names <- function(name) {
  names(formals(get(name, envir = asNamespace("deli"))))
}

# What each argument of each exported estimating equation is. Every argument of
# every one of them appears exactly once, and the tests below fail if the set
# drifts from the signatures in either direction.
ee_formal_kinds <- function() {
  c(
    # Filled by the formula interface and never forwarded through `...`.
    theta = "interface",
    X = "interface",
    # One value per row of the fitted sample. Shared across the regression,
    # GLM, and survival families.
    weights = "per_observation",
    offset = "per_observation",
    event = "per_observation",
    time = "per_observation",
    y = "per_observation",
    # Causal designs: treatment, instrument, and the nuisance and
    # counterfactual design matrices.
    A = "per_observation",
    Z = "per_observation",
    W = "per_observation",
    V = "per_observation",
    X1 = "per_observation",
    X0 = "per_observation",
    # Measurement error: mismeasured values, gold-standard values, and the
    # sample indicator.
    y_star = "per_observation",
    a = "per_observation",
    a_star = "per_observation",
    r = "per_observation",
    # Missingness and dose-response data.
    delta = "per_observation",
    q_eval = "per_observation",
    dose = "per_observation",
    response = "per_observation",
    # The outcome distribution and how the mean is linked to the design, and
    # the fixed variance power the tweedie reads from `hyperparameter`.
    distribution = "specification",
    link = "specification",
    model = "specification",
    hyperparameter = "specification",
    # The penalized-regression settings. `penalty` may carry one value per
    # coefficient, which is a value per design column rather than per row, so
    # it describes the model and applies unchanged to new data.
    penalty = "specification",
    center = "specification",
    epsilon = "specification",
    ratio = "specification",
    gamma = "specification",
    s = "specification",
    # The robust and dose-response loss functions and their tuning constant.
    loss = "specification",
    k = "specification",
    # Bounds and scalars naming a point on a fitted curve or a censoring
    # limit, none of which is indexed by observation. `truncate` is the
    # length-2 vector the registry's own documentation contrasts with a
    # censoring indicator.
    lower = "specification",
    upper = "specification",
    ed50 = "specification",
    steepness = "specification",
    truncate = "specification",
    q = "specification",
    # The pooled logistic time design and the grid it is defined on, both
    # indexed by time step rather than by observation, which is why
    # `plogit_predict()` reads them from the specification half.
    S = "specification",
    unique_times = "specification",
    # The spline specifications, one per design column.
    specifications = "specification",
    # Switches and a function-valued argument.
    force_continuous = "specification",
    log_theta = "specification",
    H_function = "specification",
    # The external coefficient a regression-calibration fit corrects, a scalar
    # rather than a column of data.
    beta = "specification"
  )
}

# A one-row-per-argument list to hand the split, so the two functions are asked
# about every name at once. The values are placeholders: the split reads names
# and nothing else.
all_ee_args <- function() {
  kinds <- ee_formal_kinds()
  stats::setNames(as.list(seq_along(kinds)), names(kinds))
}

# The data behind the two tests that go through the formula interface.
registry_data <- function() {
  data.frame(
    count = c(3, 4, 5, 6, 7, 8, 3, 4),
    x = c(1, 2, 3, 4, 5, 6, 7, 8),
    exposure = c(2, 3, 4, 5, 6, 7, 8, 9)
  )
}

# The classification table against the signatures ----------------------------

test_that("the table classifies every argument of every exported equation", {
  equations <- exported_ee_functions()
  # Fails if the enumeration stops finding the equations at all, which would
  # leave the comparison below trivially true.
  expect_gt(length(equations), 30L)

  formal_names <- unique(unlist(lapply(equations, ee_formal_names)))
  expect_setequal(names(ee_formal_kinds()), formal_names)
})

test_that("the table classifies each argument exactly once and by one kind", {
  kinds <- ee_formal_kinds()

  expect_identical(anyDuplicated(names(kinds)), 0L)
  expect_setequal(
    unique(unname(kinds)),
    c("interface", "per_observation", "specification")
  )
})

# The registry against the classification table ------------------------------

test_that("the registry names exactly the observation-indexed arguments", {
  kinds <- ee_formal_kinds()
  expect_setequal(
    per_observation_ee_args(),
    names(kinds)[kinds == "per_observation"]
  )
})

test_that("the registry names neither the parameters nor the design", {
  # The two arguments a forwarded one can never displace. `y`, which the
  # interface also fills, is on the registry for the reason the header gives and
  # the test at the end of this file demonstrates.
  expect_false("theta" %in% per_observation_ee_args())
  expect_false("X" %in% per_observation_ee_args())
})

test_that("the split routes every argument of every equation the table's way", {
  kinds <- ee_formal_kinds()
  args <- all_ee_args()

  expect_setequal(
    names(spec_obs_args(args)),
    names(kinds)[kinds == "per_observation"]
  )
  expect_setequal(
    names(spec_ee_args(args)),
    names(kinds)[kinds != "per_observation"]
  )
  # The two halves partition the arguments, so no name is recorded twice and
  # none is dropped.
  expect_setequal(
    c(names(spec_ee_args(args)), names(spec_obs_args(args))),
    names(kinds)
  )
})

# Why the classification is reachable at all ---------------------------------

test_that("a forwarded response is recorded as per-observation data", {
  # The interface passes the response positionally into the first argument that
  # is neither `theta` nor `X` nor named in `...`, so naming `y` moves the
  # response along to `weights` and the caller's `y` is forwarded like any
  # other argument. That is the one way an argument the interface also fills
  # reaches the split, and it is why `y` is on the registry.
  data <- registry_data()
  m <- m_estimate(
    count ~ x,
    data = data,
    .ee = ee_regression,
    model = "linear",
    y = data$exposure
  )
  spec <- m@model_spec

  expect_identical(spec$ee_obs_args, list(y = data$exposure))
  expect_identical(spec$ee_spec_args, list(model = "linear"))
  # The response the formula named is recorded on its own, where it always is.
  expect_identical(unname(spec$y), data$count)
})

test_that("the design and the parameters cannot reach the split at all", {
  # `theta` and `X` are passed by name, so a dot of either name makes the call
  # match one formal twice and the fit fails before anything is recorded. That
  # is what makes the registry's silence about them safe, and in particular
  # what makes it safe that `X` is the one observation-indexed design argument
  # it does not name.
  data <- registry_data()
  fit <- function(...) {
    m_estimate(
      count ~ x,
      data = data,
      .ee = ee_regression,
      model = "linear",
      init = c(0, 0),
      ...
    )
  }

  expect_error(fit(X = base::cbind(1, data$x)), "multiple actual arguments")
  expect_error(fit(theta = c(0, 0)), "multiple actual arguments")
})

test_that("the one argument two equations disagree about cannot collide", {
  # `delta` is a per-observation missingness indicator in
  # `ee_mean_sensitivity_analysis()` and a scalar effective-dose level in
  # `ee_emax_ed()` and `ee_loglogistic_ed()`. A registry keyed by name holds one
  # answer for both, and the per-observation one is what it holds. Nothing is
  # misclassified today because neither `_ed` equation can be driven by a
  # formula at all: neither takes an `X` argument for the model matrix, so
  # `check_formula_ee_signature()` refuses the fit before any argument is split.
  # Giving either of them a design would make this test the place that says so.
  expect_true("delta" %in% per_observation_ee_args())
  for (name in c("ee_emax_ed", "ee_loglogistic_ed")) {
    expect_false("X" %in% ee_formal_names(name), label = name)
  }
})
