# Abstract parent class for the deli estimators.
#
# deli_estimator holds the properties and post-estimation behavior shared by
# MEstimator and GMMEstimator. The shared generics (coefficients, variance,
# confidence intervals, tidiers, influence functions, printing, and so on)
# dispatch on this class so that both estimators inherit one definition. It is
# abstract and never instantiated directly; the two concrete estimators supply
# their own constructors and any class-specific properties.
#
# The property definitions live here and are inherited by both children.
# MEstimator adds nothing beyond them; GMMEstimator adds the weight matrix and
# the over-identification controls.
#
# summed_equations and check_summed_equations are the pair that says how the
# estimating functions are reduced to one value per equation, and both children
# take them because both reduce them: the M path for its solver and its bread,
# the GMM path for its objective and its bread. They sit together because the
# second is about the first alone, which is also why neither is an argument of
# estimate(): a fit that supplies no reduction has nothing to check, and the
# settings estimate() takes describe one solve rather than the system.
#
# model_spec is the one property no constructor takes. A fit built by hand from
# a `stacked_equations` closure has no model to describe: the estimator is
# handed an arbitrary function and cannot recover a formula, a design, or a
# response from it. Only the formula interface knows those, so
# `m_estimate.formula()` and `gmm_estimate.formula()` set the property by
# assignment after construction and every other route leaves it `NULL`. Keeping
# it out of the constructors is what lets the public signatures stay as they
# are; see `formula_model_spec()` in R/model-spec.R for what the list holds.
deli_estimator <- new_class(
  "deli_estimator",
  abstract = TRUE,
  properties = list(
    stacked_equations = class_function,
    summed_equations = new_property(
      class = NULL | class_function,
      default = NULL
    ),
    check_summed_equations = new_property(
      class = class_logical,
      default = TRUE
    ),
    init = class_numeric,
    subset = new_property(
      class = NULL | class_integer,
      default = NULL
    ),
    finite_correction = new_property(
      class = NULL | class_character,
      default = NULL
    ),
    model_spec = new_property(
      class = NULL | class_list,
      default = NULL
    ),
    n_obs = new_property(
      class = NULL | class_integer,
      default = NULL
    ),
    n_params = new_property(class_integer),
    theta = new_property(
      class = NULL | class_numeric,
      default = NULL
    ),
    bread = new_property(
      class = NULL | class_double,
      default = NULL
    ),
    meat = new_property(
      class = NULL | class_double,
      default = NULL
    ),
    variance = new_property(
      class = NULL | class_double,
      default = NULL
    ),
    asymptotic_variance = new_property(
      class = NULL | class_double,
      default = NULL
    )
  )
)
