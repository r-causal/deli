# ---- Why the three broom generics are re-exported ----------------------------
# `tidy()`, `glance()`, and `augment()` are generics owned by the generics
# package, and deli's methods for them are registered against those objects.
# Re-exporting them means `library(deli)` alone puts the generics on the search
# path, so `tidy(fit)` works without attaching generics or broom first.
#
# What is exported here is the same object as `generics::tidy()`, not a generic
# of deli's own carrying the same name. That matters: a separate generic would
# answer the bare call while splitting dispatch across two objects, so a fit
# handed to a broom-aware package would find no method. Because the binding is
# the generics package's own, the method registrations in R/tidiers.R and
# R/augment.R dispatch through it unchanged.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment
