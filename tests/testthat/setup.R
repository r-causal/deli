# Session options do not propagate to parallel test workers, so the
# warning promotion the suite is verified under is set here, where every
# worker sources it. CRAN runs keep the default warning behavior.
if (Sys.getenv("NOT_CRAN") == "true") {
  options(warn = 2)
}
