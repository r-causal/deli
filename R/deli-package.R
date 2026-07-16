#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import S7
#' @importFrom stats dnorm plogis pnorm pt qlogis qnorm quantile qt
#' @importFrom utils capture.output
## usethis namespace: end
NULL

.onLoad <- function(...) {
  S7::methods_register()
}
