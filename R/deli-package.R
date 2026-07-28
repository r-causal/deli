#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import S7
#' @importFrom stats dnorm median plogis pnorm pt qlogis qnorm quantile qt
## usethis namespace: end
NULL

.onLoad <- function(...) {
  S7::methods_register()
}
