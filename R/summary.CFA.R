#' Summary for CFA objects
#'
#' Returns a summary of the CFA results
#'
#' @param object An \code{\link{CFA}} object
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#'\donttest{
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#'
#' #esimtate CFA
#' cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = TRUE, data = wmt2)
#'
#' #summary of CFA results
#' summary(cfa.wmt)
#'}
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#' @export
#'
## S3 method for class 'CFA'
#'
#Summary function for CFA
summary.CFA <- function(object, ...) {
  cat("Summary: Confirmatory Factor Analysis:\n")
  print(object$summary)
  cat("\n FIt Measures:\n")
  print(object$fit.measures)
}
