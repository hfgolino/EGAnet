#' Summary for EGA objects
#'
#' Returns a summary of the EGA results
#'
#' @param object An \code{\link{EGA}} object
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \donttest{
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#'
#' #summary of EGA reults
#' summary(ega.wmt)
#'}
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @export
#'
## S3 method for class 'EGA'
#'
#Summary function for EGA
summary.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}
