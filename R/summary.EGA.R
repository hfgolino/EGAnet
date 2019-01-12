#'  Summary method for EGA objects.
#'
#' \code{summary} Returns a summary of the EGA results.
#'
#' @param object An EGA object
#' @param ... potentially further arguments; unused currently.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' summary(ega.wmt)
#'
#' \dontrun{
#' summary(EGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

## S3 method for class 'EGA'

summary.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}
