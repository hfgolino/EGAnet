#' Print method for EGA objects.
#'
#' Returns a summary of the EGA results
#'
#' @param x An \code{\link{EGA}} object
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#'\donttest{
#' #estimtae EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#'
#' #print EGA results
#' print(ega.wmt)
#'}
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @export
#'
## S3 method for class 'EGA'
#'
#Print EGA function
print.EGA <- function(x, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(x$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dim.variables)
}
