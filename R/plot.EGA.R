#'  Plot method for EGA objects.
#'
#' \code{plot} Plots the EGA result using \code{\link{qgraph}}
#'
#' @param object An EGA object
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' summary(ega.wmt)
#' plot(ega.wmt)
#'
#' \dontrun{
#' plot(EGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

## S3 method for class 'EGA'

plot.EGA <- function(object, layout = "spring", vsize = 6, ...) {
  groups = as.factor(object$wc)
  variable = object$glasso
  qgraph(variable, layout = layout, vsize = vsize, groups = groups, ...)
}
