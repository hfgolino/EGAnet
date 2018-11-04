#'  Plot method for CFA objects.
#'
#' \code{plot} Plots the CFA structure using \code{\link{semPlot}}
#'
#' @param object An CFA object.
#' @param layout Layout of plot (see \code{\link[semPlot]{semPaths}}).
#' Defaults to "spring"
#' @param vsize Size of objects in plot.
#' Defaults to 6
#' @param ... Arguments passed to ’semPaths’ in semPlot.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = FALSE, data = wmt2)
#' plot(cfa.wmt)
#'
#' \dontrun{
#' plot(CFA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#' @export

## S3 method for class 'CFA'

#Plot CFA:
plot.CFA <- function(object, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(object$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
           sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
           "std", cut = 0.5, ...)
}
