#' Plot method for \code{\link[EGAnet]{EGA}}objects
#'
#' Plots the \code{\link[EGAnet]{EGA}} result using \code{\link[qgraph]{qgraph}}
#'
#' @param x An \code{\link[EGAnet]{EGA}} object
#'
#' @param title Character. Title of the plot
#'
#' @param vsize An integer indicating the size of the nodes.
#' Default vsize = 6
#'
#' @param ... Arguments passed to \code{\link[qgraph]{qgraph}}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' 
#' }
#'
#' # Summary of EGA results
#' summary(ega.wmt)
#'
#' # Plot EGA network
#' plot(ega.wmt, vsize = 6, label.prop = 1)
#' 
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @export
#'
## S3 method for class 'EGA'
#'
#Plot EGA function
plot.EGA <- function(x, title = "", vsize = 6,  ...) {
    plot.ega <- qgraph::qgraph(x$network, layout = "spring", vsize = vsize, groups = as.factor(x$wc), ...)

}
