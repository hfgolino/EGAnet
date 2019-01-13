#' Plot method for EGA objects.
#'
#' Plots the EGA result using \code{\link{qgraph}}
#'
#' @param x An \code{\link[EGA]{EGA}} object
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
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' 
#' #summary of EGA results
#' summary(ega.wmt)
#'
#' #plot EGA network
#' plot(ega.wmt, vsize = 6)
#'
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @export
#'
## S3 method for class 'EGA'
#'
#Plot EGA function
plot.EGA <- function(x, title = "", vsize = 6,  ...) {
    plot.ega <- qgraph::qgraph(x$network, layout = "spring", vsize = vsize, groups = as.factor(x$wc), ...)

}
