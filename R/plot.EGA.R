#'  Plot method for EGA objects.
#'
#' \code{plot} Plots the EGA result using \code{\link{qgraph}}
#'
#' @param ega.obj An EGA object
#' @param title Character. Title of the plot
#' @param vsize An integer indicating the size of the nodes. Default vsize = 6.
#' @param ... Arguments passed to 'qgraph'.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' summary(ega.wmt)
#' plot(ega.wmt, vsize = 6)
#'
#' \dontrun{
#' plot(EGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

## S3 method for class 'EGA'

plot.EGA <- function(ega.obj, title = "", vsize = 6,  ...){
    
    plot.ega <- qgraph::qgraph(ega.obj$network, layout = "spring", vsize = vsize, groups = as.factor(ega.obj$wc), ...)

}
