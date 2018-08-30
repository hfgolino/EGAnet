#'  Plot method for bootEGA objects.
#'
#' \code{plot.bootEGA} Plots the bootEGA typical structure using \code{\link{qgraph}}
#'
#' @param bootega.obj A bootEGA object
#' @param vsize An integer indicating the size of the nodes. Default vsize = 6.
#' @param ... Arguments passed to 'qgraph'.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' \dontrun{
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 10, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "GGM",
#' type = "parametric", ncores = 4, confirm = ega.wmt$wc)
#' plot.bootEGA(boot.wmt)
#' }
#'
#' \dontrun{
#' plot(bootEGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

## S3 method for class 'bootEGA'

plot.bootEGA <- function(bootega.obj, vsize = 6,  ...){
  require(qgraph)
  qgraph(bootega.obj$typicalGraph$graph, layout = "spring",
         groups = as.factor(bootega.obj$typicalGraph$wc),
         vsize = vsize, ...)

}
