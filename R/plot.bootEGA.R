#' Plot method for \code{\link[EGAnet]{bootEGA}} objects
#'
#' Plots \code{\link[EGAnet]{bootEGA}} typical structure using \code{\link[qgraph]{qgraph}}
#'
#' @param x A \code{\link[EGAnet]{bootEGA}} object
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
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#'
#' # Estimate bootEGA
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 10, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "GGM",
#' type = "parametric", ncores = 4, confirm = ega.wmt$wc)
#' 
#' }
#'
#' # Plot bootEGA
#' plot(boot.wmt)
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @export
#'
## S3 method for class 'bootEGA'
#'
#Plot bootEGA function
plot.bootEGA <- function(x, vsize = 6,...){

  qgraph::qgraph(x$typicalGraph$graph, layout = "spring",
         groups = as.factor(x$typicalGraph$wc),
         vsize = vsize, ...)

}
