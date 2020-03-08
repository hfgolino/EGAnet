#' Plot for net.loads objects of \code{\link[EGAnet]{net.loads}} results
#'
#' Returns a plot of the net.loads results of \code{\link[EGAnet]{net.loads}} results
#'
#' @param object An \code{\link[EGAnet]{net.loads}} object
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(wmt)
#'
#' }
#'
#' # Network loadings
#' nloads <- net.loads(ega.wmt)
#' 
#' # Plot
#' plot(nloads)
#' 
#' @export
#'
## S3 method for class 'NetLoads'
#'
#Plot function for NetLoads
#Updated 05.03.2020
plot.NetLoads <- function(object, ...) {
  
  plot(object$plot)
}
