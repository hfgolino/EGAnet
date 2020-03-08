#' Summary for net.loads objects of \code{\link[EGAnet]{net.loads}} results
#'
#' Returns a summary of the net.loads results of \code{\link[EGAnet]{net.loads}} results
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
#' # Summary
#' summary(nloads)
#' 
#' @export
#'
## S3 method for class 'NetLoads'
#'
#Summary function for NetLoads
#Updated 05.03.2020
summary.NetLoads <- function(object, ...) {
  
  object$std[which(abs(object$std) <= object$MinLoad, arr.ind = TRUE)] <- ""
  
  print(object$std)
}
