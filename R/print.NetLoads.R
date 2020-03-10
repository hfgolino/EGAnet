#' Print for net.loads objects of \code{\link[EGAnet]{net.loads}} results
#'
#' Returns a summary of the net.loads results of \code{\link[EGAnet]{net.loads}} results
#'
#' @param x An \code{\link[EGAnet]{net.loads}} object
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
#' # Print
#' print(nloads)
#' 
#' @export
#'
## S3 method for class 'NetLoads'
#'
#Print function for NetLoads
#Updated 09.03.2020
print.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$MinLoad, arr.ind = TRUE)] <- ""
  
  print(x$std)
  message("Loadings <= ", x$MinLoad, " are blank")
}
