#' Print method for \code{\link[EGAnet]{dynEGA}} objects
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} objects
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} objects
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Estimate dynEGA
#' dyn.random <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' level = "population", model = "glasso")
#'
#' #Print dynEGA results
#' print(dyn.random)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA' (Level: Population)
#'
#Print dynEGA function
print.dynEGA<- function(x, ...) {
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(x$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dynEGA$dim.variables)
}
