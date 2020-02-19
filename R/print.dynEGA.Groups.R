#' Print method for \code{\link[EGAnet]{dynEGA}} objects (Group Effects)
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} objects (Group Effects)
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} objects (Group Effects)
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Estimate dynEGA
#' dyn.group <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' effects = "group", model = "glasso")
#'
#' #Print dynEGA results
#' print(dyn.group)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @importFrom graphics par
#'
#' @export
#'
## S3 method for class 'dynEGA.Groups' (Group Effects)
#'
#Print dynEGA function
print.dynEGA.Groups <- function(x, ...) {
    for(i in 1:length(x$dynEGA)){
      cat("dynEGA Results (Group Effects):\n")
      cat("Group:", names(x$dynEGA[i]))
      cat("\nNumber of Dimensions:\n")
      print(x$dynEGA[[i]]$n.dim)
      cat("\nItems per Dimension:\n")
      print(x$dynEGA[[i]]$dim.variables)
    }
}
