#' Summary for \code{\link[EGAnet]{dynEGA}} objects (Group Effects)
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} results (Group Effects)
#'
#' @param object An \code{\link[EGAnet]{dynEGA}} object (Group Effects)
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
#' #Summary of dynEGA results
#' summary(dyn.group)
#'}
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA.Groups' (Group Effects)
#'
#Summary function for dynEGA (Group Effects)
summary.dynEGA.Groups <- function(object, ...) {
  for(i in 1:length(object$dynEGA)){
    cat("dynEGA Results (Group Effects):\n")
    cat("Group:", names(object$dynEGA[i]))
    cat("\nNumber of Dimensions:\n")
    print(object$dynEGA[[i]]$n.dim)
    cat("\nItems per Dimension:\n")
    print(object$dynEGA[[i]]$dim.variables)
  }
}

