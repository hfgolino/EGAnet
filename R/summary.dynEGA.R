#' Summary for \code{\link[EGAnet]{dynEGA}} objects
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} results
#'
#' @param object An \code{\link[EGAnet]{dynEGA}} object
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
#' #Summary of dynEGA reults
#' summary(dyn.random)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA'
#'
#Summary function for dynEGA
summary.dynEGA <- function(object, ...) {
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(object$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dynEGA$dim.variables)
}

