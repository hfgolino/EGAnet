#' Summary for \code{\link[EGAnet]{dynEGA}} objects (Level: Individual)
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} results (Level: Individual)
#'
#' @param object An \code{\link[EGAnet]{dynEGA}} object (Level: Individual)
#'
#' @param ... potentially further arguments (\strong{unused currently})
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Estimate dynEGA
#' dyn.individual <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' level = "individual", model = "glasso")
#'
#' #Summary of dynEGA results
#' summary(dyn.individual)
#'}
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA.Individuals' (Level: Individual)
#'
#Summary function for dynEGA (Level: Individual)
summary.dynEGA.Individuals <- function(object, ...) {
  cat("Number of Cases (individuals): \n")
  number <- length(object$dynEGA)
  print(number)
  cat("Summary statistics (number of factors/communities): \n")
  dim <- sapply(object$dynEGA, "[[", 3)
  cat("Mean:", mean(dim), "\n")
  cat("Median:", median(dim), "\n")
  cat("Min:", min(dim), "\n")
  cat("Max:", max(dim), "\n")
}

