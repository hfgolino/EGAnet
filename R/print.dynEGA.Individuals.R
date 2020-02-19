#' Print method for \code{\link[EGAnet]{dynEGA}} objects (Fixed Effects - Intraindividual Structure)
#'
#' Returns a summary of the \code{\link[EGAnet]{dynEGA}} objects (Fixed Effects - Intraindividual Structure)
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} objects (Fixed Effects - Intraindividual Structure)
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
#' effects = "fixed", model = "glasso")
#'
#' #Print dynEGA results
#' print(dyn.individual)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA.Individuals' (Fixed Effects - Intraindividual Structure)
#'
#Print dynEGA function
print.dynEGA.Individuals <- function(x, ...) {
  cat("Number of Cases (individuals): \n")
  number <- length(x$dynEGA)
  print(number)
  cat("Summary statistics (number of factors/communities): \n")
  dim <- sapply(x$dynEGA, "[[", 3)
  cat("Mean:", mean(dim), "\n")
  cat("Median:", median(dim), "\n")
  cat("Min:", min(dim), "\n")
  cat("Max:", max(dim), "\n")
}
