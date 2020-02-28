#' Plot method for \code{\link[EGAnet]{dynEGA}}objects (Level: Population)
#'
#' Plots the \code{\link[EGAnet]{dynEGA}} result using \code{\link[qgraph]{qgraph}}
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} object (Level: Population)
#'
#' @param title Character. Title of the plot
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
#' # Estimate dynEGA
#' dyn.random <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' level = "population", model = "glasso")
#'
#' #Summary of dynEGA reults
#' summary(dyn.random)
#'
#' # Plot EGA network
#' plot(dyn.random, vsize = 6, label.prop = 1)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA'
#'
#Plot dynEGA function (Level: Population)
plot.dynEGA <- function(x, title = "", vsize = 6,  ...) {
  plot.dynEGA <- qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA$wc), ...)

}

