#' Plot method for \code{\link[EGAnet]{dynEGA}} objects (Fixed Effects)
#'
#' Plots the \code{\link[EGAnet]{dynEGA}} result using \code{\link[qgraph]{qgraph}}
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} object (Fixed Effects)
#'
#' @param title Character. Title of the plot
#'
#' @param vsize An integer indicating the size of the nodes.
#' Default vsize = 6
#'
#' @param id An integer or character indicating the ID of the individual to plot.
#'
#' @param ... Arguments passed to \code{\link[qgraph]{qgraph}}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Estimate dynEGA
#' dyn.intra <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' effects = "fixed", model = "glasso")
#'
#' #Summary of dynEGA reults
#' summary(dyn.intra)
#'
#' # Plot EGA network
#' plot(dyn.intra, vsize = 6, label.prop = 1, id = "ID1")
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA.Individuals'
#'
#Plot dynEGA function (Fixed Effects)
plot.dynEGA.Individuals <- function(x, title = "", vsize = 6,  id = NULL, ...) {
  plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)

}

