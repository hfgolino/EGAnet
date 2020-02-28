#' Plot method for \code{\link[EGAnet]{dynEGA}}objects (Level: Group)
#'
#' Plots the \code{\link[EGAnet]{dynEGA}} result using \code{\link[qgraph]{qgraph}}
#'
#' @param x An \code{\link[EGAnet]{dynEGA}} object (Level: Group)
#'
#' @param ncol Number of columns
#'
#' @param nrow Number of rows
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
#' dyn.group <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, group = 22, use.derivatives = 1,
#' level = "group", model = "glasso")
#'
#' #Summary of dynEGA reults
#' summary(dyn.group)
#'
#' # Plot EGA network
#' plot(dyn.group, vsize = 6, label.prop = 1)
#'}
#'
#' @seealso \code{\link[EGAnet]{dynEGA}} to estimate the number of dimensions in multivariate time series using dynEGA.
#'
#' @export
#'
## S3 method for class 'dynEGA.Groups'
#'
#Plot dynEGA function (Level: Group)
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", vsize = 6,  ...) {
  par(mfrow=c(nrow,ncol))
  for(i in 1:length(x$dynEGA)){
    qgraph::qgraph(x$dynEGA[[i]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[i]]$wc), ...)
    title(names(x$dynEGA)[[i]], ...)}
}

