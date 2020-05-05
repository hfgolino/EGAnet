#------------------------------------------
## S3Methods plot() // Updated 24.03.2020
#------------------------------------------

#' S3Methods for Plotting
#'
#' @name plots
#'
#' @aliases
#' plot.bootEGA
#' plot.CFA
#' plot.dynEGA
#' plot.dynEGA.Groups
#' plot.dynEGA.Individuals
#' plot.EGA
#' plot.NetLoads
#'
#' @usage
#' \method{plot}{bootEGA}(x, vsize = 6, ...)
#'
#' \method{plot}{CFA}(x, layout = "spring", vsize = 6, ...)
#'
#' \method{plot}{dynEGA}(x, title = "", vsize = 6,  ...)
#'
#' \method{plot}{dynEGA.Groups}(x, ncol, nrow, title = "", vsize = 6,  ...)
#'
#' \method{plot}{dynEGA.Individuals}(x, title = "", vsize = 6,  id = NULL, ...)
#'
#' \method{plot}{EGA}(x, title = "", vsize = 6,  ...)
#'
#' \method{plot}{NetLoads}(x, ...)
#'
#' @description Plots for \code{EGAnet} objects
#'
#' @param x Object from \code{EGAnet} package
#'
#' @param vsize Numeric.
#' Size of vertices in network plots.
#' Defaults to \code{6}
#'
#' @param layout Character.
#' Layout of plot (see \code{\link[semPlot]{semPaths}}).
#' Defaults to "spring"
#'
#' @param ncol Numeric.
#' Number of columns
#'
#' @param nrow Numeric.
#' Number of rows
#'
#' @param title Character.
#' Title of the plot.
#' Defaults to \code{""}
#'
#' @param id Numeric.
#' An integer or character indicating the ID of the individual to plot
#'
#' @param ... Arguments passed on to
#'
#' \itemize{
#'
#' \item{\code{\link[qgraph]{qgraph}}}
#' {Functions: bootEGA, dynEGA, dynEGA.Groups, dynEGA.Individuals, EGA, and net.loads}
#'
#' \item{\code{\link[semPlot]{semPaths}}}
#' {Functions: CFA}
#'
#' }
#'
#' @return Plots of \code{EGAnet} object
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics plot
#'
#Plot bootEGA----
# Updated 02.05.2020
#' @export
plot.bootEGA <- function(x, vsize = 6,...){

  qgraph::qgraph(x$typicalGraph$graph, layout = "spring",
                 groups = as.factor(x$typicalGraph$wc),
                 vsize = vsize, ...)

}

#Plot CFA----
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}

#Plot dynEGA function (Level: Group)----
# Updated 02.05.2020
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", vsize = 6,  ...) {
  par(mfrow=c(nrow,ncol))
  for(i in 1:length(x$dynEGA)){
    qgraph::qgraph(x$dynEGA[[i]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[i]]$wc), ...)
    title(names(x$dynEGA)[[i]], ...)}
}

#Plot dynEGA function (Level: Individual)----
# Updated 02.05.2020
#' @export
plot.dynEGA.Individuals <- function(x, title = "", vsize = 6,  id = NULL, ...) {
  plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)

}

#Plot dynEGA function (Level: Population)----
# Updated 02.05.2020
#' @export
plot.dynEGA <- function(x, title = "", vsize = 6,  ...) {
  plot.dynEGA <- qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA$wc), ...)

}

#Plot EGA----
# Updated 02.05.2020
#' @export
plot.EGA <- function(x, title = "", vsize = 6,  ...) {
  plot.ega <- qgraph::qgraph(x$network, layout = "spring", vsize = vsize, groups = as.factor(x$wc), ...)

}

#Plot net.loads----
# Updated 02.05.2020
#' @export
plot.NetLoads <- function(x, ...) {

  plot(x$plot)
}
