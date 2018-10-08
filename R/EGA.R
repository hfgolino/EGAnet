#'  Apply the Exploratory Graph Analysis technique
#'
#' \code{EGA} Estimates the number of dimensions of a given dataset/instrument using graphical lasso or the TMFG method and a random walk algorithm. The glasso regularization parameter
#' is set via EBIC.
#'
#' @param data A dataframe with the variables to be used in the analysis, or a correlation matrix. If the data used is a correlation matrix, the arguments *matrix* and *n* will need to be specified.
#' @param plot.EGA Logical. If TRUE, returns a plot of the network and its estimated dimensions.
#' @param model A string indicating the method to use. Current options are:
#' -\code{glasso}:
#' {Gaussian graphical model estimation using graphical LASSO with extended Bayesian information criterion to select optimal regularization parameter (default method). Using \code{\link[qgraph]{EBICglasso}} from the \code{\link[qgraph]{qgraph}} package.}
#' \code{TMFG}:
#' {Estimates a Triangulated Maximally Filtered Graph, using the function \code{\link[NetworkToolbox]{TMFG}} of the \code{\link[NetworkToolbox]{NetworkToolbox}} package}
#' @param n Integer. Sample size, if the data provided is a correlation matrix.
#' @param steps Number of steps to be used in walktrap algorithm.
#' Defaults to 4
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso", plot.EGA = TRUE)
#' summary(ega.wmt)
#' plot(ega.wmt)
#'
#' ega.intel <- EGA(data = intelligenceBattery[,8:66], model = "glasso", plot.EGA = TRUE)
#' summary(ega.intel)
#' plot(ega.intel)
#'
#' \dontrun{
#' EGA(a)
#' }
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

# EGA default function - 11/21/2017
EGA <- function(data, model = c("glasso", "TMFG"), plot.EGA = TRUE, n = NULL, steps = 4) {
  if(!require(qgraph)) {
    message("installing the 'qgraph' package")
    install.packages("qgraph")
  }

  if(!require(igraph)) {
    message("installing the 'igraph' package")
    install.packages("igraph")
  }

  if(!require(NetworkToolbox)) {
    message("installing the 'NetworkToolbox' package")
    install.packages("NetworkToolbox")
  }
  {
    if(is.null(model)){
      model = "glasso"
    }
  if(!is.matrix(data)){
    data <- as.data.frame(data)
    if(model == "glasso"){
      cor.data <- cor_auto(data)
      estimated.network <- EBICglasso(S = cor.data, n = nrow(data), lambda.min.ratio = 0.1)
    } else if(model == "TMFG"){
      cor.data <- cor(data)
      estimated.network <- TMFG(data)$A

  }
} else if(is.matrix(data)){
    cor.data <- data
    if(model == "glasso"){
      estimated.network <- EBICglasso(S = data, n = n, lambda.min.ratio = 0.1)
    } else if(model == "TMFG"){
      estimated.network <- TMFG(data)$A
      }
    }
  }

    graph <- as.igraph(qgraph(abs(estimated.network), layout = "spring", vsize = 3, DoNotPlot = TRUE))
    wc <- walktrap.community(graph, steps = steps)
    names(wc$membership) <- colnames(data)
    n.dim <- max(wc$membership)
    a <- list()
    a$n.dim <- n.dim
    a$correlation <- cor.data
    a$network <- estimated.network
    a$wc <- wc$membership
    dim.variables <- data.frame(items = colnames(data), dimension = a$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]), ]
    a$dim.variables <- dim.variables
    class(a) <- "EGA"
    if (plot.EGA == TRUE) {
        plot.ega <- qgraph(estimated.network, layout = "spring", vsize = 6, groups = as.factor(wc$membership))
      }
    return(a)
}


