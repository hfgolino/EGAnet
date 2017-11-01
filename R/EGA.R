#'  Apply the Exploratory Graph Analysis technique
#'
#' \code{EGA} Estimates the number of dimensions of a given dataset/instrument using graphical lasso and a random walk algorithm. The glasso regularization parameter
#' is set via EBIC.
#'
#' @param data A dataframe with the variables to be used in the analysis, or a correlation matrix. If the data used is a correlation matrix, the arguments *matrix* and *n* will need to be specified.
#' @param plot.EGA Logical. If TRUE, returns a plot of the network of partial correlations estimated via graphical lasso and its estimated dimensions.
#' @param matrix Logical. If TRUE, will treat the data as a correlation matrix, and a *n* (sample size) will need to be specified. Default to FALSE.
#' @param n Integer. Sample size, if the data provided is a correlation matrix.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' summary(ega.wmt)
#' plot(ega.wmt)
#'
#' ega.intel <- EGA(data = intelligenceBattery[,8:66])
#' summary(ega.intel)
#' plot(ega.intel)
#'
#' \dontrun{
#' EGA(a)
#' }
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

# EGA default function - 01/18/2017
EGA <- function(data, plot.EGA = TRUE, matrix = FALSE, n = NULL) {
  if(!require(qgraph)) {
    message("installing the 'qgraph' package")
    install.packages("qgraph")
  }

  if(!require(igraph)) {
    message("installing the 'igraph' package")
    install.packages("igraph")
  }

  if(matrix == FALSE){
    data <- as.data.frame(data)
    cor.data <- cor_auto(data)
    glasso.ebic <- EBICglasso(S = cor.data, n = nrow(data), lambda.min.ratio = 0.1)
    graph.glasso <- as.igraph(qgraph(abs(glasso.ebic), layout = "spring", vsize = 3, DoNotPlot = TRUE))
    wc <- walktrap.community(graph.glasso)
    n.dim <- max(wc$membership)

    if (plot.EGA == TRUE) {
      plot.ega <- qgraph(glasso.ebic, layout = "spring", vsize = 6, groups = as.factor(wc$membership))
    }

    a <- list()
    a$n.dim <- n.dim
    a$correlation <- cor.data
    a$glasso <- glasso.ebic
    a$wc <- wc$membership
    dim.variables <- data.frame(items = colnames(data), dimension = a$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]), ]
    a$dim.variables <- dim.variables
    class(a) <- "EGA"
    return(a)
  } else{
    cor.data <- data
    glasso.ebic <- EBICglasso(S = cor.data, n = n, lambda.min.ratio = 0.1)
    graph.glasso <- as.igraph(qgraph(abs(glasso.ebic), layout = "spring", vsize = 3, DoNotPlot = TRUE))
    wc <- walktrap.community(graph.glasso)
    n.dim <- max(wc$membership)

    if (plot.EGA == TRUE) {
      plot.ega <- qgraph(glasso.ebic, layout = "spring", vsize = 6, groups = as.factor(wc$membership))
    }

    a <- list()
    a$n.dim <- n.dim
    a$correlation <- data
    a$glasso <- glasso.ebic
    a$wc <- wc$membership
    dim.variables <- data.frame(items = colnames(data), dimension = a$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]), ]
    a$dim.variables <- dim.variables
    class(a) <- "EGA"
    return(a)
  }
}

