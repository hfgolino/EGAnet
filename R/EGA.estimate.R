#' A Wrapper Function for \code{EGA}
#'
#' Estimates the number of dimensions of a given dataset/instrument
#' using graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}})
#' method and the walktrap community detection algorithm (\code{\link[igraph]{cluster_walktrap}}).
#' The glasso regularization parameter is set via EBIC.
#'
#' @param data A dataframe with the variables to be used in the analysis or a correlation matrix.
#' If the data used is a correlation matrix, the argument \code{n} will need to be specified.
#'
#' @param n Integer.
#' Sample size, if the data provided is a correlation matrix
#'
#' @param model A string indicating the method to use.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param steps Number of steps to be used in \code{\link[igraph]{cluster_walktrap}} algorithm.
#' Defaults to 4.
#'
#' @param ... Additional arguments to be passed to \code{\link{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @return Returns a list containing:
#'
#' \item{estimated.network}{A symmetric network estimated using either the
#' \code{\link{EBICglasso.qgraph}} or \code{\link[NetworkToolbox]{TMFG}}}
#'
#' \item{wc}{A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network}
#'
#' \item{n.dim}{A scalar of how many total dimensions were identified in the network}
#'
#' \item{cor.data}{The zero-order correlation matrix}
#'
#' @examples
#'
#' \donttest{
#' #estimate EGA
#' ega.wmt <- EGA.estimate(data = wmt2[,7:24], model = "glasso")
#'
#' #estimate EGAtmfg
#' ega.wmt <- EGA.estimate(data = wmt2[,7:24], model = "TMFG")
#'
#' #estimate EGA
#' ega.intel <- EGA.estimate(data = intelligenceBattery[,8:66], model = "glasso")
#' }
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PloS one}, \emph{12(6)}, e0174035..
#' doi: \href{https://doi.org/10.1371/journal.pone.0174035}{journal.pone.0174035}
#'
#' Golino, H. F., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#' doi: \href{https://doi.org/10.1016/j.intell.2017.02.007}{j.intell.2017.02.007}
#'
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (in press).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}.
#' doi: \href{https://psyarxiv.com/gzcre/}{10.31234/osf.io/gzcre}
#'
#' @export
#'
# Estimates EGA
EGA.estimate <- function(data, n = NULL,
                         model = c("glasso", "TMFG"),
                         steps = 4, ...)
{
  # Check if model is missing
  if(missing(model))
  {model <- "glasso"
  }else{model <- match.arg(model)}

  # Check if data is correlation matrix and positive definite
  if(nrow(data) != ncol(data))
  {
    # Obtain n
    n <- nrow(data)

    # Compute correlation matrix
    cor.data <- qgraph::cor_auto(data, forcePD = TRUE)
  }else{

    # Check if positive definite
    if(any(eigen(data)$values < 0))
    {
      # Let user know
      warning("Correlation matrix is not positive definite.\nForcing positive definite matrix using Matrix::nearPD()\nResults may be unreliable")

      # Force positive definite matrix
      cor.data <- as.matrix(Matrix::nearPD(data)$mat)
    }else{cor.data <- data}
  }

  # Estimate network
  if(model == "glasso")
  {

    gamma.values <- c(0.50, 0.25, 0)

    for(j in 1:length(gamma.values))
    {
      estimated.network <- EBICglasso.qgraph(data = cor.data,
                                             n = n,
                                             lambda.min.ratio = 0.1,
                                             returnAllResults = FALSE,
                                             gamma = gamma.values[j],
                                             ...)

      if(all(NetworkToolbox::strength(estimated.network)>0))
      {
        message(paste("Network estimated with gamma = ",gamma.values[j],sep=""))
        break
      }
    }
  }else if(model == "TMFG")
  {estimated.network <- NetworkToolbox::TMFG(cor.data, ...)$A}

  # Convert to igraph
  graph <- suppressWarnings(NetworkToolbox::convert2igraph(abs(estimated.network)))

  # Check for unconnected nodes
  if(igraph::vcount(graph)!=ncol(data))
  {
    warning("Estimated network contains unconnected nodes:\n",
            paste(names(which(NetworkToolbox::strength(estimated.network)==0)), collapse = ", "))

    unconnected <- which(NetworkToolbox::strength(estimated.network)==0)
  }

  # Run walktrap
  wc <- igraph::walktrap.community(graph, steps = steps)

  # Obtain community memberships
  wc <- wc$membership
  init.wc <- as.vector(matrix(NA, nrow = 1, ncol = ncol(data)))
  init.wc[1:length(wc)] <- wc
  wc <- init.wc

  # Replace unconnected nodes with NA communities
  if(exists("unconnected"))
  {wc[unconnected] <- NA}

  names(wc) <- colnames(data)
  n.dim <- max(wc, na.rm = TRUE)

  # Return results
  res <- list()
  res$estimated.network <- estimated.network
  res$wc <- wc
  res$n.dim <- n.dim
  res$cor.data <- cor.data

  return(res)
}
