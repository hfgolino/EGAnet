#' A Sub-routine Function for \code{EGA}
#'
#' Estimates the number of dimensions of a given dataset or correlation matrix
#' using the graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}})
#' network estimation methods.
#' 
#' Two community detection algorithms, Walktrap (Pons & Latapy, 2006) and
#' Louvain (Blondel et al., 2008), are pre-programmed because of their
#' superior performance in simulation studies on psychological
#' data generated from factor models (Christensen & Golino; 2020; Golino et al., 2020).
#' Notably, any community detection algorithm from the \code{\link{igraph}}
#' can be used to estimate the number of communities (see examples).
#'
#' @param data Matrix or data frame.
#' Variables (down columns) or correlation matrix.
#' If the input is a correlation matrix,
#' then argument \code{n} (number of cases) is \strong{required}
#'
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
#'
#' @param model Character.
#' A string indicating the method to use.
#' 
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
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#' 
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#' 
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and Hudson Golino <hfg9s at virginia.edu>
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
#' \dontshow{# Fast for CRAN checks
#' # Pearson's correlation matrix
#' wmt <- cor(wmt2[,7:24])
#' 
#' # Estimate EGA
#' ega.wmt <- EGA.estimate(data = wmt, n = nrow(wmt2), model = "glasso")
#' 
#' }
#'
#' \donttest{
#' # Estimate EGA
#' ega.wmt <- EGA.estimate(data = wmt2[,7:24], model = "glasso")
#'
#' # Estimate EGAtmfg
#' ega.wmt <- EGA.estimate(data = wmt2[,7:24], model = "TMFG")
#' 
#' # Estimate EGA with Spinglass
#' ega.wmt <- EGA.estimate(data = wmt2[,7:24], model = "glasso",
#' algorithm = igraph::cluster_spinglass)
#' }
#' 
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' # Louvain algorithm \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#' doi: \href{https://doi.org/10.1088/1742-5468/2008/10/P10008}{10.1088/1742-5468/2008/10/P10008}
#' 
#' # Compared all \emph{igraph} community detections algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data \cr
#' Christensen, A. P., & Golino, H. (under review).
#' Estimating factors with psychometric networks: A Monte Carlo simulation comparing community detection algorithms.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/hz89e}{10.31234/osf.io/hz89e}
#' 
#' # Original simulation and implementation of EGA \cr
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
#' # Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#' doi: \href{https://doi.org/10.1037/met0000255}{10.1037/met0000255}
#' 
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#' doi: \href{https://doi.org/10.7155/jgaa.00185}{10.7155/jgaa.00185}
#'
#' @export
#'
# Estimates EGA
# Updated 21.10.2020
EGA.estimate <- function(data, n = NULL,
                         model = c("glasso", "TMFG"), model.args = list(),
                         algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                         corr = c("cor_auto", "pearson", "spearman"), ...)
{
  
  # Get additional arguments
  add.args <- list(...)
  
  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      paste(
        "The 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")",
        sep = ""
      )
    )

    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }
  
  #### ARGUMENTS HANDLING ####

  # Missing arguments
  
  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }

  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  # Model Arguments
  ## Check for model
  if(model == "glasso"){
    model.formals <- formals(EBICglasso.qgraph)
  }else{model.formals <- formals(NetworkToolbox::TMFG)}
  
  ## Check for input model arguments
  if(length(model.args) != 0){
    
    ### Check for matching arguments
    if(any(names(model.args) %in% names(model.formals))){
      
      model.replace.args <- model.args[na.omit(match(names(model.formals), names(model.args)))]
  
      model.formals[names(model.replace.args)] <- model.replace.args
    }
    
  }
  
  ## Remove ellipses
  if("..." %in% names(model.formals)){
    model.formals[which(names(model.formals) == "...")] <- NULL
  }
  
  # Algorithm Arguments
  ## Check for algorithm
  if(!is.function(algorithm)){
    
    if(algorithm == "walktrap"){
      algorithm.formals <- formals(igraph::cluster_walktrap)
    }else if(algorithm == "louvain"){
      algorithm.formals <- formals(igraph::cluster_louvain)
    }
    
  }else{algorithm.formals <- formals(algorithm)}
  
  ## Check for input algorithm arguments
  if(length(algorithm.args) != 0){
    
    ### Check for matching arguments
    if(any(names(algorithm.args) %in% names(algorithm.formals))){
      
      algorithm.replace.args <- algorithm.args[na.omit(match(names(algorithm.formals), names(algorithm.args)))]
      
      algorithm.formals[names(algorithm.replace.args)] <- algorithm.replace.args
    }
    
  }
  
  ## Remove ellipses
  if("..." %in% names(algorithm.formals)){
    algorithm.formals[which(names(algorithm.formals) == "...")] <- NULL
  }
  
  ## Remove weights from igraph functions' arguments
  if("weights" %in% names(algorithm.formals)){
    algorithm.formals[which(names(algorithm.formals) == "weights")] <- NULL
  }

  #### ARGUMENTS HANDLING ####

  # Check if data is correlation matrix and positive definite
  if(nrow(data) != ncol(data)){
    
    # Obtain n
    n <- nrow(data)

    # Compute correlation matrix

    cor.data <- switch(corr,
                       cor_auto = qgraph::cor_auto(data, forcePD = TRUE),
                       pearson = cor(data, use = "pairwise.complete.obs", method = "pearson"),
                       spearman = cor(data, use = "pairwise.complete.obs", method = "spearman")
    )

    # Check if positive definite
    if(any(eigen(cor.data)$values < 0)){
      
      # Let user know
      warning("Correlation matrix is not positive definite.\nForcing positive definite matrix using Matrix::nearPD()\nResults may be unreliable")

      # Force positive definite matrix
      cor.data <- as.matrix(Matrix::nearPD(cor.data, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
    }
    
  }else{

    # Check if positive definite
    if(any(eigen(data)$values < 0)){
      
      # Let user know
      warning("Correlation matrix is not positive definite.\nForcing positive definite matrix using Matrix::nearPD()\nResults may be unreliable")

      # Force positive definite matrix
      cor.data <- as.matrix(Matrix::nearPD(data, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
      
    }else{cor.data <- data}
  }
  
  #### ADDITIONAL ARGUMENTS HANDLING ####
  
  model.formals$data <- cor.data
  model.formals$n <- n
  
  #### ADDITIONAL ARGUMENTS HANDLING ####

  # Estimate network
  if(model == "glasso")
  {
    # GLASSO additional arguments
    ## Lambda
    if(!"lambda.min.ratio" %in% names(model.args)){
      model.formals$lambda.min.ratio <- 0.1
    }
    
    if(!"gamma" %in% names(model.args)){
      gamma.values <- seq(from = 0.50, to = 0, length.out = 3)
    }else{gamma.values <- model.formals$gamma}
    
    for(j in 1:length(gamma.values)){
      
      # Re-instate gamma values
      model.formals$gamma <- gamma.values[j]
      
      # Estimate network
      estimated.network <- do.call(EBICglasso.qgraph, model.formals)

      if(all(abs(NetworkToolbox::strength(estimated.network))>0)){
        
        message(paste("Network estimated with:\n",
                      " \u2022 gamma = ", gamma.values[j], "\n",
                      " \u2022 lambda.min.ratio = ", model.formals$lambda.min.ratio,
                      sep=""))
        break
      }
    }
    
  }else if(model == "TMFG"){
    estimated.network <- NetworkToolbox::TMFG(cor.data)$A
  }

  # Convert to igraph
  graph <- suppressWarnings(NetworkToolbox::convert2igraph(abs(estimated.network)))

  # Check for unconnected nodes
  if(igraph::vcount(graph)!=ncol(data)){
    
    warning("Estimated network contains unconnected nodes:\n",
            paste(names(which(NetworkToolbox::strength(estimated.network)==0)), collapse = ", "))

    unconnected <- which(NetworkToolbox::degree(estimated.network)==0)
    
  }

  # Run community detection algorithm
  algorithm.formals$graph <- graph
  
  if(!is.function(algorithm)){
    
    wc <- switch(algorithm,
                 walktrap = do.call(igraph::cluster_walktrap, as.list(algorithm.formals)),
                 louvain = do.call(igraph::cluster_louvain, as.list(algorithm.formals))
    )
    
  }else{wc <- do.call(what = algorithm, args = as.list(algorithm.formals))}

  # Obtain community memberships
  wc <- wc$membership
  init.wc <- as.vector(matrix(NA, nrow = 1, ncol = ncol(data)))
  init.wc[1:length(wc)] <- wc
  wc <- init.wc

  # Replace unconnected nodes with NA communities
  if(exists("unconnected")){
    wc[unconnected] <- NA
  }

  names(wc) <- colnames(data)
  n.dim <- max(wc, na.rm = TRUE)

  # Return results
  res <- list()
  res$network <- estimated.network
  res$wc <- wc
  res$n.dim <- n.dim
  res$cor.data <- cor.data
  
  if(model == "glasso"){
    res$gamma <- model.formals$gamma
    res$lambda <- model.formals$lambda.min.ratio
  }

  return(res)
}