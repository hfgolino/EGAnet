#' A Sub-routine Function for \code{EGA}
#'
#' Estimates the number of dimensions of a given dataset or correlation matrix
#' using the graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[EGAnet]{TMFG}})
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
#' or \code{\link[EGAnet]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#' 
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
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
#' @param consensus.iter Numeric.
#' Number of iterations to perform in consensus clustering for the Louvain algorithm
#' (see Lancichinetti & Fortunato, 2012).
#' Defaults to \code{100}
#' 
#' @param consensus.method Character.
#' What consensus clustering method should be used? 
#' Defaults to \code{"highest_modularity"}.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\strong{\code{highest_modularity}}}
#' {Uses the community solution that achieves the highest modularity
#' across iterations}
#' 
#' \item{\strong{\code{most_common}}}
#' {Uses the community solution that is found the most
#' across iterations}
#' 
#' \item{\strong{\code{iterative}}}
#' {Identifies the most common community solutions across iterations
#' and determines how often nodes appear in the same community together.
#' A threshold of 0.30 is used to set low proportions to zero.
#' This process repeats iteratively until all nodes have a proportion of
#' 1 in the community solution.
#' }
#' 
#' \item{\code{lowest_tefi}}
#' {Uses the community solution that achieves the lowest \code{\link[EGAnet]{tefi}}
#' across iterations}
#' 
#' \item{\code{most_common_tefi}}
#' {Uses the most common number of communities detected across the number
#' of iterations. After, if there is more than one solution for that number
#' of communities, then the solution with the lowest \code{\link[EGAnet]{tefi}
#' is used}}
#' 
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
#' \code{\link{EBICglasso.qgraph}} or \code{\link[EGAnet]{TMFG}}}
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
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA.estimate(data = wmt)
#'
#' # Estimate EGAtmfg
#' ega.wmt.tmfg <- EGA.estimate(data = wmt, model = "TMFG")
#'
#' # Estimate EGA with Louvain algorithm
#' ega.wmt.louvain <- EGA.estimate(data = wmt, algorithm = "louvain")
#'
#' # Estimate EGA with Spinglass algorithm
#' ega.wmt.spinglass <- EGA.estimate(
#'   data = wmt,
#'   algorithm = igraph::cluster_spinglass # any {igraph} algorithm
#' )}
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' # Louvain algorithm \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#'
#' # Compared all \emph{igraph} community detections algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data \cr
#' Christensen, A. P., & Golino, H. (under review).
#' Estimating factors with psychometric networks: A Monte Carlo simulation comparing community detection algorithms.
#' \emph{PsyArXiv}.
#'
#' # Original simulation and implementation of EGA \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#'
#' Golino, H. F., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#'
#' # Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#'
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @export
#'
# Estimates EGA
# Updated 31.01.2023
EGA.estimate <- function(
    data, n = NULL,
    corr = c("cor_auto", "pearson", "spearman"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "leiden", "louvain"), algorithm.args = list(),
    consensus.method = c(
      "highest_modularity",
      "most_common",
      "iterative",
      "lowest_tefi",
      "most_common_tefi"
    ), consensus.iter = 100, 
  ...
)
{
  # Make the data a matrix
  data <- as.matrix(data)

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
  }else{model <- tolower(match.arg(model))}

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }
  
  if(missing(consensus.method)){
    consensus.method <- "most_common_tefi"
  }else{consensus.method <- tolower(match.arg(consensus.method))}

  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- tolower(match.arg(corr))}

  # Model function
  model.FUN <- switch(
    model,
    "glasso" = EBICglasso.qgraph,
    "tmfg" = TMFG
  )
  
  # Model arguments
  model.ARGS <- obtain.arguments(
    FUN = model.FUN,
    FUN.args = model.args
  )
  
  # Algorithm function
  if(!is.function(algorithm)){
    algorithm.FUN <- switch(
      algorithm,
      "walktrap" = igraph::cluster_walktrap,
      "leiden" = igraph::cluster_leiden,
      "louvain" = igraph::cluster_louvain
    )
  }else{
    algorithm.FUN <- algorithm
  }
  
  # Algorithm arguments
  algorithm.ARGS <- obtain.arguments(
    FUN = algorithm.FUN,
    FUN.args = algorithm.args
  )

  ## Remove weights from igraph functions' arguments
  if("weights" %in% names(algorithm.ARGS)){
    algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
  }

  # Check if data are data
  if(!isSymmetric(unname(as.matrix(data)))){
    
    # Obtain n
    n <- nrow(data)
    
    # Compute correlation matrix
    correlation <- switch(
      corr,
      "cor_auto" = auto.correlate(data),
      "pearson" = cor(data, use = "pairwise.complete.obs", method = "pearson"),
      "spearman" = cor(data, use = "pairwise.complete.obs", method = "spearman")
    )
    
  }else{ # Set correlation as data
    correlation <- data
  }
  
  # Check if positive definite
  if(any(eigen(correlation)$values < 0)){
    
    # Let user know
    warning("Correlation matrix is not positive definite.\nForcing positive definite matrix using Matrix::nearPD()\nResults may be unreliable")
    
    # Force positive definite matrix
    correlation <- as.matrix(Matrix::nearPD(correlation, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
    
  }
  
  # Add arguments to model
  model.ARGS$data <- correlation
  model.ARGS$n <- n
  
  # Estimate network
  if(model == "glasso"){
    
    # Check for lambda.min.ratio
    if(!"lambda.min.ratio" %in% names(model.args)){
      model.ARGS$lambda.min.ratio <- 0.1
    }

    # Check for gamma
    if(!"gamma" %in% names(model.args)){
      gamma.values <- seq(from = 0.50, to = 0, length.out = 3)
    }else{gamma.values <- model.args$gamma}

    # Loop through gamma values
    for(j in 1:length(gamma.values)){

      # Re-instate gamma values
      model.ARGS$gamma <- gamma.values[j]

      # Estimate network
      estimated.network <- do.call(
        model.FUN, model.ARGS
      )
      
      # Check for disconnected nodes
      if(all(abs(strength(estimated.network)) > 0)){
        break
      }
      
    }

  }else if(model == "tmfg"){
    estimated.network <- TMFG(correlation)$A
    colnames(estimated.network) <- colnames(correlation)
    rownames(estimated.network) <- rownames(correlation)
  }

  # Check for unconnected nodes
  if(all(degree(estimated.network) == 0)){

    # Initialize community membership list
    wc <- list()
    wc$membership <- rep(NA, ncol(estimated.network))
    warning(
      "Estimated network contains unconnected nodes:\n",
      paste(names(which(strength(estimated.network)==0)), collapse = ", ")
    )

    unconnected <- which(degree(estimated.network) == 0)

  }else{

    if(any(degree(estimated.network) == 0)){

      warning(
        "Estimated network contains unconnected nodes:\n",
        paste(names(which(strength(estimated.network)==0)), collapse = ", ")
      )

      unconnected <- which(degree(estimated.network) == 0)

    }
    
    # Check if algorithm is a function
    if(is.function(algorithm)){
      
      # Convert to igraph
      graph <- suppressWarnings(convert2igraph(abs(estimated.network)))
      
      # Run community detection algorithm
      algorithm.ARGS$graph <- graph
      
      # Call community detection algorithm
      wc <- do.call(
        what = algorithm.FUN,
        args = algorithm.ARGS
      )
      
    }else if(tolower(algorithm) == "louvain"){
      
      # Initialize community membership list
      wc <- list()
      
      # Check for lower order results
      if("lower.louvain" %in% names(add.args)){
        
        louvain.order <- ifelse(
          add.args$lower.louvain, "lower", "higher"
        )
        
      }else{
        louvain.order <- "higher"
      }
      
      # Population community membership list
      if(consensus.method != "most_common_tefi"){
        
        consensus <- consensus_clustering(
          network = estimated.network,
          corr = correlation,
          order = louvain.order,
          consensus.iter = consensus.iter,
          resolution = algorithm.ARGS$resolution,
          type = consensus.method
        )
        wc$membership <- consensus[[consensus.method]]
        
      }else{
        
        consensus <- most_common_tefi(
          network = estimated.network,
          corr = correlation,
          order = louvain.order,
          consensus.iter = consensus.iter,
          resolution = algorithm.ARGS$resolution
        )
        wc$membership <- consensus$most_common
        
      }
      
    }else{
      
      # Convert to igraph
      graph <- suppressWarnings(convert2igraph(abs(estimated.network)))
      
      # Run community detection algorithm
      algorithm.ARGS$graph <- graph
      
      # Call community detection algorithm
      wc <- do.call(
        what = algorithm.FUN,
        args = algorithm.ARGS
      )
      
    }

  }
  
  # Obtain community memberships
  wc <- wc$membership

  # Set up missing memberships
  init.wc <- as.vector(matrix(NA, nrow = 1, ncol = ncol(data)))
  init.wc[1:length(wc)] <- wc
  wc <- init.wc

  # Replace unconnected nodes with NA communities
  if(exists("unconnected")){
    wc[unconnected] <- NA
  }
  
  # Re-index communities
  wc <- suppressWarnings(
    reindex_comm(wc)
  )
  
  # Replace singleton communities with NA
  frequencies <- table(wc)
  
  # Check for singleton communities
  if(any(frequencies == 1)){
    
    # Singleton communities
    singletons <- as.numeric(names(frequencies)[which(frequencies == 1)])
    
    # Replace singletons with NA
    wc[!is.na(match(wc, singletons))] <- NA
  }

  # Name communities
  names(wc) <- colnames(data)
  
  # Obtain dimensions
  n.dim <- suppressWarnings(length(unique(na.omit(wc))))
  
  # Set methods arguments
  methods <- list(
    model = model,
    model.args = model.ARGS,
    algorithm = algorithm,
    algorithm.args = algorithm.ARGS,
    corr = corr,
    consensus.method = consensus.method,
    consensus.iter = consensus.iter
  )

  # Return results
  res <- list(
    network = estimated.network, wc = wc,
    n.dim = n.dim, cor.data = correlation,
    Methods = methods
  )
  
  # Check for consensus
  if(exists("consensus")){
    res$consensus <- consensus
  }

  return(res)
}

