#' Applies Several Different Approaches to Detect Unidimensional Communities
#'
#' A function to apply several approaches to detect a unidimensional community in 
#' networks. There have many different approaches recently such as expanding
#' the correlation matrix to have orthogonal correlations (\code{"expand"}),
#' applying the Leading Eigenvalue community detection algorithm
#' \code{\link[igraph]{cluster_leading_eigen}} to the correlation matrix
#' (\code{"LE"}), and applying the Louvain community detection algorithm
#' \code{\link[igraph]{cluster_louvain}} to the correlation matrix (\code{"louvain"})
#'
#' @param data Numeric matrix or data frame.
#' Either data representing \emph{only} the variables of interest, or
#' a correlation matrix. Data that are not numeric will be
#' removed from the dataset
#' 
#' @param n Numeric (length = 1).
#' Sample size if \code{data} is a correlation matrix
#' 
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"} to automatically compute
#' appropriate correlations using \code{\link[EGAnet]{auto.correlate}}.
#' \code{"pearson"} and \code{"spearman"} are provide for completeness.
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"}}
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"}}
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"}}
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"}}
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"}}
#' {Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details}
#' 
#' }
#' 
#' @param uni.method Character (length = 1).
#' Unidimensional method to apply to the data.
#' Defaults to \code{"louvain"} based on recent evidence (Christensen, 2023).
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"expand"}}
#' {Expands the correlation matrix by four variables that are all
#' correlated 0.50 with each other and 0.00 to the empirical data
#' (Golino et al., 2020)}
#' 
#' \item{\code{"LE"}}
#' {Applies the Leading Eigenvalue algorithm \code{\link[igraph]{cluster_leading_eigen}}
#' to the correlation matrix (Christensen et al., 2023)}
#' 
#' \item{\code{"louvain"}}
#' {Applies the Louvain algorithm with consensus clustering
#' \code{\link[EGAnet]{community.consensus}} to the correlation matrix (Christensen, 2023)}
#' 
#' }
#' 
#' @param verbose Boolean.
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#' 
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}}, \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.consensus}}, and \code{\link[EGAnet]{community.detection}}
#' 
#' @return Returns the memberships of the community detection algorithm.
#' The memberships will output \emph{regardless} of whether the
#' network is unidimensional
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Louvain with Consensus Clustering (default)
#' community.unidimensional(wmt)
#' 
#' # Leading Eigenvalue
#' community.unidimensional(wmt, uni.method = "LE")
#'
#' # Expand
#' community.unidimensional(wmt, uni.method = "expand")
#'
#' @references
#' \strong{Expand approach} \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020).
#' Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors:
#' A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320. 
#' 
#' \strong{Leading Eigenvalue approach} \cr
#' Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (2023).
#' Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation.
#' \emph{Behavior Research Methods}.
#' 
#' \strong{Louvain approach} \cr
#' Christensen, A. P. (2023).
#' Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison.
#' \emph{PsyArXiv}.
#' 
#' @export
#'
# Compute unidimensional approaches for EGA
# Updated 16.06.2023
community.unidimensional <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    uni.method = c("expand", "LE", "louvain"),
    verbose = FALSE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", community.unidimensional)
  na.data <- set_default(na.data, "pairwise", auto.correlate)  
  model <- set_default(model, "glasso", network.estimation)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  
  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Make sure there are variable names
  data <- ensure_dimension_names(data)
  
  # Check for whether data are data or correlation matrix
  if(is_symmetric(data)){
    
    # Check for sample size
    if(is.null(n)){
      stop("A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.")
    }

    # Set data as correlation matrix
    correlation_matrix <- data
    
  }else{ # Assume 'data' is data
    
    # Check for appropriate variables
    data <- usable_data(data, verbose)
    
    # Obtain sample size
    n <- nrow(data)
    
    # Check for automatic correlations
    if(corr == "auto"){
      
      # Obtain arguments for `auto.correlate`
      auto_ARGS <- obtain_arguments(FUN = auto.correlate, FUN.args = ellipse)
      
      # Supply data and correlation method
      auto_ARGS$data <- data
      auto_ARGS$corr <- "pearson"
      
      # Obtain correlation matrix
      correlation_matrix <- do.call(
        what = auto.correlate,
        args = auto_ARGS
      )
      
    }else{
      
      # Obtain correlations using base R
      correlation_matrix <- cor(data, use = na.data, method = corr)
      
    }
    
  }
  
  # Apply unidimensional approach
  if(uni.method == "expand"){
    
    # Perform "expand" approach
    membership <- expand(
      correlation_matrix, n, model,
      verbose, ellipse
    )
    
  }else if(uni.method == "le"){
    
    # Perform Leading Eigenvalue approach (pass arguments on)
    membership <- community.detection(
      correlation_matrix, algorithm = "leading_eigen", ellipse
    )
    
  }else if(uni.method == "louvain"){
    
    # Perform Louvain with Consensus approach
    membership <- consensus_wrapper(
      correlation_matrix, verbose, ellipse
    )
    
  }
  
  # Add methods attribute
  attr(membership, "methods") <- list(
    corr = corr, model = model,
    uni.method = uni.method
  )
  
  # No S3 methods -- not intended for individual use
  
  # Return membership
  return(membership)
  
}

# Bug Checking ----
# ## Basic input
# data = wmt2[,7:24]; n = NULL;
# corr = "auto"; na.data = "pairwise";
# model = "glasso"; uni.method = "expand";
# verbose = FALSE; ellipse = list()

#' @noRd
# "Expand" Correlation approach ----
# Updated 13.06.2023
expand <- function(correlation_matrix, n, model, verbose, ellipse)
{
  
  # Number of variables
  variables <- ncol(correlation_matrix)
  
  # New total variables
  new_total <- variables + 4
  
  # Create new matrix
  expanded_matrix <- matrix(
    0, nrow = new_total, ncol = new_total
  )
  
  # Create orthogonal correlation matrix
  orthogonal_matrix <- matrix(0.50, nrow = 4, ncol = 4)
  diag(orthogonal_matrix) <- 1
  
  # Insert original correlation matrix
  expanded_matrix[
    seq_len(variables), seq_len(variables)
  ] <- correlation_matrix
  
  # Insert orthogonal correlations
  expanded_matrix[
    (variables + 1):new_total,
    (variables + 1):new_total
  ] <- orthogonal_matrix
  
  # Obtain network estimation arguments
  network_ARGS <- obtain_arguments(network.estimation, ellipse)
  
  # Set data, n, and model
  network_ARGS$data <- expanded_matrix
  network_ARGS$n <- n
  network_ARGS$model <- model
  
  # Apply network estimation method
  network <- do.call(network.estimation, network_ARGS)
  
  # Obtain community detection arguments
  community_ARGS <- obtain_arguments(community.detection, ellipse)
  
  # Set network
  community_ARGS$network <- network

  # Apply community detection algorithm
  membership <- do.call(community.detection, community_ARGS)
  
  # Remove additional variables
  membership <- membership[seq_len(variables)]
  
  # Add back names
  names(membership) <- colnames(correlation_matrix)
  
  # Return membership
  return(membership)
  
}

#' @noRd
# Wrapper for Louvain consensus ----
# Updated 15.06.2023
consensus_wrapper <- function(correlation_matrix, verbose, ellipse)
{
  
  # Obtain arguments
  consensus_ARGS <- obtain_arguments(community.consensus, ellipse)
  
  # Add network
  consensus_ARGS$network <- correlation_matrix
  
  # Set progress
  consensus_ARGS$progress <- verbose
  
  # Apply Louvain with consensus approach
  membership <- do.call(
    what = community.consensus,
    args = consensus_ARGS
  )$selected_solution
  
  # Add back names
  names(membership) <- colnames(correlation_matrix)
  
  # Return membership
  return(membership)
  
}

