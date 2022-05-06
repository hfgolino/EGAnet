#' Louvain Community Detection Algorithm
#'
#' @description Computes the Louvain community detection algorithm (Blondel et al., 2008)
#'
#' @param A Matrix or data frame.
#' A network adjacency matrix
#' 
#' @param method Character.
#' Whether modularity or \code{\link[EGAnet]{tefi}} should
#' be used to optimize communities.
#' Defaults to \code{"modularity"}
#'
#' @param resolution Numeric.
#' Resolution parameter for computing modularity.
#' Defaults to \code{1}.
#' Values smaller than 1 favor larger communities;
#' values larger than 1 favor smaller communities
#' 
#' @param corr Matrix or data frame.
#' Correlation matrix to be used when \code{method = "tefi"}
#' 
#' @details This version was adapted from the Matlab code available here:
#' https://perso.uclouvain.be/vincent.blondel/research/louvain.html. The code
#' was adjusted to mirror the results of \code{\link[igraph]{cluster_louvain}}.
#' The Louvain algorithm's results can vary depending on node ordering. In this
#' version, nodes are \strong{not} shuffled so that consistent results can be
#' achieved with the same node ordering. Results from \code{\link[igraph]{cluster_louvain}}
#' will shuffle nodes \strong{within} the function and therefore will sometimes produce
#' similar results and sometimes produce slightly different results. This version
#' is based all in R and therefore is slower than the version in \code{\link{igraph}}.
#'
#' @return Returns a list containing:
#'
#' \item{wc}{A matrix of lower to higher order community membership
#' detected in the network}
#'
#' \item{modularity}{A vector of modularity values corresponding the rows
#' of the \code{wc} matrix}
#'
#' @examples
#' # Load data
#' dep <- depression[,24:44]
#' 
#' # Estimate correlations
#' corr <- qgraph::cor_auto(dep)
#' 
#' # Estimate network
#' net <- EBICglasso.qgraph(corr, n = nrow(dep))
#'
#' # Estimate communities using modularity
#' louvain(net, method = "modularity")
#' 
#' # Estimate communities using tefi
#' louvain(net, method = "tefi", corr = corr)
#'
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s@virginia.edu>
#' 
#' @export
#'
# Louvain
# Updated 06.05.2022
louvain <- function(
    A,
    method = c("modularity", "tefi"),
    resolution = 1,
    corr = NULL
)
{
  
  # Check for missing arguments
  if(missing(method)){
    method <- "modularity"
  }else{
    method <- tolower(match.arg(method))
  }
  
  # Check for correlation matrix
  if(method == "tefi"){
    
    # Ensure matrix
    corr <- as.matrix(corr)
    
    # Ensure absolute
    corr <- abs(corr)
    
    # Missing correlation matrix
    if(missing(corr)){
      stop("Correlation matrix is necessary to compute TEFI")
    }
    
    # Symmetric matrix
    if(nrow(corr) != ncol(corr)){
      stop("Matrix input as correlation matrix does not have the same number of rows and columns")
    }
    
  }
  
  # Ensure matrix
  A <- as.matrix(A)
  
  # Set up results
  community_results <- matrix(
    0, nrow = ncol(A), ncol = ncol(A)
  )
  q_results <- numeric(ncol(A))
  
  # Set count
  count <- 1
  
  # Lower order ----
  results <- louvain_communities(
    newA = A,
    method = method,
    resolution = resolution,
    corr = corr
  )
  community_results[count,] <- results$communities
  q_results[count] <- results$modularity
  
  # Higher order ----
  
  # While loop
  while(TRUE){
    
    # Obtain new adjacency matrix with nodes merged
    newA <- make_higher_order(A, community_results[count,])
    
    # Increase count
    count <- count + 1
    
    # Obtain communities and modularity
    results <- louvain_communities(
      newA = newA,
      method = method,
      resolution = resolution,
      corr = corr,
      original_A = A,
      previous_communities = community_results[count - 1,],
      previous_modularity = q_results[count - 1]
    )
    community_results[count,] <- results$communities
    q_results[count] <- results$modularity
    
    # Check if communities match previous
    if(all(community_results[count,] == community_results[count - 1,])){
      break
    }
    
    
  }
  
  # Check results (checking for duplicates)
  duplicate_df <- as.data.frame(community_results)
  
  # Obtain duplicate indices
  dupe_ind <- duplicated(duplicate_df)
  
  # Rows for non-duplicates
  non_dupes <- as.matrix(duplicate_df[!dupe_ind,])
  q_results <- q_results[!dupe_ind]
  
  # Remove all zero rows
  keep_results <- apply(non_dupes, 1, function(x){
    ifelse(all(x == 0), FALSE, TRUE)
  })
  
  # Keep remaining
  community_results <- non_dupes[keep_results,]
  q_results <- q_results[keep_results]
  
  # Ensure matrix for community results
  if(!is.matrix(community_results)){
    
    community_results <- matrix(
      community_results, nrow = 1
    )
    
  }
  
  # Set up names
  if(!is.null(colnames(A))){
    
    # Name communities
    colnames(community_results) <- colnames(A)
    
  }
  
  # Set up results to return
  results <- list()
  results$wc <- community_results
  if(method == "modularity"){
    results$modularity <- round(q_results, 5)
  }else if(method == "tefi"){
    results$tefi <- round(q_results, 5)
  }
  
  # Return results
  return(results)
  
}
