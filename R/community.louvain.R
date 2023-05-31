#' Apply the Louvain Community Detection Algorithm
#'
#' A specific function to apply the Louvain community detection algorithm. There
#' are two versions of the algorithm available: standard and signed. The standard
#' Louvain algorithm uses the absolute values of the network. The signed Louvain algorithm
#' uses the signed values of the network, which is made possible by using a signed version
#' of modularity (see Gomez et al., 2009).
#'
#' @param network Matrix or \code{\link{igraph}} network object
#' 
#' @param method Character (length = 1).
#' Whether \code{"standard"} or \code{"signed"} algorithm should be used.
#' Defaults to \code{"standard"}
#' 
#' @param resolution Numeric (length = 1).
#' A parameter that adjusts modularity to allow the algorithm to
#' prefer smaller (\code{resolution} > 1) or larger
#' (0 < \code{resolution} < 1) communities.
#' Defaults to \code{1} (standard modularity computation).
#' Currently, this argument is only available for \code{"standard"}.
#' Future versions may allow \code{"signed"} to take advantage of
#' this parameter
#' 
#' @return Returns each level of memberships from the Louvain algorithm
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate network
#' network <- EBICglasso.qgraph(data = wmt)
#' 
#' # Compute standard Louvain
#' community.louvain(network, method = "standard")
#' 
#' # Compute signed Louvain
#' community.louvain(network, method = "signed")
#'
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' Gomez, S., Jensen, P., & Arenas, A. (2009).
#' Analysis of community structure in networks of correlated data.
#' \emph{Physical Review E}, \emph{80}(1), 016114.
#'
#' @export
#'
# Compute Louvain communities for EGA
# Updated 30.05.2023
community.louvain <- function(
    network, method = c("signed", "standard"), 
    resolution = 1
)
{
  
  # Set missing method
  if(missing(method)){
    method <- "standard"
  }else{method <- tolower(match.arg(method))}
  
  # Determine class of network
  if(is(network, "igraph")){
    
    # Convert to network matrix
    network_matrix <- igraph2matrix(network)
    
    # Check for absolute
    if(method == "standard"){
      network_matrix <- abs(network_matrix)
    }
    
    # Convert to {igraph} network (ensures absolute even if {igraph})
    igraph_network <- convert2igraph(network_matrix)
    
    
  }else{
    
    # Ensure network is matrix
    network <- as.matrix(network)
    
    # Check for absolute
    if(method == "standard"){
      network <- abs(network)
    }
    
    # Store network as network matrix
    network_matrix <- network
    
    # Convert to {igraph} network
    igraph_network <- convert2igraph(network)
    
  }
  
  # Check for names
  if(is.null(colnames(network_matrix))){
    
    # Assign names
    names(network_matrix) <- paste0(
      "V", formatC(
        x = 1:ncol(network_matrix),
        digits = digits(ncol(network_matrix)) - 1,
        flag = "0", format = "d"
      )
    )
    
  }
  
  # Obtain strength
  node_strength <- colSums(abs(network_matrix), na.rm = TRUE)
  
  # Initialize memberships as missing
  membership <- rep(NA, length(node_strength))
  
  # Determine whether all nodes are disconnected
  if(all(node_strength == 0)){
    
    # Send warning
    warning(
      "The network input is empty. All community memberships are missing."
    )
    
  }else{ # Carry on if at least one node is connected
    
    # Check if any nodes are disconnected
    if(any(node_strength == 0)){
      
      # Determine unconnected nodes
      unconnected <- node_strength == 0
      
      # Send warning
      warning(
        "The network input contains unconnected nodes:\n",
        paste(names(node_strength)[unconnected], collapse = ", ")
      )
      
    }
    
    # Algorithm function
    algorithm.FUN <- switch(
      method,
      "standard" = igraph::cluster_louvain,
      "signed" = signed.louvain
    )
    
    # Algorithm arguments
    algorithm.ARGS <- obtain.arguments(
      FUN = algorithm.FUN,
      FUN.args = list(resolution = resolution)
    )
    
    # Remove weights from igraph functions' arguments
    if("weights" %in% names(algorithm.ARGS)){
      algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
    }
    
    # Check for proper network
    if(method == "standard"){
      algorithm.ARGS[[1]] <- igraph_network
    }else if(method == "signed"){
      algorithm.ARGS[[1]] <- network_matrix
    }
    
    # Get result
    result <- do.call(algorithm.FUN, as.list(algorithm.ARGS))$memberships
    
  }
  
  # Name nodes
  colnames(result) <- colnames(network_matrix)
  
  # Return membership
  return(result)
  
}

