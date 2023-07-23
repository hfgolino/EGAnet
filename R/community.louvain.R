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
#' @param signed Boolean.
#' Whether the standard or signed algorithm should be used.
#' Defaults to \code{FALSE} or standard
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
#' community.louvain(network, signed = FALSE)
#' 
#' # Compute signed Louvain
#' community.louvain(network, signed = TRUE)
#'
#' @references
#' \strong{Louvain algorithm} \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' \strong{Signed modularity} \cr
#' Gomez, S., Jensen, P., & Arenas, A. (2009).
#' Analysis of community structure in networks of correlated data.
#' \emph{Physical Review E}, \emph{80}(1), 016114.
#'
#' @export
#'
# Compute Louvain communities for EGA
# Updated 23.07.2023
community.louvain <- function(
    network, signed = FALSE, 
    resolution = 1
)
{
  
  # Send errors (bookmark to have full input check)
  typeof_error(signed, "logical")
  
  # Get networks
  networks <- obtain_networks(network, signed)
  igraph_network <- networks$igraph_network
  network_matrix <- networks$network_matrix
  
  # Make sure there are variable names
  network_matrix <- ensure_dimension_names(network_matrix)
  
  # Obtain strength
  node_strength <- strength(network_matrix)
  
  # Initialize memberships as missing
  membership <- rep(NA, length(node_strength))
  
  # Determine unconnected nodes
  unconnected <- node_strength == 0
  
  # Determine whether all nodes are disconnected
  if(all(unconnected)){
    warning("The network input is empty. All community memberships are missing.", call. = FALSE)
  }else{ # Carry on if at least one node is connected
    
    # Check if any nodes are disconnected
    if(any(unconnected)){
      warning(
        "The network input contains unconnected nodes:\n",
        paste(names(node_strength)[unconnected], collapse = ", "),
        call. = FALSE
      )
    }
    
    # Algorithm function
    algorithm.FUN <- swiftelse(
      signed, signed.louvain, igraph::cluster_louvain
    )
    
    # Algorithm arguments
    algorithm.ARGS <- obtain_arguments(
      FUN = algorithm.FUN,
      FUN.args = list(resolution = resolution)
    )
    
    # Remove weights from igraph functions' arguments
    if("weights" %in% names(algorithm.ARGS)){
      algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
    }
    
    # Check for proper network
    if(signed){
      algorithm.ARGS[[1]] <- network_matrix
    }else{
      algorithm.ARGS[[1]] <- igraph_network
    }
    
    # Get result
    result <- do.call(algorithm.FUN, as.list(algorithm.ARGS))$memberships
    
  }
  
  # Name nodes
  dimnames(result)[[2]] <- dimnames(network_matrix)[[2]]
  
  # Return membership
  return(result)
  
}
