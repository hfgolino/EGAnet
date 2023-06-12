#' Apply a Community Detection Algorithm
#'
#' General function to apply community detection algorithms available in
#' \code{\link{igraph}}. Follows the \code{\link{EGAnet}} approach of setting
#' singleton and disconnected nodes to missing (\code{NA})
#'
#' @param network Matrix or \code{\link{igraph}} network object
#' 
#' @param algorithm Character or \code{\link{igraph}} \code{cluster_*} function.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"edge_betweenness"}}
#' {See \code{\link[igraph]{cluster_edge_betweenness}} for more details}
#' 
#' \item{\code{"fast_greedy"}}
#' {See \code{\link[igraph]{cluster_fast_greedy}} for more details}
#' 
#' \item{\code{"fluid"}}
#' {See \code{\link[igraph]{cluster_fluid_communities}} for more details}
#' 
#' \item{\code{"infomap"}}
#' {See \code{\link[igraph]{cluster_infomap}} for more details}
#' 
#' \item{\code{"label_prop"}}
#' {See \code{\link[igraph]{cluster_label_prop}} for more details}
#' 
#' \item{\code{"leading_eigen"}}
#' {See \code{\link[igraph]{cluster_leading_eigen}} for more details}
#' 
#' \item{\code{"leiden"}}
#' {See \code{\link[igraph]{cluster_leiden}} for more details}
#' 
#' \item{\code{"louvain"}}
#' {See \code{\link[igraph]{cluster_louvain}} for more details}
#' 
#' \item{\code{"optimal"}}
#' {See \code{\link[igraph]{cluster_optimal}} for more details}
#' 
#' \item{\code{"signed_louvain"}}
#' {See \code{\link[EGAnet]{signed.louvain}} for more details}
#' 
#' \item{\code{"spinglass"}}
#' {See \code{\link[EGAnet]{cluster_spinglass}} for more details}
#' 
#' \item{\code{"walktrap"}}
#' {See \code{\link[EGAnet]{cluster_walktrap}} for more details}
#' 
#' }
#'
#' @param absolute Boolean.
#' Should absolute network values be used?
#' Defaults to \code{TRUE}.
#' Caution should be used when setting to \code{FALSE}.
#' Most algorithms are not able to handle signed (negative)
#' weights
#' 
#' @param allow.singleton Boolean.
#' Whether singleton or single node communities should be allowed.
#' Defaults to \code{FALSE}.
#' When \code{FALSE}, singleton communities will be set to
#' missing (\code{NA}); otherwise, when \code{TRUE}, singleton
#' communities will be allowed
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link{igraph}}'s community detection functions
#' (see \code{algorithm} for arguments for each algorithm)
#'
#' @return Returns memberships from a community detection algorithm
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
#' # Compute Edge Betweenness
#' community.detection(network, algorithm = "edge_betweenness")
#' 
#' # Compute Fast Greedy
#' community.detection(network, algorithm = "fast_greedy")
#' 
#' # Compute Fluid
#' community.detection(
#'   network, algorithm = "fluid",
#'   no.of.communities = 2 # needs to be set
#' )
#' 
#' # Compute Infomap
#' community.detection(network, algorithm = "infomap")
#' 
#' # Compute Label Propagation
#' community.detection(network, algorithm = "label_prop")
#' 
#' # Compute Leading Eigenvalue
#' community.detection(network, algorithm = "leading_eigen")
#' 
#' # Compute Louvain
#' community.detection(network, algorithm = "louvain")
#' 
#' # Compute Optimal (identifies maximum modularity solution)
#' community.detection(network, algorithm = "optimal")
#' 
#' # Compute Signed Louvain
#' community.detection(network, algorithm = "signed_louvain")
#' 
#' # Compute Spinglass
#' community.detection(network, algorithm = "spinglass")
#' 
#' # Compute Walktrap
#' community.detection(network, algorithm = "walktrap")
#' 
#' # Example with {igraph} network
#' community.detection(
#'   convert2igraph(network), algorithm = "walktrap"
#' )
#'
#' @references
#' Csardi, G., & Nepusz, T. (2006). The igraph software package for complex network research.
#' \emph{InterJournal, Complex Systems}, 1695.
#'
#' @export
#'
# Compute communities for EGA
# Updated 10.06.2023
community.detection <- function(
    network, algorithm = c(
      "edge_betweenness", "fast_greedy",
      "fluid", "infomap", "label_prop",
      "leading_eigen", "louvain", "optimal",
      "signed_louvain", "spinglass", "walktrap"
    ),
    absolute = TRUE, allow.singleton = FALSE,
    ...
)
{
  
  # Set default algorithm if missing
  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else{algorithm <- tolower(match.arg(algorithm))}
  
  # Determine class of network
  if(is(network, "igraph")){
    
    # Convert to network matrix
    network_matrix <- igraph2matrix(network)
    
    # Check for absolute
    if(isTRUE(absolute)){
      network_matrix <- abs(network_matrix)
    }
    
    # Convert to {igraph} network (ensures absolute even if {igraph})
    igraph_network <- convert2igraph(network_matrix)
    
    
  }else{
    
    # Ensure network is matrix
    network <- as.matrix(network)
    
    # Check for absolute
    if(isTRUE(absolute)){
      network <- abs(network)
    }
    
    # Store network as network matrix
    network_matrix <- network
    
    # Convert to {igraph} network
    igraph_network <- convert2igraph(network)
    
  }
  
  # Check for names
  network_matrix <- ensure_dimension_names(network_matrix)
  
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
    
    # Determine unconnected nodes
    unconnected <- node_strength == 0
    
    # Check if any nodes are disconnected
    if(any(unconnected)){
      
      # Send warning
      warning(
        "The network input contains unconnected nodes:\n",
        paste(names(node_strength)[unconnected], collapse = ", ")
      )
      
    }
    
    # Algorithm function
    if(!is.function(algorithm)){
      
      algorithm.FUN <- switch(
        tolower(algorithm),
        "edge_betweenness" = igraph::cluster_edge_betweenness,
        "fast_greedy" = igraph::cluster_fast_greedy,
        "fluid" = igraph::cluster_fluid_communities,
        "infomap" = igraph::cluster_infomap,
        "label_prop" = igraph::cluster_label_prop,
        "leading_eigen" = igraph::cluster_leading_eigen,
        "leiden" = igraph::cluster_leiden,
        "louvain" = igraph::cluster_louvain,
        "optimal" = igraph::cluster_optimal,
        "signed_louvain" = signed.louvain,
        "spinglass" = igraph::cluster_spinglass,
        "walktrap" = igraph::cluster_walktrap
      )
      
    }else{
      
      # Set algorithm otherwise
      algorithm.FUN <- algorithm
      
    }
    
    # Obtain ellipse arguments
    ellipse <- list(...)
    
    # Algorithm arguments
    algorithm.ARGS <- obtain_arguments(
      FUN = algorithm.FUN,
      FUN.args = ellipse
    )
    
    # Check for Leading Eigenvalue (needs ARPACK)
    if(algorithm == "leading_eigen" & !"options" %in% names(ellipse)){
      algorithm.ARGS$options <- igraph::arpack_defaults
    }
    
    # Remove weights from igraph functions' arguments
    if("weights" %in% names(algorithm.ARGS)){
      algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
    }
    
    # Set up network
    if(tolower(algorithm) == "signed_louvain"){
      
      # Add {igraph} network
      algorithm.ARGS[[1]] <- network_matrix
      
      # Get result
      result <- do.call(algorithm.FUN, as.list(algorithm.ARGS))$memberships
      
      # Obtain membership (higher-order)
      membership[!unconnected] <- result[nrow(result),]
      
    }else{
      
      # Add {igraph} network
      algorithm.ARGS[[1]] <- igraph_network
      
      # Get result
      result <- do.call(algorithm.FUN, as.list(algorithm.ARGS))$membership
      
      # Obtain membership
      membership[!unconnected] <- result[!unconnected]
      
    }
  
  }
  
  # Check singleton behavior
  if(!isTRUE(allow.singleton)){
    
    # Determine whether there are any singleton communities
    membership_frequency <- table(membership)
    
    # Check for frequencies equal to one
    if(any(membership_frequency == 1)){
      
      # Identify communities
      singleton_communities <- as.numeric(
        names(membership_frequency)[
          membership_frequency == 1
        ]
      )
      
      # Set values to NA
      membership[
        membership %in% singleton_communities
      ] <- NA
      
    }
    
  }
  
  # Name nodes
  names(membership) <- colnames(network_matrix)
  
  # Return membership
  return(membership)
  
}

