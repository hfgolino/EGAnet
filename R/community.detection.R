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
#' {See \code{\link[igraph]{cluster_leiden}} for more details.
#' \emph{Note}: The Leiden algorithm will default to the
#' Constant Potts Model objective function
#' (\code{objective_function = "CPM"}). Set
#' \code{objective_function = "modularity"} to use
#' modularity instead (see examples)}
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
#' @param signed Boolean.
#' Should signed network values be used?
#' Defaults to \code{FALSE}.
#' Caution should be used when setting to \code{TRUE}.
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
#' @param membership.only Boolean.
#' Whether the memberships only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' community detection algorithm
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
#' # Compute Leading Eigenvector
#' community.detection(network, algorithm = "leading_eigen")
#' 
#' # Compute Leiden (with modularity)
#' community.detection(
#'   network, algorithm = "leiden",
#'   objective_function = "modularity"
#' )
#' 
#' # Compute Louvain
#' community.detection(network, algorithm = "louvain")
#' 
#' # Compute Optimal (identifies maximum modularity solution)
#' community.detection(network, algorithm = "optimal")
#' 
#' # Compute Signed Louvain
#' community.detection(
#'   network,
#'   algorithm = "signed_louvain",
#'   signed = TRUE
#'   # needs to know the network *itself*
#'   # should be signed
#' )
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
# Updated 12.07.2023
community.detection <- function(
    network, algorithm = c(
      "edge_betweenness", "fast_greedy",
      "fluid", "infomap", "label_prop",
      "leading_eigen", "leiden", "louvain", "optimal",
      "signed_louvain", "spinglass", "walktrap"
    ),
    signed = FALSE, allow.singleton = FALSE,
    membership.only = TRUE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  
  # Get networks
  networks <- obtain_networks(network, signed)
  igraph_network <- networks$igraph_network
  network_matrix <- networks$network_matrix
  
  # Check for names
  network_matrix <- ensure_dimension_names(network_matrix)
  
  # Obtain network dimensions
  dimensions <- dim(network_matrix)
  
  # Obtain strength
  node_strength <- strength(network_matrix)
  
  # Initialize memberships as missing
  membership <- rep(NA, dimensions[2])
  
  # Determine whether all nodes are disconnected
  if(all(node_strength == 0)){
    
    # Send warning
    warning(
      "The network input is empty. All community memberships are missing."
    )
    
    # Set algorithm arguments
    algorithm.ARGS <- list()
    
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
    if("options" %in% names(algorithm.ARGS) & !"options" %in% names(ellipse)){
      algorithm.ARGS$options <- igraph::arpack_defaults
    }
    
    # Remove weights from igraph functions' arguments
    if("weights" %in% names(algorithm.ARGS)){
      algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
    }
    
    # Set up network
    if(is.function(algorithm) || algorithm != "signed_louvain"){
      algorithm.ARGS[[1]] <- igraph_network
    }else{ # Signed Louvain
      algorithm.ARGS[[1]] <- network_matrix
    }
    
    # Get result
    result <- do.call(algorithm.FUN, as.list(algorithm.ARGS))$membership
    
    # Obtain membership
    membership[!unconnected] <- result[!unconnected]
  
  }
  
  # Check singleton behavior
  if(isFALSE(allow.singleton)){
    
    # Determine whether there are any singleton communities
    membership_frequency <- table(membership)
    
    # Singletons
    singletons <- membership_frequency == 1
    
    # Check for frequencies equal to one
    if(any(singletons)){
      
      # Identify communities
      singleton_communities <- as.numeric(
        names(membership_frequency)[singletons]
      )
      
      # Set values to NA
      membership[membership %in% singleton_communities] <- NA
      
    }
    
  }
  
  # Name nodes
  names(membership) <- dimnames(network_matrix)[[2]]
  
  # Add methods to membership attributes
  attr(membership, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm), signed = signed,
    objective_function = algorithm.ARGS$objective_function
    # `objective_function` will be NULL unless it's there!
  )
  
  # Check for whether all results should be returned
  if(isTRUE(membership.only)){
    
    # Set class
    class(membership) <- "EGA.community"
    
    # Only return membership
    return(membership)
    
  }else{
    
    # Set up results
    results <- list(
      membership = membership,
      output = result
    )
    
    # Set class
    class(results) <- "EGA.community"
    
    # Return all results
    return(results)
    
  }
  
}

# Bug Checking ----
# ## Basic input
# network = network.estimation(wmt2[,7:24], model = "glasso");
# algorithm = "walktrap"; signed = FALSE;
# allow.singleton = FALSE; membership.only = TRUE;
# ellipse = list();

#' @exportS3Method 
# S3 Print Method ----
# Updated 05.07.2023
print.EGA.community <- function(x, boot = FALSE, ...)
{
  
  # Determine whether result is a list
  if(is.list(x)){
    membership <- x$membership
  }else{
    membership <- x
  }
  
  # Determine number of communities
  communities <- unique_length(membership)
  
  # Obtain algorithm name (if available)
  algorithm <- attr(membership, "methods")$algorithm
  
  # Determine whether algorithm was a function
  if(!is.function(algorithm)){
    
    # Check for signed
    algorithm_name <- ifelse(
      attr(membership, "methods")$signed,
      paste("Signed", algorithm),
      algorithm
    )
    
    # Check for Leiden
    if(algorithm == "Leiden"){
      
      # Obtain objective function
      objective_function <- attr(membership, "methods")$objective_function
      
      # Set up algorithm name
      objective_name <- ifelse(
        is.null(objective_function),
        "CPM", objective_function
      )
      
      # Expand "CPM"
      objective_name <- ifelse(
        objective_name == "CPM",
        "Constant Potts Model", "Modularity"
      )
      
      # Finalize algorithm name
      algorithm_name <- paste(
        algorithm, "with", objective_name
      )
      
    }
    
    # Set up methods
    cat(paste0("Algorithm: "), algorithm_name)
    
    # Check for bootEGA
    if(isFALSE(boot)){
      cat("\n\n") # add breakspace
    }

  }
  
  # Check for bootEGA
  if(isFALSE(boot)){
   
    # Print communities
    cat(paste0("Number of communities: "), communities)
    cat("\n\n") # Add breakspace
    
    # Remove class and attribute for clean print
    membership <- remove_attributes(membership)
    
    # Print membership
    print(membership)
     
  }
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 05.07.2023
summary.EGA.community <- function(object, boot = FALSE, ...)
{
  print(object, boot = boot, ...) # same as print
}

#' @noRd
# Obtain appropriate algorithm name ----
# Updated 22.06.2023
obtain_algorithm_name <- function(algorithm)
{
  
  # Set algorithm names (a hash table of sorts)
  algorithm_names <- c(
    # Directly from {igraph}'s `res$algorithm` output (for functions)
    "edge betweenness" = "Edge Betweenness", "fast greedy" = "Fast-greedy",
    "fluid communities" = "Fluid", "infomap" = "Infomap",
    "label propagation" = "Label Propagation",
    "leading eigenvector" = "Leading Eigenvector",
    "leiden" = "Leiden", "multi level" = "Louvain",
    "optimal" = "Optimal", "spinglass" = "Spinglass",
    "walktrap" = "Walktrap",
    # Algorithms with different characters from {EGAnet} (for characters)
    "edge_betweenness" = "Edge Betweenness", "fast_greedy" = "Fast-greedy",
    "fluid" = "Fluid", "label_prop" = "Label Propagation",
    "leading_eigen" = "Leading Eigenvector", "louvain" = "Louvain",
    "signed_louvain" = "Louvain"
  )
  
  # Check for function
  if(is.function(algorithm)){
    
    # Obtain function code
    function_code <- capture.output(algorithm)
    
    # Determine whether algorithm is {igraph}
    if(any(grepl("namespace:igraph", function_code))){
      
      # Proceed with identifying the algorithm
      algorithm_line <- function_code[
        grepl("res\\$algorithm", function_code)
      ]
      
      # Everything between \" and \"
      algorithm_name <- trimws(
        unlist(
          regmatches(
            algorithm_line, 
            gregexpr("(?<=\")[^\"]*(?=\")", algorithm_line, perl = TRUE)
          )
        )
      )
    
    }else{
      
      # Return algorithm for what it is
      ## Not expected to be used,
      # will likely break `community.detection` function
      return(algorithm)
      
    }
    
  }else if(is.character(algorithm)){
    
    # Set algorithm name as algorithm
    algorithm_name <- algorithm
    
  }
  
  # Obtain name from hash table
  return(unname(algorithm_names[algorithm_name]))
  
}

