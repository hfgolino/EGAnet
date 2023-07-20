#' Applies the Consensus Clustering Method (\code{\link[EGAnet]{community.louvain}} only)
#'
#' Applies the consensus clustering method introduced by (Lancichinetti & Fortunato, 2012).
#' The original implementation of this method applies a community detection algorithm repeatedly
#' to the same network. With stochastic networks, the algorithm is likely to identify different
#' community solutions with many repeated applications.
#'
#' @param network Matrix or \code{\link{igraph}} network object
#' 
#' @param signed Boolean.
#' Whether the standard or signed algorithm should be used.
#' Defaults to \code{FALSE} or standard
#' 
#' @param order Character (length = 1).
#' Whether \code{"lower"} or \code{"higher"} order memberships from
#' the Louvain algorithm should be obtained for the consensus.
#' The \code{"lower"} order Louvain memberships are from the first
#' initial pass of the Louvain algorithm whereas the \code{"higher"}
#' order Louvain memberships are from the last pass of the Louvain
#' algorithm.
#' Defaults to \code{"higher"}
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
#' @param consensus.method Character (length = 1).
#' Available options for arriving at a consensus (\emph{Note}: 
#' All methods except \code{"iterative"} are considered experimental
#' until validated):
#' 
#' \itemize{
#' 
#' \item{\code{"highest_modularity"}}
#' {\strong{EXPERIMENTAL.} Selects the community solution with the highest modularity across
#' the applications. Modularity is a reasonable metric for identifying the number
#' of communities in a network but it comes with limitations (e.g., resolution limit)}
#' 
#' \item{\code{"iterative"}}
#' {The original approach proposed by Lancichinetti & Fortunato (2012). See
#' "Details" for more information}
#' 
#' \item{\code{"most_common"}}
#' {Selects the community solution that appears the most
#' frequently across the applications. The idea behind this method is that the solution
#' that appears most often will be the most likely solution for the algorithm as well
#' as most reproducible. Can be less stable as the number of nodes increase requiring
#' a larger value for \code{consensus.iter}.  This method is the \strong{default}}
#' 
#' \item{\code{"lowest_tefi"}}
#' {\strong{EXPERIMENTAL.} Selects the community solution with the lowest Total Entropy
#' Fit Index (\code{\link[EGAnet]{tefi}}) across the applications. TEFI is a reasonable metric
#' to identify the number of communities in a network based on Golino, Moulder et al. (2020)}
#' 
#' }
#' 
#' @param consensus.iter Numeric (length = 1).
#' Number of algorithm applications to the network.
#' Defaults to \code{1000}
#' 
#' @param correlation.matrix Symmetric matrix.
#' Used for computation of \code{\link[EGAnet]{tefi}}.
#' Only needed when \code{consensus.method = "tefi"}
#' 
#' @param membership.only Boolean.
#' Whether the memberships only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' community detection algorithm
#' 
#' @param ...
#' Not actually used but makes it either for general functionality
#' in the package
#' 
#' @details The goal of the consensus clustering method is to identify a stable solution across
#' algorithm applications to derive a "consensus" clustering. The standard or "iterative"
#' approach is to apply the community detection algorithm \emph{N} times. Then, a co-occurrence
#' matrix is created representing how often each pair of nodes co-occurred across the
#' applications. Based on some cut-off value (e.g., 0.30), co-occurrences below this value
#' are set to zero, forming a "new" sparse network. The procedure proceeds until all nodes
#' co-occur with all other nodes in their community (or a proportion of 1.00).
#' 
#' Variations of this procedure are also available in this package but are
#' \strong{experimental}. Use these experimental procedures with caution.
#' More work is necessary before these experimental procedures are validated
#' 
#' @return Returns a list containing...
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate correlation matrix
#' correlation.matrix <- auto.correlate(wmt)
#' 
#' # Estimate network
#' network <- EBICglasso.qgraph(data = wmt)
#' 
#' # Compute standard Louvain with highest modularity approach
#' community.consensus(
#'   network,
#'   consensus.method = "highest_modularity"
#' )
#' 
#' # Compute standard Louvain with iterative (original) approach
#' community.consensus(
#'   network,
#'   consensus.method = "iterative"
#' )
#' 
#' # Compute standard Louvain with most common approach
#' community.consensus(
#'   network,
#'   consensus.method = "most_common"
#' )
#' 
#' # Compute standard Louvain with lowest TEFI approach
#' community.consensus(
#'   network,
#'   consensus.method = "lowest_tefi",
#'   correlation.matrix = correlation.matrix
#' )
#' 
#' # Compute signed Louvain with highest modularity approach
#' community.consensus(
#'   network, signed = TRUE,
#'   consensus.method = "highest_modularity"
#' )
#' 
#' # Compute signed Louvain with iterative (original) approach
#' community.consensus(
#'   network, signed = TRUE,
#'   consensus.method = "iterative"
#' )
#' 
#' # Compute signed Louvain with most common approach
#' community.consensus(
#'   network, signed = TRUE,
#'   consensus.method = "most_common"
#' )
#' 
#' # Compute signed Louvain with lowest TEFI approach
#' community.consensus(
#'   network, signed = TRUE,
#'   consensus.method = "lowest_tefi",
#'   correlation.matrix = correlation.matrix
#' )
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
#' Lancichinetti, A., & Fortunato, S. (2012).
#' Consensus clustering in complex networks.
#' \emph{Scientific Reports}, \emph{2}(1), 1–7.
#' 
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @export
#'
# Compute consensus clustering for EGA
# Updated 14.07.2023
community.consensus <- function(
    network, signed = FALSE, 
    order = c("lower", "higher"), resolution = 1,
    consensus.method = c(
      "highest_modularity", "iterative",
      "most_common", "lowest_tefi"
    ), consensus.iter = 1000, 
    correlation.matrix = NULL,
    membership.only = TRUE,
    ...
)
{
  
  # Check for higher or lower order solution
  order <- set_default(order, "higher", community.consensus)
  consensus.method <- set_default(consensus.method, "most_common", community.consensus)
  
  # Get networks
  networks <- obtain_networks(network, signed)
  igraph_network <- networks$igraph_network
  network_matrix <- networks$network_matrix
  
  # Make sure there are variable names
  network_matrix <- ensure_dimension_names(network_matrix)
  
  # Obtain network dimensions
  dimensions <- dim(network_matrix)
  
  # Check for lowest TEFI method
  if(consensus.method == "lowest_tefi"){
    
    # Check for NULL correlation matrix
    if(is.null(correlation.matrix)){
      stop("A correlation matrix is required to compute TEFI. Please supply network's corresponding correlation matrix using the 'correlation.matrix' argument.")
    }else if(!is_symmetric(correlation.matrix)){ # Check for symmetric correlation matrix
      stop("Correlation matrix is not symmetric.")
    }else if(dim(correlation.matrix)[2] != dimensions[2]){ # Check for same number of variables
      stop("Number of variables in 'correlation.matrix' does not match number of variables in 'network'. Double check that the correlation matrix matches the dimensions of the network.")
    }
    
  }
  
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
    if(isTRUE(signed)){
      algorithm.FUN <- signed.louvain
    }else{
      algorithm.FUN <- igraph::cluster_louvain
    }
    
    # Algorithm arguments
    algorithm.ARGS <- obtain_arguments(
      FUN = algorithm.FUN,
      FUN.args = list(resolution = resolution)
    )
    
    # Remove weights from igraph functions' arguments
    if("weights" %in% names(algorithm.ARGS)){
      algorithm.ARGS["weights"] <- NULL
    }
    
    # Check for proper network
    if(isTRUE(signed)){
      algorithm.ARGS[[1]] <- network_matrix
    }else{
      algorithm.ARGS[[1]] <- igraph_network
    }
    
    # Get consensus method function
    consensus.FUN <- switch(
      consensus.method,
      "highest_modularity" = highest_modularity,
      "iterative" = iterative,
      "most_common" = most_common,
      "lowest_tefi" = lowest_tefi
    )
    
    # Set arguments
    consensus.ARGS <- list(
      FUN = algorithm.FUN,
      FUN.ARGS = as.list(algorithm.ARGS),
      order = order,
      consensus.iter = consensus.iter,
      correlation.matrix = correlation.matrix,
      dimensions = dimensions
    )
    
    # Get result
    result <- do.call(consensus.FUN, as.list(consensus.ARGS))
    
  }
  
  # Force into vector
  # `reindex_memberships` internal is in `community.detection`
  result$selected_solution <- reindex_memberships(
    force_vector(result$selected_solution)
  )
  
  # Obtain network names
  network_names <- dimnames(network_matrix)[[2]]
  
  # Ensure names
  names(result$selected_solution) <- network_names
  
  # Set methods attribute
  attr(result$selected_solution, "methods") <- list(
    algorithm = "Louvain",
    signed = signed, order = order, 
    consensus.method = consensus.method,
    consensus.iter = consensus.iter
  )
  
  # Check for membership only
  if(isTRUE(membership.only)){
    
    # Set class
    class(result$selected_solution) <- "EGA.consensus"
    
    # Return only membership
    return(result$selected_solution)
    
  }else{
    
    # Add names for returning full results
    ## Memberships
    if("memberships" %in% names(result)){
      dimnames(result$memberships)[[2]] <- network_names
    }
    
    ## Proportion Table
    if("proportion_table" %in% names(result)){
      dimnames(result$proportion_table)[[2]] <- c(network_names, "Proportion")
    }
    
    # Set class
    class(result) <- "EGA.consensus"
    
    # Return full results
    return(result)
    
  }
  
}

# Bug check ----
# network = ega.wmt$network; signed = FALSE;
# order = "higher"; resolution = 1;
# consensus.method = "most_common";
# consensus.iter = 1000;
# correlation.matrix = ega.wmt$correlation;

#' @exportS3Method 
# S3 Print Method ----
# Updated 29.06.2023
print.EGA.consensus <- function(x, ...)
{
  
  # Determine whether result is a list
  if(is.list(x)){
    membership <- x$selected_solution
  }else{
    membership <- x
  }
  
  # Obtain method
  method <- attr(membership, "methods")
  
  # Obtain consensus name
  consensus_name <- switch(
    method$consensus.method,
    "highest_modularity" = "Highest Modularity",
    "iterative" = "Iterative",
    "most_common" = "Most Common",
    "lowest_tefi" = "Lowest TEFI"
  )
  
  # Print method information
  cat(
    paste0(
      "Consensus Method: ", consensus_name,
      " (", method$consensus.iter, " iterations)",
      "\nAlgorithm: ", swiftelse(method$signed, "Signed Louvain", "Louvain"),
      "\nOrder: ", totitle(method$order)
    )
  )
  
  # Add breakspace
  cat("\n\n")
  
  # Determine number of communities
  communities <- unique_length(membership)
  
  # Print communities
  cat(paste0("Number of communities: "), communities)
  cat("\n\n") # Add breakspace
  
  # Remove attribute for clean print
  membership <- remove_attributes(membership)
  
  # Print membership
  print(membership)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 29.06.2023
summary.EGA.consensus <- function(object, ...)
{
  
  # Determine whether result is a list
  if(is.list(object)){
    membership <- object$selected_solution
  }else{
    membership <- object
  }
  
  # Obtain method
  method <- attr(membership, "methods")
  
  # Obtain consensus name
  consensus_name <- switch(
    method$consensus.method,
    "highest_modularity" = "Highest Modularity",
    "iterative" = "Iterative",
    "most_common" = "Most Common",
    "lowest_tefi" = "Lowest TEFI"
  )
  
  # Print method information
  cat(
    paste0(
      "Consensus Method: ", consensus_name,
      " (", method$consensus.iter, " iterations)",
      "\nAlgorithm: ", swiftelse(method$signed, "Signed Louvain", "Louvain"),
      "\nOrder: ", totitle(method$order)
    )
  )
  
  # Add breakspace
  cat("\n\n")
  
  # Determine number of communities
  communities <- unique_length(membership)
  
  # Print communities
  cat(paste0("Number of communities: "), communities)
  cat("\n\n") # Add breakspace
  
  # Remove attribute for clean print
  membership <- remove_attributes(membership)
  
  # Print membership
  print(membership)
  
}

#' @noRd
# Standard application method ----
# Updated 06.07.2023
consensus_application <- function(
    FUN, FUN.ARGS, consensus.iter
)
{

  # Apply algorithm
  return(
    lapply(
      seq_len(consensus.iter), do.call,
      what = FUN, args = FUN.ARGS
    )
  )
  
}

#' @noRd
# Highest modularity method ----
# Updated 01.07.2023
highest_modularity <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter, 
    correlation.matrix, # not used
    dimensions
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter
  )
  
  # Obtain modularities
  if(order == "lower"){
    modularities <- nvapply(
      communities, function(x){
        x$modularity[1]
      }
    )
  }else if(order == "higher"){
    modularities <- nvapply(
      communities, function(x){
        x$modularity[length(x$modularity)]
      }
    )
  }
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[1,]
      }, LENGTH = dimensions[2]
    ))
  }else if(order == "higher"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[dim(x$memberships)[1],]
      }, LENGTH = dimensions[2]
    ))
  }
  
  # Set up return list
  return(
    list(
      selected_solution = memberships[which.max(modularities),],
      memberships = memberships,
      modularities = modularities
    )
  )
  
}

#' @noRd
# Function to check for binary matrix ----
# Updated 01.07.2023
binary_check <- function(b_matrix)
{
  return(all(b_matrix == 0 | b_matrix == 1))
}

#' @noRd
# Iterative method ----
# Updated 06.07.2023
iterative <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter,
    correlation.matrix, # not used
    dimensions
)
{
  
  # Initialize iterations
  iterations <- 1
    
  # Loop over for consensus
  while(TRUE){
    
    # Apply algorithm
    communities <- consensus_application(
      FUN = FUN, FUN.ARGS = FUN.ARGS,
      consensus.iter = consensus.iter
    )
    
    # Obtain memberships
    if(order == "lower"){
      memberships <- t(nvapply(
        communities, function(x){
          x$memberships[1,]
        }, LENGTH = dimensions[2]
      ))
    }else if(order == "higher"){
      memberships <- t(nvapply(
        communities, function(x){
          x$memberships[dim(x$memberships)[1],]
        }, LENGTH = dimensions[2]
      ))
    }
    
    # Initialize consensus matrix
    consensus_matrix <- matrix(
      nrow = dimensions[2], ncol = dimensions[2]
    )
    
    # Loop over to get proportions
    for(i in seq_len(dimensions[2])){
      for(j in i:dimensions[2]){
        
        # Fill consensus matrix
        consensus_matrix[i,j] <- consensus_matrix[j,i] <- sum(
          memberships[,i] == memberships[,j],
          na.rm = TRUE
        )
        
      }
    }
    
    # Divide by iterations
    consensus_matrix <- consensus_matrix / consensus.iter
    
    # Threshold values
    consensus_matrix[consensus_matrix < 0.30] <- 0
    
    # Check for break
    if(binary_check(consensus_matrix)){
      break
    }
    
    # Increase iterations
    iterations <- iterations + 1
    
    # Update network (if continuing)
    if(is(FUN.ARGS[[1]], "igraph")){
      FUN.ARGS[[1]] <- convert2igraph(consensus_matrix)
    }
    
  }
  
  # Set up return list
  return(
    list(
      selected_solution = memberships[1,],
      iterations = iterations
    )
  )
  
}

#' @noRd
# Lowest TEFI method ----
# Updated 06.07.2023
lowest_tefi <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter, 
    correlation.matrix, # used in function
    dimensions
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter
  )
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[1,]
      }, LENGTH = dimensions[2]
    ))
  }else if(order == "higher"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[dim(x$memberships)[1],]
      }, LENGTH = dimensions[2]
    ))
  }
  
  # Prepare a data frame
  unique_table <- unique(memberships, MARGIN = 1)
  
  # Get dimensions
  dimensions <- dim(unique_table)
  
  # Apply TEFI over non-duplicated memberships
  tefis <- nvapply(
    seq_len(dimensions[1]), function(i){
      tefi(
        correlation.matrix, 
        unique_table[i,,drop = TRUE]
      )$VN.Entropy.Fit
    }
  )
  
  # Set up return list
  return(
    list(
      selected_solution = unique_table[which.min(tefis),],
      memberships = memberships,
      tefis = tefis,
      proportion_table = unique_table
    )
  )
  
}

#' @noRd
# Most Common method ----
# Updated 06.07.2023
most_common <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter,
    correlation.matrix, # not used
    dimensions
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter
  )
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[1,]
      }, LENGTH = dimensions[2]
    ))
  }else if(order == "higher"){
    memberships <- t(nvapply(
      communities, function(x){
        x$memberships[dim(x$memberships)[1],]
      }, LENGTH = dimensions[2]
    ))
  }
  
  # Prepare a data frame
  proportion_table <- count_table(memberships, proportion = TRUE)

  # Set up return list
  return(
    list(
      selected_solution = proportion_table[
        which.max(proportion_table$Value), seq_len(dimensions[2])
      ],
      memberships = memberships,
      proportion_table = proportion_table
    )
  )
  
}

