#' @title Applies the Consensus Clustering Method (Louvain only)
#'
#' @description Applies the consensus clustering method introduced by (Lancichinetti & Fortunato, 2012).
#' The original implementation of this method applies a community detection algorithm repeatedly
#' to the same network. With stochastic networks, the algorithm is likely to identify different
#' community solutions with many repeated applications.
#'
#' @param network Matrix or \code{\link{igraph}} network object
#' 
#' @param order Character (length = 1).
#' Defaults to \code{"higher"}.
#' Whether \code{"lower"} or \code{"higher"} order memberships from
#' the Louvain algorithm should be obtained for the consensus.
#' The \code{"lower"} order Louvain memberships are from the first
#' initial pass of the Louvain algorithm whereas the \code{"higher"}
#' order Louvain memberships are from the last pass of the Louvain
#' algorithm
#' 
#' @param resolution Numeric (length = 1).
#' A parameter that adjusts modularity to allow the algorithm to
#' prefer smaller (\code{resolution} > 1) or larger
#' (0 < \code{resolution} < 1) communities.
#' Defaults to \code{1} (standard modularity computation)
#' 
#' @param consensus.method Character (length = 1).
#' Defaults to \code{"most_common"}.
#' Available options for arriving at a consensus (\emph{Note}: 
#' All methods except \code{"iterative"} are considered experimental
#' until validated):
#' 
#' \itemize{
#' 
#' \item{\code{"highest_modularity"} --- }
#' {\strong{EXPERIMENTAL.} Selects the community solution with the highest modularity across
#' the applications. Modularity is a reasonable metric for identifying the number
#' of communities in a network but it comes with limitations (e.g., resolution limit)}
#' 
#' \item{\code{"iterative"} --- }
#' {The original approach proposed by Lancichinetti & Fortunato (2012). See
#' "Details" for more information}
#' 
#' \item{\code{"most_common"} --- }
#' {Selects the community solution that appears the most
#' frequently across the applications. The idea behind this method is that the solution
#' that appears most often will be the most likely solution for the algorithm as well
#' as most reproducible. Can be less stable as the number of nodes increase requiring
#' a larger value for \code{consensus.iter}.  This method is the \strong{default}}
#' 
#' \item{\code{"lowest_tefi"} --- }
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
#' @param allow.singleton Boolean (length = 1).
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
#' @param ...
#' Not actually used but makes it easier for general functionality
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
#' \emph{At this time, seed setting for consensus clustering is not supported}
#' 
#' @return Returns either a vector with the selected solution
#' or a list when \code{membership.only = FALSE}:
#' 
#' \item{selected_solution}{Resulting solution from the consensus method}
#' 
#' \item{memberships}{Matrix of memberships across the consensus iterations}
#' 
#' \item{proportion_table}{For methods that use frequency, a table that
#' reports those frequencies alongside their corresponding memberships}
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
#' @references
#' \strong{Louvain algorithm} \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' \strong{Consensus clustering} \cr
#' Lancichinetti, A., & Fortunato, S. (2012).
#' Consensus clustering in complex networks.
#' \emph{Scientific Reports}, \emph{2}(1), 1–7.
#' 
#' \strong{Entropy fit indices} \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @export
#'
# Compute consensus clustering for EGA ----
# Updated 06.08.2023
community.consensus <- function(
    network, 
    order = c("lower", "higher"), resolution = 1,
    consensus.method = c(
      "highest_modularity", "iterative",
      "most_common", "lowest_tefi"
    ), consensus.iter = 1000, 
    correlation.matrix = NULL,
    allow.singleton = FALSE,
    membership.only = TRUE,
    ...
)
{
  
  # Check for higher or lower order solution
  order <- set_default(order, "higher", community.consensus)
  consensus.method <- set_default(consensus.method, "most_common", community.consensus)
  
  # Arguments errors
  community.consensus_errors(
    network, resolution, consensus.iter,
    correlation.matrix, allow.singleton, membership.only
  )
  
  # Check for {igraph} network
  if(is(network, "igraph")){
    network <- igraph2matrix(network)
  }
  
  # Use network matrix for now
  network <- abs(as.matrix(network))
  
  # Make sure there are variable names
  network <- ensure_dimension_names(network)
  
  # Obtain network dimensions
  dimensions <- dim(network)
  
  # Check for lowest TEFI method
  if(consensus.method == "lowest_tefi"){
    
    # Check for NULL correlation matrix
    if(is.null(correlation.matrix)){
      stop("A correlation matrix is required to compute TEFI. Please supply network's corresponding correlation matrix using the 'correlation.matrix' argument.", call. = FALSE)
    }else if(!is_symmetric(correlation.matrix)){ # Check for symmetric correlation matrix
      stop("Correlation matrix is not symmetric.", call. = FALSE)
    }else if(dim(correlation.matrix)[2] != dimensions[2]){ # Check for same number of variables
      stop("Number of variables in 'correlation.matrix' does not match number of variables in 'network'. Double check that the correlation matrix matches the dimensions of the network.", call. = FALSE)
    }
    
  }
  
  # Obtain strength
  node_strength <- colSums(network, na.rm = TRUE)
  
  # Get number of nodes
  nodes <- length(node_strength)
  
  # Initialize memberships as missing
  membership <- rep(NA, dimensions[2])
  
  # Determine unconnected nodes
  unconnected <- node_strength == 0
  
  # Determine whether all nodes are disconnected
  if(all(unconnected)){
    
    # Send warning about empty network
    warning("The network input is empty. All community memberships are missing.", call. = FALSE)
    
    # Set up results
    result <- list(selected_solution = rep(NA, nodes))
    
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
    algorithm.FUN <- igraph::cluster_louvain
    
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
    algorithm.ARGS[[1]] <- convert2igraph(network)
    
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
  result$selected_solution <- force_vector(result$selected_solution)
  
  # Obtain network names
  network_names <- dimnames(network)[[2]]
  
  # Ensure names
  names(result$selected_solution) <- network_names
  
  # Check singleton behavior
  if(!allow.singleton){
    
    # Determine whether there are any singleton communities
    membership_frequency <- fast_table(result$selected_solution)
    
    # Singletons
    singletons <- membership_frequency == 1
    
    # Check for frequencies equal to one
    if(any(singletons)){
      
      # Identify communities
      singleton_communities <- as.numeric(
        names(membership_frequency)[singletons]
      )
      
      # Set values to NA
      result$selected_solution[
        result$selected_solution %in% singleton_communities
      ] <- NA
      
    }
    
  }
  
  # `reindex_memberships` internal is in `community.detection`
  result$selected_solution <- reindex_memberships(result$selected_solution)
  
  # Set methods attribute
  attr(result$selected_solution, "methods") <- list(
    algorithm = "Louvain", order = order, 
    consensus.method = consensus.method,
    consensus.iter = consensus.iter
  )
  
  # Check for membership only
  if(membership.only){
    
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
# network = ega.wmt$network;
# order = "higher"; resolution = 1
# consensus.method = "most_common"
# consensus.iter = 1000; membership.only = TRUE
# correlation.matrix = ega.wmt$correlation; seed = NULL

#' @noRd
# Errors ----
# Updated 06.08.2023
community.consensus_errors <- function(
    network, resolution, consensus.iter,
    correlation.matrix, allow.singleton, membership.only
) 
{
  
  # 'network' errors
  if(!is(network, "igraph")){
    object_error(network, c("matrix", "data.frame"))
  }
  
  # 'resolution' errors
  length_error(resolution, 1)
  typeof_error(resolution, "numeric")
  range_error(resolution, c(0, Inf))
  
  # 'consensus.iter' errors
  length_error(consensus.iter, 1)
  typeof_error(consensus.iter, "numeric")
  range_error(consensus.iter, c(1, Inf))
  
  # 'correlation.matrix' errors
  if(!is.null(correlation.matrix)){
    object_error(correlation.matrix, c("matrix", "data.frame"))
  }
  
  # 'allow.singleton' errors
  length_error(allow.singleton, 1)
  typeof_error(allow.singleton, "logical")
  
  # 'membership.only' errors
  length_error(membership.only, 1)
  typeof_error(membership.only, "logical")
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 29.06.2023
print.EGA.consensus <- function(x, ...)
{
  
  # Determine whether result is a list
  membership <- swiftelse(is.list(x), x$selected_solution, x)
  
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
      "\nAlgorithm: Louvain",
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
  print(object, ...) # same as print
}

#' @noRd
# Standard application method ----
# Updated 02.08.2023
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
        unique_table[i,,drop = TRUE] # drop to vector
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

