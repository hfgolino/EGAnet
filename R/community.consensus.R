#' Applies the Consensus Clustering Method (\code{\link[EGAnet]{community.louvain}} only)
#'
#' Applies the consensus clustering method introduced by (Lancichinetti & Fortunato, 2012).
#' The original implementation of this method applies a community detection algorithm repeatedly
#' to the same network. With stochastic networks, the algorithm is likely to identify different
#' community solutions with many repeated applications.
#' 
#' The goal of the consensus clustering method is to identify a stable solution across
#' algorithm applications to derive a "consensus" clustering. The standard or "iterative"
#' approach is to apply the community detection algorithm \emph{N} times. Then, a co-occurrence
#' matrix is created representing how often each pair of nodes co-occurred across the
#' applications. Based on some cut-off value (e.g., 0.30), co-occurrences below this value
#' are set to zero, forming a "new" sparse network. The procedure proceeds until all nodes
#' co-occur with all other nodes in their community (or a proportion of 1.00).
#' 
#' Variations of this procedure are also available in this package but are
#' \strong{experimental}. Use these experimental procedures with caution.
#' More work is necessary before these experimental procedures are validated.
#' Please use the \code{"iterative"} approach before using other procedures
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
#' "Details" for more information. This method is the \strong{default}}
#' 
#' \item{\code{"most_common"}}
#' {\strong{EXPERIMENTAL.} Selects the community solution that appears the most
#' frequently across the applications. The idea behind this method is that the solution
#' that appears most often will be the most likely solution for the algorithm as well
#' as most reproducible. Can be less stable as the number of nodes increase requiring
#' a larger value for \code{consensus.iter}}
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
#' @param progress Boolean.
#' Whether progress should be printed.
#' Defaults to \code{TRUE}
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
# Updated 13.06.2023
community.consensus <- function(
    network, signed = FALSE, 
    order = c("lower", "higher"), resolution = 1,
    consensus.method = c(
      "highest_modularity", "iterative",
      "most_common", "lowest_tefi"
    ), consensus.iter = 1000, 
    correlation.matrix = NULL,
    progress = TRUE
)
{
  
  # Make order lower
  if(missing(order)){
    order <- "higher"
  }else{order <- tolower(match.arg(order))}
  
  # Set missing consensus method
  if(missing(consensus.method)){
    consensus.method <- "iterative"
  }else{consensus.method <- tolower(match.arg(consensus.method))}
  
  # Check for lowest TEFI method
  if(consensus.method == "lowest_tefi"){
    
    # Check for NULL correlation matrix
    if(is.null(correlation.matrix)){
      stop("A correlation matrix is required to compute TEFI. Please supply network's corresponding correlation matrix using the 'correlation.matrix' argument.")
    }else if(!is_symmetric(correlation.matrix)){ # Check for symmetric correlation matrix
      stop("Correlation matrix is not symmetric.")
    }else if(ncol(correlation.matrix) != ncol(network)){ # Check for same number of variables
      stop("Number of variables in 'correlation.matrix' does not match number of variables in 'network'. Double check that the correlation matrix matches the dimensions of the network.")
    }
    
  }
  
  # Determine class of network
  if(is(network, "igraph")){
    
    # Convert to network matrix
    network_matrix <- igraph2matrix(network)
    
    # Check for absolute
    if(!isTRUE(signed)){
      network_matrix <- abs(network_matrix)
    }
    
    # Convert to {igraph} network (ensures absolute even if {igraph})
    igraph_network <- convert2igraph(network_matrix)
    
    
  }else{
    
    # Ensure network is matrix
    network <- as.matrix(network)
    
    # Check for absolute
    if(!isTRUE(signed)){
      network <- abs(network)
    }
    
    # Store network as network matrix
    network_matrix <- network
    
    # Convert to {igraph} network
    igraph_network <- convert2igraph(network)
    
  }
  
  # Make sure there are variable names
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
    if(!isTRUE(signed)){
      algorithm.FUN <- igraph::cluster_louvain
    }else{
      algorithm.FUN <- signed.louvain
    }
    
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
    if(!isTRUE(signed)){
      algorithm.ARGS[[1]] <- igraph_network
    }else{
      algorithm.ARGS[[1]] <- network_matrix
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
      progress = progress
    )
    
    # Get result
    result <- do.call(consensus.FUN, as.list(consensus.ARGS))
    
  }
  
  # Return results
  return(result)
  
}

# Bug check ----
# network = ega.wmt$network; signed = FALSE;
# order = "higher"; resolution = 1;
# consensus.method = "lowest_tefi";
# consensus.iter = 1000; progress = TRUE;
# correlation.matrix = ega.wmt$correlation;

#' @noRd
# Standard application method ----
# Updated 27.05.2023
consensus_application <- function(
    FUN, FUN.ARGS, 
    consensus.iter, progress
)
{
  
  # Apply algorithm
  communities <- lapply(
    1:consensus.iter, function(i){
      
      # Perform call
      output <- do.call(what = FUN, args = FUN.ARGS)
      
      # Check for progress
      if(isTRUE(progress)){
        
        # Obtain text for progress
        progress_text <- paste0(
          "\r Consensus iteration ", formatC(
            i, digits = digits(consensus.iter) - 1,
            format = "d", flag = "0"
          ), " of ", consensus.iter,
          " complete."
        )
        
        # Print progress in message text
        cat(colortext(text = progress_text, defaults = "message"))
        
        # Close cat
        cat("\r")
        
      }
      
      # Return output
      return(output)
      
    }
  )
  
  # Return communities
  return(communities)
  
}

#' @noRd
# Highest modularity method ----
# Updated 17.05.2023
highest_modularity <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter, 
    correlation.matrix, # not used
    progress
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter,
    progress = progress
  )
  
  # Obtain modularities
  if(order == "lower"){
    modularities <- sapply(
      communities, function(x){
        x$modularity[1]
      }
    )
  }else if(order == "higher"){
    modularities <- sapply(
      communities, function(x){
        x$modularity[length(x$modularity)]
      }
    )
  }
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[1,]
      }
    ))
  }else if(order == "higher"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[nrow(x$memberships),]
      }
    ))
  }
  
  # Obtain solution
  selected_solution <- memberships[which.max(modularities),]
  
  # Set up return list
  results <- list(
    selected_solution = selected_solution,
    memberships = memberships,
    modularities = modularities
  )
  
  # Return results
  return(results)
  
}

#' @noRd
# Function to check for binary matrix ----
# Updated 27.05.2023
binary_check <- function(b_matrix){
  all(b_matrix == 0 | b_matrix == 1)
}

#' @noRd
# Iterative method ----
# Updated 27.05.2023
iterative <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter,
    correlation.matrix, # not used
    progress
)
{
  
  # Initialize iterations
  iterations <- 1
    
  # Loop over for consensus
  while(TRUE){
    
    # Apply algorithm
    communities <- consensus_application(
      FUN = FUN, FUN.ARGS = FUN.ARGS,
      consensus.iter = consensus.iter,
      progress = progress
    )
    
    # Obtain memberships
    if(order == "lower"){
      memberships <- t(sapply(
        communities, function(x){
          x$memberships[1,]
        }
      ))
    }else if(order == "higher"){
      memberships <- t(sapply(
        communities, function(x){
          x$memberships[nrow(x$memberships),]
        }
      ))
    }
    
    # Obtain maximum memberships
    maximum_memberships <- max(memberships, na.rm = TRUE)
    
    # Initialize consensus matrix
    consensus_matrix <- matrix(
      0, nrow = ncol(memberships),
      ncol = ncol(memberships)
    )
    
    # Loop over to get proportions
    for(i in 1:ncol(memberships)){
      for(j in i:ncol(memberships)){
        
        # Fill consensus matrix
        consensus_matrix[i,j] <- trace(
          table(
            factor( # Factor ensures proper table
              memberships[,i],
              levels = 1:maximum_memberships
            ), 
            factor(
              memberships[,j],
              levels = 1:maximum_memberships
            )
          )
        ) / consensus.iter
        
        # Fill other side
        consensus_matrix[j,i] <- consensus_matrix[i,j]
        
      }
    }
    
    # Threshold values
    consensus_matrix[consensus_matrix < 0.30] <- 0
    
    # Check for break
    if(isTRUE(binary_check(consensus_matrix))){
      break
    }
    
    # Increase iterations
    iterations <- iterations + 1
    
    # Update network (if continuing)
    if(is(FUN.ARGS[[1]], "igraph")){
      FUN.ARGS[[1]] <- convert2igraph(consensus_matrix)
    }
    
  }
  
  # Obtain solution
  selected_solution <- memberships[1,]
  
  # Set up return list
  results <- list(
    selected_solution = selected_solution,
    iterations = iterations
  )
  
  # Return results
  return(results)
  
}

#' @noRd
# Lowest TEFI method ----
# Updated 17.05.2023
lowest_tefi <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter, 
    correlation.matrix, # used in function
    progress
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter,
    progress = progress
  )
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[1,]
      }
    ))
  }else if(order == "higher"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[nrow(x$memberships),]
      }
    ))
  }
  
  # Make data frame
  memberships_df <- as.data.frame(memberships)
  
  # Obtain duplicate indices
  dupicated_index <- duplicated(memberships_df)
  
  # Rows for non-duplicates
  non_duplicated_df <- data.frame(memberships_df[!dupicated_index,])
  
  # Apply TEFI over non-duplicated memberships
  tefis <- sapply(
    1:nrow(non_duplicated_df), function(i){
      tefi(
        abs(correlation.matrix), 
        as.vector(as.matrix(non_duplicated_df[i,]))
      )$VN.Entropy.Fit
    }
  )
  
  # Obtain solution
  selected_solution <- non_duplicated_df[which.min(tefis),]
  
  # Set up return list
  results <- list(
    selected_solution = selected_solution,
    memberships = memberships,
    tefis = tefis
  )
  
  # Return results
  return(results)
  
}

#' @noRd
# Most Common method ----
# Updated 27.05.2023
most_common <- function(
    FUN, FUN.ARGS, 
    order, consensus.iter,
    correlation.matrix, # not used
    progress
)
{
  
  # Apply algorithm
  communities <- consensus_application(
    FUN = FUN, FUN.ARGS = FUN.ARGS,
    consensus.iter = consensus.iter,
    progress = progress
  )
  
  # Obtain memberships
  if(order == "lower"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[1,]
      }
    ))
  }else if(order == "higher"){
    memberships <- t(sapply(
      communities, function(x){
        x$memberships[nrow(x$memberships),]
      }
    ))
  }
  
  # Make data frame
  memberships_df <- as.data.frame(memberships)
  
  # Obtain duplicate indices
  dupicated_index <- duplicated(memberships_df)
  
  # Rows for non-duplicates
  non_duplicated_df <- data.frame(memberships_df[!dupicated_index,])
  
  # Rows for duplicates
  duplicated_df <- data.frame(memberships_df[dupicated_index,])
  
  # Match duplicates with non-duplicates
  duplicated_counts <- table(
    match(
      data.frame(t(duplicated_df)), data.frame(t(non_duplicated_df))
    ))
  
  # Obtain counts
  counts <- rep(1, nrow(non_duplicated_df))
  counts[as.numeric(names(duplicated_counts))] <- 
    counts[as.numeric(names(duplicated_counts))] + 
    duplicated_counts
  
  # Add counts to non-duplicated memberships
  proportion_table <- data.frame(
    Proportion = counts / consensus.iter,
    non_duplicated_df
  )
  
  # Obtain solution
  selected_solution <- proportion_table[
    which.max(proportion_table$Proportion), -1
  ]
  
  # Set up return list
  results <- list(
    selected_solution = selected_solution,
    memberships = memberships,
    proportion_table = proportion_table
  )
  
  # Return results
  return(results)
  
}

