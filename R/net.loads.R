#' Network Loadings
#'
#' @description Computes the between- and within-community
#' \code{strength} of each item
#' for each community. This function uses the
#' \code{comcat} and
#' \code{stable} functions to calculate
#' the between- and within-community strength of each item, respectively.
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' A network adjacency matrix
#'
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#' 
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure.
#' Defaults to \code{\link[GPArotation]{geominQ}}.
#' For a list of rotations, see \code{\link{GPArotation}}
#' 
#' @param min.load Numeric.
#' Sets the minimum loading allowed in the standardized
#' network loading matrix. Values equal or greater than
#' the minimum loading are kept in the output. Values
#' less than the minimum loading are removed. This matrix can
#' be viewed using \code{print()} or \code{summary()}
#' Defaults to \code{0}
#'
#' @return Returns a list containing:
#'
#' \item{unstd}{A matrix of the unstandardized within- and between-community
#' strength values for each node}
#'
#' \item{std}{A matrix of the standardized within- and between-community
#' strength values for each node}
#' 
#' \item{minLoad}{The minimum loading to appear in summary of network loadings.
#' Use \code{print()} or \code{summary()} to view}
#'
#' @details Simulation studies have demonstrated that a node's strength
#' centrality is roughly equivalent to factor loadings
#' (Christensen, Golino, & Silvia, 2019; Hallquist, Wright, & Molenaar, in press).
#' Hallquist and colleagues (in press) found that node strength represented a
#' combination of dominant and cross-factor loadings. This function computes
#' each node's strength within each specified dimension, providing a rough
#' equivalent to factor loadings (including cross-loadings).
#'
#' For more details, type \code{vignette("Network_Scores")}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#'
#' # Network loadings
#' net.loads(ega.wmt)
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   ega.wmt,
#'   stats = "net.loads"
#' )}
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#'
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. M. (2019).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}, 1-25.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings
# Updated 10.07.2023
# Default = "BRM" or `net.loads` from version 1.2.3
# Experimental = new signs and cross-loading adjustment
net.loads <- function(
    A, wc, loading.method = c("BRM", "experimental"),
    rotation = "geominQ", ...
)
{
  
  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  loading.method <- set_default(loading.method, "brm", net.loads)
  rotation <- set_default(rotation, "geominq", net.loads)
  
  # Organize and extract input
  input <- organize_input(A, wc)
  A <- input$A; wc <- input$wc
  
  # Get number of nodes and their names
  nodes <- length(wc); node_names <- names(wc)
  
  # Get unique communities (`NA` is OK)
  unique_communities <- sort(unique(wc)) # put in order
  
  # Get number of communities
  communities <- length(unique_communities)
  
  # Get return order of node names and communities (without NA)
  return_node_order <- order(node_names)
  return_community_order <- order(
    unique_communities[unique_communities != "NA"]
  )
  
  # If all singleton communities, then send NA for all
  if(nodes == communities){
    
    # Send results
    return(
      list(
        unstd = unstandardized[return_node_order, return_community_order],
        std = unstandardized[return_node_order, return_community_order]
      )
    )
    
  }
  
  # Not singleton dimensons, so carry on
  
  # Check for method
  if(loading.method == "brm"){
    
    # Compute unstandardized loadings (absolute sums)
    unstandardized <- absolute_weights(A, wc, nodes, unique_communities)
    
    # Add signs to the loadings
    unstandardized <- old_add_signs(unstandardized, A, wc, unique_communities)
    
    # Before rounding occured prior to standardization (no rounding)
    standardized <- t(
      t(unstandardized) / sqrt(colSums(abs(unstandardized), na.rm = TRUE))
    )
    
    # Get descending order
    standardized <- descending_order(standardized, wc, unique_communities)
    
    # Set up results
    results <- list(
      unstd = unstandardized[dimnames(standardized)[[1]],],
      std = standardized
    )
    
    # Add attributes
    attr(results, "methods") <- list(
      loading.method = loading.method, rotation = rotation
    )
    
    # Add class
    class(results) <- "net.loads"
    
    # Return results
    return(results)
    
  }else{ # Pass on to experimental
    
    # Initialize loading matrix
    loading_matrix <- matrix(
      0, nrow = nodes, ncol = communities,
      dimnames = list(node_names, unique_communities)
    )
    
    # Initialize sign vector
    signs <- rep(1, nodes)
    names(signs) <- node_names
    
    # Populate loading matrix
    for(community in unique_communities){
      
      # Get community index
      community_index <- wc == community
    
      # Determine positive direction for dominant loadings
      target_network <- obtain_signs(A[community_index, community_index, drop = FALSE])
      
      # Compute absolute sum for dominant loadings
      loading_matrix[community_index, community] <- colSums(target_network, na.rm = TRUE)
      
      # Determine positive direction for dominant loadings
      signs[community_index] <- attr(target_network, "signs")
      
    }
    
    # Check for unidimensional structure
    if(communities > 1){
      
      # Check for any negative signs
      if(any(signs == -1)){
        
        # Make a copy of the network
        A_copy <- A
        
        # Flip signs
        A[signs == -1,] <- A[,signs == -1] <- A_copy[signs == -1,]
        
        
      }
      
      # Populate loading matrix with cross-loadings
      for(community in unique_communities){
        for(cross in unique_communities){
          
          # No need for same community loadings
          if(community != cross){
            
            # Get community index
            community_index <- wc == community
            
            # Compute algebraic sum for cross-loadings
            loading_matrix[community_index, cross] <- colSums(
              A[wc == cross, community_index, drop = FALSE], na.rm = TRUE
            )
            
          }
          
        }
      }
      
    }
    
    # Set signs
    loading_matrix <- loading_matrix * signs

    # Obtain standardized loadings
    standardized <- t(
      t(loading_matrix) / sqrt(colSums(abs(loading_matrix), na.rm = TRUE))
    )
    
    
    
    
    
  }
  
  
  # Check for singleton communities
  if(length_wc == length(unique_wc)){
    
    # Initialize results
    unstd <- matrix(NA, nrow = ncol(A), ncol = ncol(A))
    colnames(unstd) <- colnames(A)
    row.names(unstd) <- colnames(A)
    
    # Set up results
    results <- list(
      unstd = unstd,
      std = unstd
    )
    
  }else{ # Not singleton dimensions
    
    # Reorder communities
    wc <- wc[wc_order]
    
    # Reorder network
    A <- A[wc_order, wc_order]
    
    # Initialize loading matrix
    loading_matrix <- matrix(
      0, nrow = ncol(A),
      ncol = length(unique_wc)
    )
    
    # Initialize sign vector
    signs <- rep(1, ncol(A)) # start with all positive orientation
    names(signs) <- colnames(A)
    
    # Add column and row names
    row.names(loading_matrix) <- colnames(A)
    colnames(loading_matrix) <- unique_wc
    
    # Populate loading matrix
    for(dominant in unique_wc){
      
      # Obtain target portion of network
      target_network <- A[wc == dominant, wc == dominant]
      
      # Determine positive direction for dominant loadings
      sign_updated <- obtain_signs(target_network)
      
      # Update the target network
      target_network <- sign_updated$target_network
      
      # Obtain the sum
      target_sum <- sum_function(target_network)
      
      # Compute absolute sum for dominant loadings
      loading_matrix[wc == dominant, as.character(dominant)] <- 
        target_sum # / sqrt(sum(abs(target_sum)))
      
      # Determine positive direction for dominant loadings
      signs[wc == dominant] <- sign_updated$signs
      
    }
    
    # Check for cross-loadings
    if(length(unique_wc) > 1){
      
      # Initialize reversed A
      A_reversed <- A
      
      # Create duplicate of network
      if(sum.method == "signed"){
        
        if(any(signs == -1)){
          A_reversed[which(signs == -1),] <- -A[which(signs == -1),]
          A_reversed[,which(signs == -1)] <- -A[,which(signs == -1)]
        }
        
      }
      
      # Populate loading matrix
      for(dominant in unique_wc){
        for(cross in unique_wc){
          
          # Do not use dominant loadings
          if(dominant != cross){
            
            # Obtain the sum
            target_sum <- sum_function(A_reversed[wc == cross, wc == dominant])
            
            # Compute algebraic sum for cross-loadings
            loading_matrix[wc == dominant, as.character(cross)] <- 
              target_sum # / sqrt(sum(abs(target_sum)))
            
          }
          
        }
      }
      
      # Check for sum method absolute
      if(sum.method == "absolute"){
        
        # Add signs (old way)
        loading_matrix <- old.add.signs(
          comm.str = loading_matrix,
          A = A, wc = wc, dims = unique_wc
        )
        
      }
      
    }
    
    # Set signs
    loading_matrix <- loading_matrix * signs
    
    # Check for flipping orientation
    if(isTRUE(positive.orientation)){
      
      # Using signs, ensure positive orientation based
      # on most common direction
      for(dominant in unique_wc){
        
        # Determine dominant orientation
        orientation <- sum(signs[wc == dominant])
        
        # Check for negative orientation
        if(orientation <= -1){
          
          # Reverse dominant variables signs across all communities
          loading_matrix[wc == dominant,] <-
            -loading_matrix[wc == dominant,]
          
          # Check for cross-loadings
          if(length(unique_wc) > 1){
            
            # Reverse cross-loading signs on target community
            loading_matrix[wc != dominant, as.character(dominant)] <-
              -loading_matrix[wc != dominant, as.character(dominant)]
            
          }
          
        }
        
      }
      
    }
    
    # Obtain standardized loadings
    standardized <- t(
      t(loading_matrix) /
        sqrt(colSums(abs(loading_matrix)))
    )
    
    # Set up for rotation
    
    # Check for {GPArotation} and {fungible}
    # Function in `helpers-general.R`
    check_package(c("GPArotation", "fungible"))
    
    # Obtain rotation from GPArotation package
    rotation_names <- ls(asNamespace("GPArotation"))
    
    # Check if rotation exists
    rotation_names_lower <- tolower(rotation_names)
    
    # Obtain rotation arguments
    rot_arguments <- list(...)
    
    # Check if rotation exists
    if(tolower(rotation) %in% rotation_names_lower){
      
      if(tolower(rotation) != "oblimin"){
        
        # Obtain arguments
        rotation_arguments <- obtain.arguments(
          FUN = psych::faRotations, FUN.args = list(rotate = rotation)
        )
        
        # Check for arguments
        rotation_arguments$loadings <- standardized
        rotation_arguments$n.rotations <- ifelse(
          "n.rotations" %in% names(rot_arguments),
          rot_arguments$n.rotations,
          10
        )
        rotation_arguments$maxit <- ifelse(
          "maxit" %in% names(rot_arguments),
          rot_arguments$maxit,
          1000
        )
        
        # Add other arguments
        rotation_arguments <- c(
          rotation_arguments, rot_arguments[
            which(!names(rot_arguments) %in% names(rotation_arguments))
          ]
        )
        
        # Add default for "geominQ"
        if(tolower(rotation) == "geominq"){
          
          # Check for epsilon
          if(!"eps" %in% names(rot_arguments)){
            
            # Set up standard >= 4 dimensions
            eps <- 0.01
            
            # Check for 2 or 3 dimensions
            eps <- ifelse(ncol(standardized) == 2, 0.0001, eps)
            eps <- ifelse(ncol(standardized) == 3, 0.001, eps)
            
            # Set up defaults
            rotation_arguments$eps <- eps
            
          }
          
        }
        
        # Set loadings
        rotation_arguments$loadings <- as.matrix(standardized)
        
        # Obtain rotated loadings
        rotated <- do.call(
          what = psych::faRotations,
          args = as.list(rotation_arguments)
        )
        
      }else{
        
        # Obtain arguments
        rotation_arguments <- obtain.arguments(
          FUN = oblimin_rotate, FUN.args = list(...)
        )
        
        # Check for arguments
        rotation_arguments$n.rotations <- ifelse(
          "n.rotations" %in% names(rot_arguments),
          rot_arguments$n.rotations,
          10
        )
        rotation_arguments$maxit <- ifelse(
          "maxit" %in% names(rot_arguments),
          rot_arguments$maxit,
          1000
        )
        
        # Set loadings
        rotation_arguments$loadings <- as.matrix(standardized)
        
        # Obtain rotated loadings
        rotated <- do.call(
          what = oblimin_rotate,
          args = as.list(rotation_arguments)
        )
        
      }
      
      # Re-align rotated loadings
      aligned_output <- fungible::faAlign(
        F1 = as.matrix(standardized),
        F2 = as.matrix(rotated$loadings),
        Phi2 = as.matrix(rotated$Phi)
      )
      
      # Update aligned loadings
      aligned_loadings <- aligned_output$F2
      colnames(aligned_loadings) <- colnames(standardized)
      row.names(aligned_loadings) <- row.names(standardized)
      
      # Update aligned correlations
      aligned_Phi <- aligned_output$Phi2
      colnames(aligned_Phi) <- colnames(standardized)
      row.names(aligned_Phi) <- colnames(standardized)
      
      # Re-assign rotated values
      rotated$loadings <- aligned_loadings
      rotated$Phi <- aligned_Phi
      
      
    }
  }
  
  # Set up results
  results <- list(
    unstd = loading_matrix,
    std = standardized,
    rotated = rotated,
    minLoad = min.load
  )
  
  # Set class
  class(results) <- "NetLoads"
  
  return(results)
  
  
}


# Bug checking ----
# 
# set.seed(1234)
# 
# # Generate data
# sim_data <- latentFactoR::simulate_factors(
#   factors = 3,
#   variables = 10,
#   loadings = 0.60,
#   cross_loadings = 0.10,
#   correlations = 0.30,
#   sample_size = 1000,
#   variable_categories = 5,
#   skew_range = c(-1, 1)
# )
# 
# # Add wording effects (for negative loadings)
# sim_data <- latentFactoR::add_wording_effects(
#   sim_data, method = "mixed"
# )
# 
# # Estimate EGA
# ega <- EGA(sim_data$data, plot.EGA = FALSE)
# ega$wc[8] <- NA
# A = ega; loading.method = "brm"
# rotation = "geominq"

#' @noRd
# Organize input ----
# Updated 10.07.2023
organize_input <- function(A, wc)
{
  
  # Check for `EGA` object
  if(any(class(A) %in% c("EGA", "EGA.fit", "riEGA"))){
    
    # Get `EGA` object
    ega_object <- get_EGA_object(A)
    
    # Set network and memberships
    A <- ega_object$network
    wc <- ega_object$wc
    
  }else{
    
    # Produce errors for miss aligned data
    length_error(wc, dim(A)[2]) # length between network and memberships
    object_error(A, c("matrix", "data.frame")) # must be matrix or data frame
    object_error(wc, c("vector", "matrix", "data.frame")) # must be one of these
    
  }
  
  # Generally, good to proceed
  A <- as.matrix(A); wc <- force_vector(wc)
  
  # Set memberships as string
  wc <- paste(wc)
  
  # Ensure names
  A <- ensure_dimension_names(A)
  names(wc) <- dimnames(A)[[2]]
  
  # Set orders
  ordering <- order(wc)
  
  # Return ordered network and memberships
  return(
    list(A = A[ordering, ordering], wc = wc[ordering])
  )
  
}

#' @noRd
# Obtain signs ----
# Function to obtain signs on dominant community
# Updated 11.07.2023
obtain_signs <- function(target_network)
{
  
  # Initialize signs to all positive orientation
  signs <- rep(1, dim(target_network)[2])
  names(signs) <- dimnames(target_network)[[2]]
  
  # Initialize row sums and minimum value
  row_sums <- rowSums(target_network, na.rm = TRUE)
  minimum_value <- which.min(row_sums)
  
  # Set while loop
  while(sign(row_sums[minimum_value]) == -1){
    
    # Get minimum value name
    minimum_name <- names(minimum_value)
    
    # Flip variable
    target_network[minimum_name,] <- 
      target_network[,minimum_name] <-
      -target_network[,minimum_name]
    
    # Set sign as flipped
    signs[minimum_name] <- -signs[minimum_name]
    
    # Update row sums and minimum value
    row_sums <- rowSums(target_network, na.rm = TRUE)
    minimum_value <- which.min(row_sums)
    
  }
  
  # Add signs as an attribute to the target network
  attr(target_network, "signs") <- signs
  
  # Return results
  return(target_network)
  
}

#' @noRd
# Descending order ----
# Updated 11.07.2023
descending_order <- function(standardized, wc, unique_communities) 
{
  
  # Initialize order names
  order_names <- character(dim(standardized)[1])
  
  # Loop over communities
  for(community in unique_communities){
    
    # Get community index
    community_index <- wc == community
    
    # Get order
    ordering <- order(standardized[community_index, community])
    
    # Input ordering into order names
    order_names[community_index] <- dimnames(standardized)[[1]][community_index]
    
  }
  
  
  
  # Return reordered results
  return(standardized[order_names,])
  
}

#%%%%%%%%%%%%%%%%%
# BRM Legacy ----
#%%%%%%%%%%%%%%%%%

#' @noRd
## Absolute weights ("BRM") ----
# Updated 10.07.2023
absolute_weights <- function(A, wc, nodes, unique_communities)
{
  
  # Ensure network is absolute
  A <- abs(A)
  
  # Loop over communities
  return(
    nvapply(
      unique_communities, function(community){
        colSums(A[wc == community,, drop = FALSE], na.rm = TRUE)
      }, LENGTH = nodes
    )
  )
  
}

#' @noRd
## Add signs ("BRM") ----
# From CRAN version 1.2.3
# Updated 10.07.2023
old_add_signs <- function(unstandardized, A, wc, unique_communities)
{
  
  # Loop over main loadings
  for(community in unique_communities){
    
    # Get community index
    community_index <- wc == community
    
    # Get number of nodes
    node_count <- sum(community_index)
    
    # Get community sub-network
    community_network <- A[community_index, community_index, drop = FALSE]
    
    # Initialize sign matrix
    community_signs <- sign(community_network)
    
    # Initialize signs to all positive
    signs <- rep(1, node_count)
    
    # Loop over nodes
    for(node in seq_len(node_count)){
      
      # Make copy of signs
      signs_copy <- community_signs
      
      # Get current maximum sum
      current_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
      
      # Flip sign of each node
      community_signs[node,] <- -community_signs[node,]
      
      # Get new maximum sum
      new_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
      
      # Check for increase
      if(new_max > current_max){
        signs[node] <- -1 # with increase, flip sign
      }else{ # otherwise, return sign matrix to original state
        community_signs <- signs_copy
      }
      
    }
    
    # Update signs in loadings
    unstandardized[community_index, community] <-
      unstandardized[community_index, community] * signs
    
    # Sweep across community
    A[, community_index] <- sweep(
      A[, community_index, drop = FALSE], MARGIN = 2, signs, `*`
    )
    
  }
  
  # Loop over communities
  for(community1 in unique_communities){
    
    # Get first community index
    community_index1 <- wc == community1
    
    # Get number of nodes
    node_count <- sum(community_index1)
    
    # Loop over other communities
    for(community2 in unique_communities){
      
      # Check for the same community
      if(community1 != community2){
        
        # Get second community index
        community_index2 <- wc == community2
        
        # Get community sub-network
        community_network <- A[community_index1, community_index2, drop = FALSE]
        
        # Initialize sign matrix
        community_signs <- sign(community_network)
        
        # Initialize signs to all positive
        signs <- rep(1, node_count)
        
        # Loop over nodes
        for(node in seq_len(node_count)){
          
          # Make copy of signs
          signs_copy <- community_signs
          
          # Get current maximum sum
          current_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
          
          # Flip sign of each node
          community_signs[node,] <- -community_signs[node,]
          
          # Get new maximum sum
          new_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
          
          # Check for increase
          if(new_max > current_max){
            signs[node] <- -1 # with increase, flip sign
          }else{ # otherwise, return sign matrix to original state
            community_signs <- signs_copy
          }
          
        }
        
        # Update signs in loadings
        unstandardized[community_index1, community2] <-
          unstandardized[community_index1, community2] * signs
        
      }
      
    }
    
  }
  
  # Flip direction of community with main loadings
  for(community in unique_communities){
    
    # Get community indices
    community_index <- wc == community
    
    # Determine direction with sign
    if(sign(sum(unstandardized[community_index, community])) != 1){
      unstandardized[community_index,] <- 
        -unstandardized[community_index,]
    }
    
  }
  
  # Return unstandardized loadings
  return(unstandardized)
  
}

