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
# Updated 22.04.2023
# Cross-loadings and signs updated 12.04.2023
# Rotations added 20.10.2022
net.loads <- function(
    A, wc, sum.method = c("absolute", "signed"),
    positive.orientation = FALSE,
    rotation = "geominQ",
    min.load = 0, ...
)
{
  
  # Check for sum method
  if(missing(sum.method)){
    sum.method <- "signed"
  }
  
  # Select sum function
  if(sum.method == "absolute"){
    sum_function <- function(x){colSums(abs(x))}
  }else if(sum.method == "signed"){
    sum_function <- colSums
  }
  
  # Check for EGA object
  if(is(A, "EGA")){
    
    # Order
    wc_order <- match(colnames(A$network), names(A$wc))
    
    # Grab communities
    wc <- A$wc
    
    # Replace 'A' with 'EGA' network
    A <- A$network
    
  }else{
    
    # Obtain membership order
    wc_order <- order(wc)
    
  }
  
  # Ensure names
  names(wc) <- colnames(A)
  row.names(A) <- colnames(A)
  
  # Make "A" a matrix
  A <- as.matrix(A)
  
  # Check for symmtric
  if(!is_symmetric(A)){ # Function in `helpers-general.R`
    stop("Input for 'A' must be a symmetric n x n matrix.")
  }
  
  # Obtain NA community memberships
  NA_wc <- is.na(wc)
  
  # Check for any NA community memberships
  if(any(NA_wc)){
    
    # Print warning message
    ## Let's user know that these nodes have been removed
    warning(
      paste(
        "The following nodes were found to have NA community membership:",
        paste0(colnames(A)[NA_wc], collapse = ", "), "\n",
        "The nodes have been removed from the loadings matrix."
      )
    )
    
    # Remove from network
    A <- A[!NA_wc, !NA_wc]
    
    # Remove from memberships
    wc <- wc[!NA_wc]
    
  }
  
  # Obtain number of unique communities
  unique_wc <- unique(wc)
  
  # Obtain length of communities
  length_wc <- length(wc)
  
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

# set.seed(1234)
# 
# # Generate data
# sim_data <- latentFactoR::simulate_factors(
#   factors = 3,
#   variables = 3,
#   loadings = 0.60,
#   cross_loadings = 0.10,
#   correlations = 0.30,
#   sample_size = 1000,
#   variable_categories = 5,
#   skew_range = c(-1, 1)
# )
# 
# # Estimate EGA
# ega <- EGA(sim_data$data)
# A = ega; rotation = "geominQ";
# min.load = 0; rot_arguments = list();
# source("./utils-EGAnet.R")
# source("./helpers-general.R")
# source("./helpers-functions.R")
# source("./helpers-errors.R")

# Descending order ----
#' @noRd
# Function to order loadings largest to smallest
# within their respective factors
descend.ord <- function(loads, wc){
  
  # Initialize ordering vector
  ord.names <- vector("character")
  
  # Loop through dimensions
  for(i in colnames(loads)){
    ord <- order(loads[names(which(wc == i)),i], decreasing = TRUE)
    ord.names <- c(ord.names, names(which(wc == i))[ord])
  }
  
  # Reorder
  reord <- loads[ord.names,]
  
  # Check for matrix
  if(!is.matrix(reord)){
    reord <- as.matrix(reord)
  }
  
  # Make sure names
  row.names(reord) <- ord.names
  colnames(reord) <- colnames(loads)
  
  return(reord)
  
}

# Obtain signs ----
#' @noRd
# Function to obtain signs on dominant community
obtain_signs <- function(target_network)
{
  
  # Initialize signs
  signs <- rep(1, ncol(target_network)) # start with all positive orientation
  names(signs) <- colnames(target_network)
  
  # Set minimum 
  row_sums <- -1
  minimum_value <- 1
  
  # Set while loop
  while(sign(row_sums[minimum_value]) == -1){
    
    # Sum of rows
    row_sums <- rowSums(target_network, na.rm = TRUE)
    
    # Find minimum value
    minimum_value <- which.min(row_sums)
    
    # Check for negative
    if(sign(row_sums[minimum_value]) == -1){
      
      # Flip variable
      target_network[names(minimum_value),] <- -target_network[names(minimum_value),]
      target_network[,names(minimum_value)] <- -target_network[,names(minimum_value)]
      
      # Set sign as flipped
      signs[names(minimum_value)] <- -signs[names(minimum_value)]
      
    }
    
  }
  
  # Set up results
  results <- list(
    target_network = target_network,
    signs = signs
  )
  
  # Return results
  return(results)
  
}

# Obtain signs ----
#' @noRd
# Function to obtain signs on dominant community
old.add.signs <- function(comm.str, A, wc, dims, pos.manifold)
{
  
  # Set NA to "NA"
  if(any(is.na(wc))){
    wc <- ifelse(is.na(wc), "NA", wc)
  }
  
  # Loop through self
  for(i in dims){
    
    # Set minimum 
    row_sums <- -1
    minimum_value <- 1
    
    # Set while loop
    while(sign(row_sums[minimum_value]) == -1){
      
      # Sum of rows
      row_sums <- rowSums(A[wc == i, wc == i], na.rm = TRUE)
      
      # Find minimum value
      minimum_value <- which.min(row_sums)
      
      # Check for negative
      if(sign(row_sums[minimum_value]) == -1){
        
        # Flip variable
        A[names(minimum_value), wc == i] <- 
          -A[names(minimum_value), wc == i]
        A[wc == i, names(minimum_value)] <- 
          -A[wc == i, names(minimum_value)]
        
        # Add negative
        comm.str[names(minimum_value), as.character(i)] <- 
          -comm.str[names(minimum_value), as.character(i)]
        
      }
      
    }
    
  }
  
  # Check for unidimensional structure
  if(ncol(comm.str) > 1){
    
    # Set combinations
    combinations <- combn(
      dims, m = 2
    )
    
    # Loop through combinations
    for(i in 1:ncol(combinations)){
      
      # Set targets
      target1 <- combinations[1,i]
      target2 <- combinations[2,i]
      
      # Set minimum 
      row_sums <- -1
      minimum_value <- 1
      
      # Set while loop
      while(sign(row_sums[minimum_value]) == -1){
        
        # Sum of rows
        row_sums <- rowSums(A[wc == target1, wc == target2], na.rm = TRUE)
        
        # Find minimum value
        minimum_value <- which.min(row_sums)
        
        # Check for negative
        if(sign(row_sums[minimum_value]) == -1){
          
          # Flip variable
          A[names(minimum_value), wc == target2] <- 
            -A[names(minimum_value), wc == target2]
          
          # Add negative
          comm.str[names(minimum_value), as.character(target2)] <- 
            -comm.str[names(minimum_value), as.character(target2)]
          
        }
        
      }
      
    }
    
    # Loop through combinations (switches `target1` with `target2`)
    for(i in 1:ncol(combinations)){
      
      # Set targets
      target1 <- combinations[2,i]
      target2 <- combinations[1,i]
      
      # Set minimum 
      row_sums <- -1
      minimum_value <- 1
      
      # Set while loop
      while(sign(row_sums[minimum_value]) == -1){
        
        # Sum of rows
        row_sums <- rowSums(A[wc == target1, wc == target2], na.rm = TRUE)
        
        # Find minimum value
        minimum_value <- which.min(row_sums)
        
        # Check for negative
        if(sign(row_sums[minimum_value]) == -1){
          
          # Flip variable
          A[names(minimum_value), wc == target2] <- 
            -A[names(minimum_value), wc == target2]
          
          # Add negative
          comm.str[names(minimum_value), as.character(target2)] <- 
            -comm.str[names(minimum_value), as.character(target2)]
          
        }
        
      }
      
    }
    
  }
  
  # Flip dimensions (if necessary)
  # if(!pos.manifold)
  # {
  #   for(i in 1:length(dims))
  #   {
  #     wc.sign <- sign(sum(comm.str[which(wc==dims[i]),i]))
  # 
  #     if(wc.sign != 1)
  #     {comm.str[which(wc==dims[i]),] <- -comm.str[which(wc==dims[i]),]}
  #   }
  # }
  
  # res <- list()
  # res$comm.str <- comm.str
  # res$A <- A
  
  return(comm.str)
}



