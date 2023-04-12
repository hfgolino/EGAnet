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
#' @param min.load Numeric.
#' Sets the minimum loading allowed in the standardized
#' network loading matrix. Values equal or greater than
#' the minimum loading are kept in the output. Values
#' less than the minimum loading are removed. This matrix can
#' be viewed using \code{print()} or \code{summary()}
#' Defaults to \code{0}
#' 
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure. For a list of rotations,
#' see \link{GPArotation}
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
# Updated 12.04.2023
# Cross-loadings and signs updated 12.04.2023
# Rotations added 20.10.2022
net.loads <- function(
    A, wc, rotation = "oblimin",
    min.load = 0,
    ...
)
{
  
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
  
  # Ensure names in memberships
  names(wc) <- colnames(A)
  
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
    
  }else if(length(unique_wc) == 1){# Check for single community
    
    # Reorder communities
    wc <- wc[wc_order]
    
    # Reorder network
    A <- A[wc_order, wc_order]
    
    # Compute absolute sum of nodes
    absolute_sum <- colSums(abs(A))
    
    # Come back to this section
    stop("unidimensional structures are under construction...")
    
  }else{ # Not singleton dimensions, not unidimensional
    
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
    signs <- numeric(ncol(A))
    names(signs) <- colnames(A)
    
    # Add column and row names
    row.names(loading_matrix) <- colnames(A)
    colnames(loading_matrix) <- unique_wc
    
    # Populate loading matrix
    for(dominant in unique_wc){
      
      # Obtain target portion of network
      target_network <- A[wc == dominant, wc == dominant]
      
      # Compute absolute sum for dominant loadings
      loading_matrix[wc == dominant, as.character(dominant)] <- 
        colSums(abs(target_network))
      
      # Determine positive direction for dominant loadings
      ## Compute total sum of signs
      target_signs <- colSums(sign(target_network))
      
      ## Check for negative signs
      if(all(target_signs <= -1)){
        
        ## Determine dominant sign based on other variables
        signs[wc == dominant] <- sign(colSums(A[wc != dominant, wc == dominant]))
        
        
      }else if(sum(target_signs) >= 0){
        
        ## If target signs are equal to or greater than break even,
        ## then reverse the negative signs
        signs[wc == dominant][which(target_signs < 0)] <- -1 
        signs[wc == dominant][which(target_signs >= 0)] <- 1
        
      }else if(sum(target_signs) < 0){
        
        ## If target signs are less than break even,
        ## then reverse the positive signs
        signs[wc == dominant][which(target_signs < 0)] <- 1 
        signs[wc == dominant][which(target_signs >= 0)] <- -1
        
      }
      
    }
    
    # Set dominant loadings
    loading_matrix <- loading_matrix * signs
    
    # Initialize reversed A
    A_reversed <- A
    
    # Create duplicate of network
    if(any(signs == -1)){
      A_reversed[which(signs == -1),] <- -A[which(signs == -1),]
    }
    
    # Populate loading matrix
    for(dominant in unique_wc){
      for(cross in unique_wc){
        
        # Do not use dominant loadings
        if(dominant != cross){
          
          # Compute algebraic sum for cross-loadings
          loading_matrix[wc == dominant, as.character(cross)] <- 
            colSums(A_reversed[wc == cross, wc == dominant])
          
        }
        
      }
    }
    
    # Using signs, ensure positive orientation based
    # on most common direction
    for(dominant in unique_wc){
      
      # Determine dominant orientation
      orientation <- sum(signs[wc == dominant])
      
      # Check for negative orientation
      if(orientation <= -1){
        
        loading_matrix[wc == dominant,] <-
          -loading_matrix[wc == dominant,]
        
      }
      
    }
    
    # Obtain standardized loadings
    standardized <- t(
      t(loading_matrix) /
      sqrt(colSums(abs(loading_matrix)))
    )
    
    # Apply descending order
    standardized <- descend.ord(standardized, wc)
    loading_matrix <- loading_matrix[
      row.names(standardized), colnames(standardized)
    ] # orients matrix into same ordering
    
    # Set up for rotation
    
    # Check for {GPArotation} and {fungible}
    # Function in `helpers-general.R`
    check_package(c("GPArotation", "fungible"))
    
    # Obtain rotation from GPArotation package
    rotation_names <- ls(asNamespace("GPArotation"))
    
    # Check if rotation exists
    rotation <- tolower(rotation)
    rotation_names_lower <- tolower(rotation_names)
    
    # Obtain rotation arguments
    rot_arguments <- list(...)
    
    # Check if rotation exists
    if(rotation %in% rotation_names_lower){
      
      if(rotation != "oblimin"){
        
        # Obtain arguments
        rotation_arguments <- obtain.arguments(
          FUN = psych::faRotations, FUN.args = list(rotate = rotation)
        )
        
        # Check for arguments
        rotation_arguments$loadings <- std
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
        F2 = as.matrix(rotated$loadings)
      )
      
      # Update aligned loadings
      aligned_loadings <- aligned_output$F2
      colnames(aligned_loadings) <- colnames(standardized)
      
      # Rename Phi
      colnames(rotated$Phi) <- colnames(standardized)
      row.names(rotated$Phi) <- colnames(standardized)
      
      # Re-assign rotated loadings
      rotated$loadings <- aligned_loadings
    
    
    }
  }
  
  # Set up results
  results <- list(
    unstd = loading_matrix,
    std = standardized,
    rotated = rotated,
    minLoad = min.load
  )
  
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
# A = ega; rotation = "oblimin";
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







