#' Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (TEFI) of a dimensionality structure using Von Neumman's entropy when the input is a correlation matrix.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data A matrix, data frame, or correlation matrix
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{VN.Entropy.Fit}{The Entropy Fit Index using Von Neumman's entropy}
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Estimate EGA model
#' ega.wmt <- EGA(
#'   data = wmt, model = "glasso",
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#'
#' # Compute entropy indices
#' tefi(data = ega.wmt$correlation, structure = ega.wmt$wc)
#'
#' @references
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com>, and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
# Updated 26.06.2023
tefi <- function(data, structure)
{
  
  # Generic function to get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = 1, # set to 1 to avoid error
    corr = "auto", na.data = "pairwise", 
    verbose = FALSE
  )
  
  # Get absolute correlations
  data <- abs(output$correlation_matrix)
  
  # Check structure
  if(anyNA(structure)){
    
    # Determine variables that are NA
    rm.vars <- is.na(structure)

    # Send warning message
    warning(paste("Some variables did not belong to a dimension:", colnames(data)[rm.vars]))
    message("Use caution: These variables have been removed from the TEFI calculation")

    # Keep available variables
    data <- data[!rm.vars, !rm.vars]
    
  }

  # Obtain density matrix
  density_matrix <- data / dim(data)[2]
  
  # Obtain Von Neumann's entropy
  h.vn <- -trace(density_matrix %*% log(density_matrix))
  
  # Obtain communities
  communities <- unique_length(structure)
  
  # Initialize Von Neumman entropy by dimension
  h.vn.wc <- nnapply(seq_len(communities), function(community){
    
    # Get indices
    indices <- structure == community
    
    # Obtain community matrix
    community_matrix <- data[indices, indices]
    
    # Obtain community density matrix
    community_density <- community_matrix / dim(community_matrix)[2]
    
    # Return Von Neumann entropy
    return(-trace(community_density %*% log(community_density)))
    
  })
  
  # Compute values
  ## Sum of community Von Neumann
  sum.h.vn.wc <- sum(h.vn.wc, na.rm = TRUE)
  ## Mean of community Von Neumann
  mean.h.vn.wc <- sum.h.vn.wc / communities
  ## Difference between total and total community
  Hdiff <- h.vn - sum.h.vn.wc
  ## Average entropy
  mean.vn <- mean.h.vn.wc - h.vn
  
  # Set up results
  results <- data.frame(
    "VN.Entropy.Fit" = mean.vn + (Hdiff * sqrt(communities)),
    "Total.Correlation" = sum.h.vn.wc - h.vn,
    "Average.Entropy" = mean.vn
  )
  
  # Return results
  return(results)
  
}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]; ega.wmt <- EGA(data, plot.EGA = FALSE)
# data <- ega.wmt$correlation
# structure <- ega.wmt$wc