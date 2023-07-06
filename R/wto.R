#' Weighted Topological Overlap
#' 
#' @description Computes weighted topological overlap following
#' the Novick et al. (2009) definition
#'
#' @param network Symmetric matrix.
#' Input for a symmetric network matrix
#' 
#' @param signed Boolean (length = 1).
#' Whether the signed version should be used.
#' Defaults to \code{TRUE}.
#' Use \code{FALSE} for absolute values
#' 
#' @param diagonal.zero Boolean (length = 1).
#' Whether diagonal of overlap matrix should be set to zero.
#' Defaults to \code{TRUE}.
#' Use \code{FALSE} to allow overlap of a node with itself
#' 
#' @examples 
#' # Obtain network
#' network <- network.estimation(wmt2[,7:24], model = "glasso")
#' 
#' # Compute wTO
#' wto(network)
#' 
#' @references 
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' 
#' @return A symmetric matrix of weighted topological overlap
#' values between each pair of variables
#' 
#' @export
#' 
# Weighted Topological Overlap ----
# Updated 02.07.2023
wto <- function (network, signed = TRUE, diagonal.zero = TRUE)
{
  
  # Remove {EGAnet} bulk
  network <- remove_attributes(network)
  
  # Ensure network is matrix
  network <- as.matrix(network)
  
  # Get dimensions of the network
  dimensions <- dim(network)
  
  # Obtain absolute network values
  absolute_network <- abs(network)
  
  # Determine whether absolute values should constitute the network
  if(isFALSE(signed)){
    network <- absolute_network
  }
  
  # Obtain node strengths
  node_strengths <- strength(absolute_network)
  
  # Obtain variable pair minimums
  minimum_pairs <- matrix(
    c(
      rep(node_strengths, each = dimensions[2]),
      rep(node_strengths, times = dimensions[2])
    ), ncol = 2
  )
  
  # Create minimum matrix
  minimum_matrix <- matrix(
    ifelse(
      minimum_pairs[,1] < minimum_pairs[,2],
      minimum_pairs[,1], minimum_pairs[,2]
    ), 
    nrow = dimensions[2], ncol = dimensions[2]
  )
  
  # Divide numerator by denominator
  omega <- (crossprod(network) + network) / 
           (minimum_matrix + 1 - absolute_network)
  
  # Set diagonal to zero
  if(isTRUE(diagonal.zero)){
    diag(omega) <- 0
  }
  
  # Return weighted topological overlap
  return(omega)
  
}

# Bug Checking ----
## Basic input
# network = network.estimation(wmt2[,7:24], model = "glasso")
# signed = TRUE; diagonal.zero = TRUE
