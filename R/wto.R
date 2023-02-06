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
#' @param diagonal_zero Boolean (length = 1).
#' Whether diagonal of overlap matrix should be set to zero.
#' Defaults to \code{TRUE}.
#' Use \code{FALSE} to allow overlap of a node with itself
#' 
#' @examples 
#' # Obtain network
#' network <- EBICglasso.qgraph(wmt2[,7:24])
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
# Updated 03.02.2023
wto <- function (network, signed = TRUE, diagonal_zero = TRUE)
{
  
  # Ensure network is matrix
  network <- as.matrix(network)
  
  # Obtain absolute network values
  absolute_network <- abs(network)
  
  # Determine whether absolute values
  # should constitute the network
  if(!isTRUE(signed)){
    network <- absolute_network
  }
  
  # Obtain numerator
  numerator <- (network %*% t(network)) + # sum of connections
    network # connections between edge node
  
  # Obtain node strengths
  node_strengths <- colSums(absolute_network, na.rm = TRUE)
  
  # Obtain variable pair minimums
  minimum_df <- data.frame(
    node_i = rep(
      node_strengths,
      each = length(node_strengths)
    ),
    node_j = rep(
      node_strengths,
      times = length(node_strengths)
    )
  )
  
  # Obtain minimums for each pair
  minimum_vector <- apply(
    as.matrix(minimum_df), 1,
    min, na.rm = TRUE
  )
  
  # Create matrix
  minimum_matrix <- matrix(
    minimum_vector,
    nrow = nrow(network),
    ncol = ncol(network)
  )
  
  # Obtain denominator
  denominator <- minimum_matrix + 1 - abs(network)
  
  # Divide numerator by denominator
  omega <- numerator / denominator
  
  # Set diagonal to zero
  if(isTRUE(diagonal_zero)){
    diag(omega) <- 0
  }
  
  # Return weighted topological overlap
  return(omega)
  
}
