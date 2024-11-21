#' @title Weighted Topological Overlap
#' 
#' @description Computes weighted topological overlap following
#' the Novick et al. (2009) definition
#'
#' @param network Symmetric matrix or data frame.
#' A symmetric network
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
#' \strong{Original formalization} \cr
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
# About 10x faster than `wTO::wTO`
# Updated 07.08.2023
wto <- function (network, signed = TRUE, diagonal.zero = TRUE)
{
  
  # Check for errors, remove attributes, and ensure network is matrix
  network <- wto_errors(network, signed, diagonal.zero)
  
  # Get dimensions of the network
  dimensions <- dim(network)
  
  # Obtain absolute network values
  absolute_network <- abs(network)
  
  # Determine whether absolute values should constitute the network
  if(!signed){
    network <- absolute_network
  }
  
  # Obtain node strengths
  node_strengths <- colSums(absolute_network, na.rm = TRUE)
  
  # Obtain variable pair minimums
  strength_each <- rep(node_strengths, each = dimensions[2])
  strength_times <- rep(node_strengths, times = dimensions[2])
  
  # Create minimum matrix
  minimum_matrix <- matrix(
    swiftelse(
      strength_each < strength_times,
      strength_each, strength_times
    ), 
    nrow = dimensions[2], ncol = dimensions[2]
  )
  
  # Divide numerator by denominator
  omega <- (crossprod(network) + network) / 
           (minimum_matrix + 1 - absolute_network)
  
  # Set diagonal to zero
  if(diagonal.zero){
    diag(omega) <- 0
  }
  
  # Return weighted topological overlap
  return(omega)
  
}

#' @noRd
# Argument errors ----
# Updated 13.08.2023
wto_errors <- function(network, signed, diagonal.zero)
{
  
  # 'network' errors
  object_error(network, c("matrix", "data.frame", "tibble"), "wto")
  
  # 'signed' errors
  length_error(signed, 1, "wto")
  typeof_error(signed, "logical", "wto")
  
  # 'diagonal.zero' errors
  length_error(diagonal.zero, 1, "wto")
  typeof_error(diagonal.zero, "logical", "wto")
  
  # Return network without attributes and as matrix
  return(as.matrix(remove_attributes(network)))
  
}

# Bug Checking ----
## Basic input
# network = network.estimation(wmt2[,7:24], model = "glasso")
# signed = TRUE; diagonal.zero = TRUE
