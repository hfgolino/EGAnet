#' @title Convert \code{\link{igraph}} network to matrix
#'
#' @description Converts \code{\link{igraph}} network to matrix
#'
#' @param igraph_network \code{\link{igraph}} network object
#' 
#' @param diagonal Numeric (length = 1).
#' Value to be placed on the diagonal of \code{network}.
#' Defaults to \code{0}
#' 
#' @examples
#' # Convert network to {igraph}
#' igraph_network <- convert2igraph(ega.wmt$network)
#' 
#' # Convert network back to matrix
#' igraph2matrix(igraph_network)
#'
#' @return Returns a network in the \code{\link{igraph}} format
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
#' 
# Convert {igraph} network to matrix
# Updated 03.08.2023
igraph2matrix <- function (igraph_network, diagonal = 0)
{
  
  # Argument errors
  igraph2matrix_errors(igraph_network, diagonal)
  
  # Convert {igraph} network to matrix
  network <- as.matrix(
    igraph::as_adjacency_matrix(
      graph = igraph_network,
      type = "both", attr = "weight"
    )
  )
  
  # Get node names
  node_names <- igraph::vertex.attributes(igraph_network)$`FALSE`
  
  # Add back names
  dimnames(network) <- list(node_names, node_names)
  
  # Make diagonal zero
  diag(network) <- diagonal
  
  # Return network
  return(network)
  
}

#' @noRd
# Argument errors
# Updated 13.08.2023
igraph2matrix_errors <- function(igraph_network, diagonal)
{
  
  # 'igraph_network' errors
  class_error(igraph_network, "igraph", "igraph2matrix")
  
  # 'diagonal' errors
  length_error(diagonal, 1, "igraph2matrix")
  typeof_error(diagonal, "numeric", "igraph2matrix")
  range_error(diagonal, c(-1, 1), "igraph2matrix")
  
}
