#' Convert \code{\link{igraph}} network to matrix
#'
#' @description Converts \code{\link{igraph}} network to matrix
#'
#' @param igraph_network \code{\link{igraph}} network object
#' 
#' @param diagonal Numeric.
#' Value to be placed on the diagonal of \code{A}.
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
# Convert {igraph} network to matrix
# Updated 11.05.2023
igraph2matrix <- function (igraph_network, diagonal = 0)
{
  
  # Convert {igraph} network to matrix
  network <- as.matrix(
    igraph::as_adjacency_matrix(
      graph = igraph_network,
      type = "both",
      attr = "weight"
    )
  )
  
  # Add back names
  colnames(network) <- igraph::vertex.attributes(igraph_network)$`FALSE`
  row.names(network) <- colnames(network)
  
  # Make diagonal zero
  diag(network) <- diagonal
  
  # Return network
  return(network)
  
}
