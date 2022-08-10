#' Convert networks to \code{\link{igraph}}
#'
#' @description Converts networks to \code{\link{igraph}} format
#'
#' @param A Matrix or data frame.
#' \emph{N} x \emph{N} matrix where \emph{N} is the number of nodes
#' 
#' @param diagonal Numeric.
#' Value to be placed on the diagonal of \code{A}.
#' Defaults to \code{0}
#' 
#' @examples
#' convert2igraph(ega.wmt$network)
#'
#' @return Returns a network in the \code{\link{igraph}} format
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Convert network to {igraph}
# Updated 10.08.2022
convert2igraph <- function (A, diagonal = 0)
{
  # Make diagonal zero
  diag(A) <- diagonal
  
  return(
    suppressWarnings(
      igraph::graph_from_adjacency_matrix(
        as.matrix(A), weighted = TRUE, mode = "undirected",
        add.colnames = FALSE
      )
    )
  )
}
