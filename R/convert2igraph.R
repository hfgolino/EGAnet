#' Convert networks to \code{\link{igraph}}
#'
#' @description Converts networks to \code{\link{igraph}} format
#'
#' @param A Matrix or data frame.
#' \emph{N} x \emph{N} matrix where \emph{N} is the number of nodes
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
# Updated 14.07.2022
convert2igraph <- function (A)
{
  return(
    suppressWarnings(
      igraph::graph_from_adjacency_matrix(
        as.matrix(A), weighted = TRUE, mode = "undirected",
        add.colnames = FALSE
      )
    )
  )
}
