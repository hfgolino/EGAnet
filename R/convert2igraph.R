#' @title Convert networks to \code{igraph}
#'
#' @description Converts networks to \code{igraph} format
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
#' @return Returns a network in the \code{igraph} format
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Convert matrix to {igraph} network
# Updated 09.08.2023
convert2igraph <- function (A, diagonal = 0)
{

  # Argument errors (return A in case of tibble)
  A <- convert2igraph_errors(A, diagonal)

  # Convert to matrix
  A <- as.matrix(A)

  # Change diagonal (to zero)
  diag(A) <- diagonal

  # Return {igraph} network
  return(
    silent_call(
      igraph::graph_from_adjacency_matrix(
        A, weighted = TRUE, mode = "undirected",
        add.colnames = FALSE
      )
    )
  )

}

#' @noRd
# Argument errors
# Updated 13.08.2023
convert2igraph_errors <- function(A, diagonal)
{

  # 'A' errors
  object_error(A, c("matrix", "data.frame", "tibble"), "convert2igraph")

  # Check for tibble
  if(get_object_type(A) == "tibble"){
    A <- as.data.frame(A)
  }

  # 'diagonal' errors
  length_error(diagonal, 1, "convert2igraph")
  typeof_error(diagonal, "numeric", "convert2igraph")
  range_error(diagonal, c(-1, 1), "convert2igraph")

  # Return A in case of tibble
  return(A)

}


