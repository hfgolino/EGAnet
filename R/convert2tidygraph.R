#' Convert networks to \code{\link{tidygraph}}
#'
#' @description Converts networks to \code{\link{tidygraph}} format
#'
#' @param EGA.object
#' A single \code{\link{EGAnet}} object containing the outputs
#' \code{$network} and \code{$wc}
#' 
#' @examples
#' convert2tidygraph(ega.wmt)
#'
#' @return Returns a network in the \code{\link{tidygraph}} format
#' 
#' @author Dominique Makowski, Hudson Golino <hfg9s at virginia.edu>, & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Convert network to {tidygraph}
# Updated 14.07.2023
convert2tidygraph <- function (EGA.object)
{
  
  # Set up nodes
  nodes <- data.frame(
    name = names(EGA.object$wc),
    dimension = as.character(EGA.object$wc)
  )
  
  # Get number of nodes
  node_count <- length(EGA.object$wc)
  
  # Get node sequence
  node_sequence <- seq_len(node_count)
  
  # Set up edges
  edges <- data.frame(
    from = rep(node_sequence, each = node_count),
    to = rep(node_sequence, times = node_count),
    link = abs(force_vector(EGA.object$network))
  )
  
  # Return result
  return(
    list(
      nodes = nodes,
      edges = edges[edges$link != 0,]
    )
  )
  
}
