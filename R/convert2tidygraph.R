#' Convert networks to \code{tidygraph}
#'
#' @description Converts networks to \code{tidygraph} format
#'
#' @param EGA.object
#' A single \code{\link{EGAnet}} object containing the outputs
#' \code{$network} and \code{$wc}
#'
#' @examples
#' convert2tidygraph(ega.wmt)
#' convert2tidygraph(boot.wmt)
#'
#' @return Returns a network in the \code{tidygraph} format
#'
#' @author Dominique Makowski, Hudson Golino <hfg9s at virginia.edu>, & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Convert network to {tidygraph}
# Updated 14.07.2023
convert2tidygraph <- function(EGA.object)
{

  if("bootEGA" %in% class(EGA.object)) {
    # Bootstrapped EGA ---
    if(!"typicalGraph" %in% names(EGA.object)) {
      stop("The bootEGA object must contain a typicalGraph object. Set `typicalStructure=TRUE`.")
    }

    if("lower_order" %in% names(EGA.object$typicalGraph)) {
      # Hierarchical EGA ---
      lower_nodes <- .convert2tidygraph_nodes(EGA.object$typicalGraph$lower_order$wc)
      higher_nodes <- .convert2tidygraph_nodes(EGA.object$typicalGraph$higher_order$wc)
      higher_nodes$dimension <- paste0("H", higher_nodes$dimension)  # Discriminate from lower
      nodes <- rbind(lower_nodes, higher_nodes)

      lower_edges <- .convert2tidygraph_edges(EGA.object$typicalGraph$lower_order$graph)
      higher_edges <- .convert2tidygraph_edges(EGA.object$typicalGraph$higher_order$graph)
      edges <- rbind(lower_edges, higher_edges)
      edges$type <- "real"

      # Make edges from lower to higher order
      for(i in 1:nrow(lower_nodes)) {
        edges <- rbind(
          edges,
          data.frame(from = lower_nodes$name[i], to = lower_nodes$dimension[i], link = 1, type = "virtual")
          )
      }
    } else {
      # Non-hierarchical EGA ---
      nodes <- .convert2tidygraph_nodes(EGA.object$typicalGraph$wc)
      edges <- .convert2tidygraph_edges(EGA.object$typicalGraph$graph)
    }

  } else {
    # Normal EGA ---
    nodes <- .convert2tidygraph_nodes(EGA.object$wc)
    edges <- .convert2tidygraph_edges(EGA.object$network)
  }

  # Return result
  return(
    list(
      nodes = nodes,
      edges = edges[edges$link != 0,]
    )
  )
}


#' @keywords internal
#' @noRd
.convert2tidygraph_nodes <- function(wc){
  data.frame(
    name = names(wc),
    dimension = as.character(wc)
  )
}

#' @keywords internal
#' @noRd
.convert2tidygraph_edges <- function(network){
  data.frame(
    from = rep(rownames(network), times = ncol(network)),
    to = rep(colnames(network), each = nrow(network)),
    link = as.vector(network)
  )
}
