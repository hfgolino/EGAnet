#' Convert EGAnet objects to tidygraph
#'
#' Converts an EGA object to the input required by \code{tidygraph::tbl_graph()}
#' so that it can be plotted using \code{ggraph}.
#'
#' @param ega An EGA object.
#'
#' @examples
#' \dontrun{
#'   library(tidygraph)
#'   library(ggraph)
#'
#'   ega.wmt <- EGA(wmt2[,7:24], plot.EGA = FALSE)
#'
#'   x <- ega_to_tidygraph(ega.wmt)
#'
#'
#'   graph <- tidygraph::tbl_graph(nodes = x$nodes, edges = x$edges)
#'
#'   ggraph::ggraph(graph) +
#'    geom_edge_link(aes(colour = link)) +
#'    geom_node_point(aes(colour = dimension))
#' }
#' @export
ega_to_tidygraph <- function(ega) {
  # Nodes
  nodes <- ega$dim.variables
  row.names(nodes) <- NULL
  names(nodes)[names(nodes) == "items"] <- c("name")
  nodes$dimension <- as.character(nodes$dimension)

  # Edges
  edges <- expand.grid(data.frame(from = as.numeric(as.factor(nodes$name)),
                                  to = as.numeric(as.factor(nodes$name))))
  edges$link <- NA
  for(row in 1:nrow(edges)) {
    subset <- edges[row, ]
    edges[row, "link"] <- ega$network[nodes$name[subset$from], nodes$name[subset$to]]
  }
  edges <- edges[edges$link > 0, ]
  list(nodes = nodes, edges = edges)
}