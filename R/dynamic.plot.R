#'  Dynamic Plot method for EGA objects.
#'
#' \code{dynamic.plot} Plots the EGA result using \code{\link{plotly}}
#'
#' @param ega.obj An EGA object
#' @param title Character. Title of the plot
#' @param vsize An integer indicating the size of the nodes. Default vsize = 30.
#' @param opacity A numeric value indicating the opacity of the edges. Default opacity = 0.4
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.EGA = TRUE)
#' summary(ega.wmt)
#' dynamic.plot(ega.wmt, title = "", vsize = 30, opacity = 0.4)
#'
#' \dontrun{
#' dynamic.plot(EGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

## S3 method for class 'EGA'

dynamic.plot <- function(ega.obj, title = "", vsize = 30, opacity = 0.4){
  require(qgraph)
  require(igraph)
  require(plotly)

  graph.glasso <- as.igraph(qgraph(ega.obj$glasso, DoNotPlot = TRUE))
  vert <- V(graph.glasso)
  es <- as.data.frame(get.edgelist(graph.glasso))
  edge.width <- E(graph.glasso)$weight
  L <- qgraph.layout.fruchtermanreingold(edgelist = as.matrix(es),
                                         weights = edge.width, vcount = length(ega.obj$wc))
  Nv <- length(vert)
  Ne <- length(es[1]$V1)
  Xn <- L[,1]
  Yn <- L[,2]
  network <- plot_ly(x = ~Xn, y = ~Yn, mode = "markers", text = paste("Variable: ",vert$label), hoverinfo = "text",
                     color = as.factor(ega.obj$wc),
                     marker = list(size = vsize,
                                   width = 2)) %>%
    add_annotations(x = Xn,
                    y = Yn,
                    text = vert$label,
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE,
                    ax = 20,
                    ay = -40)
  edge_shapes <- list()
  for(i in 1:Ne) {
    v0 <- es[i,]$V1
    v1 <- es[i,]$V2
    edge_shape = list(opacity = opacity,
                      type = "line",
                      line = list(color = ifelse(edge.width[i]>=0, "green", "red"), width = abs(edge.width[i])*10,
                                  hoverinfo = "text", color = "black",
                                  hoverlabel = list(bgcolor = "white"),
                                  text = ~paste("R.Part.Cor.:", round(edge.width[i],3))),
                      x0 = Xn[v0],
                      y0 = Yn[v0],
                      x1 = Xn[v1],
                      y1 = Yn[v1]
    )
    edge_shapes[[i]] <- edge_shape
  }
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  plot <- layout(
    network,
    title = title,
    shapes = edge_shapes,
    xaxis = axis,
    yaxis = axis,
    legend = list(x = 100, y = 0.5)
  )
  plot
}
