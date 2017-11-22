#Methods:

#Plot EGA:
plot.EGA <- function(ega.obj, title = "", vsize = 30, opacity = 0.4){
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
  print(plot)
}



#Summary EGA:
summary.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}

#Print EGA:
print.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}

#Plot CFA:
plot.CFA <- function(object, layout = "spring", vsize = 6, ...) {
  semPaths(object$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
           sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
           "std", cut = 0.5)
}

#Summary CFA:
summary.CFA <- function(object, ...) {
  cat("Summary: Confirmatory Factor Analysis:\n")
  print(object$summary)
  cat("\n FIt Measures:\n")
  print(object$fit.measures)
}

