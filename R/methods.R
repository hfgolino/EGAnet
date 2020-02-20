# Methods:

#Summary function for dynEGA
summary.dynEGA <- function(object, ...) {
  cat("dynEGA Results (Random Effects):\n")
  cat("\nNumber of Dimensions:\n")
  print(object$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dynEGA$dim.variables)
}

#Summary function for dynEGA (Group Effects)
summary.dynEGA.Groups <- function(object, ...) {
  for(i in 1:length(object$dynEGA)){
    cat("dynEGA Results (Group Effects):\n")
    cat("Group:", names(object$dynEGA[i]))
    cat("\nNumber of Dimensions:\n")
    print(object$dynEGA[[i]]$n.dim)
    cat("\nItems per Dimension:\n")
    print(object$dynEGA[[i]]$dim.variables)
  }
}

#Summary function for dynEGA (Fixed Effects - Intraindividual Structure)
summary.dynEGA.Individuals <- function(object, ...) {
  cat("Number of Cases (individuals): \n")
  number <- length(object$dynEGA)
  print(number)
  cat("Summary statistics (number of factors/communities): \n")
  dim <- sapply(object$dynEGA, "[[", 3)
  cat("Mean:", mean(dim), "\n")
  cat("Median:", median(dim), "\n")
  cat("Min:", min(dim), "\n")
  cat("Max:", max(dim), "\n")
}

#Print dynEGA.Groups function
print.dynEGA.Groups <- function(x, ...) {
  for(i in 1:length(x$dynEGA)){
    cat("dynEGA Results (Group Effects):\n")
    cat("Group:", names(x$dynEGA[i]))
    cat("\nNumber of Dimensions:\n")
    print(x$dynEGA[[i]]$n.dim)
    cat("\nItems per Dimension:\n")
    print(x$dynEGA[[i]]$dim.variables)
  }
}

#Print dynEGA.Individuals function

print.dynEGA.Individuals <- function(x, ...) {
  cat("Number of Cases (individuals): \n")
  number <- length(x$dynEGA)
  print(number)
  cat("Summary statistics (number of factors/communities): \n")
  dim <- sapply(x$dynEGA, "[[", 3)
  cat("Mean:", mean(dim), "\n")
  cat("Median:", median(dim), "\n")
  cat("Min:", min(dim), "\n")
  cat("Max:", max(dim), "\n")
}

#Print dynEGA function
print.dynEGA<- function(x, ...) {
  cat("dynEGA Results (Random Effects):\n")
  cat("\nNumber of Dimensions:\n")
  print(x$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dynEGA$dim.variables)
}

#Plot dynEGA function (Group Effects)
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", vsize = 6,  ...) {
  par(mfrow=c(nrow,ncol))
  for(i in 1:length(x$dynEGA)){
    qgraph::qgraph(x$dynEGA[[i]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[i]]$wc), ...)
    title(names(x$dynEGA)[[i]], ...)}
}

#Plot dynEGA function (Random Effects)
plot.dynEGA <- function(x, title = "", vsize = 6,  ...) {
  plot.dynEGA <- qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA$wc), ...)

}

#Plot dynEGA function (Fixed Effects)
plot.dynEGA.Individuals <- function(x, title = "", vsize = 6,  id, ...) {
  plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)

}

# Plot EGA:

plot.EGA <- function(ega.obj, title = "", vsize = 6,  ...){
  plot.ega <- qgraph::qgraph(ega.obj$network, layout = "spring", vsize = vsize, groups = as.factor(ega.obj$wc), ...)
}


# Plot bootEGA:
plot.bootEGA <- function(bootega.obj, title = "", vsize = 6,  ...){
  qgraph::qgraph(bootega.obj$typicalGraph$graph, layout = "spring",
         groups = as.factor(bootega.obj$typicalGraph$wc),
         vsize = vsize, ...)

}

# Dynamic Plot:
dynamic.plot <- function(ega.obj, title = "", vsize = 30, opacity = 0.4){

  graph.glasso <- NetworkToolbox::convert2igraph(ega.obj$network)
  vert <- igraph::V(graph.glasso)
  es <- as.data.frame(igraph::get.edgelist(graph.glasso))
  edge.width <- igraph::E(graph.glasso)$weight
  L <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = as.matrix(es),
                                         weights = edge.width, vcount = length(ega.obj$wc))
  Nv <- length(vert)
  Ne <- length(es[1]$V1)
  Xn <- L[,1]
  Yn <- L[,2]
  network <- plotly::plot_ly(x = ~Xn, y = ~Yn, mode = "markers", text = paste("Variable: ",vert$label), hoverinfo = "text",
                     color = as.factor(ega.obj$wc),
                     marker = list(size = vsize,
                                   width = 2)) %>%
    plotly::add_annotations(x = Xn,
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
  plot <- plotly::layout(
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
  semPlot::semPaths(object$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
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


# #Summary bootEGA:
# summary.bootEGA <- function(object, ...) {
#   cat("bootEGA Results:\n")
#   cat("\nNumber of Dimensions:\n")
#   print(object$n.dim)
#   cat("\nItems per Dimension:\n")
#   print(object$dim.variables)
# }
#
# #Print EGA:
# print.EGA <- function(object, ...) {
#   cat("EGA Results:\n")
#   cat("\nNumber of Dimensions:\n")
#   print(object$n.dim)
#   cat("\nItems per Dimension:\n")
#   print(object$dim.variables)
# }
