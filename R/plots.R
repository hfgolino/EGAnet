#------------------------------------------
## S3Methods plot() // Updated 15.10.2020
#------------------------------------------

#' S3Methods for Plotting
#'
#' @name plots
#'
#' @aliases
#' plot.bootEGA
#' plot.CFA
#' plot.dynEGA
#' plot.dynEGA.Groups
#' plot.dynEGA.Individuals
#' plot.EGA
#' plot.NetLoads
#'
#' @usage
#' \method{plot}{bootEGA}(x, vsize = 6, plot = c("GGally", "qgraph"), ...)
#'
#' \method{plot}{CFA}(x, layout = "spring", vsize = 6, ...)
#'
#' \method{plot}{dynEGA}(x, title = "", vsize = 6, plot = c("GGally", "qgraph"), ...)
#'
#' \method{plot}{dynEGA.Groups}(x, ncol, nrow, title = "", vsize = 6, plot = c("GGally", "qgraph"),  ...)
#'
#' \method{plot}{dynEGA.Individuals}(x, title = "", vsize = 6,  id = NULL, plot = c("GGally", "qgraph"), ...)
#'
#' \method{plot}{EGA}(x, title = "", vsize = 6, plot = c("GGally", "qgraph"), ...)
#'
#' \method{plot}{NetLoads}(x, ...)
#'
#' @description Plots for \code{EGAnet} objects
#'
#' @param x Object from \code{EGAnet} package
#'
#' @param vsize Numeric.
#' Size of vertices in network plots.
#' Defaults to \code{6}
#'
#' @param layout Character.
#' Layout of plot (see \code{\link[semPlot]{semPaths}}).
#' Defaults to "spring"
#'
#' @param ncol Numeric.
#' Number of columns
#'
#' @param nrow Numeric.
#' Number of rows
#'
#' @param title Character.
#' Title of the plot.
#' Defaults to \code{""}
#'
#' @param id Numeric.
#' An integer or character indicating the ID of the individual to plot
#'
#' @param plot Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param ... Arguments passed on to
#'
#' \itemize{
#'
#' \item{\code{\link[qgraph]{qgraph}}}
#' {Functions: bootEGA, dynEGA, dynEGA.Groups, dynEGA.Individuals, EGA, and net.loads}
#'
#' \item{\code{\link[semPlot]{semPaths}}}
#' {Functions: CFA}
#'
#' }
#'
#' @return Plots of \code{EGAnet} object
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics plot
#'
#Plot bootEGA----
# Updated 02.05.2020
#' @export
plot.bootEGA <- function(x, vsize = 6, plot = c("GGally","qgraph"),...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  ### Plot ###
  if(plot == "qgraph"){
    qgraph::qgraph(x$typicalGraph$graph, layout = "spring",
                   groups = as.factor(x$typicalGraph$wc),
                   vsize = vsize, ...)
  }else if(plot == "GGally"){
    # weighted  network
    network1 <- network::network(x$typicalGraph$graph,
                        ignore.eval = FALSE,
                        names.eval = "weights",
                        directed = FALSE)
    network::set.vertex.attribute(network1, attrname= "Communities", value = x$typicalGraph$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$typicalGraph$graph))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(scales::rescale(as.vector(x$typicalGraph$graph),
                                                         to = c(.001, 1.75)),
                                         nrow = nrow(x$typicalGraph$graph),
                                         ncol = ncol(x$typicalGraph$graph)))

    # Layout "Spring"
    graph1 <- NetworkToolbox::convert2igraph(x$typicalGraph$graph)
    edge.list <- igraph::as_edgelist(graph1)
    layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                               weights =
                                                                 abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                               vcount = ncol(x$typicalGraph$graph))


    set.seed(1234)
    GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                           color = "Communities", edge.color = c("color"),
                           alpha = 0.7, size = 12, edge.alpha = 0.4,
                           mode =  layout.spring,
                           label.size = 5,
                           label = colnames(x$typicalGraph$graph)) + ggplot2::theme(legend.title = ggplot2::element_blank())

  }
}

#Plot CFA----
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}

#Plot dynEGA function (Level: Group)----
# Updated 15.10.2020
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", vsize = 6, plot = c("GGally","qgraph"),...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  ### Plot ###
  if(plot == "qgraph"){
  par(mfrow=c(nrow,ncol))
  for(i in 1:length(x$dynEGA)){
    qgraph::qgraph(x$dynEGA[[i]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[i]]$wc), ...)
    title(names(x$dynEGA)[[i]], ...)}
  } else if(plot == "GGally"){
    par(mfrow=c(nrow,ncol))
    for(i in 1:length(x$dynEGA)){
    # weighted  network
    network1[[i]] <- network::network(x$dynEGA[[i]]$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)

    network::set.vertex.attribute( network1[[i]] , attrname= "Communities", value = x$dynEGA[[i]]$wc)
    network::set.vertex.attribute( network1[[i]] , attrname= "Names", value = network::network.vertex.names(network1[[i]]))
    network::set.edge.attribute( network1[[i]] , "color", ifelse(network::get.edge.value(network1[[i]], "weights") > 0, "darkgreen", "red"))
    network::set.edge.value( network1[[i]] ,attrname="AbsWeights",value=abs(x$dynEGA[[i]]$network))
    network::set.edge.value(network1[[i]],attrname="ScaledWeights",
                            value=matrix(scales::rescale(as.vector(x$dynEGA[[i]]$network),
                                                         to = c(.001, 1.75)),
                                         nrow = nrow(x$dynEGA[[i]]$network),
                                         ncol = ncol(x$dynEGA[[i]]$network)))

    # Layout "Spring"
    graph1[[i]] <- NetworkToolbox::convert2igraph(x$dynEGA[[i]]$network)
    edge.list[[i]] <- igraph::as_edgelist(graph1[[i]])
    layout.spring[[i]] <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list[[i]],
                                                               weights =
                                                                 abs(igraph::E(graph1[[i]])$weight/max(abs(igraph::E(graph1[[i]])$weight)))^2,
                                                               vcount = ncol(x$dynEGA[[i]]$network))


    set.seed(1234)
    GGally::ggnet2(network1[[i]], edge.size = "ScaledWeights", palette = "Set1",
                   color = "Communities", edge.color = c("color"),
                   alpha = 0.7, size = 12, edge.alpha = 0.4,
                   mode =  layout.spring[[i]],
                   label.size = 5,
                   label = colnames(x$dynEGA[[i]]$network))+ggplot2::theme(legend.title = ggplot2::element_blank())
    }
  }
}

#Plot dynEGA function (Level: Individual)----
# Updated 10.15.2020
#' @export
plot.dynEGA.Individuals <- function(x, title = "", vsize = 6,  id = NULL, plot = c("GGally","qgraph"),...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  ### Plot ###
  if(plot == "qgraph"){
    plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)
   }else if(plot == "GGally"){
   # weighted  network
   network1 <- network::network(x$dynEGA[[id]]$network,
                                ignore.eval = FALSE,
                                names.eval = "weights",
                                directed = FALSE)

   network::set.vertex.attribute(network1, attrname= "Communities", value = x$dynEGA[[id]]$wc)
   network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
   network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
   network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$dynEGA[[id]]$network))
   network::set.edge.value(network1,attrname="ScaledWeights",
                           value=matrix(scales::rescale(as.vector(x$dynEGA[[id]]$network),
                                                        to = c(.001, 1.75)),
                                        nrow = nrow(x$dynEGA[[id]]$network),
                                        ncol = ncol(x$dynEGA[[id]]$network)))

   # Layout "Spring"
   graph1 <- NetworkToolbox::convert2igraph(x$dynEGA[[id]]$network)
   edge.list <- igraph::as_edgelist(graph1)
   layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                              weights =
                                                                abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                              vcount = ncol(x$dynEGA[[id]]$network))


   set.seed(1234)
   GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                  color = "Communities", edge.color = c("color"),
                  alpha = 0.7, size = 12, edge.alpha = 0.4,
                  mode =  layout.spring,
                  label.size = 12,
                  label = colnames(x$dynEGA[[id]]$network))+ggplot2::theme(legend.title = ggplot2::element_blank())
 }
}

#Plot dynEGA function (Level: Population)----
# Updated 10.15.2020
#' @export
plot.dynEGA <- function(x, title = "", vsize = 6,  plot = c("GGally","qgraph"),...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  ### Plot ###
  if(plot == "qgraph"){
    qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = vsize, groups = as.factor(x$dynEGA$wc), ...)
    }else if(plot == "GGally"){
  # weighted  network
  network1 <- network::network(x$dynEGA$network,
                               ignore.eval = FALSE,
                               names.eval = "weights",
                               directed = FALSE)

  network::set.vertex.attribute(network1, attrname= "Communities", value = x$dynEGA$wc)
  network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$dynEGA$network))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(scales::rescale(as.vector(x$dynEGA$network),
                                                       to = c(.001, 1.75)),
                                       nrow = nrow(x$dynEGA$network),
                                       ncol = ncol(x$dynEGA$network)))

  # Layout "Spring"
  graph1 <- NetworkToolbox::convert2igraph(x$dynEGA$network)
  edge.list <- igraph::as_edgelist(graph1)
  layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                             weights =
                                                               abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                             vcount = ncol(x$dynEGA$network))


  set.seed(1234)
  GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                 color = "Communities", edge.color = c("color"),
                 alpha = 0.7, size = 12, edge.alpha = 0.4,
                 mode =  layout.spring,
                 label.size = 5,
                 label = colnames(x$dynEGA$network))+ggplot2::theme(legend.title = ggplot2::element_blank())
  }
}

#Plot EGA----
# Updated 02.05.2020
#' @export
plot.EGA <- function(x, title = "", vsize = 6,  plot = c("GGally","qgraph"),...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  ### Plot ###
  if(plot == "qgraph"){
  plot.ega <- qgraph::qgraph(x$network, layout = "spring", vsize = vsize, groups = as.factor(x$wc), ...)
  }else if(plot == "GGally"){
    # weighted  network
    network1 <- network::network(x$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)

    network::set.vertex.attribute(network1, attrname= "Communities", value = x$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$network))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(scales::rescale(as.vector(x$network),
                                                         to = c(.001, 1.75)),
                                         nrow = nrow(x$network),
                                         ncol = ncol(x$network)))

    # Layout "Spring"
    graph1 <- NetworkToolbox::convert2igraph(x$network)
    edge.list <- igraph::as_edgelist(graph1)
    layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                               weights =
                                                                 abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                               vcount = ncol(x$network))


    set.seed(1234)
    GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                   color = "Communities", edge.color = c("color"),
                   alpha = 0.7, size = 12, edge.alpha = 0.4,
                   mode =  layout.spring,
                   label.size = 5,
                   label = colnames(x$network))+ggplot2::theme(legend.title = ggplot2::element_blank())
  }
}

#Plot net.loads----
# Updated 02.05.2020
#' @export
plot.NetLoads <- function(x, ...) {

  plot(x$plot)
}
