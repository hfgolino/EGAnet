#------------------------------------------
## S3Methods plot() // Updated 16.06.2021
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
#' \method{plot}{bootEGA}(x, plot.type = c("GGally", "qgraph"), plot.args = list(), ...)
#'
#' \method{plot}{CFA}(x, layout = "spring", vsize = 6, ...)
#'
#' \method{plot}{dynEGA}(x, title = "", plot.type = c("GGally", "qgraph"), plot.args = list(), ...)
#'
#' \method{plot}{dynEGA.Groups}(x, ncol, nrow, title = "",
#' plot.type = c("GGally", "qgraph"), plot.args = list(),  ...)
#'
#' \method{plot}{dynEGA.Individuals}(x, title = "",  id = NULL,
#' plot.type = c("GGally", "qgraph"), plot.args = list(), ...)
#'
#' \method{plot}{EGA}(x, title = "", plot.type = c("GGally", "qgraph"), plot.args = list(), ...)
#'
#' \method{plot}{NetLoads}(x, ...)
#'
#' @description Plots for \code{EGAnet} objects
#'
#' @param x Object from \code{EGAnet} package
#' 
#' @param vsize Numeric.
#' Size of vertices in \code{\link[EGAnet]{CFA}} plots.
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
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#' 
#' \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#' 
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for 
#' more details and examples}
#' 
#' }
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
# PLOTS----
# Updated 16.06.2021
#' @export
# Plot bootEGA----
# Updated 16.06.2021
plot.bootEGA <- function(x, plot.type = c("GGally","qgraph"),
                         plot.args = list(), produce = TRUE, ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$typicalGraph$graph, layout = "spring",
                               groups = as.factor(x$typicalGraph$wc),
                               vsize = plot.args$vsize, ...)
  }else if(plot.type == "GGally"){
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$typicalGraph$graph <- ifelse(abs(as.matrix(x$typicalGraph$graph)) <= .00001, 0, as.matrix(x$typicalGraph$graph))  
    
    # Reorder network and communities
    x$typicalGraph$graph <- x$typicalGraph$graph[order(x$typicalGraph$wc), order(x$typicalGraph$wc)]
    x$typicalGraph$wc <- x$typicalGraph$wc[order(x$typicalGraph$wc)]
    
    # weighted  network
    network1 <- network::network(x$typicalGraph$graph,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    
    if(exists("legend.names")){
      for(l in 1:max(x$typicalGraph$wc, na.rm = TRUE)){
        x$typicalGraph$wc[x$typicalGraph$wc == l] <- legend.names[l]
      }
    }
    
    network::set.vertex.attribute(network1, attrname= "Communities", value = x$typicalGraph$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$typicalGraph$graph))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(rescale.edges(x$typicalGraph$graph, plot.args$edge.size),
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$node.shape <- plot.args$shape
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- NULL
    plot.args$palette <- NULL
    
    lower <- abs(x$typicalGraph$graph[lower.tri(x$typicalGraph$graph)])
    non.zero <- sqrt(lower[lower != 0])
    
    plot.args$edge.alpha <- non.zero
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$typicalGraph$graph)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$node.size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$node.size/2}
    
    ega.plot <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
  }
  set.seed(NULL)
  
  name <- colnames(x$typicalGraph$graph)
  
  # Custom nodes: transparent insides and dark borders
  ega.plot <- ega.plot + 
    ggplot2::geom_point(ggplot2::aes(color = color), size = node.size,
               color = color_palette_EGA(color.palette, na.omit(x$typicalGraph$wc)),
               shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = color), size = node.size + .5,
               color = color_palette_EGA(color.palette, na.omit(x$typicalGraph$wc)),
               shape = 19, alpha = plot.args$alpha) +
    ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot dynEGA function (Level: Group)----
# Updated 16.06.2021
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", plot.type = c("GGally","qgraph"),
                               plot.args = list(), produce = TRUE, ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  
  ### Plot ###
  if(plot.type == "qgraph"){
    par(mfrow=c(nrow,ncol))
    for(i in 1:length(x$dynEGA)){
      qgraph::qgraph(x$dynEGA[[i]]$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$dynEGA[[i]]$wc), ...)
      title(names(x$dynEGA)[[i]], ...)}
  } else if(plot.type == "GGally"){
    
    plots.net <- vector("list", length = length(x$dynEGA))
    network1 <- vector("list", length = length(x$dynEGA))
    graph1 <- vector("list", length = length(x$dynEGA))
    edge.list <- vector("list", length = length(x$dynEGA))
    layout.spring <- vector("list", length = length(x$dynEGA))
    
    for(i in 1:length(x$dynEGA)){
      
      # Insignificant values (keeps ggnet2 from erroring out)
      x$dynEGA[[i]]$network <- ifelse(abs(as.matrix(x$dynEGA[[i]]$network)) <= .00001, 0, as.matrix(x$dynEGA[[i]]$network))
      
      # Reorder network and communities
      x$dynEGA[[i]]$network <- x$dynEGA[[i]]$network[order(x$dynEGA[[i]]$wc), order(x$dynEGA[[i]]$wc)]
      x$dynEGA[[i]]$wc <- x$dynEGA[[i]]$wc[order(x$dynEGA[[i]]$wc)]
      
      # weighted  network
      network1[[i]] <- network::network(x$dynEGA[[i]]$network,
                                        ignore.eval = FALSE,
                                        names.eval = "weights",
                                        directed = FALSE)
      
      network::set.vertex.attribute(network1[[i]], attrname= "Communities", value = x$dynEGA[[i]]$wc)
      network::set.vertex.attribute(network1[[i]], attrname= "Names", value = network::network.vertex.names(network1[[i]]))
      network::set.edge.attribute(network1[[i]], "color", ifelse(network::get.edge.value(network1[[i]], "weights") > 0, "darkgreen", "red"))
      network::set.edge.value(network1[[i]], attrname="AbsWeights",value=abs(x$dynEGA[[i]]$network))
      network::set.edge.value(network1[[i]],attrname="ScaledWeights",
                              value=matrix(rescale.edges(x$dynEGA[[i]]$network, default.args$edge.size),
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
      plot.args$net <- network1[[i]]
      plot.args$node.color <- "Communities"
      plot.args$node.alpha <- plot.args$alpha
      node.size <- plot.args$node.size
      plot.args$node.size <- 0
      plot.args$node.shape <- plot.args$shape
      plot.args$edge.color <- "color"
      plot.args$edge.size <- "ScaledWeights"
      plot.args$color.palette <- NULL
      plot.args$palette <- NULL
      
      lower <- abs(x$dynEGA[[i]]$network[lower.tri(x$dynEGA[[i]]$network)])
      non.zero <- sqrt(lower[lower != 0])
      
      plot.args$edge.alpha <- non.zero
      plot.args$mode <- layout.spring[[i]]
      plot.args$label <- colnames(x$dynEGA[[i]]$network)
      plot.args$node.label <- plot.args$label
      if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$node.size/2}
      if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$node.size/2}
      
      
      palette <- color_palette_EGA(color.palette, x$dynEGA[[i]]$wc)
      palette <- ifelse(is.na(palette), "white", palette)
      
      plots.net[[i]] <- suppressWarnings(
        suppressMessages(
          do.call(GGally::ggnet2, plot.args) + 
            ggplot2::theme(legend.title = ggplot2::element_blank())
        )
      )
      
      name <- colnames(x$dynEGA[[i]]$network)
      
      # Custom nodes: transparent insides and dark borders
      plots.net[[i]] <- plots.net[[i]] + 
        ggplot2::geom_point(ggplot2::aes(color = color), size = node.size,
                            color = palette,
                            shape = 1, stroke = 1.5, alpha = .8) +
        ggplot2::geom_point(ggplot2::aes(color = color), size = node.size + .5,
                            color = palette,
                            shape = 19, alpha = plot.args$alpha) +
        ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
        ggplot2::guides(
          color = ggplot2::guide_legend(override.aes = list(
            color = unique(palette[order(x$dynEGA[[i]]$wc)]),
            size = node.size,
            alpha = plot.args$alpha,
            stroke = 1.5
          ))
        )
      
    }
    group.labels <- names(x$dynEGA)
    set.seed(NULL)
    if(isTRUE(produce)){
      ggpubr::ggarrange(plotlist=plots.net, ncol = ncol, nrow = nrow, labels = group.labels, label.x = 0.3)
    }else{return(plots.net)}
  }
}

# Plot dynEGA function (Level: Individual)----
# Updated 16.06.2021
#' @export
plot.dynEGA.Individuals <- function(x, title = "",  id = NULL, plot.type = c("GGally","qgraph"),
                                    plot.args = list(), produce = TRUE, ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  ### Plot ###
  if(plot.type == "qgraph"){
    plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)
  }else if(plot.type == "GGally"){
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$dynEGA[[id]]$network <- ifelse(abs(as.matrix(x$dynEGA[[id]]$network)) <= .00001, 0, as.matrix(x$dynEGA[[id]]$network))
    
    # Reorder network and communities
    x$dynEGA[[id]]$network <- x$dynEGA[[id]]$network[order(x$dynEGA[[id]]$wc), order(x$dynEGA[[id]]$wc)]
    x$dynEGA[[id]]$wc <- x$dynEGA[[id]]$wc[order(x$dynEGA[[id]]$wc)]
    
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
                            value=matrix(rescale.edges(x$dynEGA[[id]]$network, plot.args$edge.size),
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$node.shape <- plot.args$shape
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- NULL
    plot.args$palette <- NULL
    
    lower <- abs(x$dynEGA[[id]]$network[lower.tri(x$dynEGA[[id]]$network)])
    non.zero <- sqrt(lower[lower != 0])
    
    plot.args$edge.alpha <- non.zero
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$dynEGA[[id]]$network)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$node.size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$node.size/2}
    
    palette <- color_palette_EGA(color.palette, x$dynEGA[[id]]$wc)
    palette <- ifelse(is.na(palette), "white", palette)
    
    ega.plot <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
    name <- colnames(x$dynEGA[[id]]$network)
    
    # Custom nodes: transparent insides and dark borders
    ega.plot <- ega.plot + 
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size,
                          color = palette,
                          shape = 1, stroke = 1.5, alpha = .8) +
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size + .5,
                          color = palette,
                          shape = 19, alpha = plot.args$alpha) +
      ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette[order(x$dynEGA[[id]]$wc)]),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ))
      )
    
    set.seed(NULL)
    
    if(isTRUE(produce)){
      plot(ega.plot)
    }else{return(ega.plot)}
  }
}

# Plot dynEGA function (Level: Population)----
# Updated 16.06.2021
#' @export
plot.dynEGA <- function(x, title = "", plot.type = c("GGally","qgraph"),
                        plot.args = list(), produce = TRUE, ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$dynEGA$wc), ...)
  }else if(plot.type == "GGally"){
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$dynEGA$network <- ifelse(abs(as.matrix(x$dynEGA$network)) <= .00001, 0, as.matrix(x$dynEGA$network))
    
    # Reorder network and communities
    x$dynEGA$network <- x$dynEGA$network[order(x$dynEGA$wc), order(x$dynEGA$wc)]
    x$dynEGA$wc <- x$dynEGA$wc[order(x$dynEGA$wc)]
    
    # weighted  network
    network1 <- network::network(x$dynEGA$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    if(exists("legend.names")){
      for(l in 1:x$dynEGA$n.dim){
        x$dynEGA$wc[x$dynEGA$wc == l] <- legend.names[l]
      }
    }
    
    network::set.vertex.attribute(network1, attrname= "Communities", value = x$dynEGA$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$dynEGA$network))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(rescale.edges(x$dynEGA$network, plot.args$edge.size),
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- NULL
    plot.args$palette <- NULL
    
    lower <- abs(x$dynEGA$network[lower.tri(x$dynEGA$network)])
    non.zero <- sqrt(lower[lower != 0])
    
    plot.args$edge.alpha <- non.zero
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$dynEGA$network)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}

    palette <- color_palette_EGA(color.palette, x$dynEGA$wc)
    palette <- ifelse(is.na(palette), "white", palette)
    
    ega.plot <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
    name <- colnames(x$dynEGA$network)
    
    # Custom nodes: transparent insides and dark borders
    ega.plot <- ega.plot + 
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size,
                          color = palette,
                          shape = 1, stroke = 1.5, alpha = .8) +
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size + .5,
                          color = palette,
                          shape = 19, alpha = plot.args$alpha) +
      ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette[order(x$dynEGA$wc)]),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ))
      )
    
  }
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot EGA----
# Updated 16.06.2021
#' @export
plot.EGA <- function(x,  title = "", plot.type = c("GGally","qgraph"),
                     plot.args = list(), produce = TRUE, ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$wc), ...)
  }else if(plot.type == "GGally"){
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$network <- ifelse(abs(as.matrix(x$network)) <= .00001, 0, as.matrix(x$network))
    
    if(exists("legend.names")){
      for(l in 1:length(unique(legend.names))){
        x$wc[x$wc == l] <- legend.names[l]
      }
    }
    
    # Reorder network and communities
    x$network <- x$network[order(x$wc), order(x$wc)]
    x$wc <- x$wc[order(x$wc)]
    
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
                            value=matrix(rescale.edges(x$network, plot.args$edge.size),
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$color.palette <- NULL
    plot.args$palette <- NULL
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    
    lower <- abs(x$network[lower.tri(x$network)])
    non.zero <- sqrt(lower[lower != 0])
    
    plot.args$edge.alpha <- non.zero
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$network)
    plot.args$node.label <- rep("", ncol(x$network))
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    palette <- color_palette_EGA(color.palette, x$wc)
    palette <- ifelse(is.na(palette), "white", palette)
    
    ega.plot <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
    name <- colnames(x$network)
    
    # Custom nodes: transparent insides and dark borders
    ega.plot <- ega.plot + 
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size,
                          color = palette,
                          shape = 1, stroke = 1.5, alpha = .8) +
      ggplot2::geom_point(ggplot2::aes(color = color), size = node.size + .5,
                          color = palette,
                          shape = 19, alpha = plot.args$alpha) +
      ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ))
      )
    
  }
  
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot EGA.fit----
# Updated 17.03.2021
#' @export
plot.EGA.fit <- function(x,  title = "", plot.type = c("GGally","qgraph"),
                     plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  plot.EGA(x$EGA, plot.type = plot.type, plot.args = plot.args)
}

#Plot net.loads----
# Updated 02.05.2020
#' @export
plot.NetLoads <- function(x, ...) {

  plot(x$plot)
}

#Plot CFA----
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}
