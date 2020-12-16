# Methods:

#Summary function for dynEGA
summary.dynEGA <- function(object, ...) {
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(object$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dynEGA$dim.variables)
}

#Summary function for dynEGA (Level: Group)
summary.dynEGA.Groups <- function(object, ...) {
  for(i in 1:length(object$dynEGA)){
    cat("dynEGA Results (Level: Group):\n")
    cat("Group:", names(object$dynEGA[i]))
    cat("\nNumber of Dimensions:\n")
    print(object$dynEGA[[i]]$n.dim)
    cat("\nItems per Dimension:\n")
    print(object$dynEGA[[i]]$dim.variables)
  }
}

#Summary function for dynEGA (Level: Individual - Intraindividual Structure)
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
    cat("dynEGA Results (Level: Group):\n")
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
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(x$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dynEGA$dim.variables)
}

# Plot dynEGA function (Level: Group)----
# Updated 16.12.2020
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", plot.type = c("GGally","qgraph"),
                               plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
        plot.args[which(names(plot.args) == "legend.names")] <- NULL
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
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
      plot.args$net <- network1
      plot.args$node.color <- "Communities"
      plot.args$node.alpha <- plot.args$alpha
      plot.args$node.shape <- plot.args$shape
      plot.args$node.size <- plot.args$size
      plot.args$edge.color <- "color"
      plot.args$edge.size <- "ScaledWeights"
      plot.args$color.palette <- "Set1"
      plot.args$mode <- layout.spring
      plot.args$label <- colnames(x$dynEGA[[i]]$network)
      plot.args$node.label <- plot.args$label
      if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
      if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
      
      
      plots.net[[i]] <-  do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
      
    }
    group.labels <- names(x$dynEGA)
    set.seed(NULL)
    ggpubr::ggarrange(plotlist=plots.net, ncol = ncol, nrow = nrow, labels = group.labels, label.x = 0.3)
  }
}

# Plot dynEGA function (Level: Population)----
# Updated 16.12.2020
#' @export
plot.dynEGA <- function(x, title = "", plot.type = c("GGally","qgraph"),
                        plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
        plot.args[which(names(plot.args) == "legend.names")] <- NULL
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
  }
  
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$dynEGA$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$dynEGA$wc), ...)
  }else if(plot.type == "GGally"){
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    plot.args$node.size <- plot.args$size
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- "Set1"
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$dynEGA$network)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    ega.plot <- do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
    
  }
  set.seed(NULL)
  plot(ega.plot)
}

# Plot dynEGA function (Level: Individual)----
# Updated 16.12.2020
#' @export
plot.dynEGA.Individuals <- function(x, title = "",  id = NULL, plot.type = c("GGally","qgraph"),
                                    plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
        plot.args[which(names(plot.args) == "legend.names")] <- NULL
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
  }
  
  ### Plot ###
  if(plot.type == "qgraph"){
    plot.dynEGA.Individuals <- qgraph::qgraph(x$dynEGA[[id]]$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$dynEGA[[id]]$wc), ...)
  }else if(plot.type == "GGally"){
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    plot.args$node.size <- plot.args$size
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- "Set1"
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$dynEGA[[id]]$network)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    ega.plot <- do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
    
    set.seed(NULL)
    plot(ega.plot)
  }
}

# Plot EGA----
# Updated 16.12.2020
#' @export
plot.EGA <- function(x, title = "",  plot.type = c("GGally","qgraph"),
                     plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
        plot.args[which(names(plot.args) == "legend.names")] <- NULL
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
  }
  
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$wc), ...)
  }else if(plot.type == "GGally"){
    # weighted  network
    network1 <- network::network(x$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    if(exists("legend.names")){
      for(l in 1:x$n.dim){
        x$wc[x$wc == l] <- legend.names[l]
      }
    }
    
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    plot.args$node.size <- plot.args$size
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- "Set1"
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$network)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    ega.plot <- do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
    
  }
  set.seed(NULL)
  
  plot(ega.plot)
}

# Plot bootEGA----
# Updated 16.12.2020
plot.bootEGA <- function(x, plot.type = c("GGally","qgraph"),
                         plot.args = list(), ...){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
        plot.args[which(names(plot.args) == "legend.names")] <- NULL
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
  }
  
  ### Plot ###
  if(plot.type == "qgraph"){
    ega.plot <- qgraph::qgraph(x$typicalGraph$graph, layout = "spring",
                   groups = as.factor(x$typicalGraph$wc),
                   vsize = plot.args$vsize, ...)
  }else if(plot.type == "GGally"){
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
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    plot.args$node.size <- plot.args$size
    plot.args$edge.color <- "color"
    plot.args$edge.size <- "ScaledWeights"
    plot.args$color.palette <- "Set1"
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$typicalGraph$graph)
    plot.args$node.label <- plot.args$label
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    
    ega.plot <- do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
    
  }
  set.seed(NULL)
  plot(ega.plot)
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

#Summary NetLoads
#Updated 09.03.2020
summary.NetLoads <- function(object, ...) {
  object$std[which(abs(object$std) <= object$MinLoad, arr.ind = TRUE)] <- ""
  print(object$std)
  message("Loadings <= ", object$MinLoad, " are blank")
}

#Print NetLoads
#Updated 09.03.2020
print.NetLoads <- function(object, ...) {
  object$std[which(abs(object$std) <= object$MinLoad, arr.ind = TRUE)] <- ""
  print(object$std)
  message("Loadings <= ", object$MinLoad, " are blank")
}

#Plot function for NetLoads
#Updated 05.03.2020
plot.NetLoads <- function(object, ...) {
  plot(object$plot)
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
