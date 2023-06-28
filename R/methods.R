#%%%%%%%%%%%%%%%%%%%%%%%%#
#### EGAnet S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%#

# Updated 07.06.2022 (DD.MM.YYYY)

# print() Methods ----

# Print dynEGA (Level: Population)
# Updated 13.05.2022
#' @export
print.dynEGA<- function(x, ...) {
  
  # Print communities
  cat(paste(
    "Number of communities (population-level):",
    x$dynEGA$n.dim,
    "\n\n"
  ))
  
  # Set up item placements
  item_placement <- x$dynEGA$wc
  names(item_placement) <- gsub(
    ".Ord*.", "", names(item_placement)
  )
  
  # Print item placements
  print(item_placement)
  
  # Print GLLA methods
  cat("\nGLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## Print GLLA
  print(glla.methods, quote = FALSE)
  
  # Print EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## Print EGA
  print(ega.methods, quote = FALSE)
  
}

# Print dynEGA (Level: Groups)
# Updated 13.05.2022
#' @export
print.dynEGA.Groups <- function(x, ...) {
  
  for(i in 1:(length(x$dynEGA) - 1)){
    
    # Print communities
    cat(paste(
      "Number of communities (group-level):",
      x$dynEGA[[i]]$n.dim, "\n",
      paste("Group:", names(x$dynEGA[i])),
      "\n\n"
    ))
    
    # Set up item placements
    item_placement <- x$dynEGA[[i]]$wc
    names(item_placement) <- gsub(
      ".Ord*.", "", names(item_placement)
    )
    
    # Print item placements
    print(
      item_placement
    )
    
    # Add space
    cat("\n")
    
  }
  
  # Print GLLA methods
  cat("GLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## Print GLLA
  print(glla.methods, quote = FALSE)
  
  # Print EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## Print EGA
  print(ega.methods, quote = FALSE)
  
}

# Print dynEGA (Level: Individuals)
# Updated 24.06.2022
#' @export
print.dynEGA.Individuals <- function(x, ...) {
  
  # Number of people
  cat("Number of Cases (individuals): ")
  number <- length(x$dynEGA) - 1
  cat(number, "\n")
  
  # Summary statistics
  cat("\nSummary statistics (number of communities): \n")
  dim <- unlist(
    lapply(x$dynEGA, function(y){
      y$n.dim
    })
  )
  
  ## Set up summary
  summary.methods <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(summary.methods) <- c(
    "Mean =",
    "Median =",
    "Min =",
    "Max ="
  )
  colnames(summary.methods) <- ""
  
  ## Insert values
  summary.methods["Mean =",] <- mean(dim, na.rm = TRUE)
  summary.methods["Median =",] <- median(dim, na.rm = TRUE)
  summary.methods["Min =",] <- min(dim, na.rm = TRUE)
  summary.methods["Max =",] <- max(dim, na.rm = TRUE)
  
  ## Print summary
  print(summary.methods, quote = FALSE)
  
  # Print GLLA methods
  cat("\nGLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## Print GLLA
  print(glla.methods, quote = FALSE)
  
  # Print EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## Print EGA
  print(ega.methods, quote = FALSE)
  
}

#Print Network Loadings
# Updated 13.05.2022
#' @export
print.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$minLoad, arr.ind = TRUE)] <- ""
  
  print(x$std)
  message("Loadings <= |", x$minLoad, "| are blank")
}

#Print Measurement Invariance
# Updated 10.02.2022
#' @export
print.invariance <- function(x, ...) {
  print(x$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#Print Hierarchical EGA
# Updated 16.11.2022
#' @export
print.hierEGA <- function(x, ...) {
  
  # Print lower order communities
  cat(paste(
    "Lower order communities:",
    x$lower_order$n.dim,
    "\n\n"
  ))
  print(x$lower_order$wc)
  
  # Print higher order communities
  cat(
    paste(
      "\nHigher order communities:",
      x$higher_order$EGA$n.dim,
      "\n\n"
    )
  )
  print(x$higher_order$EGA$wc)
  
  # Print methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 7, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method =",
    "Scores =",
    "Consensus Method =",
    "Consensus Iterations ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- switch(
    tolower(x$Methods$uni.method),
    "expand" = "expand correlation matrix",
    "le" = "leading eigenvalue",
    "louvain" = "louvain with consensus clustering"
  )
  methods.matrix["Scores =",] <- x$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", x$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- x$Methods$consensus.iter
  
  print(methods.matrix, quote = FALSE)
  
}

#Print Network Descriptives
# Updated 19.07.2022
#' @export
print.Descriptives <- function(x, ...)
{
  # Print weights
  cat("\nWeights:\n")
  
  ## Set up weights
  weights.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(weights.matrix) <- c(
    "Mean =",
    "Standard Deviation =",
    "Range =",
    "Density ="
  )
  colnames(weights.matrix) <- ""
  
  weights.matrix["Mean =",] <- x["Mean_weight"]
  weights.matrix["Standard Deviation =",] <- x["SD_weight"]
  weights.matrix["Range =",] <- paste(
    x["Min_weight"], "to", x["Max_weight"]
  )
  weights.matrix["Density =",] <- x["Density"]
  
  print(weights.matrix, quote = FALSE)
  
  # Print weights
  cat("\nGlobal Properties:\n")
  
  ## Set up weights
  global.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(global.matrix) <- c(
    "Average Shortest Path Length (ASPL) =",
    "Clustering Coefficient (CC) =",
    # "Small-world (Random) =",
    # "Small-world (Humphries & Gurney, 2008) =",
    "Small-world (Telesford et al., 2011) =",
    "R-squared Fit for Scale-free Network ="
  )
  colnames(global.matrix) <- ""
  
  global.matrix["Average Shortest Path Length (ASPL) =",] <- x["ASPL"]
  global.matrix["Clustering Coefficient (CC) =",] <- x["CC"]
  # global.matrix["Small-world (Random) =",] <- x["swn.rand"]
  # global.matrix["Small-world (Humphries & Gurney, 2008) =",] <- x["swn.HG"]
  global.matrix["Small-world (Telesford et al., 2011) =",] <- x["swn.TJHBL"]
  global.matrix["R-squared Fit for Scale-free Network =",] <- x["scale-free_R-sq"]
  
  print(global.matrix, quote = FALSE)
  
  ## Add interpretation
  cat("---")
  # cat("\nSmall-world (Random): 0 'not small-world' 1 'small-world' Inf")
  # cat("\nSmall-world (Humphries & Gurney, 2008): 0 'not small-world' 3 'small-world' Inf")
  cat("\nSmall-world (Telesford et al., 2011): -1 'lattice' 0 'random' 1; near 0 = small-world")
  
}

# summary() Methods ----

# summary dynEGA (Level: Population)
# Updated 13.05.2022
#' @export
summary.dynEGA <- function(object, ...) {
  
  # summary communities
  cat(paste(
    "Number of communities (population-level):",
    object$dynEGA$n.dim,
    "\n\n"
  ))
  
  # Set up item placements
  item_placement <- object$dynEGA$wc
  names(item_placement) <- gsub(
    ".Ord*.", "", names(item_placement)
  )
  
  # summary item placements
  print(item_placement)
  
  # summary GLLA methods
  cat("\nGLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- object$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- object$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- object$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  print(glla.methods, quote = FALSE)
  
  # summary EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    object$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    object$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- object$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- object$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  print(ega.methods, quote = FALSE)
  
}

# summary dynEGA (Level: Groups)
# Updated 13.05.2022
#' @export
summary.dynEGA.Groups <- function(object, ...) {
  
  for(i in 1:(length(object$dynEGA) - 1)){
    
    # summary communities
    cat(paste(
      "Number of communities (group-level):",
      object$dynEGA[[i]]$n.dim, "\n",
      paste("Group:", names(object$dynEGA[i])),
      "\n\n"
    ))
    
    # Set up item placements
    item_placement <- object$dynEGA[[i]]$wc
    names(item_placement) <- gsub(
      ".Ord*.", "", names(item_placement)
    )
    
    # summary item placements
    print(
      item_placement
    )
    
    # Add space
    cat("\n")
    
  }
  
  # summary GLLA methods
  cat("GLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- object$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- object$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- object$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  print(glla.methods, quote = FALSE)
  
  # summary EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    object$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    object$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- object$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- object$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  print(ega.methods, quote = FALSE)
  
}

# summary dynEGA (Level: Individuals)
# Updated 13.05.2022
#' @export
summary.dynEGA.Individuals <- function(object, ...) {
  
  # Number of people
  cat("Number of Cases (individuals): ")
  number <- length(object$dynEGA) - 1
  cat(number, "\n")
  
  # Summary statistics
  cat("\nSummary statistics (number of communities): \n")
  dim <- unlist(
    lapply(object$dynEGA, function(x){
      x$n.dim
    })
  )
  
  ## Set up summary
  summary.methods <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(summary.methods) <- c(
    "Mean =",
    "Median =",
    "Min =",
    "Max ="
  )
  colnames(summary.methods) <- ""
  
  ## Insert values
  summary.methods["Mean =",] <- mean(dim, na.rm = TRUE)
  summary.methods["Median =",] <- median(dim, na.rm = TRUE)
  summary.methods["Min =",] <- min(dim, na.rm = TRUE)
  summary.methods["Max =",] <- max(dim, na.rm = TRUE)
  
  ## summary summary
  print(summary.methods, quote = FALSE)
  
  # summary GLLA methods
  cat("\nGLLA Methods:\n")
  
  ## Set up GLLA methods
  glla.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(glla.methods) <- c(
    "Embedding Dimensions =",
    "Embedding Offset (tau) =",
    "Lag (delta) ="
  )
  colnames(glla.methods) <- ""
  
  ## Insert values
  glla.methods["Embedding Dimensions =",] <- object$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- object$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- object$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  print(glla.methods, quote = FALSE)
  
  # summary EGA methods
  cat("\nEGA Methods:\n")
  
  ## Set up EGA methods
  ega.methods <- matrix(
    nrow = 3, ncol = 1
  )
  row.names(ega.methods) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="
  )
  colnames(ega.methods) <- ""
  
  ## Insert values
  ega.methods["Correlations =",] <- ifelse(
    object$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    object$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- object$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- object$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  print(ega.methods, quote = FALSE)
  
}

#summary Network Loadings
# Updated 13.05.2022
#' @export
summary.NetLoads <- function(object, ...) {
  
  object$std[which(abs(object$std) <= object$minLoad, arr.ind = TRUE)] <- ""
  
  print(object$std)
  message("Loadings <= |", object$minLoad, "| are blank")
}

#summary Measurement Invariance
# Updated 10.02.2022
#' @export
summary.invariance <- function(object, ...) {
  print(object$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#summary Hierarchical EGA
# Updated 09.05.2022
#' @export
summary.hierEGA <- function(object, ...) {
  
  # Print lower order communities
  cat(paste(
    "Lower order communities:",
    object$lower_order$n.dim,
    "\n\n"
  ))
  print(object$lower_order$wc)
  
  # Print higher order communities
  cat(
    paste(
      "\nHigher order communities:",
      object$higher_order$EGA$n.dim,
      "\n\n"
    )
  )
  print(object$higher_order$EGA$wc)
  
  # Print methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 7, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method =",
    "Scores =",
    "Consensus Method =",
    "Consensus Iterations ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    object$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    object$Methods$corr
  )
  methods.matrix["Model =",] <- object$Methods$model
  methods.matrix["Algorithm =",] <- object$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- switch(
    tolower(object$Methods$uni.method),
    "expand" = "expand correlation matrix",
    "le" = "leading eigenvalue",
    "louvain" = "louvain with consensus clustering"
  )
  methods.matrix["Scores =",] <- object$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", object$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- object$Methods$consensus.iter
  
  print(methods.matrix, quote = FALSE)
  
}

#Summary Network Descriptives
# Updated 19.07.2022
#' @export
summary.Descriptives <- function(object, ...)
{
  # Print weights
  cat("\nWeights:\n")
  
  ## Set up weights
  weights.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(weights.matrix) <- c(
    "Mean =",
    "Standard Deviation =",
    "Range =",
    "Density ="
  )
  colnames(weights.matrix) <- ""
  
  weights.matrix["Mean =",] <- object["Mean_weight"]
  weights.matrix["Standard Deviation =",] <- object["SD_weight"]
  weights.matrix["Range =",] <- paste(
    object["Min_weight"], "to", object["Max_weight"]
  )
  weights.matrix["Density =",] <- object["Density"]
  
  print(weights.matrix, quote = FALSE)
  
  # Print weights
  cat("\nGlobal Properties:\n")
  
  ## Set up weights
  global.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(global.matrix) <- c(
    "Average Shortest Path Length (ASPL) =",
    "Clustering Coefficient (CC) =",
    # "Small-world (Random) =",
    # "Small-world (Humphries & Gurney, 2008) =",
    "Small-world (Telesford et al., 2011) =",
    "R-squared Fit for Scale-free Network ="
  )
  colnames(global.matrix) <- ""
  
  global.matrix["Average Shortest Path Length (ASPL) =",] <- object["ASPL"]
  global.matrix["Clustering Coefficient (CC) =",] <- object["CC"]
  # global.matrix["Small-world (Random) =",] <- x["swn.rand"]
  # global.matrix["Small-world (Humphries & Gurney, 2008) =",] <- x["swn.HG"]
  global.matrix["Small-world (Telesford et al., 2011) =",] <- object["swn.TJHBL"]
  global.matrix["R-squared Fit for Scale-free Network =",] <- object["scale-free_R-sq"]
  
  print(global.matrix, quote = FALSE)
  
  ## Add interpretation
  cat("---")
  # cat("\nSmall-world (Random): 0 'not small-world' 1 'small-world' Inf")
  # cat("\nSmall-world (Humphries & Gurney, 2008): 0 'not small-world' 3 'small-world' Inf")
  cat("\nSmall-world (Telesford et al., 2011): -1 'lattice' 0 'random' 1; near 0 = small-world")
  
}

# plot() Methods ----

# Plot bootEGA
# Updated 07.07.2022
plot.bootEGA <- function(x, title = "",
                         plot.args = list(), produce = TRUE, ...){
  
  # Organize plot arguments
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
  }
  if("legend.title" %in% names(plot.args)){
    legend.title <- plot.args$legend.title
  }
  plot.args <- GGally.args(plot.args)
  color.palette <- plot.args$color.palette
  
  ### Plot ###
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
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
  network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$typicalGraph$graph))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(rescale.edges(x$typicalGraph$graph, plot.args$edge.size),
                                       nrow = nrow(x$typicalGraph$graph),
                                       ncol = ncol(x$typicalGraph$graph)))
  
  # Layout "Spring"
  graph1 <- convert2igraph(x$typicalGraph$graph)
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
  plot.args$edge.lty <- "line"
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
  
  palette <- color_palette_EGA(color.palette, as.numeric(factor(x$typicalGraph$wc)))
  palette <- ifelse(is.na(palette), "white", palette)
  
  ega.plot <- suppressWarnings(
    suppressMessages(
      do.call(GGally::ggnet2, plot.args) + 
        ggplot2::theme(legend.title = ggplot2::element_blank())
    )
  )
  
  name <- colnames(x$typicalGraph$graph)
  
  name.split <- lapply(name, function(x){
    unlist(strsplit(x, split = " "))
  })
  
  name <- unlist(
    lapply(name.split, function(x){
      
      len <- length(x)
      
      if(len > 1){
        
        add.line <- round(len / 2)
        
        paste(
          paste(x[1:add.line], collapse = " "),
          paste(x[(add.line+1):length(x)], collapse = " "),
          sep = "\n"
        )
        
      }else{x}
      
    })
  )
  
  # Border color
  if(all(color.palette == "grayscale" |
         color.palette == "greyscale" |
         color.palette == "colorblind")){
    border.color <- ifelse(palette == "white", "white", "black")
  }else{border.color <- palette}
  
  # Custom nodes: transparent insides and dark borders
  ega.plot <- ega.plot + 
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                        color = border.color,
                        shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
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
  
  # Add legend and plot title
  if(title != ""){
    ega.plot <- ega.plot +
      ggplot2::labs(title = title)
  }
  
  if(exists("legend.title")){
    ega.plot <- ega.plot +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ), title = legend.title)
      ) +
      ggplot2::theme(
        legend.title = ggplot2::element_text()
      )
  }
  
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
  
}

# Plot dynEGA function (Level: Group)
# Updated 07.07.2022
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "",
                               plot.args = list(), produce = TRUE, ...){
  
  # Organize plot arguments
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
  }
  if("legend.title" %in% names(plot.args)){
    legend.title <- plot.args$legend.title
  }
  plot.args <- GGally.args(plot.args)
  color.palette <- plot.args$color.palette
  reset.plot.args <- plot.args
  
  if("methods" %in% tolower(names(x$dynEGA))){
    x$dynEGA <- x$dynEGA[-which(tolower(names(x$dynEGA)) == "methods")]
  }
  
  
  ### Plot ###
  plots.net <- vector("list", length = length(x$dynEGA))
  network1 <- vector("list", length = length(x$dynEGA))
  graph1 <- vector("list", length = length(x$dynEGA))
  edge.list <- vector("list", length = length(x$dynEGA))
  layout.spring <- vector("list", length = length(x$dynEGA))
  
  for(i in 1:length(x$dynEGA)){
    
    plot.args <- reset.plot.args
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$dynEGA[[i]]$network <- ifelse(abs(as.matrix(x$dynEGA[[i]]$network)) <= .00001, 0, as.matrix(x$dynEGA[[i]]$network))
    
    # Reorder network and communities
    x$dynEGA[[i]]$network <- x$dynEGA[[i]]$network[order(x$dynEGA[[i]]$wc), order(x$dynEGA[[i]]$wc)]
    x$dynEGA[[i]]$wc <- x$dynEGA[[i]]$wc[order(x$dynEGA[[i]]$wc)]
    
    # weighted  network
    network1 <- network::network(x$dynEGA[[i]]$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    network::set.vertex.attribute(network1, attrname= "Communities", value = x$dynEGA[[i]]$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
    network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
    network::set.edge.value(network1, attrname="AbsWeights",value=abs(x$dynEGA[[i]]$network))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(rescale.edges(x$dynEGA[[i]]$network, plot.args$edge.size),
                                         nrow = nrow(x$dynEGA[[i]]$network),
                                         ncol = ncol(x$dynEGA[[i]]$network)))
    
    # Layout "Spring"
    graph1[[i]] <- convert2igraph(x$dynEGA[[i]]$network)
    edge.list[[i]] <- igraph::as_edgelist(graph1[[i]])
    layout.spring[[i]] <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list[[i]],
                                                                    weights =
                                                                      abs(igraph::E(graph1[[i]])$weight/max(abs(igraph::E(graph1[[i]])$weight)))^2,
                                                                    vcount = ncol(x$dynEGA[[i]]$network))
    
    
    set.seed(1234)
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$node.shape <- plot.args$shape
    plot.args$edge.lty <- "line"
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
    
    palette <- color_palette_EGA(color.palette, as.numeric(factor(x$dynEGA[[i]]$wc)))
    palette <- ifelse(is.na(palette), "white", palette)
    
    plots.net[[i]] <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
    name <- colnames(x$dynEGA[[i]]$network)
    
    name.split <- lapply(name, function(x){
      unlist(strsplit(x, split = " "))
    })
    
    name <- unlist(
      lapply(name.split, function(x){
        
        len <- length(x)
        
        if(len > 1){
          
          add.line <- round(len / 2)
          
          paste(
            paste(x[1:add.line], collapse = " "),
            paste(x[(add.line+1):length(x)], collapse = " "),
            sep = "\n"
          )
          
        }else{x}
        
      })
    )
    
    # Border color
    if(all(color.palette == "grayscale" |
           color.palette == "greyscale" |
           color.palette == "colorblind")){
      border.color <- ifelse(palette == "white", "white", "black")
    }else{border.color <- palette}
    
    # Custom nodes: transparent insides and dark borders
    plots.net[[i]] <- plots.net[[i]] + 
      ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                          color = border.color,
                          shape = 1, stroke = 1.5, alpha = .8) +
      ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
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
    
    # Add legend and plot title
    if(title != ""){
      plots.net[[i]] <- plots.net[[i]] + 
        ggplot2::labs(title = title)
    }
    
    if(exists("legend.title")){
      plots.net[[i]] <- plots.net[[i]] + 
        ggplot2::guides(
          color = ggplot2::guide_legend(override.aes = list(
            color = unique(palette),
            size = node.size,
            alpha = plot.args$alpha,
            stroke = 1.5
          ), title = legend.title)
        ) +
        ggplot2::theme(
          legend.title = ggplot2::element_text()
        )
    }
    
  }
  group.labels <- names(x$dynEGA)
  set.seed(NULL)
  if(isTRUE(produce)){
    ggpubr::ggarrange(plotlist=plots.net, ncol = ncol, nrow = nrow, labels = group.labels, label.x = 0.3)
  }else{return(plots.net)}
}

# Plot dynEGA function (Level: Individual)
# Updated 07.07.2022
#' @export
plot.dynEGA.Individuals <- function(x, title = "",  id = NULL,
                                    plot.args = list(), produce = TRUE, ...){
  
  # Organize plot arguments
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
  }
  if("legend.title" %in% names(plot.args)){
    legend.title <- plot.args$legend.title
  }
  plot.args <- GGally.args(plot.args)
  color.palette <- plot.args$color.palette
  
  ### Plot ###
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
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
  network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$dynEGA[[id]]$network))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(rescale.edges(x$dynEGA[[id]]$network, plot.args$edge.size),
                                       nrow = nrow(x$dynEGA[[id]]$network),
                                       ncol = ncol(x$dynEGA[[id]]$network)))
  
  # Layout "Spring"
  graph1 <- convert2igraph(x$dynEGA[[id]]$network)
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
  plot.args$edge.lty <- "line"
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
  
  palette <- color_palette_EGA(color.palette, as.numeric(factor(x$dynEGA[[id]]$wc)))
  palette <- ifelse(is.na(palette), "white", palette)
  
  ega.plot <- suppressWarnings(
    suppressMessages(
      do.call(GGally::ggnet2, plot.args) + 
        ggplot2::theme(legend.title = ggplot2::element_blank())
    )
  )
  
  name <- colnames(x$dynEGA[[id]]$network)
  
  name.split <- lapply(name, function(x){
    unlist(strsplit(x, split = " "))
  })
  
  name <- unlist(
    lapply(name.split, function(x){
      
      len <- length(x)
      
      if(len > 1){
        
        add.line <- round(len / 2)
        
        paste(
          paste(x[1:add.line], collapse = " "),
          paste(x[(add.line+1):length(x)], collapse = " "),
          sep = "\n"
        )
        
      }else{x}
      
    })
  )
  
  # Border color
  if(all(color.palette == "grayscale" |
         color.palette == "greyscale" |
         color.palette == "colorblind")){
    border.color <- ifelse(palette == "white", "white", "black")
  }else{border.color <- palette}
  
  # Custom nodes: transparent insides and dark borders
  ega.plot <- ega.plot + 
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                        color = border.color,
                        shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
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
  
  # Add legend and plot title
  if(title != ""){
    ega.plot <- ega.plot +
      ggplot2::labs(title = title)
  }
  
  if(exists("legend.title")){
    ega.plot <- ega.plot +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ), title = legend.title)
      ) +
      ggplot2::theme(
        legend.title = ggplot2::element_text()
      )
  }
  
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
  
}

# Plot dynEGA function (Level: Population)
# Updated 07.07.2022
#' @export
plot.dynEGA <- function(x, title = "",
                        plot.args = list(), produce = TRUE, ...){
  
  # Organize plot arguments
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
  }
  if("legend.title" %in% names(plot.args)){
    legend.title <- plot.args$legend.title
  }
  plot.args <- GGally.args(plot.args)
  color.palette <- plot.args$color.palette
  
  
  ### Plot ###
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
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
  network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$dynEGA$network))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(rescale.edges(x$dynEGA$network, plot.args$edge.size),
                                       nrow = nrow(x$dynEGA$network),
                                       ncol = ncol(x$dynEGA$network)))
  
  # Layout "Spring"
  graph1 <- convert2igraph(x$dynEGA$network)
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
  plot.args$edge.lty <- "line"
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
  
  palette <- color_palette_EGA(color.palette, as.numeric(factor(x$dynEGA$wc)))
  palette <- ifelse(is.na(palette), "white", palette)
  
  ega.plot <- suppressWarnings(
    suppressMessages(
      do.call(GGally::ggnet2, plot.args) + 
        ggplot2::theme(legend.title = ggplot2::element_blank())
    )
  )
  
  name <- colnames(x$dynEGA$network)
  
  name.split <- lapply(name, function(x){
    unlist(strsplit(x, split = " "))
  })
  
  name <- unlist(
    lapply(name.split, function(x){
      
      len <- length(x)
      
      if(len > 1){
        
        add.line <- round(len / 2)
        
        paste(
          paste(x[1:add.line], collapse = " "),
          paste(x[(add.line+1):length(x)], collapse = " "),
          sep = "\n"
        )
        
      }else{x}
      
    })
  )
  
  # Border color
  if(all(color.palette == "grayscale" |
         color.palette == "greyscale" |
         color.palette == "colorblind")){
    border.color <- ifelse(palette == "white", "white", "black")
  }else{border.color <- palette}
  
  # Custom nodes: transparent insides and dark borders
  ega.plot <- ega.plot + 
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                        color = border.color,
                        shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
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
  
  # Add legend and plot title
  if(title != ""){
    ega.plot <- ega.plot +
      ggplot2::labs(title = title)
  }
  
  if(exists("legend.title")){
    ega.plot <- ega.plot +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette),
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ), title = legend.title)
      ) +
      ggplot2::theme(
        legend.title = ggplot2::element_text()
      )
  }
  
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot invariance
# Updated 07.07.2022
#' @export
plot.invariance <- function(
    x, title = "", labels = NULL,
    rows, columns,
    plot.args = list(), ...
)
{
  # Obtain structure
  structure <- x$memberships
  
  # Prepare EGA results for plots
  input_EGA <- lapply(x$groups$EGA, function(x){
    
    # Make class 'EGA'
    class(x) <- "EGA"
    
    # Return list
    return(x)
    
  })
  
  # Set structure
  input_EGA <- lapply(input_EGA, function(x){
    
    # Set structure
    x$wc <- structure
    
    # Return list
    return(x)
    
  })
  
  # Check for labels
  if(is.null(labels)){
    labels <- names(input_EGA)
  }
  
  # Check for rows
  if(missing(rows)){
    rows <- 1
  }
  
  # Check for columns
  if(missing(columns)){
    columns <- length(input_EGA)
  }
  
  # Set up plot arguments
  if(any(x$results$p <= .05)){
    
    # Check for "alpha" in plot.args
    if(!"alpha" %in% names(plot.args)){
      plot.args$alpha <- ifelse(
        x$results$p <= .05, .8, .2
      )
    }
    
  }
  
  # Obtain plots
  plots <- compare_EGA_plots(
    input.list = input_EGA,
    labels = labels, rows = rows,
    columns = columns,
    plot.args = plot.args
  )
  
}


#Plot CFA
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}

