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
# Updated 13.05.2022
#' @export
print.dynEGA.Individuals <- function(x, ...) {
  
  # Number of people
  cat("Number of Cases (individuals): ")
  number <- length(x$dynEGA) - 1
  cat(number, "\n")
  
  # Summary statistics
  cat("\nSummary statistics (number of communities): \n")
  dim <- sapply(x$dynEGA[-length(x$dynEGA)], "[[", 3)
  
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
  summary.methods["Mean =",] <- mean(dim)
  summary.methods["Median =",] <- median(dim)
  summary.methods["Min =",] <- min(dim)
  summary.methods["Max =",] <- max(dim)
  
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

# Print EGA
# Updated 13.05.2022
#' @export
print.EGA <- function(x, ...) {
  
  # Print lower order communities
  cat(paste(
    "Number of communities:",
    x$n.dim,
    "\n\n"
  ))
  print(x$wc)
  
  # Print methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  print(methods.matrix, quote = FALSE)
  
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
# Updated 09.05.2022
#' @export
print.hierEGA <- function(x, ...) {
  
  # Print lower order communities
  cat(paste(
    "Lower order communities:",
    x$hierarchical$lower_order$n.dim,
    "\n\n"
  ))
  print(x$hierarchical$lower_order$wc)
  
  # Print higher order communities
  cat(
    paste(
      "\nHigher order communities:",
      x$hierarchical$higher_order$EGA$n.dim,
      "\n\n"
    )
  )
  print(x$hierarchical$higher_order$EGA$wc)
  
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
    x$hierarchical$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$hierarchical$Methods$corr
  )
  methods.matrix["Model =",] <- x$hierarchical$Methods$model
  methods.matrix["Algorithm =",] <- x$hierarchical$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$hierarchical$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  methods.matrix["Scores =",] <- x$hierarchical$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", x$hierarchical$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- x$hierarchical$Methods$consensus.iter
  
  print(methods.matrix, quote = FALSE)
  
}

#Print Residual EGA
# Updated 13.05.2022
#' @export
print.riEGA <- function(x, ...) {
  
  # Print lower order communities
  cat(paste(
    "Number of communities:",
    x$EGA$n.dim,
    "\n\n"
  ))
  print(x$EGA$wc)
  
  # Print loadings if RI was necessary
  if("RI" %in% names(x)){
    
    ## Loadings
    ri_loadings <- round(as.vector(x$RI$loadings), 3)
    names(ri_loadings) <- row.names(x$RI$loadings)
    
    ## Print loadings
    cat("\nRandom-intercept loadings:\n\n")
    print(ri_loadings)
    
  }
  
  # Print methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  print(methods.matrix, quote = FALSE)
  
}

# summary() Methods ----

# summary dynEGA (Level: Population)
# Updated 13.05.2022
#' @export
summary.dynEGA<- function(x, ...) {
  
  # summary communities
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
  
  # summary item placements
  summary(item_placement)
  
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
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  summary(glla.methods, quote = FALSE)
  
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
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  summary(ega.methods, quote = FALSE)
  
}

# summary dynEGA (Level: Groups)
# Updated 13.05.2022
#' @export
summary.dynEGA.Groups <- function(x, ...) {
  
  for(i in 1:(length(x$dynEGA) - 1)){
    
    # summary communities
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
    
    # summary item placements
    summary(
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
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  summary(glla.methods, quote = FALSE)
  
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
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  summary(ega.methods, quote = FALSE)
  
}

# summary dynEGA (Level: Individuals)
# Updated 13.05.2022
#' @export
summary.dynEGA.Individuals <- function(x, ...) {
  
  # Number of people
  cat("Number of Cases (individuals): ")
  number <- length(x$dynEGA) - 1
  cat(number, "\n")
  
  # Summary statistics
  cat("\nSummary statistics (number of communities): \n")
  dim <- sapply(x$dynEGA[-length(x$dynEGA)], "[[", 3)
  
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
  summary.methods["Mean =",] <- mean(dim)
  summary.methods["Median =",] <- median(dim)
  summary.methods["Min =",] <- min(dim)
  summary.methods["Max =",] <- max(dim)
  
  ## summary summary
  summary(summary.methods, quote = FALSE)
  
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
  glla.methods["Embedding Dimensions =",] <- x$dynEGA$Methods$glla$n.embed
  glla.methods["Embedding Offset (tau) =",] <- x$dynEGA$Methods$glla$tau
  glla.methods["Lag (delta) =",] <- x$dynEGA$Methods$glla$delta
  
  ## summary GLLA
  summary(glla.methods, quote = FALSE)
  
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
    x$dynEGA$Methods$EGA$corr == "cor_auto",
    "auto (from qgraph)",
    x$dynEGA$Methods$EGA$corr
  )
  ega.methods["Model =",] <- x$dynEGA$Methods$EGA$model
  ega.methods["Algorithm =",] <- x$dynEGA$Methods$EGA$algorithm
  
  ## summary EGA
  summary(ega.methods, quote = FALSE)
  
}

# summary EGA
# Updated 13.05.2022
#' @export
summary.EGA <- function(x, ...) {
  
  # summary lower order communities
  cat(paste(
    "Number of communities:",
    x$n.dim,
    "\n\n"
  ))
  summary(x$wc)
  
  # summary methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  summary(methods.matrix, quote = FALSE)
  
}

#summary Network Loadings
# Updated 13.05.2022
#' @export
summary.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$minLoad, arr.ind = TRUE)] <- ""
  
  summary(x$std)
  message("Loadings <= |", x$minLoad, "| are blank")
}

#summary Measurement Invariance
# Updated 10.02.2022
#' @export
summary.invariance <- function(x, ...) {
  summary(x$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#summary Hierarchical EGA
# Updated 09.05.2022
#' @export
summary.hierEGA <- function(x, ...) {
  
  # summary lower order communities
  cat(paste(
    "Lower order communities:",
    x$hierarchical$lower_order$n.dim,
    "\n\n"
  ))
  summary(x$hierarchical$lower_order$wc)
  
  # summary higher order communities
  cat(
    paste(
      "\nHigher order communities:",
      x$hierarchical$higher_order$EGA$n.dim,
      "\n\n"
    )
  )
  summary(x$hierarchical$higher_order$EGA$wc)
  
  # summary methods
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
    x$hierarchical$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$hierarchical$Methods$corr
  )
  methods.matrix["Model =",] <- x$hierarchical$Methods$model
  methods.matrix["Algorithm =",] <- x$hierarchical$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$hierarchical$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  methods.matrix["Scores =",] <- x$hierarchical$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", x$hierarchical$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- x$hierarchical$Methods$consensus.iter
  
  summary(methods.matrix, quote = FALSE)
  
}

#summary Residual EGA
# Updated 13.05.2022
#' @export
summary.riEGA <- function(x, ...) {
  
  # summary lower order communities
  cat(paste(
    "Number of communities:",
    x$EGA$n.dim,
    "\n\n"
  ))
  summary(x$EGA$wc)
  
  # summary loadings if RI was necessary
  if("RI" %in% names(x)){
    
    ## Loadings
    ri_loadings <- round(as.vector(x$RI$loadings), 3)
    names(ri_loadings) <- row.names(x$RI$loadings)
    
    ## summary loadings
    cat("\nRandom-intercept loadings:\n\n")
    summary(ri_loadings)
    
  }
  
  # summary methods
  cat("\nMethods:\n")
  
  ## Set up methods
  methods.matrix <- matrix(
    nrow = 4, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    x$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  summary(methods.matrix, quote = FALSE)
  
}

# plot() Methods ----

# Plot bootEGA
# Updated 28.07.2021
plot.bootEGA <- function(x, plot.type = c("GGally","qgraph"),
                         plot.args = list(), produce = TRUE, ...){
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
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
    
  }
  
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot dynEGA function (Level: Group)
# Updated 28.07.2021
#' @export
plot.dynEGA.Groups <- function(x, ncol, nrow, title = "", plot.type = c("GGally","qgraph"),
                               plot.args = list(), produce = TRUE, ...){
  
  # Remove methods from input list
  if("Methods" %in% names(x$dynEGA)){
    x$dynEGA <- x$dynEGA[-which(names(x$dynEGA) == "Methods")]
  }
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
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
      network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
      network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
      network::set.edge.value(network1[[i]], attrname="AbsWeights",value=abs(x$dynEGA[[i]]$network))
      network::set.edge.value(network1[[i]],attrname="ScaledWeights",
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
      plot.args$net <- network1[[i]]
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
      
    }
    group.labels <- names(x$dynEGA)
    set.seed(NULL)
    if(isTRUE(produce)){
      ggpubr::ggarrange(plotlist=plots.net, ncol = ncol, nrow = nrow, labels = group.labels, label.x = 0.3)
    }else{return(plots.net)}
  }
}

# Plot dynEGA function (Level: Individual)
# Updated 28.07.2021
#' @export
plot.dynEGA.Individuals <- function(x, title = "",  id = NULL, plot.type = c("GGally","qgraph"),
                                    plot.args = list(), produce = TRUE, ...){
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
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
    
    set.seed(NULL)
    
    if(isTRUE(produce)){
      plot(ega.plot)
    }else{return(ega.plot)}
  }
}

# Plot dynEGA function (Level: Population)
# Updated 28.07.2021
#' @export
plot.dynEGA <- function(x, title = "", plot.type = c("GGally","qgraph"),
                        plot.args = list(), produce = TRUE, ...){
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
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
    
  }
  set.seed(NULL)
  
  if(isTRUE(produce)){
    plot(ega.plot)
  }else{return(ega.plot)}
}

# Plot EGA
# Updated 07.02.2022
#' @export
plot.EGA <- function(x,  title = "", plot.type = c("GGally","qgraph"),
                     plot.args = list(), produce = TRUE, ...){
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  ## Check for sum
  sum_nodes <- sum(x$network)
  
  if(sum_nodes != 0){
    
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
      network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
      network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
      network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$network))
      network::set.edge.value(network1,attrname="ScaledWeights",
                              value=matrix(rescale.edges(x$network, plot.args$edge.size),
                                           nrow = nrow(x$network),
                                           ncol = ncol(x$network)))
      
      # Layout "Spring"
      graph1 <- convert2igraph(x$network)
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
      plot.args$edge.lty <- "line"
      plot.args$edge.size <- "ScaledWeights"
      
      lower <- abs(x$network[lower.tri(x$network)])
      non.zero <- sqrt(lower[lower != 0])
      
      plot.args$edge.alpha <- non.zero
      plot.args$mode <- layout.spring
      plot.args$label <- colnames(x$network)
      plot.args$node.label <- rep("", ncol(x$network))
      if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
      if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
      
      palette <- color_palette_EGA(color.palette, as.numeric(factor(x$wc)))
      palette <- ifelse(is.na(palette), "white", palette)
      
      ega.plot <- suppressWarnings(
        suppressMessages(
          do.call(GGally::ggnet2, plot.args) + 
            ggplot2::theme(legend.title = ggplot2::element_blank())
        )
      )
      
      name <- colnames(x$network)
      
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
      
    }
    
    set.seed(NULL)
    
    if(isTRUE(produce)){
      plot(ega.plot)
    }else{return(ega.plot)}
    
  }else{
    message("Network is empty. No plot produced.")
  }
  
  
}

# Plot EGA.fit
# Updated 17.03.2021
#' @export
plot.EGA.fit <- function(x,  title = "", plot.type = c("GGally","qgraph"),
                         plot.args = list(), ...){
  
  # MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  plot.EGA(x$EGA, plot.type = plot.type, plot.args = plot.args)
}

#Plot net.loads
# Updated 02.05.2020
#' @export
plot.NetLoads <- function(x, ...) {
  
  plot(x$plot)
}

# Plot invariance
# Updated 10.02.2022
#' @export
plot.invariance <- function(
    x, title = "", labels = NULL,
    rows, columns, plot.type = c("GGally","qgraph"),
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
  
  # Check for plot type
  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{
    plot.type <- match.arg(plot.type)
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
  plots <- compare.EGA.plots(
    input_list = input_EGA,
    labels = labels, rows = rows,
    columns = columns, plot.type = plot.type,
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






