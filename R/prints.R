#-------------------------------------------
## S3Methods print() // Updated 13.05.2022
#-------------------------------------------

#' S3Methods for Printing
#' 
#' @name prints
#'
#' @aliases 
#' print.dynEGA
#' print.dynEGA.Groups
#' print.dynEGA.Individuals
#' print.EGA
#' print.NetLoads
#' print.invariance
#' print.hierEGA
#' 
#' @usage
#' \method{print}{dynEGA}(x,  ...)
#' 
#' \method{print}{dynEGA.Groups}(x, ...)
#' 
#' \method{print}{dynEGA.Individuals}(x, ...)
#' 
#' \method{print}{EGA}(x,  ...)
#' 
#' \method{print}{NetLoads}(x, ...)
#' 
#' \method{print}{invariance}(x, ...)
#' 
#' \method{print}{hierEGA}(x, ...)
#' 
#' @description Prints for \code{EGAnet} objects
#' 
#' @param x Object from \code{EGAnet} package
#' 
#' @param ... Additional arguments
#' 
#' @return Prints \code{EGAnet} object
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics par
#' 
# Print dynEGA (Level: Population)----
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

# Print dynEGA (Level: Groups)----
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

# Print dynEGA (Level: Individuals)----
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

# Print EGA----
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
    nrow = 3, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="# ,
    # "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  # methods.matrix["Unidimensional Method =",] <- ifelse(
  #   x$Methods$uni.method == "LE",
  #   "leading eigenvalue",
  #   "expand correlation matrix"
  # )
  
  print(methods.matrix, quote = FALSE)
  
}

#Print Network Loadings----
# Updated 13.05.2022
#' @export
print.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$minLoad, arr.ind = TRUE)] <- ""
  
  print(x$std)
  message("Loadings <= |", x$minLoad, "| are blank")
}

#Print Measurement Invariance----
# Updated 10.02.2022
#' @export
print.invariance <- function(x, ...) {
  print(x$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#Print Hierarchical EGA----
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
    nrow = 6, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm =",
    # "Unidimensional Method =",
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
  # methods.matrix["Unidimensional Method =",] <- ifelse(
  #     x$hierarchical$Methods$uni.method == "LE",
  #     "leading eigenvalue",
  #     "expand correlation matrix"
  #   )
  methods.matrix["Scores =",] <- x$hierarchical$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", x$hierarchical$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- x$hierarchical$Methods$consensus.iter
  
  print(methods.matrix, quote = FALSE)
  
}

#Print Residual EGA----
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
    nrow = 3, ncol = 1
  )
  row.names(methods.matrix) <- c(
    "Correlations =",
    "Model =",
    "Algorithm ="# ,
    # "Unidimensional Method ="
  )
  colnames(methods.matrix) <- ""
  
  methods.matrix["Correlations =",] <- ifelse(
    x$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    x$Methods$corr
  )
  methods.matrix["Model =",] <- x$Methods$model
  methods.matrix["Algorithm =",] <- x$Methods$algorithm
  # methods.matrix["Unidimensional Method =",] <- ifelse(
  #   x$Methods$uni.method == "LE",
  #   "leading eigenvalue",
  #   "expand correlation matrix"
  # )
  
  print(methods.matrix, quote = FALSE)
  
}

#Print Network Descriptives----
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
