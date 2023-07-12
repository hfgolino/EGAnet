#-------------------------------------------
## S3Methods summary() // Updated 13.05.2022
#-------------------------------------------

#' S3Methods for summarying
#' 
#' @name summarys
#'
#' @aliases 
#' summary.hierEGA
#' 
#' @usage 
#' \method{summary}{hierEGA}(object, ...)
#' 
#' @description summarys for \code{EGAnet} objects
#' 
#' @param object Object from \code{EGAnet} package
#' 
#' @param ... Additional arguments
#' 
#' @return summarys \code{EGAnet} object
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics par
#summary Hierarchical EGA----
# Updated 16.11.2022
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

#Summary Network Descriptives----
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