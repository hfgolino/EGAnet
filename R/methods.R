#%%%%%%%%%%%%%%%%%%%%%%%%#
#### EGAnet S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%#

# Updated 07.06.2022 (DD.MM.YYYY)

# print() Methods ----

#Print Network Loadings
# Updated 13.05.2022
#' @export
print.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$minLoad, arr.ind = TRUE)] <- ""
  
  print(x$std)
  message("Loadings <= |", x$minLoad, "| are blank")
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

#summary Network Loadings
# Updated 13.05.2022
#' @export
summary.NetLoads <- function(object, ...) {
  
  object$std[which(abs(object$std) <= object$minLoad, arr.ind = TRUE)] <- ""
  
  print(object$std)
  message("Loadings <= |", object$minLoad, "| are blank")
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

#Plot CFA
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}

