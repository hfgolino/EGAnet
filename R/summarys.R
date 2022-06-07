#-------------------------------------------
## S3Methods summary() // Updated 13.05.2022
#-------------------------------------------

#' S3Methods for summarying
#' 
#' @name summarys
#'
#' @aliases 
#' summary.dynEGA
#' summary.dynEGA.Groups
#' summary.dynEGA.Individuals
#' summary.EGA
#' summary.NetLoads
#' summary.invariance
#' summary.hierEGA
#' summary.riEGA
#' 
#' @usage
#' \method{summary}{dynEGA}(x,  ...)
#' 
#' \method{summary}{dynEGA.Groups}(x, ...)
#' 
#' \method{summary}{dynEGA.Individuals}(x, ...)
#' 
#' \method{summary}{EGA}(x,  ...)
#' 
#' \method{summary}{NetLoads}(x, ...)
#' 
#' \method{summary}{invariance}(x, ...)
#' 
#' \method{summary}{hierEGA}(x, ...)
#' 
#' \method{summary}{riEGA}(x, ...)
#' 
#' @description summarys for \code{EGAnet} objects
#' 
#' @param x Object from \code{EGAnet} package
#' 
#' @param ... Additional arguments
#' 
#' @return summarys \code{EGAnet} object
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics par
#' 
# summary dynEGA (Level: Population)----
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

# summary dynEGA (Level: Groups)----
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

# summary dynEGA (Level: Individuals)----
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

# summary EGA----
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

#summary Network Loadings----
# Updated 13.05.2022
#' @export
summary.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$minLoad, arr.ind = TRUE)] <- ""
  
  summary(x$std)
  message("Loadings <= |", x$minLoad, "| are blank")
}

#summary Measurement Invariance----
# Updated 10.02.2022
#' @export
summary.invariance <- function(x, ...) {
  summary(x$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#summary Hierarchical EGA----
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

#summary Residual EGA----
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
