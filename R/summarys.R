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
#' \method{summary}{dynEGA}(object,  ...)
#' 
#' \method{summary}{dynEGA.Groups}(object, ...)
#' 
#' \method{summary}{dynEGA.Individuals}(object, ...)
#' 
#' \method{summary}{EGA}(object,  ...)
#' 
#' \method{summary}{NetLoads}(object, ...)
#' 
#' \method{summary}{invariance}(object, ...)
#' 
#' \method{summary}{hierEGA}(object, ...)
#' 
#' \method{summary}{riEGA}(object, ...)
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
#' 
# summary dynEGA (Level: Population)----
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

# summary dynEGA (Level: Groups)----
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

# summary dynEGA (Level: Individuals)----
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

# summary EGA----
# Updated 13.05.2022
#' @export
summary.EGA <- function(object, ...) {
  
  # summary lower order communities
  cat(paste(
    "Number of communities:",
    object$n.dim,
    "\n\n"
  ))
  print(object$wc)
  
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
    object$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    object$Methods$corr
  )
  methods.matrix["Model =",] <- object$Methods$model
  methods.matrix["Algorithm =",] <- object$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    object$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  print(methods.matrix, quote = FALSE)
  
}

#summary Network Loadings----
# Updated 13.05.2022
#' @export
summary.NetLoads <- function(object, ...) {
  
  object$std[which(abs(object$std) <= object$minLoad, arr.ind = TRUE)] <- ""
  
  print(object$std)
  message("Loadings <= |", object$minLoad, "| are blank")
}

#summary Measurement Invariance----
# Updated 10.02.2022
#' @export
summary.invariance <- function(object, ...) {
  print(object$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}

#summary Hierarchical EGA----
# Updated 09.05.2022
#' @export
summary.hierEGA <- function(object, ...) {
  
  # summary lower order communities
  cat(paste(
    "Lower order communities:",
    object$hierarchical$lower_order$n.dim,
    "\n\n"
  ))
  print(object$hierarchical$lower_order$wc)
  
  # summary higher order communities
  cat(
    paste(
      "\nHigher order communities:",
      object$hierarchical$higher_order$EGA$n.dim,
      "\n\n"
    )
  )
  print(object$hierarchical$higher_order$EGA$wc)
  
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
    object$hierarchical$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    object$hierarchical$Methods$corr
  )
  methods.matrix["Model =",] <- object$hierarchical$Methods$model
  methods.matrix["Algorithm =",] <- object$hierarchical$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    object$hierarchical$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  methods.matrix["Scores =",] <- object$hierarchical$Methods$scores
  methods.matrix["Consensus Method =",] <- gsub(
    "_", " ", object$hierarchical$Methods$consensus.method
  )
  methods.matrix["Consensus Iterations =",] <- object$hierarchical$Methods$consensus.iter
  
  print(methods.matrix, quote = FALSE)
  
}

#summary Residual EGA----
# Updated 13.05.2022
#' @export
summary.riEGA <- function(object, ...) {
  
  # summary lower order communities
  cat(paste(
    "Number of communities:",
    object$EGA$n.dim,
    "\n\n"
  ))
  print(object$EGA$wc)
  
  # summary loadings if RI was necessary
  if("RI" %in% names(object)){
    
    ## Loadings
    ri_loadings <- round(as.vector(object$RI$loadings), 3)
    names(ri_loadings) <- row.names(object$RI$loadings)
    
    ## summary loadings
    cat("\nRandom-intercept loadings:\n\n")
    print(ri_loadings)
    
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
    object$Methods$corr == "cor_auto",
    "auto (from qgraph)",
    object$Methods$corr
  )
  methods.matrix["Model =",] <- object$Methods$model
  methods.matrix["Algorithm =",] <- object$Methods$algorithm
  methods.matrix["Unidimensional Method =",] <- ifelse(
    object$Methods$uni.method == "LE",
    "leading eigenvalue",
    "expand correlation matrix"
  )
  
  print(methods.matrix, quote = FALSE)
  
}
