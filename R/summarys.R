#---------------------------------------------
## S3Methods summary() // Updated 25.03.2020
#---------------------------------------------

#' S3Methods for Summaries
#' 
#' @name summarys
#'
#' @aliases 
#' summary.CFA
#' summary.dynEGA
#' summary.dynEGA.Groups
#' summary.dynEGA.Individuals
#' summary.EGA
#' summary.NetLoads
#' 
#' @usage
#' \method{summary}{CFA}(object, ...)
#' 
#' \method{summary}{dynEGA}(object, ...)
#' 
#' \method{summary}{dynEGA.Groups}(object, ...)
#' 
#' \method{summary}{dynEGA.Individuals}(object, ...)
#' 
#' \method{summary}{EGA}(object, ...)
#' 
#' \method{summary}{NetLoads}(object, ...)
#' 
#' @description Summaries for \code{EGAnet} objects
#' 
#' @param object Object from \code{EGAnet} package
#' 
#' @param ... Additional arguments
#' 
#' @return Summarizes \code{EGAnet} object
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#' 
# Summary CFA----
summary.CFA <- function(object, ...) {
  cat("Summary: Confirmatory Factor Analysis:\n")
  print(object$summary)
  cat("\n FIt Measures:\n")
  print(object$fit.measures)
}

# Summary dynEGA (Level: Population)----
summary.dynEGA <- function(object, ...) {
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(object$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dynEGA$dim.variables)
}

# Summary dynEGA (Level: Group)----
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

# Summary dynEGA (Level: Individual)----
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

# Summary EGA----
summary.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}

# Summary Network Loadings----
summary.NetLoads <- function(object, ...) {
  
  object$std[which(abs(object$std) <= object$MinLoad, arr.ind = TRUE)] <- ""
  
  print(object$std)
  message("Loadings <= ", object$MinLoad, " are blank")
}