#-------------------------------------------
## S3Methods print() // Updated 25.03.2020
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
# Updated 02.05.2020
#' @export
print.dynEGA<- function(x, ...) {
  cat("dynEGA Results (Level: Population):\n")
  cat("\nNumber of Dimensions:\n")
  print(x$dynEGA$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dynEGA$dim.variables)
}

# Print dynEGA (Level: Groups)----
# Updated 02.05.2020
#' @export
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

# Print dynEGA (Level: Individuals)----
# Updated 02.05.2020
#' @export
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

# Print EGA----
# Updated 02.05.2020
#' @export
print.EGA <- function(x, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(x$n.dim)
  cat("\nItems per Dimension:\n")
  print(x$dim.variables)
}

#Print Network Loadings----
# Updated 02.05.2020
#' @export
print.NetLoads <- function(x, ...) {
  
  x$std[which(abs(x$std) <= x$MinLoad, arr.ind = TRUE)] <- ""
  
  print(x$std)
  message("Loadings <= ", x$MinLoad, " are blank")
}

#Print Measurement Invariance----
# Updated 10.02.2022
#' @export
print.invariance <- function(x, ...) {
  print(x$results, row.names = FALSE)
  cat("---\n")
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
}
