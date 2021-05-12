#' Descriptive Statistics for Networks
#'
#' Computes descriptive statistics for network models
#'
#' @param network Matrix, data frame,
#' \code{\link[qgraph]{qgraph}}, or \code{\link[EGAnet]{EGA}} object
#'
#' @return 
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \donttest{# EGA example
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' ega.wmt <- EGA(data = wmt, plot.type = "qgraph")
#' }
#' 
#' # Compute descriptives
#' network.descriptives(ega.wmt)
#'
#' @references
#' # swn.HG \cr
#' Humphries, M. D., & Gurney, K. (2008).
#' Network 'small-world-ness': A quantitative method for determining canonical network equivalence.
#' \emph{PLoS one}, \emph{3}, e0002051
#' 
#' # swn.TJHBL \cr
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367-375
#' 
#' # scale-free_R-sq \cr
#' Langfelder, P., & Horvath, S. (2008).
#' WGCNA: an R package for weighted correlation network analysis.
#' \emph{BMC Bioinformatics}, \emph{9}, 559
#'
#' @export
#' 
#' @importFrom methods is
#'
# Network Descriptives
# Updated 12.05.2021
network.descriptives <- function(network)
{
  # Check for input
  if(isTRUE(methods::is(network, "qgraph"))){
    network <- qgraph::getWmat(network)
  }else if(isTRUE(methods::is(network, "EGA"))){
    network <- network$network
  }
  
  # Initialize descriptives matrix
  desc <- numeric(11)
  names(desc) <- c("Mean_pcor", "SD_pcor", "Min_pcor", "Max_pcor",
                   "Density", "ASPL", "CC",
                   "swn.rand", "swn.HG", "swn.TJHBL", "scale-free_R-sq")
  
  # Collect descriptives
  connectivity <- NetworkToolbox::conn(network)
  degree <- NetworkToolbox::degree(network)
  desc["Mean_pcor"] <- connectivity$mean
  desc["SD_pcor"] <- connectivity$sd
  desc["Min_pcor"] <- min(connectivity$weights)
  desc["Max_pcor"] <- max(connectivity$weights)
  desc["Density"] <- connectivity$density
  desc["ASPL"] <- NetworkToolbox::pathlengths(network)$ASPL
  desc["CC"] <- NetworkToolbox::clustcoeff(network)$CC
  desc["swn.rand"] <- NetworkToolbox::smallworldness(network, method = "rand")$swm
  desc["swn.HG"] <- NetworkToolbox::smallworldness(network, method = "HG")$swm
  desc["swn.TJHBL"] <- NetworkToolbox::smallworldness(network, method = "TJHBL")$swm
  desc["scale-free_R-sq"] <- scaleFreeFitIndex(degree, nBreaks = 10)
  
  return(round(desc, 3))
    
}