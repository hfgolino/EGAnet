#' Descriptive Statistics for Networks
#'
#' Computes descriptive statistics for network models
#'
#' @param network Matrix, data frame,
#' \code{\link[qgraph]{qgraph}}, or \code{\link[EGAnet]{EGA}} object
#'
#' @return Numeric vector including:
#' 
#' \item{Mean_weight}{The average of the edge weights in the network}
#' 
#' \item{SD_weight}{The standard deviation of the edge weights in the network}
#' 
#' \item{Min_weight}{The minimum of the edge weights in the network}
#' 
#' \item{Max_weight}{The minimum of the edge weights in the network}
#' 
#' \item{Density}{The density of the network}
#' 
#' \item{ASPL}{The average shortest path length (ASPL) of the network (computed as unweighted)}
#' 
#' \item{CC}{The clustering coefficent (CC) of the network (computed as unweighted)}
#' 
#' \item{swn.rand}{Small-worldness measure based on random networks:
#' 
#' \deqn{swn.rand = (ASPL / ASPL_random) / (CC / CC_random)}
#' 
#' \code{swn.rand} > 1 suggests the network is small-world}
#' 
#' \item{swn.HG}{Small-worldness measure based on Humphries & Gurney (2008):
#' 
#' \deqn{swn.HG = (transitivity / transitivity_random) / (ASPL / ASPL_random)}
#' 
#' \code{swn.HG} > 1 suggests the network is small-world}
#' 
#' \item{swn.TJHBL}{Small-worldness measure based on Telesford, Joyce, Hayasaka, Burdette, & Laurienti (2011):
#' 
#' \deqn{swn.TJHBL = (ASPL_random / ASPL) - (CC / CC_lattice)}
#' 
#' \code{swn.TJHBL} near 0 suggests the network is small-world,
#' positive values suggest more random network characteristics,
#' negative values suggest more lattice network characteristics}
#' 
#' \item{scale-free_R-sq}{The R-squared fit of whether the degree distribution
#' follows the power-law (many small degrees, few large degrees)}
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
  names(desc) <- c("Mean_weight", "SD_weight", "Min_weight", "Max_weight",
                   "Density", "ASPL", "CC",
                   "swn.rand", "swn.HG", "swn.TJHBL", "scale-free_R-sq")
  
  # Collect descriptives
  connectivity <- NetworkToolbox::conn(network)
  degree <- NetworkToolbox::degree(network)
  desc["Mean_weight"] <- connectivity$mean
  desc["SD_weight"] <- connectivity$sd
  desc["Min_weight"] <- min(connectivity$weights)
  desc["Max_weight"] <- max(connectivity$weights)
  desc["Density"] <- connectivity$density
  desc["ASPL"] <- NetworkToolbox::pathlengths(network)$ASPL
  desc["CC"] <- NetworkToolbox::clustcoeff(network)$CC
  desc["swn.rand"] <- NetworkToolbox::smallworldness(network[degree != 0,degree != 0], method = "rand")$swm
  desc["swn.HG"] <- NetworkToolbox::smallworldness(network[degree != 0,degree != 0], method = "HG")$swm
  desc["swn.TJHBL"] <- NetworkToolbox::smallworldness(network[degree != 0,degree != 0], method = "TJHBL")$swm
  desc["scale-free_R-sq"] <- as.matrix(EGAnet:::scaleFreeFitIndex(degree, nBreaks = 10)["Rsquared.SFT"])
  
  return(round(desc, 3))
    
}
