#' @title Computes the Signed Louvain Community Detection Algorithm
#' 
#' @description This function implementations the Louvain algorithm (Blondel et al., 2008)
#' that attempts to maximize signed modularity (Gomez et al., 2009). When there
#' are no negative values in the network or absolute values are used, this function
#' is equivalent to \code{\link[igraph]{cluster_louvain}}
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
#' 
#' @param resolution Numeric (length = 1).
#' 
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results
#' Set for reproducible results.
#' See https://github.com/hfgolino/EGAnet/wiki/Reproducibility-and-PRNG
#' for more details on random number generation in \code{\link{EGAnet}}.
#' 
#' Louvain is a special case for reproducibility. It is standard practice
#' to shuffle the node order in the Louvain algorithm. By setting a seed,
#' node order will still be shuffled but the shuffle will be predictable
#' because of the seed.
#' 
#' Shuffling can happen multiple times in the Louvain algorithm meaning
#' that a single seed isn't enough; however, the current implementation
#' uses the current seed and increments the seed by \code{1} with each pass
#' 
#' @return Returns a list:
#' 
#' \item{membership}{Highest level of multi-level community matrix}
#' 
#' \item{memberships}{Multi-level community matrix}
#' 
#' \item{modularity}{Modularity values for each level}
#' 
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate network
#' network <- network.estimation(wmt, model = "glasso")
#' 
#' # Estimate signed Louvain
#' signed.louvain(network)
#' 
#' @references
#' \strong{Louvain algorithm} \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' \strong{Signed modularity} \cr
#' Gomez, S., Jensen, P., & Arenas, A. (2009).
#' Analysis of community structure in networks of correlated data.
#' \emph{Physical Review E}, \emph{80}(1), 016114.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Signed Louvain communities
# Updated 01.08.2023
signed.louvain <- function(network, resolution = 1, seed = NULL)
{
  
  # <= 85 nodes, faster or equal to {igraph}
  # > 85 nodes, {igraph} is faster
  
  # Argument errors
  signed.louvain_errors(network, resolution, seed)

  # Ensure data is a matrix
  network <- remove_attributes(as.matrix(network))
  
  # Call from C
  output <- .Call(
    "r_signed_louvain",
    network, resolution,
    swiftelse(is.null(seed), 0, seed),
    PACKAGE = "EGAnet"
  )
  
  # Get dimensions of memberships
  dimensions <- dim(output$memberships)
  
  # Get dimension names of network
  network_names <- dimnames(network)[[2]]

  # Check for variable names
  if(!is.null(network_names)){
    
    # Add names to output
    dimnames(output$memberships)[[2]] <- network_names
    names(output$modularity) <- seq_len(dimensions[1])
    
  }
  
  # Add highest level
  output$membership <- output$memberships[dimensions[1],]
  
  # Reorder output
  return(output[c("membership", "memberships", "modularity")])
  
  
}

#' @noRd
# Argument errors ----
# Updated 01.08.2023
signed.louvain_errors <- function(network, resolution, seed)
{
 
  # 'network' errors
  object_error(network, c("matrix", "data.frame"))
  
  # 'resolution' errors
  typeof_error(resolution, "numeric")
  length_error(resolution, 1)
  range_error(resolution, c(0, Inf))
  
  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1)
    typeof_error(seed, "numeric")
    range_error(seed,  c(0, as.double(.Machine$integer.max)))
  }
  
}
