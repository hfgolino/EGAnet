#' Computes the Signed Louvain Community Detection Algorithm
#' 
#' This function implementations the Louvain algorithm (Blondel et al., 2008)
#' that attempts to maximize signed modularity (Gomez et al., 2009)
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
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
# Updated 25.06.2023
signed.louvain <- function(network)
{
  
  # Ensure data is a matrix
  network <- as.matrix(network)
  
  # Call from C
  output <- .Call(
    "r_signed_louvain", network, PACKAGE = "EGAnet"
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
  output <- output[c("membership", "memberships", "modularity")]

  # Return
  return(output)
  
  
}
