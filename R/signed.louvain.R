#' Computes the Signed Louvain Community Detection Algorithm
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
#' 
#' @return Returns a list:
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
#' network <- EBICglasso.qgraph(data = wmt)
#' 
#' # Estimate signed Louvain
#' signed.louvain(network)
#' 
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
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
# Updated 09.06.2023
signed.louvain <- function(network)
{
  
  # Ensure data is a matrix
  network <- as.matrix(network)
  
  # Call from C
  output <- .Call(
    "r_signed_louvain",
    network,
    PACKAGE = "EGAnet"
  )

  # Check for variable names
  if(!is.null(colnames(network))){
    
    # Add names to output
    colnames(output$memberships) <- colnames(network)
    names(output$modularity) <- 1:nrow(output$memberships)
    
  }

  # Return
  return(output)
  
  
}
