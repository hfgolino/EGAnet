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
# Updated 15.06.2023
signed.louvain <- function(network)
{
  
  # Ensure data is a matrix
  network <- as.matrix(network)
  
  # Call from C
  output <- .Call(
    "r_signed_louvain", network, PACKAGE = "EGAnet"
  )

  # Check for variable names
  if(!is.null(colnames(network))){
    
    # Add names to output
    colnames(output$memberships) <- colnames(network)
    names(output$modularity) <- nrow_sequence(output$memberships)
    
  }
  
  # Add highest level
  output$membership <- output$memberships[nrow(output$memberships),]
  
  # Reorder output
  output <- output[c("membership", "memberships", "modularity")]

  # Return
  return(output)
  
  
}
