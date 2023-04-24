#' Computes the signed Louvain community detection algorithm
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
#' 
#' @return Returns a list:
#' 
#' \item{communities}{Multi-level community matrix}
#' 
#' \item{modularities}{Modularity values for each level}
#' 
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#' 
#' # Estimate signed Louvain
#' signed_louvain(ega.wmt$network)
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Obtain signed louvain communities
# Updated 24.04.2023
signed_louvain <- function(network)
{
  
  # Ensure data is an double matrix
  network <- apply(as.matrix(network), 2, as.double)
  
  # Call from C
  output <- .Call(
    "r_signed_louvain",
    network,
    PACKAGE = "EGAnet"
  )
  
  # Check for variable names
  if(!is.null(colnames(network))){
    
    # Add names to output
    row.names(output$communities) <- colnames(network)
    names(output$modularities) <- 1:ncol(output$communities)
    
    # Transpose communities
    output$communities <- t(output$communities)
    
  }

  # Return
  return(output)
  
  
}
