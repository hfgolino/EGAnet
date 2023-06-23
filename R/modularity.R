#' Computes the (Signed) Modularity Statistic
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
#' 
#' @param memberships Numeric (length = \code{ncol(network)}).
#' A numeric vector of integer values corresponding to 
#' each node's community membership
#' 
#' @param signed Boolean.
#' Whether signed or absolute modularity should be computed.
#' The most common modularity metric is defined by positive values only. 
#' Gomez et al. (2009) introduced a signed version of modularity that
#' will discount modularity for edges with negative values. This property
#' isn't always desired for psychometric networks. If \code{TRUE}, then
#' this signed modularity metric will be computed. If \code{FALSE}, then
#' the absolute value of the edges in the network (using \code{abs}) will
#' be used to compute modularity.
#' Defaults to \code{FALSE}
#'  
#' @return Returns the modularity statistic
#' 
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate EGA
#' ega.wmt <- EGA(wmt, model = "glasso")
#' 
#' # Compute standard (absolute values) modularity
#' modularity(
#'   network = ega.wmt$network,
#'   memberships = ega.wmt$wc,
#'   signed = FALSE
#' )
#' 
#' # Compute signed modularity
#' modularity(
#'   network = ega.wmt$network,
#'   memberships = ega.wmt$wc,
#'   signed = TRUE
#' )
#' 
#' @references
#' Gomez, S., Jensen, P., & Arenas, A. (2009).
#' Analysis of community structure in networks of correlated data.
#' \emph{Physical Review E}, \emph{80}(1), 016114.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Modularity statistic
# Updated 23.06.2023
modularity <- function(network, memberships, signed = FALSE)
{
  
  # Ensure data is a matrix
  network <- as.matrix(network)
  
  # Obtain dimensions
  dimensions <- dim(network)
  
  # Ensure names
  network <- ensure_dimension_names(network)
  
  # Ensure memberships is a vector
  memberships <- force_vector(memberships)
  
  # Membership length
  membership_length <- length(memberships)
  
  # Ensure membership length equals nodes
  if(dimensions[2] != membership_length){
    stop(
      paste0(
        "Number of nodes in the 'network' (`ncol` = ",
        dimensions[2],
        ") does not equal the number of nodes in the 'memberships' (`length` = ",
        membership_length, ")"
      )
    )
  }
  
  # Apply network names to memberships
  names(memberships) <- colnames(network)
  
  # Check for absolute
  if(isFALSE(signed)){
    network <- abs(network)
  }
  
  # Obtain for missing memberships
  remove_nodes <- is.na(memberships)
  
  # Set network and memberships
  network <- network[!remove_nodes, !remove_nodes]
  memberships <- as.integer(memberships[!remove_nodes])
  
  # Check for any missing
  if(any(remove_nodes)){
    
    # Set missing node names
    missing_nodes <- colnames(network)[remove_nodes]
    
    # Push warning
    warning(
      paste0(
        "Nodes were missing values in 'memberships'. Modularity ",
        "was computed without these nodes: ",
        missing_nodes
      )
    )
    
  }
  
  # Call from C
  output <- .Call(
    "r_signed_modularity",
    network, memberships,
    PACKAGE = "EGAnet"
  )

  # Return
  return(output)
  
  
}