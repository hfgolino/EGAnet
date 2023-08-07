#' @title Computes the (Signed) Modularity Statistic
#'
#' @description Computes (signed) modularity statistic
#' given a network and community structure. Allows the
#' resolution parameter to be set 
#'
#' @param network Matrix or data frame.
#' A symmetric matrix representing a network
#' 
#' @param memberships Numeric (length = \code{ncol(network)}).
#' A numeric vector of integer values corresponding to 
#' each node's community membership
#' 
#' @param resolution Numeric (length = 1).
#' A parameter that adjusts modularity to
#' prefer smaller (\code{resolution} > 1) or larger
#' (0 < \code{resolution} < 1) communities.
#' Defaults to \code{1} (standard modularity computation)
#' 
#' @param signed Boolean (length = 1).
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
#' # 0.1697952
#' 
#' # Compute signed modularity
#' modularity(
#'   network = ega.wmt$network,
#'   memberships = ega.wmt$wc,
#'   signed = TRUE
#' )
#' # 0.1701946
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
# Updated 07.08.2023
modularity <- function(network, memberships, resolution = 1, signed = FALSE)
{
  
  # Argument errors (returns 'memberships' as a vector)
  memberships <- modularity_errors(
    network, memberships, resolution, signed
  )
  
  # Ensure data is a matrix
  network <- as.matrix(network)
  
  # Obtain dimensions
  dimensions <- dim(network)
  
  # Ensure names
  network <- ensure_dimension_names(network)

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
      ), call. = FALSE
    )
  }
  
  # Get network names
  network_names <- dimnames(network)[[2]]
  
  # Apply network names to memberships
  names(memberships) <- network_names
  
  # Check for absolute
  if(!signed){
    network <- abs(network)
  }
  
  # Obtain for missing memberships
  remove_nodes <- is.na(memberships)
  
  # Check for any missing
  if(any(remove_nodes)){
    
    # Set missing node names
    missing_nodes <- network_names[remove_nodes]
    
    # Push warning
    warning(
      paste0(
        "Nodes were missing values in 'memberships'. Modularity ",
        "was computed without these nodes: ",
        missing_nodes
      ), call. = FALSE
    )
    
  }
  
  # Call from C
  return(
    .Call(
      "r_signed_modularity",
      network[!remove_nodes, !remove_nodes], 
      as.integer(memberships[!remove_nodes]),
      resolution,
      PACKAGE = "EGAnet"
    )
  )
  
}

#' @noRd
# Argument errors ----
# Updated 04.08.2023
modularity_errors <- function(
    network, memberships, resolution, signed
)
{
  
  # 'network' errors
  object_error(network, c("matrix", "data.frame"))
  
  # 'memberships' errors
  object_error(memberships, c("vector", "matrix", "data.frame"))
  memberships <- force_vector(memberships)
  length_error(memberships, dim(network)[2])
  
  # 'resolution' errors
  length_error(resolution, 1)
  typeof_error(resolution, "numeric")
  range_error(resolution, c(0, Inf))
  
  # 'signed' errors
  length_error(signed, 1)
  typeof_error(signed, "logical")
  
  # Return memberships as a vector
  return(memberships)
  
}
  
