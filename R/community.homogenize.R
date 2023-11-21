#' @title Homogenize Community Memberships
#'
#' @description Memberships from community detection algorithms do not always
#' align numerically. This function seeks to homogenize 
#' community memberships between a target membership (the 
#' membership to homogenize toward) and one or more other 
#' memberships. This function is the core of the 
#' \code{\link[EGAnet]{dimensionStability}} and 
#' \code{\link[EGAnet]{itemStability}} functions
#'
#' @param target.membership Vector, matrix, or data frame.
#' The target memberships that all other memberships input into
#' \code{convert.membership} should be homogenize \strong{toward}
#'
#' @param convert.membership Vector, matrix, or data frame.
#' Either a vector of memberships the same length as
#' \code{target.membership} or a matrix or data frame of many
#' membership solutions with either across rows or down columns the same
#' length as \code{target.membership} (this function will automatically
#' determine this orientation for you with precedence given solutions
#' \emph{across rows})
#' 
#' @return Returns a vector or matrix the length or size of
#' \code{convert.membership} with memberships homogenized toward
#' \code{target.membership}
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Get network
#' network <- network.estimation(wmt2[,7:24])
#' 
#' # Apply Walktrap
#' network_walktrap <- community.detection(
#'   network, algorithm = "walktrap"
#' )
#' 
#' # Apply Louvain
#' network_louvain <- community.detection(
#'   network, algorithm = "louvain"
#' )
#' 
#' # Homogenize toward Walktrap
#' community.homogenize(network_walktrap, network_louvain)
#'
#' @references
#' \strong{Original implementation of bootEGA} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#' 
#' @export
#'
# Make memberships as homogeneous as possible ----
# Updated 19.11.2023
community.homogenize <- function(target.membership, convert.membership)
{
  
  # Send error if memberships are not a vector, matrix, or data frame
  object_error(target.membership, c("vector", "matrix", "data.frame"), "community.homogenize")
  object_error(convert.membership, c("vector", "matrix", "data.frame"), "community.homogenize")
  
  # Ensure target membership is a vector
  target.membership <- force_vector(target.membership)

  # Ensure conversion matrix is a case-by-row matrix
  convert.membership <- matrixize_conversion(
    convert.membership, length(target.membership)
  )
  
  # Set index ID
  numeric_ID <- nrow_sequence(convert.membership)

  # Add unique vector ID
  names(numeric_ID) <- vector2factor(convert.membership)
  
  # To speed up process, get unique solutions only
  unique_homogenized <- row_apply(
    convert.membership[numeric_ID[unique(names(numeric_ID))],, drop = FALSE],
    single_homogenize, target.membership = target.membership
  )
  
  # Return the ordered homogenized memberships
  return(unique_homogenized[as.numeric(names(numeric_ID)),, drop = FALSE])
  
}

#' @noRd
# Make conversion membership matrix ----
# Updated 06.07.2023
matrixize_conversion <- function(convert.membership, target_length)
{
  
  # Check for whether conversion membership is a vector
  if(get_object_type(convert.membership) == "vector"){
    
    # Make sure length equals target membership
    length_error(convert.membership, target_length, "community.homogenize")
    
    # Convert to matrix
    return(matrix(convert.membership, nrow = 1))
    
  }else{
    
    # Ensure matrix
    convert.membership <- as.matrix(convert.membership)
    
    # Get dimensions of conversion matrix
    dimensions <- dim(convert.membership)
    
    # Determine which dimensions align with target membership
    if(# edge case: square matrix, then assume case-by-row matrix already
      dimensions[1] == target_length &
      dimensions[2] != target_length
    ){ # Transpose and return
      return(t(convert.membership))
    }else{ # Otherwise, just return
      return(convert.membership)
    }
    
  }
  
}

#' @noRd
# Get Rand values ----
# Updated 23.07.2023
get_rand <- function(convert_keep, target_keep)
{
  
  # This branching is ugly... but it's faster
  # than computing the Rand index for every community
  # (about 50% faster than without it)
  return(
    nvapply(
      seq_len(max(convert_keep)), function(community){
  
        # Get indices
        community_index <- convert_keep == community
        
        # Total indices
        total_indices <- sum(community_index)
        
        # Get unique lengths
        target_length <- unique_length(target_keep[community_index])
        convert_length <- unique_length(convert_keep[community_index])
        
        # Branch for different cases before Rand index
        if(total_indices < 3){
          return(0) # Doublets go last
        }else if(
          target_length == total_indices || convert_length == total_indices
        ){ # All singleton, throw zero
          return(0) 
        }else if((target_length * convert_length) == 1){ 
          return(1) # Both all one value, then Rand index = 1
        }else{ 
          return( # Otherwise, return Rand index
            igraph::compare(
              target_keep[community_index],
              convert_keep[community_index],
              method = "rand"
            )
          )
        }
        
      }
    )
  )
  
}

#' @noRd
# Core Homogenize Function ----
# Updated 19.11.2023
single_homogenize <- function(target.membership, convert_single)
{

  # Handle NAs
  keep_memberships <- !(is.na(target.membership) | is.na(convert_single))
  
  # Get memberships to keep
  target_keep <- target.membership[keep_memberships]

  # Set up conversion as unique values
  convert_keep <- reindex_memberships(convert_single[keep_memberships])

  # Create table between target and conversion
  conversion_table <- table(convert_keep, target_keep)
  
  # Store conversion order
  convert_order <- order(
    get_rand(convert_keep, target_keep), # first order by Rand index
    rowSums(conversion_table), # then order by number of nodes in community
    decreasing = TRUE
  )
  
  # Order table by communities based first on Rand, then on frequencies
  # After, get maximum value for each conversion
  conversion_table <- max.col(conversion_table[convert_order,], "first")
  
  # Set duplicates
  conversion_duplicates <- duplicated(conversion_table)
  
  # For duplicates, replace with additional communities
  conversion_table[conversion_duplicates] <- 
    max(conversion_table) + seq_len(sum(conversion_duplicates))
  
  # Create key by adding names to conversion table
  names(conversion_table) <- convert_order
  
  # Replace values in original membership
  convert_single[keep_memberships] <- conversion_table[as.character(convert_keep)]

  # Get conversion based on conversion table
  return(convert_single)
  
}

