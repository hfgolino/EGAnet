#' @title Frobenius Norm (Similarity)
#'
#' @description Computes the Frobenius Norm (Ulitzsch et al., 2023)
#'
#' @param network1 Matrix or data frame.
#' Network to be compared
#'
#' @param network2 Matrix or data frame.
#' Second network to be compared
#'
#' @examples
#' # Obtain wmt2 data
#' wmt <- wmt2[,7:24]
#'
#' # Set seed (for reproducibility)
#' set.seed(1234)
#'
#' # Split data
#' split1 <- sample(
#'   1:nrow(wmt), floor(nrow(wmt) / 2)
#' )
#' split2 <- setdiff(1:nrow(wmt), split1)
#'
#' # Obtain split data
#' data1 <- wmt[split1,]
#' data2 <- wmt[split2,]
#'
#' # Perform EBICglasso
#' glas1 <- EBICglasso.qgraph(data1)
#' glas2 <- EBICglasso.qgraph(data2)
#'
#' # Frobenius norm
#' frobenius(glas1, glas2)
#' # 0.7070395
#'
#' @return Returns Frobenius Norm
#'
#' @references
#' \strong{Simulation Study} \cr
#' Ulitzsch, E., Khanna, S., Rhemtulla, M., & Domingue, B. W. (2023).
#' A graph theory based similarity metric enables comparison of subpopulation psychometric networks
#' \emph{Psychological Methods}.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
#'
# Frobenius Norm
# Updated 10.07.2024
frobenius <- function(network1, network2)
{

  # Argument errors (send back networks in case of tibble)
  error_return <- frobenius_errors(network1, network2)

  # Return similarity
  return(
    1 / (
      1 +
      sqrt(sum((error_return$network1 - error_return$network2)^2, na.rm = TRUE)) /
      sqrt(dim(error_return$network1)[2] / 2)
    )
  )

}

#' @noRd
# Argument errors ----
# Updated 13.08.2023
frobenius_errors <- function(network1, network2)
{

  # 'network1' errors
  object_error(network1, c("matrix", "data.frame", "tibble"), "jsd")

  # Check for tibble
  if(!is(network1, "matrix")){
    network1 <- as.matrix(network1)
  }

  # 'network2' errors
  object_error(network2, c("matrix", "data.frame"), "jsd")

  # Check for tibble
  if(!is(network2, "matrix")){
    network2 <- as.matrix(network2)
  }

  # Return networks
  return(list(network1 = network1, network2 = network2))

}