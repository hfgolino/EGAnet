#' @title Time-delay Embedding
#'
#' @description Reorganizes a single observed time series into an embedded matrix. The embedded 
#' matrix is constructed with replicates of an individual time series that are offset from 
#' each other in time. The function requires two parameters, one that specifies the number 
#' of observations to be used (i.e., the number of embedded dimensions) and the other that 
#' specifies the number of observations to offset successive embeddings
#'
#' @param x Numeric vector.
#' An observed time series to be reorganized into a time-delayed embedded matrix.
#'
#' @param E Numeric (length = 1).
#' Number of embedded dimensions or the number of observations to 
#' be used. \code{E = 5}, for example, will generate a matrix with 
#' five columns corresponding to five consecutive observations across
#' each row of the embedded matrix
#'
#' @param tau Numeric (length = 1).
#' Number of observations to offset successive embeddings. 
#' A tau of one uses adjacent observations.
#' Default is \code{tau = 1}
#'
#' @return Returns a numeric matrix
#'
#' @examples
#' # A time series with 8 time points
#' time_series <- 49:56
#' 
#' # Time series embedding
#' Embed(time_series, E = 5, tau = 1)
#'
#' @references
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14}, 367-386.
#'
#' @author Pascal Deboeck <pascal.deboeck at psych.utah.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# Embedding matrix ----
# Updated 31.07.2023
Embed <- function(x, E, tau) {
  
  # Pre-compute series length - E * tau
  # Uses argument errors to extract length
  # Either length of `x` is returned or
  # there is an argument error
  sl_E_tau <- Embed_errors(x, E, tau) - E * tau

  # Sequence along time series
  return(
    nvapply(
      seq_len(E), function(i){
        x[(1 + (i - 1) * tau):(sl_E_tau + i * tau)]
      }, LENGTH = sl_E_tau + tau
    )
  )
  
}

#' @noRd
# Argument errors ----
# Updated 31.07.2023
Embed_errors <- function(x, E, tau)
{
  
  # 'x' errors
  object_error(x, "vector")
  typeof_error(x, "numeric")
  
  # Get length of `x`
  x_length <- length(x)
  
  # 'E' errors
  length_error(E, 1)
  typeof_error(E, "numeric")
  range_error(E, c(1, x_length))
  
  # 'tau' errors
  length_error(tau, 1)
  typeof_error(tau, "numeric")
  range_error(tau, c(1, x_length))
  
  # Return 'x' length
  return(x_length)
  
}

# Bug checking ----
## Basic input
# x <- 49:56; E = 4; tau = 1

# # Original function for comparison
# Embed <- function(x,E,tau) {
#   len <- length(x)
#   out <- x[1:(len-(E*tau)+tau)]
#   for(i in 2:E) { out <- cbind(out,x[(1+((i-1)*tau)):(len-(E*tau)+(i*tau))]) }
#   return(out)
# }
