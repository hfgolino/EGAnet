#' Time-delay Embedding
#'
#' @description Reorganizes an individual’s observed time
#' series into an embedded matrix. The embedded matrix is constructed with replicates of an
#' individual time series that are offset from each other in time. The function requires
#' two parameters, one that specifies the number of observations to be used (i.e. the number of
#' embedded dimensions) and the other that specifies the number of observations to offset successive embeddings.
#'
#' @param x Vector.
#' An observed time series to be reorganized into a time-delayed embedded matrix.
#'
#' @param E Integer.
#' Number of embedded dimensions or the number of observations to be used. For example,
#' an \code{"E = 5"} will generate a matrix with five columns, meaning that five consecutive observations are used to create each row of the embedded matrix.
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings. A tau of one uses adjacent observations.
#' Default is \code{"tau = 1"}.
#'
#' @return Returns a matrix containing the embedded matrix.
#'
#' @examples
#' # A time series with 8 time points
#' tseries <- 49:56
#' embed.tseries <- Embed(tseries, E = 4, tau = 1)
#'
#'
#' @references
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14}, 367-386.
#'
#' @author Pascal Deboeck <pascal.deboeck at psych.utah.edu>
#'
#' @export
# Embed
# Updated 06.07.2023
Embed <- function(x, E, tau) {
  
  # Pre-compute series length - E * tau
  sl_E_tau <- length(x) - E * tau

  # Sequence along time series
  return(
    nvapply(
      seq_len(E), function(i){
        x[(1 + (i - 1) * tau):(sl_E_tau + i * tau)]
      }, LENGTH = sl_E_tau + tau
    )
  )
  
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
