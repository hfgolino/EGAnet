#' Generalized Local Linear Approximation
#'
#' @description Estimates the derivatives of a time series using generalized local linear approximation (GLLA).
#' GLLA is a filtering method for estimating derivatives from data that uses time delay embedding and a variant of Savitzky-Golay filtering to accomplish the task.
#'
#' @param x Vector.
#' An observed time series.
#'
#' @param n.embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function).
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings in the \code{\link[EGAnet]{Embed}} function. A tau of one uses adjacent observations.
#' Default is \code{"tau = 1"}.
#'
#' @param delta Integer.
#' The time between successive observations in the time series.
#' Default is \code{"delta = 1"}.
#'
#' @param order Integer.
#' The maximum order of the derivative to be estimated. For example,
#' \code{"order = 2"} will return a matrix with three columns with the estimates
#' of the observed scores and the first and second derivative for each row of the embedded
#' matrix (i.e. the reorganization of the time series implemented via
#' the \code{\link[EGAnet]{Embed}} function).
#'
#' @return Returns a matrix containing n columns, in which n is one plus the maximum order of the
#' derivatives to be estimated via generalized local linear approximation.
#'
#' @examples
#' # A time series with 8 time points
#' tseries <- 49:56
#' deriv.tseries <- glla(tseries, n.embed = 4, tau = 1, delta = 1, order = 2)
#'
#' @references
#'
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#'
#'
#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
# Generalized local linear approximation
# Updated 28.06.2023
glla <- function(x, n.embed, tau, delta, order){
  
  # Get time-delay embedding
  embedding <- Embed(x = x, E = n.embed, tau = tau)
  
  # Estimate derivatives
  ## Formally, `embedding %*% L %*% (solve(t(L) %*% L))`
  ## See helper function `glla_setup` below
  derivative_estimates <- embedding %*% glla_setup(
    n.embed = n.embed, tau = tau,
    delta = delta, order = order
  )
  
  # Add names
  if(order != 0){
    dimnames(derivative_estimates)[[2]] <- c(
      "Obs", paste0("DerivOrd", seq_len(order))
    )
  }else{
    dimnames(derivative_estimates)[[2]] <- "Obs"
  }
  
  # Return derivative estimates
  return(derivative_estimates)
  
}

# Bug checking ----
## Basic input
# x <- 49:56; n.embed = 4; tau = 1
# delta = 1; order = 2

# # Original function for comparison
# glla <- function(x, n.embed, tau, delta, order){
#   X <- Embed(x, E = n.embed, tau = tau)
#   
#   # Weights:
#   v <- 1:n.embed
#   mv <- mean(v)
#   order.der <- seq(from = 0, to = order)
#   L <- matrix(NA, nrow = n.embed, ncol = length(order.der))
#   for(i in 1:length(order.der)){
#     L[,i] <- ((tau*delta*v-mv)^order.der[[i]])/factorial(order.der[[i]])
#   }
#   
#   # Estimate the derivatives
#   Y <- X%*%L%*%solve(t(L)%*%L)
#   colnames(Y) <- c("Obs", paste0("DerivOrd", 1:order))
#   return(Y)
# }

#' @noRd
# GLLA Setup ----
# The purpose of this function is to avoid multiple repetitive
# computations of the same exact matrices when performing
# sample-wide dynamic EGA. By pre-computing the L matrix,
# embeddings can be multiplied by the same L matrix without
# the need to compute it for every single participant
# Updated 28.06.2023
glla_setup <- function(n.embed, tau, delta, order)
{
  
  # Set up weights
  embed_sequence <- seq_len(n.embed)
  
  # Get mean of embedding sequence
  mean_sequence <- mean(embed_sequence)
  
  # Derivative order
  derivative_order <- 0:order
  
  # Pre-compute tau * delta * embedding sequence - mean sequence
  embedding_value <- tau * delta * embed_sequence - mean_sequence
  
  # Get L matrix
  L <- nnapply(derivative_order, function(derivative){
    embedding_value^derivative / factorial(derivative)
  }, LENGTH = n.embed)
  
  # Return L matrix
  return(L %*% (solve(t(L) %*% L)))
  
}
