#' @title Generalized Local Linear Approximation
#'
#' @description Estimates the derivatives of a time series using generalized 
#' local linear approximation (GLLA). GLLA is a filtering method for 
#' estimating derivatives from data that uses time delay embedding and a 
#' variant of Savitzky-Golay filtering to accomplish the task.
#'
#' @param x Numeric vector.
#' An observed time series
#'
#' @param n.embed Numeric (length = 1).
#' Number of embedded dimensions (the number of observations 
#' to be used in the \code{\link[EGAnet]{Embed}} function)
#'
#' @param tau Numeric (length = 1).
#' Number of observations to offset successive embeddings in 
#' the \code{\link[EGAnet]{Embed}} function. A \code{tau} of one 
#' uses adjacent observations.
#' Default is \code{1}
#'
#' @param delta Numeric (length = 1).
#' The time between successive observations in the time series.
#' Default is \code{1}
#'
#' @param order Numeric (length = 1).
#' The maximum order of the derivative to be estimated. For example,
#' \code{"order = 2"} will return a matrix with three columns with the estimates
#' of the observed scores and the first and second derivative for each row of the embedded
#' matrix (i.e. the reorganization of the time series implemented via
#' the \code{\link[EGAnet]{Embed}} function)
#'
#' @return Returns a matrix containing \emph{n} columns in which \emph{n}
#' is one plus the maximum order of the derivatives to be estimated via 
#' generalized local linear approximation
#'
#' @examples
#' # A time series with 8 time points
#' tseries <- 49:56
#' deriv.tseries <- glla(tseries, n.embed = 4, tau = 1, delta = 1, order = 2)
#'
#' @references
#' \strong{GLLA implementation} \cr
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#'
#' \strong{Filtering procedure} \cr
#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#' 
# Generalized local linear approximation ----
# Updated 03.08.2023
glla <- function(x, n.embed, tau, delta, order){
  
  # Arguments errors
  glla_errors(x, n.embed, tau, delta, order)
  
  # Estimate derivatives
  ## Formally, `embedding %*% L %*% (solve(t(L) %*% L))`
  ## See helper function `glla_setup` below
  derivative_estimates <- Embed(x = x, E = n.embed, tau = tau) %*% 
                          glla_setup(n.embed, tau, delta, order)
  
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

#' @noRd
# GLLA errors
# Updated 03.08.2023
glla_errors <- function(x, n.embed, tau, delta, order)
{
  
  # 'x' errors
  object_error(x, c("vector", "matrix"))
  typeof_error(x, "numeric")
  
  # 'n.embed' errors
  typeof_error(n.embed, "numeric")
  length_error(n.embed, 1)
  
  # 'tau' errors
  typeof_error(tau, "numeric")
  length_error(tau, 1)
  
  # 'delta' errors
  typeof_error(delta, "numeric")
  length_error(delta, 1)
  
  # 'order' errors
  typeof_error(order, "numeric")
  length_error(order, 1)
  
  
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
# the need to compute it for every single participant x variable
# Updated 07.07.2023
glla_setup <- function(n.embed, tau, delta, order)
{
  
  # Set up weights
  embed_sequence <- seq_len(n.embed)
  
  # Pre-compute tau * delta * embedding sequence - mean sequence
  embedding_value <- tau * delta * embed_sequence - mean(embed_sequence)
  
  # Get L matrix
  L <- nvapply(0:order, function(derivative) {
    embedding_value^derivative / gamma(derivative + 1)
    # `gamma` is the .Primitive for `factorial`
    # same as `factorial` under the hood
  }, LENGTH = n.embed)
  
  # Return L matrix
  # L %*% (solve(t(L) %*% L))
  return(L %*% solve(crossprod(L)))
  
}