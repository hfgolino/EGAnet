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
#'
#' # A time series with 8 time points
#' tseries <- 49:56
#' deriv.tseries <- glla(tseries, n.embed = 4, tau = 1, delta = 1, order = 2)
#'
#'
#' @references
#'
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#' \doi{10.1037/a0016622}
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#' \doi{10.1037/a0016622}
#'
#'
#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#' \doi{10.1021/ac60214a047}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
# glla
# Updated 22.05.2021
#'
glla <- function(x, n.embed, tau, delta, order){
  X <- Embed(x, E = n.embed, tau = tau)

  # Weights:
  v <- 1:n.embed
  mv <- mean(v)
  order.der <- seq(from = 0, to = order)
  L <- matrix(NA, nrow = n.embed, ncol = length(order.der))
  for(i in 1:length(order.der)){
    L[,i] <- ((tau*delta*v-mv)^order.der[[i]])/factorial(order.der[[i]])
  }

  # Estimate the derivatives
  Y <- X%*%L%*%solve(t(L)%*%L)
  colnames(Y) <- c("Obs", paste0("DerivOrd", 1:order))
  return(Y)
}
#----
