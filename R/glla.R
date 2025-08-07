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
#' @param na.derivative Character (length = 1).
#' How should missing data in the embeddings be handled?
#' Available options (see Boker et al. (2018) for more details):
#'
#' \itemize{
#'
#' \item \code{"none"} (default) --- does nothing and leaves \code{NA}s in data
#'
#' \item \code{"kalman"} --- uses Kalman smoothing (\code{\link[stats]{KalmanSmooth}}) with
#' structural time series models (\code{\link[stats]{StructTS}}) to impute missing values.
#' This approach models the underlying temporal dependencies (trend, seasonality, autocorrelation)
#' to generate estimates for missing observations while preserving the original time scale.
#' More computationally intensive than the other methods but typically provides the
#' most accurate imputation by respecting the stochastic properties of the time series
#'
#' \item \code{"rowwise"} --- adjusts time interval with respect to each embedding ensuring
#' time intervals are adaptive to the missing data (tends to be more accurate than \code{"none"})
#'
#' \item \code{"skipover"} --- "skips over" missing data and treats the non-missing points
#' as continuous points in time (note that the time scale shifts to the "per mean time interval,"
#' which is different and \emph{larger} than the original scale)
#'
#' }
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
#' \strong{Missing Data} \cr
#' Boker, S. M., Tiberio, S. S., & Moulder, R. G. (2018).
#' Robustness of time delay embedding to sampling interval misspecification.
#' In K. van Montfort, J. H. L. Oud, & M. C. Voelkle (Eds.),
#' \emph{Continuous Time Modeling in the Behavioral and Related Sciences} (pp. 239–258).
#' Springer International Publishing.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Generalized local linear approximation ----
# Updated 07.08.2025
glla <- function(x, n.embed, tau, delta, order, na.derivative = c("none", "kalman", "rowwise", "skipover"))
{

  # Set default
  na.derivative <- set_default(na.derivative, "none", glla)

  # Arguments errors
  glla_errors(x, n.embed, tau, delta, order)

  # Estimate derivatives
  ## Formally, `embedding %*% L %*% (solve(t(L) %*% L))`
  ## See helper function `glla_setup` below
  if(na.derivative != "none"){

    # Handle missing derivatives
    derivative_estimates <- switch(
      na.derivative,
      "kalman" = Embed(x = impute_kalman(x), E = n.embed, tau = tau) %*% glla_setup(n.embed, tau, delta, order),
      "skipover" = no_correction(x = x, n.embed = n.embed, tau = tau, delta = delta, order = order),
      "rowwise" = rowwise_correction(x = x, n.embed = n.embed, tau = tau, delta = delta, order = order)
    )

  }else{

    # Do derivatives
    derivative_estimates <- Embed(x = x, E = n.embed, tau = tau) %*%
                            glla_setup(n.embed, tau, delta, order)

  }

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
# Updated 13.08.2023
glla_errors <- function(x, n.embed, tau, delta, order)
{

  # 'x' errors
  object_error(x, c("vector", "matrix"), "glla")
  typeof_error(x, "numeric", "glla")

  # 'n.embed' errors
  typeof_error(n.embed, "numeric", "glla")
  length_error(n.embed, 1, "glla")

  # 'tau' errors
  typeof_error(tau, "numeric", "glla")
  length_error(tau, 1, "glla")

  # 'delta' errors
  typeof_error(delta, "numeric", "glla")
  length_error(delta, 1, "glla")

  # 'order' errors
  typeof_error(order, "numeric", "glla")
  length_error(order, 1, "glla")


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

#' @noRd
# Kalman Smoothing ----
# Updated 07.08.2025
impute_kalman <- function(x)
{

  # Make copy
  x_copy <- x

  # Missing indices
  index <- is.na(x)

  # If no missing values, return original
  if(!anyNA(x)){
    return(x)
  }

  # Check for first NA
  if(index[1]){
    x[1] <- x[which.min(is.na(x))]
  }

  # Fit StructTS model
  mod <- silent_call(try(StructTS(x)$model0, silent = TRUE))

  # If error, then just return the time series
  if(is(mod, "try-error")){
    return(x_copy)
  }

  # Extract smoothed values ONLY for missing observations
  x[index] <- KalmanSmooth(x, mod)$smooth[index,, drop = FALSE] %*% as.matrix(mod$Z)

  # Return imputed values
  return(x)

}

#' @noRd
# No correction ----
# Updated 06.08.2025
no_correction <- function(x, n.embed, tau, delta, order)
{

  # Time series length
  ts_length <- length(x)
  times <- seq_len(ts_length)

  # Set names
  names(x) <- paste0("V", seq_len(ts_length))

  # Rows and names
  rows <- (ts_length - n.embed) + 1

  # Set derivative matrix
  derivative_matrix <- matrix(
    NA, nrow = rows, ncol = order + 1,
    dimnames = list(names(x)[seq_len(rows)], NULL)
  )

  # Get embedding
  embedding <- Embed(na.omit(x), n.embed, tau) %*% glla_setup(n.embed, tau, delta, order)

  # Insert into derivative matrix
  derivative_matrix[row.names(embedding),] <- embedding

  # Return derivative matrix
  return(derivative_matrix)

}

#' @noRd
# Row-wise correction ----
# Updated 06.08.2025
rowwise_correction <- function(x, n.embed, tau, delta, order)
{

  # Get length of derivatives
  full_length <- length(x)

  # Set row sequence
  row_sequence <- seq_len(full_length)

  # Set row indices
  row_index <- Embed(row_sequence, n.embed, tau)

  # Check for NA indices
  available_index <- !is.na(x)

  # Set embedding
  embedding <- Embed(x[available_index], n.embed, tau)

  # Set time
  time <- Embed(seq_along(x)[available_index], n.embed, tau)

  # Match ending time to row index
  derivative_row <- match(time[,1], row_index[,1])

  # Set derivative matrix
  derivative_matrix <- matrix(NA, nrow = nrow(row_index), ncol = order + 1)

  # Compute intervals
  intervals <- apply(time, 1, function(x){sum(range(x)) / 2})

  # Get length of derivatives
  dx_length <- length(intervals)

  # Loop over to compute derivatives
  for(i in seq_len(dx_length)){

    # Compute loading matrix
    L <- matrix(1, nrow = n.embed, ncol = 1)

    # Add derivatives
    if(order == 1){
      L <- cbind(L, time[i,] - intervals[i])
    }else{

      # Compute center
      center <- time[i,] - intervals[i]

      # Add second order derivatives
      L <- cbind(L, center, center^2)

    }

    # Compute derivative
    derivative_matrix[derivative_row[i],] <- embedding[i,] %*% (L %*% solve(crossprod(L)))

  }

  # Return derivative matrix
  return(derivative_matrix)

}
