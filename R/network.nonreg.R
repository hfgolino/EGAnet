#' @title Iterative Threshold Non-regularized Network Estimation
#'
#' @description A relatively simple approach to non-regularized network estimation.
#' Works by incrementally thresholding values in the inverse covariance (precision) matrix
#' and computing some criterion. After the initial pass, the selected precision matrix
#' undergoes further fitting using \code{"L-BFGS-B"} optimization to minimize the criterion
#' (uses \code{\link[stats]{optim}})
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param n Numeric (length = 1).
#' Sample size \strong{must} be provided if \code{data} provided is a correlation matrix
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)
#'
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
#' Arguments can be passed along to the function
#'
#' \item \code{"cosine"} --- Uses \code{\link[EGAnet]{cosine}} to compute cosine similarity
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
#'
#' }
#'
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#'
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param criterion Character (length = 1).
#' Which criterion should be minimized?
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"AIC"} --- Akaike information criterion
#'
#' \item \code{"BIC"} --- Bayesian information criterion
#'
#' \item \code{"EBIC"} --- extended Bayesian information criterion
#'
#' }
#'
#' Defaults to \code{"BIC"}
#'
#' @param gamma Numeric (length = 1)
#' EBIC tuning parameter.
#' Defaults to \code{0.50} and is generally a good choice.
#' Setting to \code{0} will cause regular BIC to be used
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to \code{\link[EGAnet]{auto.correlate}}
#'
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and
#' Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain default (BIC)
#' nonreg_graph <- network.nonreg(data = wmt)
#'
#' # Obtain AIC
#' nonreg_aic_graph <- network.nonreg(data = wmt, criterion = "AIC")
#'
#' # Obtain EBIC
#' nonreg_ebic_graph <- network.nonreg(data = wmt, criterion = "EBIC")
#'
#' @noRd
#'
# Non-regularization using optimization ----
# Updated 28.12.2024
network.nonreg <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "cosine", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    criterion = c("AIC", "BIC", "EBIC"),
    gamma = 0.50, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network.nonreg)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  criterion <- set_default(criterion, "bic", network.nonreg)

  # Argument errors (return data in case of tibble)
  data <- nonreg_errors(data, n, gamma, verbose, ...)

  # Get dimensions of the data
  dimensions <- dim(data)

  # Get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n,
    corr = corr, na.data = na.data,
    verbose = verbose, needs_usable = FALSE, # skips usable data check
    ...
  )

  # Get outputs
  data <- output$data; n <- output$n
  S <- output$correlation_matrix

  # Obtain inverse and make a copy
  original_K <- K <- solve(S)

  # Pre-compute values
  half_n <- n / 2
  log_n <- log(n)
  lower_triangle <- lower.tri(K)

  # Get absolute K
  absolute <- abs(K)
  sorted <- sort(absolute[lower.tri(absolute)])

  # Create diagonal of inverse variances
  diagonal_K <- diag(diag(original_K))

  # Obtain fit function
  fit_COST <- switch(
    criterion,
    "aic" = aic_cost,
    "bic" = bic_cost,
    "ebic" = ebic_cost
  )

  # Compute fit values
  fits <- nvapply(sorted, function(value){

    # Set K back to original
    K <- original_K

    # Set K values to zero
    K[absolute <= value] <- 0

    # Get lower diagonal of K
    lower_K <- K[lower_triangle]

    # Set zeros
    zeros <- lower_K != 0

    # Compute BIC
    fit_COST(
      nonzero = lower_K[zeros], diagonal_K = diagonal_K,
      S = S, half_n = half_n, param = dimensions[2], log_n = log_n,
      lower_triangle = lower_triangle, zeros = zeros, gamma = gamma
    )

  })

  # Set K back to original
  K <- original_K

  # Set K values to zero
  K[absolute <= sorted[[which.min(fits)]]] <- 0

  # Get lower diagonal of K
  lower_K <- K[lower_triangle]

  # Set zeros
  zeros <- lower_K != 0

  # Get diagonal of absolute
  absolute_diag <- diag(absolute)

  # Obtain updated parameters
  updated <- optim(
    fn = fit_COST, par = lower_K[zeros], gr = nonreg_gradient,
    S = S, diagonal_K = diagonal_K,
    half_n = half_n, param = dimensions[2], log_n = log_n,
    lower_triangle = lower_triangle, zeros = zeros, gamma = gamma,
    method = "L-BFGS-B", lower = -absolute_diag,
    upper = absolute_diag
  )

  # Place in estimate
  diagonal_K[lower_triangle][zeros] <- updated$par
  diagonal_K <- t(diagonal_K)
  diagonal_K[lower_triangle][zeros] <- updated$par

  # Convert to zero-order correlations
  R <- solve(diagonal_K); diag(R) <- 1

  # Convert to partial correlations
  W <- -cov2cor(diagonal_K); diag(W) <- 0

  # Return results
  return(
    list(
      network = W, K = diagonal_K, R = R,
      correlation = S, criterion = criterion,
      value = updated$value, gamma = gamma,
      convergence = updated$convergence
    )
  )

}

#' @noRd
# Errors ----
# Updated 28.12.2024
nonreg_errors <- function(data, n, gamma, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "nonreg.network")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "nonreg.network")
    typeof_error(n, "numeric", "nonreg.network")
  }

  # 'verbose' errors
  length_error(verbose, 1, "nonreg.network")
  typeof_error(verbose, "logical", "nonreg.network")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

#' @noRd
# AIC cost function ----
# Updated 28.12.2024
aic_cost <- function(nonzero, diagonal_K, S, half_n, param, log_n, lower_triangle, zeros, gamma)
{

  # Set up K
  diagonal_K[lower_triangle][zeros] <- nonzero
  diagonal_K <- t(diagonal_K)
  diagonal_K[lower_triangle][zeros] <- nonzero

  # Obtain edge count
  E <- edge_count(diagonal_K, param)

  # Pre-compute SK
  SK <- S %*% diagonal_K

  # Return AIC
  return(-2 * (half_n * (log(det(diagonal_K)) - sum(diag(SK)))) + 2 * E)

}

#' @noRd
# BIC cost function ----
# Updated 28.12.2024
bic_cost <- function(nonzero, diagonal_K, S, half_n, param, log_n, lower_triangle, zeros, gamma)
{

  # Set up K
  diagonal_K[lower_triangle][zeros] <- nonzero
  diagonal_K <- t(diagonal_K)
  diagonal_K[lower_triangle][zeros] <- nonzero

  # Obtain edge count
  E <- edge_count(diagonal_K, param)

  # Pre-compute SK
  SK <- S %*% diagonal_K

  # Return BIC
  return(-2 * (half_n * (log(det(diagonal_K)) - sum(diag(SK)))) + E + log_n)

}

#' @noRd
# EBIC cost function ----
# Updated 28.12.2024
ebic_cost <- function(nonzero, diagonal_K, S, half_n, param, log_n, lower_triangle, zeros, gamma)
{

  # Set up K
  diagonal_K[lower_triangle][zeros] <- nonzero
  diagonal_K <- t(diagonal_K)
  diagonal_K[lower_triangle][zeros] <- nonzero

  # Obtain edge count
  E <- edge_count(diagonal_K, param)

  # Pre-compute SK
  SK <- S %*% diagonal_K

  # Return EBIC
  return(-2 * (half_n * (log(det(diagonal_K)) - sum(diag(SK)))) + E + log_n + 4 * E * gamma * log(param))

}

#' @noRd
# Non-regularized gradient function ----
# *Should* be generally OK for each cost function
# Updated 28.12.2024
nonreg_gradient <- function(nonzero, diagonal_K, S, half_n, param, log_n, lower_triangle, zeros, gamma)
{

  # Set up K
  diagonal_K[lower_triangle][zeros] <- nonzero
  diagonal_K <- t(diagonal_K)
  diagonal_K[lower_triangle][zeros] <- nonzero

  # Return gradient
  return(-2 * (solve(diagonal_K) - S)[lower_triangle][zeros] / half_n)

}
