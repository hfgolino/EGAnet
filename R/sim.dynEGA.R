#' sim.dynEGA Data
#'
#' A simulated (multivariate time series) data with 24 variables,
#' 100 individual observations, 50 time points per individual and 
#' 2 groups of individuals
#' 
#' Data were generated using the \code{\link[EGAnet]{simDFM}} function
#' with the following arguments:
#' 
#' \strong{Group 1}
#' 
#' \code{simDFM(
#'   variab = 8, timep = 50, nfact = 3,
#'   error = 0.125, dfm = "DAFS",
#'   loadings = 0.60, autoreg = 0.80,
#'   crossreg = 0.00, var.shock = 0.36,
#'   cov.shock = 0.18, variation = TRUE
#' )}
#' 
#' \strong{Group 2}
#' 
#' \code{simDFM(
#'   variab = 6, timep = 50, nfact = 4,
#'   error = 0.125, dfm = "DAFS",
#'   loadings = 0.60, autoreg = 0.80,
#'   crossreg = 0.00, var.shock = 0.36,
#'   cov.shock = 0.18, variation = TRUE
#' )}
#'
#' @name sim.dynEGA
#'
#' @docType data
#'
#' @usage data(sim.dynEGA)
#'
#' @format A 5000 x 26 multivariate time series
#'
#' @keywords datasets
#'
#' @examples
#' data("sim.dynEGA")
#'
NULL
# Updated 31.10.2023
