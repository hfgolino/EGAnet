#' sim.dynEGA Data
#'
#' A simulated (multivariate time series) data with 24 variables,
#' 100 individual observations, 50 time points per individual and 
#' 2 groups of individuals.
#' 
#' Data were generated using the \code{\link[EGAnet]{simDFM}} function
#' with the following arguments:
#' 
#' \code{simDFM(
#'   variab = 12, timep = 50, nfact = 2,
#'   error = 0.175, dfm = "DAFS",
#'   loadings = 0.70, autoreg = 0.50,
#'   crossreg = 0.10, var.shock = 0.09,
#'   cov.shock = 0.30, variation = TRUE
#' )}
#' 
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
# Updated 04.09.2022
