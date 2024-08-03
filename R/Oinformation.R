#' @title O-information
#'
#' @description Computes the O-information of a dataset
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param base Numeric (length = 1).
#' Base to use for entropy.
#' Defaults to \code{exp(1)} or \code{2.718282}
#'
#' @return Returns a list containing:
#'
#' \item{Total.Cor}{The total correlation of the dataset}
#'
#' \item{Dual.Total.Cor}{The dual total correlation of the dataset}
#'
#' \item{O.Information}{The O-information of the dataset}
#'
#' \item{Normalized.Total.Cor}{Normalized total correlation}
#'
#' \item{Normalized.Dual.Total.Cor}{Normalized dual total correlation}
#'
#' \item{Normalized.O.Information}{Normalized O-information}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Compute O-information
#' Oinformation(wmt2[,7:24])
#'
#' @references
#' \strong{Formalization of O-information} \cr
#' Crutchfield, J. P. (1994). The calculi of emergence: Computation, dynamics and induction.
#' \emph{Physica D: Nonlinear Phenomena}, \emph{75}(1-3), 11-54.
#'
#' \strong{Applied implementation} \cr
#' Marinazzo, D., Van Roozendaal, J., Rosas, F. E., Stella, M., Comolatti, R., Colenbier, N., Stramaglia, S., & Rosseel, Y. (2024).
#' An information-theoretic approach to build hypergraphs in psychometrics.
#' \emph{Behavior Research Methods}, 1-23.
#'
#' @export
#'
# O-information ----
# Updated 03.08.2024
Oinformation <- function(data, base = 2.718282)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "totalCor")

  # 'base' errors
  typeof_error(base, "numeric", "Oinformation")
  length_error(base, 1, "Oinformation")

  # Ensure data is a matrix
  data <- as.matrix(usable_data(data, verbose = TRUE))

  # Obtain values
  tc <- totalCor(data, base = base)
  dtc <- dualTotalCor(data, base = base)

  # Return O-information
  return(
    list(
      Total.Cor = tc$Total.Cor,
      Dual.Total.Cor = dtc$Dual.Total.Cor,
      O.Information = tc$Total.Cor - dtc$Dual.Total.Cor,
      Normalized.Total.Cor = tc$Normalized,
      Normalized.Dual.Total.Cor = dtc$Normalized,
      Normalized.O.Information = tc$Normalized - dtc$Normalized
    )
  )

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]