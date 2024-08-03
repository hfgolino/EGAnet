#' @title Total Correlation
#'
#' @description Computes the total correlation of a dataset
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
#' \item{Ind.Entropies}{Individual entropies for each variable}
#'
#' \item{Joint.Entropy}{The joint entropy of the dataset}
#'
#' \item{Total.Cor}{The total correlation of the dataset}
#'
#' \item{Normalized}{Total correlation divided by the sum of the individual entropies minus the maximum of the individual entropies}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Compute total correlation
#' totalCor(wmt2[,7:24])
#'
#' @references
#' \strong{Formalization of total correlation} \cr
#' Watanabe, S. (1960).
#' Information theoretical analysis of multivariate correlation.
#' \emph{IBM Journal of Research and Development} \emph{4}, 66-82.
#'
#' \strong{Applied implementation} \cr
#' Felix, L. M., Mansur-Alves, M., Teles, M., Jamison, L., & Golino, H. (2021).
#' Longitudinal impact and effects of booster sessions in a cognitive training program for healthy older adults.
#' \emph{Archives of Gerontology and Geriatrics}, \emph{94}, 104337.
#'
#' @export
#'
# Total Correlation ----
# Updated 03.08.2024
totalCor <- function(data, base = 2.718282)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "totalCor")

  # 'base' errors
  typeof_error(base, "numeric", "totalCor")
  length_error(base, 1, "totalCor")

  # Ensure data is a matrix
  data <- as.matrix(usable_data(data, verbose = TRUE))

  # Get data dimensions
  dimensions <- dim(data)

  # Get number of bins
  bins <- floor(sqrt(dimensions[1L] / 5))

  # Set bin length
  bin_length <- bins + 1L

  # Get bin cuts
  bin_cuts <- lapply(seq_len(dimensions[2L]), function(variable){

    # Get range
    data_range <- range(data[,variable], na.rm = TRUE)

    # Return cuts
    return(
      cut(
        data[,variable],
        breaks = seq.int(data_range[1L], data_range[2L], length.out = bin_length),
        include.lowest = TRUE
      )
    )

  })

  # Get bin frequencies
  bin_frequencies <- nvapply(bin_cuts, fast_table, LENGTH = bins) / dimensions[1]

  # Get entropies
  H <- nvapply(seq_len(dimensions[2L]), function(variable){

    # Get non-zero frequencies
    bin_non_zero <- bin_frequencies[bin_frequencies[,variable] > 0, variable]

    # Return entropy
    return(entropy(bin_non_zero, base = base))

  })

  # Get joint entropy of positive non-zero values
  H_joint <- joint_entropy(bin_cuts, base = base) # function in `dualTotatCor.R`
  H_sum <- sum(H, na.rm = TRUE)
  Total.Cor <- H_sum - H_joint

  # Return list
  return(
    list(
      Ind.Entropies = H,
      Joint.Entropy = H_joint,
      Total.Cor = Total.Cor,
      Normalized = Total.Cor / (H_sum - max(H, na.rm = TRUE))
    )
  )

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]