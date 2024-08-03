#' @title Dual Total Correlation
#'
#' @description Computes the dual total correlation of a dataset
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
#' \item{Cond.Entropies}{Conditional entropies}
#'
#' \item{Joint.Entropy}{The joint entropy of the dataset}
#'
#' \item{Dual.Total.Cor}{The dual total correlation of the dataset}
#'
#' \item{Normalized}{Dual total correlation divided by the joint entropy}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Compute dual total correlation
#' dualTotalCor(wmt2[,7:24])
#'
#' @references
#' \strong{Formalization of dual total correlation used} \cr
#' Varley, T. F., Pope, M., Faskowitz, J., & Sporns, O. (2023).
#' Multivariate information theory uncovers synergistic subsystems of the human cerebral cortex.
#' \emph{Communications Biology}, \emph{6}(1), 451.
#'
#' @export
#'
# Dual Total Correlation ----
# Updated 03.08.2024
dualTotalCor <- function(data, base = 2.718282)
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

  # Variable sequence
  variable_sequence <- seq_len(dimensions[2L])

  # Get number of bins
  bins <- floor(sqrt(dimensions[1L] / 5))

  # Set bin length
  bin_length <- bins + 1L

  # Get bin cuts
  bin_cuts <- lapply(variable_sequence, function(variable){

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
  H <- nvapply(variable_sequence, function(variable){

    # Get non-zero frequencies
    bin_non_zero <- bin_frequencies[bin_frequencies[,variable] > 0, variable]

    # Return entropy
    return(entropy(bin_non_zero, base = base))

  })

  # Get joint entropy of positive non-zero values
  H_joint <- joint_entropy(bin_cuts, base)

  # Get joint entropy without each variable
  H_individual_joint <- nvapply(
    variable_sequence, function(i){
      return(joint_entropy(bin_cuts[-i], base))
    }
  )

  # Get conditional entropies
  H_conditional <- H_joint - H_individual_joint

  # Get dual total correlation
  Dual.Total.Cor <- H_joint - sum(H_conditional, na.rm = TRUE)

  # Return list
  return(
    list(
      Ind.Entropies = H,
      Cond.Entropies = H_conditional,
      Joint.Entropy = H_joint,
      Dual.Total.Cor = Dual.Total.Cor,
      Normalized = Dual.Total.Cor / H_joint
    )
  )

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]


#' @noRd
# Get joint entropy
# Updated 03.08.2024
joint_entropy <- function(bin_cuts, base)
{

  # Get joint frequency table
  joint_frequency <- count_table(
    do.call(cbind, bin_cuts), proportion = TRUE
  )$Value

  # Return joint entropy
  return(entropy(joint_frequency[joint_frequency > 0], base = base))

}

