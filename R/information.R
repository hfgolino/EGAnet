#' @title Information Theory Metrics
#'
#' @description A general function to compute several different information theory metrics
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param base Numeric (length = 1).
#' Base of logarithm to use for entropy.
#' Common options include:
#'
#' \itemize{
#'
#' \item \code{2} --- bits
#'
#' \item \code{2.718282} --- nats
#'
#' \item \code{10} --- bans
#'
#' }
#'
#' Defaults to \code{exp(1)} or \code{2.718282}
#'
#' @param bins Numeric (length = 1).
#' Number of bins if data are not discrete.
#' Defaults to \code{floor(sqrt(nrow(data) / 5))}
#'
#' @param statistic Character.
#' Information theory statistics to compute.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"entropy"} --- Shannon's entropy (Shannon, 1948) for each variable in \code{data}.
#' Values range from \code{0} to \code{log(k)} where \code{k} is the number of categories for the variable
#'
#' \item \code{"joint.entropy"} --- shared uncertainty over all variables in \code{data}.
#' Values range from the maximum of the individual entropies to the sum of individual entropies
#'
#' \item \code{"conditional.entropy"} --- uncertainty remaining after considering all other
#' variables in \code{data}. Values range from \code{0} to the individual entropy of the
#' conditioned variable
#'
#' \item \code{"total.correlation"} --- generalization of mutual information to more than
#' two variables (Watanabe, 1960). Quantifies the redundancy of information in \code{data}.
#' Values range from \code{0} to the sum of individual entropies minus the maximum of the
#' individual entropies
#'
#' \item \code{"dual.total.correlation"} --- "shared randomness" or total uncertainty remaining in
#' the \code{data} (Han, 1978). Values range from \code{0} to joint entropy
#'
#' \item \code{"o.information"} --- quantifies the extent to which the \code{data} is represented
#' by lower-order (\code{> 0}; redundancy) or higher-order (\code{< 0}; synergy) constraint
#' (Crutchfield, 1994)
#'
#' }
#'
#' By default, all statistics are computed
#'
#' @return Returns list containing \emph{only} requested \code{statistic}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # All measures
#' information(wmt2[,7:24])
#'
#' # One measures
#' information(wmt2[,7:24], statistic = "joint.entropy")
#'
#' @references
#' \strong{Shannon's entropy} \cr
#' Shannon, C. E. (1948). A mathematical theory of communication.
#' \emph{The Bell System Technical Journal}, \emph{27}(3), 379-423.
#'
#' \strong{Formalization of total correlation} \cr
#' Watanabe, S. (1960).
#' Information theoretical analysis of multivariate correlation.
#' \emph{IBM Journal of Research and Development} \emph{4}, 66-82.
#'
#' \strong{Applied implementation of total correlation} \cr
#' Felix, L. M., Mansur-Alves, M., Teles, M., Jamison, L., & Golino, H. (2021).
#' Longitudinal impact and effects of booster sessions in a cognitive training program for healthy older adults.
#' \emph{Archives of Gerontology and Geriatrics}, \emph{94}, 104337.
#'
#' \strong{Formalization of dual total correlation} \cr
#' Te Sun, H. (1978).
#' Nonnegative entropy measures of multivariate symmetric correlations.
#' \emph{Information and Control}, \emph{36}, 133-156.
#'
#' \strong{Formalization of O-information} \cr
#' Crutchfield, J. P. (1994). The calculi of emergence: Computation, dynamics and induction.
#' \emph{Physica D: Nonlinear Phenomena}, \emph{75}(1-3), 11-54.
#'
#' \strong{Applied implementation of O-information} \cr
#' Marinazzo, D., Van Roozendaal, J., Rosas, F. E., Stella, M., Comolatti, R., Colenbier, N., Stramaglia, S., & Rosseel, Y. (2024).
#' An information-theoretic approach to build hypergraphs in psychometrics.
#' \emph{Behavior Research Methods}, 1-23.
#'
#' @export
#'
# Information function ----
# Updated 06.08.2024
information <- function(
    data, base = 2.718282, bins = floor(sqrt(nrow(data) / 5)),
    statistic = c(
      "entropy", "joint.entropy",
      "conditional.entropy", "total.correlation",
      "dual.total.correlation", "o.information"
    )
)
{

  # Argument errors (return data in case of tibble)
  data <- information_errors(data, base, bins)

  # Check for missing arguments (argument, default, function)
  statistic <- set_default(statistic, "entropy", information, several.ok = TRUE)

  # Initialize entropy list
  entropy_list <- vector("list", length(statistic))
  names(entropy_list) <- statistic

  # Get data dimensions
  dimensions <- dim(data)

  # Obtain variable sequence
  variable_sequence <- seq_len(dimensions[2L])

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
  entropy_list$entropy <- nvapply(variable_sequence, function(variable){

    # Get non-zero frequencies
    bin_non_zero <- bin_frequencies[bin_frequencies[,variable] > 0, variable]

    # Return entropy
    return(entropy(bin_non_zero, base = base))

  })

  # All other statistics need joint entropy
  if(length(statistic[statistic != "entropy"]) != 0){
    entropy_list$joint.entropy <- joint_entropy(bin_cuts, base)
  }

  # Check for statistics that need conditional entropy
  if(
    any(c("conditional.entropy", "dual.total.correlation", "o.information") %in% statistic)
  ){

    # Get joint entropy without each variable
    H_individual_joint <- nvapply(
      variable_sequence, function(i){
        return(joint_entropy(bin_cuts[-i], base))
      }
    )

    # Get conditional entropies
    entropy_list$conditional.entropy <- entropy_list$joint.entropy - H_individual_joint

  }

  # Check for statistics that need total correlation
  if(any(c("total.correlation", "o.information") %in% statistic)){
    entropy_list$total.correlation <- sum(entropy_list$entropy, na.rm = TRUE) - entropy_list$joint.entropy
  }

  # Check for statistics that need dual total correlation
  if(any(c("dual.total.correlation", "o.information") %in% statistic)){
    entropy_list$dual.total.correlation <- entropy_list$joint.entropy -
                                           sum(entropy_list$conditional.entropy, na.rm = TRUE)
  }

  # Check for O-information
  if("o.information" %in% statistic){
    entropy_list$o.information <- entropy_list$total.correlation - entropy_list$dual.total.correlation
  }

  # Return list
  return(entropy_list[statistic])

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]

#' @noRd
# Errors ----
# Updated 04.08.2024
information_errors <- function(data, base, bins)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "information")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'base' errors
  typeof_error(base, "numeric", "information")
  length_error(base, 1, "information")

  # 'bins' errors
  typeof_error(bins, "numeric", "information")
  length_error(bins, 1, "information")

  # Ensure data is a matrix
  return(as.matrix(usable_data(data, verbose = TRUE)))

}

#' @noRd
# Basic entropy calculation ----
# Returns entropy value
# Updated 29.06.2023
entropy <- function(values, base = 2.718282)
{
  return(-sum(values * log(values, base = base), na.rm = TRUE))
}

#' @noRd
# Get joint entropy ----
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

#' @noRd
# Get mutual information ----
# Updated 04.08.2024
mutual_information <- function(Hx, Hy, Hxy)
{
  return(Hx + Hy - Hxy)
}
