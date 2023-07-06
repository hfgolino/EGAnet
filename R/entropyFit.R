#' Entropy Fit Index
#'
#' @description Computes the fit of a dimensionality structure using empirical entropy.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data Matrix or data frame.
#' Contains variables to be used in the analysis
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Total.Correlation.MM}{Miller-Madow correction for the total correlation of the dataset}
#'
#' \item{Entropy.Fit}{The Entropy Fit Index}
#'
#' \item{Entropy.Fit.MM}{Miller-Madow correction for the Entropy Fit Index}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA model
#' ega.wmt <- EGA(data = wmt)}
#'
#' # Compute entropy indices
#' entropyFit(data = wmt, structure = ega.wmt$wc)
#'
#' @references
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link[EGAnet]{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com> and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
# Entropy Fit Index
# Updated 06.07.2023
entropyFit <- function (data, structure)
{
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Get data dimensions
  dimensions <- dim(data)
  
  # Get number of communities
  communities <- unique_length(structure)
  
  # Get number of bins
  bins <- floor(sqrt(dimensions[1] / 5))
  
  # Obtain summed data
  if(communities == dimensions[2L]){
    summed_data <- data # scores are already summed
  }else{
    
    # Get sums by community
    summed_data <- nvapply(seq_len(communities), function(community){
      rowSums(data[,structure == community, drop = FALSE], na.rm = TRUE)  
    }, LENGTH = dimensions[1L])
    
  }
  
  # Set bin length
  bin_length <- bins + 1
  
  # Get bin cuts
  bin_cuts <- lapply(seq_len(communities), function(community){
    
    # Get range
    data_range <- range(summed_data[,community], na.rm = TRUE)
    
    # Return cuts
    return(
      cut(
        summed_data[,community], 
        breaks = seq.int(data_range[1L], data_range[2L], length.out = bin_length),
        include.lowest = TRUE
      )
    )
    
  })
  
  # Get frequencies
  bin_frequencies <- nvapply(bin_cuts, table, LENGTH = bins) / dimensions[1L]
  
  # Get entropies
  H <- nvapply(seq_len(communities), function(community){
    
    # Get non-zero frequencies
    bin_non_zero <- bin_frequencies[bin_frequencies[,community] > 0, community]
    
    # Return entropy
    return(entropy(bin_non_zero))
    
  })
  
  # Get joint frequency table
  joint_frequency <- count_table(
    do.call(cbind, bin_cuts)
  )$Value / dimensions[1L]
  
  # Get non-zero frequencies
  joint_non_zero <- joint_frequency[joint_frequency > 0]
  
  # Get joint entropy
  H_joint <- entropy(joint_non_zero)
  
  # Get maximum sums
  max_sum <- rowSums(data, na.rm = TRUE)
  
  # Obtain range
  max_range <- range(max_sum, na.rm = TRUE)
  
  # Count the cuts
  max_frequency <- count_table(
    cut(
      max_sum, 
      breaks = seq.int(max_range[1], max_range[2], length.out = bin_length), 
      include.lowest = TRUE
    )
  )$Value / dimensions[1L]
  
  # Get non-zero frequencies
  max_non_zero <- max_frequency[max_frequency > 0]
  
  # Get maximum entropy
  H_max <- entropy(max_non_zero)
  
  # Miller-Madow Bias Correction (for individual communities)
  ## Pre-compute denominator for corrections
  MM_denominator <- 2L * dimensions[1L]
  ## Entropy
  MM_non_zero <- colSums(bin_frequencies != 0, na.rm = TRUE)
  MM_H <- H + (MM_non_zero - 1) / MM_denominator
  ## Joint Entropy
  MM_joint_non_zero <- sum(joint_frequency != 0, na.rm = TRUE)
  MM_H_joint <- H_joint + (MM_joint_non_zero - 1) / MM_denominator
  
  # Compute mean entropy
  H_mean <- mean(H, na.rm = TRUE)
  
  # Pre-compute denominator for entropy fit measures
  EF_denominator <- (H_max - H_mean) * sqrt(communities)
  
  # Set up data frame
  return(
    fast.data.frame(
      c(
        sum(H, na.rm = TRUE) - H_joint,
        sum(MM_H, na.rm = TRUE) - MM_H_joint,
        H_mean - H_joint + EF_denominator,
        mean(MM_H, na.rm = TRUE) - MM_H_joint + EF_denominator,
        H_mean - H_joint
      ), ncol = 5L,
      colnames = c(
        "Total.Correlation", "Total.Correlation.MM",
        "Entropy.Fit", "Entropy.Fit.MM",
        "Average.Entropy"
      )
    )
  )
  
}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]; ega.wmt <- EGA(data, plot.EGA = FALSE)
# structure <- ega.wmt$wc
