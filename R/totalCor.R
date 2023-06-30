#' Total Correlation
#'
#' Computes the total correlation of a dataset
#'
#' @param data Matrix or data frame.
#' Variables to be used in the analysis
#' 
#' @return Returns a list containing:
#' 
#' \item{Ind.Entropies}{Individual entropies for each variable}
#' 
#' \item{Joint.Entropy}{The joint entropy of the dataset}
#' 
#' \item{Total.Cor}{The total correlation of the dataset}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Compute total correlation
#' totalCor(wmt2[,7:24])
#' 
#' @references 
#' Watanabe, S. (1960).
#' Information theoretical analysis of multivariate correlation.
#' \emph{IBM Journal of Research and Development} \emph{4}, 66-82.
#' 
#' # Implementation
#' Felix, L. M., Mansur-Alves, M., Teles, M., Jamison, L., & Golino, H. (2021).
#' Longitudinal impact and effects of booster sessions in a cognitive training program for healthy older adults.
#' \emph{Archives of Gerontology and Geriatrics}, \emph{94}, 104337.
#' 
#' @export
#'
# Total Correlation
# Updated 29.06.2023
totalCor <- function(data)
{
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Get data dimensions
  dimensions <- dim(data)
  
  # Get number of bins
  bins <- floor(sqrt(dimensions[1] / 5))
  
  # Set bin length
  bin_length <- bins + 1
  
  # Get bin cuts
  bin_cuts <- lapply(seq_len(dimensions[2]), function(variable){
    
    # Get range
    data_range <- range(data[,variable], na.rm = TRUE)
    
    # Return cuts
    return(
      cut(
        data[,variable],
        breaks = seq.int(data_range[1], data_range[2], length.out = bin_length),
        include.lowest = TRUE
      )
    )
    
  })
  
  # Get bin sums
  bin_sums <- nnapply(bin_cuts, function(x){
    table(x)
  }, LENGTH = bins)
  
  # Get frequencies
  bin_frequencies <- bin_sums / dimensions[1]
  
  # Get entropies
  H <- nnapply(seq_len(dimensions[2]), function(variable){
    
    # Get non-zero frequencies
    bin_non_zero <- bin_frequencies[bin_frequencies[,variable] > 0, variable]
    
    # Return entropy
    return(entropy(bin_non_zero))
    
  })
  
  # Get joint frequency table
  joint_frequency <- count_table(
    do.call(cbind, bin_cuts)
  )$Value / dimensions[1]
  
  # Get non-zero frequencies
  joint_non_zero <- joint_frequency[joint_frequency > 0]
  
  # Get joint entropy
  H_joint <- entropy(joint_non_zero)
  
  # Return list
  return(
    list(
      Ind.Entropies = H,
      Joint.Entropy = H_joint,
      Total.Cor = sum(H, na.rm = TRUE) - H_joint
    )
  )

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]