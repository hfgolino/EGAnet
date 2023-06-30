#' Total Correlation Matrix
#'
#' Computes the pairwise total correlation for a dataset
#'
#' @param data Matrix or data frame.
#' Variables to be used in the analysis
#'
#' @return Returns a square matrix with pairwise total correlations
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Compute total correlation matrix
#' totalCorMat(wmt2[,7:24])
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
# Updated 28.06.2023
totalCorMat <- function(data){
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Get dimensions of data
  dimensions <- dim(data)
  
  # Ensure variable names
  data <- ensure_dimension_names(data)
  
  # Get variable names
  variable_names <- dimnames(data)[[2]]
  
  # Initialize total correlation matrix
  total_correlation <- matrix(
    nrow = dimensions[2], ncol = dimensions[2],
    dimnames = list(variable_names, variable_names)
  )
  
  # Fill matrix
  for(i in seq_len(dimensions[2])){
    
    # Loop over other variables
    for(j in i:dimensions[2]){
      
      # Fill both sides of the matrix
      total_correlation[i,j] <- total_correlation[j,i] <- totalCor(data[,c(i, j)])$Total.Cor
      
    }
    
  }
  
  # Return total correlation matrix
  return(total_correlation)
}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]

