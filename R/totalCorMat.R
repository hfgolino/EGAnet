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
#' # Compute total correlation
#' totalCorMat(wmt2[,7:24])
#'
#' @references
#' Watanabe, S. (1960).
#' Information theoretical analysis of multivariate correlation.
#' \emph{IBM Journal of Research and Development} \emph{4}, 66-82.
#'
#'
#' @export
#'
# Total COrrelation
# Updated 14.02.2021
totalCorMat <- function(data){
  TotCorMat <- matrix(NA,nrow=ncol(data),ncol=ncol(data))
  TotCorMat <- matrix(apply(expand.grid( x = 1:ncol(data), y = 1:ncol(data)), 1,
                              function(x){ return(totalCor( data[ , c(x[1], x[2])])$Total.Cor)}),
                      ncol = ncol(data))
  colnames(TotCorMat) <- colnames(data)
  rownames(TotCorMat) <- colnames(data)
  return(TotCorMat)
}



