#' Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (TEFI) of a dimensionality structure using Von Neumman's entropy when the input is a correlation matrix.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data A dataframe or correlation matrix
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{VN.Entropy.Fit}{The Entropy Fit Index using Von Neumman's entropy}
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA model
#' ega.wmt <- EGA(data = wmt, model = "glasso")
#'
#' }
#'
#' # Compute entropy indices
#' tefi(data = ega.wmt$correlation, structure = ega.wmt$wc)
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' 
#' @references 
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' doi: \href{https://doi.org/10.31234/osf.io/mtka2}{10.31234/osf.io/mtka2}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com>, and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
# Updated 21.10.2020
tefi <- function(data, structure){
  if(ncol(data)!=nrow(data)){
    data <- qgraph::cor_auto(data)
  }
  
  if(any(is.na(structure)))
  {
    rm.vars <- which(is.na(structure))
    
    warning(paste("Some variables did not belong to a dimension:", colnames(data)[rm.vars]))
    message("Use caution: These variables have been removed from the TEFI calculation")
    
    data <- data[-rm.vars, -rm.vars]
  }
  
  data <- abs(data)
  cor1 <- data/ncol(data)
  h.vn <- -matrixcalc::matrix.trace(cor1%*%(log(cor1)))
  
  n <- max(structure)
  cor.fact <- vector("list")
  eigen.fact <- vector("list")
  l.eigen.fact <- vector("list")
  h.vn.fact <- vector("list")
  for(i in 1:n){
    cor.fact[[i]] <- data[which(structure==unique(structure)[i]),which(structure==unique(structure)[i])]
    cor.fact[[i]] <- cor.fact[[i]]/ncol(cor.fact[[i]])
    h.vn.fact[[i]] <- -matrixcalc::matrix.trace(cor.fact[[i]]%*%log(cor.fact[[i]]))
  }
  
  h.vn.fact2 <- unlist(h.vn.fact)
  
  # Difference between Max the sum of the factor entropies:
  Hdiff <- h.vn-sum(h.vn.fact2)
  results <- data.frame(matrix(NA, nrow = 1, ncol = 3))
  colnames(results) <- c("VN.Entropy.Fit", "Total.Correlation","Average.Entropy")
  results$VN.Entropy.Fit <- (mean(h.vn.fact2)-h.vn)+(Hdiff*(sqrt(n)))
  results$Total.Correlation <- sum(h.vn.fact2)-h.vn
  results$Average.Entropy <- mean(h.vn.fact2)-h.vn
  return(results)
}
#----
