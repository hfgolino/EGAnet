#' Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (TEFI) of a dimensionality structure using Von Neumman's entropy when the input is a correlation matrix.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data A dataset or a correlation matrix
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{Entropy.Fit}{The Entropy Fit Index using Von Neumman's entropy}
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
#' tefi(data = wmt, structure = ega.wmt$wc)
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
tefi <- function(data, structure){
  if(!is.matrix(data)){
    cor1 <- qgraph::cor_auto(data)/ncol(data)
    h.vn <- -matrixcalc::matrix.trace(cor1%*%log(cor1))

    n <- max(structure)
    cor.fact <- vector("list")
    eigen.fact <- vector("list")
    l.eigen.fact <- vector("list")
    h.vn.fact <- vector("list")
    for(i in 1:n){
      cor.fact[[i]] <- qgraph::cor_auto(data[,which(structure==unique(structure)[i])])
      cor.fact[[i]] <- cor.fact[[i]]/ncol(cor.fact[[i]])
      h.vn.fact[[i]] <- -matrixcalc::matrix.trace(cor.fact[[i]]%*%log(cor.fact[[i]]))
    }

  } else{
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
