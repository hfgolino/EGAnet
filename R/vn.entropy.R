#' Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit of a dimensionality structure using Von Neumman's entropy when the input is a correlation matrix.
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
#' \donttest{
#' #estimate EGA model
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#'
#' #compute entropy indices
#' vn.entropy(data = wmt2[,7:24], structure = ega.wmt$wc)
#'}
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Entropy Fit Index
# VN Entropy Function (for correlation matrices)
vn.entropy <- function(data, structure){

  uniq <- unique(structure)
  num.comm <- structure

  if(!is.matrix(data)){
    cor1 <- qgraph::cor_auto(data)/ncol(data)
    eigen1 <- eigen(cor1)$values
    h.vn <- -sum(eigen1*log(eigen1))

    n <- max(structure)
    cor.fact <- vector("list")
    eigen.fact <- vector("list")
    l.eigen.fact <- vector("list")
    h.vn.fact <- vector("list")
    for(i in 1:n){
      cor.fact[[i]] <- qgraph::cor_auto(data[,which(structure==unique(structure)[i])])
      cor.fact[[i]] <- cor.fact[[i]]/ncol(cor.fact[[i]])
      eigen.fact[[i]] <- eigen(cor.fact[[i]])$values
      l.eigen.fact[[i]] <- eigen.fact[[i]]*log(eigen.fact[[i]])
      h.vn.fact[[i]] <- -sum(l.eigen.fact[[i]])
    }

    # Joint entropy using Kronecker product of a list of matrices
    cor.joint <- vector("list", n)
    for(i in 1:n){
      cor.joint[[i]] <- cor1[which(num.comm==uniq[i]),which(num.comm==uniq[i])]/table(num.comm)[[i]]
    }
    kronecker.product <- klin::klin.klist(cor.joint)
    eigen.kronecker <- eigen(kronecker.product)$values
    h.vn.joint <- -sum(eigen.kronecker*log(eigen.kronecker))

  } else{
    cor1 <- data/ncol(data)
    eigen1 <- eigen(cor1)$values
    h.vn <- -sum(eigen1*log(eigen1))

    n <- max(structure)
    cor.fact <- vector("list")
    eigen.fact <- vector("list")
    l.eigen.fact <- vector("list")
    h.vn.fact <- vector("list")
    for(i in 1:n){
      cor.fact[[i]] <- data[which(structure==unique(structure)[i]),which(structure==unique(structure)[i])]
      cor.fact[[i]] <- cor.fact[[i]]/ncol(cor.fact[[i]])
      eigen.fact[[i]] <- eigen(cor.fact[[i]])$values
      l.eigen.fact[[i]] <- eigen.fact[[i]]*log(eigen.fact[[i]])
      h.vn.fact[[i]] <- -sum(l.eigen.fact[[i]])
    }
    # Joint entropy using Kronecker product of a list of matrices
    cor.joint <- vector("list", n)
    for(i in 1:n){
      cor.joint[[i]] <- cor1[which(num.comm==uniq[i]),which(num.comm==uniq[i])]/table(num.comm)[[i]]
    }
    kronecker.product <- klin::klin.klist(cor.joint)
    eigen.kronecker <- eigen(kronecker.product)$values
    h.vn.joint <- -sum(eigen.kronecker*log(eigen.kronecker))
  }

  h.vn.fact2 <- unlist(h.vn.fact)

  # Difference between Max the sum of the factor entropies:
  Hdiff <- h.vn-mean(h.vn.fact2)

  results <- data.frame(matrix(NA, nrow = 1, ncol = 4))
  colnames(results) <- c("Entropy.Fit","Entropy.Fit2", "Total.Correlation","Average.Entropy")
  results$Entropy.Fit <- ((mean(h.vn.fact2)-h.vn.joint))-((Hdiff-(mean(h.vn.fact2))*sqrt(n)))
  results$Total.Correlation <- sum(h.vn.fact2)-h.vn.joint
  results$Average.Entropy <- mean(h.vn.fact2)-h.vn.joint
  return(results)
}
