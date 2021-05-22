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
#' 
#' @export
#'
# Total Correlation
# Updated 22.05.2021
totalCor <- function(data){
  #number of dimensions
  n <- ncol(data)

  #initialize entropy vector
  H <- vector("numeric", length = n)
  bins <- floor(sqrt(nrow(data) / 5))
  seque <- matrix(NA, nrow = bins + 1, ncol = n)
  sums <- matrix(NA, nrow = nrow(data), ncol = n)
  bin.sums <- vector("list", n)
  bin.sums2 <- matrix(NA, nrow = bins, ncol = n)
  Freq <- matrix(NA, nrow = bins, ncol = n)

  #compute empirical entropy for each community or item
  for(i in 1:n)
  {
    seque[,i] <- seq(from = range(data[,i], na.rm = TRUE)[1], to = range(data[,i], na.rm = TRUE)[2], length.out = bins + 1)
    bin.sums[[i]] <- table(cut(data[,i], breaks = seque[,i], include.lowest = TRUE))
    bin.sums2[,i] <- as.vector(unlist(bin.sums[[i]]))
    Freq[,i] <- bin.sums2[,i]/sum(bin.sums2[,i])
    H[i] <- -sum(ifelse(Freq[,i]>0,Freq[,i] * log(Freq[,i]),0))
  }

  # Joint Entropy:
  bin.sums3 <- data.frame(matrix(NA, nrow = nrow(data), ncol = n))
  for(i in 1:n){
    bin.sums3[,i] <- cut(data[,i], breaks = seque[,i], include.lowest = TRUE)
  }
  joint.table <- plyr::count(bin.sums3)$freq

  freq.joint <- joint.table / sum(joint.table)
  joint.entropy <- -sum(ifelse(freq.joint >0,freq.joint * log(freq.joint),0))

  results <- vector("list")
  results$Ind.Entropies <- H
  results$Joint.Entropy <- joint.entropy
  results$Total.Cor <- sum(H)-joint.entropy
  return(results)
}
