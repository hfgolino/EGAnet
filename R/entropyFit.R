#' Fit of Factor Model using Empirical Entropy
#' @description Uses the \code{\link[entropy]{entropy.empirical}} function
#' from the \code{\link{entropy}} package to compute the entropy of a
#' factor model. Mutual information (\code{\link[entropy]{mi.plugin}})
#' is used to derive associations between factors. The final \code{entropyFit}
#' value is entropy + mutual information.
#' Lower values suggest better fit with the data
#'
#'
#' @param data A dataset
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#' @param bins An integer, indicating thhe number of bins to use.
#'
#' @return Returns the values of entropy per factor, the mean empirical entropy (i.e., the average entropy across factors),
#' the minimum joint entropy (i.e. joint entropy for all variables) and the adjusted entropy (ratio between joint entropy and
#' the minumum joint entropy).
#'
#'
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#'
#' entropyFit(data = wmt2[,7:24], structure = ega.wmt$wc, bins = 10)
#'
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
entropyFit <- function (data, structure, bins = 10)
{
  if(all(range(data)==c(0,1)))
  {data <- ifelse(data==1,2,1)}

  #convert structure to number if necessary
  if(is.character(structure))
  {
    uni <- unique(structure)
    num.comm <- structure

    for(i in 1:length(uni))
    {num.comm[which(num.comm==uniq[i])] <- i}

  } else {num.comm <- structure}

  #number of dimensions
  n <- max(num.comm)
  #communities sorted low to high
  uniq <- sort(unique(num.comm))

  #initialize entropy vector
  H <- vector("numeric",length=n)
  sums <- matrix(NA,nrow=nrow(data),ncol=n)
  bin.sums <- matrix(NA,nrow=nrow(data),ncol=n)
  Freq[,i] <- matrix(NA,nrow=nrow(data),ncol=n)
  seque <- matrix(NA,nrow=bins+1,ncol=n)

  #compute empirical entropy for each community
  for(i in 1:n)
  {
    sums[,i] <- rowSums(data[,which(num.comm==uniq[i])])
    seque[,i] <- seq(from = range(sums[,i])[1], to = range(sums[,i])[2], length.out = bins +
                       1)
    bin.sums[,i] <- cut(sums[,i], breaks = seque[,i], include.lowest = TRUE)
    Freq[,i] <- bin.sums[,i]/sum(bin.sums[,i])
    H[i] <- -sum(Freq[,i] * log(Freq[,i]))

  }

  # Joint Entropy:
  joint.table <- table(NA)
  for(i in 1:n){
    joint.table = table(cut(sums[,i], breaks = seque[,i], include.lowest = TRUE))
  }

  freq.joint <- joint.table/sum(joint.table)
  joint.entropy <- -sum(freq.joint[,i] * log(freq.joint[,i]))


  # Minumum Joint Entropy:
  seque2 <- matrix(NA,nrow=bins+1,ncol=ncol(data))
  joint.table2 <- table(NA)
  for(i in 1:ncol(data)){
    seque2[,i] <- seq(from = range(data[,i])[1], to = range(data[,i])[2], length.out = bins +
                        1)
    joint.table2 = table(cut(data[,i], breaks = seque2[,i], include.lowest = TRUE))
  }

  freq.min.joint <- joint.table2/sum(joint.table2)
  minimum.entropy <- -sum(freq.min.joint[,i] * log(freq.min.joint[,i]))

  #compute mean emprirical entropy
  #(empirical entropy per dimension)
  ent <- mean(H)

  result <- list()
  result$Ind.Entropy <- H
  result$Mean.Entropy <- ent
  result$Joint.Entropy <- joint.entropy
  result$Min.JointEntropy <- minimum.entropy
  result$Adj.Entropy <- joint.entropy/minimum.entropy
  return(result)
}
#----
