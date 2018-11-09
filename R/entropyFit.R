#' Fit of Factor Model using Empirical Entropy
#' @description Compute compute the entropy fit of a given structure.
#' Lower values of the adjusted entropy fit and the adjusted entropy fit
#' ranging from 0 to 1 suggest better fit with the data
#'
#' @param data A dataset
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#'
#' @return Returns the values of entropy per factor, the mean empirical entropy (i.e., the average entropy across factors),
#' the joint entropy (i.e. joint entropy for all factors), a modified empirical entropy per factor (that uses the frequencies of the sumscores),
#' a mean modified empirical entropy. It also returns an adjusted entropy (ratio between the mean empirical entropy and
#' the joint entropy) and an adjusted entropy varying from 0 to 1.
#'
#'
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#'
#' entropyFit(data = wmt2[,7:24], structure = ega.wmt$wc)
#'
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
entropyFit <- function (data, structure)
{
  require(plyr)
  require(infotheo)

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

## Traditional Entropy:

  #number of dimensions
  n <- max(num.comm)
  #communities sorted low to high
  uniq <- sort(unique(num.comm))

  #initialize entropy vector
  H <- vector("numeric",length=n)
  bins <- floor(sqrt(nrow(data)/5))
  seque <- matrix(NA,nrow=bins+1,ncol=n)
  sums <- matrix(NA,nrow=nrow(data),ncol=n)
  bin.sums <- vector("list", n)
  bin.sums2 <- matrix(NA, nrow=bins, ncol = n)
  Freq <- matrix(NA,nrow=bins,ncol=n)

  #compute empirical entropy for each community
  for(i in 1:n)
  {
    sums[,i] <- rowSums(data[,which(num.comm==uniq[i])])
    seque[,i] <- seq(from = range(sums[,i])[1], to = range(sums[,i])[2], length.out = bins+1)
    bin.sums[[i]] <- table(cut(sums[,i], breaks = seque[,i], include.lowest = TRUE))
    bin.sums2[,i] <- as.vector(unlist(bin.sums[[i]]))
    Freq[,i] <- bin.sums2[,i]/sum(bin.sums2[,i])
    H[i] <- -sum(ifelse(Freq[,i]>0,Freq[,i] * log(Freq[,i]),0))
  }

  # Joint Entropy:

  bin.sums3 <- data.frame(matrix(NA, nrow = nrow(data), ncol = n))
  joint.table <- vector()
  for(i in 1:n){
    bin.sums3[,i] <- cut(sums[,i], breaks = seque[,i], include.lowest = TRUE)
    joint.table = plyr::count(bin.sums3)$freq
  }

  freq.joint <- joint.table/sum(joint.table)
  joint.entropy <- -sum(ifelse(freq.joint >0,freq.joint * log(freq.joint),0))

# Modified Entropy:

  #initialize entropy vector
  Hmod <- vector("numeric",length=n)
  sums.mod <- matrix(NA,nrow=nrow(data),ncol=n)
  Freq.mod <- matrix(NA,nrow=nrow(data),ncol=n)

  #compute empirical entropy for each community
  for(i in 1:n)
  {
    sums.mod[,i] <- rowSums(data[,which(num.comm==uniq[i])])
    Freq.mod[,i] <- sums.mod[,i]/sum(sums.mod[,i])
    Hmod[i] <- -sum(ifelse(Freq.mod[,i]>0,Freq.mod[,i] * log(Freq.mod[,i]),0))
  }


  #compute mean emprirical entropy
  #(empirical entropy per dimension)
  ent <- mean(H)
  mod.ent <- mean(Hmod)

  result <- list()
  result$Ind.Entropy <- H
  result$Mean.Entropy <- ent
  result$Joint.Entropy <- joint.entropy
  result$Ind.Mod.Entropy <- Hmod
  result$Mean.Mod.Entropy <- mod.ent
  result$Adj.Entropy <- ent-joint.entropy
  result$Adj.Entropy2 <- exp(ent/joint.entropy)/(1+exp(ent/joint.entropy))
  result$Adj.Entropy3 <- mean(H-joint.entropy)
  result$Adj.Entropy4 <- 1-(exp(mean(H-joint.entropy))/(1+exp(mean(H-joint.entropy))))
  return(result)
}
#----
