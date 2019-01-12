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
#' the joint entropy (i.e. joint entropy for all factors), the entropy fit index (EFI), total correlation and modified versions
#' of all these indices using the Miller Madow correction.
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

  #compute empirical entropy for each community or item
  for(i in 1:n)
  {
    if(n != ncol(data)){
      sums[,i] <- rowSums(data[,which(num.comm==uniq[i])])
    } else{
        sums[,i] <- data[,i]
      }
    seque[,i] <- seq(from = range(sums[,i])[1], to = range(sums[,i])[2], length.out = bins+1)
    bin.sums[[i]] <- table(cut(sums[,i], breaks = seque[,i], include.lowest = TRUE))
    bin.sums2[,i] <- as.vector(unlist(bin.sums[[i]]))
    Freq[,i] <- bin.sums2[,i]/sum(bin.sums2[,i])
    H[i] <- -sum(ifelse(Freq[,i]>0,Freq[,i] * log(Freq[,i]),0))
  }

  # Joint Entropy:

  bin.sums3 <- data.frame(matrix(NA, nrow = nrow(data), ncol = n))
  joint.table <- vector("numeric")
  for(i in 1:n){
    bin.sums3[,i] <- cut(sums[,i], breaks = seque[,i], include.lowest = TRUE)
    joint.table = plyr::count(bin.sums3)$freq
  }

  freq.joint <- joint.table/sum(joint.table)
  joint.entropy <- -sum(ifelse(freq.joint >0,freq.joint * log(freq.joint),0))

    # Maximum Entropy:
    sums.max <- vector("numeric")
    sums.max <- rowSums(data)
    joint.table.max <- vector("numeric")
    seque.min <- seq(from = range(sums.max)[1], to = range(sums.max)[2], length.out = bins+1)
    bin.sums.min <- cut(sums.max, breaks = seque.min, include.lowest = TRUE)
    joint.table.max = plyr::count(bin.sums.min)$freq

    freq.joint.max <- joint.table.max/sum(joint.table.max)
    Hmax <- -sum(ifelse(freq.joint.max >0,freq.joint.max * log(freq.joint.max),0))

  # # Miller-Madow Bias Correction:
  # # Individual Factors:
   non.zero.bins1 <- vector("numeric",length=n)
   H.miller.madow <- vector("numeric",length=n)
   for(i in 1:n){
     non.zero.bins1[i] <- length(bin.sums2[bin.sums2[,i]!=0,i])
     H.miller.madow[i] <- H[i]+((non.zero.bins1[i]-1)/(2*(nrow(data))))
     }

  # Joint Entropy with Miller-Madow Bias Correction:

   non.zero.bins.joint <- length(joint.table[joint.table!=0])
   joint.miller.madow <- joint.entropy+((non.zero.bins.joint-1)/(2*(nrow(data))))


  #compute mean emprirical entropy
  #(empirical entropy per dimension)
  ent <- mean(H)

  result <- list()
  result$Ind.Entropy <- H
  result$Mean.Entropy <- ent
  result$Joint.Entropy <- joint.entropy
  result$H.Miller.Madow <- H.miller.madow
  result$Mean.Entropy.MM <- mean(H.miller.madow)
  result$Joint.Miller.Madow <- joint.miller.madow
  result$Total.Correlation <- sum(H)-joint.entropy
  result$Total.Correlation.MM <- sum(H.miller.madow)-joint.miller.madow
  result$Entropy.Fit <- (ent-joint.entropy)+((Hmax-ent)*sqrt(n))
  result$Entropy.Fit.MM <- (mean(H.miller.madow)-joint.miller.madow)+((Hmax-mean(H.miller.madow))*(sqrt(n)))
  result$Average.Entropy <- mean(H)-joint.entropy
  return(result)
}
#----
