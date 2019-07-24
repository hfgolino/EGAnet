#' Entropy Fit Index
#'
#' @description Computes the fit of a dimensionality structure using empirical entropy.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data Matrix or data frame.
#' Contains variables to be used in the analysis
#'
#' @param structure A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Total.Correlation.MM}{Miller-Madow correction for the total correlation of the dataset}
#'
#' \item{Entropy.Fit}{The Entropy Fit Index}
#'
#' \item{Entropy.Fit.MM}{Miller-Madow correction for the Entropy Fit Index}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' @examples
#'
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
#' entropyFit(data = wmt, structure = ega.wmt$wc)
#'
#' @references
#' Golino, H. F., Moulder, R., Shi, D., Christensen, A. P., Neito, M. D., Nesselroade, J. R., & Boker, S. M. (under review)
#' Entropy Fit Index: A new fit measure for assessing the structure and dimensionality of multiple latent variables.
#' Retrieved from: https://www.researchgate.net/profile/Hudson_Golino/publication/333753928_Entropy_Fit_Index_A_New_Fit_Measure_for_Assessing_the_Structure_and_Dimensionality_of_Multiple_Latent_Variables/
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link[EGAnet]{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com> and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
#Entropy Fit Index
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

  result <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  colnames(result) <- c("Total.Correlation", "Total.Correlation.MM","Entropy.Fit",
                        "Entropy.Fit.MM", "Average.Entropy")
  result$Total.Correlation <- sum(H)-joint.entropy
  result$Total.Correlation.MM <- sum(H.miller.madow)-joint.miller.madow
  result$Entropy.Fit <- (ent-joint.entropy)+((Hmax-ent)*(sqrt(n)))
  result$Entropy.Fit.MM <- (mean(H.miller.madow)-joint.miller.madow)+((Hmax-ent)*(sqrt(n)))
  result$Average.Entropy <- mean(H)-joint.entropy
  return(result)
}
#----
