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
#' @param comm A vector with community numbers or labels for each item.
#' Can be theoretical factors or communities detected by \code{\link{EGA}}
#' 
#' @return Returns the empirical entropy per factor (i.e., the average entropy across factors)
#' plus mutual information
#' 
#' @examples
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' 
#' entropyFit(data = wmt2[,7:24], comm = ega.wmt$wc)
#' 
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' 
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
entropyFit <- function (data, comm)
{
    if(all(range(data)==c(0,1)))
    {data <- ifelse(data==1,2,1)}
    
    #convert comm to number if necessary
    if(is.character(comm))
    {
        uni <- unique(comm)
        num.comm <- comm
        
        for(i in 1:length(uni))
        {num.comm[which(num.comm==uniq[i])] <- i}
        
    } else {num.comm <- comm}
    
    #number of communities
    n <- max(num.comm)
    #communities sorted low to high
    uniq <- sort(unique(num.comm))
    
    #initialize entropy vector
    H <- vector("numeric",length=n)
    sums <- matrix(NA,nrow=nrow(data),ncol=n)
    
    #compute empirical entropy each community
    for(i in 1:n)
    {
        len <- length(which(num.comm==uniq[i]))
        
        if(len==1)
        {sums[,i] <- data[,which(num.comm==uniq[i])]
        }else{sums[,i] <- rowSums(data[,which(num.comm==uniq[i])])}
        
        H[i] <- entropy::entropy.empirical(sums[,i])
    }
    
    #mutual information
    if(n!=1)
    {mi <- entropy::mi.empirical(sums)}
    
    #compute mean emprirical entropy
    #(empirical entropy per community)
    ent <- mean(H)
    
    result <- list()
    result$ind <- H
    if(n!=1)
    {
        result$mut <- mi
        result$full <- (ent + mi)
    }else{result$full <- ent}
    
    return(result)
}
#----