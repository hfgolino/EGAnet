#' EGA Optimal Model Fit using Entropy
#' 
#' @description Estimates the best fitting model using \code{\link[EGA]{EGA}}.
#' The number of steps in the \code{\link[igraph]{cluster_walktrap}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link[EGA]{entropyFit}}.
#'
#' @param data A dataset
#'
#' @param model A string indicating the method to use.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\code{\strong{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#' 
#' \item{\code{\strong{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#' 
#' }
#' 
#' @param steps Range of steps to be used in the model selection.
#' Defaults from 3 to 8 steps (based on Pons & Latapy, 2006)
#'
#' @return Returns a list containing:
#' 
#' \item{EGA}{The \code{\link[EGA]{EGA}} output for the best fitting model}
#' 
#' \item{steps}{The number of steps used in the best fitting model from
#' the \code{\link[igraph]{cluster_walktrap}} algorithm}
#' 
#' \item{EntropyFit}{The Entropy Fit Index for the unique solutions given the range of steps
#' (vector names represent the number of steps)}
#' 
#' \item{Lowest.EntropyFit}{The lowest value for the Entropy Fit Index}
#'
#' @examples
#' #estimate normal EGAtmfg
#' tmfg <- EGA(data = wmt2[,7:24], model = "TMFG")
#' 
#' #estimate optimal EGAtmfg
#' tmfg.opt <- EGA.fit(data = wmt2[,7:24], model = "TMFG")
#'
#' #estimate Entropy Fit Index
#' entropyFit(data = wmt2[,7:24], structure = tmfg.opt$wc)$Entropy.Fit
#' entropyFit(data = wmt2[,7:24], structure = tmfg$wc)$Entropy.Fit
#'
#'
#' @references
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#' doi:\href{https://doi.org/10.7155/jgaa.00185}{10.7155/jgaa.00185}
#'
#' @seealso \code{\link[EGA]{bootEGA}} to investigate the stability of EGA's estimation via bootstrap,
#' \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA, 
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
EGA.fit <- function (data, model = c("glasso","TMFG"),
                     steps = c(3,4,5,6,7,8))
{
    if(missing(model))
    {model <- "glasso"
    }else{model <- match.arg(model)}
    
    if(missing(steps))
    {steps <- c(3,4,5,6,7,8)
    }else{steps <- steps}
    
    num <- length(steps)
    
    mods <- list()
    dims <- matrix(NA, nrow = ncol(data), ncol = num)
    
    #Generate walktrap models
    for(i in 1:num)
    {
        message(paste("Estimating EGA model",i,"of",num,sep=" "))
        mods[[as.character(steps[i])]] <- EGA(data = data,
                                              model = model,
                                              steps = steps[i],
                                              plot.EGA = FALSE)
        
        dims[,i] <- mods[[as.character(steps[i])]]$wc
    }
    
    colnames(dims) <- as.character(steps)
    
    #check for unique number of dimensions
    uniq.dim <- vector("numeric",length=(ncol(dims)))
    
    for(i in 2:(ncol(dims)))
    {
        uniq.dim[i] <- igraph::compare(dims[,i-1],dims[,i],method = "nmi")
        names(uniq.dim) <- paste(steps,sep="")
    }
    
    uniq <- unique(as.matrix(uniq.dim))
    
    step <- steps[which(uniq!=1)]
    
    len <- length(step)
    
    best.fit <- list()
    
    #if all models are the same
    if(len==1)
    {
        best.fit$EGA <- mods[[1]]
        best.fit$steps <- 4
        message("All EGA models are identical.")
    }else{
        
        ent.vec <- vector("numeric",length=len)
        
        for(i in 1:len)
        {ent.vec[i] <- entropyFit(data, mods[[as.character(step[i])]]$wc)$Entropy.Fit}
        
        names(ent.vec) <- step
        
        best.fit$EGA <- mods[[as.character(step[which(ent.vec==min(ent.vec))])]]
        best.fit$steps <- step[which(ent.vec==min(ent.vec))]
        best.fit$EntropyFit <- ent.vec
        best.fit$Lowest.EntropyFit <- ent.vec[which(ent.vec==min(ent.vec))]
    }
    
    return(best.fit)
    
}
