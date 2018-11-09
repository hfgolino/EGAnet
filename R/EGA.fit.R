#' EGA Optimal Model Fit
#' @description Estimates the best fitting model using \code{\link{EGA}}.
#' The number of steps in the \code{\link[igraph]{walktrap.community}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link{entropyFit}}.
#'
#' @param data A dataset
#'
#' @param model A string indicating the method to use. Current options are:
#' -\code{glasso}:
#' {Gaussian graphical model estimation using graphical LASSO with extended Bayesian information criterion to select optimal regularization parameter (default method). Using \code{\link[qgraph]{EBICglasso}} from the \code{\link[qgraph]{qgraph}} package.}
#' \code{TMFG}:
#' {Estimates a Triangulated Maximally Filtered Graph, using the function \code{\link[NetworkToolbox]{TMFG}} of the \code{\link[NetworkToolbox]{NetworkToolbox}} package}
#'
#' @param steps Range of steps to be used in the model selection.
#' Defaults from 3 to 8 steps (based on Pons & Latapy, 2006)
#'
#' @return Returns a list containing the \code{\link{EGA}} output
#' and the number of steps found to provide the optimal solution
#'
#' @examples
#' tmfg <- EGA.fit(data = wmt2[,7:24], model = "TMFG")
#'
#' entropyFit(data = wmt2[,7:24], structure = tmfg$wc)
#'
#' entropyFit(data = wmt2[,7:24], structure = EGA(data = wmt2[,7:24], model = "TMFG")$wc)
#'
#' @references
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#' doi:\href{https://doi.org/10.7155/jgaa.00185}{10.7155/jgaa.00185}
#'
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
EGA.fit <- function (data, model = c("glasso","TMFG"),
                     steps = c(3,4,5,6,7,8))
{
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
    uniq.dim <- vector("numeric",length=(ncol(dims)-1))

    for(i in 1:(ncol(dims)-1))
    {
        uniq.dim[i] <- igraph::compare(dims[,i],dims[,i+1],method = "nmi")
        names(uniq.dim)[i] <- paste(steps[i+1],sep="")
    }

    uniq <- unique(as.matrix(uniq.dim))

    if(uniq[1,1]!=1)
    {step <- c(3,as.numeric(row.names(uniq)))
    }else{step <- as.numeric(row.names(uniq))}

    len <- length(step)

    #if all models are the same
    if(len==1)
    {
        best.fit <- mods[[1]]
        best.fit$steps <- 4
        message("All EGA models are identical.")
    }else{

        ent.vec <- vector("numeric",length=len)

        for(i in 1:len)
        {ent.vec[i] <- entropyFit(data, mods[[as.character(step[i])]]$wc)$Adj.Entropy}

        best.fit <- mods[[as.character(step[which(ent.vec==min(ent.vec))])]]
        best.fit$steps <- step[which(ent.vec==min(ent.vec))]
        best.fit$EntropyFit <- ent.vec
        best.fit$Lowest.EntropyFit <- ent.vec[which(ent.vec==min(ent.vec))]
    }

    return(best.fit)

}
