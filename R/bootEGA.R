#' Dimension Stability Analysis of \code{\link[EGAnet]{EGA}}
#'
#' \code{bootEGA} Estimates the number of dimensions of \emph{n} bootstraps
#' using the empirical (partial) correlation matrix (parametric) or resampling from
#' the empirical dataset (non-parametric). It also estimates a typical
#' median network structure, which is formed by the median or mean pairwise (partial)
#' correlations over the \emph{n} bootstraps.
#'
#' @param data Matrix or data frame.
#' Includes the variables to be used in the \code{bootEGA} analysis
#'
#' @param n Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{500} is recommended
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"glasso"}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}}}
#'
#' \item{\strong{\code{"TMFG"}}}
#' {Estimates a Triangulated Maximally Filtered Graph.
#' See \code{\link[NetworkToolbox]{TMFG}}}
#'
#' }
#'
#' @param type Character.
#' A string indicating the type of bootstrap to use.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"parametric"}}}
#' {Generates \code{n} new datasets (multivariate normal random distributions) based on the
#' original dataset, via the \code{\link[mvtnorm]{Mvnorm}} function of the mvtnorm package}
#'
#' \item{\strong{\code{"resampling"}}}
#' {Generates n random subsamples of the original data}
#'
#' }
#'
#' @param typicalStructure Boolean.
#' If \code{TRUE}, returns the typical network of partial correlations
#' (estimated via graphical lasso or via TMFG) and estimates its dimensions.
#' The "typical network" is the median of all pairwise correlations over the \emph{n} bootstraps.
#' Defaults to \code{TRUE}
#'
#' @param plot.typicalStructure Boolean.
#' If \code{TRUE}, returns a plot of the typical network (partial correlations),
#' which is the median of all pairwise correlations over the \emph{n} bootstraps,
#' and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{4}.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#' 
#' @param ... Additional arguments to be passed to \code{\link{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @return Returns a list containing:
#'
#' \item{n}{Number of replica samples in bootstrap}
#'
#' \item{boot.ndim}{Number of dimensions identified in each replica sample}
#'
#' \item{boot.wc}{Item allocation for each replica sample}
#'
#' \item{bootGraphs}{Networks of each replica sample}
#'
#' \item{summary.table}{Summary table containing number of replica samples, median,
#' standard deviation, standard error, and 95\% confidence intervals}
#'
#' \item{frequency}{Proportion of times the number of dimensions was identified
#' (e.g., .85 of 1,000 = 850 times that specific number of dimensions was found)}
#'
#' \item{EGA}{Output of the original \code{\link[EGAnet]{EGA}} results}
#'
#' \item{typicalGraph}{A list containing:
#'
#' \itemize{
#'
#' \item{\strong{\code{graph}}}
#' {Network matrix of the median network structure}
#'
#' \item{\strong{\code{typical.dim.variables}}}
#' {An ordered matrix of item allocation}
#'
#' \item{\strong{\code{wc}}}
#' {Item allocation of the median network}
#'
#'     }
#' }
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#'
#' # bootEGA glasso example
#' boot.wmt <- bootEGA(data = wmt, n = 500, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso", type = "parametric", ncores = 4)
#' }
#'
#' # Load data
#' intwl <- intelligenceBattery[,8:66]
#'
#' \dontrun{
#' # bootEGA TMFG example
#' boot.intwl <- bootEGA(data = intelligenceBattery[,8:66], n = 500, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "TMFG", type = "parametric", ncores = 4)
#'
#' }
#'
#' @references
#' Christensen, A. P., & Golino, H. F. (2019).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/9deay}{10.31234/osf.io/9deay}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom stats cov median sd qt
#'
#' @export
#'
# Bootstrap EGA
bootEGA <- function(data, n,
                    model = c("glasso", "TMFG"), type = c("parametric", "resampling"),
                    typicalStructure = TRUE, plot.typicalStructure = TRUE, ncores = 4, ...) {
    
    #number of cases
    cases <- nrow(data)
    
    #set inverse covariance matrix for parametric approach
    if(type=="parametric")  # Use a parametric approach:
    {
        if(model=="glasso")
        {
            g <- -EBICglasso.qgraph(qgraph::cor_auto(data), n = cases, lambda.min.ratio = 0.1, returnAllResults = FALSE, ...)
            diag(g) <- 1
        }else if(model=="TMFG")
        {
            g <- -NetworkToolbox::LoGo(data, normal = TRUE, partial=TRUE, ...)
            diag(g) <- 1
        }
    }
    
    #initialize data list
    datalist <- list()
    
    #initialize count
    count <- 0
    
    #let user know data generation has started
    message("\nGenerating data...", appendLF = FALSE)
    
    repeat{
        
        #increase count
        count <- count + 1
        
        #generate data
        if(type == "parametric")
        {datalist[[count]] <- mvtnorm::rmvnorm(cases, sigma = corpcor::pseudoinverse(g))
        }else if(type == "resampling")
        {datalist[[count]] <- data[sample(1:cases, replace=TRUE),]}
        
        #break out of repeat
        if(count == n)
        {break}
    }
    
    #let user know data generation has ended
    message("done", appendLF = TRUE)
    
    #initialize correlation matrix list
    corlist <- list()
    
    #let user know data generation has started
    message("\nComputing correlation matrices...\n", appendLF = FALSE)
    
    #Parallel processing
    cl <- parallel::makeCluster(ncores)
    
    #Export variables
    parallel::clusterExport(cl = cl,
                            varlist = c("datalist", "corlist", "cases", ...),
                            envir=environment())
    
    #Compute correlation matrices
    corlist <- pbapply::pblapply(X = datalist, cl = cl,
                                 FUN = qgraph::cor_auto)
    
    #let user know data generation has started
    message("Estimating networks...\n", appendLF = FALSE)
    
    #Estimate networks
    if(model == "glasso")
    {
        boots <- pbapply::pblapply(X = corlist, cl = cl,
                                   FUN = EBICglasso.qgraph,
                                   n = cases,
                                   lambda.min.ratio = 0.1,
                                   returnAllResults = FALSE,
                                   ...)
    }else if(model == "TMFG")
    {
        boots <- pbapply::pblapply(X = corlist, cl = cl,
                                   FUN = NetworkToolbox::TMFG,
                                   normal = TRUE,
                                   ...)
        
        for(i in 1:n)
        {boots[[i]] <- boots[[i]]$A}
    }
    
    parallel::stopCluster(cl)
    
    #let user know results are being computed
    message("Computing results...", appendLF = FALSE)
    
    bootGraphs <- vector("list", n)
    for (i in 1:n) {
        bootGraphs[[i]] <- boots[[i]]
        colnames(bootGraphs[[i]]) <- colnames(data)
        rownames(bootGraphs[[i]]) <- colnames(data)
    }
    boot.igraph <- vector("list", n)
    for (l in 1:n) {
        boot.igraph[[l]] <- NetworkToolbox::convert2igraph(abs(bootGraphs[[l]]))
    }
    boot.wc <- vector("list", n)
    for (m in 1:n) {
        boot.wc[[m]] <- igraph::walktrap.community(boot.igraph[[m]])
    }
    boot.ndim <- matrix(NA, nrow = n, ncol = 2)
    for (m in 1:n) {
        boot.ndim[m, 2] <- max(boot.wc[[m]]$membership)
    }
    
    colnames(boot.ndim) <- c("Boot.Number", "N.Dim")
    
    boot.ndim[, 1] <- seq_len(n)
    if (typicalStructure == TRUE) {
        if(model=="glasso")
        {typical.Structure <- apply(simplify2array(bootGraphs),1:2, median)
        }else if(model=="TMFG")
        {typical.Structure <- apply(simplify2array(bootGraphs),1:2, mean)}
        typical.igraph <- NetworkToolbox::convert2igraph(abs(typical.Structure))
        typical.wc <- igraph::walktrap.community(typical.igraph)
        typical.ndim <- max(typical.wc$membership)
        dim.variables <- data.frame(items = colnames(data), dimension = typical.wc$membership)
    }
    if (plot.typicalStructure == TRUE) {
        plot.typical.ega <- qgraph::qgraph(typical.Structure, layout = "spring",
                                           vsize = 6, groups = as.factor(typical.wc$membership))
    }
    Median <- median(boot.ndim[, 2])
    se.boot <- 1.253 * sd(boot.ndim[, 2])
    ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
    ci <- se.boot * ciMult
    summary.table <- data.frame(n.Boots = n, median.dim = Median,
                                SE.dim = se.boot, CI.dim = ci,
                                Lower = Median - ci, Upper = Median + ci)
    
    #compute frequency
    dim.range <- range(boot.ndim[,2])
    lik <- matrix(0, nrow = diff(dim.range)+1, ncol = 2)
    colnames(lik) <- c("# of Factors", "Frequency")
    count <- 0
    
    for(i in seq(from=min(dim.range),to=max(dim.range),by=1))
    {
        count <- count + 1
        lik[count,1] <- i
        lik[count,2] <- length(which(boot.ndim[,2]==i))/n
    }
    
    #let user know results have been computed
    message("done", appendLF = TRUE)
    
    result <- list()
    result$n <- n
    result$boot.ndim <- boot.ndim
    result$boot.wc <- boot.wc
    result$bootGraphs <- bootGraphs
    result$summary.table <- summary.table
    result$frequency <- lik
    result$EGA <- suppressMessages(suppressWarnings(EGA(data = data, model = model, plot.EGA = FALSE)))
    
    # Typical structure
    if (typicalStructure == TRUE) {
        typicalGraph <- list()
        typicalGraph$graph <- typical.Structure
        typicalGraph$typical.dim.variables <- dim.variables[order(dim.variables[,2]), ]
        typicalGraph$wc <- typical.wc$membership
        result$typicalGraph <- typicalGraph
    }
    
    class(result) <- "bootEGA"
    return(result)
}
#----
