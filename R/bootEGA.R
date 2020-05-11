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
#' @param algorithm A string indicating the algorithm to use.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#' 
#' \item{\strong{\code{louvain}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_louvain}}}
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
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
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
#' standard deviation, standard error, 95\% confidence intervals, and quantiles (lower = 2.5\% and upper = 97.5\%)}
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
#' @importFrom stats cov median sd qt quantile
#'
#' @export
#'
# Bootstrap EGA
# Updated 11.05.2020
bootEGA <- function(data, n,
                    model = c("glasso", "TMFG"), algorithm = c("walktrap", "louvain"),
                    type = c("parametric", "resampling"),
                    typicalStructure = TRUE, plot.typicalStructure = TRUE, ncores, ...) {
  
  #### MISSING ARGUMENTS HANDLING ####
  
  if(missing(model))
  {model <- "glasso"
  }else{model <- match.arg(model)}
  
  if(missing(algorithm))
  {algorithm <- "walktrap"
  }else{algorithm <- match.arg(algorithm)}
  
  if(missing(type))
  {type <- "parametric"
  }else{type <- match.arg(type)}
  
  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}
  
  #### MISSING ARGUMENTS HANDLING ####
  
  #number of cases
  cases <- nrow(data)
  
  #set inverse covariance matrix for parametric approach
  if(type=="parametric")  # Use a parametric approach:
  {
    g <- -EGA(data, n = cases, model = model, algorithm = algorithm, ...)$network
    diag(g) <- 1
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
  
  #let user know data generation has started
  message("Estimating networks...\n", appendLF = FALSE)
  
  #Parallel processing
  cl <- parallel::makeCluster(ncores)
  
  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "cases", "model", "algorithm", ...),
                          envir=environment())
  
  #Estimate networks
  boots <- pbapply::pblapply(X = datalist, cl = cl,
                             FUN = EGA,
                             model = model,
                             algorithm = algorithm,
                             n = cases,
                             ...)
  
  #Stop parallel processing
  parallel::stopCluster(cl)
  
  #Let user know results are being computed
  message("Computing results...", appendLF = FALSE)
  
  #Obtain networks
  bootGraphs <- lapply(boots, function(x){
    net <- x$network
    colnames(net) <- colnames(data)
    row.names(net) <- colnames(data)
    return(net)
  })
  
  #Obtain communities
  boot.wc <- lapply(boots, function(x){
    wc <- x$wc
    names(wc) <- colnames(data)
    return(wc)
  })
  
  #Obtain number of dimensions
  boot.ndim <- matrix(NA, nrow = n, ncol = 2)
  ndims <- lapply(boots, function(x){x$n.dim})
  boot.ndim[,2] <- unlist(ndims)
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")
  boot.ndim[,1] <- seq_len(n)
  
  if (typicalStructure == TRUE) {
    if(model=="glasso")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, median)
    }else if(model=="TMFG")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, mean)}
    typical.igraph <- NetworkToolbox::convert2igraph(abs(typical.Structure))
    
    typical.wc <- switch(algorithm,
                         walktrap = igraph::cluster_walktrap(typical.igraph),
                         louvain = igraph::cluster_louvain(typical.igraph)
                         )
    
    typical.ndim <- max(typical.wc$membership)
    dim.variables <- data.frame(items = colnames(data), dimension = typical.wc$membership)
  }
  if (plot.typicalStructure == TRUE) {
    plot.typical.ega <- qgraph::qgraph(typical.Structure, layout = "spring",
                                       vsize = 6, groups = as.factor(typical.wc$membership))
  }
  Median <- median(boot.ndim[, 2], na.rm = TRUE)
  se.boot <- sd(boot.ndim[, 2], na.rm = TRUE)
  ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
  ci <- se.boot * ciMult
  quant <- quantile(boot.ndim[,2], c(.025, .975), na.rm = TRUE)
  summary.table <- data.frame(n.Boots = n, median.dim = Median,
                              SE.dim = se.boot, CI.dim = ci,
                              Lower = Median - ci, Upper = Median + ci,
                              Lower.Quantile = quant[1], Upper.Quantile = quant[2])
  row.names(summary.table) <- NULL
  
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
