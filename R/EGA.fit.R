#' \code{\link[EGAnet]{EGA}} Optimal Model Fit using the Total Entropy Fit Index  (\code{\link[EGAnet]{tefi}})
#'
#' @description Estimates the best fitting model using \code{\link[EGAnet]{EGA}}.
#' The number of steps in the \code{\link[igraph]{cluster_walktrap}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link[EGAnet]{tefi}}.
#'
#' @param data Matrix or data frame.
#' Dataset or correlation matrix
#' 
#' @param n Integer.
#' Sample size (if the data provided is a correlation matrix)
#' 
#' @param uni.method Character.
#' What unidimensionality method should be used? 
#' Defaults to \code{"LE"}.
#' Current options are:
#' 
#' \itemize{
#'
#' \item{\strong{\code{expand}}}
#' {Expands the correlation matrix with four variables correlated .50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This is the method used in the Golino et al. (2020)
#' \emph{Psychological Methods} simulation.}
#'
#' \item{\strong{\code{LE}}}
#' {Applies the leading eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the leading eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' }
#' 
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}
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
#' See \code{\link[EGAnet]{TMFG}}}
#'
#' }
#' 
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Defaults to \code{"walktrap"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}} or \code{\link[igraph]{cluster_leiden}}.
#' Options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{steps}}}
#' {Number of steps used in the Walktrap algorithm. Defaults to \code{c(3:8)}}
#'
#' \item{\strong{\code{leiden}}}
#' {Resolution parameter used in the Leiden algorithm. Defaults to \code{seq(0, 2, .001)}.
#' Higher values lead to smaller communities, lower values lead to larger communities}
#'
#' }
#'
#' @return Returns a list containing:
#'
#' \item{EGA}{The \code{\link[EGAnet]{EGA}} output for the best fitting model}
#'
#' \item{steps}{The number of steps used in the best fitting model from
#' the \code{\link[igraph]{cluster_walktrap}} algorithm}
#' 
#' \item{resolution_parameter}{The resolution parameter used in the best fitting model from
#' the \code{\link[igraph]{cluster_leiden}} algorithm}
#'
#' \item{EntropyFit}{The \code{\link[EGAnet]{tefi}} Index for the unique solutions given the range of steps
#' (vector names represent the number of steps)}
#'
#' \item{Lowest.EntropyFit}{The lowest value for the \code{\link[EGAnet]{tefi}} Index}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt, 
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Estimate optimal EGA
#' fit.wmt <- EGA.fit(data = wmt)
#' 
#' # Plot optimal fit
#' plot(fit.wmt$EGA)
#'
#' # Estimate CFAs
#' cfa.ega <- CFA(ega.wmt, estimator = "WLSMV", data = wmt)
#' cfa.fit <- CFA(fit.wmt$EGA, estimator = "WLSMV", data = wmt)
#'
#' # Compare CFAs
#' lavaan::lavTestLRT(
#'   cfa.ega$fit, cfa.fit$fit,
#'   method = "satorra.bentler.2001"
#' )}
#'
#' @references
#' # Entropy fit measures \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (in press).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' 
#' # Simulation for EGA.fit \cr
#' Jamison, L., Christensen, A. P., & Golino, H. (under review).
#' Optimizing Walktrap's community detection in networks using the Total Entropy Fit Index.
#' \emph{PsyArXiv}.
#' 
#' # Leiden algorithm \cr
#' Traag, V. A., Waltman, L., & Van Eck, N. J. (2019).
#' From Louvain to Leiden: guaranteeing well-connected communities.
#' \emph{Scientific Reports}, \emph{9}(1), 1-12.
#' 
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @seealso \code{\link[EGAnet]{bootEGA}} to investigate the stability of EGA's estimation via bootstrap,
#' \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA,
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# EGA fit
# Updated 11.22.2021
EGA.fit <- function (
  data, n = NULL, uni.method = c("expand", "LE"),
  corr = c("cor_auto", "pearson", "spearman"),
  model = c("glasso","TMFG"),
  algorithm = c("leiden", "walktrap"),
  algorithm.args = list(
    steps = c(3:8),
    resolution_parameter = seq(0, 2, .001)
  )
)
{
  if(missing(uni.method)){
    uni.method <- "LE"
  }else{uni.method <- match.arg(uni.method)}
  
  # Check if uni.method = "LE" has been used
  if(uni.method == "LE"){
    # Give change warning
    warning(
      paste(
        "Previous versions of EGAnet (<= 0.9.8) checked unidimensionality using",
        styletext('uni.method = "expand"', defaults = "underline"),
        "as the default"
      )
    )
  }else if(uni.method == "expand"){
    # Give change warning
    warning(
      paste(
        "Newer evidence suggests that",
        styletext('uni.method = "LE"', defaults = "underline"),
        'is more accurate than uni.method = "expand" (see Christensen, Garrido, & Golino, 2021 in references).',
        '\n\nIt\'s recommended to use uni.method = "LE"'
      )
    )
  }
  
  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else{algorithm <- tolower(match.arg(algorithm))}
  
  # Check for algorithm arguments
  if(algorithm == "walktrap"){
    if(!"steps" %in% names(algorithm.args)){
      algorithm.args <- list()
      algorithm.args$steps <- c(3:8)
    }
  }else if(algorithm == "leiden"){
    if(!"resolution_parameter" %in% names(algorithm.args)){
      algorithm.args <- list()
      algorithm.args$resolution_parameter <- seq(0, 2, .001)
    }
  }
  
  #Speed up process with data
  if(nrow(data) != ncol(data)){
    n <- nrow(data)
    data <- switch(corr,
                   "cor_auto" = auto.correlate(data),
                   "pearson" = cor(data, use = "pairwise.complete.obs"),
                   "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
    )
  }
  
  # Check algorithm
  if(algorithm == "walktrap"){
    
    best.fit <- list()
    
    steps <- algorithm.args$steps
    
    num <- length(steps)
    
    mods <- list()
    dims <- matrix(NA, nrow = ncol(data), ncol = num)
    
    #Generate walktrap models
    for(i in 1:num)
    {
      message(paste("Estimating EGA -- Walktrap model",i,"of",num,sep=" "))
      mods[[as.character(steps[i])]] <- EGA(data = data,
                                            n = n,
                                            model = model,
                                            model.args = list(steps = steps[i]),
                                            algorithm = "walktrap",
                                            plot.EGA = FALSE)
      
      dims[,i] <- mods[[as.character(steps[i])]]$wc
    }
    
    colnames(dims) <- as.character(steps)
    
    #remove solutions with missing dimensions
    rm.cols <- which(apply(apply(dims, 2, is.na), 2, any))
    
    if(length(rm.cols) != 0)
    {
      dims <- dims[,-rm.cols]
      steps <- steps[-rm.cols]
    }
    
    #check for unique number of dimensions
    step <- as.numeric(colnames(dims)[which(!duplicated(homogenize.membership(dims[,1], dims), MARGIN = 2))])
    
    len <- length(step)
    
    #if all models are the same
    if(len==1)
    {
      best.fit$EGA <- mods[[1]]
      best.fit$steps <- 4
      Sys.sleep(1)
      message("\nAll EGA models are identical.")
      Sys.sleep(1)
    }else{
      
      ent.vec <- vector("numeric",length=len)
      
      for(i in 1:len)
      {ent.vec[i] <- tefi(abs(mods[[as.character(step[i])]]$correlation), mods[[as.character(step[i])]]$wc)$VN.Entropy.Fit}
      
      names(ent.vec) <- step
      
      best.fit$EGA <- mods[as.character(step[which(ent.vec==min(ent.vec))])]
      best.fit$steps <- step[which(ent.vec==min(ent.vec))]
      best.fit$EntropyFit <- ent.vec
      best.fit$Lowest.EntropyFit <- ent.vec[which(ent.vec==min(ent.vec))]
    }
    
    # Get information for EGA Methods section
    args <- list()
    
    args$model <- model
    args$algorithm <- "walktrap"
    args$steps <- range(steps)
    args$entropy <- best.fit$Lowest.EntropyFit
    args$solutions <- best.fit$EntropyFit
    
    best.fit$Methods <- args
    
  }else if(algorithm == "leiden"){
    
    best.fit <- list()
    
    resolution_parameter <- algorithm.args$resolution_parameter
    
    num <- length(resolution_parameter)
    
    mods <- list()
    dims <- matrix(NA, nrow = ncol(data), ncol = num)
    
    #Generate Leiden models
    ## Estimate EGA
    ega <- suppressMessages(
      suppressWarnings(
        EGA(data = data,
            n = n,
            model = model,
            plot.EGA = FALSE)
      )
    )
    
    ## Obtain network
    net <- ega$network
    
    ## Convert to igraph
    g <- convert2igraph(abs(net))
    
    ## Estimate Leiden results
    results <- list()
    
    for(i in 1:length(resolution_parameter)){
      
      mods[[as.character(resolution_parameter[i])]] <- igraph::cluster_leiden(
        g,
        resolution_parameter = resolution_parameter[i]
        
      )
      
      dims[,i] <- mods[[as.character(resolution_parameter[i])]]$membership
      
    }
    
    colnames(dims) <- as.character(resolution_parameter)
    
    #remove solutions with missing dimensions
    rm.cols <- which(apply(apply(dims, 2, is.na), 2, any))
    
    if(length(rm.cols) != 0)
    {
      dims <- dims[,-rm.cols]
      resolution_parameter <- resolution_parameter[-rm.cols]
    }
    
    #check for unique number of dimensions
    resolution <- as.numeric(colnames(dims)[which(!duplicated(homogenize.membership(dims[,1], dims), MARGIN = 2))])
    
    len <- length(resolution)
    
    #if all models are the same
    if(len==1)
    {
      best.fit$EGA <- mods[[1]]
      best.fit$resolution <- resolution_parameter[[1]]
      Sys.sleep(1)
      message("\nAll EGA models are identical.")
      Sys.sleep(1)
    }else{
      
      ent.vec <- vector("numeric",length=len)
      
      for(i in 1:len)
      {ent.vec[i] <- tefi(abs(ega$correlation), mods[[as.character(resolution[i])]]$membership)$VN.Entropy.Fit}
      
      names(ent.vec) <- resolution
      
      best.fit$EGA <- suppressMessages(
        suppressWarnings(
          EGA(data = data,
              n = n,
              model = model,
              algorithm = igraph::cluster_leiden,
              algorithm.args = list(
                resolution_parameter = resolution[which(ent.vec==min(ent.vec))]
              ),
              plot.EGA = FALSE)
        )
      )
      best.fit$resolution_parameter <- resolution[which(ent.vec==min(ent.vec))]
      best.fit$EntropyFit <- ent.vec
      best.fit$Lowest.EntropyFit <- ent.vec[which(ent.vec==min(ent.vec))]
    }
    
    # Get information for EGA Methods section
    args <- list()
    
    args$model <- model
    args$algorithm <- "leiden"
    args$steps <- range(resolution)
    args$entropy <- best.fit$Lowest.EntropyFit
    args$solutions <- best.fit$EntropyFit
    
    best.fit$Methods <- args
    
  }
      
  class(best.fit) <- "EGA.fit"
      
  return(best.fit)
}
#----
