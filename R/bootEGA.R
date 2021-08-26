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
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
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
#' @param iter Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{500} is recommended
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
#' original dataset, via the \code{\link[MASS]{mvrnorm}} function}
#'
#' \item{\strong{\code{"resampling"}}}
#' {Generates n random subsamples of the original data}
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
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
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
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
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
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#' 
#'  \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#' 
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for 
#' more details and examples}
#' 
#' }
#'
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{iter}{Number of replica samples in bootstrap}
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
#'    }
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \donttest{# bootEGA glasso example
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' boot.wmt <- bootEGA(data = wmt, iter = 500, plot.type = "qgraph",
#' type = "parametric", ncores = 2)
#'
#' # bootEGA TMFG example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, model = "TMFG",
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#'
#' # bootEGA Louvain example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, algorithm = "louvain",
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#'
#' # bootEGA Spinglass example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, model = "TMFG", plot.type = "qgraph",
#' algorithm = igraph::cluster_spinglass, type = "parametric", ncores = 2)
#' }
#'
#' # Load data
#' intwl <- intelligenceBattery[,8:66]
#'
#' \donttest{# Another bootEGA example
#' boot.intwl <- bootEGA(data = intwl, iter = 500,
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#' }
#'
#' @references
#' # Original implementation of bootEGA \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}.
#' \doi{10.31234/osf.io/9deay}
#'
#' # Structural consistency (see \code{\link[EGAnet]{dimensionStability}}) \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#' \doi{10.1002/per.2265}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom stats cov median sd qt quantile
#'
#' @export
#'
# Bootstrap EGA
# Updated 05.08.2021
bootEGA <- function(data, n = NULL, uni.method = c("expand", "LE"), iter,
                    type = c("parametric", "resampling"),
                    corr = c("cor_auto", "pearson", "spearman"),
                    model = c("glasso", "TMFG"), model.args = list(),
                    algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                    typicalStructure = TRUE, plot.typicalStructure = TRUE,
                    plot.type = c("GGally", "qgraph"),
                    plot.args = list(), ncores, ...) {

  #### DEPRECATED ARGUMENTS ####

  # Get additional arguments
  add.args <- list(...)

  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){

    # Give deprecation warning
    warning(
      paste(
        "The 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")",
        sep = ""
      )
    )

    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }
  
  # Check if uni has been input as an argument
  if("uni" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      "The 'uni' argument has been deprecated in all EGA functions."
    )
  }

  #### DEPRECATED ARGUMENTS ####
  
  # Message function
  message(styletext(styletext("\nBootstrap Exploratory Graph Analysis\n", defaults = "underline"), defaults = "bold"))

  #### MISSING ARGUMENTS HANDLING ####
  
  if(missing(uni.method)){
    uni.method <- "LE"
  }else{uni.method <- match.arg(uni.method)}
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}

  if(missing(model)){
    model <- "glasso"
  }else{
    model <- match.arg(model)
  }

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- match.arg(algorithm)
  }

  if(missing(type)){
    type <- "parametric"
  }else{
    type <- match.arg(type)
  }

  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }

  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{
    plot.type <- match.arg(plot.type)
  }

  ## Check for input plot arguments
  color.palette <- "polychrome"

  #### MISSING ARGUMENTS HANDLING ####
  
  # Let user know setting
  message(paste(" \u2022 type = ", type, "\n",
                " \u2022 iterations = ", iter, "\n",
                " \u2022 model = ", model, "\n",
                " \u2022 algorithm = ",
                gsub(
                  "igraph", "",
                  gsub(
                    "::", "",
                    gsub(
                      "cluster_", "",
                      paste(substitute(algorithm), collapse = "")
                    )
                  )
                ),
                "\n",
                " \u2022 correlation = ", corr, "\n",
                " \u2022 unidimensional check = ", ifelse(
                  uni.method == "LE",
                  "leading eigenvalue",
                  "correlation matrix expansion"
                ), "\n",
                sep=""))

  #number of cases
  if(is.null(n)){
    
    if(isSymmetric(as.matrix(data))){
      stop("The argument 'n' is missing for a symmetric matrix")
    }else{
      cases <- nrow(data)
    }
    
  }else{
    cases <- n
  }
  
  #empirical EGA
  empirical.EGA <- suppressMessages(suppressWarnings(EGA(data = data, n = cases, uni.method = uni.method, corr = corr,
                                                         model = model, model.args = model.args,
                                                         algorithm = algorithm, algorith.args = algorithm.args,
                                                         plot.EGA = FALSE)))

  #set inverse covariance matrix for parametric approach
  if(type == "parametric"){  # Use a parametric approach
    
    ## Compute correlation matrix
    cor.data <- empirical.EGA$correlation    

    # Generating data will be continuous
    corr.method <- "pearson"
    
  }else if(type == "resampling"){
    
    # Check if matrix is symmetric
    if(isSymmetric(data)){
      warning("The argument 'data' is symmetric and therefore treated as a correlation matrix. Parametric bootstrap will be used instead")
      type <- "parametric"
      corr.method <- "pearson"
    }else{
      corr.method <- corr
    }
    
  }

  #initialize data list
  datalist <- list()

  #initialize count
  count <- 0

  #let user know data generation has started
  message("Generating data...", appendLF = FALSE)

  repeat{

    #increase count
    count <- count + 1

    #generate data
    if(type == "parametric"){

      datalist[[count]] <- MASS::mvrnorm(cases, mu = rep(0, ncol(cor.data)), Sigma = cor.data)

    }else if(type == "resampling"){

      datalist[[count]] <- data[sample(1:cases, replace=TRUE),]

    }

    #break out of repeat
    if(count == iter)
    {break}
  }

  #let user know data generation has ended
  message("done", appendLF = TRUE)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "uni.method", "cases", "corr",
                                      "model", "model.args",
                                      "algorithm", "algorithm.args"),
                          envir=environment())

  #let user know data generation has started
  message("Estimating EGA networks...\n", appendLF = FALSE)

  #Estimate networks
  boots <- pbapply::pblapply(
    X = datalist, cl = cl,
    FUN = EGA,
    uni.method = uni.method, corr = corr.method,
    model = model, model.args = model.args,
    algorithm = algorithm, algorith.args = algorithm.args,
    plot.EGA = FALSE
  )

  parallel::stopCluster(cl)

  #let user know results are being computed
  message("Computing results...\n")

  #get networks
  bootGraphs <- lapply(boots, function(x, col.names){
    net <- x$network
    colnames(net) <- col.names
    row.names(net) <- col.names
    return(net)
  }, col.names = colnames(data))

  #get community membership
  boot.wc <- lapply(boots, function(x, col.names){
    wc <- x$wc
    names(wc) <- col.names
    return(wc)
  }, col.names = colnames(data))

  #get dimensions
  boot.ndim <- matrix(NA, nrow = iter, ncol = 2)
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")

  boot.ndim[,1] <- seq_len(iter)
  boot.ndim[,2] <- unlist(
    lapply(boots, function(x){
      x$n.dim
    })
  )

  if (typicalStructure){

    typical.Structure <- switch(model,
                                "glasso" = apply(simplify2array(bootGraphs),1:2, median),
                                "TMFG" = apply(simplify2array(bootGraphs),1:2, mean)
                         )

    # Sub-routine to following EGA approach (handles undimensional structures)
    typical.wc <- suppressWarnings(
      suppressMessages(
        
        typicalStructure.network(A = typical.Structure, corr = corr,
                                 model = model, model.args = model.args,
                                 n = cases, uni.method = uni.method, algorithm = algorithm,
                                 algorithm.args = algorithm.args)
        
      )
    )

    typical.ndim <- length(na.omit(unique(typical.wc)))
    
    if(typical.ndim == 1){typical.wc[1:length(typical.wc)] <- 1}
    
    dim.variables <- data.frame(items = colnames(data), dimension = typical.wc)
  }
  
  Median <- median(boot.ndim[, 2], na.rm = TRUE)
  se.boot <- sd(boot.ndim[, 2], na.rm = TRUE)
  ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
  ci <- se.boot * ciMult
  quant <- quantile(boot.ndim[,2], c(.025, .975), na.rm = TRUE)
  summary.table <- data.frame(n.Boots = iter, median.dim = Median,
                              SE.dim = se.boot, CI.dim = ci,
                              Lower.CI = Median - ci, Upper.CI = Median + ci,
                              Lower.Quantile = quant[1], Upper.Quantile = quant[2])
  row.names(summary.table) <- NULL

  #compute frequency
  dim.range <- range(boot.ndim[,2], na.rm = TRUE)
  lik <- matrix(0, nrow = diff(dim.range)+1, ncol = 2)
  colnames(lik) <- c("# of Factors", "Frequency")
  count <- 0

  for(i in seq(from=min(dim.range),to=max(dim.range),by=1))
  {
    count <- count + 1
    lik[count,1] <- i
    lik[count,2] <- length(which(boot.ndim[,2]==i))/iter
  }

  result <- list()
  result$iter <- iter
  result$type <- type
  result$boot.ndim <- boot.ndim
  result$boot.wc <- boot.wc
  result$bootGraphs <- bootGraphs
  result$summary.table <- summary.table
  result$frequency <- lik
  result$EGA <- empirical.EGA

  # Typical structure
  if (typicalStructure) {
    
    typicalGraph <- list()
    typicalGraph$graph <- typical.Structure
    typicalGraph$typical.dim.variables <- dim.variables[order(dim.variables[,2]), ]
    typicalGraph$wc <- typical.wc
    result$typicalGraph <- typicalGraph
    
  }
  
  # Add plot arguments (for itemStability)
  result$color.palette <- color.palette

  class(result) <- "bootEGA"
  
  if(typicalStructure & plot.typicalStructure){
    result$plot.typical.ega <- plot(result)
  }

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

  return(result)
}
#----
