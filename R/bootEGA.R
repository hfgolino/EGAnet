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
#' Defaults to \code{"louvain"}.
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
#' {Applies the Leading Eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' \item{\strong{\code{louvain}}}
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix using a resolution parameter = 0.95.
#' If the number of dimensions is 1, then the Louvain solution is used; otherwise,
#' regular EGA is used. This method was validated in the Christensen (2022) simulation.}
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
#' @param seed Numeric.
#' Seed to reproduce results. Defaults to \code{1234}. For random results, set to \code{NULL}
#'
#' @param corr Character.
#' Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' 
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
#' @param EGA.type Character.
#' Type of EGA model to use.
#' 
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\code{\link[EGAnet]{EGA}}}
#' {Uses standard exploratory graph analysis}
#' 
#' \item{\code{\link[EGAnet]{EGA.fit}}}
#' {Uses \code{\link[EGAnet]{tefi}} to determine best fit of
#' \code{\link[EGAnet]{EGA}}}
#' 
#' \item{\code{\link[EGAnet]{hierEGA}}}
#' {Uses hierarhical exploratory graph analysis}
#' 
#' \item{\code{\link[EGAnet]{riEGA}}}
#' {Uses random-intercept exploratory graph analysis}
#' 
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
#' or \code{\link[EGAnet]{TMFG}}
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
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}.
#' Defaults to \code{objective_function = "modularity"}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#' 
#' @param consensus.iter Numeric.
#' Number of iterations to perform in consensus clustering for the Louvain algorithm
#' (see Lancichinetti & Fortunato, 2012).
#' Defaults to \code{100}
#' 
#' @param consensus.method Character.
#' What consensus clustering method should be used? 
#' Defaults to \code{"highest_modularity"}.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\strong{\code{highest_modularity}}}
#' {Uses the community solution that achieves the highest modularity
#' across iterations}
#' 
#' \item{\strong{\code{most_common}}}
#' {Uses the community solution that is found the most
#' across iterations}
#' 
#' \item{\strong{\code{iterative}}}
#' {Identifies the most common community solutions across iterations
#' and determines how often nodes appear in the same community together.
#' A threshold of 0.30 is used to set low proportions to zero.
#' This process repeats iteratively until all nodes have a proportion of
#' 1 in the community solution.
#' }
#' 
#' \item{\code{lowest_tefi}}
#' {Uses the community solution that achieves the lowest \code{\link[EGAnet]{tefi}}
#' across iterations}
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
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' See \code{\link[GGally]{ggnet2}} for
#' full list of arguments:
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
#' @param progress Boolean.
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' For Windows, \code{FALSE} is about 2x faster
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
#' \donttest{# Standard EGA example
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Produce Methods section
#' methods.section(boot.wmt)
#' 
#' # Louvain example
#' boot.wmt.louvain <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   algorithm = "louvain",
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Spinglass example
#' boot.wmt.spinglass <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   algorithm = igraph::cluster_spinglass, # use any function from {igraph}
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )
#'
#' # EGA fit example
#' boot.wmt.fit <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   EGA.type = "EGA.fit",
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Hierarchical EGA example
#' boot.wmt.hier <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   EGA.type = "hierEGA",
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Random-intercept EGA example
#' boot.wmt.ri <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   EGA.type = "riEGA",
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )}
#'
#' @references
#' # Original implementation of bootEGA \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#'
#' # Structural consistency (see \code{\link[EGAnet]{dimensionStability}}) \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom stats cov median sd qt quantile
#' @importFrom methods formalArgs
#'
#' @export
#'
# Bootstrap EGA
# Updated 28.08.2022
bootEGA <- function(
    data, n = NULL, uni.method = c("expand", "LE", "louvain"), iter,
    type = c("parametric", "resampling"), seed = 1234,
    corr = c("cor_auto", "pearson", "spearman"),
    EGA.type = c("EGA", "EGA.fit", "hierEGA", "riEGA"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "leiden", "louvain"), algorithm.args = list(),
    consensus.method = c(
      "highest_modularity",
      "most_common",
      "iterative",
      "lowest_tefi"
    ), consensus.iter = 100,
    typicalStructure = TRUE, plot.typicalStructure = TRUE,
    plot.args = list(), ncores, progress = TRUE, ...
) 
{
  
  # Make data a matrix
  data <- as.matrix(data)
  
  #### DEPRECATED ARGUMENTS
  
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
  
  #### MISSING ARGUMENTS HANDLING
  
  if(missing(uni.method)){
    uni.method <- "LE"
  }else{uni.method <- match.arg(uni.method)}
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  if(missing(EGA.type)){
    EGA.type <- "EGA"
  }else{EGA.type <- match.arg(EGA.type)}
  
  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}
  
  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- match.arg(algorithm)
  }else{
    algorithm <- algorithm
  }
  
  if(missing(type)){
    type <- "parametric"
  }else{type <- match.arg(type)}
  
  if(missing(consensus.method)){
    consensus.method <- "most_common"
  }else{consensus.method <- tolower(match.arg(consensus.method))}
  
  
  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }
  
  ## Check for input plot arguments
  if("color.palette" %in% names(plot.args)){
    color.palette <- plot.args$color.palette
  }else{color.palette <- "polychrome"}
  
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
  
  # Obtain EGA function
  ega_function <- switch(
    tolower(EGA.type),
    "ega" = EGA,
    "ega.fit" = EGA.fit,
    "hierega" = hierEGA,
    "riega" = riEGA
  )
  
  # EGA calls
  ega_args <- list(
    data = data, n = cases, uni.method = uni.method,
    corr = corr, model = model, model.args = model.args,
    algorithm = algorithm, algorithm.args = algorithm.args,
    consensus.method = consensus.method,
    consensus.iter = consensus.iter,
    plot.EGA = FALSE
  )
  
  # Remove calls that are not in formal arguments
  ega_args <- ega_args[
    !is.na(
      match(names(ega_args), formalArgs(ega_function))
    )
  ]
  
  # Estimate empirical EGA (run quietly)
  empirical_EGA <- suppressWarnings(
    suppressMessages(
      do.call(
        ega_function, ega_args
      )
    )
  )
  
  # EGA output
  if("EGA" %in% names(empirical_EGA)){
    # Obtain EGA output from EGA option 
    ega_output <- empirical_EGA$EGA
  }else if(!is(empirical_EGA, "EGA")){
    
    # Base on EGA type
    if(tolower(EGA.type) == "hierega"){
      
      # Obtain lower order output
      ega_output <- empirical_EGA$hierarchical$lower_order
      
    }
    
  }else{
    # Assume output is from standard EGA
    ega_output <- empirical_EGA
  }
  
  
  #set inverse covariance matrix for parametric approach
  if(type == "parametric"){  # Use a parametric approach
    
    ## Compute correlation matrix
    cor.data <- ega_output$correlation
    
    # Generating data will be continuous
    corr.method <- "pearson"
    
  }else if(type == "resampling"){
    
    # Check if matrix is symmetric
    if(isSymmetric(as.matrix(data))){
      warning("The argument 'data' is symmetric and therefore treated as a correlation matrix. Parametric bootstrap will be used instead")
      type <- "parametric"
      corr.method <- "pearson"
    }else{
      corr.method <- corr
    }
    
  }
  
  # Set seed
  set.seed(seed)
  
  if(is.null(seed)){
    warning("Results are unique. Set the 'seed' argument for reproducible results (see examples)")
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
      
      datalist[[count]] <- MASS_mvrnorm(cases, mu = rep(0, ncol(cor.data)), Sigma = cor.data)
      
    }else if(type == "resampling"){
      
      datalist[[count]] <- data[sample(1:cases, replace=TRUE),]
      
    }
    
    #break out of repeat
    if(count == iter)
    {break}
  }
  
  #let user know data generation has ended
  message("done", appendLF = TRUE)
  
  # Perform bootstrap using parallel processing
  # See in `utils-EGAnet`
  boots <- parallel_process(
    datalist = datalist,
    iter = iter,
    progress = progress,
    FUN = ega_function,
    FUN_args = ega_args,
    ncores = ncores
  )
  
  # Bootstrap output
  if("EGA" %in% names(boots[[1]])){
    # Obtain EGA output from EGA option 
    boot_output <- lapply(boots, function(x){
      x$EGA
    })
  }else if(!is(boots[[1]], "EGA")){
    
    # Base on EGA type
    if(tolower(EGA.type) == "hierega"){
      
      # Obtain lower order output
      boot_output <- lapply(boots, function(x){
        x$hierarchical$lower_order
      })
      boot_output_higher <- lapply(boots, function(x){
        x$hierarchical$higher_order$EGA
      })
      
    }
    
  }else{
    # Assume output is from standard EGA
    boot_output <- boots
  }
  
  #let user know results are being computed
  message("Computing results...\n")
  
  #get networks
  bootGraphs <- lapply(boot_output, function(x, col.names){
    net <- x$network
    colnames(net) <- col.names
    row.names(net) <- col.names
    return(net)
  }, col.names = colnames(data))
  
  #get community membership
  boot.wc <- lapply(boot_output, function(x, col.names){
    wc <- x$wc
    names(wc) <- col.names
    return(wc)
  }, col.names = colnames(data))
  
  #get dimensions
  boot.ndim <- matrix(NA, nrow = iter, ncol = 2)
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")
  
  boot.ndim[,1] <- seq_len(iter)
  boot.ndim[,2] <- unlist(
    lapply(boot_output, function(x){
      x$n.dim
    })
  )
  
  if (typicalStructure){
    
    typical.Structure <- switch(
      model,
      "glasso" = apply(simplify2array(bootGraphs),1:2, median),
      "TMFG" = apply(simplify2array(bootGraphs),1:2, mean)
      
    )
    
    # Sub-routine to following EGA approach (handles unidimensional structures)
    typical.wc <- suppressWarnings(
      suppressMessages(
        
        typicalStructure.network(
          A = typical.Structure, corr = corr,
          model = model, model.args = model.args,
          n = cases, uni.method = uni.method, algorithm = algorithm,
          algorithm.args = algorithm.args,
          consensus.method = consensus.method,
          consensus.iter = consensus.iter
        )
        
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
  
  for(i in seq(from=min(dim.range),to=max(dim.range),by=1)){
    count <- count + 1
    lik[count,1] <- i
    lik[count,2] <- length(which(boot.ndim[,2]==i))/iter
  }
  
  # Higher order EGA
  if(tolower(EGA.type) == "hierega"){
    
    #get networks
    bootGraphs_higher <- lapply(boot_output_higher, function(x){
      net <- x$network
      return(net)
    })
    
    #get community membership
    boot.wc_higher <- lapply(boot_output_higher, function(x){
      wc <- x$wc
      return(wc)
    })
    
    #get dimensions
    boot.ndim_higher <- matrix(NA, nrow = iter, ncol = 2)
    colnames(boot.ndim_higher) <- c("Boot.Number", "N.Dim")
    
    boot.ndim_higher[,1] <- seq_len(iter)
    boot.ndim_higher[,2] <- unlist(
      lapply(boot_output_higher, function(x){
        x$n.dim
      })
    )
    
    if (typicalStructure){
      
      typical.Structure_higher <- switch(model,
                                         "glasso" = apply(simplify2array(bootGraphs_higher),1:2, median),
                                         "TMFG" = apply(simplify2array(bootGraphs_higher),1:2, mean)
      )
      
      # Sub-routine to following EGA approach (handles undimensional structures)
      typical.wc_higher <- suppressWarnings(
        suppressMessages(
          
          typicalStructure.network(
            A = typical.Structure_higher, corr = corr,
            model = model, model.args = model.args,
            n = cases, uni.method = uni.method, algorithm = algorithm,
            algorithm.args = algorithm.args,
            consensus.method = consensus.method,
            consensus.iter = consensus.iter
          )
          
        )
      )
      
      typical.ndim_higher <- length(na.omit(unique(typical.wc_higher)))
      
      if(typical.ndim_higher == 1){typical.wc_higher[1:length(typical.wc_higher)] <- 1}
      
      dim.variables_higher <- data.frame(items = names(typical.wc_higher), dimension = typical.wc_higher)
    }
    
    Median_higher <- median(boot.ndim_higher[, 2], na.rm = TRUE)
    se.boot_higher <- sd(boot.ndim_higher[, 2], na.rm = TRUE)
    ciMult_higher <- qt(0.95/2 + 0.5, nrow(boot.ndim_higher) - 1)
    ci_higher <- se.boot_higher * ciMult_higher
    quant_higher <- quantile(boot.ndim_higher[,2], c(.025, .975), na.rm = TRUE)
    summary.table_higher <- data.frame(n.Boots = iter, median.dim = Median_higher,
                                       SE.dim = se.boot_higher, CI.dim = ci_higher,
                                       Lower.CI = Median_higher - ci_higher, Upper.CI = Median_higher + ci_higher,
                                       Lower.Quantile = quant_higher[1], Upper.Quantile = quant_higher[2])
    row.names(summary.table_higher) <- NULL
    
    #compute frequency
    dim.range_higher <- range(boot.ndim_higher[,2], na.rm = TRUE)
    lik_higher <- matrix(0, nrow = diff(dim.range_higher)+1, ncol = 2)
    colnames(lik_higher) <- c("# of Factors", "Frequency")
    count <- 0
    
    for(i in seq(from=min(dim.range_higher),to=max(dim.range_higher),by=1)){
      count <- count + 1
      lik_higher[count,1] <- i
      lik_higher[count,2] <- length(which(boot.ndim_higher[,2]==i))/iter
    }
    
  }
  
  # Reset seed
  set.seed(NULL)
  
  # Set up result list
  if(tolower(EGA.type) == "hierega"){
    result <- list(
      iter = iter, type = type, boot.ndim = boot.ndim,
      boot.wc = boot.wc, bootGraphs = bootGraphs,
      summary.table = summary.table, frequency = lik,
      EGA = empirical_EGA, EGA.type = EGA.type
    )
  }else{
    result <- list(
      iter = iter, type = type, boot.ndim = boot.ndim,
      boot.wc = boot.wc, bootGraphs = bootGraphs,
      summary.table = summary.table, frequency = lik,
      EGA = ega_output, EGA.type = EGA.type
    )
  }
  
  # Typical structure
  if (typicalStructure) {
    
    typicalGraph <- list(
      graph = typical.Structure,
      typical.dim.variables = dim.variables[order(dim.variables[,2]), ],
      wc = typical.wc
    )
    
    result$typicalGraph <- typicalGraph
    
  }
  
  # higher order
  if(tolower(EGA.type) == "hierega"){
    
    # Set up result list
    result_higher <- list(
      iter = iter, type = type, boot.ndim = boot.ndim_higher,
      boot.wc = boot.wc_higher, bootGraphs = bootGraphs_higher,
      summary.table = summary.table_higher, frequency = lik_higher
    )
    
    # Typical structure
    if (typicalStructure) {
      
      typicalGraph_higher <- list(
        graph = typical.Structure_higher,
        typical.dim.variables = dim.variables_higher[order(dim.variables_higher[,2]), ],
        wc = typical.wc_higher
      )
      
      result_higher$typicalGraph <- typicalGraph_higher
      
    }
    
  }
  
  # Add plot arguments (for itemStability)
  result$color.palette <- color.palette
  
  class(result) <- "bootEGA"
  
  if(tolower(EGA.type) == "hierega"){
    class(result_higher) <- "bootEGA"
  }
  
  if(typicalStructure & plot.typicalStructure){
    
    if(tolower(EGA.type) != "hierega"){
      result$plot.typical.ega <- plot(
        result,
        plot.args = plot.args
      )
    }else{
      
      # Set up output
      hier_plot <- suppressMessages(
        suppressWarnings(
          suppressPackageStartupMessages(
            ggpubr::ggarrange(
              plot(
                result,
                plot.args = plot.args,
                produce = FALSE
              ), # plot lower-order
              plot(
                result_higher,
                plot.args = plot.args,
                produce = FALSE
              ), # plot higher-order
              labels = c("Lower Order", "Higher Order")
            )
          )
        )
      )
      
      # Output plots
      result$plot.typical.ega <- suppressMessages(
        suppressWarnings(
          suppressPackageStartupMessages(
            plot(hier_plot)
          )
        )
      )
      
    }
    
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
  
  # Set results for higher order EGA
  if(tolower(EGA.type) == "hierega"){
    
    result <- list(
      result_lower = result,
      result_higher = result_higher
    )
    
  }
  
  return(result)
  
}
#----