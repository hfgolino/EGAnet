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
#' \item{\code{most_common_tefi}}
#' {Uses the most common number of communities detected across the number
#' of iterations. After, if there is more than one solution for that number
#' of communities, then the solution with the lowest \code{\link[EGAnet]{tefi}
#' is used}}
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
#' \dontrun{
#' # Standard EGA example
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 500,
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Produce Methods section
#' methods.section(boot.wmt)
#' 
#' # Louvain example
#' boot.wmt.louvain <- bootEGA(
#'   data = wmt, iter = 500,
#'   algorithm = "louvain",
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Spinglass example
#' boot.wmt.spinglass <- bootEGA(
#'   data = wmt, iter = 500,
#'   algorithm = igraph::cluster_spinglass, # use any function from {igraph}
#'   type = "parametric", ncores = 2
#' )
#'
#' # EGA fit example
#' boot.wmt.fit <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "EGA.fit",
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Hierarchical EGA example
#' boot.wmt.hier <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "hierEGA",
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Random-intercept EGA example
#' boot.wmt.ri <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "riEGA",
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
# Updated 30.06.2023
bootEGA <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    iter = 500, type = c("parametric", "resampling"),
    ncores, EGA.type = c("EGA", "EGA.fit", "hierEGA", "riEGA"),
    typicalStructure = TRUE, plot.typicalStructure = TRUE,
    verbose = TRUE, seed = 1234,
    ...
) 
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- ifelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  type <- set_default(type, "parametric", bootEGA)
  EGA.type <- set_default(EGA.type, "ega", bootEGA)
  
  # Set cores
  if(missing(ncores)){ncores <- round(parallel::detectCores() / 2)}
  
  # `EGA.estimate` will handle legacy arguments and data processing 
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Get data dimensions
  dimensions <- dim(data)
  
  # Ensure data has names
  data <- ensure_dimension_names(data)
  
  # Get variable names
  variable_names <- dimnames(data)[[2]]

  # NEED HANDLING FOR DIFFERENT TYPES ONLY EGA RIGHT NOW

  # Obtain EGA function
  ega_function <- switch(
    EGA.type,
    "ega" = EGA,
    "ega.fit" = EGA.fit,
    "hierega" = hierEGA,
    "riega" = riEGA
  )
  
  # Estimate empirical EGA
  empirical_EGA <- ega_function(
    data = data, n = n, corr = corr,
    na.data = na.data, model = model,
    algorithm = algorithm, uni.method = uni.method,
    plot.EGA = FALSE, verbose = FALSE, ...
  )
  # If "n" is NULL and a correlation matrix is supplied,
  # then an error will be thrown within `EGA.estimate`,
  # so no need to check for "n"
  
  # Generate data
  bootstrap_data <- reproducible_bootstrap(
    data = data, samples = iter, cases = empirical_EGA$n,
    mu = rep(0, dimensions[2]), Sigma = empirical_EGA$correlation,
    seed = seed, type = type
  )
  
  # Perform bootstrap using parallel processing
  boots <- parallel_process(
    # Parallel processing arguments
    iterations = iter,
    datalist = bootstrap_data,
    ncores = ncores, progress = verbose,
    # Standard EGA arguments
    FUN = ega_function, n = n, corr = corr,
    na.data = na.data, model = model,
    algorithm = algorithm, uni.method = uni.method,
    plot.EGA = FALSE, verbose = FALSE,
    ...
  )
  
  # NEED HANDLING FOR DIFFERENT TYPES ONLY EGA RIGHT NOW
  # EGA.fit = empirical_EGA$EGA
  # riEGA = empirical_EGA$EGA
  # hierEGA = will have "higher" and "lower"
  
  # Get results
  results <- prepare_bootEGA_results(boots)
  
  # Organize results
  result <- list(
    iter = iter, type = type, boot.ndim = results$boot_n.dim,
    boot.wc = results$boot_memberships, bootGraphs = results$boot_networks,
    summary.table = results$summary_table, frequency = results$frequencies,
    EGA = empirical_EGA, EGA.type = EGA.type
  )
  
  # No attributes needed (all information is contained in output)
  
  # Set class
  class(result) <- "bootEGA"
  
  # Check for typical structure results
  if(isTRUE(typicalStructure)){
    
    # Obtain results
    result$typicalGraph <- estimate_typicalStructure(
      data, results, empirical_EGA, verbose, ...
    )
    
    # Check for plot
    if(isTRUE(plot.typicalStructure)){
      
      # Get plot
      result$plot.typical.ega <- plot(result)
      
      # Actually send plot
      plot(result$plot.typical.ega)
      
    }
    
  }
  
  # Return result
  return(result)
  
}

# Bug checking ----
# DATA
# data = wmt2[,7:24]; n = NULL; corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"; uni.method = "louvain"
# iter = 100; type = "parametric"; ncores = 8; EGA.type = "EGA"
# typicalStructure = TRUE; plot.typicalStructure = TRUE;
# verbose = TRUE; seed = 1234
# r_sample_seeds <- EGAnet:::r_sample_seeds
# r_sample_with_replacement <- EGAnet:::r_sample_with_replacement
# r_sample_without_replacemetn <- EGAnet:::r_sample_without_replacement
# Need above functions for testing!

#' @noRd
# Typical network and memberships ----
# Updated 30.06.2023
estimate_typicalStructure <- function(
    data, results, empirical_EGA, verbose, ...
)
{

  # Get ellipse arguments
  ellipse <- list(...)
  
  # Attributes from empirical EGA
  model_attributes <- attr(empirical_EGA$network, "methods")
  algorithm_attributes <- attr(empirical_EGA$wc, "methods")
  unidimensional_attributes <- attr(empirical_EGA, "unidimensional")
  
  # Set model, algorithm and unidimensional method
  model <- tolower(model_attributes$model)
  algorithm <- tolower(algorithm_attributes$algorithm)
  uni.method <- tolower(unidimensional_attributes$uni.method)
  
  # Get network
  network <- switch(
    model,
    "bggm" = apply(simplify2array(results$boot_networks),1:2, median),
    "glasso" = apply(simplify2array(results$boot_networks),1:2, median),
    "tmfg" = apply(simplify2array(results$boot_networks),1:2, mean)
  )
  
  # Check for function or non-Louvain method
  if(
    is.function(algorithm) ||
    !algorithm %in% c("louvain", "signed_louvain")
  ){
    
    # Apply non-Louvain method
    wc <- do.call(
      what = community.detection,
      args = c(
        list(
          network = network, algorithm = algorithm,
          membership.only = TRUE
        ),
        ellipse # pass on ellipse
      )
    )
    
  }else{ # for Louvain, use consensus clustering
    
    # Check for consensus method
    ellipse$consensus.method <- ifelse(
      "consensus.method" %in% names(ellipse),
      ellipse$consensus.method, "most_common" # default
    )
    
    # Check for consensus iterations
    ellipse$consensus.iter <- ifelse(
      "consensus.iter" %in% names(ellipse),
      ellipse$consensus.iter, 1000 # default
    )
    
    # Apply consensus clustering
    wc <- do.call(
      what = community.consensus,
      args = c(
        list(
          network = network,
          signed = algorithm == "signed_louvain",
          membership.only = TRUE,
          verbose = FALSE
        ),
        ellipse # pass on ellipse
      )
    )
    
  }
  
  # Obtain arguments for model
  model_ARGS <- switch(
    model,
    "bggm" = c(
      obtain_arguments(BGGM::estimate, model_attributes),
      obtain_arguments(BGGM:::select.estimate, model_attributes)
    ),
    "glasso" = obtain_arguments(EBICglasso.qgraph, model_attributes),
    "tmfg" = obtain_arguments(TMFG, model_attributes)
  )
  
  # Make adjustments for each model (removes extraneous arguments)
  model_ARGS <- adjust_model_arguments(model_ARGS)
  # Found in `EGA` internals
  
  # Set up arguments for unidimensional
  unidimensional_ARGS <- list( # standard arguments
    data = data, n = empirical_EGA$n, 
    corr = model_attributes$corr, 
    na.data = model_attributes$na.data,
    model = model, uni.method = uni.method,
    verbose = FALSE
  )
  
  # `data` at this point will be data or correlation matrix
  # For non-BGGM network estimation, OK to use correlation matrix
  if(model_attributes$model != "bggm"){
    unidimensional_ARGS$data <- empirical_EGA$correlation
  }
  
  # Additional arguments for model
  unidimensional_ARGS <- c(unidimensional_ARGS, model_ARGS)
  
  # Third, obtain the unidimensional result
  unidimensional_result <- do.call(
    what = community.unidimensional,
    args = unidimensional_ARGS
  )
  
  # Determine result
  if(unique_length(unidimensional_result) == 1){
    wc <- unidimensional_result
  }
  
  # Obtain number of dimensions
  n.dim <- unique_length(wc)
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- fast.data.frame(
    c(dimnames(empirical_EGA$network)[[2]], wc),
    nrow = dim(empirical_EGA$network)[2], ncol = 2,
    colnames = c("items", "dimension")
  )
  
  # Dimension variables data frame by dimension
  dim.variables <- dim.variables[
    order(dim.variables$dimension),
  ]
  
  # Return results
  return(
    list(
      graph = network,
      typical.dim.variables = dim.variables,
      wc = wc, n.dim = n.dim
    )
  )

}

#' @noRd
# Prepare `bootEGA` results ----
# Self-contained to work on `EGA` bootstraps
# Updated 29.06.2023
prepare_bootEGA_results <- function(boot_object)
{
  
  # Get number of iterations
  iter <- length(boot_object)
  
  # Get network
  boot_networks <- lapply(boot_object, function(x){
    return(x$network)
  })
  
  # Get memberships
  boot_memberships <- t(nnapply(boot_object, function(x){
    return(x$wc)
  }, LENGTH = dim(boot_networks[[1]])[2]))
  
  # Get bootstrap dimensions
  boot_n.dim <- nnapply(boot_object, function(x){x$n.dim})
  
  # Set up data frame
  boot_dimensions <- fast.data.frame(
    c(seq_len(iter), boot_n.dim),
    nrow = iter, ncol = 2,
    colnames = c("Boot.Number", "N.Dim")
  )
  
  # Pre-compute standard deviation and confidence interval
  median_boot <- median(boot_n.dim, na.rm = TRUE)
  se_boot <- sd(boot_n.dim, na.rm = TRUE)
  ci_boot <- se_boot * qt(0.95 / 2 + 0.5, iter - 1)
  quantile_boot <- quantile(boot_n.dim, c(0.025, 0.975), na.rm = TRUE)
  
  # Set up descriptive statistics
  summary_table <- fast.data.frame(
    c(
      iter, median_boot, se_boot, ci_boot, 
      median_boot - ci_boot, median_boot + ci_boot,
      quantile_boot
    ), ncol = 8,
    colnames = c(
      "n.Boots", "median.dim", "SE.dim", "CI.dim",
      "Lower.CI", "Upper.CI", "Lower.Quantile", "Upper.Quantile"
    ), row.names = ""
  )
  
  # Compute frequency
  frequencies <- count_table(boot_n.dim, proportion = TRUE)
  dimnames(frequencies)[[2]] <- c("# of Factors", "Frequency")
  
  # Return results
  return(
    list(
      boot_networks = boot_networks,
      boot_memberships = boot_memberships,
      boot_n.dim = boot_n.dim,
      summary_table = summary_table,
      frequencies = frequencies
    )
  )
  
  
}

