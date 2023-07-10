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
# Updated 10.07.2023
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
    verbose = TRUE, # seed = 1234,
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
    "riega" = riEGA,
    "hierega" = stop("'EGA.type = \"hierEGA\"' is not yet supported. Support coming soon...", call. = FALSE),
    stop(
      paste0(
        "EGA.type = \"", EGA.type, "\" is not supported. Please use one of the default options: ", 
        paste0("\"", as.character(formals(bootEGA)$EGA.type)[-1], "\"", collapse = ", ")
      )
    )
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
  
  # Obtain EGA output
  empirical_EGA_output <- get_EGA_object(empirical_EGA)
  
  # Generate data
  bootstrap_data <- reproducible_bootstrap(
    data = data, samples = iter, cases = empirical_EGA_output$n,
    mu = rep(0, dimensions[2]), Sigma = empirical_EGA_output$correlation,
    type = type # , seed = seed
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
  
  # Obtain bootstrap EGA output
  bootstrap_EGA_output <- lapply(boots, get_EGA_object)
  
  # Get results
  results <- prepare_bootEGA_results(bootstrap_EGA_output, iter)
  
  # Add additional results
  results[c("type", "EGA", "EGA.type")] <- list(
    type, empirical_EGA, EGA.type
  )
  
  # No attributes needed (all information is contained in output)
  
  # Set class
  class(results) <- "bootEGA"
  
  # Check for typical structure results
  if(isTRUE(typicalStructure)){
    
    # Obtain results
    results$typicalGraph <- estimate_typicalStructure(
      data, results, verbose, ...
    )
    
    # Check for plot
    if(isTRUE(plot.typicalStructure)){
      
      # Get plot
      results$plot.typical.ega <- plot(results, ...)
      
      # Actually send plot
      silent_plot(results$plot.typical.ega)
      
    }
    
  }
  
  # Return result
  return(results)
  
}

# Bug checking ----
# DATA
# data = wmt2[,7:24]; n = NULL; corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"; uni.method = "louvain"
# iter = 100; type = "parametric"; ncores = 8; EGA.type = "EGA"
# typicalStructure = TRUE; plot.typicalStructure = FALSE;
# verbose = TRUE; seed = 1234
# r_sample_seeds <- EGAnet:::r_sample_seeds
# r_sample_with_replacement <- EGAnet:::r_sample_with_replacement
# r_sample_without_replacement <- EGAnet:::r_sample_without_replacement
# Need above functions for testing!

#' @exportS3Method 
# S3 Print Method ----
# Updated 06.07.2023
print.bootEGA <- function(x, ...)
{
  
  # Ensure proper EGA object
  ega_object <- get_EGA_object(x)
  
  # Print network information
  send_network_methods(ega_object$network, boot = TRUE)
  
  # Add line break
  cat("\n")
  
  # Print community detection
  print(ega_object$wc, boot = TRUE)
  
  # Add line break
  cat("\n")
  
  # Do not print unidimensional for `EGA.fit`
  if(x$EGA.type != "ega.fit"){
    
    # Get unidimensional attributes
    unidimensional_attributes <- attr(ega_object, "unidimensional")
    
    # Obtain unidimensional method
    unidimensional_method <- switch(
      unidimensional_attributes$uni.method,
      "expand" = "Expand",
      "le" = "Leading Eigenvector",
      "louvain" = "Louvain"
    )
    
    # Set up unidimensional print
    if(unidimensional_method == "Louvain"){
      
      # Set up consensus attributes
      consensus_attributes <- unidimensional_attributes$consensus
      
      # Obtain consensus name
      consensus_name <- switch(
        consensus_attributes$consensus.method,
        "highest_modularity" = "Highest Modularity",
        "iterative" = "Iterative",
        "most_common" = "Most Common",
        "lowest_tefi" = "Lowest TEFI"
      )
      
      # Update unidimensional method text
      unidimensional_method <- paste0(
        unidimensional_method, " (", consensus_name,
        " for ", consensus_attributes$consensus.iter,
        " iterations)"
      )
      
    }
    
    # Print unidimensional
    cat("Unidimensional Method: ", unidimensional_method)
    
  }
  
  # Add break space
  cat("\n\n----\n\n")
  
  # Set proper EGA type name
  ega_type <- switch(
    x$EGA.type,
    "ega" = "EGA",
    "ega.fit" = "EGA.fit",
    "hierega" = "hierEGA",
    "riega" = "riEGA"
  )
  
  # Print EGA type
  cat(paste0("EGA Type: ", ega_type), "\n")
  
  # Set up methods
  cat(
    paste0(
      "Bootstrap Samples: ",
      x$iter, " (", totitle(x$type), ")\n"
    )
  )
  
  # Print frequency table
  frequency_df <- as.data.frame(
    do.call(rbind, lapply(x$frequency, as.character))
  )
  # Adjust dimension names (`dimnames` doesn't work)
  colnames(frequency_df) <- NULL
  row.names(frequency_df) <- c("", "Frequency: ")
  # Finally, print
  print(frequency_df)
  
  # Print summary table
  cat(
    paste0(
      "\nMedian dimensions: ", x$summary.table$median.dim,
      " [", round(x$summary.table$Lower.CI, 2), ", ",
      round(x$summary.table$Upper.CI, 2), "] 95% CI"
    )
  )

}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 05.07.2023
summary.bootEGA <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 05.07.2023
plot.bootEGA <- function(x, ...)
{
  
  # Return plot
  single_plot(
    network = x$typicalGraph$graph,
    wc = x$typicalGraph$wc,
    ...
  )
  
}

#' @noRd
# Prepare `bootEGA` results ----
# Self-contained to work on `EGA` bootstraps
# Updated 06.07.2023
prepare_bootEGA_results <- function(boot_object, iter)
{
  
  # Get networks
  boot_networks <- lapply(boot_object, function(x){x$network})
  
  # Get memberships
  boot_memberships <- t(nvapply(boot_object, function(x){x$wc}, LENGTH = dim(boot_networks[[1]])[2]))
  
  # Get bootstrap dimensions
  boot_n.dim <- nvapply(boot_object, function(x){x$n.dim})
  
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
  
  # Set up descriptive statistics
  summary_table <- fast.data.frame(
    c(
      iter, median_boot, se_boot, ci_boot, 
      median_boot - ci_boot, median_boot + ci_boot,
      quantile(boot_n.dim, c(0.025, 0.975), na.rm = TRUE)
    ), ncol = 8,
    colnames = c(
      "n.Boots", "median.dim", "SE.dim", "CI.dim",
      "Lower.CI", "Upper.CI", "Lower.Quantile", "Upper.Quantile"
    ), row.names = ""
  )
  
  # Compute frequency
  frequencies <- count_table(boot_n.dim, proportion = TRUE)
  dimnames(frequencies)[[2]] <- c("# of Factors", "Frequency")
  
  # Return results (keep legacy names)
  return(
    list(
      iter = iter,
      bootGraphs = boot_networks,
      boot.wc = boot_memberships,
      boot.ndim = boot_n.dim,
      summary.table = summary_table,
      frequency = frequencies
    )
  )
  
}

#' @noRd
# Typical Walktrap `EGA.fit` ----
# Updated 07.07.2023
typical_walktrap_fit <- function(network, dimensions, ellipse)
{
  
  # Check for parameter search space
  if(!"steps" %in% names(ellipse)){
    steps <- 3:8 # default
  }else{
    
    # Set steps
    steps <- ellipse$steps
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "steps"]
    
  }
  
  # Get the length of the steps
  step_length <- length(steps)
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    steps, function(step){
      do.call(
        what = community.detection,
        args = c(
          list( # Necessary call
            network = network,
            algorithm = "walktrap",
            steps = step
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    steps, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = steps
    )
  )
  
}

#' @noRd
# Typical Leiden `EGA.fit` ----
# Updated 07.07.2023
typical_leiden_fit <- function(network, dimensions, ellipse)
{
  
  # Check for objective function
  if(!"objective_function" %in% names(ellipse)){
    objective_function <- "CPM"
    # default for {igraph} and `community.detection`
  }else{
    
    # Set objective function
    objective_function <- ellipse$objective_function
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "objective_function"]
    
  }
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    
    # Switch based on objective function
    if(objective_function == "CPM"){
      resolution_parameter <- seq.int(0, max(abs(network)), 0.01) # default 
    }else if(objective_function == "modularity"){
      resolution_parameter <- seq.int(0, 2, 0.05) # default
    }
    
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Get the length of the resolution parameter
  resolution_parameter_length <- length(resolution_parameter)
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    resolution_parameter, function(resolution_parameter){
      do.call(
        what = community.detection,
        args = c(
          list( # Necessary call
            network = network,
            algorithm = "leiden",
            resolution_parameter = resolution_parameter,
            objective_function = objective_function
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    resolution_parameter, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = resolution_parameter
    )
  )
  
}

#' @noRd
# Typical Louvain `EGA.fit` ----
# Updated 07.07.2023
typical_louvain_fit <- function(network, dimensions, ellipse)
{
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    resolution_parameter <- seq.int(0, 2, 0.05) # default
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Get the length of the resolution parameter
  resolution_parameter_length <- length(resolution_parameter)
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    resolution_parameter, function(resolution_parameter){
      do.call(
        what = community.consensus,
        args = c(
          list( # Necessary call
            network = network,
            resolution = resolution_parameter
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    resolution_parameter, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = resolution_parameter
    )
  )
  
}

#' @noRd
# Typical `EGA.fit` memberships ----
# Needs to be handled consistently
# Updated 07.07.2023
estimate_typical_EGA.fit <- function(results, ellipse)
{
  
  # Get proper EGA object
  ega_object <- get_EGA_object(results)
  
  # Attributes from empirical EGA
  model_attributes <- attr(ega_object$network, "methods")
  algorithm_attributes <- attr(ega_object$wc, "methods")
  
  # Set model, algorithm and unidimensional method
  model <- tolower(model_attributes$model)
  algorithm <- tolower(algorithm_attributes$algorithm)
  
  # Get network
  network <- switch(
    model,
    "bggm" = symmetric_matrix_lapply(results$bootGraphs, median),
    "glasso" = symmetric_matrix_lapply(results$bootGraphs, median),
    "tmfg" = symmetric_matrix_lapply(results$bootGraphs, mean)
  )
  
  # Make sure proper names are there
  dimnames(network) <- dimnames(ega_object$network)
  
  # Get network dimensions
  dimensions <- dim(network)
  
  # Branch for fit function
  fit_FUN <- switch(
    algorithm,
    "walktrap" = typical_walktrap_fit,
    "leiden" = typical_leiden_fit,
    "louvain" = typical_louvain_fit
  )
  
  # Get fit result
  fit_result <- fit_FUN(network, dimensions, ellipse)
  
  # Set objective function (will be NULL for Louvain and Walktrap)
  objective_function <- fit_result$objective_function
  
  # Obtain only unique solutions
  search_unique <- unique_solutions(fit_result$search_matrix)
  
  # Determine best fitting solution
  fit_values <- apply(
    search_unique, 1, function(membership){
      tefi(ega_object$correlation, membership)
    }$VN.Entropy.Fit
  )
  
  # Determine best fit index
  best_index <- which.min(fit_values)
  
  # Obtain best solution
  best_solution <- search_unique[best_index,]
  
  # Obtain signed argument
  if(!"signed" %in% names(ellipse)){
    signed <- FALSE # default
  }else{signed <- ellipse$signed}
  
  # Add methods to membership attributes
  attr(best_solution, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm),
    # `obtain_algorithm_name` is with `community.detection`
    signed = signed, objective_function = objective_function
  )
  
  # Set class (proper `EGA.estimate` printing)
  class(best_solution) <- "EGA.community"
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- fast.data.frame(
    c(dimnames(ega_object$network)[[2]], best_solution),
    nrow = dimensions[2], ncol = 2,
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
      wc = best_solution, n.dim = unique_length(best_solution),
      EntropyFit = fit_values,
      Lowest.EntropyFit = fit_values[best_index],
      parameter.space = fit_result$parameters
    )
  )
  
}

#' @noRd
# Typical network and memberships ----
# Updated 06.07.2023
estimate_typicalStructure <- function(
    data, results, verbose, ...
)
{

  # Get ellipse arguments
  ellipse <- list(...)
  
  # If results are from `EGA.fit`, handle separately
  if(results$EGA.type == "ega.fit"){
    return(estimate_typical_EGA.fit(results, ellipse))
  }else{ # Get proper EGA object
    ega_object <- get_EGA_object(results)
  }
  
  # Attributes from empirical EGA
  model_attributes <- attr(ega_object$network, "methods")
  algorithm_attributes <- attr(ega_object$wc, "methods")
  unidimensional_attributes <- attr(ega_object, "unidimensional")
  
  # Set model, algorithm and unidimensional method
  model <- tolower(model_attributes$model)
  algorithm <- tolower(algorithm_attributes$algorithm)
  uni.method <- tolower(unidimensional_attributes$uni.method)
  
  # Get network
  network <- switch(
    model,
    "bggm" = symmetric_matrix_lapply(results$bootGraphs, median),
    "glasso" = symmetric_matrix_lapply(results$bootGraphs, median),
    "tmfg" = symmetric_matrix_lapply(results$bootGraphs, mean)
  )
  
  # Make sure proper names are there
  dimnames(network) <- dimnames(ega_object$network)
  
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
    if(!"consensus.method" %in% names(ellipse)){
      ellipse$consensus.method <- "most_common" # default
    }
    
    # Check for consensus iterations
    if(!"consensus.iter" %in% names(ellipse)){
      ellipse$consensus.iter <- 1000 # default
    }

    # Apply consensus clustering
    wc <- do.call(
      what = community.consensus,
      args = c(
        list(
          network = network,
          signed = algorithm == "signed_louvain",
          membership.only = TRUE
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
      overwrite_arguments(
        # defaults for `BGGM:::select.estimate`
        list(cred = 0.95, alternative = "two.sided"), model_attributes
      )
    ),
    "glasso" = obtain_arguments(EBICglasso.qgraph, model_attributes),
    "tmfg" = obtain_arguments(TMFG, model_attributes)
  )
  
  # Make adjustments for each model (removes extraneous arguments)
  model_ARGS <- adjust_model_arguments(model_ARGS)
  # Found in `EGA` internals
  
  # Set up arguments for unidimensional
  unidimensional_ARGS <- list( # standard arguments
    data = data, n = ega_object$n, 
    corr = model_attributes$corr, 
    na.data = model_attributes$na.data,
    model = model, uni.method = uni.method,
    verbose = FALSE
  )
  
  # `data` at this point will be data or correlation matrix
  # For non-BGGM network estimation, OK to use correlation matrix
  if(model_attributes$model != "bggm"){
    unidimensional_ARGS$data <- ega_object$correlation
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
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- fast.data.frame(
    c(dimnames(ega_object$network)[[2]], wc),
    nrow = dim(ega_object$network)[2], ncol = 2,
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
      wc = wc, n.dim = unique_length(wc)
    )
  )

}

