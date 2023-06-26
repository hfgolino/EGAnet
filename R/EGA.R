#' Applies the Exploratory Graph Analysis technique
#'
#' Estimates the number of dimensions of a given dataset or correlation matrix
#' using the graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[EGAnet]{TMFG}})
#' network estimation methods.
#'
#' Two community detection algorithms, Walktrap (Pons & Latapy, 2006) and
#' Louvain (Blondel et al., 2008), are pre-programmed because of their
#' superior performance in simulation studies on psychological
#' data generated from factor models (Christensen & Golino; 2020; Golino et al., 2020).
#' Notably, any community detection algorithm from the \code{\link{igraph}}
#' can be used to estimate the number of communities (see examples).
#'
#' @param data Matrix or data frame.
#' Variables (down columns) or correlation matrix.
#' If the input is a correlation matrix,
#' then argument \code{n} (number of cases) is \strong{required}
#'
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
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
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter}
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
#' \item{\code{most_common_tefi}}
#' {Uses the most common number of communities detected across the number
#' of iterations. After, if there is more than one solution for that number
#' of communities, then the solution with the lowest \code{\link[EGAnet]{tefi}
#' is used}}
#' 
#' }
#'
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
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
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>, Maria Dolores Nieto <acinodam at gmail.com> and Luis E. Garrido <garrido.luiseduardo at gmail.com>
#'
#' @return Returns a list containing:
#'
#' \item{network}{A symmetric network estimated using either the
#' \code{\link{EBICglasso.qgraph}} or \code{\link[EGAnet]{TMFG}}}
#'
#' \item{wc}{A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network}
#'
#' \item{n.dim}{A scalar of how many total dimensions were identified in the network}
#'
#' \item{cor.data}{The zero-order correlation matrix}
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Summary statistics
#' summary(ega.wmt)
#' 
#' # Produce Methods section
#' methods.section(ega.wmt)
#' 
#' # Estimate EGAtmfg
#' ega.wmt.tmfg <- EGA(
#'   data = wmt, model = "TMFG"
#' )
#'
#' # Estimate EGA with Louvain algorithm
#' ega.wmt.louvain <- EGA(
#'   data = wmt, algorithm = "louvain"
#' )
#' 
#' # Estimate EGA with Leiden algorithm
#' ega.wmt.leiden <- EGA(
#'   data = wmt, algorithm = "leiden"
#' )
#'
#' # Estimate EGA with Spinglass algorithm
#' ega.wmt.spinglass <- EGA(
#'   data = wmt,
#'   algorithm = igraph::cluster_spinglass
#' )}
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' # Louvain algorithm \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#'
#' # Comprehensive unidimensionality simulation \cr
#' Christensen, A. P. (2022).
#' Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison.
#' \emph{PsyArXiv}.
#'
#' # Compared all \emph{igraph} community detections algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data \cr
#' # Also implements the Leading Eigenvalue unidimensional method \cr
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2021).
#' Comparing community detection algorithms in psychological data: A Monte Carlo simulation.
#' \emph{PsyArXiv}.
#'
#' # Original simulation and implementation of EGA \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#'
#' Golino, H. F., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#'
#' # Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#'
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
# Main EGA function
# Updated 23.06.2023
EGA <- function (
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    plot.EGA = TRUE, verbose = FALSE,
    ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- ifelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  
  # Ensure data has names
  data <- ensure_dimension_names(data)
  
  # First, obtain the multidimensional result
  ## Proper correlations are automatically computed in `EGA.estimate`
  ## Network estimation settings can be passed to unidimensional result
  multidimensional_result <- EGA.estimate(
    data = data, n = n, corr = corr, na.data = na.data,
    model = model, algorithm = algorithm,
    verbose = verbose, ...
  )
  
  # Store model attributes
  model_attributes <- attr(multidimensional_result$network, "methods")
  
  # Obtain arguments for model
  model_ARGS <- switch(
    model_attributes$model,
    "bggm" = c(
      obtain_arguments(BGGM::estimate, model_attributes),
      obtain_arguments(BGGM:::select.estimate, model_attributes)
    ),
    "glasso" = obtain_arguments(EBICglasso.qgraph, model_attributes),
    "tmfg" = obtain_arguments(TMFG, model_attributes)
  )
  
  # Make adjustments for each model (removes extraneous arguments)
  model_ARGS <- adjust_model_arguments(model_ARGS)
  
  # Set up arguments for unidimensional
  unidimensional_ARGS <- list( # standard arguments
    data = data, n = n, corr = corr, na.data = na.data,
    model = model, uni.method = uni.method,
    verbose = verbose
  )
  
  # `data` at this point will be data or correlation matrix
  # For non-BGGM network estimation, OK to use correlation matrix
  if(model_attributes$model != "bggm"){
    unidimensional_ARGS$data <- multidimensional_result$cor.data
    unidimensional_ARGS$n <- multidimensional_result$n
  }
  
  # Additional arguments for model
  unidimensional_ARGS <- c(unidimensional_ARGS, model_ARGS)
  
  # Third, obtain the unidimensional result
  unidimensional_result <- do.call(
    what = community.unidimensional,
    args = unidimensional_ARGS
  )
  
  # Unidimensional?
  unidimensional <- unique_length(unidimensional_result) == 1
  
  # Determine result
  if(unidimensional){
    multidimensional_result$wc <- unidimensional_result
  }
  
  # Obtain number of dimensions
  multidimensional_result$n.dim <- unique_length(multidimensional_result$wc)
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- data.frame(
    items = dimnames(data)[[2]],
    dimension = as.vector(multidimensional_result$wc)
  )
  
  # Dimension variables data frame by dimension
  multidimensional_result$dim.variables <- dim.variables[
    order(dim.variables$dimension),
  ]
  
  # For legacy, change names of output
  ## Change correlation matrix of `cor.data` to `correlation`
  names(multidimensional_result)[
    names(multidimensional_result) == "cor.data"
  ] <- "correlation"
  
  # Add unidimensional attributes
  attr(multidimensional_result, "unidimensional") <- list(
    uni.method = uni.method, unidimensional = unidimensional
  )
  
  # Check for Louvain method (for consensus information)
  if(uni.method == "louvain"){
    attr(multidimensional_result, "unidimensional")$consensus <-
      attr(unidimensional_result, "methods")
  }
  
  # Set class
  class(multidimensional_result) <- "EGA"
  
  # Check for plot
  if(isTRUE(plot.EGA)){
    
    # Set up plot
    multidimensional_result$Plot.EGA <- plot(multidimensional_result)
    
    # Actually send the plot
    plot(multidimensional_result$Plot.EGA)
    
  }

  # Return EGA
  return(multidimensional_result)
  
}

# Bug checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL; corr = "auto"
# na.data = "pairwise"; model = "glasso"; algorithm = "walktrap"
# uni.method = "louvain"; plot.EGA = TRUE; verbose = FALSE

#' @exportS3Method 
# S3 Print Method ----
# Updated 22.06.2023
print.EGA <- function(x, ...)
{
  
  # Print network estimation
  print(x$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Print community detection
  print(x$wc)
  
  # Add break space
  cat("\n----\n\n")
  
  # Get unidimensional attributes
  unidimensional_attributes <- attr(x, "unidimensional")
  
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
  cat(
    paste0(
      "Unidimensional Method: ", unidimensional_method, "\n",
      "Unidimensional: ", ifelse(
        unidimensional_attributes$unidimensional,
        "Yes", "No"
      )
    )
  )
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 22.06.2023
summary.EGA <- function(object, ...)
{
  
  # Print network estimation
  print(object$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Print community detection
  print(object$wc)
  
  # Add break space
  cat("\n----\n\n")
  
  # Get unidimensional attributes
  unidimensional_attributes <- attr(object, "unidimensional")
  
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
  cat(
    paste0(
      "Unidimensional Method: ", unidimensional_method, "\n",
      "Unidimensional: ", ifelse(
        unidimensional_attributes$unidimensional,
        "Yes", "No"
      )
    )
  )
  
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 22.06.2023
plot.EGA <- function(x, ...)
{
  
  # Return plot
  single_plot(
    network = x$network,
    wc = x$wc,
    ...
  )
  
}


#' @noRd
# Cleaning adjusts model arguments
# Updated 23.06.2023
adjust_model_arguments <- function(model_ARGS)
{
  
  # Arguments to set to NULL
  null_arguments <- c(
    # All already exist in `unidimensional_ARGS`
    "data", "n", "corr", "na.data", "model", "verbose",
    # EBICglasso-specific arguments
    "penalizeMatrix",
    # BGGM-specific arguments
    "Y", "object"
  )
  
  # Set intersecting arguments to NULL
  model_ARGS[
    intersect(names(model_ARGS), null_arguments)
  ] <- NULL
  
  # Return model arguments
  return(model_ARGS)
  
}







