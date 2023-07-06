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
#' @param data Numeric matrix or data frame.
#' Either data representing \emph{only} the variables of interest or
#' a correlation matrix. Data that are not numeric will be
#' removed from the dataset
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} is a correlation matrix
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"} to automatically compute
#' appropriate correlations using \code{\link[EGAnet]{auto.correlate}}.
#' \code{"pearson"} and \code{"spearman"} are provide for completeness.
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"}}
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"}}
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"}}
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"}}
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"}}
#' {Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details}
#' 
#' }
#' 
#' @param algorithm Character or \code{\link{igraph}} \code{cluster_*} function.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#' 
#' \itemize{
#'
#' \item{\code{"leiden"}}
#' {See \code{\link[igraph]{cluster_leiden}} for more details}
#' 
#' \item{\code{"louvain"}}
#' {By default, \code{"louvain"} will implement the non-signed version
#' of the Louvain algorithm using the consensus clustering method 
#' (see \code{\link[EGAnet]{community.consensus}} for more information). 
#' This function will implement \code{consensus.method = "most_common"}
#' and \code{consensus.iter = 1000} unless specified otherwise}
#' 
#' \item{\code{"walktrap"}}
#' {This algorithm is the default. See \code{\link[EGAnet]{cluster_walktrap}} for more details}
#' 
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
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#' 
#' @param verbose Boolean.
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}}, \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}}, and
#' \code{\link[EGAnet]{community.consensus}}
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
# Updated 05.07.2023
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
  dim.variables <- fast.data.frame(
    c(dimnames(data)[[2]], as.vector(multidimensional_result$wc)),
    nrow = length(multidimensional_result$wc), ncol = 2,
    colnames = c("items", "dimension")
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
    multidimensional_result$Plot.EGA <- plot(
      multidimensional_result, ...
    )
    
    # Actually send the plot
    silent_plot(multidimensional_result$Plot.EGA)
    
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
# Updated 05.07.2023
summary.EGA <- function(object, ...)
{
  print(object, ...) # same as print
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







