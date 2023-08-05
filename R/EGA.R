#' @title Exploratory Graph Analysis
#'
#' @description Estimates the number of communities (dimensions) of
#' a dataset or correlation matrix using a network estimation method
#' (Golino & Epskamp, 2017; Golino et al., 2020). After, a community
#' detection algorithm is applied (Christensen et al., 2023) for
#' multidimensional data. A unidimensional check is also applied
#' based on findings from Golino et al. (2020) and Christensen (2023) 
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} provided is a correlation matrix
#' 
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"auto"} --- }
#' {Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)}
#' 
#' \item{\code{"pearson"} --- }
#' {Pearson's correlation is computed for all variables regardless of
#' categories}
#' 
#' \item{\code{"spearman"} --- }
#' {Spearman's rank-order correlation is computed for all variables
#' regardless of categories}
#' 
#' }
#' 
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"} --- }
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"} --- }
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"} --- }
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"} --- }
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"} --- }
#' {Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details}
#' 
#' }
#' 
#' @param algorithm Character or 
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#' 
#' \itemize{
#'
#' \item{\code{"leiden"} --- }
#' {See \code{\link[igraph]{cluster_leiden}} for more details}
#' 
#' \item{\code{"louvain"} --- }
#' {By default, \code{"louvain"} will implement the Louvain algorithm using 
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}} 
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000} 
#' unless specified otherwise}
#' 
#' \item{\code{"walktrap"} --- }
#' {See \code{\link[igraph]{cluster_walktrap}} for more details}
#' 
#' }
#'
#' @param uni.method Character (length = 1).
#' What unidimensionality method should be used? 
#' Defaults to \code{"louvain"}.
#' Available options:
#' 
#' \itemize{
#'
#' \item{\code{"expand"} --- }
#' {Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation}
#'
#' \item{\code{"LE"} --- }
#' {Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023) 
#' \emph{Behavior Research Methods} simulation}
#' 
#' \item{\code{"louvain"} --- }
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1, 
#' then the Louvain solution is used; otherwise, regular EGA is used. 
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}}
#' 
#' }
#'
#' @param plot.EGA Boolean (length = 1).
#' Defaults to \code{TRUE}.
#' Whether the plot should be returned with the results.
#' Set to \code{FALSE} for no plot
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}}, 
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}}, 
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{community.unidimensional}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>, Maria Dolores Nieto <acinodam at gmail.com> and Luis E. Garrido <garrido.luiseduardo at gmail.com>
#'
#' @return Returns a list containing:
#'
#' \item{network}{A matrix containing a network estimated using 
#' \code{link[EGAnet]{network.estimation}}}
#'
#' \item{wc}{A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network}
#'
#' \item{n.dim}{A scalar of how many total dimensions were identified in the network}
#'
#' \item{correlation}{The zero-order correlation matrix}
#' 
#' \item{n}{Number of cases in \code{data}}
#' 
#' \item{dim.variables}{An ordered matrix of item allocation}
#' 
#' \item{TEFI}{\code{link[EGAnet]{tefi}} for the estimated structure}
#' 
#' \item{plot.EGA}{Plot output if \code{plot.EGA = TRUE}}
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Print results
#' print(ega.wmt)
#' 
#' # Estimate EGAtmfg
#' ega.wmt.tmfg <- EGA(
#'   data = wmt, model = "TMFG",
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Estimate EGA with Louvain algorithm
#' ega.wmt.louvain <- EGA(
#'   data = wmt, algorithm = "louvain",
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#' 
#' # Estimate EGA with an {igraph} function (Fast-greedy)
#' ega.wmt.greedy <- EGA(
#'   data = wmt,
#'   algorithm = igraph::cluster_fast_greedy,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' @references
#' \strong{Original simulation and implementation of EGA} \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#' 
#' \strong{Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data} \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#' 
#' \strong{Compared all \emph{igraph} community detection algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data} \cr
#' \strong{Also implements the Leading Eigenvalue unidimensional method} \cr
#' Christensen, A. P., Garrido, L. E., Pena, K. G., & Golino, H. (2023).
#' Comparing community detection algorithms in psychological data: A Monte Carlo simulation.
#' \emph{Behavior Research Methods}.
#' 
#' \strong{Comprehensive unidimensionality simulation} \cr
#' Christensen, A. P. (2023).
#' Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison.
#' \emph{PsyArXiv}.
#'
#' \strong{Compared all} \code{\link{igraph}} \strong{community detection algorithms, simulation with continuous and polytomous data} \cr
#' Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (2023).
#' Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation.
#' \emph{Behavior Research Methods}.
#' 
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
# EGA ----
# Updated 02.08.2023
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
  corr <- swiftelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  
  # Argument errors
  EGA_errors(data, n, plot.EGA, verbose)
  
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
  
  # Determine `select` arguments
  ## Uses `overwrite_arguments` instead of
  ## `obtain_arguments` because `BGGM::select`
  ## is an S3method. It's possible to use
  ## `BGGM:::select.estimate` but CRAN gets
  ## mad about triple colons
  bggm_select_ARGS <- list( # defaults for `BGGM:::select.estimate`
    cred = 0.95, alternative = "two.sided"
  )
  
  # Obtain arguments for model
  model_ARGS <- switch(
    model_attributes$model,
    "bggm" = c(
      obtain_arguments(BGGM::estimate, model_attributes),
      overwrite_arguments(bggm_select_ARGS, model_attributes)
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
  
  # Add TEFI to the result
  multidimensional_result$TEFI <- tefi(multidimensional_result)$VN.Entropy.Fit
  
  # Check for plot
  if(plot.EGA){
    
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

#' @noRd
# Errors ----
# Updated 27.07.2023
EGA_errors <- function(data, n, plot.EGA, verbose)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame"))
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1)
    typeof_error(n, "numeric")
  }
  
  # 'plot.EGA' errors
  length_error(plot.EGA, 1)
  typeof_error(plot.EGA, "logical")
  
  # 'verbose' errors
  length_error(verbose, 1)
  typeof_error(verbose, "logical")
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 02.08.2023
print.EGA <- function(x, ...)
{
  
  # Make network have S3 class
  class(x$network) <- "EGA.network"
  
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
  if(
    unidimensional_method == "Louvain" &
    "consensus.iter" %in% names(unidimensional_attributes$consensus)
  ){
    
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
      "Unidimensional: ", swiftelse(
        unidimensional_attributes$unidimensional,
        "Yes", "No"
      )
    )
  )
  
  # Check for "TEFI" in output
  if("TEFI" %in% names(x)){
    
    # Add break space
    cat("\n\n----\n\n")
    
    # Print TEFI
    cat(paste0("TEFI: ", round(x$TEFI, 3)))
    
  }
  
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
