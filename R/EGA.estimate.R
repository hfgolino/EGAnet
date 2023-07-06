#' Estimates \code{EGA} for Multidimensional Structures
#'
#' A basic function to estimate \code{EGA} for multidimensional structures.
#' This function does not include the unidimensional check and it does not
#' plot the results. This function can be used as a streamlined approach
#' for quick \code{EGA} estimation when unidimensionality or visualization
#' is not a priority
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
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @return Returns a list containing:
#'
#' \item{estimated.network}{A symmetric network estimated using either the
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
#' # Estimate EGA
#' ega.wmt <- EGA.estimate(data = wmt)
#' 
#' # Estimate EGA with TMFG
#' ega.wmt.tmfg <- EGA.estimate(data = wmt, model = "TMFG")
#' 
#' # Estimate Bayesian EGA (BEGA)
#' bega.wmt <- EGA.estimate(
#'   data = wmt, model = "BGGM",
#'   analytic = TRUE # faster example for CRAN
#' )
#' 
#' # Estimate EGA with Spinglass algorithm
#' ega.wmt.spinglass <- EGA.estimate(
#'   data = wmt,
#'   algorithm = igraph::cluster_spinglass # any {igraph} algorithm
#' )
#'
#' @references
#' \strong{Original simulation and implementation of EGA} \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#' 
#' \strong{Introduced unidimensional checks, simulation with continuous and dichotomous data} \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#' 
#' 
#' \strong{Compared all} \code{\link{igraph}} \strong{community detection algorithms, simulation with continuous and polytomous data} \cr
#' Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (2023).
#' Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation.
#' \emph{Behavior Research Methods}.
#'
#' @export
#'
# Estimates multidimensional EGA only (no automatic plots)
# Updated 02.07.2023
EGA.estimate <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    verbose = FALSE, ...
)
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- ifelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)

  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Handle legacy arguments (`model.args` and `algorithm.args`)
  ellipse <- legacy_EGA_args(ellipse)
  
  # Ensure data has names
  data <- ensure_dimension_names(data)
  
  # First, get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n, 
    corr = corr, na.data = na.data, 
    verbose = verbose, ...
  )
  
  # Get outputs
  data <- output$data; n <- output$n
  correlation_matrix <- output$correlation_matrix
  
  # Model branching checks (in order of most likely to be used)
  # GLASSO = can handle correlation matrix but needs iterative gamma procedure
  # BGGM = needs original data and not correlation matrix
  # TMFG = can handle correlation matrix
  if(model == "glasso"){
    
    # Use wrapper to clean up iterative gamma procedure
    network <- glasso_wrapper(
      data = correlation_matrix, n = n, corr = corr,
      na.data = na.data, model = model, network.only = TRUE,
      verbose = verbose, ellipse = ellipse
    )
    
  }else if(model == "bggm"){
    
    # Check for correlation input
    if(is_symmetric(data)){
      stop("A symmetric matrix was provided in the 'data' argument. For 'model = \"BGGM\"', the original data is required.")
    }
    
    # Estimate network
    network <- do.call(
      what = network.estimation,
      args = c(
        list( # functions passed into this function
          data = data, n = n, corr = corr, na.data = na.data,
          model = model, network.only = TRUE, verbose = verbose
        ),
        ellipse # pass on ellipse
      )
    )
    
  }else if(model == "tmfg"){
    
    # Estimate network (BGGM *or* TMFG)
    network <- do.call(
      what = network.estimation,
      args = c(
        list( # functions passed into this function
          data = correlation_matrix, n = n, corr = corr, na.data = na.data,
          model = model, network.only = TRUE, verbose = verbose
        ),
        ellipse # pass on ellipse
      )
    )

  }
  
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
          signed =  algorithm == "signed_louvain",
          membership.only = TRUE
        ),
        ellipse # pass on ellipse
      )
    )

  }
  
  # Set up results
  results <- list(
    network = network, wc = wc,
    n.dim = unique_length(wc),
    cor.data = correlation_matrix,
    n = n
  )
  
  # Set class (attributes are stored in `network` and `wc`)
  class(results) <- "EGA.estimate"
  
  # Return results
  return(results)
  
}

# Bug checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL
# corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = igraph::cluster_spinglass
# verbose = FALSE; ellipse = list()

#' @exportS3Method 
# S3 Print Method ----
# Updated 22.06.2023
print.EGA.estimate <- function(x, ...)
{
  
  # Print network estimation
  print(x$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Print community detection
  print(x$wc)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 22.06.2023
summary.EGA.estimate <- function(object, ...)
{
  
  # Print network estimation
  summary(object$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Print community detection
  summary(object$wc)
  
  
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 22.06.2023
plot.EGA.estimate <- function(x, ...)
{
  
  # Return plot
  single_plot(
    network = x$network,
    wc = x$wc,
    ...
  )
  
}

#' @noRd
# Wrapper for GLASSO ----
# Updated 02.07.2023
glasso_wrapper <- function(
    network, data, n, corr, na.data,
    model, network.only, verbose,
    ellipse
)
{

  # Remove `gamma` from ellipse
  if("gamma" %in% names(ellipse)){
    ellipse <- ellipse[-which(names(ellipse) == "gamma")]
  }
  
  # Set gamma 
  gamma <- 0.50
  
  # Check for any disconnected nodes
  while(TRUE){
    
    # Re-estimate network
    network <- do.call(
      what = network.estimation,
      args = c(
        list( # functions passed into this function
          data = data, n = n, corr = corr, na.data = na.data,
          model = model, network.only = TRUE, verbose = verbose,
          gamma = gamma
        ),
        ellipse # pass on ellipse
      )
    )
    
    # Check for disconnected nodes
    if(any(strength(network) == 0) & gamma != 0){
      gamma <- gamma - 0.25 # decrease gamma
    }else{
      break # all nodes are connected or gamma equals zero
    }
    
  }
  
  # Return network
  return(network)
  
}

