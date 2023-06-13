#' Apply a Network Estimation Method
#'
#' General function to apply network estimation methods in \code{\link{EGAnet}}
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
#' @param method Character (length = 1).
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
#' @param network.only Boolean.
#' Whether the network only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' network estimation method
#' 
#' @param verbose Boolean.
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}} and the different
#' network estimation methods (see \code{method} for
#' method specific details)
#'
#' @return Returns a matrix populated with a network from the input data
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # BGGM with analytic solution
#' bggm_network <- network.estimation(
#'   data = wmt, method = "BGGM",
#'   analytic = TRUE # faster example for CRAN
#' )
#' 
#' # EBICglasso
#' glasso_network <- network.estimation(
#'   data = wmt, method = "glasso"
#' )
#' 
#' # TMFG
#' tmfg_network <- network.estimation(
#'   data = wmt, method = "TMFG"
#' )
#'
#' @references
#' # Graphical Least Absolute Shrinkage and Selection Operator (GLASSO)
#' 
#' # GLASSO with Extended Bayesian Information Criterion (EBICglasso)
#' 
#' # Bayesian Gaussian Graphical Model (BGGM)
#'
#' # Triangulated Maximally Filtered Graph (TMFG)
#'
#' @export
#'
# Compute networks for EGA
# Updated 13.06.2023
network.estimation <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    method = c("BGGM", "glasso", "TMFG"),
    network.only = TRUE,
    verbose = FALSE,
    ...
)
{
  
  # Missing arguments
  ## Correlation method
  if(missing(corr)){
    corr <- "auto"
  }else{corr <- tolower(match.arg(corr))}
  ## Missing data method
  if(missing(na.data)){
    na.data <- "pairwise"
  }else{na.data <- tolower(match.arg(na.data))}
  ## Network method
  if(missing(method)){ # No default!
    stop("No network estimation 'method' was selected. Please select either \"BGGM\", \"glasso\", or \"TMFG\"")
  }else{method <- tolower(match.arg(method))}
  # Makes 'method' *lowercase* for parsimonious handling later
  
  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Make sure there are variable names
  data <- ensure_dimension_names(data)
  
  # Check for whether data are data or correlation matrix
  if(is_symmetric(data)){
    
    # Check for sample size
    if(is.null(n)){
      stop("A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.")
    }
    
    # Check for 'BGGM' network estimation (can't be correlation matrix)
    if(method == "bggm"){
      stop("A symmetric matrix was provided in the 'data' argument. 'method = \"BGGM\"' is not compatiable with a correlation matrix. Please use the original data instead")
    }
    
    # Set data as correlation matrix
    correlation_matrix <- data
    
  }else{ # Assume 'data' is data
    
    # Check for appropriate variables
    data <- usable_data(data, verbose)
    
    # Obtain sample size
    n <- nrow(data)
    
  }
  
  # For 'method = "BGGM"', take a different path...
  if(method == "bggm"){
  
    # Obtain arguments for `BGGM`
    ## See "helpers-general.R" for `silent_load`
    bggm_estimate_ARGS <- obtain_arguments(
      silent_load(BGGM::estimate), 
      FUN.args = ellipse
    )
    
    # Check for override of 'type'
    if(!"type" %in% names(ellipse)){
      
      # 'type' was not supplied; apply appropriate type for user
      bggm_estimate_ARGS$type <- find_BGGM_type(data, ellipse)

    }
    
    # Send 'data' and 'verbose' to arguments
    bggm_estimate_ARGS$Y <- data
    bggm_estimate_ARGS$progress <- verbose
    
    # Perform BGGM estimate
    bggm_output <- do.call(
      what = BGGM::estimate,
      args = bggm_estimate_ARGS
    )
    
    # Determine `select` arguments
    bggm_select_ARGS <- obtain_arguments(BGGM:::select.estimate, FUN.args = ellipse)
    
    # Send 'bggm_output' to arguments
    bggm_select_ARGS$object <- bggm_output
    
    # Perform BGGM select
    bggm_select <- do.call(
      what = BGGM::select,
      args = bggm_select_ARGS
    )
    
    # Obtain network (regardless of 'network.only')
    estimated_network <- bggm_select$pcor_adj
    
  }else{ # Non-BGGM methods
    
    # Check for automatic correlations
    if(corr == "auto"){
      
      # Add arguments to 'ellipse'
      ellipse$corr <- "pearson"; ellipse$na.data <- na.data;
      
      # Obtain arguments for `auto.correlate`
      auto_ARGS <- obtain_arguments(FUN = auto.correlate, FUN.args = ellipse)
      
      # Supply data
      auto_ARGS$data <- data
      
      # Obtain correlation matrix
      correlation_matrix <- do.call(
        what = auto.correlate,
        args = auto_ARGS
      )
      
    }else{
      
      # Obtain correlations using base R
      correlation_matrix <- cor(data, use = na.data, method = corr)
      
    }
    
    # Obtain estimation method function
    estimation_FUN <- switch(
      method,
      "glasso" = EBICglasso.qgraph,
      "tmfg" = TMFG
    )
    
    # Obtain estimation method arguments
    estimation_ARGS <- obtain_arguments(estimation_FUN, ellipse)
    
    # Set data, sample size, output, and verbose
    estimation_ARGS$data <- correlation_matrix
    estimation_ARGS$n <- n
    estimation_ARGS$returnAllResults <- !network.only
    estimation_ARGS$verbose <- verbose
    
    # Perform estimation
    estimation_OUTPUT <- do.call(
      what = estimation_FUN,
      args = estimation_ARGS
    )

    # Obtain estimated network
    if(!isTRUE(network.only)){
      
      # Switch with method
      estimated_network <- estimation_OUTPUT$network
      
    }else{
      
      # Output is network only
      estimated_network <- estimation_OUTPUT
      
    }
    
  }
  
  # At the end, ensure network is named
  colnames(estimated_network) <- colnames(data)
  row.names(estimated_network) <- colnames(data)
  
  # Set up for return
  if(isTRUE(network.only)){ # only return network
    return(estimated_network)
  }else{ # sort out output to send back
    
    # BGGM or other method
    if(method == "bggm"){
      
      # Set up results
      results <- list(
        estimated_network = estimated_network,
        bggm_estimate = bggm_output,
        bggm_select = bggm_select
      )
      
    }else{
      
      # Set up results
      results <- list(
        estimated_network = estimated_network,
        method_output = estimation_OUTPUT
      )
      
    }
    
    # Return results
    return(results)
  
  }

}

# Bug Checking ----
## Basic input
# data = wmt2; n = NULL;
# corr = "auto"; method = "bggm";
# na.data = "pairwise"; network.only = TRUE;
# verbose = FALSE; ellipse = list();

# Function to find 'type' argument for `BGGM` ----
#' @noRd
# Updated 09.06.2023
find_BGGM_type <- function(data, ellipse)
{
  
  # Obtain categories of the data
  categories <- data_categories(data)
  
  # Obtain unique categories
  unique_categories <- unique(categories)
  
  # If all the categories are the same,
  # then either binary or ordinal
  if(isTRUE(length(unique_categories) == 1)){
    
    # Perform switch
    type <- ifelse(unique_categories == 2, "binary", "ordinal")
    
  }else{ # Not all are the same, then mixed, continuous, or not handled
    
    # If missing 'ordinal.categories', then make it so! (it's not)
    if(!"ordinal.categories" %in% names(ellipse)){
      ellipse$ordinal.categories <- 7 # same as `auto.correlate`
    }
    
    # Determine whether continuous or mixed is appropriate
    if(all(categories > ellipse$ordinal.categories)){
      type <- "continuous"
    }else if(all(categories >= 2) & all(categories <= ellipse$ordinal.categories)){
      type <- "mixed"
    }else{ # not handled... :(
      stop(
        paste0(
          "The `BGGM::estimate` function cannot handle mixed continuous and ordinal data (see documentation `?BGGM::estimate`). ",
          "The range of categories detected in 'data' was between ",
          paste0(range(categories), collapse = " and "), ". \n\n",
          "If you believe this error is a mistake, add the argument 'ordinal.categories' and ",
          "set the argument to either `", min(categories) - 1, "` to treat all variables as ",
          "continuous or `", max(categories), "` to treat all variables as ordinal"
        )
      )
    }
    
  }
  
  # Return type
  return(type)
  
}

