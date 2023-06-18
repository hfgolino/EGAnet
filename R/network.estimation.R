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
#' network estimation methods (see \code{model} for
#' model specific details)
#'
#' @return Returns a matrix populated with a network from the input data
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # EBICglasso (default for EGA functions)
#' glasso_network <- network.estimation(
#'   data = wmt, model = "glasso"
#' )
#' 
#' # Plot network
#' plot(glasso_network)
#' 
#' # BGGM with analytic solution
#' bggm_network <- network.estimation(
#'   data = wmt, model = "BGGM",
#'   analytic = TRUE # faster example for CRAN
#' )
#' 
#' # TMFG
#' tmfg_network <- network.estimation(
#'   data = wmt, model = "TMFG"
#' )
#'
#' @references
#' \strong{Graphical Least Absolute Shrinkage and Selection Operator (GLASSO)} \cr
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008). 
#' Sparse inverse covariance estimation with the graphical lasso. 
#' \emph{Biostatistics}, \emph{9}(3), 432–441.
#' 
#' \strong{GLASSO with Extended Bayesian Information Criterion (EBICglasso)} \cr
#' Epskamp, S., & Fried, E. I. (2018).
#' A tutorial on regularized partial correlation networks.
#' \emph{Psychological Methods}, \emph{23}(4), 617–634.
#' 
#' \strong{Bayesian Gaussian Graphical Model (BGGM)} \cr
#' Williams, D. R. (2021).
#' Bayesian estimation for Gaussian graphical models: Structure learning, predictability, and network comparisons.
#' \emph{Multivariate Behavioral Research}, \emph{56}(2), 336–352. 
#'
#' \strong{Triangulated Maximally Filtered Graph (TMFG)} \cr
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}, 161-178.
#'
#' @export
#'
# Compute networks for EGA
# Updated 16.06.2023
network.estimation <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    network.only = TRUE,
    verbose = FALSE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", network.estimation)
  na.data <- set_default(na.data, "pairwise", network.estimation)
  model <- set_default(model, "glasso", network.estimation)
  
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
    if(model == "bggm"){
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
  
  # For 'model = "BGGM"', take a different path...
  if(model == "bggm"){
  
    # Obtain arguments for `BGGM`
    ## See "helpers.R" for `silent_load`
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
    bggm_select_ARGS <- obtain_arguments(
      BGGM:::select.estimate, FUN.args = ellipse
    )
    
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
    
    # Check for whether `correlation_matrix` exists
    if(!exists("correlation_matrix")){
      
      # Check for automatic correlations
      if(corr == "auto"){
        
        # Obtain correlation matrix
        correlation_matrix <- auto.correlate(
          data = data, corr = "pearson",
          na.data = na.data, verbose = verbose,
          ...
        )
        
      }else{
        
        # Obtain correlations using base R
        correlation_matrix <- cor(data, use = na.data, method = corr)
        
      }
      
    }
    
    # Obtain estimation method function
    estimation_FUN <- switch(
      model,
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
    if(isFALSE(network.only)){
      
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
  
  # Add class
  class(estimated_network) <- "EGA.network"
  
  # Add methods attribute for BGGM
  ## Methods for GLASSO and TMFG are already there
  if(model == "bggm"){
    attr(estimated_network, "methods") <- list(
      type = bggm_estimate_ARGS$type,
      analytic = bggm_estimate_ARGS$analytic,
      prior_sd = bggm_estimate_ARGS$prior_sd,
      iter = bggm_estimate_ARGS$iter,
      cred = bggm_select_ARGS$cred,
      alternative = bggm_select_ARGS$alternative
    )
  }
  
  # Add model to attributes
  attr(estimated_network, "methods")$model <- model
  
  # Set up for return
  if(isTRUE(network.only)){
    return(estimated_network) # only return network
  }else{ # sort out output to send back
    
    # BGGM or other model
    if(model == "bggm"){
      
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
        output = estimation_OUTPUT
      )
      
    }
    
    # Return results
    return(results)
  
  }

}

# Bug Checking ----
## Basic input
# data = wmt2; n = NULL;
# corr = "auto"; model = "glasso";
# na.data = "pairwise"; network.only = TRUE;
# verbose = FALSE; ellipse = list();

#' @noRd
# Send Network Methods for S3 ----
# Updated 16.06.2023
send_network_methods <- function(estimated_network)
{
  
  # Obtain methods
  methods <- attr(estimated_network, "methods")
  
  # Set model
  model <- methods$model
  
  # Send output text based on model
  if(model == "bggm"){ # BGGM
    
    # Send model
    cat(
      paste0(
        "Model: ", toupper(methods$model),
        " (", methods$type, ")"
      )
    )
    
    # Send prior information
    cat(
      paste0(
        "\nPrior SD: ", methods$prior_sd,
        " (", methods$iter, " iterations)"
      )
    )
  
    # Send statistical information
    cat(
      paste0(
        "\nCredible Interval: ", methods$cred,
        " (", gsub("\\.", "-", methods$alternative), ")"
      )
    )
    
  }else if(model == "glasso"){ # GLASSO
    
    # Obtain model selection
    model.selection <- methods$model.selection
    
    # Determine whether EBIC was used
    if(model.selection == "ebic"){
      model.selection_text <- paste0(
        " (", toupper(model.selection), 
        " with gamma = ", methods$gamma, ")"
      )
    }else if(model.selection == "jsd"){
      model.selection_text <- paste0(
        " (", toupper(model.selection), ")"
      )
    }
    
    # Send model
    cat(
      paste0(
        "Model: ", toupper(methods$model),
        model.selection_text
      )
    )
    
    # Send correlations
    cat(paste0("\nCorrelations: ", methods$corr))
    
    # Send lambda information
    cat(
      paste0(
        "\nLambda: ", methods$lambda,
        " (n = ", methods$nlambda, ", ratio = ",
        methods$lambda.min.ratio, ")"
      )
    )
    
  }else if(model == "tmfg"){ # TMFG
    
    # Send model
    cat(
      paste0(
        "Model: ", toupper(methods$model)
      )
    )
    
    # Send correlations
    cat(
      paste0(
        "\nCorrelations: ", methods$corr,
        " (", ifelse(methods$partial, "partial", "zero-order"),
        ")"
      )
    )
    
  }
  
  # Add breakspace
  cat("\n\n")

}

#' @exportS3Method 
# S3 Print Method ----
# Updated 14.06.2023
print.EGA.network <- function(x, ...)
{
  
  # Determine whether result is a list
  if(is.list(x)){
    network <- x$estimated_network
  }else{
    network <- x
  }
  
  # Print methods information
  send_network_methods(network)
  
  # Obtain lower triangle
  lower_triangle <- network[lower.tri(network)]
  
  # Non-zero edges
  non_zero_edges <- lower_triangle[lower_triangle != 0]
  
  # Determine number and density of edges
  edges <- length(non_zero_edges)
  edge_density <- edges / length(lower_triangle)
  
  # Obtain summary statistics
  average_weight <- mean(non_zero_edges, na.rm = TRUE)
  sd_weight <- sd(non_zero_edges, na.rm = TRUE)
  range_weight <- range(non_zero_edges, na.rm = TRUE)
  
  # Print information about edges
  cat(paste0("Number of nodes: ", dim(network)[2], "\n"))
  cat(paste0("Number of edges: ", edges, "\n"))
  cat(paste0("Edge density: ", format_decimal(edge_density, 3)))
  
  # Add breakspace
  cat("\n\n")
  
  # Print information about weights
  cat("Non-zero edge weights: ")
  print_df <- data.frame(
    c("M", format_decimal(average_weight, 3)),
    c("SD", format_decimal(sd_weight, 3)),
    c("Min", format_decimal(range_weight[1], 3)),
    c("Max", format_decimal(range_weight[2], 3))
  )
  no_name_print(print_df)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 14.06.2023
summary.EGA.network <- function(object, ...)
{
  
  # Determine whether result is a list
  if(is.list(object)){
    network <- object$estimated_network
  }else{
    network <- object
  }
  
  # Print methods information
  send_network_methods(network)
  
  # Obtain lower triangle
  lower_triangle <- network[lower.tri(network)]
  
  # Non-zero edges
  non_zero_edges <- lower_triangle[lower_triangle != 0]
  
  # Determine number and density of edges
  edges <- length(non_zero_edges)
  edge_density <- edges / length(lower_triangle)
  
  # Obtain summary statistics
  average_weight <- mean(non_zero_edges, na.rm = TRUE)
  sd_weight <- sd(non_zero_edges, na.rm = TRUE)
  range_weight <- range(non_zero_edges, na.rm = TRUE)
  
  # Print information about edges
  cat(paste0("Number of nodes: ", dim(network)[2], "\n"))
  cat(paste0("Number of edges: ", edges, "\n"))
  cat(paste0("Edge density: ", format_decimal(edge_density, 3)))
  
  # Add breakspace
  cat("\n\n")
  
  # Print information about weights
  cat("Non-zero edge weights: ")
  print_df <- data.frame(
    c("M", format_decimal(average_weight, 3)),
    c("SD", format_decimal(sd_weight, 3)),
    c("Min", format_decimal(range_weight[1], 3)),
    c("Max", format_decimal(range_weight[2], 3))
  )
  no_name_print(print_df)
  
  
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 15.06.2023
plot.EGA.network <- function(x, ...)
{
  
  # Determine whether result is a list
  if(is.list(x)){
    network <- x$estimated_network
  }else{
    network <- x
  }
  
  # Return plot
  basic_plot_setup(network = network, ...)
  
}

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

