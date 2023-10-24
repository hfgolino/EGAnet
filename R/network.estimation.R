#' @title Apply a Network Estimation Method
#'
#' @description General function to apply network estimation methods in \code{\link{EGAnet}}
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
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
#' \item{\code{"cor_auto"} --- }
#' {Uses \code{\link[qgraph]{cor_auto}} to compute correlations. Arguments
#' can be passed along to the function}
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
#' @param network.only Boolean (length = 1).
#' Whether the network only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' network estimation method
#' 
#' @param verbose Boolean (length = 1).
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
# Compute networks for EGA ----
# Updated 24.10.2023
network.estimation <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    network.only = TRUE,
    verbose = FALSE,
    ...
)
{
  
  # Argument errors (return data in case of tibble)
  data <- network.estimation_errors(data, n, network.only, verbose, ...)
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", network.estimation)
  na.data <- set_default(na.data, "pairwise", network.estimation)
  model <- set_default(model, "glasso", network.estimation)
  
  # Check for {BGGM}
  if(model == "bggm"){
    stop("Due to CRAN check issues, `model = \"BGGM\"` is not available at the moment.")
  }
  
  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Make sure there are variable names
  data <- ensure_dimension_names(data)
  
  # First, compute correlation if necessary
  if(model != "bggm"){
    
    # Generic function to get necessary inputs
    output <- obtain_sample_correlations(
      data = data, n = n, corr = corr, 
      na.data = na.data, verbose = verbose, 
      needs_usable = FALSE, # skips usable data check
      ...
    )
    
    # Get correlations and sample size
    correlation_matrix <- output$correlation_matrix; n <- output$n
    
  }else if(is_symmetric(data)){
    
    # 'BGGM' network estimation can't be correlation matrix
    stop(
      "A symmetric matrix was provided in the 'data' argument. 'model = \"BGGM\"' is not compatiable with a correlation matrix. Please use the original data instead",
      call. = FALSE
    )
    
  }
  
  # For 'model = "BGGM"', take a different path...
  if(model == "bggm"){
  
    # # Obtain arguments for `BGGM`
    # bggm_estimate_ARGS <- obtain_arguments(
    #   silent_load(BGGM::estimate), 
    #   FUN.args = ellipse 
    # )
    # 
    # # Check for override of 'type'
    # bggm_estimate_ARGS$type <- find_BGGM_type(data, ellipse)
    # 
    # # Send 'data' and 'verbose' to arguments
    # bggm_estimate_ARGS$Y <- data
    # bggm_estimate_ARGS$progress <- verbose
    # 
    # # Perform BGGM estimate
    # bggm_output <- do.call(
    #   what = BGGM::estimate,
    #   args = bggm_estimate_ARGS
    # )
    # 
    # # Determine `select` arguments
    # ## Uses `overwrite_arguments` instead of
    # ## `obtain_arguments` because `BGGM::select`
    # ## is an S3method. It's possible to use
    # ## `BGGM:::select.estimate` but CRAN gets
    # ## mad about triple colons
    # bggm_select_ARGS <- overwrite_arguments(
    #   list( # defaults for `BGGM:::select.estimate`
    #     cred = 0.95, alternative = "two.sided"
    #   ), ARGS = ellipse
    # )
    # 
    # # Send 'bggm_output' to arguments
    # bggm_select_ARGS$object <- bggm_output
    # 
    # # Perform BGGM select
    # bggm_select <- do.call(
    #   what = BGGM::select,
    #   args = bggm_select_ARGS
    # )
    # 
    # # Obtain network (regardless of 'network.only')
    # estimated_network <- bggm_select$pcor_adj
    
  }else{ # Non-BGGM methods
    
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
    estimation_ARGS$needs_usable <- FALSE # skips usable data check

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
  
  # Transfer variable names
  estimated_network <- transfer_names(data, estimated_network)
  
  # Add methods attribute for BGGM
  ## Methods for GLASSO and TMFG are already there
  # if(model == "bggm"){
  #   attr(estimated_network, "methods") <- c(
  #     bggm_estimate_ARGS[c("type", "analytic", "prior_sd", "iter")],
  #     bggm_select_ARGS[c("cred", "alternative")]
  #   )
  # }
  
  # Add "model", "corr", and "na.data" to attributes
  attr(estimated_network, "methods")[
    c("model", "corr", "na.data")
  ] <- c(model, corr, na.data)
  
  # Set up for return
  if(network.only){
    return(estimated_network) # only return network
  }else{ # sort out output to send back
    
    # BGGM or other model
    if(model == "bggm"){
      
      # # Set up results
      # return(
      #   list(
      #     estimated_network = estimated_network,
      #     output = list(
      #       bggm_estimate = bggm_output,
      #       bggm_select = bggm_select[names(bggm_select) != "object"]
      #     )
      #   )
      # )
      
    }else{
      
      # Set up results
      return(
        list(
          estimated_network = estimated_network,
          output = estimation_OUTPUT
        )
      )
      
    }
  
  }

}

# Bug Checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL;
# corr = "auto"; model = "glasso";
# na.data = "pairwise"; network.only = TRUE;
# verbose = FALSE; ellipse = list();

#' @noRd
# Argument errors ----
# Updated 07.09.2023
network.estimation_errors <- function(data, n, network.only, verbose, ...)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.estimation")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network.estimation")
    typeof_error(n, "numeric", "network.estimation")
  }
  
  # 'network.only' errors
  length_error(network.only, 1, "network.estimation")
  typeof_error(network.only, "logical", "network.estimation")
  
  # 'verbose' errors
  length_error(verbose, 1, "network.estimation")
  typeof_error(verbose, "logical", "network.estimation")
  
  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }
  
  # Return usable data in case of tibble
  return(data)
  
}

#' @noRd
# Send Network Methods for S3 ----
# Updated 07.07.2023
send_network_methods <- function(estimated_network, boot = FALSE)
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
    
    # Send whether analytic solution
    cat(
      paste0(
        "\nAnalytic: ", 
        swiftelse(methods$analytic, "Yes", "No")
      )
    )
    
  }else if(model == "glasso"){ # GLASSO
    
    # Obtain model selection
    model.selection <- methods$model.selection
    
    # Determine whether EBIC was used
    if(isFALSE(boot)){
      
      # Add gamma
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
      
    }else{
      
      # Don't print gamma
      if(model.selection == "ebic"){
        model.selection_text <- paste0(
          " (", toupper(model.selection), ")"
        )
      }else if(model.selection == "jsd"){
        model.selection_text <- paste0(
          " (", toupper(model.selection), ")"
        )
      }
      
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
    
    # Check for bootEGA
    if(isFALSE(boot)){
      
      # Send lambda information
      cat(
        paste0(
          "\nLambda: ", methods$lambda,
          " (n = ", methods$nlambda, ", ratio = ",
          methods$lambda.min.ratio, ")"
        )
      )
      
    }
    
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
        " (", swiftelse(methods$partial, "partial", "zero-order"),
        ")"
      )
    )
    
  }
  
  # Check for bootEGA
  if(isFALSE(boot)){
    cat("\n\n") # Add breakspace
  }

}

#' @exportS3Method 
# S3 Print Method ----
# Updated 06.08.2023
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
  
  # Obtain summary statistics (check for no edges)
  if(edges != 0){
    average_weight <- mean(non_zero_edges, na.rm = TRUE)
    sd_weight <- sd(non_zero_edges, na.rm = TRUE)
    range_weight <- range(non_zero_edges, na.rm = TRUE)
  }else{
    average_weight <- sd_weight <- 0
    range_weight <- numeric(2)
  }
  
  
  # Print information about edges
  cat(paste0("Number of nodes: ", dim(network)[2], "\n"))
  cat(paste0("Number of edges: ", edges, "\n"))
  cat(paste0("Edge density: ", format_decimal(edge_density, 3)))
  
  # Add breakspace
  cat("\n\n")
  
  # Print information about weights
  cat("Non-zero edge weights: \n")
  print_df <- fast.data.frame(
    c(
      format_decimal(average_weight, 3),
      format_decimal(sd_weight, 3),
      format_decimal(range_weight[1], 3),
      format_decimal(range_weight[2], 3)
    ), ncol = 4,
    colnames = c("M", "SD", "Min", "Max")
  )
  print(print_df, quote = FALSE, row.names = FALSE)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 05.07.2023
summary.EGA.network <- function(object, ...)
{
  print(object, ...) # same as print
}

# Function to find 'type' argument for `BGGM` ----
#' @noRd
# Updated 01.07.2023
find_BGGM_type <- function(data, ellipse)
{
  
  # Check if "type" already exists
  if("type" %in% names(ellipse)){
    
    # To avoid conflict "type" with `bootEGA`
    if(
      tolower(ellipse$type) %in% 
      c("binary", "ordinal", "continuous", "mixed")
    ){
      return(ellipse$type)
    }
    
  }
  
  # Obtain categories of the data
  categories <- data_categories(data)
  
  # Obtain unique categories
  unique_categories <- unique_length(categories)
  
  # If all the categories are the same,
  # then either binary or ordinal
  if(isTRUE(unique_categories == 1)){
    
    # Perform switch
    # (since all categories are the same
    # just use the first index)
    type <- unname(swiftelse(categories[1], "binary", "ordinal"))
    
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
      
      # Get range
      category_range <- range(categories)
      
      # Send error
      stop(
        paste0(
          "The `BGGM::estimate` function cannot handle mixed continuous and ordinal data (see documentation `?BGGM::estimate`). ",
          "The range of categories detected in 'data' was between ",
          paste0(category_range, collapse = " and "), ". \n\n",
          "If you believe this error is a mistake, add the argument 'ordinal.categories' and ",
          "set the argument to either `", category_range[1] - 1, "` to treat all variables as ",
          "continuous or `", category_range[2], "` to treat all variables as ordinal"
        ), call. = FALSE
      )
    }
    
  }
  
  # Return type
  return(type)
  
}

