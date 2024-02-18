#' @title Approaches to Detect Unidimensional Communities
#'
#' @description A function to apply several approaches to detect a unidimensional community in
#' networks. There have many different approaches recently such as expanding
#' the correlation matrix to have orthogonal correlations (\code{"expand"}),
#' applying the Leading Eigenvalue community detection algorithm
#' \code{\link[igraph]{cluster_leading_eigen}} to the correlation matrix
#' (\code{"LE"}), and applying the Louvain community detection algorithm
#' \code{\link[igraph]{cluster_louvain}} to the correlation matrix (\code{"louvain"}).
#' Not necessarily intended for individual use -- it's better to use \code{\link[EGAnet]{EGA}}
#'
#' @param data Matrix or data frame.
#' Should consist only of variables that are desired to be in analysis
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
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)
#'
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
#' Arguments can be passed along to the function
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
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
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"BGGM"} --- Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{?BGGM::estimate} for more details
#'
#' \item \code{"glasso"} --- Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details
#'
#' \item \code{"nonreg"} --- Computes the Maximum Likelihood non-regularized
#' approach. See \code{\link[EGAnet]{ggm_inference.GGMnonreg}} for more details
#'
#' \item \code{"TMFG"} --- Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details
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
#' \item \code{"expand"} --- Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation
#'
#' \item \code{"LE"} --- Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023)
#' \emph{Behavior Research Methods} simulation
#'
#' \item \code{"louvain"} --- Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Louvain solution is used; otherwise, regular EGA is used.
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}
#'
#' }
#'
#' @param verbose Boolean.
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{community.detection}}
#'
#' @return Returns the memberships of the community detection algorithm.
#' The memberships will output \emph{regardless} of whether the
#' network is unidimensional
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Louvain (default)
#' community.unidimensional(wmt)
#'
#' # Louvain with consensus clustering
#' community.unidimensional(wmt, consensus.iter = 1000)
#'
#' # Leading Eigenvector
#' community.unidimensional(wmt, uni.method = "LE")
#'
#' # Expand
#' community.unidimensional(wmt, uni.method = "expand")
#'
#' @references
#' \strong{Expand approach} \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020).
#' Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors:
#' A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#'
#' \strong{Leading Eigenvector approach} \cr
#' Christensen, A. P., Garrido, L. E., Guerra-Pena, K., & Golino, H. (2023).
#' Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation.
#' \emph{Behavior Research Methods}.
#'
#' \strong{Louvain approach} \cr
#' Christensen, A. P. (2023).
#' Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison.
#' \emph{PsyArXiv}.
#'
#' @export
#'
# Compute unidimensional approaches for EGA
# Updated 18.02.2024
community.unidimensional <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "nonreg", "TMFG"),
    uni.method = c("expand", "LE", "louvain"),
    verbose = FALSE,
    ...
)
{

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", community.unidimensional)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)

  # Argument errors (return data in case of tibble)
  data <- community.unidimensional_errors(data, n, verbose, ...)

  # Make sure there are variable names
  data <- ensure_dimension_names(data)

  # Generic function to get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n, corr = corr,
    na.data = na.data, verbose = verbose,
    needs_usable = FALSE, # skips usable data check
    ...
  )

  # Check for incompatible method combinations
  if(model %in% c("bggm", "nonreg") && uni.method == "expand"){

    # Return unidimensional approach
    # No S3 methods -- not intended for individual use
    return(expand_data(output$data, output$n, model, list(...)))

  }else{

    # Return unidimensional approach
    # No S3 methods -- not intended for individual use
    return(
      switch( # Ordered by most common usage
        uni.method,
        "louvain" = consensus_wrapper(output$correlation_matrix, verbose, list(...)),
        "le" = community.detection(output$correlation_matrix, algorithm = "leading_eigen", ...),
        "expand" = expand(output$correlation_matrix, output$n, model, verbose, list(...))
      )
    )

  }

}

# Bug Checking ----
# ## Basic input
# data = wmt2[,7:24]; n = NULL;
# corr = "auto"; na.data = "pairwise";
# model = "glasso"; uni.method = "expand";
# verbose = FALSE; ellipse = list()

#' @noRd
# Errors ----
# Updated 07.09.2023
community.unidimensional_errors <- function(data, n, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "community.unidimensional")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "community.unidimensional")
    typeof_error(n, "numeric", "community.unidimensional")
  }

  # 'verbose' errors
  length_error(verbose, 1, "community.unidimensional")
  typeof_error(verbose, "logical", "community.unidimensional")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

#' @noRd
# "Expand" Correlation approach ----
# Updated 12.07.2023
expand <- function(correlation_matrix, n, model, verbose, ellipse)
{

  # Number of variables
  variables <- dim(correlation_matrix)[2]

  # New total variables
  new_total <- variables + 4

  # Create new matrix
  expanded_matrix <- matrix(0, nrow = new_total, ncol = new_total)

  # Create orthogonal correlation matrix
  orthogonal_matrix <- matrix(0.50, nrow = 4, ncol = 4)
  diag(orthogonal_matrix) <- 1

  # Set original dimensions
  original_dimensions <- seq_len(variables)

  # Insert original correlation matrix
  expanded_matrix[
    original_dimensions, original_dimensions
  ] <- correlation_matrix

  # Set expanded dimensions
  expanded_dimensions <- (variables + 1):new_total

  # Insert orthogonal correlations
  expanded_matrix[
    expanded_dimensions, expanded_dimensions
  ] <- orthogonal_matrix

  # Obtain estimation method function
  estimation_FUN <- switch(
    model,
    "glasso" = EBICglasso.qgraph,
    "tmfg" = TMFG
  )

  # Obtain estimation method arguments
  network_ARGS <- obtain_arguments(estimation_FUN, ellipse)

  # Set data
  if(model == "tmfg"){ # Use correlation matrix
    network_ARGS$data <- correlation_matrix
  }else{ # Normal approach
    network_ARGS$data <- expanded_matrix
  }

  # Set data, sample size, output, and verbose
  network_ARGS$n <- n
  network_ARGS$returnAllResults <- FALSE

  # Apply network estimation method
  network <- do.call(estimation_FUN, network_ARGS)

  # Obtain community detection arguments
  community_ARGS <- obtain_arguments(community.detection, ellipse)

  # Set network
  if(model == "tmfg"){

    # Set network into expanded matrix
    expanded_matrix[original_dimensions, original_dimensions] <- network

    # Use expanded matrix as network
    community_ARGS$network <- expanded_matrix

  }else{ # Normal approach
    community_ARGS$network <- network
  }

  # Apply community detection algorithm
  membership <- do.call(community.detection, community_ARGS)

  # Remove additional variables
  membership <- membership[original_dimensions]

  # Add back names
  names(membership) <- dimnames(correlation_matrix)[[2]]

  # Return membership
  return(membership)

}

#' @noRd
# "Expand" Data approach ----
# Updated 18.02.2024
expand_data <- function(data, n, model, ellipse)
{

  # Set Cholesky based on a population correlation matrix
  # of all r's = 0.50 (i.e., loadings = 0.70)
  cholesky <- matrix(
    c(
      1, 0.5000000, 0.5000000, 0.5000000,
      0, 0.8660254, 0.2886751, 0.2886751,
      0, 0.0000000, 0.8164966, 0.2041241,
      0, 0.0000000, 0.0000000, 0.7905694
    ), nrow = 4, ncol = 4, byrow = TRUE
  )

  # Generate data
  simulated_data <- MASS_mvrnorm_quick(
    seed = NULL, p = 4, np = 4 * n, diag(4)
  ) %*% cholesky

  # Get median categories of original data
  original_categories <- median(data_categories(data), na.rm = TRUE)

  # Check for need to categories
  if(original_categories <= 6){

    # Categorize the data
    simulated_data <- expand_categorize(
      simulated_data, original_categories
    )

  }

  # Add variable names
  dimnames(simulated_data)[[2]] <- paste0("sim_V", 1:4)

  # Combine data
  combined_data <- cbind(simulated_data, data)

  # Ensure 'verbose' is FALSE
  ellipse$verbose <- FALSE

  # Apply BGGM or non-regularized
  unidimensional_output <- do.call(
    what = EGA.estimate,
    args = c(
      list(
        data = combined_data,
        n = n, model = model
      ),
      ellipse
    )
  )

  # Return memberships
  return(reindex_memberships(unidimensional_output$wc[-c(1:4)]))


}

#' @noRd
# Categorization function adapted from {latentFactoR}
# Updated 13.10.2023
expand_categorize <- function(data, categories)
{

  # Skew is always zero
  skew_values <- switch(
    as.character(categories),
    "2" = 0,
    "3" = c(-0.4307, 0.4307),
    "4" = c(-0.6745, 0.0000, 0.6745),
    "5" = c(-0.8416, -0.2533, 0.2534, 0.8416),
    "6" = c(-0.9674, -0.4307, 0.0000, 0.4307, 0.9674)
  )

  # Categorize biased data with updated thresholds
  for(i in (length(skew_values) + 1):1){

    # First category
    if(i == 1){
      data[data < skew_values[i]] <- i
    }else if(i == length(skew_values) + 1){ # Last category
      data[data >= skew_values[i-1]] <- i
    }else{ # Middle category
      data[data >= skew_values[i-1] & data < skew_values[i]] <- i
    }

  }

  # Return categorized data
  return(data)

}

#' @noRd
# Wrapper for Louvain consensus ----
# Updated 23.07.2023
consensus_wrapper <- function(correlation_matrix, verbose, ellipse)
{

  # Check for consensus method
  if(any(c("consensus.method", "consensus.iter") %in% names(ellipse))){

    # Obtain arguments
    consensus_ARGS <- obtain_arguments(community.consensus, ellipse)

    # Set arguments
    consensus_ARGS[
      c("network", "consensus.method", "membership.only")
    ] <- list(correlation_matrix, "most_common", TRUE)

    # Apply Louvain with consensus approach
    membership <- do.call(
      what = community.consensus,
      args = consensus_ARGS
    )

  }else{ # Use single-shot Louvain
    membership <- community.detection(correlation_matrix, algorithm = "louvain")
  }

  # Add back names
  names(membership) <- dimnames(correlation_matrix)[[2]]

  # Return membership
  return(membership)

}
