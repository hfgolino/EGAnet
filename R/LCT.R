#' @title Loadings Comparison Test
#'
#' @description An algorithm to identify whether data were generated from a
#' factor or network model using factor and network loadings.
#' The algorithm uses heuristics based on theory and simulation. These
#' heuristics were then submitted to several deep learning neural networks
#' with 240,000 samples per model with varying parameters.
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
#' \item \code{"TMFG"} --- Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details
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
#' \item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details
#' 
#' \item \code{"louvain"} --- By default, \code{"louvain"} will implement the Louvain algorithm using 
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}} 
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000} 
#' unless specified otherwise
#' 
#' \item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details
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
#' @param iter Numeric (length = 1).
#' Number of replicate samples to be drawn from a multivariate
#' normal distribution (uses \code{MASS::mvrnorm}).
#' Defaults to \code{100} (recommended)
#' 
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results.
#' Set for reproducible results.
#' See \href{https://github.com/hfgolino/EGAnet/wiki/Reproducibility-and-PRNG}{Reproducibility and PRNG}
#' for more details on random number generation in \code{\link{EGAnet}}
#' 
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments that can be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{EGA}}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen at gmail.com>
#'
#' @return Returns a list containing:
#' 
#' \item{empirical}{Prediction of model based on empirical dataset only}
#' 
#' \item{bootstrap}{Prediction of model based on means of the loadings across
#' the bootstrap replicate samples}
#' 
#' \item{proportion}{Proportions of models suggested across bootstraps}
#'
#' @examples
#' # Get data
#' data <- psych::bfi[,1:25]
#' 
#' \dontrun{# Compute LCT
#' ## Factor model
#' LCT(data)} 
#' 
#' @references
#' \strong{Model training and validation} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Factor or network model? Predictions from neural networks.
#' \emph{Journal of Behavioral Data Science}, \emph{1}(1), 85-126.
#' 
#' @export
#'
# Loadings Comparison Test ----
# Updated 05.09.2023
LCT <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    iter = 100, seed = NULL, verbose = TRUE, ...
)
{
  
  # Store random state (if there is one)
  store_state()
  
  # Argument errors (return data in case of tibble)
  data <- LCT_errors(data, n, iter, verbose)
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", LCT)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  
  # Ensure data has names
  data <- ensure_dimension_names(data)
  
  # First, get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n, 
    corr = corr, na.data = na.data, 
    verbose = verbose, ...
  )
  
  # Get parameters for parametric bootstrap
  mvrnorm_parameters <- mvrnorm_precompute(
    cases = output$n,
    Sigma = output$correlation
  )
  
  # Initialize factor and network loading matrices
  fl_proportions <- nl_proportions <- matrix(0, nrow = iter, ncol = 5)
  
  # Set up progress bar
  if(verbose){pb <- txtProgressBar(max = iter, style = 3)}
  
  # Check for seed
  if(!is.null(seed)){
    seeds <- reproducible_seeds(iter, seed)
  }else{ 
    
    # Set all seeds to zero (or random)
    seeds <- rep(0, iter)
    
    # Send message about NULL seed
    message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results")
    
  }
  
  # Perform loop
  for(iteration in seq_len(iter)){
    
    # Generate data
    generated_data <- reproducible_bootstrap(
      seed = seeds[iteration], type = "parametric",
      mvrnorm_parameters = mvrnorm_parameters
    )
    
    # Estimate EGA
    ega <- EGA(generated_data, plot.EGA = FALSE, ...)
    
    # Get correlation matrix
    generated_correlation <- ega$correlation
    
    # Remove variables with missing membership
    keep_membership <- !is.na(ega$wc)
    ega$wc <- ega$wc[keep_membership]
    ega$network <- ega$network[keep_membership, keep_membership]
    generated_correlation <- generated_correlation[keep_membership, keep_membership]
    
    # Estimate network loadings
    network_loadings <- abs(net.loads(ega, loading.method = "BRM")$std)
    
    # Get proportions
    nl_proportions[iteration,] <- c(
      mean(network_loadings >= 0.15, na.rm = TRUE),
      mean(network_loadings >= 0.25, na.rm = TRUE),
      mean(network_loadings >= 0.35, na.rm = TRUE),
      assigned_proportions(network_loadings, ega$wc, "network"),
      unassigned_proportions(network_loadings, ega$wc, "network")
    )
    
    # Get dimension sequence
    dimension_sequence <- seq_len(ega$n.dim)
    
    # Estimate factor loadings
    factor_loadings <- silent_call(
      abs(
        psych::fa(
          generated_correlation, nfactors = ega$n.dim, n.obs = output$n
        )$loadings[,dimension_sequence]
      )
    )
    
    # Get "memberships"
    factor_wc <- max.col(factor_loadings, ties.method = "first")
    dimnames(factor_loadings)[[2]] <- dimension_sequence
    names(factor_wc) <- dimnames(factor_loadings)[[1]]
    
    # Get proportions
    fl_proportions[iteration,] <- c(
      mean(factor_loadings >= 0.40, na.rm = TRUE),
      mean(factor_loadings >= 0.55, na.rm = TRUE),
      mean(factor_loadings >= 0.70, na.rm = TRUE),
      assigned_proportions(factor_loadings, factor_wc, "factor"),
      unassigned_proportions(factor_loadings, factor_wc, "factor")
    )
    
    # Update progress bar
    if(verbose){setTxtProgressBar(pb, iteration)}
    
  }
  
  # Close progress bar
  if(verbose){close(pb)}
  
  # Combine proportions into a matrix
  loadings_matrix <- cbind(nl_proportions, fl_proportions)
  
  # Set NA values to zero
  loadings_matrix[] <- swiftelse(is.na(loadings_matrix), 0, loadings_matrix)

  # Collapse across matrix
  bootstrap_predictions <- apply(na.omit(loadings_matrix), 1, dnn.predict)
  
  # Get table
  bootstrap_predictions <- fast_table(bootstrap_predictions) / iter
  
  # Set new table
  new_table <- c("1" = 0, "2" = 0)
  
  # Set names
  new_table[names(bootstrap_predictions)] <- bootstrap_predictions
  names(new_table) <- c("Factor", "Network")
  
  # Restore random state (if there is one)
  restore_state()

  # Return results
  return(
    list(
      # Single-shot
      empirical = swiftelse(
        dnn.predict(loadings_matrix[1,]) == 1,
        "Factor", "Network"
      ),
      # Bootstrap
      bootstrap = swiftelse(
        dnn.predict(colMeans(loadings_matrix, na.rm = TRUE)) == 1,
        "Factor", "Network"
      ),
      # Proportion
      proportion = new_table
    )
  )
  
}

#' @noRd
# Errors ----
# Updated 19.08.2023
LCT_errors <- function(data, n, iter, verbose)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "LCT")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "LCT")
    typeof_error(n, "numeric", "LCT")
  }
  
  # 'iter' errors
  length_error(iter, 1, "LCT")
  typeof_error(iter, "numeric", "LCT")
  range_error(iter, c(1, Inf), "LCT")
  
  # 'verbose' errors
  length_error(verbose, 1, "LCT")
  typeof_error(verbose, "logical", "LCT")
  
  # Return usable data in case of tibble
  return(usable_data(data, verbose))
  
}

#' @noRd
# Assigned proportions ----
# Updated 05.08.2023
assigned_proportions <- function(loadings, memberships, type)
{
  
  # Get cut-off value
  cut_off <- swiftelse(type == "network", 0.15, 0.40)
  
  # Return assigned proportions
  return(
    mean(
      ulapply(unique(memberships), function(community){
        
        # Get community index
        community_index <- memberships == community
        
        # Return assigned proportion
        return(loadings[names(memberships)[community_index], as.character(community)] >= cut_off)
        
      }),
      na.rm = TRUE
    )
  )
  
}

#' @noRd
# Unassigned proportions ----
# Updated 05.08.2023
unassigned_proportions <- function(loadings, memberships, type)
{
  
  # Get cut-off value
  cut_off <- swiftelse(type == "network", 0.15, 0.40)
  
  # Return assigned proportions
  return(
    mean(
      ulapply(unique(memberships), function(community){
        
        # Get community index
        community_index <- memberships != community
        
        # Return assigned proportion
        return(loadings[names(memberships)[community_index], as.character(community)] >= cut_off)
        
      }),
      na.rm = TRUE
    )
  )
  
}
