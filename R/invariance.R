#' @title Measurement Invariance of \code{\link[EGAnet]{EGA}} Structure
#'
#' @description Estimates configural invariance using \code{\link[EGAnet]{bootEGA}}
#' on all data (across groups) first. After configural variance is established,
#' then metric invariance is tested using the community structure that established
#' configural invariance (see \strong{Details} for more information on this process)
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param groups Numeric or character vector (length = \code{nrow(data)}).
#' Group membership corresponding to each case in data
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}.
#' If supplied, then configural invariance check is skipped (i.e., configural
#' invariance is assumed based on the given structure)
#'
#' @param iter Numeric (length = 1).
#' Number of iterations to perform for the permutation.
#' Defaults to \code{500} (recommended)
#'
#' @param configural.threshold Numeric (length = 1).
#' Value to use a threshold in \code{\link[EGAnet]{itemStability}} to determine
#' which items should be removed during configural invariance (see \strong{Details}
#' for more information).
#' Defaults to \code{0.70} (recommended)
#'
#' @param configural.type Character (length = 1).
#' Type of bootstrap to use for configural invariance in \code{\link[EGAnet]{bootEGA}}.
#' Defaults to \code{"parametric"}
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
#' @param ncores Numeric (length = 1).
#' Number of cores to use in computing results.
#' Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then type: \code{parallel::detectCores()}
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
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{EGA}},
#' \code{\link[EGAnet]{bootEGA}}, and
#' \code{\link[EGAnet]{net.loads}}
#'
#' @details In traditional psychometrics, measurement invariance is performed in
#' sequential testing from more flexible (more free parameters) to more rigid
#' (fewer free parameters) structures. Measurement invariance in network
#' psychometrics is no different.
#'
#' \strong{Configural Invariance}
#'
#' To establish configural invariance, the data are collapsed across groups
#' and a common sample structure is identified used \code{\link[EGAnet]{bootEGA}}
#' and \code{\link[EGAnet]{itemStability}}. If some variables have a replication
#' less than 0.70 in their assigned dimension, then they are considered unstable
#' and therefore not invariant. These variables are removed and this process
#' is repeated until all items are considered stable (replication values greater
#' than 0.70) or there are no variables left. If configural invariance cannot be
#' established, then the last run of results are returned and metric invariance
#' is not tested (because configural invariance is not met). Importantly, if any
#' variables \emph{are} removed, then configural invariance is not met for the
#' original structure. Any removal would suggest only partial configural invariance
#' is met.
#'
#' \strong{Metric Invariance}
#'
#' The variables that remain after configural invariance are submitted to metric
#' invariance. First, each group estimates a network and then network loadings
#' (\code{\link[EGAnet]{net.loads}}) are computed using the assigned
#' community memberships (determined during configural invariance). Then,
#' the difference between the assigned loadings of the groups is computed. This
#' difference represents the empirical values. Second, the group memberships
#' are permutated and networks are estimated based on the these permutated
#' groups for \code{iter} times. Then, network loadings are computed and
#' the difference between the assigned loadings of the group is computed, resulting
#' in a null distribution. The empirical difference is then compared against
#' the null distribution using a two-tailed \emph{p}-value based on the number
#' of null distribution differences that are greater and less than the empirical
#' differences for each variable. Both uncorrected and false discovery rate
#' corrected \emph{p}-values are returned in the results. Uncorrected \emph{p}-values
#' are flagged for significance along with the direction of group differences.
#'
#' \strong{Three or More Groups}
#' 
#' When there are 3 or more groups, the function performs metric invariance testing by comparing all possible pairs of groups. Specifically:
#' \itemize{
#' \item \emph{Pairwise Comparisons}: The function generates all possible unique group pairings and computes the differences in network loadings for each pair. The same community structure, derived from configural invariance or provided by the user, is used for all groups.
#' \item \emph{Permutation Testing}: For each group pair, permutation tests are conducted to assess the statistical significance of the observed differences in loadings. P-values are calculated based on the proportion of permuted differences that are greater than or equal to the observed difference.
#' \item \emph{Result Compilation}: The function compiles the results for each pair, including both uncorrected and FDR-corrected (Benjamini-Hochberg) p-values, and the direction of differences. It returns a summary of the findings for all pairwise comparisons.
#' }
#' This approach allows for a detailed examination of metric invariance across multiple groups, ensuring that all potential differences are thoroughly assessed while maintaining the ability to identify specific group differences.
#'
#' For more details, see Jamison, Golino, and Christensen (2023)
#'
#' @return Returns a list containing:
#'
#' \item{configural.results}{\code{\link[EGAnet]{bootEGA}} results from
#' the final run that produced configural invariance. This output will be
#' output on the final run of unsuccessful configural invariance runs}
#'
#' \item{memberships}{Original memberships provided in \code{structure}
#' or from \code{\link[EGAnet]{EGA}} if \code{structure = NULL}}
#'
#' \item{EGA}{Original \code{\link[EGAnet]{EGA}} results for the full sample}
#'
#' \item{groups}{A list containing:
#'
#' \itemize{
#'
#' \item \code{\link[EGAnet]{EGA}} --- \code{\link[EGAnet]{EGA}} results for each group
#'
#' \item \code{loadings} --- Network loadings (\code{\link[EGAnet]{net.loads}}) for each group
#'
#' \item \code{loadingsDifference} --- Difference between the dominant loadings of each group
#'
#' }
#'
#' }
#'
#' \item{permutation}{A list containing:
#'
#' \itemize{
#'
#' \item \code{groups} --- Permutated groups acorss iterations
#'
#' \item \code{loadings} --- Network loadings (\code{\link[EGAnet]{net.loads}}) for each group for each permutation
#'
#' \item \code{loadingsDifference} --- Difference between the dominant loadings of each group for each permutation
#'
#' }
#'
#' }
#'
#' \item{results}{Data frame of the results (which are printed)}
#'
#' @author Laura Jamison <lj5yn@virginia.edu>,
#' Hudson F. Golino <hfg9s at virginia.edu>, and
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#'
#' @references
#' \strong{Original implementation} \cr
#' Jamison, L., Golino, H., & Christensen, A. P. (2023).
#' Metric invariance in exploratory graph analysis via permutation testing.
#' \emph{PsyArXiv}.
#'
#' @examples
#' # Load data
#' wmt <- wmt2[-1,7:24]
#'
#' # Groups
#' groups <- rep(1:2, each = nrow(wmt) / 2)
#'
#' \dontrun{
#' # Measurement invariance
#' results <- invariance(wmt, groups, ncores = 2)
#'
#' # Plot with uncorrected alpha = 0.05
#' plot(results, p_type = "p", p_value = 0.05)
#'
#' # Plot with BH-corrected alpha = 0.10
#' plot(results, p_type = "p_BH", p_value = 0.10)}
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
#'
# Measurement Invariance
# Updated 10.08.2024
invariance <- function(
    data, groups, structure = NULL,
    iter = 500, configural.threshold = 0.70,
    configural.type = c("parametric", "resampling"),
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    ncores, seed = NULL, verbose = TRUE,
    search = TRUE,
    ...
)
{
  
  # Store random state (if there is one)
  store_state()
  
  # Check for missing arguments
  configural.type <- set_default(configural.type, "parametric", invariance)
  corr <- set_default(corr, "auto", invariance)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Argument errors (returns 'data' and 'groups')
  error_return <- invariance_errors(
    data, groups, iter, configural.threshold,
    ncores, seed, verbose
  )
  
  # Get data and groups
  data <- error_return$data; groups <- error_return$groups
  
  
  # Ensure data has variable names and get dimensions
  data <- ensure_dimension_names(data)
  dimensions <- dim(data)
  dimension_names <- dimnames(data)
  
  # Get unique groups
  unique_groups <- na.omit(unique(groups))
  num_iterations <- length(unique_groups)
  # Generate all possible group pairings
  group_pairs <- combn(unique_groups, 2, simplify = FALSE)
  num_iterations <- length(group_pairs)
  
  # Send message about total number of iterations
  message(paste("Metric invariance testing will be performed over", num_iterations, "group pair comparisons."))
  
  
  # If structure is supplied, then skip configural invariance
  if(is.null(structure)){
    
    # Send message
    message("Testing configural invariance...")
    
    # Perform configural invariance
    configural_results <- configural(
      data = data, iter = iter, structure = structure,
      configural.threshold = configural.threshold,
      configural.type = configural.type, corr = corr,
      na.data = na.data, model = model, algorithm = algorithm,
      uni.method = uni.method, ncores = ncores, seed = seed,
      verbose = verbose, ...
    )
    
    # Check for configural invariance
    if(!configural_results$configural_flag){
      
      # Send message
      message("\nConfigural invariance was not found. Terminating invariance testing...")
      
      # Return configural invariance results
      return(configural_results)
      
    }
    
    # Configural invariance was found, continue with metric
    
    # Send message
    message(paste(
      "\nConfigural invariance was found with",
      length(configural_results$stable_items), "variables\n"
    ))
    
    # Update data
    data <- configural_results$data
    
    # Update original EGA
    original_EGA <- configural_results$boot_object$EGA
    
    # Set structure based on original EGA
    structure <- original_EGA$wc
    
  } else {
    # Process structure if supplied
    structure <- remove_attributes(structure)
    structure <- force_vector(structure)
    names(structure) <- dimension_names[[2]]
  }
  
  # Get community names
  community_names <- as.character(unique(structure))
  
  # Send message about continuing on with metric invariance
  message("Testing metric invariance...")
  
  # Perform EGA and loadings calculations for all groups
  group_ega <- lapply(unique_groups, function(group){
    EGA(
      data = data[groups == group,],
      corr = corr, model = model,
      algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, ...
    )
  })
  names(group_ega) <- unique_groups
  
  group_loadings <- lapply(group_ega, function(x){
    loadings <- as.matrix(
      net.loads(A = x$network, wc = structure, ...)$std
    )
    return(loadings[dimension_names[[2]], community_names, drop = FALSE])
  })
  
  # Get seeds
  seeds <- reproducible_seeds(iter, seed)
  
  # Perform pairwise comparisons
  pairwise_results <- lapply(group_pairs, function(pair) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    # Original difference
    original_difference <- group_loadings[[as.character(group1)]] - group_loadings[[as.character(group2)]]
    
    # Obtain original assigned difference
    original_assigned_difference <- ulapply(
      community_names, function(community){
        original_difference[structure == community, community]
      }
    )
    
    # Ensure same order as original data
    original_assigned_difference <- original_assigned_difference[dimension_names[[2]]]
    
    # Permutate groups for this pair
    perm_groups <- lapply(
      seeds, function(seed_value){
        shuffle(groups[groups %in% pair], seed = seed_value)
      }
    )
    
    # Perform permutation estimation of loadings
    permutated_loadings <- parallel_process(
      iterations = iter,
      datalist = perm_groups,
      function(
    permutation, pair, structure,
    data, corr, model, algorithm, uni.method,
    ...
      ){
        # Estimate loadings
        return(
          lapply(pair, function(group){
            network <- EGA(
              data = data[permutation == group,],
              corr = corr, model = model,
              algorithm = algorithm, uni.method = uni.method,
              plot.EGA = FALSE, ...
            )$network
            
            loadings <- as.matrix(
              net.loads(A = network, wc = structure, ...)$std
            )
            
            return(loadings)
          })
        )
      },
    pair = pair, structure = structure,
    data = data[groups %in% pair,], corr = corr, model = model,
    algorithm = algorithm, uni.method = uni.method, ...,
    ncores = ncores, progress = verbose
    )
    
    # Compute differences (ensure same ordering)
    difference_list <- lapply(permutated_loadings, function(x){
      x[[1]][dimension_names[[2]], community_names, drop = FALSE] -
        x[[2]][dimension_names[[2]], community_names, drop = FALSE]
    })
    
    # Obtain assigned loadings only
    assigned_list <- lapply(difference_list, function(one_difference){
      differences <- ulapply(
        community_names, function(community){
          one_difference[structure == community, community]
        }
      )
      return(differences[dimension_names[[2]]])
    })
    
    # Create results
    permutation_counts <- lapply(assigned_list, function(x){
      abs(x) >= abs(original_assigned_difference)
    })
    
    # Replace the first permutation with all TRUE (original differences)
    permutation_counts[[1]] <- rep(TRUE, dimensions[2])
    
    # Compute p-values
    p_value <- rowMeans(
      do.call(cbind, permutation_counts), na.rm = TRUE
    )
    
    # Results data frame
    results_df <- data.frame(
      Membership = remove_attributes(structure),
      Difference = round(original_assigned_difference, 3),
      p = round(p_value, 3),
      p_BH = round(p.adjust(p_value, method = "BH"), 3)
    )
    
    # Order by dimension
    results_df <- results_df[order(results_df$Membership),]
    
    # Add significance
    sig <- swiftelse(results_df$p <= 0.10, ".", "")
    sig <- swiftelse(results_df$p <= 0.05, "*", sig)
    sig <- swiftelse(results_df$p <= 0.01, "**", sig)
    results_df$sig <- swiftelse(results_df$p <= 0.001, "***", sig)
    
    # Add direction
    direction <- paste(
      pair[1],
      swiftelse(sign(results_df$Difference) == 1, ">", "<"),
      pair[2]
    )
    results_df$Direction <- swiftelse(results_df$p <= 0.05, direction, "")
    
    return(list(
      pair = pair,
      results = results_df
    ))
  })
  
  # Compile results
  all_results <- list(
    configural.results = if(exists("configural_results")) configural_results else NULL,
    memberships = structure,
    EGA = if(exists("original_EGA")) original_EGA else NULL,
    groups = list(
      EGA = group_ega,
      loadings = group_loadings
    ),
    pairwise_results = pairwise_results
  )
  
  class(all_results) <- "invariance"
  
  return(all_results)
}

# Bug checking ----
## Basic input
# data = wmt2[-1, 7:24]; groups = rep(1:2, each = nrow(data) / 2)
# iter = 500; structure = NULL; configural.threshold = 0.70
# configural.type = "parametric"; corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"; uni.method = "louvain"
# ncores = 8; verbose = TRUE

#' @noRd
# Errors ----
# Updated 19.08.2023
invariance_errors <- function(
    data, groups, iter, configural.threshold,
    ncores, seed, verbose
)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "invariance")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'groups' errors
  object_error(groups, c("vector", "matrix", "data.frame"), "invariance")
  groups <- force_vector(groups)
  length_error(groups, dim(data)[1], "invariance")
  
  # 'iter' errors
  length_error(iter, 1, "invariance")
  typeof_error(iter, "numeric", "invariance")
  range_error(iter, c(1, Inf), "invariance")
  
  # 'configural.threshold' errors
  length_error(configural.threshold, 1, "invariance")
  typeof_error(configural.threshold, "numeric", "invariance")
  range_error(configural.threshold, c(0, 1), "invariance")
  
  # 'ncores' errors
  length_error(ncores, 1, "invariance")
  typeof_error(ncores, "numeric", "invariance")
  range_error(ncores, c(1, parallel::detectCores()), "invariance")
  
  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "invariance")
    typeof_error(seed, "numeric", "invariance")
    range_error(seed,  c(0, Inf), "invariance")
  }
  
  # 'verbose' errors
  length_error(verbose, 1, "invariance")
  typeof_error(verbose, "logical", "invariance")
  
  # Return usable data and groups
  return(list(data = usable_data(data, verbose), groups = groups))
  
}

#' @exportS3Method
# S3 Print Method ----
# Updated 10.08.2024
# Updated print method
print.invariance <- function(x, ...) {
  cat("Invariance analysis results:\n\n")
  
  for (result in x$pairwise_results) {
    cat("Comparison:", paste(result$pair, collapse = " vs "), "\n")
    print(result$results)
    cat("\n----\n")
  }
  
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1\n")
}

#' @exportS3Method
# S3 Summary Method ----
# Updated 10.08.2024
# Updated summary method
summary.invariance <- function(object, ...) {
  cat("Summary of invariance analysis:\n\n")
  
  cat("Number of groups:", length(object$groups$EGA), "\n")
  cat("Number of pairwise comparisons:", length(object$pairwise_results), "\n\n")
  
  for (result in object$pairwise_results) {
    cat("Comparison:", paste(result$pair, collapse = " vs "), "\n")
    cat("Number of noninvariant items (p < 0.05):", sum(result$results$p < 0.05), "\n")
    cat("Number of noninvariant items (p_BH < 0.05):", sum(result$results$p_BH < 0.05), "\n\n")
  }
  
  cat("Use print() for detailed results.\n")
}

#' @noRd
# Set group comparison plots ----
# Updated 10.08.2024
group_setup <- function(
    EGA_object, plot_ARGS,
    nodes, noninvariant,
    ...
)
{
  # Get ellipse arguments for defaults
  ellipse <- list(...)
  
  # Set edge size
  if(!"edge.size" %in% names(ellipse)){
    ellipse$edge.size <- 8 # default in `basic_plot_setup`
  }
  
  # Make copy of EGA object
  EGA_object_copy <- EGA_object
  
  # Ensure noninvariant is a numeric vector
  noninvariant <- as.numeric(as.vector(unlist(noninvariant)))
  
  # Set noninvariant edges to 1
  for(i in seq_len(nodes)){
    # Check for noninvariance
    if(noninvariant[i] == 1){
      # Get community index
      community_index <- EGA_object$wc == EGA_object$wc[i]
      
      # Target node
      target_node <- EGA_object$network[i, community_index]
      
      # Replace non-zero edges with 1
      EGA_object_copy$network[i, community_index] <-
        EGA_object_copy$network[community_index, i] <-
        ifelse(target_node != 0, 1, target_node)
    }
  }
  
  # Set up plot for edges
  edge_plot <- basic_plot_setup(
    network = EGA_object_copy$network,
    wc = EGA_object_copy$wc, ...,
    arguments = TRUE
  )
  
  # Based on significance...
  ## Change alpha (ensure numeric)
  plot_ARGS$node.alpha <- ifelse(noninvariant == 1, 0.75, 0.25)
  
  ## Change edge type
  plot_ARGS$edge.lty[
    edge_plot$ARGS$edge.size == ellipse$edge.size
  ] <- "dashed"
  
  # Set up EGA arguments
  plot_ARGS$net <- NULL
  plot_ARGS$network <- EGA_object$network
  plot_ARGS$wc <- EGA_object$wc
  
  # Return plot arguments
  return(plot_ARGS)
}


#' @exportS3Method
# S3 Plot Method ----
# Updated 10.08.2024
plot.invariance <- function(x, p_type = c("p", "p_BH"), p_value = 0.05, ...) {
  p_type <- match.arg(p_type)
  
  # Check for appropriate p-value range
  if (p_value <= 0 || p_value >= 1) {
    stop("p_value must be between 0 and 1")
  }
  
  # Get all unique groups
  all_groups <- unique(unlist(lapply(x$pairwise_results, function(r) r$pair)))
  
  # Combine all pairwise results
  all_results <- do.call(rbind, lapply(x$pairwise_results, function(r) {
    cbind(r$results, Group1 = r$pair[1], Group2 = r$pair[2])
  }))
  
  # Determine overall noninvariant items (as numeric)
  noninvariant <- sapply(unique(all_results$Item), function(item) {
    as.numeric(any(all_results$Item == item & all_results[, p_type] <= p_value))
  })
  
  # Get number of nodes
  nodes <- length(noninvariant)
  
  # Create a list to store all group plots
  group_plots <- list()
  
  for (group in all_groups) {
    # Set up group plot
    group_plot <- basic_plot_setup(
      network = x$groups$EGA[[group]]$network,
      wc = x$memberships,
      ...,
      arguments = TRUE
    )
    
    # Update plot for noninvariant nodes
    group_plot <- do.call(
      what = basic_plot_setup,
      args = group_setup(
        EGA_object = x$groups$EGA[[group]],
        plot_ARGS = group_plot$ARGS,
        nodes = nodes,
        noninvariant = noninvariant
      )
    )
    
    # Set up p-value title
    if (p_type == "p") {
      title <- bquote(paste("Group ", .(group), " (", italic(p), " ", .(ifelse(p_value < 0.5, "<", ">")), " ", .(p_value), ")"))
    } else {
      title <- bquote(paste("Group ", .(group), " (", italic(p)[adj.], " ", .(ifelse(p_value < 0.5, "<", ">")), " ", .(p_value), ")"))
    }
    
    # Update legend guide
    group_plot <- group_plot +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          title = title,
          title.position = "top",
          override.aes = list(
            alpha = c(0.25, 0.75),
            size = group_plot$ARGS$node.size,
            stroke = 1.5
          )
        )
      ) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
      )
    
    group_plots[[group]] <- group_plot
  }
  
  # Calculate the number of rows and columns for the grid
  n_groups <- length(all_groups)
  n_cols <- ceiling(sqrt(n_groups))
  n_rows <- ceiling(n_groups / n_cols)
  
  # Arrange all group plots in a grid
  ggpubr::ggarrange(
    plotlist = group_plots,
    ncol = n_cols,
    nrow = n_rows,
    common.legend = TRUE,
    legend = "bottom"
  )
}

#' @noRd
# Configural invariance ----
# Updated 13.04.2024
configural <- function(
    data, iter, structure, configural.threshold,
    configural.type, corr, na.data, model,
    algorithm, uni.method, ncores, seed, verbose,
    ...
){
  
  # Initialize check for items and stability values
  items <- dim(data)[2]
  stability <- FALSE
  
  # Perform `while` loop
  while(items > 0 & !stability){
    
    # Perform bootEGA
    boot <- bootEGA(
      data = data, corr = corr, na.data = na.data,
      model = model, algorithm = algorithm,
      uni.method = uni.method, iter = iter,
      type = configural.type, ncores = ncores,
      EGA.type = "EGA", typicalStructure = FALSE,
      plot.itemStability = FALSE,
      plot.typicalStructure = FALSE,
      seed = seed, verbose = verbose,
      clear = TRUE, suppress = TRUE, # additional internal arguments to `bootEGA`
      ...
    )
    
    # Get stable items
    stable_items <- boot$stability$item.stability$item.stability$empirical.dimensions >= configural.threshold
    
    # Perform checks
    stability <- all(stable_items)
    
    # Number of stable items
    items <- sum(stable_items)
    
    # Update data
    data <- data[,stable_items]
    
  }
  
  # Send results
  return(
    list(
      data = data,
      stable_items = stable_items,
      boot_object = boot,
      item_stability = boot$stability$item.stability,
      configural_flag = items > 0
    )
  )
  
}


