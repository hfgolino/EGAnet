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
#' @param groups Numeric or character vector (length = \code{ncol(data)}).
#' Group membership corresponding to each case in data
#' 
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}
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
#' {See \code{\link[EGAnet]{cluster_walktrap}} for more details}
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
#' At this time, only two groups are supported. There is a method proposed to
#' test three or more groups in Jamison, Golino, and Christensen (2023) but
#' this approach has not been thoroughly vetted and validated. Future versions
#' of the package will provide support for three or more groups once there is
#' an established consensus for best practice.
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
#' \item{\code{\link[EGAnet]{EGA}} --- }
#' {\code{\link[EGAnet]{EGA}} results for each group}
#' 
#' \item{\code{loadings} --- }
#' {Network loadings (\code{\link[EGAnet]{net.loads}}) for each group}
#' 
#' \item{\code{loadingsDifference} --- }
#' {Difference between the dominant loadings of each group}
#' 
#' }
#' 
#' }
#' 
#' \item{permutation}{A list containing:
#' 
#' \itemize{
#' 
#' \item{\code{groups} --- }
#' {Permutated groups acorss iterations}
#' 
#' \item{\code{loadings} --- }
#' {Network loadings (\code{\link[EGAnet]{net.loads}}) for each group for each permutation}
#' 
#' \item{\code{loadingsDifference} --- }
#' {Difference between the dominant loadings of each group for each permutation}
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
# Updated 03.08.2023
invariance <- function(
    data, groups, structure = NULL,
    iter = 500, configural.threshold = 0.70,
    configural.type = c("parametric", "resampling"),
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    ncores, seed = NULL, verbose = TRUE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
  configural.type <- set_default(configural.type, "parametric", invariance)
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- swiftelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Send error if 'groups' is not a vector, matrix, or data frame
  object_error(groups, c("vector", "matrix", "data.frame"))
  groups <- force_vector(groups)
  
  # Ensure data has variable names
  data <- ensure_dimension_names(data)
  
  # Get dimensions and dimension names of the data
  original_dimensions <- dim(data)
  original_dimension_names <- dimnames(data)
  
  # Send error if 'data' and 'groups' differ in length
  if(original_dimensions[1] != length(groups)){
    stop(
      "Number of cases in 'data' do not match the length of 'groups'. Please check that these numbers match: `nrow(data) == length(groups)`",
      call. = FALSE
    )
  }
  
  # Get ellipse arguments
  ellipse <- list(...)
  
  # Check for seed
  if(is.null(seed)){ # Send message about seed
    message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results\n")
  }
  
  # Check for legacy argument 'memberships
  if("memberships" %in% names(ellipse)){
    structure <- ellipse$memberships
  }
  
  # Get unique groups
  unique_groups <- na.omit(unique(groups))
  
  # Send warning about only two groups (for now)
  if(length(unique_groups) > 2){
    
    # Send warning
    warning(
      "More than two groups is not yet supported. Using the first two groups...",
      call. = FALSE
    )
    
    # Update unique groups
    unique_groups <- unique_groups[c(1L, 2L)]
    
    # Keep indices
    keep_index <- groups %in% unique_groups
    
    # Update groups
    groups <- groups[keep_index]
    
    # Update data
    data <- data[keep_index,]
    
    # Update data dimensions
    dimensions <- dim(data)
    
  }
  
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
    length(configural_results$stable_items), "variables", 
    "\n\nTesting metric invariance..."
  ))

  # Update data
  data <- configural_results$data
  
  # Update dimensions and dimension names of the data
  dimensions <- dim(data)
  dimension_names <- dimnames(data)
  
  # Update original EGA
  original_EGA <- configural_results$boot_object$EGA
  
  # Send error if 'structure' is not a vector, matrix, or data frame
  if(!is.null(structure)){
    
    # Remove attributes
    structure <- remove_attributes(structure)
    
    # If not `NULL`, then check for error in object type
    object_error(structure, c("vector", "matrix", "data.frame"))
    
    # Make sure 'structure' is a vector
    structure <- force_vector(structure)
    
    # Make sure 'structure' has names
    names(structure) <- original_dimension_names[[2]]
    
    # Reduce based on stable items
    structure <- structure[names(configural_results$stable_items)]
    
  }else{ # Set structure based on original `EGA`
    structure <- original_EGA$wc
  }

  # Get community names
  community_names <- as.character(unique(structure))
  
  # Obtain original group EGAs
  group_ega <- lapply(unique_groups, function(group){
    
    # Return `EGA`
    EGA(
      data = data[groups == group,], 
      corr = corr, model = model,
      algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, ...
    )
    
  })
  
  # Rename list
  names(group_ega) <- unique_groups
  
  # Original network loadings
  group_loadings <- lapply(group_ega, function(x){
    
    # Obtain loadings
    loadings <- as.matrix(
      net.loads(A = x$network, wc = structure, ...)$std
    )
    
    # Reorder and return loadings
    return(loadings[dimension_names[[2]], community_names])
    
  })
  
  # Original difference
  original_difference <- group_loadings[[1]] - group_loadings[[2]]
  
  # Obtain original assigned difference
  original_assigned_difference <- ulapply(
    community_names, function(community){
      original_difference[structure == community, community]
    }
  )
  
  # Ensure same order as original data
  original_assigned_difference <- original_assigned_difference[dimension_names[[2]]]
  
  # Get seeds
  seeds <- reproducible_seeds(iter, seed)
  
  # Permutate groups
  perm_groups <- lapply(
    seeds, function(seed_value){
      shuffle(groups, seed = seed_value)
    }
  )
  
  # Perform permutation estimation of loadings
  permutated_loadings <- parallel_process(
    iterations = iter,
    datalist = perm_groups,
    function(
      permutation, unique_groups, structure,
      data, corr, model, algorithm, uni.method,
      ...
    ){
      
      # Estimate loadings
      return( # By groups
        lapply(unique_groups, function(group){
          
          
          # Get network
          network <- EGA(
            data = data[permutation == group,], 
            corr = corr, model = model,
            algorithm = algorithm, uni.method = uni.method,
            plot.EGA = FALSE, ...
          )$network
          
          # Obtain loadings
          loadings <- as.matrix(
            net.loads(A = network, wc = structure)$std
          )
          
          # Return loadings
          return(loadings)
          
        })
      )
      
    }, 
    # Make sure the additional objects get in
    unique_groups = unique_groups, structure = structure,
    data = data, corr = corr, model = model,
    algorithm = algorithm, uni.method = uni.method, ...,
    ncores = ncores, progress = verbose
  )
  
  # Compute differences (ensure same ordering)
  difference_list <- lapply(permutated_loadings, function(x){
    x[[1]][dimension_names[[2]], community_names] -
      x[[2]][dimension_names[[2]], community_names]
  })
  
  # Obtain assigned loadings only
  assigned_list <- lapply(difference_list, function(one_difference){
    
    # Get differences
    differences <- ulapply(
      community_names, function(community){
        one_difference[structure == community, community]
      }
    )
    
    # Ensure same order as original data
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
    unique_groups[1],
    swiftelse(sign(results_df$Difference) == 1, ">", "<"),
    unique_groups[2]
  )
  results_df$Direction <- swiftelse(results_df$p <= 0.05, direction, "")
  
  # Results list
  results <- list(
    configural.results = configural_results,
    memberships = structure,
    EGA = original_EGA,
    groups = list(
      EGA = group_ega,
      loadings = group_loadings,
      loadingsDifference = original_assigned_difference
    ),
    permutation = list(
      groups = perm_groups,
      loadings = permutated_loadings,
      loadingsDifference = assigned_list
    ),
    results = results_df
  )
  
  # Add class
  class(results) <- "invariance"
  
  # Return results
  return(results)
  
}

# Bug checking ----
## Basic input
# data = wmt2[-1, 7:24]; groups = rep(1:2, each = nrow(data) / 2)
# iter = 500; structure = NULL; configural.threshold = 0.70
# configural.type = "parametric"; corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"; uni.method = "louvain"
# ncores = 8; verbose = TRUE

#' @exportS3Method 
# S3 Print Method ----
# Updated 10.07.2023
print.invariance <- function(x, ...) {
  
  # Print results "as-is"
  print(x$results)
  cat("----\n") # Add breakspace and significance code
  cat("Signif. code: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1")
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 10.07.2023
summary.invariance <- function(object, ...) {
  print(object, ...) # same as print
}

#' @noRd
# Set group comparison plots ----
# Updated 24.07.2023
group_setup <- function(
    EGA_object, plot_ARGS,
    nodes, noninvariant,
    ...
)
{
  
  # Get ellipse arguments for defaults
  ellipse <- list(...)
  
  # Set edge size
  if(!"edge.size" %in% ellipse){
    ellipse$edge.size <- 8 # default in `basic_plot_setup`
  }
  
  # Make copy of EGA object
  EGA_object_copy <- EGA_object
  
  # Set noninvariant edges to 1
  for(i in seq_len(nodes)){
    
    # Check for noninvariance
    if(noninvariant[i]){
      
      # Get community index
      community_index <- EGA_object$wc == EGA_object$wc[i]
      
      # Target node
      target_node <- EGA_object$network[i, community_index]
      
      # Replace non-zero edges with 1
      EGA_object_copy$network[i, community_index] <-
        EGA_object_copy$network[community_index, i] <-
        swiftelse(target_node != 0, 1, target_node)
      
    }
    
  }
  
  # Set up plot for edges
  edge_plot <- basic_plot_setup(
    network = EGA_object_copy$network, 
    wc = EGA_object_copy$wc, ..., 
    arguments = TRUE
  )

  # Based on significance...
  ## Change alpha
  plot_ARGS$node.alpha <- swiftelse(noninvariant, 0.75, 0.25)
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
# Updated 04.08.2023
plot.invariance <- function(x, p_type = c("p", "p_BH"), p_value = 0.05, ...)
{
  
  # Set default for p-value type
  p_type <- swiftelse(missing(p_type), "p", match.arg(p_type))
  
  # Check for appropriate p-value range
  range_error(p_value, c(0, 1))
  
  # Ensure same memberships
  x$groups$EGA[[1]]$wc <- x$EGA$wc
  x$groups$EGA[[2]]$wc <- x$EGA$wc
  
  # Obtain noninvariant items
  noninvariant <- x$results[names(x$EGA$wc), p_type] <= p_value
  
  # Get number of nodes
  nodes <- length(noninvariant)

  # Set up first group plot
  first_group <- basic_plot_setup(
    network = x$groups$EGA[[1]]$network,
    wc = x$groups$EGA[[1]]$wc,  ...,
    arguments = TRUE
  )
  
  # Get plot arguments
  second_ARGS <- first_group$ARGS
  
  # Remove some arguments from `first_ARGS`
  ## Essentially, the same call but allows some freedom
  second_ARGS[c(
    "net", "edge.alpha", "edge.color", "edge.lty", "edge.size"
  )] <- NULL
  
  # Add network and memberships
  second_ARGS$network <- x$groups$EGA[[2]]$network
  second_ARGS$wc <- x$groups$EGA[[2]]$wc
  second_ARGS$arguments <- TRUE
  
  # Set up second group plot
  second_group <- do.call(basic_plot_setup, second_ARGS)
  
  # Get updated plots for each group
  ## First group
  first_group <- do.call(
    what = basic_plot_setup,
    args = group_setup(
      EGA_object = x$groups$EGA[[1]],
      plot_ARGS = first_group$ARGS,
      nodes = nodes, noninvariant = noninvariant
    )
  )
  ## Second group
  second_group <- do.call(
    what = basic_plot_setup,
    args = group_setup(
      EGA_object = x$groups$EGA[[2]],
      plot_ARGS = second_group$ARGS,
      nodes = nodes, noninvariant = noninvariant
    )
  )
  
  # Update `alpha` guide
  first_group$guides$colour$override.aes$alpha <- 0.25
  second_group$guides$colour$override.aes$alpha <- 0.75
  
  # Set up p-value title
  if(p_type == "p") {
    invariant_title <- bquote(
      paste("Invariant (", italic(p), " > ", .(p_value), ")")
    )
    noninvariant_title <- bquote(
      paste("Noninvariant (", italic(p), " < ", .(p_value), ")")
    )
  }else{
    invariant_title <- bquote(
      paste("Invariant (", italic(p)[adj.], " > ", .(p_value), ")")
    )
    noninvariant_title <- bquote(
      paste("Noninvariant (", italic(p)[adj.], " < ", .(p_value), ")")
    )
  }
  
  # Update `title` guide
  first_group$guides$colour$title <- invariant_title
    
  second_group$guides$colour$title <- noninvariant_title
  
  # Update `title.position` guide
  first_group$guides$colour$title.position <- "top"
  second_group$guides$colour$title.position <- "top"
  
  # Adjust size and position
  first_group <- first_group +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  second_group <- second_group +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )

  # Arrange plots
  ggpubr::ggarrange(
    first_group, second_group,
    ncol = 2, nrow = 1,
    labels = names(x$groups$EGA),
    legend = "bottom",
    common.legend = FALSE
  )
  
}

#' @noRd
# Configural invariance ----
# Updated 03.08.2023
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
      plot.typicalStructure = FALSE,
      seed = seed, verbose = verbose, 
      clear = TRUE, suppress = TRUE, # additional internal arguments to `bootEGA`
      ...
    )
      
    # Perform itemStability
    item_stability <- itemStability(boot, IS.plot = FALSE, structure = structure)
    
    # Stable items
    stable_items <- item_stability$item.stability$empirical.dimensions >= configural.threshold
    
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
      item_stability = item_stability,
      configural_flag = items > 0
    )
  )
  
}








