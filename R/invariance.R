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
#' When there are 3 or more groups, the function performs metric invariance testing by comparing
#' all possible pairs of groups. Specifically:
#'
#' \itemize{
#'
#' \item \emph{Pairwise Comparisons}: The function generates all possible unique group pairings
#' and computes the differences in network loadings for each pair. The same community structure,
#' derived from configural invariance or provided by the user, is used for all groups.
#'
#' \item \emph{Permutation Testing}: For each group pair, permutation tests are conducted to
#' assess the statistical significance of the observed differences in loadings. \emph{p}-values are
#' calculated based on the proportion of permuted differences that are greater than or equal to
#' the observed difference.
#'
#' \item \emph{Result Compilation}: The function compiles the results for each pair including
#' both uncorrected (\code{p}) and FDR-corrected (Benjamini-Hochberg; \code{p_BH}) \emph{p}-values,
#' and the direction of differences. It returns a summary of the findings for all pairwise comparisons.
#'
#' }
#'
#' This approach allows for a detailed examination of metric invariance across multiple groups,
#' ensuring that all potential differences are thoroughly assessed while maintaining the ability
#' to identify specific group differences.
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
#' Jamison, L., Christensen, A. P., & Golino, H. F. (2024).
#' Metric invariance in exploratory graph analysis via permutation testing.
#' \emph{Methodology}, \emph{20}(2), 144-186.
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
# Updated 13.08.2024
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

  # Get unique groups (factored)
  unique_factors <- na.omit(unique(groups))

  # Set groups as factors
  groups <- as.numeric(factor(groups, levels = unique_factors))

  # Get unique groups (numeric)
  unique_groups <- na.omit(unique(groups))

  # Generate all possible group pairings
  group_pairs <- combn(unique_groups, 2, simplify = FALSE)
  pairs_length <- length(group_pairs)

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

    # Update dimension names
    dimension_names <- dimnames(data)

    # Update original EGA
    original_EGA <- configural_results$boot_object$EGA

    # Set structure based on original EGA
    structure <- original_EGA$wc

  }else{

    # Process structure if supplied
    structure <- remove_attributes(structure)
    structure <- force_vector(structure)
    names(structure) <- dimension_names[[2]]

  }

  # Get community names
  community_names <- as.character(unique(structure))

  # Send message about continuing on with metric invariance
  message("Testing metric invariance...")

  # Perform EGA for all groups
  group_ega <- lapply(unique_groups, function(group){
    EGA(
      data = data[groups == group,],
      corr = corr, model = model,
      algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, ...
    )
  }); names(group_ega) <- unique_factors # add names

  # Calculate loadings for all groups
  group_loadings <- lapply(group_ega, function(x){
    loadings <- as.matrix(
      net.loads(A = x$network, wc = structure, ...)$std
    )
    return(loadings[dimension_names[[2]], community_names, drop = FALSE])
  })

  # Get seeds
  seeds <- reproducible_seeds(iter, seed)

  # Get pairwise differences
  original_differences <- lapply(group_pairs, function(pair){

    # Original difference
    original_difference <- group_loadings[[pair[1]]] -
      group_loadings[[pair[2]]]

    # Obtain original assigned difference
    original_assigned_difference <- ulapply(
      community_names, function(community){
        original_difference[structure == community, community]
      }
    )

    # Ensure same order as original data
    return(original_assigned_difference[dimension_names[[2]]])

  }); names(original_differences) <- lapply(
    group_pairs, function(x){
      paste0(unique_factors[x[1]], "-", unique_factors[x[2]])
    }
  )

  # Permutate groups for this pair
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
            net.loads(A = network, wc = structure, ...)$std
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
  difference_list <- lapply(group_pairs, function(pair){

    # Loop over permutations
    permutated_differences <- lapply(permutated_loadings, function(x){

      x[[pair[1]]][dimension_names[[2]], community_names, drop = FALSE] -
      x[[pair[2]]][dimension_names[[2]], community_names, drop = FALSE]

    })

    # Obtain assigned loadings only
    assigned_list <- lapply(permutated_differences, function(one_difference){
      differences <- ulapply(
        community_names, function(community){
          one_difference[structure == community, community]
        }
      )
      return(differences[dimension_names[[2]]])
    })

    # Ensure same order as original data
    return(do.call(cbind, assigned_list))

  }); names(difference_list) <- names(original_differences)

  # Set up pairwise results
  results_list <- lapply(seq_len(pairs_length), function(i){

    # Replace the first permutation with all TRUE (original differences)
    permutation_counts <- cbind(
      TRUE, abs(difference_list[[i]]) >= abs(original_differences[[i]])
    )[,-2]

    # Compute p-values
    p_value <- rowMeans(permutation_counts, na.rm = TRUE)

    # Results data frame
    results_df <- data.frame(
      Membership = remove_attributes(structure),
      Difference = round(original_differences[[i]], 3),
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
      unique_factors[group_pairs[[i]][1]],
      swiftelse(sign(results_df$Difference) == 1, ">", "<"),
      unique_factors[group_pairs[[i]][2]]
    )
    results_df$Direction <- swiftelse(results_df$p <= 0.05, direction, "")

    # Return data frames
    return(results_df)

  }); names(results_list) <- names(original_differences)

  # Results list
  results <- list(
    configural.results = swiftelse(
      exists("configural_results"),
      configural_results, NULL
    ),
    memberships = structure,
    EGA = swiftelse(
      exists("original_EGA"),
      original_EGA, NULL
    ),
    groups = list(
      EGA = group_ega,
      loadings = group_loadings,
      loadingsDifference = original_differences,
      unique_groups = unique_factors
    ),
    permutation = list(
      groups = perm_groups,
      loadings = permutated_loadings,
      loadingsDifference = difference_list
    ),
    results = results_list
  )

  # Add class
  class(results) <- "invariance"

  # Restore random state (if there is one)
  restore_state()

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
# Updated 13.08.2024
# Updated print method
print.invariance <- function(x, pairs = list(), ...)
{

  # Print title
  cat(styletext("Invariance Results", defaults = "bold"), "\n")

  # Get number of possible pairs
  total_pairs <- length(x$results)
  pairs_sequence <- seq_len(total_pairs)

  # Create combination of pairs
  possible_pairs <- combn(names(x$groups$EGA), 2, simplify = FALSE)

  # Check for pairs
  input_pairs <- length(pairs) == 0

  # Check for missing pairs
  if(input_pairs){
    pairs <- possible_pairs
  }

  # Match up pairs
  combined <- do.call(rbind, c(pairs, possible_pairs))
  pairs <- which(duplicated(combined)) - length(pairs)

  # Loop over to print results
  for(pair in pairs){

    # Print groups
    cat(
      styletext(
        paste(
          "\nComparison:", gsub("-", " vs ", names(x$results)[pair])
        ), defaults = "underline"
      ), "\n"
    )

    # Print "as-is" results
    print(x$results[[pair]])

    # Add breakspace
    cat("----\n")

  }

  # Print significance codes
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 'n.s.' 1\n\n")

  # Let user know about print
  if(input_pairs && total_pairs > 1){
    cat("Use the argument 'pairs = list()' for individual paired results using `c()` inside `list()` (for multiple pairs, use `c()` inside `list()` for each pair).\n")
  }

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 13.08.2024
summary.invariance <- function(object, ...)
{

  # Print title
  cat(styletext("Summary of Invariance Results", defaults = "bold"), "\n\n")

  # Print descriptives
  cat("Number of groups:", length(object$groups$EGA), "\n")
  cat("Number of pairwise comparisons:", length(object$results), "\n\n")

  # Loop over for p-value descriptives
  for(pair in seq_along(object$results)){

    # Print groups
    cat(
      styletext(
        paste(
          "\nComparison:", gsub("-", " vs ", names(object$results)[pair])
        ), defaults = "underline"
      ), "\n"
    )

    # Print noninvariant items
    cat(
      "Number of noninvariant items (p < 0.05):",
      sum(object$results[[pair]]$p < 0.05),
      paste0(
        "(",
        format_decimal(mean(object$results[[pair]]$p < 0.05) * 100, 1),
        "%)\n"
      )
    )
    cat(
      "Number of noninvariant items (p_BH < 0.05):",
      sum(object$results[[pair]]$p_BH < 0.05),
      paste0(
        "(",
        format_decimal(mean(object$results[[pair]]$p_BH < 0.05) * 100, 1),
        "%)\n\n"
      )
    )

  }

  # Let user know about print
  cat("Use `print()` for more detailed results.\n")

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
# Updated 13.08.2024
plot.invariance <- function(x, pairs = list(), p_type = c("p", "p_BH"), p_value = 0.05, ...)
{

  # Obtain unique groups
  unique_factors <- x$groups$unique_groups

  # Set default for p-value type
  p_type <- swiftelse(missing(p_type), "p", match.arg(p_type))

  # Check for appropriate p-value range
  range_error(p_value, c(0, 1), "plot.invariance")

  # Get number of possible pairs
  total_pairs <- length(x$results)
  pairs_sequence <- seq_len(total_pairs)

  # Create combination of pairs
  possible_pairs <- combn(names(x$groups$EGA), 2, simplify = FALSE)

  # Check for pairs
  input_pairs <- length(pairs) == 0

  # Check for missing pairs
  if(input_pairs){
    pairs <- possible_pairs
  }else if(length(pairs) > 6){ # Don't do more than 6 comparisons (equivalent to 4 groups)

    # Get first six comparisons
    pairs <- pairs[1:6]

    # Send warning message
    warning(
      paste(
        "Due to the complexity of plotting more than 6 pairwise plots,",
        "only the combinations involving the first six comparisons will be used:\n",
        paste(
          lapply(
            pairs, function(x){
              paste(unique_factors[x[1]], "vs", unique_factors[x[2]])
            }
          ), collapse = ", "
        )
      )
    )

  }

  # Match up pairs
  combined <- do.call(rbind, c(pairs, possible_pairs))
  pair_index <- which(duplicated(combined)) - length(pairs)

  # Get number of groups
  group_names <- names(x$groups$EGA)
  total_groups <- length(group_names)
  group_sequence <- seq_len(total_groups)
  total_pairs <- length(pair_index)

  # Ensure same memberships
  x$groups$EGA <- lapply(group_sequence, function(i){

    # Replace memberships with 'structure'
    x$groups$EGA[[i]]$wc <- x$memberships

    # Return entire object
    return(x$groups$EGA[[i]])

  }); names(x$groups$EGA) <- group_names

  # Obtain noninvariant items
  noninvariant <- lapply(x$results, function(result){
    result[names(x$memberships), p_type] <= p_value
  })

  # Get number of nodes
  nodes <- length(noninvariant[[1]])

  # Set up first group plot
  first_group <- basic_plot_setup(
    network = x$groups$EGA[[
      possible_pairs[[pair_index[1]]][1]
    ]]$network,
    wc = x$groups$EGA[[
      possible_pairs[[pair_index[1]]][1]
    ]]$wc,  ...,
    arguments = TRUE
  )

  # Get plot arguments
  second_ARGS <- first_group$ARGS

  # Remove some arguments from `first_ARGS`
  ## Essentially, the same call but allows some freedom
  second_ARGS[c(
    "net", "edge.alpha", "edge.color", "edge.lty", "edge.size"
  )] <- NULL
  first_ARGS <- second_ARGS

  # Set up p-value title
  if(p_type == "p"){
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

  # Determine best final arrangement
  rows <- switch(
    as.character(total_pairs),
    "1" = 1, "2" = 2, "3" = 3,
    "4" = 2, "5" = 3, "6" = 3
  )
  columns <- switch(
    as.character(total_pairs),
    "1" = 1, "2" = 1, "3" = 1,
    "4" = 2, "5" = 2, "6" = 2
  )

  # Set up for first indices
  first_indices <- seq_len(rows)

  # Set up pairwise plots
  pairwise_plots <- lapply(first_indices, function(index){

    # Add network and memberships
    first_ARGS$network <- x$groups$EGA[[
      possible_pairs[[pair_index[[index]]]][1]
    ]]$network
    first_ARGS$wc <- x$groups$EGA[[
      possible_pairs[[pair_index[[index]]]][1]
    ]]$wc
    first_ARGS$arguments <- TRUE
    second_ARGS$network <- x$groups$EGA[[
      possible_pairs[[pair_index[[index]]]][2]
    ]]$network
    second_ARGS$wc <- x$groups$EGA[[
      possible_pairs[[pair_index[[index]]]][2]
    ]]$wc
    second_ARGS$arguments <- TRUE

    # Set up group plots
    first_group <- do.call(basic_plot_setup, first_ARGS)
    second_group <- do.call(basic_plot_setup, second_ARGS)

    # Get updated plots for each group
    ## First group
    first_group <- do.call(
      what = basic_plot_setup,
      args = group_setup(
        EGA_object = x$groups$EGA[[
          possible_pairs[[pair_index[[index]]]][1]
        ]],
        plot_ARGS = first_group$ARGS,
        nodes = nodes, noninvariant = noninvariant[[pair_index[[index]]]]
      )
    )
    ## Second group
    second_group <- do.call(
      what = basic_plot_setup,
      args = group_setup(
        EGA_object = x$groups$EGA[[
          possible_pairs[[pair_index[[index]]]][2]
        ]],
        plot_ARGS = second_group$ARGS,
        nodes = nodes, noninvariant = noninvariant[[pair_index[[index]]]]
      )
    )

    # Check for last group
    if(index == rows){

      # Update legend guide
      first_group <- first_group +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            title = invariant_title,
            title.position = "top",
            override.aes = list(
              alpha = 0.25, size = second_ARGS$node.size,
              stroke = 1.5
            )
          )
        )
      second_group <- second_group +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            title = noninvariant_title,
            title.position = "top",
            override.aes = list(
              alpha = 0.75, size = second_ARGS$node.size,
              stroke = 1.5
            )
          )
        )

      # Adjust size and position
      first_group <- first_group +
        ggplot2::theme(
          legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
        )
      second_group <- second_group +
        ggplot2::theme(
          legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
        )

      # Return plot
      return(
        ggpubr::ggarrange(
          first_group, second_group,
          ncol = 2, nrow = 1,
          labels = possible_pairs[[pair_index[[index]]]],
          legend = "bottom",
          common.legend = FALSE
        )
      )

    }else{

      # Remove legends
      first_group <- first_group + ggplot2::theme(
        legend.position = "none"
      )
      second_group <- second_group + ggplot2::theme(
        legend.position = "none"
      )

      # Return plot
      return(
        ggpubr::ggarrange(
          first_group, second_group,
          ncol = 2, nrow = 1,
          labels = possible_pairs[[pair_index[[index]]]],
          legend = "none"
        )
      )

    }

  })

  # Arrange final plots
  final_plots <- ggpubr::ggarrange(
    plotlist = pairwise_plots,
    ncol = 1, nrow = rows
  )

  # Check for more
  if(total_pairs > 3){

    # Set up for second indices
    second_indices <- (rows + 1):total_pairs

    # Set up pairwise plots
    pairwise_plots <- lapply(second_indices, function(index){

      # Add network and memberships
      first_ARGS$network <- x$groups$EGA[[
        possible_pairs[[pair_index[[index]]]][1]
      ]]$network
      first_ARGS$wc <- x$groups$EGA[[
        possible_pairs[[pair_index[[index]]]][1]
      ]]$wc
      first_ARGS$arguments <- TRUE
      second_ARGS$network <- x$groups$EGA[[
        possible_pairs[[pair_index[[index]]]][2]
      ]]$network
      second_ARGS$wc <- x$groups$EGA[[
        possible_pairs[[pair_index[[index]]]][2]
      ]]$wc
      second_ARGS$arguments <- TRUE

      # Set up group plots
      first_group <- do.call(basic_plot_setup, first_ARGS)
      second_group <- do.call(basic_plot_setup, second_ARGS)

      # Get updated plots for each group
      ## First group
      first_group <- do.call(
        what = basic_plot_setup,
        args = group_setup(
          EGA_object = x$groups$EGA[[
            possible_pairs[[pair_index[[index]]]][1]
          ]],
          plot_ARGS = first_group$ARGS,
          nodes = nodes, noninvariant = noninvariant[[pair_index[[index]]]]
        )
      )
      ## Second group
      second_group <- do.call(
        what = basic_plot_setup,
        args = group_setup(
          EGA_object = x$groups$EGA[[
            possible_pairs[[pair_index[[index]]]][2]
          ]],
          plot_ARGS = second_group$ARGS,
          nodes = nodes, noninvariant = noninvariant[[pair_index[[index]]]]
        )
      )

      # Check for last group
      if(index == total_pairs){

        # Update legend guide
        first_group <- first_group +
          ggplot2::guides(
            colour = ggplot2::guide_legend(
              title = invariant_title,
              title.position = "top",
              override.aes = list(
                alpha = 0.25, size = second_ARGS$node.size,
                stroke = 1.5
              )
            )
          )
        second_group <- second_group +
          ggplot2::guides(
            colour = ggplot2::guide_legend(
              title = noninvariant_title,
              title.position = "top",
              override.aes = list(
                alpha = 0.75, size = second_ARGS$node.size,
                stroke = 1.5
              )
            )
          )

        # Adjust size and position
        first_group <- first_group +
          ggplot2::theme(
            legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
          )
        second_group <- second_group +
          ggplot2::theme(
            legend.title = ggplot2::element_text(size = 12, hjust = 0.5)
          )

        # Return plot
        return(
          ggpubr::ggarrange(
            first_group, second_group,
            ncol = 2, nrow = 1,
            labels = possible_pairs[[pair_index[[index]]]],
            legend = "bottom",
            common.legend = FALSE
          )
        )

      }else{

        # Remove legends
        first_group <- first_group + ggplot2::theme(
          legend.position = "none"
        )
        second_group <- second_group + ggplot2::theme(
          legend.position = "none"
        )

        # Return plot
        return(
          ggpubr::ggarrange(
            first_group, second_group,
            ncol = 2, nrow = 1,
            labels = possible_pairs[[pair_index[[index]]]],
            legend = "none"
          )
        )

      }

    })

    # Arrange final plots
    final_plots <- ggpubr::ggarrange(
      final_plots, ggpubr::ggarrange(
        plotlist = pairwise_plots,
        ncol = 1, nrow = rows
      ), ncol = columns, nrow = 1
    )

  }

  # Let user know about plot
  if(input_pairs && total_pairs > 1){
    cat("Use the argument 'pairs = list()' for individual paired results using `c()` inside `list()` (for multiple pairs, use `c()` inside `list()` for each pair).\n")
  }

  # Return final plots
  return(final_plots)

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


