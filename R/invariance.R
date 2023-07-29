#' Measurement Invariance of \code{\link[EGAnet]{EGA}} Structure
#'
#' @description Estimates metric invariance of \code{\link[EGAnet]{EGA}} or specified structure
#' 
#' @param data Matrix or data frame.
#' Variables to be used in the analysis
#' 
#' @param groups Vector.
#' Group membership corresponding to each case in data
#' 
#' @param iter Numeric. 
#' Number of iterations to perform for the permutation.
#' Defaults to \code{500}
#'
#' @param memberships Vector. 
#' Node membership for each community or factor. 
#' Defaults to \code{NULL}.
#' When \code{NULL}, \code{\link[EGAnet]{EGA}} is used to compute node memberships
#'
#' @param type Character.
#' Type of measurement invariance to estimate.
#' Only includes \code{"loadings"} at the moment
#' 
#' @param model Character.
#' A string indicating the method to use.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[EGAnet]{TMFG}}. By default, \code{gamma} is set to 0 in 
#' \code{\link[EGAnet]{EBICglasso.qgraph}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#' 
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#' 
#' @param uni.method Character.
#' What unidimensionality method should be used? 
#' Defaults to \code{"louvain"}.
#' Current options are:
#' 
#' \itemize{
#'
#' \item{\strong{\code{expand}}}
#' {Expands the correlation matrix with four variables correlated .50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This is the method used in the Golino et al. (2020)
#' \emph{Psychological Methods} simulation.}
#'
#' \item{\strong{\code{LE}}}
#' {Applies the Leading Eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' \item{\strong{\code{louvain}}}
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix using a resolution parameter = 0.95.
#' If the number of dimensions is 1, then the Louvain solution is used; otherwise,
#' regular EGA is used. This method was validated in the Christensen (2022) simulation.}
#' 
#' }
#' 
#' @param consensus.iter Numeric.
#' Number of iterations to perform in consensus clustering for the Louvain algorithm
#' (see Lancichinetti & Fortunato, 2012).
#' Defaults to \code{100}
#' 
#' @param consensus.method Character.
#' What consensus clustering method should be used? 
#' Defaults to \code{"highest_modularity"}.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\strong{\code{highest_modularity}}}
#' {Uses the community solution that achieves the highest modularity
#' across iterations}
#' 
#' \item{\strong{\code{most_common}}}
#' {Uses the community solution that is found the most
#' across iterations}
#' 
#' \item{\strong{\code{iterative}}}
#' {Identifies the most common community solutions across iterations
#' and determines how often nodes appear in the same community together.
#' A threshold of 0.30 is used to set low proportions to zero.
#' This process repeats iteratively until all nodes have a proportion of
#' 1 in the community solution.
#' }
#' 
#' \item{\code{lowest_tefi}}
#' {Uses the community solution that achieves the lowest \code{\link[EGAnet]{tefi}}
#' across iterations}
#' 
#' }
#' 
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#' 
#' @param progress Boolean.
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' For Windows, \code{FALSE} is about 2x faster
#'
#' @return Returns a list containing:
#' 
#' \item{memberships}{Original memberships provided in \code{memberships}
#' or from \code{\link[EGAnet]{EGA}} if \code{NULL}}
#'
#' \item{EGA}{Original \code{\link[EGAnet]{EGA}} results for the sample}
#' 
#' \item{groups}{
#' 
#' \itemize{
#' 
#' \item{\code{\link[EGAnet]{EGA}}}
#' {EGA results for each group}
#' 
#' \item{\code{loadings}}
#' {Network loadings for each group}
#' 
#' \item{\code{loadingsDifference}}
#' {Difference between the dominant loadings of each group}
#' 
#' }
#' 
#' }
#' 
#' \item{permutation}{
#' 
#' \itemize{
#' 
#' \item{\code{groups}}
#' {Permutated groups acorss iterations}
#' 
#' \item{\code{loadings}}
#' {Loadings for each group for each permutation}
#' 
#' \item{\code{loadingsDifference}}
#' {Difference between the dominant loadings of each group for each permutation}
#' 
#' }
#' 
#' }
#' 
#' \item{results}{Data frame of the results (which are printed)}
#'
#' @author Laura Jamison <lj5yn@virginia.edu>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>, and 
#' Hudson F. Golino <hfg9s at virginia.edu>
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
#' results <- invariance(wmt, groups, ncores = 2)}
#' 
#' @export
#'
# Measurement Invariance
# Updated 25.07.2023
invariance <- function(
    data, groups, iter = 500, 
    structure = NULL, type = c("loadings"),
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    ncores, verbose = TRUE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
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
  
  # Get dimensions and dimension names of the data
  dimensions <- dim(data)
  dimension_names <- dimnames(data)
  
  # Send error if 'data' and 'groups' differ in length
  if(dimensions[1] != length(groups)){
    stop(
      "Number of cases in 'data' do not match the length of 'groups'. Please check that these numbers match: `nrow(data) == length(groups)`",
      call. = FALSE
    )
  }
  
  # Obtain original EGA
  original_EGA <- EGA(
    # Standard arguments
    data = data, corr = corr, model = model,
    algorithm = algorithm, uni.method = uni.method,
    plot.EGA = FALSE, ...
    # Legacy arguments 'model.args' and 'algorithm.args'
    # are handled in `EGA`
  )
  
  # Get ellipse arguments
  ellipse <- list(...)
  
  # Check for legacy argument 'memberships
  if("memberships" %in% names(ellipse)){
    structure <- ellipse$memberships
  }
  
  # Send error if 'structure' is not a vector, matrix, or data frame
  if(!is.null(structure)){ # If not NULL
    object_error(structure, c("vector", "matrix", "data.frame"))
    structure <- force_vector(structure)
  }else{ # Set structure based on original `EGA`
    structure <- original_EGA$wc
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
      net.loads(A = x$network, wc = structure)$std
    )
    
    # Reorder and return loadings
    return(loadings[dimension_names[[2]],, drop = FALSE])
    
  })
  
  # Ensure loadings are in the same order
  group_loadings[[2]] <- group_loadings[[2]][
    , dimnames(group_loadings[[1]])[[2]]
  ]
  
  # Original difference
  original_difference <- group_loadings[[1]] - group_loadings[[2]]
  
  # Get names of dimensions for loadings
  community_names <- dimnames(original_difference)[[2]]
  
  # Obtain original assigned difference
  original_assigned_difference <- ulapply(
    community_names, function(community){
      original_difference[structure == community, community]
    }
  )
  
  # Ensure same order as original data
  original_assigned_difference <- original_assigned_difference[
    dimension_names[[2]]
  ]
  
  # Permutate groups
  perm_groups <- lapply(seq_len(iter), function(i){
    shuffle(groups)
  })
  
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
# set.seed(1234)
# data = wmt2[,7:24]; dimensions = dim(data)
# dim_sequence = seq_len(dimensions[1])
# split1 = sample(dim_sequence, floor(dimensions[1] / 2))
# split2 = setdiff(dim_sequence, split1)
# group = c(rep(1, length(split1)), rep(2, length(split2)))
# data = data[c(split1, split2),]; iter = 500; type = "loadings"
# corr = "auto"; na.data = "pairwise"; model = "glasso"
# algorithm = "walktrap"; uni.method = "louvain"
# ncores = 8; verbose = FALSE; ellipse = list()

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
# Updated 28.07.2023
plot.invariance <- function(
    x, p_type = c("p", "p_BH"), p_value = 0.05, ...
)
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
    wc = x$groups$EGA[[1]]$wc, # ...,
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
  
  # Update `title` guide
  first_group$guides$colour$title <- "Invariant"
  second_group$guides$colour$title <- "Noninvariant"
  
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










