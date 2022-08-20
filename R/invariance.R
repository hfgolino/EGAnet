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
#' or \code{\link[EGAnet]{TMFG}}
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
#' @param ... Arguments passed to \code{\link[EGAnet]{EGA}}
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
#' \donttest{
#' # Measurement invariance
#' results <- invariance(wmt, groups, ncores = 2)}
#' 
#' @export
#'
# Measurement Invariance
# Updated 20.08.2022
invariance <- function(
  data, groups, iter = 500, 
  memberships = NULL,
  type = c("loadings"),
  corr = c("cor_auto", "pearson", "spearman"),
  uni.method = c("expand", "LE", "louvain"),
  model = c("glasso", "TMFG"), model.args = list(),
  algorithm = c("walktrap", "leiden", "louvain"), algorithm.args = list(),
  consensus.method = c(
    "highest_modularity",
    "most_common",
    "iterative",
    "lowest_tefi"
  ), consensus.iter = 100, 
  ncores
)
{
  # Number of processing cores
  if(missing(ncores)){
    ncores <- round(parallel::detectCores() / 2, 0)
  }
  
  # Missing arguments
  
  if(missing(model)){
    model <- "glasso"
  }else{model <- tolower(match.arg(model))}
  
  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }
  
  if(missing(consensus.method)){
    consensus.method <- "most_common"
  }else{consensus.method <- tolower(match.arg(consensus.method))}
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- tolower(match.arg(corr))}
  
  # Model function
  model.FUN <- switch(
    model,
    "glasso" = EBICglasso.qgraph,
    "tmfg" = TMFG
  )
  
  # Model arguments
  model.ARGS <- obtain.arguments(
    FUN = model.FUN,
    FUN.args = model.args
  )
  
  # Force gamma = 0
  model.ARGS$gamma <- 0
  
  # Algorithm function
  if(!is.function(algorithm)){
    algorithm.FUN <- switch(
      algorithm,
      "walktrap" = igraph::cluster_walktrap,
      "leiden" = igraph::cluster_leiden,
      "louvain" = igraph::cluster_louvain
    )
  }else{
    algorithm.FUN <- algorithm
  }
  
  # Algorithm arguments
  algorithm.ARGS <- obtain.arguments(
    FUN = algorithm.FUN,
    FUN.args = algorithm.args
  )
  
  ## Remove weights from igraph functions' arguments
  if("weights" %in% names(algorithm.ARGS)){
    algorithm.ARGS[which(names(algorithm.ARGS) == "weights")] <- NULL
  }
  
  # Make sure data and groups match
  if(nrow(data) != length(groups)){
    stop("Number of cases in 'data' do not match the length of 'groups'. Please check that these numbers match: `nrow(data) == length(groups)`")
  }
  
  # Add data to EGA arguments
  ega_args$data <- data
  
  # Ensure class is list!
  class(ega_args) <- "list"
  
  # Set EGA arguments
  ega_args$corr <- corr
  ega_args$model <- model
  ega_args$model.args <- model.ARGS
  ega_args$algorithm <- algorithm
  ega_args$algorithm.args <- algorithm.ARGS
  ega_args$uni.method <- uni.method
  ega_args$consensus.method <- consensus.method
  ega_args$plot.EGA <- FALSE
  
  # Estimate original EGA
  original_EGA <- suppressWarnings(
    suppressMessages(
      do.call(
        what = EGA,
        args = ega_args
      )
    )
  )
  
  # Obtain memberships if missing
  if(is.null(memberships)){
    memberships <- original_EGA$wc
  }
  
  # Obtain unique groups
  unique_groups <- na.omit(unique(groups))
  
  # Obtain original group EGAs
  group_ega <- lapply(seq_along(unique_groups), function(i){
    
    # Obtain group data
    ega_args$data <- data[groups == unique_groups[i],]
    
    # Obtain network
    suppressWarnings(
      suppressMessages(
        do.call(
          what = EGA,
          args = ega_args
        )
      )
    )
    
  })
  
  # Rename list
  names(group_ega) <- unique_groups
  
  # Original network loadings
  group_loadings <- lapply(group_ega, function(x){
    
    # Obtain loadings
    loadings <- as.matrix(net.loads(
      A = x$network, wc = memberships
    )$std)
    
    # Reorder loadings
    loadings <- loadings[colnames(data),]
    
    # Check for vector
    if(is.vector(loadings)){
      loadings <- matrix(loadings, ncol = 1)
      colnames(loadings) <- 1
      row.names(loadings) <- colnames(data)
    }
    
    # Return loadings
    return(loadings)
    
  })
  
  ## Reorder order loadings to match group 1
  if(ncol(group_loadings[[1]]) != 1){
    group_loadings <- lapply(group_loadings, function(x){
      x[,colnames(group_loadings[[1]])]
    })
  }
  
  # Original difference
  original_difference <- group_loadings[[1]] - group_loadings[[2]]
  
  # Obtain original dominant difference
  original_dominant_difference <- unlist(
    lapply(seq_along(memberships), function(i){
      original_difference[i,as.character(memberships[i])]
    })
  )
  names(original_dominant_difference) <- colnames(data)
  
  # Permutate groups
  perm_groups <- lapply(1:iter, function(i){
    sample(groups, size = length(groups), replace = FALSE)
  })
  
  # Message for estimating permutated loadings
  message("Performing permutations...")
  
  # Set up parallelization
  cl <- parallel::makeCluster(ncores)
  
  # Export variables (only necessary for testing)
  # Comment out for package
  parallel::clusterExport(
    cl = cl,
    varlist = c(
      "EGA", "ega_args",
      "memberships",
      "net.loads",
      "unique_groups",
      "perm_groups",
      "data"
    )
  )
  
  # Obtain permutated loadings
  loadings_list <- pbapply::pblapply(
    X = seq_along(perm_groups),
    FUN = function(i){
      
      # Estimate loadings
      loadings_groups <- lapply(seq_along(unique_groups), function(j){
        
        # Insert permutated data
        ega_args$data <- data[which(perm_groups[[i]] == unique_groups[j]),]
        
        # Obtain network
        network <- do.call(
          what = EGA,
          args = ega_args
        )$network
        
        # Obtain loadings
        loadings <- net.loads(A = network, wc = memberships)$std
        
        # Reorder loadings
        loadings <- loadings[colnames(data),]
        
        # Check for vector
        if(is.vector(loadings)){
          loadings <- matrix(loadings, ncol = 1)
          colnames(loadings) <- 1
          row.names(loadings) <- colnames(data)
        }
        
        # Return loadings
        return(loadings)
        
      })
      
      # Name groups
      names(loadings_groups) <- unique_groups
      
      # Return EGA groups
      return(loadings_groups)
      
    },
    cl = cl
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  # Compute differences
  difference_list <- lapply(loadings_list, function(x){
    
    # Ensure same ordering of communities
    x[[2]] <- x[[2]][,colnames(x[[1]])]
    
    # Obtain difference between groups
    difference_loadings <- x[[1]] - x[[2]]
    
    # Return difference loadings
    return(difference_loadings)
    
  })
  
  # Obtain dominant loadings only
  dominant_list <- lapply(difference_list, function(x){
    
    # Obtain differences for dominant loadings
    dominant_difference <- unlist(
      lapply(seq_along(memberships), function(i){
        x[i,as.character(memberships[i])]
      })
    )
    
    # Rename
    names(dominant_difference) <- row.names(x)
    
    # Return dominant differences
    return(dominant_difference)
    
  })
  
  # Create results
  permutation_counts <- lapply(dominant_list, function(x){
    abs(x) >= abs(original_dominant_difference)
  })
  
  ## Simplify to matrix
  permutation_counts <- simplify2array(permutation_counts)
  
  ## Add a column of TRUE (original difference)
  permutation_counts[,1] <- rep(TRUE, nrow(permutation_counts))
  
  ## p-value
  p_value <- rowMeans(permutation_counts, na.rm = TRUE)
  
  # Results data frame
  results_df <- data.frame(
    Node = colnames(data),
    Membership = memberships,
    Difference = original_dominant_difference,
    p = p_value
  )
  
  # Order by dimension
  results_df <- results_df[order(results_df$Membership),]
  
  # Add significance
  sig <- ifelse(results_df$p <= .10, ".", "n.s.")
  sig <- ifelse(results_df$p <= .05, "*", sig)
  sig <- ifelse(results_df$p<= .01, "**", sig)
  sig <- ifelse(results_df$p <= .001, "***", sig)
  results_df$sig <- sig
  
  # Results list
  results <- list()
  results$memberships <- memberships
  results$EGA <- original_EGA
  results$groups$EGA <- group_ega
  results$groups$loadings <- group_loadings
  results$groups$loadingsDifference <- original_dominant_difference
  results$permutation$groups <- perm_groups
  results$permutation$loadings <- loadings_list
  results$permutation$loadingsDifference <- dominant_list
  results$results <- results_df
  
  # Add class
  class(results) <- "invariance"
  
  # Return results
  return(results)
  
}
