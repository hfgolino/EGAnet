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
#' @param memberships Vector. 
#' Node membership for each community or factor. 
#' Defaults to \code{NULL}.
#' When \code{NULL}, \code{\link[EGAnet]{EGA}} is used to compute node memberships
#'
#' @param type Character.
#' Type of measurement invariance to estimate.
#' Only includes \code{"loadings"} at the moment
#' 
#' @param iter Numeric. 
#' Number of iterations to perform for the permutation.
#' Defaults to \code{500}
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
#' \dontshow{# Fast for CRAN
#' # Measurement invariance
#' results <- invariance(wmt, groups, iter = 1, memberships = ega.wmt$wc, ncores = 2)
#' }
#'
#' \donttest{
#' # Measurement invariance
#' results <- invariance(wmt, groups, ncores = 2)
#' }
#' 
#' @export
#'
# Measurement Invariance
# Updated 22.04.2022
invariance <- function(
  data, groups, 
  memberships = NULL,
  type = c("loadings"),
  iter = 500, ncores, ...
)
{
  # Number of processing cores
  if(missing(ncores)){
    ncores <- round(parallel::detectCores() / 2, 0)
  }
  
  # Obtain additional arguments
  add_args <- list(...)
  
  # Obtain default EGA arguments
  ega_args <- formals(EGA.estimate)
  
  # Set defaults for EGA
  ega_args$model <- "glasso" # network estimation
  ega_args$algorithm <- "walktrap" # community detection
  ega_args$corr <- "cor_auto" # correlation estimation
  ega_args$verbose <- FALSE # leave out verbose
  ega_args$model.args <- list(gamma = 0) # Set gamma to zero to maximize similarities
  ega_args$... <- NULL
  
  # Check for additional arguments for EGA
  if(length(add_args) != 0){
    
    # Match arguments
    arg_names <- names(ega_args)[na.omit(match(names(add_args), names(ega_args)))]
    
    # Input arguments
    ega_args[arg_names] <- add_args[arg_names]
    
  }
  
  # Make sure data and groups match
  if(nrow(data) != length(groups)){
    stop("Number of cases in 'data' do not match the length of 'groups'. Please check that these numbers match: `nrow(data) == length(groups)`")
  }
  
  # Add data to EGA arguments
  ega_args$data <- data
  
  # Ensure class is list!
  class(ega_args) <- "list"
  
  # Estimate original EGA
  original_EGA <- suppressWarnings(
    suppressMessages(
      do.call(
        what = EGA.estimate,
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
          what = EGA.estimate,
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
    loadings <- net.loads(
      A = x$network, wc = memberships
    )$std
    
    # Reorder loadings
    loadings <- loadings[colnames(data),]
    
    # Return loadings
    return(loadings)
    
  })
  
  ## Reorder order loadings to match group 1
  group_loadings <- lapply(group_loadings, function(x){
    x[,colnames(group_loadings[[1]])]
  })
  
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
  
  # # Export variables (only necessary for testing)
  # # Comment out for package
  # parallel::clusterExport(
  #   cl = cl,
  #   varlist = c(
  #     "EGA.estimate",
  #     "net.loads",
  #     "unique_groups",
  #     "perm_groups"
  #   )
  # )
  
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
          what = EGA.estimate,
          args = ega_args
        )$network
        
        # Obtain loadings
        loadings <- net.loads(A = network, wc = memberships)$std
        
        # Reorder
        loadings <- loadings[colnames(data),]
        
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