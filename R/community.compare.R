#' @title Compares Community Detection Solutions Using Permutation
#'
#' @description A permutation implementation to determine statistical
#' significance of whether the community comparison measure is different
#' from zero
#'
#' @param base Character or numeric vector.
#' A vector of characters or numbers that are treated as the 
#' baseline communities
#' 
#' @param comparison Character or numeric vector (length = \code{length(base)}).
#' A vector of characters or numbers that are treated as the 
#' baseline communities
#' 
#' @param method Character (length = 1).
#' Comparison metrics from \code{\link[igraph]{compare}}.
#' Defaults to \code{"adjusted.rand"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item \code{"vi"} --- Variation of information (Meila, 2003)
#' 
#' \item \code{"nmi"} --- Normalized mutual information (Danon et al., 2003)
#' 
#' \item \code{"split.join"} --- Split-join distance (Dongen, 2000)
#' 
#' \item \code{"rand"} --- Rand index (Rand, 1971)
#' 
#' \item \code{"adjusted.rand"} --- adjusted Rand index (Hubert & Arabie, 1985; Steinley, 2004)
#' 
#' }
#' 
#' @param iter Numeric (length = 1).
#' Number of permutations to perform.
#' Defaults to \code{1000} (recommended)
#' 
#' @return Returns data frame containing method used (\code{Method}), empirical or observed
#' value (\code{Empirical}), and p-value based on the permutation test (\code{p.value})
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate network
#' network <- EBICglasso.qgraph(data = wmt)
#' 
#' # Compute Edge Betweenness
#' edge_between <- community.detection(network, algorithm = "edge_betweenness")
#' 
#' # Compute Fast Greedy
#' fast_greedy <- community.detection(network, algorithm = "fast_greedy")
#' 
#' # Perform permutation test
#' community.compare(edge_between, fast_greedy)
#'
#' @references
#' \strong{Implementation of Permutation Test} \cr
#' Qannari, E. M., Courcoux, P., & Faye, P. (2014). 
#' Significance test of the adjusted Rand index. Application to the free sorting task. 
#' \emph{Food Quality and Preference}, \emph{32}, 93–97.
#'
#' \strong{Variation of Information} \cr
#' Meila, M. (2003, August).
#' Comparing clusterings by the variation of information. 
#' In \emph{Learning Theory and Kernel Machines: 16th Annual Conference on Learning Theory and 7th Kernel Workshop}, 
#' COLT/Kernel 2003, Washington, DC, USA, August 24-27, 2003. Proceedings (pp. 173-187). Berlin, DE: Springer Berlin Heidelberg.
#' 
#' \strong{Normalized Mutual Information} \cr
#' Danon, L., Diaz-Guilera, A., Duch, J., & Arenas, A. (2005). 
#' Comparing community structure identification. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2005}(09), P09008.
#' 
#' \strong{Split-join Distance} \cr
#' Dongen, S. (2000). 
#' Performance criteria for graph clustering and Markov cluster experiments. 
#' \emph{CWI (Centre for Mathematics and Computer Science)}.
#' 
#' \strong{Rand Index} \cr
#' Rand, W. M. (1971). 
#' Objective criteria for the evaluation of clustering methods. 
#' \emph{Journal of the American Statistical Association}, \emph{66}(336), 846-850.
#' 
#' \strong{Adjusted Rand Index} \cr
#' Hubert, L., & Arabie, P. (1985). 
#' Comparing partitions. 
#' \emph{Journal of Classification}, \emph{2}, 193-218.
#'
#' Steinley, D. (2004). 
#' Properties of the Hubert-Arabie adjusted rand index. 
#' \emph{Psychological Methods}, \emph{9}(3), 386.
#'
#' @export
#'
# Perform permutations for clusters ----
# Updated 03.12.2023
community.compare <- function(
    base, comparison,
    method = c(
      "vi", "nmi", "split.join", "rand", "adjusted.rand"
    ), iter = 1000, verbose = TRUE
)
{
  
  # Check for missing arguments (argument, default, function)
  method <- set_default(method, "adjusted.rand", community.compare)
  
  # Argument errors (returns usable solutions)
  output <- community.compare_errors(base, comparison, iter, verbose)
  
  # Assign updated 'base' and 'comparison' values
  base <- output$base; comparison <- output$comparison
  
  # Obtain empirical value
  empirical_value <- igraph::compare(base, comparison, method = method)
  
  # Obtain sequence of elements
  element_sequence <- seq_len(output$elements)
  
  # Perform permutations
  permutation_values <- nvapply(
    seq_len(iter - 1), function(i){
      
      # Compute comparison
      igraph::compare(
        base[shuffle(element_sequence)],
        comparison[shuffle(element_sequence)],
        method = method
      )
      
    }
  )
  
  # Replace first permutation value with empirical value
  # permutation_values[1] <- empirical_value
  
  # Return results
  return(
    data.frame(
      Method = switch(
        method,
        "vi" = "Variation of Information", 
        "nmi" = "Normalized Mutual Information", 
        "split.join" = "Split Join", 
        "rand" = "Rand Index", 
        "adjusted.rand" = "Adjusted Rand Index"
      ),
      Empirical = empirical_value,
      p.value = mean(c(TRUE, permutation_values > empirical_value), na.rm = TRUE)
    )
  )
  
}

# Bug Checking ----
# ## Basic input
# # Estimate network
# network <- EBICglasso.qgraph(data = wmt2[,7:24])
# base <- community.detection(network, algorithm = "edge_betweenness")
# comparison <- community.detection(network, algorithm = "fast_greedy")

#' @noRd
# Errors ----
# Updated 03.12.2023
community.compare_errors <- function(base, comparison, iter, verbose)
{
  
  # 'base' errors
  object_error(base, c("vector", "matrix", "data.frame"), "community.compare")
  
  # Ensure vector
  base <- as.vector(base)
  
  # Continue 'base' errors
  typeof_error(base, c("numeric", "character"), "community.compare")
  
  # Determine whether conversion to numeric is necessary
  if(is.character(base)){
    base <- as.numeric(factor(base))
  }
  
  # Get length of base
  base_length <- length(base)
  
  # 'comparison' errors
  object_error(comparison, c("vector", "matrix", "data.frame"), "community.compare")
  
  # Ensure vector
  comparison <- as.vector(comparison)
  
  # Continue 'comparison' errors
  length_error(comparison, base_length, "community.compare")
  typeof_error(comparison, c("numeric", "character"), "community.compare")
  
  # Determine whether conversion to numeric is necessary
  if(is.character(comparison)){
    comparison <- as.numeric(factor(comparison))
  }
  
  # 'iter' errors
  length_error(iter, 1, "community.compare")
  typeof_error(iter, "numeric", "community.compare")
  
  # 'verbose' errors
  length_error(verbose, 1, "community.compare")
  typeof_error(verbose, "logical", "community.compare")
  
  # Error checks complete
  
  # Check for available values
  base_available <- !is.na(base)
  comparison_available <- !is.na(comparison)
  available_values <- base_available & comparison_available 
  na_values <- !available_values
  
  # Error on all NAs
  if(all(na_values)){
    
    # Determine which is all NAs
    all_NA <- swiftelse(all(!base_available), "base", "comparison")
    
    # Set error
    stop(
      paste0("Input into '", all_NA, "' had all NA values after conversion to numeric. Double check your input"),
      call. = FALSE
    )
    
  }
  
  # Set warning on variables with NA 
  if(verbose && any(na_values)){
    
    # Send warning
    warning(
      paste0(
        "Some values were NA. These indices were removed: ",
        paste0(which(na_values), collapse = ", ")
      ), call. = FALSE
    )
    
  }
  
  # Return 'base' and 'comparison'
  return(
    list(
      base = base[available_values],
      comparison = comparison[available_values],
      elements = sum(available_values)
    )
  )
  
}
