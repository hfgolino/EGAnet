#' @title Bootstrap Test for the Ergodicity Information Index
#'
#' @description Tests the Ergodicity Information Index obtained in the empirical sample 
#' with a distribution of EII obtained by bootstrap sampling
#' (see \strong{Details} for the procedure)
#'
#' @param dynEGA.object A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object. If a \code{\link[EGAnet]{dynEGA}}
#' object, then \code{level = c("individual", "population")} is required
#'
#' @param EII A \code{\link[EGAnet]{ergoInfo}} object 
#' used to estimate the Empirical Ergodicity Information Index 
#' or the estimated value of EII estimated using the \code{\link[EGAnet]{ergoInfo}} 
#' function. Inherits \code{use} from \code{\link[EGAnet]{ergoInfo}}.
#' If no \code{\link[EGAnet]{ergoInfo}} object is provided, then it is estimated
#' 
#' @param use Character (length = 1).
#' A string indicating what network element will be used
#' to compute the algorithm complexity, the list of edges or the weights of the network.
#' Defaults to \code{use = "unweighted"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item \code{"edge.list"} --- Calculates the algorithm complexity using the list of edges
#'
#' \item \code{"unweighted"} --- Calculates the algorithm complexity using the binary weights of the encoded prime 
#' transformed network. 0 = edge absent and 1 = edge present
#' 
#' \item \code{"weighted"} --- Calculates the algorithm complexity using the weights of encoded prime-weight transformed network
#' 
#' }
#' 
#' @param ordering Character (length = 1).
#' Changes ordering of edge list.
#' \code{"row"} goes across the rows;
#' \code{"column"} goes down the columns.
#' Defaults to \code{"row"}
#' 
#' @param shuffles Numeric.
#' Number of shuffles used to compute the Kolmogorov complexity.
#' Defaults to \code{5000}
#'
#' @param iter Numeric (length = 1).
#' Number of replica samples to generate from the bootstrap analysis.
#' Defaults to \code{200} (recommended)
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
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#' 
#' @details In traditional bootstrap sampling, individual participants are resampled with 
#' replacement from the empirical sample. This process is time consuming when carried out 
#' across \emph{v} number of variables, \emph{n} number of participants,
#' \emph{t} number of time points, and \emph{i} number of iterations.
#' 
#' The approach applied in \code{boot.ergoInfo}, is to obtain a sampling distribution 
#' of EII values as if all participants in the data have the population network structure. 
#' To mirror the sample of individuals, we generate as many population variants as there are participants 
#' in the empirical sample. With the new sample containing the population network variants, we compute EII 
#' with the original population network as the population network and the population network variants as the individuals. 
#' We repeat this process for \emph{X} iterations (e.g., 200). This approach creates a sampling distribution of EII that would 
#' be expected when the individuals in the population are deviations on the population structure -- that is, much of the 
#' population structure is retained but with different variations of noise in each individual. If the empirical EII is 
#' significant different than the generated distribution, then there is significant information lost when representing the sample 
#' as an aggregate, population network; otherwise, the system is determine to be ergodic and the sample can adequately be 
#' represented with the population network.
#' 
#' How to interpret the results: the result of \code{boot.ergoInfo} is a sampling distribution of EII values that would be expected if the process 
#' was ergodic (null distribution). If the empirical EII value is not significantly different from the null distribution, then  the empirical data can 
#' be expected to be generated from an ergodic process and the population structure is  sufficient to describe all 
#' individuals. If the empirical EII value is significantly different from the null distribution, 
#' then the empirical data cannot be described by the population structure -- significant 
#' information is lost when collapsing across to the population structure.
#'
#' @examples
#' # Obtain simulated data
#' sim.data <- sim.dynEGA
#' 
#' \dontrun{
#' # Dynamic EGA individual and population structures
#' dyn1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA[,-26], n.embed = 5, tau = 1,
#'   delta = 1, id = 25, use.derivatives = 1,
#'   model = "glasso", ncores = 2, corr = "pearson"
#' )
#'
#' # Empirical Ergodicity Information Index
#' eii1 <- ergoInfo(dynEGA.object = dyn1, use = "unweighted")
#'
#' # Bootstrap Test for Ergodicity Information Index
#' testing.ergoinfo <- boot.ergoInfo(
#'   dynEGA.object = dyn1, EII = eii1,
#'   ncores = 2, use = "unweighted"
#' )
#' 
#' # Plot result
#' plot(testing.ergoinfo)
#' 
#' # Example using `dynEGA`
#' dyn2 <- dynEGA(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, use.derivatives = 1, ncores = 2,
#'   level = c("individual", "population")
#' )
#' 
#' # Empirical Ergodicity Information Index
#' eii2 <- ergoInfo(dynEGA.object = dyn2, use = "unweighted")
#' 
#' # Bootstrap Test for Ergodicity Information Index
#' testing.ergoinfo2 <- boot.ergoInfo(
#'   dynEGA.object = dyn2, EII = eii2,
#'   ncores = 2
#' )
#' 
#' # Plot result
#' plot(testing.ergoinfo2)}
#'
#' @return Returns a list containing:
#' 
#' \item{empirical.ergoInfo}{Empirical Ergodicity Information Index}
#'
#' \item{boot.ergoInfo}{The values of the Ergodicity Information Index obtained in the bootstrap}
#'
#' \item{p.value}{The two-sided \emph{p}-value of the bootstrap test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII
#' with small variation in the population structure}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the bootstrap distribution of EII.}
#'
#' \item{interpretation}{How you can interpret the result of the test in plain English}
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @references 
#' \strong{Original Implementation} \cr
#' Golino, H., Nesselroade, J., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#' 
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
# Bootstrap Test for the Ergodicity Information Index
# Updated 07.11.2023
boot.ergoInfo <- function(
    dynEGA.object, EII, 
    use = c("edge.list", "unweighted", "weighted"),
    ordering = c("row", "column"), shuffles = 5000,
    iter = 200, ncores, verbose = TRUE
){
  
  # Send experimental message (for now)
  experimental("boot.ergoInfo")
  
  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "unweighted", ergoInfo)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Argument errors
  boot.ergoInfo_errors(dynEGA.object, iter, ncores, verbose)
  
  # Check for EII
  if(missing(EII)){ # If missing, then compute it
    EII <- ergoInfo(dynEGA.object, use = use, shuffles = shuffles)$EII
  }else if(is(EII, "EII")){
    
    # Get attributes
    use <- attr(EII, "methods")$use
    ordering <- attr(EII, "methods")$ordering
    shuffles <- attr(EII, "methods")$shuffles
    EII <- EII$EII # Save empirical EII for last
    
  }
  
  # Get proper objects (if not, send an error)
  # Function found in `ergoInfo`
  dynEGA.object <- get_dynEGA_object(dynEGA.object)
  
  # Only use necessary data (saves memory!)
  population_network <- dynEGA.object$population$network
  n_dimensions <- dynEGA.object$population$n.dim
  individual_sequence <- seq_along(dynEGA.object$individual)
  
  # Get lower triangle
  lower_triangle <- lower.tri(population_network)
  
  # Get rewire estimates
  rewire_estimates <- rewire_estimate(
    base = population_network,
    network_list = lapply(
      dynEGA.object$individual, function(x){x$network}
    )
  )
  
  # Determine range
  rewire_range <- range(rewire_estimates, na.rm = TRUE) / 2
  
  # Remove `dynEGA.object` from memory
  rm(dynEGA.object); clear_memory()
  
  # Get rewired networks
  rewired_EII <- parallel_process(
    iterations = iter,
    FUN = function(
      iteration,
      # dynEGA Arguments
      population_network = population_network,
      n_dimensions = n_dimensions,
      individual_sequence = individual_sequence,
      # EII Arguments
      use = use, ordering = ordering, shuffles = shuffles
    ){
      
      # Initialize `dynEGA` object structure
      rewired_dynEGA <- list(
        dynEGA = list(
          population = list(
            network = population_network,
            n.dim = n_dimensions
          ),
          individual = lapply( # Return as list named "network"
            individual_sequence, function(x){
              list(
                network = igraph_rewire(
                  network = population_network,
                  prob = runif_xoshiro(
                    1, min = rewire_range[1], 
                    max = rewire_range[2]
                  ),
                  noise = 0.05
                )
              )
            }
          )
        )
      )
      
      # Set class
      class(rewired_dynEGA) <- "dynEGA"
      
      # Return EII
      return(
        ergoInfo(
          rewired_dynEGA, use = use,
          ordering = ordering, shuffles = shuffles
        )
      )
      
    },
    # dynEGA Arguments
    population_network = population_network,
    n_dimensions = n_dimensions,
    individual_sequence = individual_sequence,
    # EII Arguments
    use = use, ordering = ordering, shuffles = shuffles,
    # Parallelization settings
    ncores = ncores, progress = verbose
  )
  
  # Get EII values
  EII_values <- nvapply(rewired_EII, function(x){x$EII})
  
  # Get indices for empirical EII greater than or equal to rewired values
  greater_than <- EII >= EII_values
  
  # Get p-value
  p_value = 2 * min(
    (sum(greater_than) + 1) / (iter + 1),
    (sum(EII <= EII_values) + 1) / (iter + 1)
  )
  
  # Get effect direction
  effect_direction <- swiftelse(
    p_value > 0.05, "n.s.",
    swiftelse(mean(greater_than) > 0.50, "greater", "less")
  )
  
  # Set up results to return
  results <- list(
    empirical.ergoInfo = EII,
    boot.ergoInfo = EII_values,
    p.value = p_value,
    effect = effect_direction,
    interpretation = switch(
      effect_direction,
      "n.s." = "The empirical EII was not different from values that would be expected if the process was ergodic, meaning the empirical data can be expected to be generated from an ergodic process and the population structure is sufficient to describe all individuals",
      "less" = "The empirical EII was significantly different from values that would be expected if the process was ergodic, meaning the empirical data cannot be described by the population structure -- significant information is lost when collapsing across to the population structure",
      "greater" = "The empirical EII was significantly different from values that would be expected if the process was ergodic, meaning the empirical data cannot be described by the population structure -- significant information is lost when collapsing across to the population structure"
    )
  )
  
  # Add "methods" attribute
  attr(results, "methods") <- list(
    use = use, ordering = ordering, shuffles = shuffles
  )
  
  # Set class
  class(results) <- "boot.ergoInfo"
  
  # Return results
  return(results)
  
}

#' @noRd
# Errors ----
# Updated 13.08.2023
boot.ergoInfo_errors <- function(dynEGA.object, iter, ncores, verbose)
{
  
  # 'dynEGA.object' errors ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "boot.ergoInfo")
  }
  
  # 'iter' errors
  length_error(iter, 1, "boot.ergoInfo")
  typeof_error(iter, "numeric", "boot.ergoInfo")
  range_error(iter, c(1, Inf), "boot.ergoInfo")
  
  # 'ncores' errors
  length_error(ncores, 1, "boot.ergoInfo")
  typeof_error(ncores, "numeric", "boot.ergoInfo")
  range_error(ncores, c(1, parallel::detectCores()), "boot.ergoInfo")
  
  # 'verbose' errors
  length_error(verbose, 1, "boot.ergoInfo")
  typeof_error(verbose, "logical", "boot.ergoInfo")
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 07.11.2023
print.boot.ergoInfo <- function(x, ...)
{
  
  # Print lower order
  cat(
    styletext(
      text = styletext(
        text =  "Empirical EII\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print EII method
  cat(
    paste0(
      "EII Method: ",
      switch(
        attr(x, "methods")$use,
        "edge.list" = "Edge List",
        "unweighted" = "Unweighted",
        "weighted" = "Weighted"
      ), "\n",
      "Ordering: ", totitle(attr(x, "methods")$ordering), "\n",
      "Shuffles: ", attr(x, "methods")$shuffles, "\n"
    )
  )
  
  # Print EII value
  cat("EII: ", round(x$empirical.ergoInfo, 4))
  
  # Add breakspace
  cat("\n\n")
  
  # Print higher order
  cat(
    styletext(
      text = styletext(
        text =  "Bootstrap EII\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print descriptives
  cat(
    paste0(
      "Mean = ", round(mean(x$boot.ergoInfo, na.rm = TRUE), 4),
      " (SD = ", round(sd(x$boot.ergoInfo, na.rm = TRUE), 4), ")",
      "\np-value = ", round(x$p.value, 4), " (", x$effect, ")",
      "\nErgodic: ", swiftelse(x$effect == "n.s.", "Yes", "No")
    )
  )
  
  cat("\n\n")
  
  # Print interpretation
  cat("Interpretation:\n", x$interpretation)

}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 26.07.2023
summary.boot.ergoInfo <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 09.07.2023
plot.boot.ergoInfo <- function(x, ...)
{
  
  # Send plot
  silent_plot(
    ggpubr::gghistogram(
      data = data.frame(EII = x$boot.ergoInfo),
      x = "EII", add = "mean", fill = "#00AFBB",
      color = "black", rug = TRUE, ylab = "Frequency",
      xlab = "Ergodicity Information Index",
      ...
    ) +
      ggplot2::geom_vline(
        xintercept = x$empirical.ergoInfo, color = "#00AFBB", linetype = "dotted" 
      )
  )
  
}

#' @noRd
# Estimate of rewiring
# Updated 11.07.2023
rewire_estimate <- function(base, network_list)
{
  
  # Get lower triangle based on base
  lower_triangle <- lower.tri(base)
  
  # Get base lower triangle
  base_lower <- base[lower_triangle]
  
  # Get binarized base
  base_lower[base_lower != 0] <- 1
  
  # Get indices that equal 1
  base_edge <- base_lower == 1
  
  # Get edges
  edges <- sum(base_edge)
  
  # Loop over network list
  return(
    nvapply(
      network_list, function(x){
        
        # Get lower triangle
        x <- x[lower_triangle]
        
        # Binarize network
        x[x != 0] <- 1
        
        # Compute proportion
        return(1 - sum(base_edge & x == 1) / edges)
        
      }
    )
  )
  
}
