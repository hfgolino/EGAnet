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
#' Defaults to \code{use = "edge.list"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\code{"edge.list"} --- }
#' {Calculates the algorithm complexity using the list of edges.}
#'
#' \item{\code{"unweighted"} --- }
#' {Calculates the algorithm complexity using the binary weights of the network.
#' 0 = edge absent and 1 = edge present}
#' 
#' \item{\code{"weighted"} --- }
#' {Calculates the algorithm complexity using the weights of encoded prime-weight transformed network
#' where the weights are multiplied by 10 and rounded (to force integers)}
#' 
#' }
#'
#' @param iter Numeric (length = 1).
#' Number of replica samples to generate from the bootstrap analysis.
#' Defaults to \code{100} (recommended)
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
#' A more efficient process, the approach applied here, is to obtain a sampling distribution 
#' of EII values as if all participants in the data have the population network structure. 
#' Sampling is not perfect and therefore random noise is added to the edges of the population 
#' structure to simulate sampling variability. This noise follows a random uniform distribution
#' ranging from -0.10 to 0.10. In addition, a proportion of edges are rewired to allow for 
#' slight variations on the population structure. The proportion of nodes that are rewired is 
#' sampled from a random uniform distribution between 0.20 to 0.40. This process is carried out 
#' for each participant resulting in \emph{n} variations of the population structure. 
#' Afterward, EII is computed. This process is carried out for \emph{i} iterations (e.g., 100).
#' 
#' The result is a sampling distribution of EII values that would be expected if the process 
#' was ergodic. If the empirical EII value is significantly less than the distribution or 
#' not significantly different, then  the empirical data can be expected to be generated 
#' from an ergodic process and the population structure is  sufficient to describe all 
#' individuals. If the empirical EII value is significantly greater than the distribution, 
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
#' eii1 <- ergoInfo(dynEGA.object = dyn1, use = "edge.list")
#'
#' # Bootstrap Test for Ergodicity Information Index
#' testing.ergoinfo <- boot.ergoInfo(
#'   dynEGA.object = dyn1, EII = eii1,
#'   ncores = 2
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
#' eii2 <- ergoInfo(dynEGA.object = dyn2, use = "edge.list")
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
# Updated 03.08.2023
boot.ergoInfo <- function(
    dynEGA.object, EII, 
    use = c("edge.list", "unweighted", "weighted"),
    iter = 100, ncores, verbose = TRUE
){
  
  # Send experimental message (for now)
  experimental("boot.ergoInfo")
  
  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "edge.list", ergoInfo)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Argument errors
  boot.ergoInfo_errors(dynEGA.object, iter, ncores, verbose)
  
  # Check for EII
  if(missing(EII)){ # If missing, then compute it
    EII <- ergoInfo(dynEGA.object, use = use)$EII # , seed = 0)$EII
  }else if(is(EII, "EII")){
    use <- attr(EII, "methods")$use; EII <- EII$EII
  }
  
  # Get proper objects (if not, send an error)
  # Function found in `ergoInfo`
  dynega_objects <- get_dynEGA_object(dynEGA.object)
  
  # Replace individual networks with population networks
  individual_networks <- lapply(
    dynega_objects$individual, function(x){dynega_objects$population$network}
  )
  
  # Get lower triangle indices (avoids repeated computation)
  lower_triangle <- lower.tri(dynega_objects$population$network)
  
  # Get rewired networks
  rewired_networks <- lapply(
    seq_len(iter), function(iteration){
      
      # Initialize `dynEGA` object structure
      rewired_dynEGA <- list(
        dynEGA = list(
          population = list(
            network = dynega_objects$population$network,
            n.dim = dynega_objects$population$n.dim
          ),
          individual = lapply( # Return as list named "network"
            individual_networks, function(x){
              list(
                network = rewire(
                  network = x, min = 0.20, max = 0.40,
                  noise = 0.10, lower_triangle = lower_triangle
                )
              )
            }
          )
        )
      )
      
      # Set class
      class(rewired_dynEGA) <- "dynEGA"
      
      # Return rewired networks
      return(rewired_dynEGA)
      
    }
  )
  
  # Perform parallelization
  rewired_EII <- parallel_process(
    iterations = iter,
    datalist = rewired_networks,
    ergoInfo, use = use,
    ncores = ncores,
    progress = verbose
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
      "n.s." = "The empirical EII was not different from what would be expected from random variation in the population structure, meaning non-significant information is lost when aggregating the results into a single, population network.",
      "less" = "The empirical EII was less than what would be expected from random variation in the population structure, meaning non-significant information is lost when aggregating the results into a single, population network.",
      "greater" = "The empirical EII was greater than what would be expected from random variation in the population structure, meaning significant information is lost when aggregating the results into a single, population network."
    )
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
# Updated 26.07.2023
print.boot.ergoInfo <- function(x, ...)
{
  
  # Message about print support
  message("No print support yet")
  
  # Print x
  x
  
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
# Rewire networks ----
# About 10x faster than previous implementation
# Updated 30.07.2023
rewire <- function(
    network, min = 0.20, max = 0.40,
    noise = 0.10, lower_triangle
)
{
  
  # Work only with the lower triangle
  lower_network <- network[lower_triangle]
  
  # Get non-zero edges
  non_zero_edges <- which(lower_network != 0)
  
  # Number of edges
  edges <- length(non_zero_edges)
  
  # Add noise
  if(!is.null(noise)){
    
    # Only add to existing edges
    lower_network[non_zero_edges] <- 
      lower_network[non_zero_edges] + runif(edges, -noise, noise)
    
  }
  
  # Number of edges to rewire
  rewire_edges <- floor(edges * runif(1, min, max))
  
  # Get rewiring indices
  rewire_index <- shuffle(non_zero_edges, size = rewire_edges)
 
  # Get replacement indices
  replace_index <- shuffle(seq_along(lower_network), size = rewire_edges)
 
  # Make a copy of the lower network
  lower_network_original <- lower_network
  
  # Replace values
  lower_network[rewire_index] <- lower_network_original[replace_index]
  lower_network[replace_index] <- lower_network_original[rewire_index]
  
  # Get nodes in original network
  nodes <- dim(network)[2]
  
  # Initialize a new network
  new_network <- matrix(
    0, nrow = nodes, ncol = nodes,
    dimnames = dimnames(network)
  )
  
  # Replace values
  new_network[lower_triangle] <- lower_network
  new_network <- t(new_network)
  new_network[lower_triangle] <- lower_network
  
  # Return the rewired network
  return(new_network)
  
}


