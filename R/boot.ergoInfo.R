#' Bootstrap Test for the Ergodicity Information Index
#'
#' @description Tests the Ergodicity Information Index obtained in the empirical sample with a distribution of EII 
#' obtained by bootstrap sampling. In traditional bootstrap sampling, individual participants are resampled with
#' replacement from the empirical sample. This process is time consuming when carried out across \emph{v} number
#' of variables, \emph{n} number of participants, \emph{t} number of time points, and \emph{i} number of iterations.
#' 
#' A more efficient process, the approach applied here, is to obtain a sampling distribution of EII values as if
#' all participants in the data have the population network structure. Sampling is not perfect and therefore
#' random noise is added to the edges of the population structure to simulate sampling variability. This noise
#' follows a random uniform distribution ranging from -0.10 to 0.10. In addition, a proportion of edges are
#' rewired to allow for slight variations on the population structure. The proportion of nodes that are rewired
#' is sampled from a random uniform distribution between 0.20 to 0.40. This process is carried out for each
#' participant resulting in \emph{n} variations of the population structure. Afterward, EII is computed. This
#' process is carried out for \emph{i} iterations (e.g., 100).
#' 
#' The result is a sampling distribution of EII values that would be expected if the process was ergodic. If
#' the empirical EII value is significantly less than the distribution or not significantly different, then 
#' the empirical data can be expected to be generated from an ergodic process and the population structure is 
#' sufficient to describe all individuals. If the empirical EII value is significantly greater than the distribution,
#' then the empirical data cannot be described by the population structure -- significant information is lost when
#' collapsing across to the population structure.
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object that is used to match the arguments of the EII object.
#'
#' @param EII A \code{\link[EGAnet]{ergoInfo}} object, used to estimate the Empirical Ergodicity Information Index, or the estimated value of EII estimated
#' using the \code{\link[EGAnet]{ergoInfo}} function. Inherits \code{use} from \code{\link[EGAnet]{ergoInfo}}
#'
#' @param iter Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{100} is recommended
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#' 
#' @param progress Boolean.
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' For Windows, \code{FALSE} is about 2x faster
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
#' eii1 <- ergoInfo(dynEGA.object = dyn1, use = "weighted")
#'
#' # Bootstrap Test for Ergodicity Information Index
#' testing.ergoinfo <- boot.ergoInfo(
#'   dynEGA.object = dyn1, EII = eii1,
#'   ncores = 2
#' )}
#'
#' @return Returns a list containing:
#'
#' \item{boot.ergoInfo}{The values of the Ergodicity Information Index obtained in the bootstrap}
#'
#' \item{p.value}{The two-sided *p*-value of the bootstrap test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII
#' with small variation in the population structure}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the bootstrap distribution of EII.}
#'
#' \item{interpretation}{How you can interpret the result of the test in plain English}
#'
#' \item{plot.dist}{Histogram of the bootstrapped ergodicity information index}
#' 
#' \item{methods}{Methods to report for print/summary S3methods and automated Methods section}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @references 
#' Golino, H., Nesselroade, J., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#'
#' @export
# Bootstrap Test for the Ergodicity Information Index
# Updated 10.07.2023
boot.ergoInfo <- function(
    dynEGA.object, EII, 
    use = c("edge.list", "unweighted"),
    iter = 100, ncores, # seed = 1234,
    verbose = TRUE
){
  
  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "edge.list", ergoInfo)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Check for appropriate class ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA")
  }
  
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
    dynega_objects$individual,
    function(x){dynega_objects$population$network}
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
  effect_direction <- ifelse(
    p_value > 0.05, "n.s.",
    ifelse(mean(greater_than) > 0.50, "greater", "less")
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

#' @exportS3Method 
# S3 Print Method ----
# Updated 09.07.2023
print.boot.ergoInfo <- function(x, ...)
{
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 09.07.2023
summary.boot.ergoInfo <- function(object, ...)
{
  
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
# Updated 09.07.2023
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
  rewire_index <- sample(non_zero_edges, rewire_edges, replace = FALSE)
  
  # Get replacement indices
  replace_index <- sample(setdiff(seq_len(edges), non_zero_edges), rewire_edges, replace = FALSE)
  
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


