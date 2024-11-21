#' @title Bootstrap Test for the Ergodicity Information Index
#'
#' @description Tests the Ergodicity Information Index obtained in the
#' empirical sample with a distribution of EII obtained by a variant of
#' bootstrap sampling (see \strong{Details} for the procedure)
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
#' @param shuffles Numeric.
#' Number of shuffles used to compute the Kolmogorov complexity.
#' Defaults to \code{5000}
#'
#' @param iter Numeric (length = 1).
#' Number of replica samples to generate from the bootstrap analysis.
#' Defaults to \code{100} (\code{1000} for robustness)
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
#' @details In traditional bootstrap sampling, individual participants are resampled
#' with replacement from the empirical sample. This process is time consuming
#' when carried out across \emph{v} number of variables, \emph{n} number of
#' participants, \emph{t} number of time points, and \emph{i} number of iterations.
#' Instead, \code{boot.ergoInfo} uses the premise of an ergodic process to
#' establish more efficient test that works directly on the sample's networks.
#'
#' With an ergodic process, the expectation is that all individuals will have
#' a systematic relationship with the population. Destroying this relationship
#' should result in a significant loss of information. Following this conjecture,
#' \code{boot.ergoInfo} shuffles a random subset of edges that exist in the
#' \strong{population} that is \emph{equal} to the number of shared edges
#' it has with an individual. An individual's unique edges remain the same,
#' controlling for their unique information. The result is a replicate individual
#' that contains the same total number of edges as the actual individual but
#' its shared information with the population has been scrambled.
#'
#' This process is repeated over each individual to create a replicate sample
#' and is repeated for \emph{X} iterations (e.g., 100). This approach creates
#' a sampling distribution that represents the expected information between
#' the population and individuals when a random process generates the shared
#' information between them. If the shared information between the population
#' and individuals in the empirical sample is sufficiently meaningful, then
#' this process should result in significant information loss.
#'
#' How to interpret the results: the result of \code{boot.ergoInfo} is a sampling
#' distribution of EII values that would be expected if the process was random
#' (null distribution). If the empirical EII value is \emph{greater than} or
#' not significantly different from the null distribution, then the empirical
#' data can be expected to be generated from an nonergodic process and the
#' population structure is not sufficient to describe all individuals. If the
#' empirical EII value is significantly \emph{lower than} the null distribution,
#' then the empirical data can be described by the population structure -- the
#' population structure is sufficient to describe all individuals.
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
#' The null hypothesis is that the empirical Ergodicity Information index is equal to or greater than the expected value of the EII
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
#' Golino, H., Nesselroade, J. R., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
# Bootstrap Test for the Ergodicity Information Index
# Updated 27.07.2024
boot.ergoInfo <- function(
    dynEGA.object, EII,
    use = c("edge.list", "unweighted", "weighted"),
    shuffles = 5000, iter = 100, ncores, verbose = TRUE
){

  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "unweighted", ergoInfo)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}

  # Argument errors
  boot.ergoInfo_errors(dynEGA.object, shuffles, iter, ncores, verbose)

  # Check for EII
  if(missing(EII)){ # If missing, then compute it
    EII <- ergoInfo(dynEGA.object, use = use, shuffles = shuffles)$EII
  }else if(is(EII, "EII")){

    # Get attributes
    use <- attr(EII, "methods")$use
    shuffles <- attr(EII, "methods")$shuffles
    EII <- EII$EII # Save empirical EII for last

  }

  # Get proper objects (if not, send an error)
  # Function found in `ergoInfo`
  dynEGA.object <- get_dynEGA_object(dynEGA.object)

  # Only use necessary data (saves memory!)
  population_network <- dynEGA.object$population$network
  n_dimensions <- dynEGA.object$population$n.dim
  individual_networks <- lapply(dynEGA.object$individual, function(x){x$network})

  # Remove `dynEGA.object` from memory
  rm(dynEGA.object); clear_memory()

  # Set up scramble dynamic EGA object (more efficient)
  scramble_dynEGA <- list(
    dynEGA = list(
      population = list(network = population_network, n.dim = n_dimensions),
      individual = vector("list", length = length(individual_networks))
    )
  )

  # Set class
  class(scramble_dynEGA) <- "dynEGA"

  # Get scrambled EII
  EII_values <- unlist(
    parallel_process(
      iterations = iter,
      FUN = function(iteration){

        # Population individuals
        scramble_dynEGA$dynEGA$individual <- lapply(
          individual_networks, function(individual_network){
            list(
              network = network_scramble(
                population_network, individual_network
              )
            )
          }
        )

        # Return EII
        return(ergoInfo(scramble_dynEGA, use = use, shuffles = shuffles)$EII)

      },
      # Parallelization settings
      ncores = ncores, progress = verbose
    )
  )

  # Get empirical EII greater than or equal to bootstrap values
  p_value <- mean(EII >= c(EII, EII_values))

  # Set up results to return
  results <- list(
    empirical.ergoInfo = EII,
    boot.ergoInfo = EII_values,
    p.value = p_value,
    interpretation = swiftelse(
      p_value > 0.05,
      paste(
        "The empirical EII was not different from values that would be expected",
        "if the process was random, meaning the empirical data cannot be described",
        "by the population structure -- significant information is lost when",
        "collapsing across to the population structure."
      ),
      paste(
        "The empirical EII was significantly less than values that would be",
        "expected if the process was random, meaning the empirical data can be",
        "expected to be generated by an ergodic process and the population",
        "structure is sufficient to describe all individuals."
      )
    )

  )

  # Add "methods" attribute
  attr(results, "methods") <- list(use = use, shuffles = shuffles, iter = iter)

  # Set class
  class(results) <- "boot.ergoInfo"

  # Return results
  return(results)

}

# Error checking ----
# # Estimate dynamic EGA
# dyn1 <- dynEGA.ind.pop(
#   data = sim.dynEGA[sim.dynEGA$Group == 2,-26], n.embed = 5, tau = 1,
#   delta = 1, id = 25, use.derivatives = 1,
#   model = "glasso", ncores = 8, corr = "pearson"
# )
#
# # Estimate Ergodicity Information Index
# eii1 <- ergoInfo(dynEGA.object = dyn1, use = "unweighted", shuffles = 100)
#
# # Set parameters
# dynEGA.object = dyn1; EII = eii1; use = "unweighted"
# ordering = "row"; shuffles = 100; iter = 100
# ncores = 8; verbose = TRUE

#' @noRd
# Errors ----
# Updated 12.11.2023
boot.ergoInfo_errors <- function(dynEGA.object, shuffle, iter, ncores, verbose)
{

  # 'dynEGA.object' errors ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "boot.ergoInfo")
  }

  # 'shuffle' errors
  length_error(shuffle, 1, "boot.ergoInfo")
  typeof_error(shuffle, "numeric", "boot.ergoInfo")
  range_error(shuffle, c(1, Inf), "boot.ergoInfo")

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
# Updated 27.07.2024
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
      "Iterations: ", attr(x, "methods")$iter,
      "\nMean = ", round(mean(x$boot.ergoInfo, na.rm = TRUE), 4),
      " (SD = ", round(sd(x$boot.ergoInfo, na.rm = TRUE), 4), ")",
      "\np-value = ", round(x$p.value, 4),
      "\nErgodic: ", swiftelse(x$p.value < 0.05, "Yes", "No")
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


