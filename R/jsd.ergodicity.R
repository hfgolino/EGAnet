#' Jensen-Shannon Distance Test for Ergodicity
#'
#' @description UPDATE
#'
#' @param dynEGA.object A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object. If a \code{\link[EGAnet]{dynEGA}}
#' object, then \code{level = c("individual", "population")} is required
#' 
#' @param method Character.
#' Method to compute Jensen-Shannon Distance.
#' Defaults to \code{"spectral"}.
#' Options:
#' 
#' \itemize{
#' 
#' \item{\code{"kld"}}
#' {Uses Kullback-Leibler Divergence}
#' 
#' \item{\code{"spectral"}}
#' {Uses eigenvalues of combinatiorial Laplacian matrix to compute
#' Von Neumann entropy}
#' 
#' }
#'
#' @examples
#' \donttest{
#' # Dynamic EGA individual and population structures
#' dyn1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA[,-26], n.embed = 5, tau = 1,
#'   delta = 1, id = 25, use.derivatives = 1,
#'   model = "glasso", ncores = 2, corr = "pearson"
#' )
#'
#' # JSD Ergodicity Test
#' testing.ergoinfo <- jsd.ergodicity(
#'   dynEGA.object = dyn1
#' )}
#' 
#' @return Returns a list containing:
#'
#' \item{empirical.JSD}{The empirical JSD values of the individual networks compared 
#' against the population network}
#'
#' \item{rewired.JSD}{The JSD values of the rewired population networks compared 
#' against the population network}
#' 
#' \item{t.test}{The full result of the \code{\link[stats]{t.test}} between
#' the \code{empirical.JSD} and \code{rewired.JSD}}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the bootstrap distribution of EII.}
#'
#' \item{interpretation}{How you can interpret the result of the test in plain English}
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @references 
#' Golino, H., Nesselroade, J., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#'
#' @noRd
#' 
# JSD Ergodicity Test
# Updated 01.11.2023
jsd.ergodicity <- function(dynEGA.object, method = c("kld", "spectral")){
  
  # Send experimental message (for now)
  experimental("jsd.ergodicity")
  
  # Check for missing arguments (argument, default, function)
  method <- set_default(method, "spectral", jsd.ergodicity)
  
  # Argument errors
  jsd.ergodicity_errors(dynEGA.object, method)
  
  # Get proper objects (if not, send an error)
  # Function found in `ergoInfo`
  dynega_objects <- get_dynEGA_object(dynEGA.object)
  
  # Set population network
  population_network <- dynega_objects$population$network
  
  # Get individual networks
  individual_networks <- lapply(
    dynega_objects$individual, function(x){x$network}
  )
  
  # Compute empirical JSD between individuals and population
  empirical_JSD <- comparison_spectral_JSD(
    base = population_network,
    network_list = individual_networks
  )
  
  
  # Replace individual networks with population networks
  # with similar densities
  rewired_networks <- lapply(
    dynega_objects$individual, function(x){
  
      return(
        igraph_rewire(
          network = population_network,
          prob = runif_xoshiro(1, min = 0.05, max = 0.15),
          noise = 0.05
        )
      )
    
    }
  )
  
  # Get rewired JSD between individuals and population
  rewired_JSD <- comparison_spectral_JSD(
    base = dynega_objects$population$network,
    network_list = rewired_networks
  )
  
  # Compare against empirical networks
  result <- t.test(empirical_JSD, rewired_JSD, var.equal = FALSE)
  
  # Get effect direction
  effect_direction <- swiftelse(
    result$p.value > 0.05, "n.s.",
    swiftelse(sign(result$statistic) == 1, "greater", "less")
  )
  
  # Set up results to return
  results <- list(
    empirical.JSD = empirical_JSD,
    rewired.JSD = rewired_JSD,
    t.test = result,
    effect = effect_direction,
    interpretation = switch(
      effect_direction,
      "n.s." = "The empirical JSD was not different from random variation in the population structure, meaning significant information is lost when aggregating the results into a single, population network.",
      "less" = "The empirical JSD was less than what would be expected from random variation in the population structure, meaning non-significant information is lost when aggregating the results into a single, population network.",
      "greater" = "The empirical JSD was greater than what would be expected from random variation in the population structure, meaning significant information is lost when aggregating the results into a single, population network."
    )
  )
  
  # Add "methods" attribute
  attr(results, "methods") <- list(method = method)
  
  # Set class
  class(results) <- "jsd.ergodicity"
  
  # Return results
  return(results)
  
}

#' @noRd
# Errors ----
# Updated 01.11.2023
jsd.ergodicity_errors <- function(dynEGA.object, method)
{
  
  # 'dynEGA.object' errors ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "jsd.ergodicity")
  }
  
  # 'method' errors
  length_error(method, 1, "jsd.ergodicity")
  typeof_error(method, "character", "jsd.ergodicity")
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 01.11.2023
print.jsd.ergodicity <- function(x, ...)
{
  
  # Print lower order
  cat(
    styletext(
      text = styletext(
        text =  "Empirical JSD\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print JSD method
  cat("JSD Method: ", attr(x, "methods")$method, "\n")
  
  # Print JSD value
  cat(
    paste0(
      "Mean = ", round(mean(x$empirical.JSD, na.rm = TRUE), 3), 
      " (SD = ", round(sd(x$empirical.JSD, na.rm = TRUE), 3), ")"
    )
  )
  
  # Add breakspace
  cat("\n\n")
  
  # Print higher order
  cat(
    styletext(
      text = styletext(
        text =  "Rewired JSD\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print descriptives
  cat(
    paste0(
      "Mean = ", round(mean(x$rewired.JSD, na.rm = TRUE), 3),
      " (SD = ", round(sd(x$rewired.JSD, na.rm = TRUE), 3), ")",
      "\np-value = ", round(x$t.test$p.value, 4), " (", x$effect, ")",
      "\nErgodic: ", swiftelse(x$effect == "less", "Yes", "No")
    )
  )
  
  cat("\n\n")
  
  # Print interpretation
  cat("Interpretation:\n", x$interpretation)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 01.11.2023
summary.jsd.ergodicity <- function(object, ...)
{
  print(object, ...) # same as print
}