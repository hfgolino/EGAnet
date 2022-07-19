#' Dimension Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, this function
#' computes the stability of dimensions. This is computed by assessing the proportion of
#' times the original dimension is exactly replicated in across bootstrap samples
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of dimStability
#'
#' @return Returns a list containing:
#'
#' \item{dimension.stability}{A list containing:
#'
#' \itemize{
#'
#' \item{\strong{\code{structural.consistency}}}
#' {The proportion of times that each empirical \code{\link[EGAnet]{EGA}} dimension
#' \emph{exactly} replicates across the \code{\link[EGAnet]{bootEGA}} samples}
#'
#' \item{\strong{\code{average.item.stability}}}
#' {The average item stability in each empirical \code{\link[EGAnet]{EGA}} dimension}
#'   }
#' }
#'
#' \item{item.stability}{Results from \code{\link[EGAnet]{itemStability}}}
#'
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \donttest{# Estimate bootstrap EGA
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 100, # recommended 500
#'   plot.typicalStructure = FALSE, # No plot for CRAN checks
#'   type = "parametric", ncores = 2
#' )}
#'
#' # Estimate stability statistics
#' res <- dimensionStability(boot.wmt)
#' res$dimension.stability
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#'
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Dimension Stability function
# Updated 18.07.2022
# Revamp 27.02.2021
dimensionStability <- function(bootega.obj, ...)
{
  
  # Check for 'hierEGA' + 'bootEGA'
  if("result_lower" %in% names(bootega.obj)){
    
    # Message user for lower order analysis
    message("Performing dimension stability analysis on lower order...", appendLF = FALSE)
    
    # Set up lower and higher order for item stability function
    higher_order_EGA<- bootega.obj$result_lower$EGA$hierarchical$higher_order$EGA
    lower_order_EGA <- bootega.obj$result_lower$EGA$hierarchical$lower_order
    bootega.obj$result_lower$EGA <- lower_order_EGA
    
    # Perform dimension stability on lower order dimensions
    lower_ds <- suppressMessages(
      dimensionStability(
        bootega.obj$result_lower, IS.plot = FALSE, structure = NULL
      )
    )
    
    # Message user lower order is done
    message("done.")
    
    # Message user for higher order analysis
    message("Performing dimension stability analysis on higher order...", appendLF = FALSE)
    
    # Set up higher order for item stability function
    bootega.obj$result_higher$EGA <- higher_order_EGA
    bootega.obj$result_higher$color.palette <- bootega.obj$result_lower$color.palette
    
    # Perform item stability on higher order dimensions
    higher_ds <- suppressMessages(
      dimensionStability(
        bootega.obj$result_higher, IS.plot = FALSE, structure = NULL
      )
    )
    
    # Message user higher order is done
    message("done.")

    # Set up results
    results <- list(
      lower_order = lower_ds,
      higher_order = higher_ds
    )
    
    # Return result
    return(results)
    
  }
  
  # Check for 'bootEGA' object
  if(is(bootega.obj) != "bootEGA"){
    stop("Input for 'bootega.obj' is not a 'bootEGA' object")
  }

  # Get additional arguments
  add.args <- list(...)

  # Check if 'orig.wc' has been input as an argument
  if("orig.wc" %in% names(add.args)){

    # Give deprecation warning
    warning(
      "The 'orig.wc' argument has been deprecated in dimensionStability.\n\nInstead, the empirical EGA estimated in bootEGA's results is used"
    )
  }

  # Compute item stability
  stability.items <- itemStability(bootega.obj)

  # Compute dimension stability ----
  ## Grab empirical membership from itemStability output
  empirical.membership <- stability.items$membership$empirical

  ## Grab unique membership from itemStability output
  unique.membership <- stability.items$membership$unique

  ## Grab bootstrap membership from itemStability output
  bootstrap.membership <- stability.items$membership$bootstrap

  ## Number of dimensions
  total.dimensions <- length(unique.membership)

  ## Initialize dimension stability vector
  stability.dimensions <- numeric(total.dimensions)

  ## Name dimensions
  names(stability.dimensions) <- unique.membership

  ## Order dimensions
  stability.dimensions <- stability.dimensions[order(names(stability.dimensions))]

  # Loop through dimensions
  for(i in 1:total.dimensions){

    # Target items
    target.items <- which(paste(empirical.membership) == paste(unique.membership[i]))

    # Get dimension stability
    target.dimension.stability <- apply(bootstrap.membership, # bootstrap membership
                                        2, # across columns
                                        function(x, target.dimension){

                                          # all items equal empirical dimension
                                          all(paste(x[target.items]) == target.dimension)

                                        }, target.dimension = paste(unique.membership[i]))

    # Input mean of dimension stability (round to 3 decimal places)
    stability.dimensions[paste(unique.membership[i])] <- round(mean(target.dimension.stability, na.rm = TRUE), 3)

  }

  # Compute average item stability ----
  average.item.stability <- stability.dimensions

  # Obtain item stabilities in empirical dimensions
  empirical.dimensions <- stability.items$item.stability$empirical.dimensions

  # Loop through dimensions
  for(i in 1:total.dimensions){

    # Target items
    target.items <- which(paste(empirical.membership) == paste(unique.membership[i]))

    # Input mean of item stability (round to 3 decimal places)
    average.item.stability[paste(unique.membership[i])] <- round(mean(empirical.dimensions[target.items], na.rm = TRUE), 3)

  }

  # Initialize and output results
  results <- list()
  results$dimension.stability <- list()
  results$dimension.stability$structural.consistency <- stability.dimensions
  results$dimension.stability$average.item.stability <- average.item.stability
  results$item.stability <- stability.items

  return(results)
}
#----
