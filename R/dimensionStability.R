#' @title Dimension Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, 
#' this function computes the stability of dimensions. Stability is 
#' computed by assessing the proportion of times the 
#' original dimension is exactly replicated in across bootstrap samples
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param IS.plot Boolean (length = 1).
#' Should the plot be produced for \code{item.replication}?
#' Defaults to \code{TRUE}
#'
#' @param structure Numeric (length = number of variables).
#' A theoretical or pre-defined structure.
#' Defaults to \code{NULL} or the empirical \code{\link[EGAnet]{EGA}}
#' result in the \code{bootega.obj}
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link[EGAnet]{itemStability}}
#'
#' @return Returns a list containing:
#'
#' \item{dimension.stability}{A list containing:
#'
#' \itemize{
#'
#' \item{\code{structural.consistency} --- }
#' {The proportion of times that each empirical \code{\link[EGAnet]{EGA}} dimension
#' \emph{exactly} replicates across the \code{\link[EGAnet]{bootEGA}} samples}
#'
#' \item{\code{average.item.stability} --- }
#' {The average item stability in each empirical \code{\link[EGAnet]{EGA}} dimension}
#' 
#' }
#' 
#' }
#'
#' \item{item.stability}{Results from \code{\link[EGAnet]{itemStability}}}
#'
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate bootstrap EGA
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 500,
#'   type = "parametric", ncores = 2
#' )}
#'
#' # Estimate stability statistics
#' dimensionStability(boot.wmt)
#'
#' @references
#' \strong{Original implementation of bootEGA} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#'
#' \strong{Conceptual introduction} \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Dimension Stability function ----
# Updated 24.07.2023
dimensionStability <- function(bootega.obj, IS.plot = TRUE, structure = NULL, ...)
{
  
  # Compute item stability (all error checks occur in `itemStability`)
  item_stability <- itemStability(bootega.obj, IS.plot, structure, ...)
  
  # Check for hierarchical EGA
  if("lower_order" %in% names(item_stability)){
    
    # Set up results list
    results <- list(
      lower_order = dimensionStability_core(item_stability$lower_order),
      higher_order = dimensionStability_core(item_stability$higher_order),
      item.stability = item_stability # redundant but cheaper print
    )
    
  }else{ # Otherwise, just run core
    results <- dimensionStability_core(item_stability)
  }
  
  # Ensure class (needed for `hierEGA` S3)
  class(results) <- "dimensionStability"
  
  # Return results
  return(results)
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 24.07.2023
print.dimensionStability <- function(x, ...)
{
  
  # First, print item stability
  print(x$item.stability)
  
  # Add breakspace
  cat("\n----\n\n")
  
  # Print structural consistency
  cat("Structural Consistency:\n\n")
  
  # Then, branch for `hierEGA`
  if("lower_order" %in% names(x)){
    
    # Print level
    cat(
      styletext(
        text = styletext(
          text =  "Lower Order\n\n", 
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )
    
    # Print lower order
    print(x$lower_order$dimension.stability$structural.consistency)
    
    # Print level
    cat(
      styletext(
        text = styletext(
          text =  "\n\nHigher Order\n\n", 
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )
    
    # Print higher order
    print(x$higher_order$dimension.stability$structural.consistency)
    
  }else{
    print(x$dimension.stability$structural.consistency)
  }
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 12.07.2023
summary.dimensionStability <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# Dimension stability core ----
# Main function -- separated to handle `hierEGA`
# Updated 24.07.2023
dimensionStability_core <- function(item_stability_object)
{
  
  # Obtain structure (convert to string for NAs)
  structure <- paste(item_stability_object$membership$structure)
  
  # Get unique structure
  unique_structure <- unique(structure)
  
  # Order unique structure
  unique_structure <- unique_structure[order(unique_structure)]
  
  # Compute dimension stability
  dimension_stability <- nvapply(
    unique_structure, function(community){
      
      # Across items in community, compute stability for dimension
      return(
        mean(
          lvapply(
            as.data.frame(
              t(item_stability_object$membership$bootstrap[,structure == community])
            ), function(row){all(row == community, na.rm = TRUE)}
          ), na.rm = TRUE
        )
      )
      
    }
  )
  
  # Compute average item stability
  average_item_stability <- nvapply(
    unique_structure, function(community){
      mean(
        item_stability_object$item.stability$empirical.dimensions[structure == community],
        na.rm = TRUE
      )
    }
  )
  
  # Set ordering
  ordering <- order(as.numeric(names(dimension_stability)))
  
  # Set up results
  results <- list(
    dimension.stability = list(
      structural.consistency = dimension_stability[ordering],
      average.item.stability = average_item_stability[ordering]
    ),
    item.stability = item_stability_object
  )
  
  # Add class
  class(results) <- "dimensionStability"
  
  # Return results
  return(results)
  
}
