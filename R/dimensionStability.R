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
#' \dontrun{
#' # Estimate bootstrap EGA
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 500,
#'   type = "parametric", ncores = 2
#' )}
#'
#' # Estimate stability statistics
#' res <- dimensionStability(boot.wmt)
#' res$dimension.stability
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   boot.wmt,
#'   stats = "dimensionStability"
#' )}
#' 
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
# Updated 12.07.2023
dimensionStability <- function(bootega.obj, IS.plot = TRUE, structure = NULL, ...)
{
  
  # Compute item stability
  item_stability <- itemStability(bootega.obj, IS.plot, structure, ...)

  # Obtain structure (convert to string for NAs)
  structure <- paste(item_stability$membership$structure)
  
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
              t(item_stability$membership$bootstrap[,structure == community])
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
        item_stability$item.stability$empirical.dimensions[structure == community],
        na.rm = TRUE
      )
    }
  )
  
  # Set up results
  results <- list(
    dimension.stability = list(
      structural.consistency = dimension_stability,
      average.item.stability = average_item_stability
    ),
    item.stability = item_stability
  )
  
  # Add class
  class(results) <- "dimensionStability"
  
  
  # Return results
  return(results)
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 12.07.2023
print.dimensionStability <- function(x, ...)
{
  
  # First, print item stability
  print(x$item.stability)
  
  # Add breakspace
  cat("\n----\n\n")
  
  # Print structural consistency
  cat("Structural Consistency:\n\n")
  print(x$dimension.stability$structural.consistency)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 12.07.2023
summary.dimensionStability <- function(object, ...)
{
  print(object, ...) # same as print
}