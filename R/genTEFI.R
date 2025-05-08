#' @title Generalized Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (Generalized TEFI) of a hierarchical or correlated bifactor
#' dimensionality structure (or \code{\link{hierEGA}} objects) using Von Neumman's entropy
#' when the input is a correlation matrix. Lower values suggest better fit of a structure to the data
#'
#' @param data Matrix, data frame, or \code{\link[EGAnet]{hierEGA}} object.
#' Can be raw data or correlation matrix
#'
#' @param structure List (length = levels).
#' A list containing the hierarchical structure.
#' Each list element corresponds to increasing levels (1 = first level, 2 = second level, etc.).
#' The length of the first element (first level) should be the same as the number of variables
#' in \code{data}. Each level after should either be the number of variables (\code{ncol(data)}) or
#' the maximum number of dimensions from the preceding level (\code{max(previous_dimensions)})
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{TRUE} to see all messages and warnings for every
#' function call.
#' Set to \code{FALSE} to ignore messages and warnings
#'
#' @return Returns a \code{levels + 1} columns data frame of the Generalized Total Entropy
#' Fit Index using Von Neumman's entropy (\code{VN.Entropy.Fit}) (first column) as well as
#' each individual's levels entropy
#'
#' @examples
#' # Example using network scores
#' opt.hier <- hierEGA(
#'   data = optimism, scores = "network",
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Compute the Generalized Total Entropy Fit Index
#' genTEFI(opt.hier)
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
# Updated 08.05.2025
genTEFI <- function(data, structure = NULL, verbose = TRUE)
{

  # Argument errors (return data in case of tibble)
  data <- genTEFI_errors(data, structure, verbose)

  # Perform TEFI
  return(tefi(data, structure, verbose))

}

#' @noRd
# Argument errors
# Updated 08.05.2025
genTEFI_errors <- function(data, structure, verbose)
{

  # 'data' errors
  if(any(!grepl("EGA", class(data)))){

    # Check for appropriate data
    object_error(data, c("matrix", "data.frame", "tibble"), "genTEFI")

    # Check for tibble
    if(get_object_type(data) == "tibble"){
      data <- as.data.frame(data)
    }

  }else{
    class_error(data, "hierEGA", "genTEFI")
  }

  # 'structure' errors
  if(!is.null(structure)){
    object_error(structure, "list", "genTEFI")
  }

  # 'verbose' errors
  length_error(verbose, 1, "genTEFI")
  typeof_error(verbose, "logical", "genTEFI")

  # Return data in case of tibble
  return(data)

}


