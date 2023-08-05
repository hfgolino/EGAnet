#' @title Generalized Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (Generalized TEFI) of a hierarchical or correlated bifactor 
#' dimensionality structure (or \code{\link{hierEGA}} objects) using Von Neumman's entropy 
#' when the input is a correlation matrix. Lower values suggest better fit of a structure to the data
#'
#' @param data Matrix, data frame, or \code{\link[EGAnet]{hierEGA}} object.
#' Can be raw data or correlation matrix
#' 
#' @param structure For high-order and correlated bifactor structures,
#' \code{structure} should be a list containing:
#' 
#' \itemize{
#' 
#' \item{\code{lower_order} --- }
#' {A vector (length = \code{ncol(data)}) representing the first-order structure 
#' (numbers or labels for each item in each first-order factor or community)}
#'
#' \item{\code{higher_order} --- }
#' {A vector (length = \code{ncol(data)} or number of \code{lower_order} communities)representing 
#' the second-order structure (numbers or labels for each item in each second-order 
#' factor or community)}
#' 
#' }
#'
#' @return Returns a three-column data frame of the Generalized Total Entropy 
#' Fit Index using Von Neumman's entropy (\code{VN.Entropy.Fit}) (first column), as well as
#' \code{Lower_Order_VN} - TEFI for the first-order factors (second column), and
#' \code{High_Order_VN}, the equivalent for the second-order factors.
#'
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
# Updated 04.08.2023
genTEFI <- function(data, structure = NULL)
{
  
  # Argument errors
  genTEFI_errors(data, structure)
  
  # Perform TEFI
  return(tefi(data, structure))

}

#' @noRd
# Argument errors
# Updated 04.08.2023
genTEFI_errors <- function(data, structure)
{
  
  # 'data' errors
  if(any(!grepl("EGA", class(data)))){
    object_error(data, c("matrix", "data.frame"))
  }else{
    class_error(data, "hierEGA")
  }
  
  # 'structure' errors
  if(!is.null(structure)){
    object_error(structure, "list")
  }
  
}


