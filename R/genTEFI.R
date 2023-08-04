#' @title Generalized Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (Generalized TEFI) of a hierarchical or correlated bifactor dimensionality structure (or \code{\link{hierEGA}} objects) using Von Neumman's entropy when the input is a correlation matrix.
#' Lower values suggest better fit of a structure to the data.
#'
#' @param data A matrix, data frame, or correlation matrix
#'
#' @param structure For high-order and correlated bifactor structures,
#' \code{structure} should be a list containing:
#' 
#' \itemize{
#' 
#' \item{\code{lower_order}}
#' {A vector representing the first-order structure (numbers or labels for each item in each first-order factor or community).}
#'
#' \item{\code{higher_order}}
#' {A vector representing the second-order structure (numbers or labels for each item in each second-order factor or community).}
#' 
#' }
#'
#' @return Returns a list containing:
#'
#' \item{VN.Entropy.Fit}{The Generalized Total Entropy Fit Index using Von Neumman's entropy}
#'
#'
#' @examples
#' 
#' # Load data
#' data <- optimism
#' 
#' # Not run: 
#' # hierEGA example
#'  opt.hier<- hierEGA(
#'  data = optimism,
#'  algorithm = "louvain")
#'  
#'  # Create a list with the lower and higher order structures:
#'  hier.structure <- vector("list")
#'  hier.structure$lower_order <- opt.hier$dim.variables$lower
#'  hier.structure$higher_order <- opt.hier$dim.variables$higher
#'  
#'  # Compute the Generalized Total Entropy Fit Index
#'  gen.tefi.opt <- tefi(opt.hier$lower_order$correlation, structure = hier.structure)
#'  
#'
#' @references
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com>, and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
# Updated 01.08.2023
genTEFI <- function(data, structure = NULL)
{
  # Check if data is an EGA object
  if(is(data, "EGA")){
    
    stop(
      paste(
        "Data is an EGA object",
        "For genTEFI, an hierEGA object should be used instead."
      ),
      call. = FALSE
    )    
  } else{
  
  # Get correlation matrix based on hierEGA
  if(is(data, "hierEGA")){

   gen.tefi <- tefi(data, structure)
    
    }else{

      # Determine if the structure is a list
      if(get_object_type(structure) == "list"){
        gen.tefi <- tefi(data, structure)
} 
      else{
        stop(
          paste(
            "Input to 'structure' was provided but did not match expected input or 'data' is
            not an hierEGA object",
            "For hierarchical structures, 'structure' should be a list with elements",
            "\"lower_order\" and \"higher_order\""
          ),
          call. = FALSE
        )
        
      }
    }
  }
}
