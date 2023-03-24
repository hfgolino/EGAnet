#' Automatic correlations based on \code{\link[qgraph]{cor_auto}}
#'
#' This wrapper is based on \code{\link[qgraph]{cor_auto}}. There
#' are some minor adjustments that make this function simpler and to
#' function within EGAnet as desired. \code{NA} values are not treated
#' as categories
#'
#' @param data Matrix or data frame.
#' Should consist only of variables that are desired to be correlated
#'
#' @param ordinal.categories Numeric (length = 1).
#' Number of categories before a variable is considered continuous.
#' Defaults to \code{7} categories before \code{8} is considered continuous
#' 
#' @param forcePD Boolean (length = 1).
#' Whether positive definite matrix should be enforced.
#' Defaults to \code{TRUE}
#' 
#' @param missing Character (length = 1).
#' Corresponds to \code{\link[lavaan]{lavCor}}'s \code{missing} argument
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages should be printed.
#' Defaults to \code{FALSE}
#' 
#' @author Sacha Epskamp <mail@sachaepskamp.com> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Obtain correlations
#' wmt_corr <- auto.correlate(wmt)
#'
#' @references
#' Sacha Epskamp, Angelique O. J. Cramer, Lourens J. Waldorp, Verena D. Schmittmann, Denny Borsboom (2012). qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software, 48(4), 1-18.
#' 
#' Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/.
#' 
#' Douglas Bates and Martin Maechler (2014). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.1-3. http://CRAN.R-project.org/package=Matrix
#'
#' @export
#'
# Automatic correlations (based on {qgraph})
# Updated 24.03.2023
auto.correlate <- function(
    data, # Matrix or data frame
    ordinal.categories = 7, # consider ordinal up to 7 categories
    forcePD = TRUE, # ensure result is positive definite
    missing = "pairwise", # use available pairwise values
    verbose = FALSE # don't print messages
)
{
  
  # Convert to data frame
  data <- as.data.frame(data)
  
  # Assumes some preprocessing has taken place
  # to only select for appropriate variables
  
  # Obtain the number of categories for each variables
  # See `helpers-general.R`
  categories <- data_categories(data)
  
  # Determine categorical variables
  categorical_variables <- which(categories <= ordinal.categories)
  
  # Check if any are categorical
  if(length(categorical_variables) != 0){
    
    # Set categorical variables to be ordered
    for(ordinal in categorical_variables){
      
      # Set variables to be ordered
      data[,ordinal] <- ordered(data[,ordinal])
      
    }
    
  }
  
  ### START COMPUTING CORRELATIONS ###
  #provide needed arguments for lavcor
  if(missing == "fiml"){
    
    #fiml needs ml, TRUE and fit to have estimation in object
    lavobject <- suppressWarnings(lavaan::lavCor(data, missing = missing, se = "standard", meanstructure = TRUE, estimator = "ML", output = "fit"))
    #compute correlation matrix from covariance matrix
    CorMat <- cov2cor(lavaan::inspect(lavobject, "est")$theta)
    class(CorMat) <- "matrix"
  }else{
    #use defaults for other options
    meanstructure <- FALSE
    estimator <- "two.step"
    CorMat <- suppressWarnings(lavaan::lavCor(data, missing = missing, meanstructure = meanstructure, estimator = estimator))
    class(CorMat) <- "matrix"
  }
  
  # Check for positive definite:
  if(isTRUE(forcePD) & !all(eigen(CorMat)$values > 0))  {
    warning("Correlation matrix is not positive definite. Finding nearest positive definite matrix")
    CorMat <- as.matrix(Matrix::nearPD(CorMat, corr = TRUE, ensureSymmetry = TRUE, keepDiag = TRUE)$mat)
  }
  
  return(CorMat)

}