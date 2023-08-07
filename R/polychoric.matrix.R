#' @title Computes Polychoric Correlations
#' 
#' @description A fast implementation of polychoric correlations in C.
#' Uses the Beasley-Springer-Moro algorithm (Boro & Springer, 1977; Moro, 1995)
#' to estimate the inverse univariate normal CDF, the Drezner-Wesolosky 
#' approximation (Drezner & Wesolosky, 1990) to estimate the bivariate normal
#' CDF, and Brent's method (Brent, 2013) for optimization of rho
#'
#' @param data Matrix or data frame.
#' A dataset with all ordinal values
#' (rows = cases, columns = variables)
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"} --- }
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"} --- }
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param empty.method Character (length = 1).
#' Method for empty cell correction.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"none"} --- }
#' {Adds no value (\code{empty.value = "none"})
#' to the empirical joint frequency table between two variables}
#' 
#' \item{\code{"zero"} --- }
#' {Adds \code{empty.value} to the cells with zero
#' in the joint frequency table between two variables}
#' 
#' \item{\code{"all"} --- }
#' {Adds \code{empty.value} to all
#' in the joint frequency table between two variables}
#' 
#' }
#' 
#' @param empty.value Character (length = 1).
#' Value to add to the joint frequency table cells.
#' Accepts numeric values between 0 and 1 or
#' specific methods:
#' 
#' \itemize{
#' 
#' \item{\code{"none"} --- }
#' {Adds no value (\code{0}) to the empirical joint
#' frequency table between two variables}
#' 
#' \item{\code{"point_five"} --- }
#' {Adds \code{0.5} to the cells defined by \code{empty.method}}
#' 
#' \item{\code{"one_over"} --- }
#' {Adds \code{1 / n} where \emph{n} equals the number of cells
#' based on \code{empty.method}. For \code{empty.method = "zero"},
#' \emph{n} equals the number of \emph{zero} cells
#' }
#' 
#' }
#' 
#' @return Returns a polychoric correlation matrix
#' 
#' @examples
#' # Load data (ensure matrix for missing data example)
#' wmt <- as.matrix(wmt2[,7:24])
#' 
#' # Compute polychoric correlation matrix
#' correlations <- polychoric.matrix(wmt)
#' 
#' # Randomly assign missing data
#' wmt[sample(1:length(wmt), 1000)] <- NA
#' 
#' # Compute polychoric correlation matrix with pairwise method
#' na_correlations <- polychoric.matrix(wmt, na.data = "pairwise")
#' 
#' @references 
#' \strong{Beasley-Moro-Springer algorithm} \cr
#' Beasley, J. D., & Springer, S. G. (1977).
#' Algorithm AS 111: The percentage points of the normal distribution.
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, \emph{26}(1), 118-121.
#' 
#' Moro, B. (1995).
#' The full monte.
#' \emph{Risk 8 (February)}, 57-58.
#' 
#' \strong{Brent optimization} \cr
#' Brent, R. P. (2013). 
#' Algorithms for minimization without derivatives.
#' Mineola, NY: Dover Publications, Inc.
#' 
#' \strong{Drezner-Wesolosky bivariate normal approximation} \cr
#' Drezner, Z., & Wesolowsky, G. O. (1990).
#' On the computation of the bivariate normal integral.
#' \emph{Journal of Statistical Computation and Simulation}, \emph{35}(1-2), 101-107.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Compute polychoric correlation matrix
# Updated 07.08.2023
polychoric.matrix <- function(
    data, na.data = c("pairwise", "listwise"),
    empty.method = c("none", "zero", "all"),
    empty.value = c("none", "point_five", "one_over")
)
{
  
  # Set default arguments if missing
  na.data <- set_default(na.data, "pairwise", polychoric.matrix)
  empty.method <- set_default(empty.method, "none", polychoric.matrix)
  if(missing(empty.value)){empty.value <- "none"}
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Get dimensions of data
  dimensions <- as.integer(dim(data))
  
  # Check for missing data
  if(na.data == "pairwise"){
    data[is.na(data)] <- 99 # "pairwise" is performed in C
  }else if(na.data == "listwise"){
    data <- na.omit(data) # no performance difference with C
  }
  
  # Set up 'empty.method' and 'empty.value' for C
  if(empty.method == "none"){
    empty.method <- 0L # Set no value
    empty.value <- 0 # Set no value
  }else{
    
    # Set 'empty.method'
    empty.method <- swiftelse(empty.method == "zero", 1L, 2L)
    
    # Set 'empty.value'
    if(is.character(empty.value)){
      empty.value <- swiftelse(empty.value == "point_five", 0.50, 2)
    }
    
  } 
  
  # Call from C
  correlations <- matrix(
    .Call(
      "r_polychoric_correlation_matrix",
      as.integer(data),
      empty.method, empty.value,
      dimensions[1], dimensions[2],
      PACKAGE = "EGAnet"
    ), nrow = dimensions[2], ncol = dimensions[2]
  )

  # Transfer variable names
  return(transfer_names(data, correlations))
  
}
