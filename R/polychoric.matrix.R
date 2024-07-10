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
#' (rows = cases, columns = variables).
#' Data are required to be between \code{0} and \code{11}.
#' Proper adjustments should be made prior to analysis (e.g.,
#' scales from -3 to 3 in increments of 1 should be shifted
#' by added 4 to all values)
#'
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param empty.method Character (length = 1).
#' Method for empty cell correction.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"none"} --- Adds no value (\code{empty.value = "none"})
#' to the empirical joint frequency table between two variables
#'
#' \item \code{"zero"} --- Adds \code{empty.value} to the cells with zero
#' in the joint frequency table between two variables
#'
#' \item \code{"all"} --- Adds \code{empty.value} to all
#' in the joint frequency table between two variables
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
#' \item \code{"none"} --- Adds no value (\code{0}) to the empirical joint
#' frequency table between two variables
#'
#' \item \code{"point_five"} --- Adds \code{0.5} to the cells defined by \code{empty.method}
#'
#' \item \code{"one_over"} --- Adds \code{1 / n} where \emph{n} equals the number of cells
#' based on \code{empty.method}. For \code{empty.method = "zero"},
#' \emph{n} equals the number of \emph{zero} cells
#'
#' }
#'
#' @param ... Not used but made available for easier
#' argument passing
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
#' # Compute polychoric correlation matrix
#' # with pairwise missing
#' pairwise_correlations <- polychoric.matrix(
#'   wmt, na.data = "pairwise"
#' )
#'
#' # Compute polychoric correlation matrix
#' # with listwise missing
#' pairwise_correlations <- polychoric.matrix(
#'   wmt, na.data = "listwise"
#' )
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
#' \strong{Drezner-Wesolowsky bivariate normal approximation} \cr
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
# Updated 10.07.2024
polychoric.matrix <- function(
    data, na.data = c("pairwise", "listwise"),
    empty.method = c("none", "zero", "all"),
    empty.value = c("none", "point_five", "one_over"),
    ...
)
{

  # Set default arguments if missing
  na.data <- set_default(na.data, "pairwise", polychoric.matrix)
  empty.method <- set_default(empty.method, "none", polychoric.matrix)
  if(missing(empty.value)){empty.value <- "none"}

  # Check for need to check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose = TRUE)
  }

  # Ensure data is a matrix
  data <- as.matrix(data)

  # Argument errors (try to return ordinal data)
  data <- polychoric.matrix_errors(data)

  # Check for missing data
  if(na.data == "pairwise"){
    data[is.na(data)] <- 99 # "pairwise" is performed in C
  }else if(na.data == "listwise"){
    data <- na.omit(data) # no performance difference with C
  }

  # Get dimensions of data
  dimensions <- as.integer(dim(data))

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

  # Determine standard deviations
  sds <- apply(data, 2, sd, na.rm = TRUE)

  # Get zeros
  zero_sd <- sds == 0

  # Check for any zeros
  if(any(zero_sd)){

    # Set correlations to NA
    correlations[zero_sd,] <- NA
    correlations[,zero_sd] <- NA

    # Set diagonal to one
    diag(correlations) <- 1

    # Send warning
    warning(
      paste0(
        "The following variables had zero standard deviation:\n",
        paste(dimnames(data)[[2]][zero_sd], collapse = ", ")
      )
    )

  }

  # Transfer variable names
  return(transfer_names(data, correlations))

}

#' @noRd
# Argument errors ----
# Updated 04.09.2023
polychoric.matrix_errors <- function(data)
{

  # Attempt to convert all data to be ordinal
  data <- column_apply(data, continuous2categorical)

  # 'data' errors
  range_error(data, c(0, 11), "polychoric.matrix")

  # Return data
  return(data)

}

#' @noRd
# Convert continuous to categorical data ----
# If continuous data is passed but has few values,
# then it is treated as categorical
# In order to pass these non-integer values to C,
# these values need to be converted
# Updated 04.09.2023
continuous2categorical <- function(values)
{

  # Get non-NA indices
  non_NA <- !is.na(values)

  # Get non-NA values
  actual_values <- values[non_NA]

  # Return values into actual values
  values[non_NA] <- swiftelse(
    all(actual_values == as.integer(actual_values)),
    actual_values, # all are integer, so return values "as-is"
    as.numeric(factor(actual_values))
    # not all are integer, so convert to integer
  )

  # Check for conversion
  return(values)

}



