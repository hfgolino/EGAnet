#' @title Cosine similarity
#'
#' @description Computes cosine similarity
#'
#' @param x Numeric vector, matrix, or data frame.
#' If \code{nrow(x) > 1}, then \code{x} will be treated as a matrix
#' to compute an \emph{n} by \emph{n} similarity matrix (\code{y} will not be used!)
#'
#' @param y Numeric vector, matrix, or data frame.
#' Only used if \code{x} is a single variable.
#' Used to compute similarity between one variable and \emph{n} other variables.
#' Defaults to \code{NULL}
#'
#' @param ...
#' Not actually used but makes it easier for general functionality
#' in the package
#'
#' @details On missing values: \code{0} will be used to replace missing values.
#' When using (matrix) multiplication, the \code{0} value cancels out the
#' product rendering the missing value as "not counting" in the sums
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain cosines
#' wmt_cosine <- cosine(wmt)
#'
#' @export
#'
# Compute cosine similarity ----
# Updated 16.03.2024
cosine <- function(x, y = NULL, ...)
{

  # Get formatted output from error checks
  output <- cosine_errors(x, y)

  # Determine dimensions
  x_dimensions <- dim(output$x)

  # Set all missing values to 0
  output$x[is.na(output$x)] <- 0

  # Assume if `x` has more than one column,
  # then perform matrix calculations
  if(x_dimensions[2] > 1){

    # Obtain dimension names
    dimnames(x)

    # Obtain cross-product
    L <- crossprod(output$x)

    # Obtain result
    result <- L / sqrt(tcrossprod(diag(L)))

    # Get dimension names
    dimension_names <- dimnames(x)[[2]]

    # Set up dimension names
    dimnames(result) <- list(dimension_names, dimension_names)

  }else{

    # Otherwise, proceed with `x` and `y`
    y_dimensions <- dim(output$y)

    # Set all missing values to 0
    output$y[is.na(output$y)] <- 0

    # Obtain result
    result <- matrix(
      nvapply(
        seq_len(y_dimensions[2]), function(i){
          single_cosine(x, output$y[,i])
        }
      ), nrow = 1
    )

  }

  # Return values
  return(result)

}


#' @noRd
# Argument errors ----
# Updated 21.09.2024
cosine_errors <- function(x, y)
{

  # 'x' errors
  object_error(x, c("vector", "matrix", "data.frame", "tibble"), "cosine")

  # Ensure matrix
  if(!is(x, "matrix")){
    x <- as.matrix(x)
  }

  # Check for NULL
  if(!is.null(y)){

    # 'y' errors
    object_error(y, c("vector", "matrix", "data.frame", "tibble"), "cosine")

    # Ensure matrix
    if(!is(y, "matrix")){
      y <- as.matrix(y)
    }

  }

  # Return data
  return(list(x = x, y = y))

}

#' @noRd
# Single cosine ----
# Updated 21.09.2024
single_cosine <- function(x, y)
{
  return(sum(x * y) / sqrt(sum(x^2) * sum(y^2)))
}
