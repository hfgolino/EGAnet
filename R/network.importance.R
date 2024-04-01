#' @title Network (Node) Importance
#'
#' @description Computes node-wise (variable) importance by permutating each
#' variable in the dataset to determine the average effect on prediction power
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Must be raw data (no correlation matrix)
#'
#' @param network Matrix or data frame.
#' A partial correlation network associated with \code{data}
#'
#' @param ordinal.categories Numeric (length = 1).
#' \emph{Up to} the number of categories \emph{before} a variable is considered continuous.
#' Defaults to \code{7} categories before \code{8} is considered continuous
#'
#' @param iter Numeric (length = 1).
#' Number of shuffles to perform for each node.
#' Defaults to \code{1}
#'
#' @return A list containing:
#'
#' \itemize{
#'
#' \item results --- average of the average change in predictive power over all iterations
#'
#' \item full --- a list of the average change in predictive power for each iteration
#'
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # EBICglasso (default for EGA functions)
#' glasso_network <- network.estimation(data = wmt, model = "glasso")
#'
#' # Check importance
#' network.importance(data = wmt, network = glasso_network)
#'
#' @seealso \code{\link[EGAnet]{network.predictabilty}} for the underlying function
#' and \code{\link[EGAnet]{network.generalizability}} for additional generalizability functionality
#'
#' @export
#'
# Compute node-wise (variable) importance ----
# Updated 01.04.2024
network.importance <- function(data, network, ordinal.categories = 7, iter = 1)
{

  # Argument errors (return data in case of tibble)
  output <- network.importance_errors(data, network, ordinal.categories)

  # Put everything into a matrix
  original.data <- as.matrix(output$data)
  network <- as.matrix(output$network)

  # Get node names
  node_names <- dimnames(network)[[2]]

  # Get dimensions
  dimensions <- dim(original.data)
  row_sequence <- seq_len(dimensions[1])
  column_sequence <- seq_len(dimensions[2])

  # Get iteration sequence
  iter_sequence <- seq_len(iter)

  # Get shuffle
  shuffled <- lapply(
    iter_sequence, function(iteration){
      shuffle(row_sequence)
    }
  )

  # Get original result
  original_result <- network.predictability(
    network = network, original.data = original.data,
    newdata = original.data, ordinal.categories = ordinal.categories
  )

  # Get flags
  flags <- attr(original_result$results, which = "flags")

  # Check for categories
  if(any(flags$categorical)){

    # Ensure categories start at 1
    one_start_list <- ensure_one_start(original.data, flags, dimensions[1])

    # Sort out data
    original.data <- one_start_list$original.data
    categorical_factors <- one_start_list$categorical_factors

  }

  # Get original means and standard deviations
  original_means <- colMeans(original.data, na.rm = TRUE)
  original_sds <- nvapply(
    column_sequence, function(i){sd(original.data[,i], na.rm = TRUE)}
  )

  # Set up factors
  original_factors <- lapply(
    column_sequence, function(i){

      # Check for categories
      if(flags$categorical[[i]]){

        # Set factors for data
        factor( # ensures proper tabling for accuracy
          original.data[,i], levels = seq.int(1, max(original.data[,i]), 1)
        )

      }else{ # Return data otherwise
        original.data[,i]
      }

    }
  )

  # Initialize new data
  newdata <- matrix(
    nvapply(column_sequence, function(i){
      (original.data[,i] - original_means[i]) / original_sds[i]
    }, LENGTH = dimensions[1]),
    ncol = dimensions[2],
    dimnames = list(NULL, node_names)
  )

  # Create copy for shuffled data
  shuffled_data <- newdata

  # Loop over shuffle indices
  result_list <- lapply(
    iter_sequence, function(iteration){

      # Get predictability for each node
      run <- lapply(column_sequence, function(node){

        # Shuffle node
        shuffled_data[,node] <- newdata[shuffled[[iteration]], node]

        # Get predictions
        predictions <- missing_matrix_multiply(shuffled_data, original_result$betas)

        # Loop over variables
        for(i in column_sequence){

          # Check for categories
          if(flags$categorical[[i]]){

            # Assign categories to each observation
            predictions[,i] <- as.numeric(
              cut(
                x = predictions[,i], breaks = c(
                  -Inf, handle_thresholds(original_factors[[i]]), Inf
                ) # `handle_thresholds` will properly adjust for missing thresholds
              )
            )

          }

        }

        # Obtain results
        results <- setup_results(
          predictions, original.data, flags, original_result$betas,
          node_names, dimensions, column_sequence
        )

        # Get differences
        differences <- original_result$results - results$results

        # Get indices
        indices <- rowSums(differences, na.rm = TRUE) != 0

        # Check for indices
        return(
          swiftelse(
            any(indices),
            colMeans(differences[indices,], na.rm = TRUE),
            structure(c(0, 0), names = dimnames(differences)[[2]])
          )
        )

      })

      # Return combined run
      return(as.data.frame(do.call(rbind, run)))

    }
  )

  # Set up results data frame
  results <- fast.data.frame(
    colMeans(do.call(rbind, lapply(result_list, unlist)), na.rm = TRUE),
    nrow = dimensions[2], ncol = 2,
    rownames = node_names, colnames = dimnames(result_list[[1]])[[2]]
  )

  # Set results
  results <- list(
    results = results,
    full = result_list
  )

  # Attach categories to results
  attr(results$results, "flags") <- flags

  # Set class
  class(results) <- "importance"

  # Return results
  return(results)

}

# Bug Checking ----
## Basic input
# data = wmt2[,7:24]; ordinal.categories = 7
# network = network.estimation(data, method = "glasso")

#' @noRd
# Argument errors ----
# Updated 30.03.2024
network.importance_errors <- function(data, network, ordinal.categories)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.importance")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'network' errors
  object_error(network, c("matrix", "data.frame", "tibble"), "network.importance")

  # Check for tibble
  if(get_object_type(network) == "tibble"){
    network <- as.data.frame(network)
  }

  # 'ordinal.categories' errors
  length_error(ordinal.categories, 1, "network.importance")
  typeof_error(ordinal.categories, "numeric", "network.importance")
  range_error(ordinal.categories, c(2, 11), "network.importance")

  # Return usable data (in case of tibble)
  return(list(network = network, data = data))

}

#' @exportS3Method
# S3 Print Method ----
# Updated 19.02.2024
print.importance <- function(x, ...)
{

  # Get flags
  flags <- attr(x$results, "flags")

  # Set up full flags
  full_flags <- lapply(flags, any)

  # Get node names
  node_names <- row.names(x$results)

  # Set up category prints
  ## Dichotomous
  if(full_flags$dichotomous){
    cat("Dichotomous\n\n")
    print(
      t(x$results[flags$dichotomous, c("Accuracy", "Kappa"), drop = FALSE]),
      quote = FALSE, digits = 3
    )
  }
  ## Polytomous
  if(full_flags$polytomous){

    ## Check for break from dichotomous
    if(full_flags$dichotomous){
      cat("\n----\n\n")
    }
    cat("Polytomous\n\n")
    print(
      t(x$results[flags$polytomous, c("Linear Kappa", "Krippendorff's Alpha"), drop = FALSE]),
      quote = FALSE, digits = 3
    )
  }
  ## Continuous
  if(full_flags$continuous){

    ## Check for break from dichotomous
    if(full_flags$categorical){
      cat("\n----\n\n")
    }
    cat("Continuous\n\n")
    print(
      t(x$results[flags$continuous, c("R2", "RMSE"), drop = FALSE]),
      quote = FALSE, digits = 3
    )
  }

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 11.02.2024
summary.importance <- function(object, ...)
{
  print(object, ...) # same as `print`
}