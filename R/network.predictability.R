#' @title Predict New Data based on Network
#'
#' @description General function to compute a network's predictive power on new data,
#' following Haslbeck and Waldorp (2018) and Williams and Rodriguez (2022)
#'
#' This implementation is different from the \code{predictability} in the \code{mgm} package
#' (Haslbeck), which is based on (regularized) regression. This implementation uses
#' the network directly, converting the partial correlations into an implied
#' precision (inverse covariance) matrix. See \strong{Details} for more information
#'
#' @param network Matrix or data frame.
#' A partial correlation network
#'
#' @param original.data Matrix or data frame.
#' Must consist only of variables to be used to estimate the \code{network}.
#' See \strong{Examples}
#'
#' @param newdata Matrix or data frame.
#' Must consist of the same variables in the same order as \code{original.data}.
#' See \strong{Examples}
#'
#' @param ordinal.categories Numeric (length = 1).
#' \emph{Up to} the number of categories \emph{before} a variable is considered continuous.
#' Defaults to \code{7} categories before \code{8} is considered continuous
#'
#' @return Returns a list containing:
#'
#' \item{predictions}{Predicted values of \code{newdata} based on the \code{network}}
#'
#' \item{betas}{Beta coefficients derived from the \code{network}}
#'
#' \item{results}{Performance metrics for each variable in \code{newdata}}
#'
#' @details This implementation of network predictability proceeds in several steps
#' with important assumptions:
#'
#' 1. Network was estimated using (partial) correlations (not regression like the
#' \code{mgm} package!)
#'
#' 2. Original data that was used to estimate the network in 1. is necessary to
#' apply the original scaling to the new data
#'
#' 3. (Linear) regression-like coefficients are obtained by reserve engineering the
#' inverse covariance matrix using the network's partial correlations (i.e.,
#' by setting the diagonal of the network to -1 and computing the inverse
#' of the opposite signed partial correlation matrix; see \code{EGAnet:::pcor2inv})
#'
#' 4. Predicted values are obtained by matrix multiplying the new data with these
#' coefficients
#'
#' 5. \strong{Dichotomous and polytomous} data are given categorical values based
#' on the \strong{original data's} thresholds and these thresholds are used to
#' convert the continuous predicted values into their corresponding categorical values
#'
#' 6. Evaluation metrics:
#'
#' \itemize{
#'
#' \item dichotomous --- Accuracy or the percent correctly predicted for the 0s and 1s
#'
#' \item polytomous --- Accuracy based on the correctly predicting the ordinal category exactly
#' (i.e., 1 = 1, 2, = 2, etc.) and a weighted accuracy such that absolute distance of the
#' predicted value from the actual value (e.g., |prediction - actual| = 1) is used
#' as the power of 0.5. This weighted approach provides an overall distance in terms of
#' accuracy where each predicted value away from the actual value is given a harsher
#' penalty (absolute difference = accuracy value): 0 = 1.000, 1 = 0.500, 2 = 0.2500,
#' 3 = 0.1250, 4 = 0.0625, etc.
#'
#' \item continuous --- R-sqaured and root mean square error
#'
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Set seed (to reproduce results)
#' set.seed(42)
#'
#' # Split data
#' training <- sample(
#'   1:nrow(wmt), round(nrow(wmt) * 0.80) # 80/20 split
#' )
#'
#' # Set splits
#' wmt_train <- wmt[training,]
#' wmt_test <- wmt[-training,]
#'
#' # EBICglasso (default for EGA functions)
#' glasso_network <- network.estimation(
#'   data = wmt_train, model = "glasso"
#' )
#'
#' # Check predictability
#' network.predictability(
#'   network = glasso_network, original.data = wmt_train,
#'   newdata = wmt_test
#' )
#'
#' @references
#' \strong{Original Implementation of Node Predictability} \cr
#' Haslbeck, J. M., & Waldorp, L. J. (2018).
#' How well do network models predict observations? On the importance of predictability in network models.
#' \emph{Behavior Research Methods}, \emph{50}(2), 853–861.
#'
#' \strong{Derivation of Regression Coefficients Used (Formula 3)} \cr
#' Williams, D. R., & Rodriguez, J. E. (2022).
#' Why overfitting is not (usually) a problem in partial correlation networks.
#' \emph{Psychological Methods}, \emph{27}(5), 822–840.
#'
#' @export
#'
# Predict new data based on network ----
# Updated 11.02.2024
network.predictability <- function(network, original.data, newdata, ordinal.categories = 7)
{

  # Check for 'newdata'
  if(missing(newdata)){newdata <- original.data}

  # Argument errors (return data in case of tibble)
  output <- network.predictability_errors(network, original.data, newdata, ordinal.categories)

  # Put everything into a matrix
  network <- as.matrix(output$network)
  original.data <- as.matrix(output$original.data)
  newdata <- as.matrix(output$newdata)

  # Get data categories
  categories <- data_categories(original.data)

  # Set flags
  flags <- list(
    dichotomous = categories == 2,
    polytomous = categories > 2 & categories <= ordinal.categories
  )

  # Set categorical/continuous flag
  flags$categorical <- flags$dichotomous | flags$polytomous

  # Get the inverse variances
  inverse_variances <- diag(pcor2inv(network))

  # Get dimensions
  dimensions <- dim(newdata)

  # Get dimension sequence
  dim_sequence <- seq_len(dimensions[2])

  # Get node names
  node_names <- colnames(network)

  # Initialize betas
  betas <- matrix(
    0, nrow = dimensions[2], ncol = dimensions[2],
    dimnames = list(node_names, node_names)
  )

  # Compute betas
  for(i in dim_sequence){
    for(j in i:dimensions[2]){

      # Populate betas
      betas[i,j] <- network[i,j] * sqrt(inverse_variances[i] / inverse_variances[j])
      betas[j,i] <- network[j,i] * sqrt(inverse_variances[j] / inverse_variances[i])

    }
  }

  # Initialize scaled new data
  newdata_scaled <- newdata

  # Obtain means and standard deviations
  original_means <- colMeans(original.data, na.rm = TRUE)
  original_sds <- nvapply(dim_sequence, function(i){sd(original.data[,i])})

  # Scale from original data
  for(i in dim_sequence){
    newdata_scaled[,i] <- (newdata[,i] - original_means[i]) / original_sds[i]
  }

  # Get predictions
  predictions <- newdata_scaled %*% betas

  # Loop over variables
  for(i in dim_sequence){

    # Check for categories
    if(flags$categorical[[i]]){

        # Get thresholds from original data
        thresholds <- obtain_thresholds(original.data[,i])

        # Assign categories to each observation
        predictions[,i] <- as.numeric(
          cut(x = predictions[,i], breaks = c(-Inf, thresholds, Inf))
        )

        # Check for lowest category
        minimum_value <- min(original.data[,i])

        # Re-adjust minimum category to 1
        if(minimum_value <= 0){
          newdata[,i] <- newdata[,i] + (abs(minimum_value) + 1)
        }

    }

  }

  # Obtain results
  results <- setup_results(
    predictions, original.data, newdata, flags, betas,
    node_names, dimensions, dim_sequence
  )

  # Attach flags to results
  attr(results$results, "flags") <- flags

  # Set class
  class(results) <- "predictability"

  # Return results
  return(results)

}

# Bug Checking ----
## Basic input
# wmt <- NetworkToolbox::neoOpen; set.seed(42); training <- sample(1:nrow(wmt), round(nrow(wmt) * 0.80))
# original.data <- wmt[training,]; newdata <- wmt[-training,]
# network <- network.estimation(original.data, model = "glasso")
# ordinal.categories = 7

#' @noRd
# Argument errors ----
# Updated 10.02.2024
network.predictability_errors <- function(network, original.data, newdata, ordinal.categories)
{

  # 'network' errors
  object_error(network, c("matrix", "data.frame", "tibble"), "network.predict")

  # Check for tibble
  if(get_object_type(network) == "tibble"){
    network <- as.data.frame(network)
  }

  # 'original.data' errors
  object_error(original.data, c("matrix", "data.frame", "tibble"), "network.predict")

  # Check for tibble
  if(get_object_type(original.data) == "tibble"){
    original.data <- as.data.frame(original.data)
  }

  # 'newdata' errors
  object_error(newdata, c("matrix", "data.frame", "tibble"), "network.predict")

  # Check for tibble
  if(get_object_type(newdata) == "tibble"){
    newdata <- as.data.frame(newdata)
  }

  # 'ordinal.categories' errors
  length_error(ordinal.categories, 1, "network.predict")
  typeof_error(ordinal.categories, "numeric", "network.predict")
  range_error(ordinal.categories, c(2, 11), "network.predict")

  # Return usable data (in case of tibble)
  return(list(network = network, original.data = original.data, newdata = newdata))

}

#' @exportS3Method
# S3 Print Method ----
# Updated 11.02.2024
print.predictability <- function(x, ...)
{

  # Get flags
  flags <- attr(x$results, "flags")

  # Add continuous flag
  flags$continuous <- !flags$categorical

  # Set up full flags
  full_flags <- lapply(flags, any)

  # Get node names
  node_names <- row.names(x$results)

  # Set up category prints
  ## Dichotomous
  if(full_flags$dichotomous){
    cat("Dichotomous\n\n")
    print(
      t(x$results[flags$dichotomous, "Accuracy", drop = FALSE]),
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
      t(x$results[flags$polytomous, c("Accuracy", "Weighted"), drop = FALSE]),
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
summary.predictability <- function(object, ...)
{
  print(object, ...) # same as `print`
}

#' @noRd
# Set up results ----
# Updated 11.02.2024
setup_results <- function(
    predictions, original.data, newdata, flags, betas,
    node_names, dimensions, dim_sequence
)
{

  # Check for all dichotomous
  if(all(flags$dichotomous)){

    # Initialize results
    results <- fast.data.frame(
      0, nrow = dimensions[2], ncol = 1,
      rownames = node_names, colnames = "Accuracy"
    )

    # Loop over and populate results
    for(i in dim_sequence){

      # Categorical measures
      results[i,] <- data_accuracy(predictions[,i], newdata[,i])[["accuracy"]]

      # Ensure predictions start as same values
      predictions[,i] <- predictions[,i] - (min(predictions[,i]) - min(original.data[,i]))

    }

  }else if(all(flags$polytomous)){ # Check for all polytomous

    # Initialize results
    results <- fast.data.frame(
      0, nrow = dimensions[2], ncol = 2,
      rownames = node_names,
      colnames = c("Accuracy", "Weighted")
    )

    # Loop over and populate results
    for(i in dim_sequence){

      # Categorical measures
      results[i,] <- data_accuracy(
        predictions[,i], newdata[,i]
      )[c("accuracy", "weighted")]

      # Ensure predictions start as same values
      predictions[,i] <- predictions[,i] - (min(predictions[,i]) - min(original.data[,i]))

    }

  }else if(all(!flags$categorical)){ # Check for all continuous


    # Initialize results
    results <- fast.data.frame(
      0, nrow = dimensions[2], ncol = 2,
      rownames = node_names,
      colnames = c("R2", "RMSE")
    )

    # Loop over and populate results
    for(i in dim_sequence){

      # Continuous measures
      results[i,] <- c(
        data_r_squared(predictions[,i], newdata[,i]),
        data_rmse(predictions[,i], newdata[,i])
      )

    }

  }else{ # Data are mixed

    # Initialize matrix
    results <- fast.data.frame(
      0, nrow = dimensions[2], ncol = 4,
      rownames = node_names,
      colnames = c("Accuracy", "Weighted", "R2", "RMSE")
    )

    # Based on categories, compute appropriate predictions
    for(i in dim_sequence){

      # Check for categories
      if(flags$dichotomous[[i]]){

        # Set result
        results[i,] <- c(
          data_accuracy(predictions[,i], newdata[,i])[["accuracy"]],
          NA, NA, NA
        )

        # Ensure categorical predictions start as same values
        predictions[,i] <- predictions[,i] - (min(predictions[,i]) - min(original.data[,i]))

      }else if(flags$polytomous[[i]]){

        # Set result
        results[i,] <- c(
          data_accuracy(predictions[,i], newdata[,i])[c("accuracy", "weighted")],
          NA, NA
        )

        # Ensure categorical predictions start as same values
        predictions[,i] <- predictions[,i] - (min(predictions[,i]) - min(original.data[,i]))

      }else{

        # Set result
        results[i,] <- c(
          NA, NA, data_r_squared(predictions[,i], newdata[,i]),
          data_rmse(predictions[,i], newdata[,i])
        )

      }

    }

  }

  # Return final results
  return(
    list(
      predictions = predictions,
      betas = betas,
      results = results
    )
  )
}
