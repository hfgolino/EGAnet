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
# Updated 16.02.2024
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

  # Get dimensions
  dimensions <- dim(newdata)

  # Get dimension sequence
  dim_sequence <- seq_len(dimensions[2])

  # Get node names
  node_names <- dimnames(network)[[2]]

  # Get data categories
  categories <- data_categories(rbind(original.data, newdata))

  # Get ranges
  ranges <- nvapply(
    seq_along(node_names), function(i){
      range(original.data[,i], newdata[,i], na.rm = TRUE)
    }, LENGTH = 2
  )

  # Set flags
  flags <- list(
    dichotomous = categories == 2,
    polytomous = categories > 2 & categories <= ordinal.categories
  )

  # Set categorical/continuous flag
  flags$categorical <- flags$dichotomous | flags$polytomous
  flags$continuous <- !flags$categorical

  # Get the inverse variances (use absolute for less than ideal matrices)
  inverse_variances <- abs(diag(pcor2inv(network)))

  # Get betas
  betas <- network * sqrt(outer(inverse_variances, inverse_variances, FUN = "/"))

  # Obtain means and standard deviations
  original_means <- colMeans(original.data, na.rm = TRUE)
  original_sds <- nvapply(
    dim_sequence, function(i){sd(original.data[,i], na.rm = TRUE)}
  )

  # Scale from original data
  newdata_scaled <- matrix(
    nvapply(dim_sequence, function(i){
      (newdata[,i] - original_means[i]) / original_sds[i]
    }, LENGTH = dimensions[1]),
    ncol = dimensions[2],
    dimnames = list(NULL, node_names)
  )

  # Get predictions
  predictions <- newdata_scaled %*% betas
  adjusted_predictions <- predictions

  # Loop over variables
  for(i in dim_sequence){

    # Check for categories
    if(flags$categorical[[i]]){

      # Set factors for data
      factored_data <- factor(
        original.data[,i], levels = seq.int(ranges[1,i], ranges[2,i], 1)
      )

      # Handle thresholds
      thresholds <- handle_thresholds(factored_data)

      # Assign categories to each observation
      predictions[,i] <- as.numeric(
        cut(x = predictions[,i], breaks = c(-Inf, thresholds, Inf))
      )

      # Check for lowest category
      minimum_value <- ranges[1,i]

      # Re-adjust minimum category to 1 for new data
      if(minimum_value <= 0){
        newdata[,i] <- newdata[,i] + (abs(minimum_value) + 1)
      }

      # Set adjusted predictions (for returning)
      adjusted_predictions[,i] <- predictions[,i] + (minimum_value - 1)

    }

  }

  # Obtain results
  results <- setup_results(
    predictions, adjusted_predictions, original.data, newdata,
    flags, betas, node_names, dimensions, dim_sequence
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
# Updated 12.02.2024
print.predictability <- function(x, ...)
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
# Handle thresholds ----
# Updated 14.02.2024
handle_thresholds <- function(factored_data)
{

  # Get thresholds from original data
  thresholds <- obtain_thresholds(factored_data)

  # Detect infinities (with signs)
  infinities <- is.infinite(thresholds) * sign(thresholds)

  # Check for infinities
  if(any(abs(infinities) == 1)){

    # Convert negative infinities
    negative_infinities <- infinities == -1
    thresholds[negative_infinities] <- -1000 + (1:sum(negative_infinities))

    # Convert positive infinities
    positive_infinities <- infinities == 1
    thresholds[positive_infinities] <- 1000 - (sum(positive_infinities):1)

  }

  # Determine identical thresholds
  while(any(colSums(outer(thresholds, thresholds, FUN = "==")) > 1)){

    # Loop over thresholds and increase from the highest threshold
    for(i in length(thresholds):2){

      # Check for identical thresholds (provide slight nudge)
      if(thresholds[i] == thresholds[i - 1]){
        thresholds[i] <- thresholds[i] + 1e-07
      }

    }

  }

  # Return thresholds
  return(thresholds)

}

#' @noRd
# Set up results ----
# Updated 12.02.2024
setup_results <- function(
    predictions, adjusted_predictions, original.data, newdata,
    flags, betas, node_names, dimensions, dim_sequence
)
{

  # Check for all dichotomous
  if(all(flags$dichotomous)){

    # Get results
    results <- as.data.frame(
      matrix(
        nvapply(
          dim_sequence, function(i){
            categorical_accuracy(predictions[,i], newdata[,i])[["accuracy"]]
          }
        ), dimnames = list(node_names, "Accuracy")
      )
    )

  }else if(all(flags$polytomous)){ # Check for all polytomous

    # Get results
    results <- as.data.frame(
      matrix(
        nvapply(
          dim_sequence, function(i){
            categorical_accuracy(predictions[,i], newdata[,i])[c("accuracy", "weighted")]
          }, LENGTH = 2
        ), ncol = 2, byrow = TRUE,
        dimnames = list(node_names, c("Accuracy", "Weighted"))
      )
    )

    for(i in dim_sequence){
      categorical_accuracy(predictions[,i], newdata[,i])[c("accuracy", "weighted")]
    }

  }else if(all(flags$continuous)){ # Check for all continuous


    # Get results
    results <- as.data.frame(
      matrix(
        nvapply(
          dim_sequence, function(i){
            continuous_accuracy(predictions[,i], newdata[,i])
          }, LENGTH = 2
        ), ncol = 2, byrow = TRUE,
        dimnames = list(node_names, c("R2", "RMSE"))
      )
    )

  }else{ # Data are mixed

    # Initialize matrix
    results <- fast.data.frame(
      NA, nrow = dimensions[2], ncol = 4,
      rownames = node_names,
      colnames = c("Accuracy", "Weighted", "R2", "RMSE")
    )

    # Check for dichotomous data
    if(any(flags$dichotomous)){

      # Get results
      results$Accuracy[flags$dichotomous] <- t(nvapply(
        dim_sequence[flags$dichotomous], function(i){
          categorical_accuracy(predictions[,i], newdata[,i])[["accuracy"]]
        }
      ))

    }

    # Check for polytomous data
    if(any(flags$polytomous)){

      # Get results
      results[flags$polytomous, c("Accuracy", "Weighted")] <- t(
        nvapply(
          dim_sequence[flags$polytomous], function(i){
            categorical_accuracy(predictions[,i], newdata[,i])[c("accuracy", "weighted")]
          }, LENGTH = 2
        )
      )

    }

    # Check for polytomous data
    if(any(flags$continuous)){

      # Get results
      results[flags$continuous, c("R2", "RMSE")] <- t(
        nvapply(
          dim_sequence[flags$continuous], function(i){
            continuous_accuracy(predictions[,i], newdata[,i])
          }, LENGTH = 2
        )
      )

    }

  }

  # Return final results
  return(
    list(
      predictions = adjusted_predictions,
      betas = betas,
      results = results
    )
  )

}

