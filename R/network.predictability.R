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
#' \item dichotomous --- \code{"Accuracy"} or the percent correctly predicted for the 0s and 1s
#' and \code{"Kappa"} or Cohen's Kappa (see cite)
#'
#' \item polytomous --- \code{"Linear Kappa"} or linearly weighted Kappa and
#' \code{"Krippendorff's alpha"} (see cite)
#'
#' \item continuous --- R-squared (\code{"R2"}) and root mean square error (\code{"RMSE"})
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
#' \strong{Cohen's Kappa} \cr
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#' \emph{Educational and Psychological Measurement}, \emph{20}(1), 37-46.
#'
#' Cohen, J. (1968). Weighted kappa: nominal scale agreement provision for scaled disagreement or partial credit.
#' \emph{Psychological Bulletin}, \emph{70}(4), 213-220.
#'
#' \strong{Krippendorff's alpha} \cr
#' Krippendorff, K. (2013).
#' Content analysis: An introduction to its methodology (3rd ed.).
#' Thousand Oaks, CA: Sage.
#'
#' @export
#'
# Predict new data based on network ----
# Updated 26.03.2024
network.predictability <- function(network, original.data, newdata, ordinal.categories = 7)
{

  # Check for 'newdata'
  if(missing(newdata)){newdata <- original.data}

  # Argument errors (return data in case of tibble)
  output <- network.predictability_errors(
    network, original.data, newdata, ordinal.categories
  )

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

  # Combine original and new data
  combined <- rbind(original.data, newdata)

  # Get sample size for the original data
  original_n <- dim(combined)[1] - dimensions[1]

  # Get data categories
  categories <- data_categories(combined)

  # Set flags
  flags <- list(
    dichotomous = categories == 2,
    polytomous = categories > 2 & categories <= ordinal.categories
  )

  # Set categorical/continuous flag
  flags$categorical <- flags$dichotomous | flags$polytomous
  flags$continuous <- !flags$categorical

  # Check for categories
  if(any(flags$categorical)){

    # Ensure categories start at 1
    one_start_list <- ensure_one_start(combined, flags, original_n)

    # Sort out data
    original.data <- one_start_list$original.data
    newdata <- matrix(one_start_list$newdata, nrow = dimensions[1], ncol = dimensions[2])
    categorical_factors <- one_start_list$categorical_factors

  }

  # Get the inverse variances of original data
  inverse_variances <- diag(solve(auto.correlate(original.data)))
  # Although it's probably more appropriate to use the network-implied
  # invariance (co)variances, there are some network estimation methods
  # that result in non-positive definite or near singular matrices with
  # negative inverse variances

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
  predictions <- missing_matrix_multiply(newdata_scaled, betas)

  # Loop over variables
  for(i in dim_sequence){

    # Check for categories
    if(flags$categorical[[i]]){

      # Set factors for data
      factored_data <- factor( # ensures proper tabling for accuracy
        original.data[,i], levels = seq.int(1, max(original.data[,i]), 1)
      )

      # Assign categories to each observation
      predictions[,i] <- as.numeric(
        cut(
          x = predictions[,i], breaks = c(
            -Inf, handle_thresholds(factored_data), Inf
          ) # `handle_thresholds` will properly adjust for missing thresholds
        )
      )

    }

  }

  # Obtain results
  results <- setup_results(
    predictions, newdata, flags, betas,
    node_names, dimensions, dim_sequence
  )

  # Check for categorical data
  if(any(flags$categorical)){

    # Get categorical predictions
    categorical_predictions <- results$predictions[, flags$categorical, drop = FALSE]

    # Convert predicted categories back into their original sequence
    for(i in ncol_sequence(categorical_predictions)){
      categorical_predictions[,i] <- categorical_factors[[i]][categorical_predictions[,i]]
    }

    # Return to results
    results$predictions[,flags$categorical] <- categorical_predictions[,, drop = FALSE]

  }

  # Attach categories to results
  attr(results$results, "flags") <- flags

  # Set class
  class(results) <- "predictability"

  # Return results
  return(results)

}

# Bug Checking ----
## Basic input
# wmt <- NetworkToolbox::neoOpen; set.seed(42); training <- sample(1:nrow(wmt), round(nrow(wmt) * 0.80))
# original.data <- as.matrix(wmt[training,]); newdata <- as.matrix(wmt[-training,])
# network <- network.estimation(original.data, model = "glasso")
# ordinal.categories = 7
# # Missing data
# original.data[sample(1:length(original.data), 1000)] <- NA
# newdata[sample(1:length(newdata), 1000)] <- NA

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
# Updated 19.02.2024
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
summary.predictability <- function(object, ...)
{
  print(object, ...) # same as `print`
}

#' @noRd
# Ensure that categorical data start at one ----
# Updated 26.02.2024
ensure_one_start <- function(combined, flags, original_n)
{

  # Convert categories
  categorical_data <- combined[, flags$categorical, drop = FALSE]

  # Get categorical sequence
  categorical_sequence <- ncol_sequence(categorical_data)

  # Initialize factors
  categorical_factors <- lapply(
    categorical_sequence, function(i){
      return(as.numeric(levels(factor(categorical_data[,i]))))
    }
  )

  # Get minimum values
  minimum_values <- nvapply(as.data.frame(categorical_data), min)

  # Set starting values to 1 for all categorical data
  for(i in categorical_sequence){

    # Check for categorical values *not* equal to one
    if(minimum_values[i] != 1){

      # Replace data
      categorical_data[,i] <- swiftelse(
        minimum_values[i] > 1,
        categorical_data[,i] - (minimum_values[i] - (minimum_values[i] - 1)),
        categorical_data[,i] + (abs(minimum_values[i]) - (abs(minimum_values[i]) - 1))
      )

    }

  }

  # Get original data indices
  original_index <- seq_len(original_n)

  # Return results
  return(
    list(
      original.data = categorical_data[original_index,, drop = FALSE],
      newdata = categorical_data[-original_index,, drop = FALSE],
      categorical_factors = categorical_factors
    )
  )

}

#' @noRd
# Missing data matrix multiplication ----
# Updated 17.02.2024
missing_matrix_multiply <- function(X, Y)
{

  # For this function, only X is expected to have missing data

  # Check for missing
  if(anyNA(X)){

    # Determine NA indices
    indices <- is.na(X)

    # Set NAs to zero
    X[indices] <- 0

    # Perform matrix multiplication
    outcome <- X %*% Y

    # Set outcomes to NA
    outcome[indices] <- NA

    # Return outcome
    return(outcome)

  }else{ # no missing return normal
    return(X %*% Y)
  }

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
# Updated 02.03.2024
setup_results <- function(
    predictions, newdata, flags, betas,
    node_names, dimensions, dim_sequence
)
{

  # Check for all dichotomous
  if(all(flags$dichotomous)){

    # Get results
    results <- as.data.frame(
      matrix(
        nvapply(
          dim_sequence, function(i){
            binary_accuracy(predictions[,i], newdata[,i])[c("accuracy", "kappa")]
          }, LENGTH = 2
        ), ncol = 2, byrow = TRUE,
        dimnames = list(node_names, c("Accuracy", "Kappa"))
      )
    )

  }else if(all(flags$polytomous)){ # Check for all polytomous

    # Get results
    results <- as.data.frame(
      matrix(
        nvapply(
          dim_sequence, function(i){
            ordinal_accuracy(predictions[,i], newdata[,i])[c("linear_kappa", "kripp_alpha")]
          }, LENGTH = 2
        ), ncol = 2, byrow = TRUE,
        dimnames = list(node_names, c("Linear Kappa", "Krippendorff's Alpha"))
      )
    )

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
      NA, nrow = dimensions[2], ncol = 6,
      rownames = node_names,
      colnames = c(
        "Accuracy", "Kappa",
        "Linear Kappa", "Krippendorff's Alpha",
        "R2", "RMSE"
      )
    )

    # Check for dichotomous data
    if(any(flags$dichotomous)){

      # Get results
      results[flags$dichotomous, c("Accuracy", "Kappa")] <- t(
        nvapply(
          dim_sequence[flags$dichotomous], function(i){
            binary_accuracy(predictions[,i], newdata[,i])[c("accuracy", "kappa")]
          }, LENGTH = 2
        )
      )

    }

    # Check for polytomous data
    if(any(flags$polytomous)){

      # Get results
      results[flags$polytomous, c("Linear Kappa", "Krippendorff's Alpha")] <- t(
        nvapply(
          dim_sequence[flags$polytomous], function(i){
            ordinal_accuracy(predictions[,i], newdata[,i])[c("linear_kappa", "kripp_alpha")]
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
      predictions = predictions,
      betas = betas,
      results = results
    )
  )

}

#' @noRd
# Ensure that categorical data start at one ----
# Updated 26.02.2024
ensure_one_start <- function(combined, flags, original_n)
{

  # Convert categories
  categorical_data <- combined[, flags$categorical, drop = FALSE]

  # Get categorical sequence
  categorical_sequence <- ncol_sequence(categorical_data)

  # Initialize factors
  categorical_factors <- lapply(
    categorical_sequence, function(i){
      return(as.numeric(levels(factor(categorical_data[,i]))))
    }
  )

  # Get minimum values
  minimum_values <- nvapply(as.data.frame(categorical_data), min)

  # Set starting values to 1 for all categorical data
  for(i in categorical_sequence){

    # Check for categorical values *not* equal to one
    if(minimum_values[i] != 1){

      # Replace data
      categorical_data[,i] <- swiftelse(
        minimum_values[i] > 1,
        categorical_data[,i] - (minimum_values[i] - (minimum_values[i] - 1)),
        categorical_data[,i] + (abs(minimum_values[i]) - (abs(minimum_values[i]) - 1))
      )

    }

  }

  # Get original data indices
  original_index <- seq_len(original_n)

  # Return results
  return(
    list(
      original.data = categorical_data[original_index,],
      newdata = categorical_data[-original_index,],
      categorical_factors = categorical_factors
    )
  )

}