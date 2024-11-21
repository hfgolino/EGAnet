#' @title Estimate the Generalizability of Network
#'
#' @description General function to compute a network's predictive power on new data,
#' following Haslbeck and Waldorp (2018) and Williams and Rodriguez (2022) and using
#' generalizability methods of data splitting, \emph{k}-folds cross-validation,
#' and leave-one-out cross-validation
#'
#' Uses \code{\link[EGAnet]{network.predictability}} as the basis to then perform
#' generalizability methods over
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param method Character (length = 1).
#' Generalizability method.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"split"} --- Performs train/test data split on the data using
#' \code{number} to adjust the size of the \strong{training} split
#'
#' \item \code{"cv"} --- (default) Performs \emph{k}-folds cross-validation using
#' \code{number} to adjust the number of folds (e.g., 5 = 80/20 splits; 10 = 90/10 splits)
#'
#' \item \code{"loocv"} --- Performs leave-one-out cross-validation. Leave-one-out
#' has a tendency to \strong{overestimate} the generalizability of the model and
#' is not recommended (\emph{k}-folds cross-validation should be preferred)
#'
#' }
#'
#' @param number Numeric (length = 1).
#' Parameter to adjust the \code{method} argument. Ranges 0-1
#' for \code{method = "split"} and 1-N for \code{method = "cv"}.
#' Defaults to \code{0.80} and \code{5}, respectively
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)
#'
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
#' Arguments can be passed along to the function
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
#'
#' }
#'
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
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
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"BGGM"} --- Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{?BGGM::estimate} for more details
#'
#' \item \code{"glasso"} --- Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details
#'
#' \item \code{"TMFG"} --- Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details
#'
#' }
#'
#' @param algorithm Character or
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#'
#' \itemize{
#'
#' \item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details
#'
#' \item \code{"louvain"} --- By default, \code{"louvain"} will implement the Louvain algorithm using
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}}
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000}
#' unless specified otherwise
#'
#' \item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details
#'
#' }
#'
#' @param uni.method Character (length = 1).
#' What unidimensionality method should be used?
#' Defaults to \code{"louvain"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"expand"} --- Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation
#'
#' \item \code{"LE"} --- Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023)
#' \emph{Behavior Research Methods} simulation
#'
#' \item \code{"louvain"} --- Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Louvain solution is used; otherwise, regular EGA is used.
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}
#'
#' }
#'
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results.
#' Set for reproducible results.
#' See \href{https://r-ega.net/articles/reproducibility-prng.html}{Reproducibility and PRNG}
#' for more details on random number generation in \code{\link{EGAnet}}
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{community.unidimensional}}
#'
#' @return Returns a list containing:
#'
#' \item{node}{Node-wise metrics output from \code{\link[EGAnet]{network.predictability}}}
#'
#' \item{community}{Community-wise metrics output from \code{\link[EGAnet]{tefi}}}
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
#' # Data splitting
#' # network.generalizability(
#' #  data = wmt2[,7:24], method = "split",
#' #  number = 0.80 # 80/20 training/testing
#' # )
#'
#' # k-folds cross-validation
#' # network.generalizability(
#' #  data = wmt2[,7:24], method = "cv",
#' #  number = 5 # 5-fold cross-validation
#' # )
#'
#' \dontrun{
#' # Leave-one-out cross-validation
#' network.generalizability(
#'   data = wmt2[,7:24], method = "loocv"
#' )}
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
#' @noRd
#'
# Perform generalizability analysis ----
# Updated 26.02.2024
network.generalizability <- function(
    data,
    # generalizability arguments
    method = c("split", "cv", "loocv"), number,
    # EGA arguments
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"), model = c("BGGM", "glasso", "TMFG"),
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    seed = NULL, ...
)
{

  # Store random state (if there is one)
  store_state()

  # Get ellipse arguments
  ellipse <- list(needs_usable = FALSE, ...)

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in (keeping non-function choices for `cor_auto`)
  method <- set_default(method, "cv", network.generalizability)
  corr <- set_default(corr, "auto", EGA)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", EGA)

  # Set number
  if(missing(number)){
    number <- switch(method, "split" = 0.80, "cv" = 5, "loocv" = 0)
  }

  # Argument errors (return data in case of tibble)
  output <- network.generalizability_errors(data, method, number, seed, ...)

  # Get necessary output (ensure matrix and names for data)
  data <- ensure_dimension_names(as.matrix(output$data))
  dimensions <- output$dimensions
  row_sequence <- seq_len(dimensions[1])

  # Check for external suppression (from `invariance`)
  if(!"suppress" %in% names(ellipse) || !ellipse$suppress){
    message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results")
  }

  # Perform analyses
  if(method == "split"){

    # Prepare indices
    shuffled_indices <- shuffle(row_sequence, round(dimensions[1] * number), seed = seed)

    # Obtain results
    train_results <- EGA(
      data = data[shuffled_indices,], corr = corr, na.data = na.data,
      model = model, algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, ...
    )

    # Get predictions
    metric_summary <- network.predictability(
      network = train_results$network, original.data = data[shuffled_indices,],
      newdata = data[-shuffled_indices,], ordinal.categories = swiftelse(
        "ordinal.categories" %in% names(ellipse), ellipse$ordinal.categories, 7
      )
    )

    # Get community results
    community_summary <- c(
      TEFI = tefi(
        data = data[-shuffled_indices,], structure = train_results$wc,
        verbose = FALSE # ignore positive definite
      )$VN.Entropy.Fit
    )

  }else if(method == "cv"){

    # Prepare indices
    shuffled_indices <- shuffle(row_sequence, dimensions[1], seed = seed)

    # Set shuffled splits
    end <- floor(
      seq.int(round(dimensions[1] / number), dimensions[1], length.out = number)
    )
    start <- c(1, end[-number] - 1)

    # Initialize lists
    train_list <- test_list <- vector("list", length = number)

    # Loop over lists
    for(i in seq_len(number)){

      # Get testing indices
      testing <- start[i]:end[i]

      # Set up lists
      train_list[[i]] <- data[shuffled_indices[-testing],]
      test_list[[i]] <- data[shuffled_indices[testing],]

    }

    # Obtain results
    train_results <- lapply(
      train_list, function(train){
        EGA(
          data = train, corr = corr, na.data = na.data,
          model = model, algorithm = algorithm, uni.method = uni.method,
          plot.EGA = FALSE, ...
        )
      }
    )

    # Get sequence for cross-validation
    cv_sequence <- seq_len(number)

    # Get predictions
    prediction_results <- lapply(
      cv_sequence, function(i){

        # Get predictability
        network.predictability(
          network = train_results[[i]]$network, original.data = train_list[[i]],
          newdata = test_list[[i]], ordinal.categories = swiftelse(
            "ordinal.categories" %in% names(ellipse), ellipse$ordinal.categories, 7
          )
        )

      }
    )

    # Get community results
    community_results <- nvapply(
      cv_sequence, function(i){

        # Get TEFI
        tefi(
          data = test_list[[i]], structure = train_results[[i]]$wc,
          verbose = FALSE # ignore positive definite
        )$VN.Entropy.Fit

      }
    )

    # Get metrics
    metric_list <- lapply(prediction_results, function(x){x$results})

    # Combine metrics
    combined_matrix <- do.call(cbind, metric_list)

    # Get metric names
    metrics_names <- dimnames(combined_matrix)[[2]]

    # Get unique metrics
    unique_metrics <- unique(metrics_names)

    # Get summaries
    metric_summary <- lapply(unique_metrics, function(metric){

      # Get selected metric
      selected_metric <- combined_matrix[,metric == metrics_names]

      # Set up return
      return(
        data.frame(
          Mean = rowMeans(selected_metric, na.rm = TRUE),
          SD = apply(selected_metric, 1, sd, na.rm = TRUE),
          Median = apply(selected_metric, 1, median, na.rm = TRUE)
        )
      )

    })

    # Set names
    names(metric_summary) <- unique_metrics

    # Set up community summary
    community_summary <- fast.data.frame(
      c(
        mean(community_results, na.rm = TRUE),
        sd(community_results, na.rm = TRUE),
        median(community_results, na.rm = TRUE)
      ), nrow = 1, ncol = 3,
      rownames = "TEFI", colnames = c("Mean", "SD", "Median")
    )

  }else if(method == "loocv"){

    # Get predictions
    prediction_results <- lapply(
      row_sequence, function(train){

        # Get EGA result
        ega_train <- EGA(
          data = data[-train,], corr = corr, na.data = na.data,
          model = model, algorithm = algorithm, uni.method = uni.method,
          plot.EGA = FALSE, ...
        )

        # Get predictability
        return(
          network.predictability(
            network = ega_train$network, original.data = data[-train,],
            newdata = data[train,,drop = FALSE], ordinal.categories = swiftelse(
              "ordinal.categories" %in% names(ellipse), ellipse$ordinal.categories, 7
            )
          )
        )

      }
    )

    # Get predictions (rather than metrics)
    prediction_matrix <- do.call(
      rbind, lapply(prediction_results, function(x){x$predictions})
    )

    # Get flags
    flags <- attr(prediction_results[[1]]$results, "flags")

    # Check for categories
    if(any(flags$categorical)){

      # Set up for combined
      combined <- rbind(data, prediction_matrix)

      # Get original sample size
      original_n <- dimensions[1]

      # Ensure categories start at 1
      one_start_list <- ensure_one_start(combined, flags, original_n)

      # Sort out data
      original.data <- one_start_list$original.data
      newdata <- one_start_list$newdata

    }

    # Set up as if at the end of `network.predictability`
    metric_summary <- setup_results(
      predictions = newdata, newdata = original.data,
      betas = NULL, node_names = dimnames(data)[[2]],
      dim_sequence = seq_len(dimensions[2])
    )

    # Attach categories to results
    attr(metric_summary$results, "flags") <- flags

    # Set class
    class(metric_summary) <- "predictability"

    # Set community summary to `NULL`
    community_summary <- NULL

  }

  # Set up results
  results <- list(
    node = metric_summary,
    community = community_summary
  )

  # Set attributes
  attr(results, "methods") <- list(method = method, number = number)

  # Set class
  class(results) <- "generalizability"

  # Restore random state (if there is one)
  restore_state()

  # Return results
  return(results)

}

# Bug Checking ----
## Basic input
# data = wmt2[,7:24]; method = "cv"; number = 5
# corr = "auto"; na.data = "pairwise"; model = "glasso"
# algorithm = "walktrap"; uni.method = "louvain"
# ellipse <- list(); seed = NULL

#' @noRd
# Argument errors ----
# Updated 13.02.2024
network.generalizability_errors <- function(data, method, number, seed, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "EGA")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Get dimensions of data
  dimensions <- dim(data)

  # 'number' errors
  length_error(number, 1, "network.generalizability")
  typeof_error(number, "numeric", "network.generalizability")
  switch(
    method,
    "split" = range_error(number, c(0, 1), "network.generalizability"),
    "cv" = range_error(number, c(1, dimensions[1]), "network.generalizability")
  )

  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "network.generalizability")
    typeof_error(seed, "numeric", "network.generalizability")
    range_error(seed,  c(0, Inf), "network.generalizability")
  }

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, TRUE)
  }

  # Return data in case of tibble
  return(list(data = data, dimensions = dimensions))


}

#' @exportS3Method
# S3 Print Method ----
# Updated 18.02.2024
print.generalizability <- function(x, ...)
{

  # Get methods
  methods <- attributes(x)$methods

  # Switch based on method
  if(methods$method == "split"){

    # Print method
    cat(
      styletext(
        text = styletext(
          text = paste0(
            "Data Split (", methods$number * 100, "/", (1 - methods$number) * 100, ")\n\n"
          ),
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )

    # Print `network.predictability` results
    print(x$node)

    # Add break
    cat("\n----\n\n")

    # Print community results
    cat("Community Metrics\n\n")
    print(x$community, quote = FALSE)

  }else if(methods$method == "cv"){

    # Print method
    cat(
      styletext(
        text = styletext(
          text = paste0(
            methods$number, "-fold Cross-validation\n\n"
          ),
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )

    # Print `network.predictability` results
    cat("Node Metrics\n\n")

    # Extract means
    metric_means <- lapply(x$node, function(x){x$Mean})

    # Combine into matrix
    metric_matrix <- round(
      matrix(
        do.call(cbind, metric_means),
        nrow = length(metric_means), byrow = TRUE,
        dimnames = list(names(x$node), row.names(x$node[[1]]))
      ), 3
    )

    # Switch out NAs with ""
    metric_matrix[is.na(metric_matrix)] <- ""

    # Print node metrics
    print(metric_matrix, quote = FALSE)

    # Add break
    cat("\n----\n\n")

    # Print community results
    cat("Community Metrics\n\n")
    print(x$community, quote = FALSE)

  }else if(methods$method == "loocv"){

    # Print method
    cat(
      styletext(
        text = styletext(
          text =  "Leave-one-out Cross-validation\n\n",
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )

    # Print `network.predictability` results
    print(x$node)

  }

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 13.02.2024
summary.generalizability <- function(object, ...)
{
  print(object, ...) # same as `print`
}
