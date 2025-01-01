#' @title Non-convex Exponential Regularization Penalty GLASSO (NEXT)
#'
#' @description The graphical least absolute shrinkage and selection operator with
#' a non-convex exponential regularization penalty
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param n Numeric (length = 1).
#' Sample size \strong{must} be provided if \code{data} provided is a correlation matrix
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
#' \item \code{"cosine"} --- Uses \code{\link[EGAnet]{cosine}} to compute cosine similarity
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
#' @param gamma Numeric (length = 1).
#' Adjusts the shape of the penalty.
#' Defaults to \code{0.25}.
#' For greater sensitivity (at the cost of specificity),
#' \code{0.10} is recommended. For greater specificity
#' (at the cost of sensitivity), \code{1} is recommended
#'
#' @param nlambda Numeric (length = 1).
#' Number of lambda values to test.
#' Defaults to \code{100}
#'
#' @param lambda.min.ratio Numeric (length = 1).
#' Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.001}
#'
#' @param fast Boolean (length = 1).
#' Whether the \code{\link[glassoFast]{glassoFast}} version should be used
#' to estimate the GLASSO.
#' Defaults to \code{TRUE}.
#'
#' The fast results \emph{may} differ by less than floating point of the original
#' GLASSO implemented by \code{\link[glasso]{glasso}} and should not impact reproducibility much (set to \code{FALSE} if concerned)
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to \code{\link[EGAnet]{auto.correlate}}
#'
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and
#' Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain network
#' RW_network <- network.next(data = wmt)
#'
#' @export
#'
# Apply NEXT regularization ----
# Updated 01.01.2025
network.next <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "cosine", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    gamma = 0.25, nlambda = 100, lambda.min.ratio = 0.001,
    fast = TRUE, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network.next)
  na.data <- set_default(na.data, "pairwise", network.next)

  # Argument errors (return data in case of tibble)
  data <- network.next_errors(data, n, gamma, nlambda, lambda.min.ratio, fast, verbose, ...)

  # Get dimensions of the data
  dimensions <- dim(data)

  # Get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n,
    corr = corr, na.data = na.data,
    verbose = verbose, needs_usable = FALSE, # skips usable data check
    ...
  )

  # Get outputs
  data <- output$data; n <- output$n
  S <- output$correlation_matrix

  # Simplify source for fewer computations (minimal improvement)
  S_zero_diagonal <- S - diag(dimensions[2]) # makes diagonal zero
  lambda.max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
  lambda.min <- lambda.min.ratio * lambda.max
  lambda <- exp(seq.int(log(lambda.min), log(lambda.max), length.out = nlambda))

  # Obtain lambda sequence
  lambda_sequence <- seq_len(nlambda)

  # Obtain precision matrix
  K <- solve(S)

  # Obtain lambda matrices
  lambda_list <- lapply(lambda, function(value){
    next_derivative(K = K, lambda = value, gamma = gamma)
  })

  # Obtain GLASSO function
  glasso_FUN <- swiftelse(fast, glassoFast::glassoFast, glasso::glasso)

  # Get function arguments
  glasso_ARGS <- obtain_arguments(
    glasso_FUN, FUN.args = list(...)
  )

  # Supply correlation matrix
  glasso_ARGS[[1]] <- S

  # Get GLASSO output
  glasso_list <- lapply(lambda_list, function(lambda_matrix){

    # Set lambda matrix
    glasso_ARGS$rho <- lambda_matrix

    # Estimate
    return(do.call(what = glasso_FUN, args = glasso_ARGS))

  })

  # Pre-compute half of n
  half_n <- n / 2

  # Log-likelihood
  lik <- nvapply(glasso_list, function(element){
    logGaus(S, element$wi, half_n)
  })

  # Compute edges
  E <- nvapply(glasso_list, function(element){
    edge_count(element$wi, dimensions[2], FALSE)
  })

  # BIC (vectorized solution; ~9x faster)
  BIC <- -2 * lik + E * log(n)

  # Optimal
  opt <- which.min(BIC)

  # Get R
  R <- glasso_list[[opt]]$w

  # Get W
  W <- wi2net(glasso_list[[opt]]$wi)
  dimnames(R) <- dimnames(W) <- dimnames(S)

  # Return results
  return(
    list(
      network = W, K = glasso_list[[opt]]$wi, R = R,
      correlation = S, BIC = BIC[[opt]], gamma = gamma
    )
  )

}

#' @noRd
# Errors ----
# Updated 01.01.2025
network.next_errors <- function(data, n, gamma, nlambda, lambda.min.ratio, fast, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.next")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network.next")
    typeof_error(n, "numeric", "network.next")
  }

  # 'gamma' errors
  length_error(gamma, 1, "network.next")
  typeof_error(gamma, "numeric", "network.next")
  range_error(gamma, c(0, Inf), "network.next")

  # 'nlambda' errors
  length_error(nlambda, 1, "network.next")
  typeof_error(nlambda, "numeric", "network.next")
  range_error(nlambda, c(1, Inf), "network.next")

  # 'lambda.min.ratio' errors
  length_error(lambda.min.ratio, 1, "network.next")
  typeof_error(lambda.min.ratio, "numeric", "network.next")
  # range_error(lambda.min.ratio, c(1e-06, 1), "network.next")

  # 'fast' errors
  length_error(fast, 1, "network.next")
  typeof_error(fast, "logical", "network.next")

  # 'verbose' errors
  length_error(verbose, 1, "network.next")
  typeof_error(verbose, "logical", "network.next")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

#' @noRd
# NEXT derivative ----
# Updated 01.01.2025
next_derivative <- function(K, lambda, gamma = 0.25){
  return(exp((-gamma * abs(K)^gamma) / lambda))
}

#' @noRd
# NEXT penalty ----
# Updated 01.01.2025
next_penalty <- function(Theta, lambda, gamma = 0.25){
  return((1 - exp((-gamma * abs(Theta)^gamma) / lambda)) * (lambda / gamma))
}