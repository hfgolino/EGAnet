#' @title Entropy Walk Penalized GLASSO
#'
#' @description The graphical least absolute shrinkage and selection operator with
#' adaptive penalization based on the random walk transition matrix of the empirical
#' partial correlation matrix
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
#' EBIC tuning parameter.
#' Defaults to \code{0} and is generally a good choice.
#' Setting to \code{0} will cause regular BIC to be used
#'
#' @param steps Numeric (length = 1).
#' Number of steps in the random walk process.
#' Defaults to \code{1}.
#' Setting steps higher will tend toward denser networks
#' (\strong{not} recommended to change)
#'
#' @param nlambda Numeric (length = 1).
#' Number of lambda values to test.
#' Defaults to \code{100}
#'
#' @param lambda.min.ratio Numeric (length = 1).
#' Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.1}.
#' \strong{NOTE} \code{qgraph} sets the default to \code{0.01}
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
#' RW_network <- network.entropyWalk(data = wmt)
#'
#' @export
#'
# Estimate based on entropy walk ----
# Updated 30.12.2024
network.entropyWalk <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "cosine", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    gamma = 0.00, steps = 1, nlambda = 100,
    lambda.min.ratio = 0.1, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network.entropyWalk)
  na.data <- set_default(na.data, "pairwise", network.entropyWalk)

  # Argument errors (return data in case of tibble)
  data <- network.entropyWalk_errors(data, n, gamma, steps, nlambda, verbose, ...)

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

  # Obtain absolute values on the inverse covariance matrix
  absolute <- abs(solve(S))

  # Solve for transition matrix and compute cross-product
  transition <- solve(diag(rowSums(absolute)), absolute)

  # Check for more than one step
  if(steps > 1){
    transition <- Reduce(`%*%`, replicate(n = steps, transition, simplify = FALSE))
  }

  # Set up transition entropy matrix (store copy)
  transition <- (transition + t(transition)) / 2
  Tentropy <- -transition * log(transition)

  # Get glasso outputs
  glasso_list <- lapply(lambda, function(value){

    # Set entropies greater than lambda to one
    Tentropy[Tentropy > value] <- 1

    # Estimate GLASSO
    glasso::glasso(s = S, rho = value * (1 - Tentropy), penalize.diagonal = FALSE, trace = 0)

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

  # EBIC (vectorized solution; ~9x faster)
  EBICs <- -2 * lik + E * log(n) + 4 * E * gamma * log(dimensions[2])

  # Optimal
  opt <- which.min(EBICs)

  # Get R
  R <- glasso_list[[opt]]$w

  # Get W
  W <- wi2net(glasso_list[[opt]]$wi)
  dimnames(R) <- dimnames(W) <- dimnames(S)

  # Return results
  return(
    list(
      network = W, K = glasso_list[[opt]]$wi, R = R,
      correlation = S, EBIC = EBICs[[opt]], gamma = gamma
    )
  )

}

#' @noRd
# Errors ----
# Updated 28.12.2024
network.entropyWalk_errors <- function(data, n, gamma, steps, nlambda, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.entropyWalk")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network.entropyWalk")
    typeof_error(n, "numeric", "network.entropyWalk")
  }

  # 'gamma' errors
  length_error(gamma, 1, "network.entropyWalk")
  typeof_error(gamma, "numeric", "network.entropyWalk")
  range_error(gamma, c(0, Inf), "network.entropyWalk")

  # 'steps' errors
  length_error(steps, 1, "network.entropyWalk")
  typeof_error(steps, "numeric", "network.entropyWalk")
  range_error(steps, c(0, Inf), "network.entropyWalk")

  # 'nlambda' errors
  length_error(nlambda, 1, "network.entropyWalk")
  typeof_error(nlambda, "numeric", "network.entropyWalk")
  range_error(nlambda, c(1, Inf), "network.entropyWalk")

  # 'verbose' errors
  length_error(verbose, 1, "network.entropyWalk")
  typeof_error(verbose, "logical", "network.entropyWalk")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}
