#' @title GLASSO with Non-convex Penalties
#'
#' @description The graphical least absolute shrinkage and selection operator with
#' a non-convex regularization penalties
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
#' @param penalty Character (length = 1).
#' iPOT (inverse power of two) = 1
#' LGP (lambda-gamma power) = 1
#' POP (plus one pareto) = 10
#' SPOT (sigmoid power of two) = 3
#'
#' @param gamma Numeric (length = 1).
#' Adjusts the shape of the penalty
#'
#' @param lambda Numeric (length = 1).
#' Adjusts the initial penalty provided to the non-convex penalty function
#'
#' @param nlambda Numeric (length = 1).
#' Number of lambda values to test.
#' Defaults to \code{100}
#'
#' @param lambda.min.ratio Numeric (length = 1).
#' Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.01}
#'
#' @param penalize.diagonal Boolean (length = 1).
#' Should the diagonal be penalized?
#' Defaults to \code{FALSE}
#'
#' @param optimize.over Character (length = 1).
#' Whether optimization of lambda, gamma, both, or no hyperparamters should be performed.
#' Defaults to \code{"none"} or no optimization
#'
#' @param ic Character (length = 1).
#' What information criterion should be used for model selection?
#' Available options include:
#'
#' \itemize{
#'
#' \item \code{"AIC"} --- Akaike's information criterion: \eqn{-2L + 2E}
#'
#' \item \code{"AICc"} --- AIC corrected: \eqn{AIC + \frac{2E^2 + 2E}{n - E - 1}}
#'
#' \item \code{"BIC"} --- Bayesian information criterion: \eqn{-2L + E \cdot \log{(n)}}
#'
#' \item \code{"EBIC"} --- Extended BIC: \eqn{BIC + 4E \cdot \gamma \cdot \log{(E)}}
#'
#' }
#'
#' Term definitions:
#'
#' \itemize{
#'
#' \eqn{n} --- sample size
#'
#' \eqn{p} --- number of variables
#'
#' \eqn{E} --- edges
#'
#' \eqn{S} --- empirical correlation matrix
#'
#' \eqn{K} --- estimated inverse covariance matrix (network)
#'
#' \eqn{L = \frac{n}{2} \cdot \log \text{det} K - \sum_{i=1}^p (SK)_{ii}}
#'
#' }
#'
#' Defaults to \code{"BIC"}
#'
#' @param ebic.gamma Numeric (length = 1)
#' Value to set gamma parameter in EBIC (see above).
#' Defaults to \code{0.50}
#'
#' \emph{Only used if \code{ic = "EBIC"}}
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
#' awe_network <- network.nonconvex(data = wmt)
#'
#' @export
#'
# Apply non-convex regularization ----
# Updated 12.01.2025
network.nonconvex <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "cosine", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    penalty = c("iPOT", "LGP", "POP", "SPOT"),
    gamma = NULL, lambda = NULL, nlambda = 50, lambda.min.ratio = 0.01,
    penalize.diagonal = TRUE, optimize.over = c("none", "lambda", "both"),
    ic = c("AIC", "AICc", "BIC", "EBIC"), ebic.gamma = 0.50,
    fast = TRUE, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network.nonconvex)
  na.data <- set_default(na.data, "pairwise", network.nonconvex)
  penalty <- set_default(penalty, "spot", network.nonconvex)
  optimize.over <- set_default(optimize.over, "none", network.nonconvex)
  ic <- set_default(ic, "bic", network.nonconvex)

  # Argument errors (return data in case of tibble)
  data <- network.nonconvex_errors(
    data, n, gamma, nlambda, lambda.min.ratio,
    penalize.diagonal, ebic.gamma, fast, verbose, ...
  )

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

  # Get number of variables
  nodes <- dim(S)[2]

  # Obtain precision matrix
  K <- solve(S)

  # Obtain GLASSO function
  glasso_FUN <- swiftelse(fast, glassoFast::glassoFast, glasso::glasso)

  # Get function arguments
  glasso_ARGS <- obtain_arguments(glasso_FUN, FUN.args = list(...))

  # Supply correlation matrix
  glasso_ARGS[[1]] <- S

  # Get derivative function
  derivative_FUN <- switch(
    penalty,
    "ipot" = ipot_derivative,
    "lgp" = lgp_derivative,
    "pop" = pop_derivative,
    "spot" = spot_derivative
  )

  # Initialize lambda matrix
  lambda_matrix <- matrix(0, nrow = nodes, ncol = nodes)

  # Check for optimization
  if(optimize.over == "none"){

    # Simplify source for fewer computations (minimal improvement)
    S_zero_diagonal <- S - diag(nodes) # makes diagonal zero
    lambda.max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
    lambda.min <- lambda.min.ratio * lambda.max
    lambda <- exp(seq.int(log(lambda.min), log(lambda.max), length.out = nlambda))

    # Obtain lambda sequence
    lambda_sequence <- seq_len(nlambda)

    # Check for gamma
    if(is.null(gamma)){

      # Set defaults
      gamma <- switch(
        penalty,
        "ipot" = 5,
        "lgp" = 5,
        "pop" = 4,
        "spot" = 3
      )

    }

    # Obtain lambda matrices
    lambda_list <- lapply(lambda, function(value){

      # Obtain lambda matrix
      lambda_matrix[] <- derivative_FUN(K = K, lambda = value, gamma = gamma)

      # Check for diagonal penalization
      if(!penalize.diagonal){
        diag(lambda_matrix) <- 0
      }

      # Return lambda matrix
      return(lambda_matrix)

    })

    # Get GLASSO output
    glasso_list <- lapply(lambda_list, function(lambda_matrix){

      # Set lambda matrix
      glasso_ARGS$rho <- lambda_matrix

      # Estimate
      return(do.call(what = glasso_FUN, args = glasso_ARGS))

    })

    # Compute ICs
    ICs <- nvapply(glasso_list, function(element){
      information_crtierion(
        S = S, K = element$wi, n = n, nodes = nodes,
        ic = ic, ebic.gamma = ebic.gamma
      )
    })

    # Optimal value
    optimal <- which.min(ICs)

    # Get R
    R <- glasso_list[[optimal]]$w

    # Get W
    W <- wi2net(glasso_list[[optimal]]$wi)
    dimnames(R) <- dimnames(W) <- dimnames(S)

    # Return results
    return(
      list(
        network = W, K = glasso_list[[optimal]]$wi, R = R,
        penalty = penalty, lambda = lambda[[optimal]], gamma = gamma,
        correlation = S, criterion = ic, IC = ICs[[optimal]]
      )
    )


  }else if(optimize.over == "lambda"){ # Same range between 0 and 1

    # Check for gamma
    if(is.null(gamma)){

      # Set defaults
      gamma <- switch(
        penalty,
        "ipot" = 5,
        "lgp" = 5,
        "pop" = 4,
        "spot" = 3
      )

    }

    # Optimize for lambda
    optimized_lambda <- optimize(
      f = lambda_optimize, interval = c(0, 1),
      gamma = gamma, K = K, S = S,
      derivative_FUN = derivative_FUN,
      glasso_FUN = glasso_FUN, glasso_ARGS = glasso_ARGS,
      lambda_matrix = lambda_matrix, penalize.diagonal = penalize.diagonal,
      ic = ic, n = n, nodes = nodes, ebic.gamma = ebic.gamma
    )

    # Obtain lambda matrix
    lambda_matrix[] <- derivative_FUN(
      K = K, lambda = optimized_lambda$minimum, gamma = gamma
    )

    # Check for diagonal penalization
    if(!penalize.diagonal){
      diag(lambda_matrix) <- 0
    }

    # Set lambda matrix
    glasso_ARGS$rho <- lambda_matrix

    # Estimate output
    output <- do.call(what = glasso_FUN, args = glasso_ARGS)

    # Get R
    R <- output$w

    # Get W
    W <- wi2net(output$wi)
    dimnames(R) <- dimnames(W) <- dimnames(S)

    # Return results
    return(
      list(
        network = W, K = output$wi, R = R,
        penalty = penalty, lambda = optimized_lambda$minimum,
        gamma = gamma, criterion = ic,
        IC = optimized_lambda$objective, correlation = S
      )
    )


  }else if(optimize.over == "gamma"){

    # Set bounds for gamma
    bounds <- switch(
      penalty,
      "ipot" = c(3, 7),
      "lgp" = c(1, 10),
      "pop" = c(1, 10),
      "spot" = c(1, 5)
    )

    # Optimize for gamma
    optimized_gamma <- optimize(
      f = gamma_optimize, interval = bounds,
      lambda = swiftelse(is.null(lambda), 0.10, lambda),
      K = K, S = S, derivative_FUN = derivative_FUN,
      glasso_FUN = glasso_FUN, glasso_ARGS = glasso_ARGS,
      lambda_matrix = lambda_matrix, penalize.diagonal = penalize.diagonal,
      ic = ic, n = n, nodes = nodes, ebic.gamma = ebic.gamma
    )

    # Obtain lambda matrix
    lambda_matrix[] <- derivative_FUN(
      K = K, lambda = lambda, gamma = optimized_gamma$minimum
    )

    # Check for diagonal penalization
    if(!penalize.diagonal){
      diag(lambda_matrix) <- 0
    }

    # Set lambda matrix
    glasso_ARGS$rho <- lambda_matrix

    # Estimate output
    output <- do.call(what = glasso_FUN, args = glasso_ARGS)

    # Get R
    R <- output$w

    # Get W
    W <- wi2net(output$wi)
    dimnames(R) <- dimnames(W) <- dimnames(S)

    # Return results
    return(
      list(
        network = W, K = output$wi, R = R,
        penalty = penalty, lambda = lambda,
        gamma = optimized_gamma$minimum, criterion = ic,
        IC = optimized_gamma$objective, correlation = S
      )
    )

  }else{ # Both at this point

    # Set bounds for gamma
    bounds <- switch(
      penalty,
      "ipot" = c(3, 7),
      "lgp" = c(1, 10),
      "pop" = c(1, 10),
      "spot" = c(1, 5)
    )

    # Perform optimization
    optimized <- DEoptim::DEoptim(
      # Parameters and function
      fn = penalty_optimize,
      # Optimization arguments
      K = K, S = S,
      derivative_FUN = derivative_FUN,
      glasso_FUN = glasso_FUN,
      glasso_ARGS = glasso_ARGS,
      lambda_matrix = lambda_matrix,
      penalize.diagonal = penalize.diagonal,
      ic = ic, n = n, nodes = nodes, ebic.gamma = ebic.gamma,
      # Optimization parameters
      lower = c(0, bounds[1]), upper = c(1, bounds[2]),
      control = DEoptim::DEoptim.control(
        NP = 50, itermax = 200, parallelType = "none",
        trace = FALSE
      )
    )

    # Obtain lambda matrix
    lambda_matrix[] <- derivative_FUN(
      K = K, lambda = optimized$optim$bestmem[[1]],
      gamma = optimized$optim$bestmem[[2]]
    )

    # Check for diagonal penalization
    if(!penalize.diagonal){
      diag(lambda_matrix) <- 0
    }

    # Set lambda matrix
    glasso_ARGS$rho <- lambda_matrix

    # Estimate output
    output <- do.call(what = glasso_FUN, args = glasso_ARGS)

    # Get R
    R <- output$w

    # Get W
    W <- wi2net(output$wi)
    dimnames(R) <- dimnames(W) <- dimnames(S)

    # Return results
    return(
      list(
        network = W, K = output$wi, R = R,
        penalty = penalty,
        lambda = optimized$optim$bestmem[[1]],
        gamma = optimized$optim$bestmem[[2]],
        correlation = S, criterion = ic,
        IC = optimized$optim$bestval
      )
    )

  }

}

#' @noRd
# Errors ----
# Updated 06.01.2025
network.nonconvex_errors <- function(
    data, n, gamma, nlambda, lambda.min.ratio,
    penalize.diagonal, ebic.gamma, fast, verbose, ...
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.nonconvex")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network.nonconvex")
    typeof_error(n, "numeric", "network.nonconvex")
  }

  # 'gamma' errors
  if(!is.null(gamma)){
    length_error(gamma, 1, "network.nonconvex")
    typeof_error(gamma, "numeric", "network.nonconvex")
    range_error(gamma, c(0, Inf), "network.nonconvex")
  }

  # 'nlambda' errors
  length_error(nlambda, 1, "network.nonconvex")
  typeof_error(nlambda, "numeric", "network.nonconvex")
  range_error(nlambda, c(1, Inf), "network.nonconvex")

  # 'lambda.min.ratio' errors
  length_error(lambda.min.ratio, 1, "network.nonconvex")
  typeof_error(lambda.min.ratio, "numeric", "network.nonconvex")
  range_error(lambda.min.ratio, c(0, 1), "network.nonconvex")

  # 'penalize.diagonal' errors
  length_error(penalize.diagonal, 1, "network.nonconvex")
  typeof_error(penalize.diagonal, "logical", "network.nonconvex")

  # 'ebic.gamma' errors
  length_error(ebic.gamma, 1, "network.nonconvex")
  typeof_error(ebic.gamma, "numeric", "network.nonconvex")
  range_error(ebic.gamma, c(0, Inf), "network.nonconvex")

  # 'fast' errors
  length_error(fast, 1, "network.nonconvex")
  typeof_error(fast, "logical", "network.nonconvex")

  # 'verbose' errors
  length_error(verbose, 1, "network.nonconvex")
  typeof_error(verbose, "logical", "network.nonconvex")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

# DERIVATIVES AND PENALTIES ----

# iPOT derivative ----
# Updated 12.01.2025
ipot_derivative <- function(K, lambda, gamma = 5)
{

  # iPOT value
  ipot_value <- 2^(-gamma)

  # Return derivative
  return(swiftelse(gamma == 0, lambda, lambda * ipot_value * abs(K)^(ipot_value - 1)))

}

#' @noRd
# iPOT penalty ----
# Updated 12.01.2025
ipot_penalty <- function(K, lambda, gamma = 5)
{
  return(lambda * abs(K)^(2^(-gamma)))
}

# LGP-norm derivative ----
# Updated 12.01.2025
lgp_derivative <- function(K, lambda, gamma = 5)
{
  return(lambda^2 * abs(K)^((lambda - gamma) / gamma) / gamma)
}

#' @noRd
# LGP-norm penalty ----
# Updated 12.01.2025
lgp_penalty <- function(K, lambda, gamma = 5)
{
  return(lambda * abs(K)^(lambda / gamma))
}

# POP derivative ----
# Updated 12.01.2025
pop_derivative <- function(K, lambda, gamma = 4)
{

  # Obtain absolute of K
  K <- abs(K) + 1

  # Return lambdas
  return((lambda * gamma) / K^(gamma + 1))

}

#' @noRd
# POP penalty ----
# Updated 12.01.2025
pop_penalty <- function(K, lambda, gamma = 4)
{

  # Obtain absolute of K
  K <- abs(K) + 1

  # Return lambdas
  return(lambda * (1 - (1 / K)^gamma))

}

# SPOT derivative ----
# Updated 12.01.2025
spot_derivative <- function(K, lambda, gamma = 3)
{

  # Obtain exponent
  exponent <- exp(-abs(K) * 2^gamma)

  # Return lambdas
  return(lambda * 2^(gamma + 1) * exponent / (exponent + 1)^2)

}

#' @noRd
# SPOT penalty ----
# Updated 12.01.2025
spot_penalty <- function(K, lambda, gamma = 3)
{
  return(2 * lambda / (1 + exp(-abs(Theta) * 2^gamma)) - lambda)
}

# The TANH penalty is equivalent to SPOT such that
# 2^(gamma - 1) == SPOT
# where
# tanh(x / 2) == 2 * sigmoid(x) - 1

# # TANH derivative
# # Updated 02.01.2025
# tanh_derivative <- function(K, lambda, gamma)
# {
#
#   # Set inner component
#   K <- abs(K * gamma)
#
#   # Return lambdas
#   return((4 * lambda * gamma) / (exp(K) + exp(-K))^2)
#
# }
#
# # TANH penalty
# # Updated 02.01.2025
# tanh_penalty <- function(K, lambda, gamma)
# {
#
#   # Set inner, positive, and negative components
#   K <- abs(K)
#   pexp <- exp(K * gamma)
#   nexp <- exp(-K * gamma)
#
#   # Return lambdas
#   return(lambda * (pexp - nexp) / (pexp + nexp))
#
# }

# OPTIMIZATION FUNCTIONS ----

#' @noRd
# lambda optimization function ----
# Updated 06.01.2025
lambda_optimize <- function(
    lambda, gamma, K, S, derivative_FUN,
    glasso_FUN, glasso_ARGS,
    lambda_matrix, penalize.diagonal,
    ic, n, nodes, ebic.gamma
)
{

  # Obtain lambda matrix
  lambda_matrix[] <- derivative_FUN(K = K, lambda = lambda, gamma = gamma)

  # Check for diagonal penalization
  if(!penalize.diagonal){
    diag(lambda_matrix) <- 0
  }

  # Obtain lambda matrix
  glasso_ARGS$rho <- lambda_matrix

  # Obtain network
  network <- do.call(what = glasso_FUN, args = glasso_ARGS)$wi

  # Compute criterion
  IC <- information_crtierion(S, network, n, nodes, ic, ebic.gamma)

  # Check for NA
  return(swiftelse(is.na(IC) | is.infinite(IC), Inf, IC))

}

#' @noRd
# gamma optimization function ----
# Updated 06.01.2025
gamma_optimize <- function(
    gamma, lambda, K, S, derivative_FUN,
    glasso_FUN, glasso_ARGS,
    lambda_matrix, penalize.diagonal,
    ic, n, nodes, ebic.gamma
)
{

  # Obtain lambda matrix
  lambda_matrix[] <- derivative_FUN(K = K, lambda = lambda, gamma = gamma)

  # Check for diagonal penalization
  if(!penalize.diagonal){
    diag(lambda_matrix) <- 0
  }

  # Obtain lambda matrix
  glasso_ARGS$rho <- lambda_matrix

  # Obtain network
  network <- do.call(what = glasso_FUN, args = glasso_ARGS)$wi

  # Compute criterion
  IC <- information_crtierion(S, network, n, nodes, ic, ebic.gamma)

  # Check for NA
  return(swiftelse(is.na(IC) | is.infinite(IC), Inf, IC))

}

#' @noRd
# Penalty optimization function ----
# Updated 06.01.2025
penalty_optimize <- function(
    params, K, S, derivative_FUN,
    glasso_FUN, glasso_ARGS,
    lambda_matrix, penalize.diagonal,
    ic, n, nodes, ebic.gamma
)
{

  # Obtain lambda matrix
  lambda_matrix[] <- derivative_FUN(K = K, lambda = params[1], gamma = params[2])

  # Check for diagonal penalization
  if(!penalize.diagonal){
    diag(lambda_matrix) <- 0
  }

  # Obtain lambda matrix
  glasso_ARGS$rho <- lambda_matrix

  # Obtain network
  network <- do.call(what = glasso_FUN, args = glasso_ARGS)$wi

  # Compute criterion
  IC <- information_crtierion(S, network, n, nodes, ic, ebic.gamma)

  # Check for NA
  return(swiftelse(is.na(IC) | is.infinite(IC), Inf, IC))

}

# INFORMATION CRITERION ----

#' @noRd
# Information criterion ----
# Updated 06.01.2024
information_crtierion <- function(S, K, n, nodes, ic, ebic.gamma)
{

  # Compute Gaussian likelihood (minus two for convenience)
  L <- -2 * (n / 2) * (log(det(K)) - sum(diag(S %*% K)))

  # Get parameters (edges)
  E <- edge_count(K, nodes)

  # Return information criterion
  return(
    switch(
      ic,
      "aic" = L + 2 * E,
      "aicc" = L + 2 * E + (2 * E^2 + 2 * E) / (n - E - 1),
      "bic" = L + E * log(n),
      "ebic" = L + E * log(n) + 4 * E * ebic.gamma * log(nodes)
    )
  )

}