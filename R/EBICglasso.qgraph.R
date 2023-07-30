#' \code{\link[qgraph]{EBICglasso}} from \code{\link{qgraph}} 1.4.4
#'
#' This function uses the \code{\link[glasso]{glasso}} package
#' (Friedman, Hastie and Tibshirani, 2011) to compute a
#' sparse gaussian graphical model with the graphical lasso
#' (Friedman, Hastie & Tibshirani, 2008).
#' The tuning parameter is chosen using the Extended Bayesian Information criterion
#' (EBIC) described by Foygel & Drton (2010).
#'
#' @param data Matrix or data frame.
#' Should consist only of variables that are desired to be in analysis
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} provided is a correlation matrix
#'
#' @param gamma Numeric (length = 1)
#' EBIC tuning parameter.
#' Defaults to \code{0.50} and is generally a good choice.
#' Setting to \code{0} will cause regular BIC to be used
#'
#' @param penalize.diagonal Boolean (length = 1).
#' Should the diagonal be penalized?
#' Defaults to \code{FALSE}
#'
#' @param nlambda Numeric (length = 1).
#' Number of lambda values to test.
#' Defaults to \code{100}
#'
#' @param lambda.min.ratio Numeric (length = 1).
#' Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.1}.
#' \strong{NOTE} \code{\link{qgraph}} sets the default to \code{0.01}
#'
#' @param returnAllResults Boolean (length = 1).
#' Whether all results should be returned.
#' Defaults to \code{FALSE} (network only).
#' Set to \code{TRUE} to access \code{\link[glasso]{glassopath}} output
#'
#' @param penalizeMatrix Boolean matrix.
#' Optional logical matrix to indicate which elements are penalized
#'
#' @param countDiagonal Boolean (length = 1).
#' Should diagonal be counted in EBIC computation?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to mimic \code{\link{qgraph}} < 1.3 behavior (not recommended!)
#'
#' @param refit Boolean (length = 1).
#' Should the optimal graph be refitted without LASSO regularization?
#' Defaults to \code{FALSE}
#' 
#' @param model.selection Character (length = 1).
#' How lambda should be selected within GLASSO.
#' Defaults to \code{"EBIC"}.
#' \code{"JSD"} is experimental and should not be used otherwise
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Arguments sent to \code{\link[glasso]{glasso}}
#'
#' @details The glasso is run for 100 values of the tuning parameter logarithmically
#' spaced between the maximal value of the tuning parameter at which all edges are zero,
#' lambda_max, and lambda_max/100. For each of these graphs the EBIC is computed and
#' the graph with the best EBIC is selected. The partial correlation matrix
#' is computed using \code{\link[qgraph]{wi2net}} and returned.
#'
#' @return A partial correlation matrix
#'
#' @references
#' \strong{Instantiation of GLASSO} \cr
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#' Sparse inverse covariance estimation with the graphical lasso.
#' \emph{Biostatistics}, \emph{9}, 432-441.
#' 
#' \strong{glasso + EBIC} \cr
#' Foygel, R., & Drton, M. (2010).
#' Extended Bayesian information criteria for Gaussian graphical models.
#' \emph{In Advances in neural information processing systems} (pp. 604-612).
#' 
#' \strong{glasso package} \cr
#' Friedman, J., Hastie, T., & Tibshirani, R. (2011).
#' glasso: Graphical lasso-estimation of Gaussian graphical models.
#' R package version 1.7.
#' 
#' \strong{Tutorial on EBICglasso} \cr
#' Epskamp, S., & Fried, E. I. (2018).
#' A tutorial on regularized partial correlation networks.
#' \emph{Psychological Methods}, \emph{23}(4), 617–634.
#'
#' @author Sacha Epskamp; for maintanence,
#' Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen at gmail.com>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' # Compute graph with tuning = 0 (BIC)
#' BICgraph <- EBICglasso.qgraph(data = wmt, gamma = 0)
#'
#' # Compute graph with tuning = 0.5 (EBIC)
#' EBICgraph <- EBICglasso.qgraph(data = wmt, gamma = 0.5)
#'
#' @export
#'
# Computes optimal glasso network based on EBIC ----
# Updated 30.07.2023
EBICglasso.qgraph <- function(
    data, # Sample covariance matrix
    n = NULL,
    gamma = 0.5,
    penalize.diagonal = FALSE, # Penalize diagonal?
    nlambda = 100,
    lambda.min.ratio = 0.1,
    returnAllResults = FALSE, # If true, returns a list
    penalizeMatrix, # Optional logical matrix to indicate which elements are penalized
    countDiagonal = FALSE, # Set to TRUE to get old qgraph behavior: conting diagonal elements as parameters in EBIC computation. This is not correct, but is included to replicate older analyses
    refit = FALSE, # If TRUE, network structure is taken and non-penalized version is computed.
    model.selection = c("EBIC", "JSD"),
    verbose = FALSE,
    ... # glasso arguments
)
{
  
  # Determine model selection
  model.selection <- set_default(model.selection, "ebic", EBICglasso.qgraph)
  
  # Argument errors
  EBICglasso.qgraph_errors(
    data, n, gamma, penalize.diagonal, nlambda,
    returnAllResults, countDiagonal, refit, verbose
  )
  
  # Obtain dimensions
  dimensions <- dim(data)
  
  # Codes originally implemented by Sacha Epskamp in his qgraph package version 1.4.4.
  # Selects optimal lambda based on EBIC for given covariance matrix.
  # EBIC is computed as in Foygel, R., & Drton, M. (2010, November). Extended Bayesian Information Criteria for Gaussian Graphical Models. In NIPS (pp. 604-612). Chicago
  
  # Generic function to get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n, corr = "auto", 
    na.data = "pairwise", verbose = verbose, ...
  )
  
  # Get correlations and sample size
  S <- output$correlation_matrix; n <- output$n
  
  # # Compute lambda sequence (code taken from huge package):
  # lambda.max <- max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  # lambda.min <- lambda.min.ratio * lambda.max
  # lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  
  # Simplify source for fewer computations (minimal improvement)
  S_zero_diagonal <- S - diag(dimensions[2]) # makes diagonal zero
  lambda.max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
  lambda.min <- lambda.min.ratio * lambda.max
  lambda <- exp(seq.int(log(lambda.min), log(lambda.max), length.out = nlambda))
  
  # Obtain lambda sequence
  lambda_sequence <- seq_len(nlambda)
  
  # Perform GLASSO path
  if(missing(penalizeMatrix)){
    glas_path <- glasso::glassopath(S, lambda, trace = 0, penalize.diagonal = penalize.diagonal, ...)
  }else{
    
    # Set up array dimensions
    new_array <- array(0, c(dimensions[2], dimensions[2], nlambda))
    
    # Initialize path to be similar to `glassopath` output
    glas_path <- list(
      w = new_array,
      wi = new_array,
      rholist = lambda
    )
    
    # Loop over lambdas
    for (i in lambda_sequence){
      res <- glasso::glasso(S, penalizeMatrix * lambda[i], trace = 0, penalize.diagonal = penalize.diagonal, ...)
      glas_path$w[,,i] <- res$w
      glas_path$wi[,,i] <- res$wi
    }
    
  }
  
  # Determine model selection criterion
  if(model.selection == "ebic"){
    
    # Log-likelihood
    lik <- nvapply(lambda_sequence, function(i){
      logGaus(S, glas_path$wi[,,i], n)
    })
    
    # Compute edges
    E <- nvapply(lambda_sequence, function(i){
      edge_count(glas_path$wi[,,i], dimensions[2], countDiagonal)
    })

    # EBIC (vectorized solution; ~9x faster)
    EBICs <- -2 * lik + E * log(n) + 4 * E * gamma * log(dimensions[2])
    
    # Maintained for legacy (replaced by vectorization above)
    # EBICs <- sapply(seq_along(lambda),function(i){
    #   EBIC(S, glas_path$wi[,,i], n, gamma, countDiagonal = countDiagonal)
    # })
    
    # Optimal
    opt <- which.min(EBICs)
    
  }else if(model.selection == "jsd"){
    
    # JSD
    JSDs <- nvapply(lambda_sequence,function(i){
      
      # Try (might be error)
      res <- try(
        jsd(S, glas_path$wi[,,i]),
        silent = TRUE
      )
      
      # Check for error
      return(swiftelse(is(res, "try-error"), NA, res))
      
    })
    
    # Optimal
    opt <- which.min(JSDs)
    
  }
  
  # Return network:
  net <- wi2net(glas_path$wi[,,opt])
  net <- transfer_names(S, net)
  
  # Check empty network:
  if(all(net == 0) & isTRUE(verbose)){
    message("An empty network was selected to be the best fitting network. Possibly set 'lambda.min.ratio' higher to search more sparse networks. You can also change the 'gamma' parameter to improve sensitivity (at the cost of specificity).")
  }
  
  # Check for whether to refit:
  if(isTRUE(refit)){
    if(isTRUE(verbose)){
      message("Refitting network without LASSO regularization")
    }
    glassoRes <- silent_call(glasso::glasso(S, 0, zero = which(net == 0 & upper.tri(net), arr.ind=TRUE), trace = 0, penalize.diagonal=penalize.diagonal, ...))
    net <- wi2net(glassoRes$wi)
    net <- transfer_names(S, net)
    optwi <- glassoRes$wi
  } else {
    optwi <- glas_path$wi[,,opt]
  }
  
  # Set methods in attributes
  attr(net, "methods") <- list(
    corr = "auto",
    model.selection = model.selection,
    lambda = lambda[opt], gamma = gamma,
    lambda.min.ratio = lambda.min.ratio,
    nlambda = nlambda, criterion = swiftelse(
      model.selection == "ebic", EBICs[opt], JSDs[opt]
    )
  )
  
  # Return
  if(isTRUE(returnAllResults)){
    
    # General result structure
    result <- list(
      results = glas_path, optnet = net,
      lambda = lambda, optwi = optwi
    )
    
    # Check for model selection
    if(model.selection == "ebic"){
      
      # Add EBICs and Log-likelihoods
      results$ebic <- EBICs; results$loglik <- lik;
      
    }else if(model.selection == "jsd"){
      
      # Add JSDs
      results$jsd <- JSDs
      
    }
    
    # Return results
    return(results)
    
  }else{
    return(net) # only return network
  }
  
}

# Bug Checking ----
# ## Basic input
# data = wmt2[,7:24]; n = NULL;
# gamma = 0.5; penalize.diagonal = FALSE;
# nlambda = 100; lambda.min.ratio = 0.1;
# returnAllResults = FALSE;
# countDiagonal = FALSE; refit = FALSE;
# model.selection = "ebic"

#' @noRd
# Errors ----
# Updated 30.07.2023
EBICglasso.qgraph_errors <- function(
    data, n, gamma, penalize.diagonal, nlambda,
    returnAllResults, countDiagonal, refit, verbose
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame"))
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1)
    typeof_error(n, "numeric")
  }
  
  # 'gamma' errors
  length_error(gamma, 1)
  typeof_error(gamma, "numeric")
  range_error(gamma, c(0, Inf))
  
  # 'penalize.diagonal' errors
  length_error(penalize.diagonal, 1)
  typeof_error(penalize.diagonal, "logical")
  
  # 'nlambda' errors
  length_error(nlambda, 1)
  typeof_error(nlambda, "numeric")
  range_error(nlambda, c(1, Inf))
  
  # 'returnAllResults' errors
  length_error(returnAllResults, 1)
  typeof_error(returnAllResults, "logical")
  
  # 'countDiagonal' errors
  length_error(countDiagonal, 1)
  typeof_error(countDiagonal, "logical")
  
  # 'refit' errors
  length_error(refit, 1)
  typeof_error(refit, "logical")
  
  # 'verbose' errors
  length_error(verbose, 1)
  typeof_error(verbose, "logical")
  
}

#' @noRd
# Log-likelihood ----
# According to huge??? : source comment
# Updated 03.07.2023
logGaus <- function(S, K, n)
{
  
  # Simply computes the Gaussian log likelihood given sample covariance and estimate of precision:
  
  # Original:
  # logGaus <- function(S,K,n)
  # {
  #   SK = S %*% K
  #   tr = function(A) sum(diag(A))
  #   n/2 * (log(det(K)) - tr(SK))
  # }
  
  # From source 
  
  return(n / 2 * (log(det(K)) - trace(crossprod(K, S))))
}

#' @noRd
# Extended Bayesian Information Criterion ----
# Here for legacy (vectorization applied in function)
# Updated 18.06.2023
EBIC <- function(S, K, n, p, gamma = 0.5, E, countDiagonal = FALSE)
{
  
  # Obtain likelihood
  L <- logGaus(S, K, n)
  
  # Determine if number of edges is missing
  ## Computes edges and avoids check 
  E <- sum(K[lower.tri(K, diag = countDiagonal)] != 0)
  
  # Number of nodes
  p <- ncol(K)
  
  # Return EBIC
  return(
    -2 * L + E * log(n) + 4 * E * gamma * log(p)
  )
  
}

#' @noRd
# Converts covariance to correlation matrix ----
# Updated 10.06.2023
wi2net <- function(x)
{
  # Get correlation matrix
  x <- -stats::cov2cor(x)
  
  # Set diagonal to zero
  diag(x) <- 0
  
  # Ensure matrix is symmetric
  x <- as.matrix(Matrix::forceSymmetric(x))
  
  # Return correlation matrix
  return(x)
  
}

