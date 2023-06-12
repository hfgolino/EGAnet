#' \code{\link[qgraph]{EBICglasso}} from \code{\link{qgraph}} 1.4.4
#'
#' This function uses the \code{\link[glasso]{glasso}} package
#' (Friedman, Hastie and Tibshirani, 2011) to compute a
#' sparse gaussian graphical model with the graphical lasso
#' (Friedman, Hastie & Tibshirani, 2008).
#' The tuning parameter is chosen using the Extended Bayesian Information criterium
#' (EBIC) described by Foygel & Drton (2010).
#'
#' @param data Matrix or data frame.
#' Either data (cases by variables) or correlation matrix.
#' If inputting a correlation matrix, then \code{n} must be
#' set
#'
#' @param n Numeric (length = 1).
#' Sample size for when a correlation matrix is input into \code{data}.
#' Defaults to \code{NULL}
#'
#' @param gamma EBIC tuning parameter. 0.5 is generally a good choice.
#' Setting to zero will cause regular BIC to be used.
#'
#' @param penalize.diagonal Should the diagonal be penalized?
#'
#' @param nlambda Number of lambda values to test.
#'
#' @param lambda.min.ratio Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.1}.
#' \strong{NOTE} \code{\link{qgraph}} sets the default to \code{0.01}
#'
#' @param returnAllResults Boolean.
#' Whether all results should be returned.
#' Defaults to \code{FALSE} (network only).
#' Set to \code{TRUE} to access \code{\link[glasso]{glassopath}} output
#'
#' @param penalizeMatrix Optional logical matrix to indicate which elements are penalized
#'
#' @param countDiagonal     Should diagonal be counted in EBIC computation?
#' Defaults to \code{FALSE}. Set to \code{TRUE} to mimic qgraph < 1.3 behavior (not recommended!).
#'
#' @param refit Logical, should the optimal graph be refitted without LASSO regularization?
#' Defaults to \code{FALSE}.
#' 
#' @param model.selection Character.
#' How lambda should be selected within GLASSO
#' 
#' @param verbose Boolean.
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
#' # Instantiation of GLASSO \cr
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#' Sparse inverse covariance estimation with the graphical lasso.
#' \emph{Biostatistics}, \emph{9}, 432-441.
#' 
#' # Tutorial on EBICglasso
#' Epskamp, S., & Fried, E. I. (2018).
#' A tutorial on regularized partial correlation networks.
#' \emph{Psychological Methods}, \emph{23}(4), 617–634.
#' 
#' # glasso package \cr
#' Friedman, J., Hastie, T., & Tibshirani, R. (2011).
#' glasso: Graphical lasso-estimation of Gaussian graphical models.
#' R package version 1.7.
#' 
#' # glasso + EBIC \cr
#' Foygel, R., & Drton, M. (2010).
#' Extended Bayesian information criteria for Gaussian graphical models.
#' \emph{In Advances in neural information processing systems} (pp. 604-612).
#'
#' @author Sacha Epskamp <mail@sachaepskamp.com>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' # Compute graph with tuning = 0 (BIC)
#' BICgraph <- EBICglasso.qgraph(
#'   data = wmt, gamma = 0
#' )
#'
#' # Compute graph with tuning = 0.5 (EBIC)
#' EBICgraph <- EBICglasso.qgraph(
#'   data = wmt, gamma = 0.5
#' )
#'
#' @export
#'
# Computes optimal glasso network based on EBIC ----
# Updated 12.06.2023
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
  if(missing(model.selection)){
    model.selection <- "ebic"
  }else{model.selection <- tolower(match.arg(model.selection))}
  
  # Codes originally implemented by Sacha Epskamp in his qgraph package version 1.4.4.
  # Selects optimal lambda based on EBIC for given covariance matrix.
  # EBIC is computed as in Foygel, R., & Drton, M. (2010, November). Extended Bayesian Information Criteria for Gaussian Graphical Models. In NIPS (pp. 604-612). Chicago
  
  # Check for whether data are data or correlation matrix
  if(is_symmetric(data)){
    
    # Check for sample size
    if(is.null(n)){
      stop("A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.")
    }
    
    # Set correlation matrix
    S <- data
    
  }else{ # Assume 'data' is data
    
    # Set correlation matrix
    S <- auto.correlate(data)
    
    # Obtain sample size
    n <- nrow(data)
    
  }
  
  # # Compute lambda sequence (code taken from huge package):
  # lambda.max <- max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  # lambda.min <- lambda.min.ratio * lambda.max
  # lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  
  # Simplify source for fewer computations (minimal improvement)
  S_zero_diagonal <- S - diag(nrow(S)) # makes diagonal zero
  lambda.max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
  lambda.min <- lambda.min.ratio * lambda.max
  lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  
  # Perform GLASSO path
  if(missing(penalizeMatrix)){
    glas_path <- glasso::glassopath(S, lambda, trace = 0, penalize.diagonal = penalize.diagonal, ...)
  }else{
    
    # Initialize path to be similar to `glassopath` output
    glas_path <- list(
      w = array(0, c(ncol(S), ncol(S), length(lambda))),
      wi = array(0, c(ncol(S), ncol(S), length(lambda))),
      rholist = lambda
    )
    
    # Loop over lambdas
    for (i in 1:nlambda){
      res <- glasso::glasso(S, penalizeMatrix * lambda[i], trace = 0, penalize.diagonal = penalize.diagonal, ...)
      glas_path$w[,,i] <- res$w
      glas_path$wi[,,i] <- res$wi
    }
    
  }
  
  # Determine model selection criterion
  if(model.selection == "ebic"){
    
    # Log-likelihood
    lik <- sapply(seq_along(lambda),function(i){
      logGaus(S, glas_path$wi[,,i], n)
    })
    
    # EBIC
    EBICs <- sapply(seq_along(lambda),function(i){
      EBIC(S, glas_path$wi[,,i], n, gamma, countDiagonal = countDiagonal)
    })
    
    # Optimal
    opt <- which.min(EBICs)
    
  }else if(model.selection == "jsd"){
    
    # JSD
    JSDs <- sapply(seq_along(lambda),function(i){
      
      # Try (might be error)
      res <- try(
        jsd(S, glas_path$wi[,,i]),
        silent = TRUE
      )
      
      # Check for error
      return(ifelse(is(res, "try-error"), NA, res))
      
    })
    
    # Optimal
    opt <- which.min(JSDs)
    
  }
  
  # Return network:
  net <- wi2net(glas_path$wi[,,opt])
  colnames(net) <- rownames(net) <- colnames(S)
  
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
    colnames(net) <- rownames(net) <- colnames(S)
    optwi <- glassoRes$wi
  } else {
    optwi <- glas_path$wi[,,opt]
  }
  
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
      
    }else if(model.selection == "jsd")
      
      # Add JSDs
      results$jsd <- JSDs
    
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
# According to huge??? : source comment
# Updated 10.06.2023
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
  
  return(
    n / 2 * (log(det(K)) - trace(S %*% K))
  )
}

#' @noRd
# Extended Bayesian Information Criterion
# Updated 12.06.2023
EBIC <- function(S, K, n, gamma = 0.5, E, countDiagonal = FALSE)
{
  
  # Obtain likelihood
  L <- logGaus(S, K, n)
  
  # Determine if number of edges is missing
  ## Computes edges and avoids check 
  E <- sum(K[lower.tri(K, diag = countDiagonal)] != 0)
  
  # Number of variables
  p <- ncol(K)
  
  # Return EBIC
  return(
    -2 * L + E * log(n) + 4 * E * gamma * log(p)
  )
  
}

#' @noRd
# Converts covariance to correlation matrix
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

