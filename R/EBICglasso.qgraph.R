#'  EBICglasso from qgraph 1.4.4
#'
#' This function uses the \code{\link[glasso]{glasso}} package
#' (Friedman, Hastie and Tibshirani, 2011) to compute a
#' sparse gaussian graphical model with the graphical lasso
#' (Friedman, Hastie & Tibshirani, 2008).
#' The tuning parameter is chosen using the Extended Bayesian Information criterium
#' (EBIC) described by Foygel & Drton (2010).
#'
#' @param S A covariance or correlation matrix
#'
#' @param n Sample size used in computing \code{S}
#'
#' @param gamma EBIC tuning parameter. 0.5 is generally a good choice.
#' Setting to zero will cause regular BIC to be used.
#'
#' @param penalize.diagonal Should the diagonal be penalized?
#'
#' @param nlambda Number of lambda values to test.
#'
#' @param lambda.min.ratio Ratio of lowest lambda value compared to maximal lambda
#'
#' @param returnAllResults   If \code{TRUE} this function does not
#' return a network but the results of the entire glasso path.
#'
#' @param checkPD
#' If \code{TRUE}, the function will check if \code{S} is positive definite
#' and return an error if not. It is not advised to use a
#' non-positive definite matrix as input as (a) that can not
#' be a covariance matrix and (b) glasso can hang if the input is not positive definite.
#'
#' @param penalizeMatrix Optional logical matrix to indicate which elements are penalized
#'
#' @param countDiagonal     Should diagonal be counted in EBIC computation?
#' Defaults to \code{FALSE}. Set to \code{TRUE} to mimic qgraph < 1.3 behavior (not recommended!).
#'
#' @param refit Logical, should the optimal graph be refitted without LASSO regularization?
#' Defaults to \code{FALSE}.
#'
#' @param ... Arguments sent to \code{\link[glasso]{glasso}}
#'
#' @details The glasso is run for 100 values of the tuning parameter logarithmically
#' spaced between the maximal value of the tuning parameter at which all edges are zero,
#' lambda_max, and lambda_max/100. For each of these graphs the EBIC is computed and
#' the graph with the best EBIC is selected. The partial correlation matrix
#' is computed using \code{\link{wi2net}} and returned.
#'
#' @return A partial correlation matrix
#'
#' @references
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#' Sparse inverse covariance estimation with the graphical lasso.
#' \emph{Biostatistics}, \emph{9}, 432-441.
#' doi: \href{https://doi.org/10.1093/biostatistics/kxm045}{10.1093/biostatistics/kxm045}
#'
#' #glasso package
#' Jerome Friedman, Trevor Hastie and Rob Tibshirani (2011).
#' glasso: Graphical lasso-estimation of Gaussian graphical models.
#' R package version 1.7.
#' \url{http://CRAN.R-project.org/package=glasso}
#'
#' Foygel, R., & Drton, M. (2010).
#' Extended Bayesian information criteria for Gaussian graphical models.
#' In Advances in neural information processing systems (pp. 604-612).
#' \url{http://papers.nips.cc/paper/4087-extended-bayesian-information-criteria-for-gaussian-graphical-models}
#'
#' #psych package
#' Revelle, W. (2014) psych: Procedures for Personality and Psychological Research,
#' Northwestern University, Evanston, Illinois, USA.
#' R package version 1.4.4.
#' \url{http://CRAN.R-project.org/package=psych}
#'
#' #Matrix package
#' Douglas Bates and Martin Maechler (2014).
#' Matrix: Sparse and Dense Matrix Classes and Methods.
#' R package version 1.1-3.
#' \url{http://CRAN.R-project.org/package=Matrix}
#'
#' @author Sacha Epskamp <mail@sachaepskamp.com>
#'
#' @examples
#' \dontrun{
#' ### Using bfi dataset from psych ###
#' library("psych")
#' data(bfi)
#'
#' # Compute correlations:
#' CorMat <- cor_auto(bfi[,1:25])
#'
#' # Compute graph with tuning = 0 (BIC):
#' BICgraph <- EBICglasso.qgraph(CorMat, nrow(bfi), 0)
#'
#' # Compute graph with tuning = 0.5 (EBIC)
#' EBICgraph <- EBICglasso.qgraph(CorMat, nrow(bfi), 0.5)
#' }
#'
#' @export
#'
# Computes optimal glasso network based on EBIC:
EBICglasso.qgraph <- function(
  S, # Sample covariance matrix
  n, # Sample size
  gamma = 0.5,
  penalize.diagonal = FALSE, # Penalize diagonal?
  nlambda = 100,
  lambda.min.ratio = 0.01,
  returnAllResults = FALSE, # If true, returns a list
  checkPD = TRUE, # Checks if matrix is positive definite and stops if not
  penalizeMatrix, # Optional logical matrix to indicate which elements are penalized
  countDiagonal = FALSE, # Set to TRUE to get old qgraph behavior: conting diagonal elements as parameters in EBIC computation. This is not correct, but is included to replicate older analyses
  refit = FALSE, # If TRUE, network structure is taken and non-penalized version is computed.
  ... # glasso arguments
) {

    # Codes originally implemented by Sacha Epskamp in his qgraph package version 1.4.4.
    # Selects optimal lamba based on EBIC for given covariance matrix.
    # EBIC is computed as in Foygel, R., & Drton, M. (2010, November). Extended Bayesian Information Criteria for Gaussian Graphical Models. In NIPS (pp. 604-612). Chicago

    # Simply computes the Gaussian log likelihood given sample covariance and estimate of precision:

    # Original:
    # logGaus <- function(S,K,n)
    # {
    #   SK = S %*% K
    #   tr = function(A) sum(diag(A))
    #   n/2 * (log(det(K)) - tr(SK))
    # }

    ## According to huge???
    logGaus <- function(S,K,n)
    {
        KS = K %*% S
        tr = function(A) sum(diag(A))
        return(n/2 * (log(det(K)) - tr(KS))  )
    }

    # Computes the EBIC:
    EBIC <- function(S,K,n,gamma = 0.5,E,countDiagonal=FALSE)
    {
        #   browser()
        L <- logGaus(S, K, n)
        if (missing(E)){
            E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
        }
        p <- nrow(K)

        # return EBIC:
        -2 * L + E * log(n) + 4 * E * gamma * log(p)
    }

    # Computes partial correlation matrix given precision matrix:
    wi2net <- function(x)
    {
        x <- -stats::cov2cor(x)
        diag(x) <- 0
        x <- Matrix::forceSymmetric(x)
        return(x)
    }

  if (checkPD){
    if (any(eigen(S)$values < 0)) stop("'S' is not positive definite")
  }

  # Standardize cov matrix:
  S <- stats::cov2cor(S)

  # Compute lambda sequence (code taken from huge package):
  lambda.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  lambda.min = lambda.min.ratio*lambda.max
  lambda = exp(seq(log(lambda.min), log(lambda.max), length = nlambda))

  # Run glasso path:
  if (missing(penalizeMatrix)){
    glas_path <- glasso::glassopath(S, lambda, trace = 0, penalize.diagonal=penalize.diagonal, ...)
  }else{
    glas_path <- list(
      w = array(0, c(ncol(S), ncol(S), length(lambda))),
      wi = array(0, c(ncol(S), ncol(S), length(lambda))),
      rholist = lambda
    )

    for (i in 1:nlambda){
      res <- glasso::glasso(S, penalizeMatrix * lambda[i], trace = 0, penalize.diagonal=penalize.diagonal, ...)
      glas_path$w[,,i] <- res$w
      glas_path$wi[,,i] <- res$wi
    }
  }


  # Compute EBICs:
  #     EBICs <- apply(glas_path$wi,3,function(C){
  #       EBIC(S, C, n, gamma)
  #     })

  lik <- sapply(seq_along(lambda),function(i){
    logGaus(S, glas_path$wi[,,i], n)
  })

  EBICs <- sapply(seq_along(lambda),function(i){
    EBIC(S, glas_path$wi[,,i], n, gamma, countDiagonal=countDiagonal)
  })

  # Smallest EBIC:
  opt <- which.min(EBICs)

  # Check if rho is smallest:
  if (opt == 1){
    warning("Network with lowest lambda selected as best network. Try setting 'lambda.min.ratio' lower.")
  }

  # Return network:
  net <- as.matrix(Matrix::forceSymmetric(wi2net(glas_path$wi[,,opt])))
  colnames(net) <- rownames(net) <- colnames(S)

  # Check empty network:
  if (all(net == 0)){
    message("An empty network was selected to be the best fitting network. Possibly set 'lambda.min.ratio' higher to search more sparse networks. You can also change the 'gamma' parameter to improve sensitivity (at the cost of specificity).")
  }

  # Refit network:
  # Refit:
  if (refit){
    message("Refitting network without LASSO regularization")
    glassoRes <- suppressWarnings(glasso::glasso(S, 0, zero = which(net == 0 & upper.tri(net), arr.ind=TRUE), trace = 0, penalize.diagonal=penalize.diagonal, ...))
    net <- as.matrix(Matrix::forceSymmetric(wi2net(glassoRes$wi)))
    colnames(net) <- rownames(net) <- colnames(S)
    optwi <- glassoRes$wi
  } else {
    optwi <- glas_path$wi[,,opt]
  }

  # Return
  if (returnAllResults){
    return(list(
      results = glas_path,
      ebic = EBICs,
      loglik = lik,
      optnet = net,
      lambda = lambda,
      optwi = optwi
    ))
  } else return(net)
}
