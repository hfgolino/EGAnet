#' @title Regularized Networks with Convex and Non-convex Penalties
#'
#' @description A general function to estimate Gaussian graphical models using
#' regularization penalties. All non-convex penalties are implemented using
#' the Local Linear Approximation (LLA: Fan & Li, 2001; Zou & Li, 2008)
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
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"atan"} --- Arctangent (Wang & Zhu, 2016)
#' \deqn{\lambda \cdot (\gamma + 2 \pi) \cdot \arctan(\frac{|x|}{\gamma})}
#'
#' \item \code{"bridge"} --- Bridge (Fu, 1998)
#' \deqn{\lambda \cdot |x|^\gamma}
#'
#' \item \code{"cauchy"} ---- Cauchy
#' \deqn{\lambda \cdot \frac{1}{\pi} \cdot \arctan{\frac{|x|}{\gamma}} + 0.5}
#'
#' \item \code{"exp"} --- EXP (Wang, Fan, & Zhu, 2018)
#' \deqn{\lambda \cdot (1 - e^{-\frac{|x|}{\gamma}})}
#'
#' \item \code{"l1"} --- LASSO (Tibshirani, 1996)
#' \deqn{\lambda \cdot |x|}
#'
#' \item \code{"l2"} --- Ridge (Hoerl & Kennard, 1970)
#' \deqn{\lambda \cdot x^2}
#'
#' \item \code{"mcp"} --- Minimax Concave Penalty (Zhang, 2010)
#' \deqn{
#' P(x; \lambda, \gamma) =
#' \begin{cases}
#' \lambda |x| - \frac{x^2}{2\gamma} & \text{if } |x| \leq \gamma\lambda \\
#' \frac{\gamma \lambda^2}{2} & \text{if } |x| > \gamma\lambda
#' \end{cases}
#' }
#'
#' \item \code{"scad"} --- Smoothly Clipped Absolute Deviation (Fan & Li, 2001)
#' \deqn{
#' P(x; \lambda, \gamma) =
#' \begin{cases}
#' \lambda |x| & \text{if } |x| \leq \lambda \\
#' -\frac{|x|^2 - 2\gamma\lambda|x| + \lambda^2}{2(\gamma - 1)} & \text{if } \lambda < |x| \leq \gamma\lambda \\
#' \frac{(\gamma + 1)\lambda^2}{2} & \text{if } |x| > \gamma\lambda
#' \end{cases}
#' }
#'
#' \item \code{"weibull"} --- Data-adaptive Weibull
#' \deqn{\lambda \cdot (1 - e^{\large(-\frac{|x|}{\gamma}\large)^k})}
#'
#' }
#'
#' @param adaptive Boolean (length = 1).
#' Whether data-adaptive parameters should be used.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to apply data-adaptive parameters
#' based on the empirical partial correlation matrix.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"cauchy"} = uses half of the interquartile range of the absolute empirical partial correlations (Bloch, 1966)
#'
#' \item \code{"exp"} =  uses median of distribution for the scale parameter (\eqn{\frac{\log{(2)}}{\lambda}})
#'
#' \item \code{"weibull"} = uses MLE estimate of shape parameter and median of distribution for the scale parameter (\eqn{\lambda \cdot (\log{(2)})^{1/k} })
#'
#' }
#'
#' @param gamma Numeric (length = 1).
#' Adjusts the shape of the penalty.
#' Defaults:
#'
#' \itemize{
#'
#' \item \code{"atan"} = 0.01
#'
#' \item \code{"bridge"} = 1
#'
#' \item \code{"cauchy"} = 0.01
#'
#' \item \code{"exp"} = 0.01
#'
#' \item \code{"mcp"} = 3
#'
#' \item \code{"scad"} = 3.7
#'
#' \item \code{"weibull"} = 0.01
#'
#' }
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
#' @param optimize.lambda Boolean (length = 1).
#' Whether optimization of lambda should be performed.
#' Defaults to \code{FALSE} or grid search over lambda.
#' If \code{TRUE}, then \code{\link[stats]{optimize}} is used
#' to find the optimal lambda
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
#' \item \code{"BIC0"} --- Bayesian information criterion not (Dicker et al., 2013): \eqn{\log{\large(\frac{D}{n - E}\large)} + \large(\frac{\log{(n)}}{n}\large) \cdot E}
#'
#' \item \code{"EBIC"} --- Extended BIC: \eqn{BIC + 4E \cdot \gamma \cdot \log{(E)}}
#'
#' \item \code{"MBIC"} --- Modified Bayesian information criterion (Wang et al., 2018):  \eqn{\log{\large(\frac{D}{n - E}\large)} + \large(\frac{\log{(n)} \cdot E}{n}\large) \cdot \log{(\log{(p)}})}
#'
#' }
#'
#' Term definitions:
#'
#' \itemize{
#'
#' \item \eqn{n} --- sample size
#'
#' \item \eqn{p} --- number of variables
#'
#' \item \eqn{E} --- edges
#'
#' \item \eqn{S} --- empirical correlation matrix
#'
#' \item \eqn{K} --- estimated inverse covariance matrix (network)
#'
#' \item \eqn{L = \frac{n}{2} \cdot \log \text{det} K - \sum_{i=1}^p (SK)_{ii}}
#'
#' \item \eqn{D = n \cdot \sum_{i=1}^p (SK)_{ii} - \log \text{det} K}
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
#' @param LLA Boolean (length = 1).
#' Should Local Linear Approximation be used to find optimal minimum?
#' Defaults to \code{FALSE} or a single-pass approximation, which can be
#' significantly faster (Zou & Li, 2008).
#' Set to \code{TRUE} to find global minimum based on convergence (\code{LLA.threshold})
#'
#' @param LLA.threshold Numeric (length = 1).
#' When performing the Local Linear Approximation, the maximum threshold
#' until convergence is met.
#' Defaults to \code{1e-04}
#'
#' @param LLA.iter Numeric (length = 1).
#' Maximum number of iterations to perform to reach convergence.
#' Defaults to \code{100}
#'
#' @param network.only Boolean (length = 1).
#' Whether the network only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' network estimation method
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to \code{\link[EGAnet]{auto.correlate}}
#'
#' @return A network matrix
#'
#' @author Alexander P. Christensen <alexpaulchristensen at gmail.com> and
#' Hudson Golino <hfg9s at virginia.edu>
#'
#' @references
#'
#' \strong{Half IQR for \eqn{\gamma} in Cauchy} \cr
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1970).
#' Continuous univariate distributions (Vol. 1).
#' New York, NY: John Wiley & Sons.
#'
#' \strong{BIC0} \cr
#' Dicker, L., Huang, B., & Lin, X. (2013).
#' Variable selection and estimation with the seamless-L0 penalty.
#' \emph{Statistica Sinica}, \emph{23}(2), 929--962.
#'
#' \strong{SCAD penalty and Local Linear Approximation} \cr
#' Fan, J., & Li, R. (2001).
#' Variable selection via nonconcave penalized likelihood and its oracle properties.
#' \emph{Journal of the American Statistical Association}, \emph{96}(456), 1348--1360.
#'
#' \strong{Bridge penalty} \cr
#' Fu, W. J. (1998).
#' Penalized regressions: The bridge versus the lasso.
#' \emph{Journal of Computational and Graphical Statistics}, \emph{7}(3), 397--416.
#'
#' \strong{L2 penalty} \cr
#' Hoerl, A. E., & Kennard, R. W. (1970).
#' Ridge regression: Biased estimation for nonorthogonal problems.
#' \emph{Technometrics}, \emph{12}(1), 55--67.
#'
#' \strong{L1 penalty} \cr
#' Tibshirani, R. (1996).
#' Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \emph{58}(1), 267--288.
#'
#' \strong{EXP penalty} \cr
#' Wang, Y., Fan, Q., & Zhu, L. (2018).
#' Variable selection and estimation using a continuous approximation to the L0 penalty.
#' \emph{Annals of the Institute of Statistical Mathematics}, \emph{70}(1), 191--214.
#'
#' \strong{Atan penalty} \cr
#' Wang, Y., & Zhu, L. (2016).
#' Variable selection and parameter estimation with the Atan regularization method.
#' \emph{Journal of Probability and Statistics}, \emph{2016}, 1--12.
#'
#' \strong{Original simulation in psychometric networks} \cr
#' Williams, D. R. (2020).
#' Beyond lasso: A survey of nonconvex regularization in Gaussian graphical models.
#' \emph{PsyArXiv}.
#'
#' \strong{MCP penalty} \cr
#' Zhang, C.-H. (2010).
#' Nearly unbiased variable selection under minimax concave penalty.
#' \emph{Annals of Statistics}, \emph{38}(2), 894--942.
#'
#' \strong{One-step Local Linear Approximation} \cr
#' Zou, H., & Li, R. (2008).
#' One-step sparse estimates in nonconcave penalized likelihood models.
#' \emph{Annals of Statistics}, \emph{36}(4), 1509--1533.
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain network
#' l1_network <- network.regularization(data = wmt)
#'
#' # Obtain Atan network
#' atan_network <- network.regularization(data = wmt, penalty = "atan")
#'
#' # Obtain data-adaptive EXP network
#' exp_network <- network.regularization(data = wmt, penalty = "exp")
#'
#' @export
#'
# Apply non-convex regularization ----
# Updated 12.01.2026
network.regularization <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "cosine", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    penalty = c("atan", "bridge", "cauchy", "exp", "l1", "l2", "mcp", "scad", "weibull"),
    adaptive = FALSE, gamma = NULL, lambda = NULL, nlambda = 50, lambda.min.ratio = 0.01,
    penalize.diagonal = TRUE, optimize.lambda = FALSE,
    ic = c("AIC", "AICc", "BIC", "BIC0", "EBIC", "MBIC"), ebic.gamma = 0.50,
    fast = TRUE, LLA = FALSE, LLA.threshold = 1e-04, LLA.iter = 100,
    network.only = TRUE, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network.regularization)
  na.data <- set_default(na.data, "pairwise", network.regularization)
  penalty <- set_default(penalty, "l1", network.regularization)
  ic <- set_default(ic, "bic", network.regularization)

  # Argument errors (return data in case of tibble)
  data <- network.regularization_errors(
    data, n, adaptive, gamma, nlambda, lambda.min.ratio,
    penalize.diagonal, optimize.lambda, ebic.gamma,
    fast, LLA, LLA.threshold, network.only, verbose, ...
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
    "atan" = atan_derivative,
    "bridge" = bridge_derivative,
    "cauchy" = cauchy_derivative,
    "exp" = exp_derivative,
    "l1" = l1_derivative,
    "l2" = l2_derivative,
    "mcp" = mcp_derivative,
    "scad" = scad_derivative,
    "weibull" = weibull_derivative
  )

  # Set shape (only used for Weibull)
  shape <- 1

  # Set gamma (set ahead of time for messaging on adaptive)
  if(is.null(gamma)){

    # Set defaults
    gamma <- switch(
      penalty,
      "atan" = 0.01,
      "bridge" = 1,
      "cauchy" = 0.01,
      "exp" = 0.01,
      "mcp" = 3,
      "scad" = 3.7,
      "weibull" = 0.01
    )

  }

  # Initialize adaptive lambda
  adaptive_lambda <- FALSE

  # Check whether penalty is adaptive option
  adaptive_option <- c("cauchy", "exp", "weibull")

  # Check for adaptive and penalty is adaptive option
  if(adaptive){

    # Check whether penalty is adaptive
    if(penalty %in% adaptive_option){

      # Send adaptive lambda
      adaptive_lambda <- TRUE

      # Set lower triangle
      lower_triangle <- lower.tri(S)

      # Set partial correlations
      P <- cor2pcor(S); lower_P <- abs(P[lower_triangle])

      if(penalty == "cauchy"){

        # Set gamma parameter
        gamma <- diff(quantile(lower_P, probs = c(0.25, 0.75))) / 2

      }else if(penalty == "exp"){

        # Obtain median of distribution
        gamma <- 0.6931472 * sum(lower_P) / sum(lower_triangle)
        # pre-computes log(2) = 0.6931472

      }else if(penalty == "weibull"){

        # Obtain Weibull estimates
        estimates <- weibull_mle(lower_P)

        # Set parameters
        shape <- min(estimates[["shape"]], 1) # cap at EXP
        gamma <- estimates[["scale"]] * 0.6931472^(1 / shape) # median
        # pre-computes log(2) = 0.6931472

      }

    }else{

      # Set message that penalty is not adaptive
      message(
        paste0(
          "The \"", penalty,
          "\" penalty is not adaptive.\n\nUsing its default, `gamma = ", gamma, "`"
        )
      )

    }

  }

  # Set LLA to FALSE for exponential penalties
  if(penalty %in% adaptive_option){
    LLA <- FALSE
  }

  # Initialize lambda matrix
  lambda_matrix <- matrix(0, nrow = nodes, ncol = nodes)

  # Check for optimization
  if(optimize.lambda){

    # Optimize for lambda
    optimized_lambda <- optimize(
      f = lambda_optimize, interval = c(0, swiftelse(penalty == "l2", 10, 1)),
      gamma = gamma, K = K, S = S,
      derivative_FUN = derivative_FUN,
      glasso_FUN = glasso_FUN, glasso_ARGS = glasso_ARGS,
      lambda_matrix = lambda_matrix, penalize.diagonal = penalize.diagonal,
      ic = ic, n = n, nodes = nodes, ebic.gamma = ebic.gamma, shape = shape
    )

    # Obtain lambda matrix
    lambda_matrix[] <- abs(derivative_FUN(
      x = K, lambda = optimized_lambda$minimum, gamma = gamma, shape = shape
    ))

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
    if(network.only){
      return(W)
    }else{
      return(
        list(
          network = W, K = output$wi, R = R,
          penalty = penalty, lambda = optimized_lambda$minimum,
          gamma = gamma, criterion = ic,
          IC = optimized_lambda$objective, correlation = S
        )
      )
    }

  }else{

    # Simplify source for fewer computations (minimal improvement)
    S_zero_diagonal <- S - diag(nodes) # makes diagonal zero
    lambda.max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
    lambda.max <- lambda.max / ifelse(adaptive_lambda, log10(n), 1) # adapt with sample size
    lambda.min <- lambda.min.ratio * lambda.max
    lambda <- exp(seq.int(log(lambda.min), log(lambda.max), length.out = nlambda))

    # Obtain lambda sequence
    lambda_sequence <- seq_len(nlambda)

    # Obtain lambda matrices
    lambda_list <- lapply(lambda, function(value){

      # Obtain lambda matrix
      lambda_matrix[] <- abs(derivative_FUN(x = K, lambda = value, gamma = gamma, shape = shape))

      # Check for diagonal penalization
      if(!penalize.diagonal){
        diag(lambda_matrix) <- 0
      }

      # Attach value
      attr(lambda_matrix, "value") <- value

      # Return lambda matrix
      return(lambda_matrix)

    })

    # Check for LLA
    if(LLA){

      # Get GLASSO output
      glasso_list <- lapply(lambda_list, function(lambda_matrix){

        # Obtain lambda value
        value <- attributes(lambda_matrix)$value

        # Set lambda matrix
        glasso_ARGS$rho <- lambda_matrix

        # Obtain estimate
        estimate <- do.call(what = glasso_FUN, args = glasso_ARGS)

        # Obtain new K
        new_K <- estimate$wi

        # Set convergence and iterations
        convergence <- Inf; iterations <- 0

        # Loop over to convergence
        while((convergence > LLA.threshold) & (iterations < LLA.iter)){

          # Set old K
          old_K <- new_K

          # Obtain lambda matrix
          lambda_matrix[] <- abs(derivative_FUN(x = old_K, lambda = value, gamma = gamma, shape = shape))

          # Check for diagonal penalization
          if(!penalize.diagonal){
            diag(lambda_matrix) <- 0
          }

          # Set lambda matrix
          glasso_ARGS$rho <- lambda_matrix

          # Obtain estimate
          estimate <- do.call(what = glasso_FUN, args = glasso_ARGS)

          # Obtain new K
          new_K <- estimate$wi

          # Increase iterations
          iterations <- iterations + 1

          # Compute convergence
          convergence <- mean(abs(new_K - old_K))

        }

        # Estimate
        return(estimate)

      })

    }else{

      # Get GLASSO output
      glasso_list <- lapply(lambda_list, function(lambda_matrix){

        # Set lambda matrix
        glasso_ARGS$rho <- lambda_matrix

        # Estimate
        return(do.call(what = glasso_FUN, args = glasso_ARGS))

      })

    }

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
    if(network.only){
      return(W)
    }else{
      return(
        list(
          network = W, K = glasso_list[[optimal]]$wi, R = R,
          penalty = penalty, lambda = lambda[[optimal]], gamma = gamma,
          correlation = S, criterion = ic, IC = ICs[[optimal]]
        )
      )
    }

  }

}

# Bug checking ----
# data = wmt2[,7:24]; n = NULL; corr = "auto"
# na.data = "pairwise"; penalty = "l1"; adaptive = FALSE
# gamma = NULL; lambda = NULL; nlambda = 50
# lambda.min.ratio = 0.01; penalize.diagonal = TRUE
# optimize.lambda = FALSE; ic = "BIC"; network.only = TRUE
# ebic.gamma = 0.5; fast = TRUE; verbose = FALSE
# LLA = FALSE; LLA.threshold = 1e-04; LLA.iter = 100


#' @noRd
# Errors ----
# Updated 12.01.2026
network.regularization_errors <- function(
    data, n, adaptive, gamma, nlambda, lambda.min.ratio,
    penalize.diagonal, optimize.lambda, ebic.gamma,
    fast, LLA, LLA.threshold, network.only, verbose, ...
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network.regularization")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network.regularization")
    typeof_error(n, "numeric", "network.regularization")
  }

  # 'adaptive' errors
  length_error(adaptive, 1, "network.regularization")
  typeof_error(adaptive, "logical", "network.regularization")

  # 'gamma' errors
  if(!is.null(gamma)){
    length_error(gamma, 1, "network.regularization")
    typeof_error(gamma, "numeric", "network.regularization")
    range_error(gamma, c(0, Inf), "network.regularization")
  }

  # 'nlambda' errors
  length_error(nlambda, 1, "network.regularization")
  typeof_error(nlambda, "numeric", "network.regularization")
  range_error(nlambda, c(1, Inf), "network.regularization")

  # 'lambda.min.ratio' errors
  length_error(lambda.min.ratio, 1, "network.regularization")
  typeof_error(lambda.min.ratio, "numeric", "network.regularization")
  range_error(lambda.min.ratio, c(0, 1), "network.regularization")

  # 'penalize.diagonal' errors
  length_error(penalize.diagonal, 1, "network.regularization")
  typeof_error(penalize.diagonal, "logical", "network.regularization")

  # 'optimize.lambda' errors
  length_error(optimize.lambda, 1, "network.regularization")
  typeof_error(optimize.lambda, "logical", "network.regularization")

  # 'ebic.gamma' errors
  length_error(ebic.gamma, 1, "network.regularization")
  typeof_error(ebic.gamma, "numeric", "network.regularization")
  range_error(ebic.gamma, c(0, Inf), "network.regularization")

  # 'fast' errors
  length_error(fast, 1, "network.regularization")
  typeof_error(fast, "logical", "network.regularization")

  # 'LLA' errors
  length_error(LLA, 1, "network.regularization")
  typeof_error(LLA, "logical", "network.regularization")

  # 'LLA.threshold' errors
  if(LLA){
    length_error(LLA.threshold, 1, "network.regularization")
    typeof_error(LLA.threshold, "numeric", "network.regularization")
    range_error(LLA.threshold, c(-Inf, 0.10), "network.regularization")
  }

  # 'network.only' errors
  length_error(network.only, 1, "network.regularization")
  typeof_error(network.only, "logical", "network.regularization")

  # 'verbose' errors
  length_error(verbose, 1, "network.regularization")
  typeof_error(verbose, "logical", "network.regularization")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

# OPTIMIZATION FUNCTIONS ----

#' @noRd
# MLE Weibull Parameters ----
# Updated 10.01.2026
weibull_mle <- function(x)
{

  # Set up MLE for shape
  shape_mle <- function(k, x, n)
  {

    # Pre-compute reused values
    x_k <- x^k
    log_x <- log(x)

    # Return MLE estimate
    return(sum(x_k * log_x) / sum(x_k) - 1 / k - sum(log_x) / n)

  }

  # Obtain MLE estimate
  shape <- uniroot(
    f = shape_mle, interval = c(0.01, 20),
    x = x, n = length(x)
  )$root

  # Return parameters
  return(c(shape = shape, scale = mean(x^shape)^(1 / shape)))

}

#' @noRd
# lambda optimization function ----
# Updated 24.11.2025
lambda_optimize <- function(
    lambda, gamma, K, S, derivative_FUN,
    glasso_FUN, glasso_ARGS,
    lambda_matrix, penalize.diagonal,
    ic, n, nodes, ebic.gamma, shape
)
{

  # Obtain lambda matrix
  lambda_matrix[] <- abs(
    derivative_FUN(x = K, lambda = lambda, gamma = gamma, shape = shape)
  )

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
  lambda_matrix[] <- derivative_FUN(x = K, lambda = lambda, gamma = gamma)

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
  lambda_matrix[] <- derivative_FUN(x = K, lambda = params[1], gamma = params[2])

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
# Updated 10.01.2026
information_crtierion <- function(S, K, n, nodes, ic, ebic.gamma)
{

  # Compute Gaussian likelihood (minus two for convenience)
  L <- swiftelse(
    ic %in% c("bic0", "mbic"),
    n * sum(diag(S %*% K)) - log(det(K)), # use deviance
    -2 * (n / 2) * (log(det(K)) - sum(diag(S %*% K))) # use log-likelihood
  )

  # Set diagonal of K to zero
  diag(K) <- 0

  # Get parameters (edges)
  E <- edge_count(K, nodes)

  # Ensure that there is enough degrees of freedom
  df <- n - E

  # Return information criterion
  return(
    switch(
      ic,
      "aic" = L + 2 * E,
      "aicc" = L + 2 * E + (2 * E^2 + 2 * E) / (n - E - 1),
      "bic" = L + E * log(n),
      "ebic" = L + E * log(n) + 4 * E * ebic.gamma * log(nodes),
      "bic0" = swiftelse( # see https://www.jstor.org/stable/24310368
        df > 0, log(L / df) + (log(n) * E / n), -Inf
      ),
      "mbic" = swiftelse( # see https://doi.org/10.1007%2Fs10463-016-0588-3
        df > 0, log(L / df) + (log(n) * E / n) * log(log(nodes)), -Inf
      )
    )
  )

}
