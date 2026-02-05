#' @title Re-fit Network
#'
#' @description Refits a network with unregularized partial correlations
#' based on some known graph (perhaps estimated with a regularization method)
#' following Hastie, Tibshirani, and Friedman (2008)
#'
#' @param S Matrix.
#' Covariance or correlation matrix
#'
#' @param A Adjacency matrix.
#' Unweighted network structure where \code{1} is an edge present
#' and \code{0} is an edge absent
#'
#' @param method Character (length = 1).
#' Whether to use the \code{\link[glasso]{glasso}} method without
#' penalization or the HTF (Haste, Tibshirani, & Friedman, 2008) method.
#' Defaults to \code{"glasso"}, which tends to be more robust
#'
#' @param tol Numeric (length = 1).
#' Tolerance for convergence. The algorithm stops when the maximum
#' absolute change in covariance matrix elements between iterations
#' is less than \code{tol}. Defaults to \code{1e-06}
#'
#' @param max.iter Numeric (length = 100).
#' Maximum number of iterations to achieve tolerance before stopping
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#'
#' @return Returns a list containing:
#'
#' \item{network}{Estimated network}
#'
#' \item{W}{Estimated covariance matrix}
#'
#' \item{Theta}{Estimated inverse covariance matrix}
#'
#' \item{iterations}{Number of iterations to converge (or maximum if it did not)}
#'
#' \item{converged}{Whether the algorithm converged}
#'
#' @examples
#' # Obtain data
#' wmt <- depression[,24:44]
#'
#' # Obtain correlation matrix
#' wmt_R <- auto.correlate(wmt)
#'
#' # Estimate network
#' wmt_network <- network.estimation(wmt_R, n = nrow(wmt))
#'
#' # Obtain adjacency
#' wmt_A <- wmt_network
#' wmt_A[] <- ifelse(wmt_A != 0, 1, 0)
#'
#' # Obtain unregularized estimate
#' wmt_unreg <- known.graph(S = wmt_R, A = wmt_A)
#'
#' @references
#' \strong{HTF Implementation on p. 631--634} \cr
#' Hastie, T., Tibshirani, R., & Friedman, J. (2008).
#' The elements of statistical learning: Data mining, inference, and prediction (2nd ed.).
#' New York, NY: Springer.
#'
#' @export
# Known graph ----
# Updated 05.02.2026
known.graph <- function(S, A, method = c("glasso", "HTF"), tol = 1e-06, max.iter = 100)
{

  # Check for method
  method <- set_default(method, "glasso", known.graph)

  # Run based on method
  if(method == "glasso"){

    # Obtain indices
    zeros <- which(A == 0, arr.ind = TRUE)
    zeros <- zeros[zeros[,1] <= zeros[,2],]

    # Obtain GLASSO output
    output <- silent_call(
      glasso::glasso(
        s = S, rho = 0, zero = zeros, thr = tol,
        maxit = max.iter, trace = 0, penalize.diagonal = FALSE
      )
    )

    # Set results
    W <- output$w
    Theta <- output$wi
    iteration <- output$niter
    converged <- output$del < tol

  }else{

    # Obtain number of nodes
    nodes <- dim(S)[2]

    # Get node sequence
    node_sequence <- seq_len(nodes)

    # Initialize matrix
    W <- S

    # Initialize zero betas
    zero_beta <- matrix(0, nrow = nodes - 1, ncol = 1)

    # Obtain non-zero edge list
    non_zero_list <- lapply(node_sequence, function(i){
      A[-i, i] != 0
    })

    # Initialize convergence criteria
    iteration <- 0; difference <- Inf

    # Loop until criterion is met or max iterations is reached
    while((difference > tol) && (iteration < max.iter)){

      # Increase iterations
      iteration <- iteration + 1

      # Set old W
      W_old <- W

      # Iterate over each node
      for(j in node_sequence){

        # Identify edges in the adjacency matrix for node j
        edges <- non_zero_list[[j]]

        # Initialize beta
        beta <- zero_beta

        # Partition for W11
        W11 <- W[-j, -j, drop = FALSE]

        # Extract the corresponding rows/columns in W11 and elements of S12
        if(any(edges)){

          # Compute S12
          S12 <- S[-j, j, drop = FALSE]

          # Solve for beta* = (W11*)^{-1} * S12*
          beta[edges,] <- solve(
            W11[edges, edges, drop = FALSE], S12[edges, , drop = FALSE]
          )

          # Update W12 = W11 * beta
          W12_new <- W11 %*% beta

          # Update the corresponding parts of W
          W[j, -j] <- W[-j, j] <- W12_new
          W[j, j] <- S[j, j] - sum(beta * W12_new)

        }

      }

      # Compute maximum absolute difference in W for convergence check
      difference <- max(abs(W - W_old))

    }

    # Obtain Theta and convergence
    Theta <- solve(W)
    converged <- difference <= tol

  }

  # Obtain P
  P <- -cov2cor(Theta); diag(P) <- 0

  # Return all output
  return(
    list(
      network = P, W = W, Theta = Theta,
      iterations = iteration,
      converged = converged

    )
  )

}
