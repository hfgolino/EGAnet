#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Optimization functions for EGM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%
## SRMR ----
#%%%%%%%%%%%

#' @noRd
# SRMR
# Updated 21.07.2025
srmr <- function(base, comparison, power = 2)
{

  # Obtain lower triangle
  lower_triangle <- lower.tri(base)

  # Return SRMR
  return(sqrt(mean((base[lower_triangle] - comparison[lower_triangle])^power)))

}

#' @noRd
# SRMR cost
# Updated 30.07.2025
srmr_cost <- function(loadings_vector, R, rows, n, v, lambda, ...)
{

  # Set regularization penalty
  penalty <- lambda * sum(loadings_vector^2)

  # Assemble loading matrix
  loadings_matrix <- matrix(loadings_vector, ncol = rows)

  # Convert to correlation matrix
  S <- tcrossprod(loadings_matrix) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  diag_S <- 1 / diag(S) # save diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations


  # Return cost
  return(
    swiftelse(
      is_positive_definite(implied_R), # check for positive definite
      sqrt(mean((implied_R - R)^2)) + penalty,
      1e10 # return horrible value if not positive definite
    )
  )

}

#' @noRd
# SRMR gradient
# Updated 30.07.2025
srmr_gradient <- function(loadings_vector, R, rows, n, v, lambda, ...)
{

  # Set regularization penalty
  penalty <- lambda * 2 * loadings_vector

  # Assemble loading matrix
  loadings_matrix <- matrix(loadings_vector, ncol = rows)

  # Convert to correlation matrix
  S <- tcrossprod(loadings_matrix) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  diag_S <- 1 / diag(S) # save diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations
  error <- implied_R - R

  # Compute error
  dError <- error / (v^2 * sqrt(mean(error^2)))

  # Gradient for S
  dS <- D %*% dError %*% D

  # Return gradient
  return(as.vector(t(2 * crossprod(loadings_matrix, dS))) + penalty)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Gaussian log-likelihood ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Log-likelihood only
# Updated 06.10.2024
log_likelihood <- function(n, p, R, S, type = c("partial", "zero"))
{

  # Set default to zero-order
  if(missing(type)){
    type <- "zero"
  }else{type <- match.arg(type)}

  # Return
  return(
    swiftelse(
      type == "zero",
      -(n / 2) * (p * log(2 * pi) + log(det(R)) + sum(diag(S %*% solve(R)))),
      (n / 2) * (log(det(R)) - sum(diag((S %*% R)))) - (n * p / 2) * log(2 * pi)
    )
  )

}

#' @noRd
# Log-likelihood cost
# Updated 30.07.2025
logLik_cost <- function(loadings_vector, R, rows, n, v, lambda, ...)
{

  # Set regularization penalty
  penalty <- lambda * sum(loadings_vector^2)

  # Assemble loading matrix
  loadings_matrix <- matrix(loadings_vector, ncol = rows)

  # Convert to correlation matrix
  S <- tcrossprod(loadings_matrix) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  diag_S <- 1 / diag(S) # save diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations

  # Return cost
  return(
    swiftelse(
      is_positive_definite(implied_R), # check for positive definite
      ((n / 2) * (v * log(2 * pi) + log(det(implied_R)) + sum(diag(R %*% solve(implied_R))))) + penalty,
      1e10 # return horrible value if not positive definite
    )
  )

}

#' @noRd
# Log-likelihood gradient
# Updated 30.07.2025
logLik_gradient <- function(loadings_vector, R, rows, n, v, lambda, ...)
{

  # Set regularization penalty
  penalty <- lambda * 2 * loadings_vector

  # Assemble loading matrix
  loadings_matrix <- matrix(loadings_vector, ncol = rows)

  # Convert to correlation matrix
  S <- tcrossprod(loadings_matrix) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  diag_S <- 1 / diag(S) # save diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations
  K <- solve(implied_R) # precision matrix

  # Compute error
  dError <- (n / 2) * (K - K %*% R %*% K)

  # Gradient for S
  dS <- D %*% dError %*% D

  # Return gradient
  return(as.vector(t(2 * crossprod(loadings_matrix, dS))) + penalty)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM optimization ----
# Updated 30.07.2025
egm_optimize <- function(
    loadings_vector, zeros,
    R, rows, n, v, lambda, opt,
    iterations = 10000, ...
)
{

  # Get optimization functions
  cost <- switch(opt, "loglik" = logLik_cost, "srmr" = srmr_cost)
  gradient <- switch(opt, "loglik" = logLik_gradient, "srmr" = srmr_gradient)

  # Set bounds for Hessian
  bounds <- 0.70 * zeros

  # Obtain result
  result <- try(
    silent_call(
      nlminb(
        start = loadings_vector,
        objective = cost, gradient = gradient,
        R = R, rows = rows, n = n, v = v, lambda = lambda,
        lower = -bounds, upper = bounds,
        control = list(
          eval.max = iterations, iter.max = iterations,
          step.min = 1e-12, step.max = 0.10
        )
      )
    ), silent = TRUE
  )

  # Check for issues, if none, then obtain hessian
  if(!is(result, "try-error")){

    # Attach hessian
    result$hessian <- optimHess(
      par = result$par, fn = cost, gr = gradient,
      R = R, rows = rows, n = n, v = v, lambda = lambda,
      lower = -bounds, upper = bounds
    )

  }

  # Return result
  return(result)

}

#' @noRd
# Hessian optimization ----
# Updated 20.06.2025
hessian_optimize <- function(
    lambda, loadings_vector, zeros,
    R, rows, n, v, opt, iterations = 10000, ...
)
{

  # Optimize over loadings
  result <- try(
    egm_optimize(
      loadings_vector = loadings_vector, zeros = zeros,
      R = R, rows = rows, n = n, v = v, opt = opt,
      lambda = lambda, iterations = iterations, ...
    ), silent = TRUE
  )

  # Get error flag
  error_flag <- is(result, "try-error")

  # Check for error
  if(is(result, "try-error")){
    return(1e12)
  }else{

    # Obtain hessian eigenvalues
    hessian_eigenvalue <- min(matrix_eigenvalues(result$hessian))
    scaled_eigenvalue <- sqrt(abs(hessian_eigenvalue))

    # Set hessian penalty
    hessian_penalty <- 0

    # Set scaling of penalty based on the objective function
    scaling <- 10^(digits(result$objective) - 1)

    # Update penalty
    if(hessian_eigenvalue < 0){
      hessian_penalty <- scaled_eigenvalue * scaling
    }else if(hessian_eigenvalue > 1){
      hessian_penalty <- scaled_eigenvalue * scaling * 0.10
    }else if(hessian_eigenvalue > 0.10){
      hessian_penalty <- scaled_eigenvalue * scaling * 0.01
    }else if(hessian_eigenvalue > 0.01){
      hessian_penalty <- scaled_eigenvalue * scaling * 0.001
    }

    # Add hessian to objective
    return(result$objective + hessian_penalty)
  }

}

#' @noRd
# Loading optimization ----
# Updated 22.07.2025
loadings_optimization <- function(
    iter = 10, loadings_vector, zeros, R,
    communities, data_dimensions, opt
){

  # Initial best parameters
  best_result <- list(par = loadings_vector)
  best_eigenvalue <- Inf; negative_eigenvalue <- -Inf

  # Set lambda minimum and maximum
  loadings_lambda_min <- 0
  loadings_lambda_max <- 10

  # Seek out positive eigenvalue
  for(i in seq_len(iter)){

    # Optimize for best quality solution (quick passes)
    lambda <- optimize(
      f = hessian_optimize, interval = c(loadings_lambda_min, loadings_lambda_max),
      loadings_vector = best_result$par, zeros = zeros,
      R = R, rows = communities, n = data_dimensions[1],
      v = data_dimensions[2], opt = opt, tol = 1e-03, iterations = 100
    )

    # Optimize over loadings
    result <- try(
      egm_optimize(
        loadings_vector = best_result$par, zeros = zeros,
        R = R, rows = communities, n = data_dimensions[1],
        v = data_dimensions[2], opt = opt, lambda = lambda$minimum
      ), silent = TRUE
    )

    # Check bad results
    if(is(result, "try-error")){
      return(NULL)
    }

    # Store eigenvalue
    current_eigenvalue <- min(matrix_eigenvalues(result$hessian))

    # Check for positive eigenvalue
    if(current_eigenvalue > 0){

      # Check if eigenvalues are better than previous
      if(current_eigenvalue < best_eigenvalue){

        # Update best eigenvalue
        best_eigenvalue <- current_eigenvalue

        # Update best results
        best_result <- result

      }

      # Check if best result is at target
      if(best_eigenvalue < 0.01){
        break
      }

      # Shrink lambda range
      interval_range <- (loadings_lambda_max - loadings_lambda_min) * 0.25

      # Set minimum and maximum
      loadings_lambda_min <- max(loadings_lambda_min, lambda$minimum - interval_range)
      loadings_lambda_max <- min(loadings_lambda_max, lambda$minimum + interval_range)

    }else{ # Negative eigenvalue

      # Check if eigenvalues are better than previous
      if(current_eigenvalue > negative_eigenvalue){

        # Keep results moving toward positive eigenvalues

        # Update best eigenvalue
        negative_eigenvalue <- current_eigenvalue

        # Update best results
        best_result <- result

        # No need to update lambda minimum or maximum

      }

    }

  }

  # Add lambda to best result
  best_result$lambda <- lambda$minimum

  # Return best result
  return(best_result)

}

#%%%%%%%%%%%%%%%#
#### NETWORK ####
#%%%%%%%%%%%%%%%#

#%%%%%%%%%%%
## SRMR ----
#%%%%%%%%%%%

#' @noRd
# SRMR network cost
# Updated 31.07.2025
srmr_network_cost <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  S <- solve(-network) # covariance matrix
  D <- diag(sqrt(1 / diag(S))) # standardization
  implied_R <- D %*% S %*% D # implied correlations

  # Return cost
  return(
    swiftelse(
      is(implied_R, "try-error") || anyNA(implied_R) || !is_positive_definite(implied_R),
      1e10, # return horrible value if not positive definite
      sqrt(mean((implied_R - R)^2))
    )
  )

}

#' @noRd
# SRMR network gradient
# Updated 30.07.2025
srmr_network_gradient <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  S <- solve(-network) # covariance matrix
  diag_S <- 1 / diag(S) # store diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations

  # Compute error
  dError <- 2 * ((implied_R - R) / sum(lower_triangle))

  # Gradient for S
  dS <- D %*% dError %*% D

  # Diagonal corrections
  diag(dS) <- diag(dS) - diag_S * rowSums(dError * implied_R)

  # Return gradient
  return(2 * (S %*% dS %*% S)[lower_triangle][zeros])

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Function ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM network optimization
# Updated 01.08.2025
egm_network_optimize <- function(
    network_vector, R, n, v,
    lower_triangle, zeros,
    max_iterations = 10000, ...
)
{

  # Set bounds
  bounds <- rep(1, length(network_vector))

  return(
    silent_call(
      nlminb(
        start = network_vector,
        objective = srmr_network_cost,
        gradient = srmr_network_gradient,
        R = R, n = n, v = v, lower_triangle = lower_triangle, zeros = zeros,
        lower = -bounds, upper = bounds,
        control = list(
          eval.max = max_iterations, iter.max = max_iterations,
          step.min = .Machine$double.eps, step.max = 1,
          abs.tol = .Machine$double.eps, rel.tol = .Machine$double.eps,
          x.tol = .Machine$double.eps, xf.tol = .Machine$double.eps,
          sing.tol = .Machine$double.eps
        )
      )
    )
  )

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### COMMUNITY CORRELATIONS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Mean absolute error cost
# Updated 31.07.2025
mae_community_cost <- function(
    between_loadings, between_indices, simple_structure,
    loading_structure, target_correlations,
    membership, total_variables, ...
)
{

  # Set between-community loadings
  loading_structure[between_indices] <- between_loadings

  # Convert to correlation matrix
  S <- tcrossprod(loading_structure) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  diag_S <- 1 / diag(S) # save diagonal
  D <- diag(sqrt(diag_S)) # standardization
  implied_R <- D %*% S %*% D # implied correlations

  # Get covariance and convert to correlations
  S_corr <- crossprod(simple_structure, D %*% S %*% D) %*% simple_structure
  D_corr <- diag(sqrt(1 / diag(S_corr)))

  # Check expected network
  network <- expected_network(loading_structure, membership, total_variables)
  network_R <- try(pcor2cor(network), silent = TRUE)

  # Return MAE
  return(
    swiftelse(
      is(network_R, "try-error") || anyNA(network_R) || !is_positive_definite(network_R),
      1e10, # send horrible cost
      mean(abs(D_corr %*% S_corr %*% D_corr - target_correlations))
    )
  )

}

#' @noRd
# Mean absolute error gradient
# Updated 31.07.2025
mae_community_gradient <- function(
    between_loadings, between_indices, simple_structure,
    loading_structure, target_correlations,
    membership, total_variables, ...
)
{

  # Set between-community loadings
  loading_structure[between_indices] <- between_loadings

  # Convert to correlation matrix
  S <- tcrossprod(loading_structure) # compute covariances
  diag(S) <- sqrt(diag(S)) # interdependence
  D <- diag(sqrt(1 / diag(S))) # standardization
  implied_R <- D %*% S %*% D # implied correlations

  # Get covariance and convert to correlations
  S_corr <- crossprod(simple_structure, D %*% S %*% D) %*% simple_structure
  diag_D_corr <- 1 / diag(S_corr)
  D_corr <- diag(sqrt(diag_D_corr))
  R_corr <- D_corr %*% S_corr %*% D_corr

  # Error of gradient
  differences <- R_corr - target_correlations
  dError <- sign(differences) / length(differences)
  dError[differences == 0] <- 0

  # Standardization gradient
  dS_corr <- D_corr %*% dError %*% D_corr

  # Add correction
  diag(dS_corr) <- diag(dS_corr) * diag_D_corr / 2

  # Implied correlations gradient
  dImplied_R <- simple_structure %*% tcrossprod(dS_corr, simple_structure)

  # Initial standardization
  dS <- D %*% dImplied_R %*% D

  # Return gradient
  return(as.vector(t(crossprod(loading_structure, dS)))[between_indices])

}

#' @noRd
# EGM community correlation optimization
# Updated 31.07.2025
egm_correlation_optimize <- function(
    loading_structure, correlations, between_indices,
    membership, total_variables, ...
)
{

  # Initialize simple structure
  simple_structure <- loading_structure

  # Bounds (can't be larger than assigned loadings)
  max_assigned <- apply(simple_structure, 1, function(x){max(abs(x))}) - 0.001
  bounds <- as.vector(max_assigned * between_indices)
  bounds <- bounds[between_indices]

  # Optimize for loadings
  cross_loadings <- silent_call(
    nlminb(
      start = runif_xoshiro(sum(between_indices), min = -0.10, max = 0.10),
      objective = mae_community_cost,
      gradient = mae_community_gradient,
      between_indices = between_indices,
      simple_structure = simple_structure,
      loading_structure = loading_structure,
      target_correlations = correlations,
      membership = membership,
      total_variables = total_variables,
      lower = -bounds, upper = bounds,
      control = list(
        eval.max = 10000, iter.max = 10000,
        abs.tol = .Machine$double.eps,
        rel.tol = .Machine$double.eps,
        x.tol = .Machine$double.eps
      )
    )
  )

  # Add back cross-loadings
  loading_structure[between_indices] <- cross_loadings$par

  # Return loadings and exact community correlations
  return(
    list(
      loadings = loading_structure,
      correlations = community_correlations(simple_structure, loading_structure)
    )
  )

}