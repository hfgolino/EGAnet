#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Optimization functions for EGM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Obtain implied correlations ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Obtain implied correlations
# Updated 04.04.2025
obtain_implied <- function(loadings_vector, rows)
{

  # Assemble loading matrix
  loadings_matrix <- matrix(loadings_vector, ncol = rows)

  # Compute partial correlation
  implied_P <- tcrossprod(loadings_matrix)

  # Compute interdependence
  diag(implied_P) <- sqrt(diag(implied_P))

  # Compute matrix I
  I <- diag(sqrt(1 / diag(implied_P)))

  # Compute implied correlations
  implied_R <- I %*% implied_P %*% I

  # Attach loadings matrix
  attr(implied_R, which = "calculations") <- list(
    loadings_matrix = loadings_matrix, implied_P = implied_P, I = I
  )

  # Get implied R
  return(implied_R)

}

#%%%%%%%%%%%%%%%%#
#### LOADINGS ####
#%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%
## SRMR ----
#%%%%%%%%%%%

#' @noRd
# SRMR
# Updated 24.09.2024
srmr <- function(base, comparison)
{

  # Obtain lower triangle
  lower_triangle <- lower.tri(base)

  # Return SRMR
  return(sqrt(mean((base[lower_triangle] - comparison[lower_triangle])^2)))

}

#' @noRd
# SRMR cost
# Updated 21.06.2025
srmr_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    lambda, ...
)
{

  # Set l2 cost
  l2_cost <- lambda * sum(loadings_vector^2)

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Check for positive definite
    return(srmr(R, implied_R) + sum(differences > 0) + l2_cost)

  }else{ # Without constraints, send it
    return(srmr(R, obtain_implied(loadings_vector, rows)) + l2_cost)
  }

}

#' @noRd
# SRMR gradient
# Updated 21.06.2025
srmr_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    lambda, ...
)
{

  # Set l2 gradient
  l2_gradient <- 2 * lambda * loadings_vector

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Compute error
    dError <- 2 * (implied_R - R) / sum(lower_triangle)
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector( # (2x leads to fewer iterations)
        t(crossprod(2 * loadings_matrix, I %*% dError %*% I)) + (differences > 0) +
          l2_gradient
      )
    )

  }else{ # Without constraints, send it

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Compute error
    dError <- 2 * (implied_R - R) / sum(lower_triangle)
    I <- attributes(implied_R)$calculations$I

    # Return gradient (2x leads to fewer iterations)
    return(
      as.vector(
        t(crossprod(2 * attributes(implied_R)$calculations$loadings, I %*% dError %*% I)) +
          l2_gradient
      )
    )

  }

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
# Updated 20.06.2025
logLik_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    lambda, ...
)
{

  # Set regularization penalty
  penalty <- lambda * sum(loadings_vector^2)

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Check for positive definite
    return(
      -log_likelihood(n, v, implied_R, R, type = "zero") + sum(differences > 0) + penalty
    )

  }else{

    # Obtain implied correlation matrix
    implied_R <- obtain_implied(loadings_vector, rows)

    # Without constraints, send it
    return(
      swiftelse(
        is_positive_definite(implied_R),
        -log_likelihood(n, v, implied_R, R, type = "zero") + penalty,
        1e10
      )
    )

  }

}

#' @noRd
# Log-likelihood gradient
# Updated 20.06.2025
logLik_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    lambda, ...
)
{

  # Set regularization penalty
  penalty <- lambda * 2 * loadings_vector

  # Check for constraint
  if(constrained){

    # Get implied and inverse R
    implied_R <- obtain_implied(loadings_vector, rows)
    inverse_R <- solve(implied_R)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Compute error (remove negative to not have to add back later)
    dError <- ((n/2) * (inverse_R - inverse_R %*% R %*% inverse_R))
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector(t(crossprod(loadings_matrix, I %*% dError %*% I)) + (differences > 0)) +
        penalty
    )

  }else{ # Without constraints, send it

    # Get implied and inverse R
    implied_R <- obtain_implied(loadings_vector, rows)
    inverse_R <- solve(implied_R)

    # Compute error
    dError <- ((n/2) * (inverse_R - inverse_R %*% R %*% inverse_R))
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector(t(crossprod(attributes(implied_R)$calculations$loadings, I %*% dError %*% I)))
      + penalty
    )

  }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM optimization ----
# Updated 20.06.2025
egm_optimize <- function(
    loadings_vector, zeros,
    R, loading_structure, rows, n, v,
    constrained, lower_triangle, lambda,
    opt, iterations = 10000, ...
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
        R = R, loading_structure = loading_structure,
        rows = rows, n = n, v = v, constrained = constrained,
        lower_triangle = lower_triangle, lambda = lambda,
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
      R = R, loading_structure = loading_structure,
      rows = rows, n = n, v = v, constrained = constrained,
      lower_triangle = lower_triangle, lambda = lambda,
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
    R, loading_structure, rows, n, v,
    constrained, lower_triangle, norm, opt,
    iterations = 10000, ...
)
{

  # Optimize over loadings
  result <- try(
    egm_optimize(
      loadings_vector = loadings_vector, zeros = zeros,
      R = R, loading_structure = loading_structure,
      rows = rows, n = n, v = v, constrained = constrained,
      lower_triangle = lower_triangle,
      lambda = lambda, opt = opt,
      iterations = iterations, ...
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
# Updated 13.07.2025
loadings_optimization <- function(
    iter = 10, loadings_vector, zeros, R, loading_structure,
    communities, data_dimensions, lower_triangle, opt
){

  # Initial best parameters
  best_result <- list(par = loadings_vector)
  best_eigenvalue <- Inf; negative_eigenvalue <- -Inf

  # Set lambda minimum and maximum
  loadings_lambda_min <- 0
  loadings_lambda_max <- 1

  # Seek out positive eigenvalue
  for(i in seq_len(10)){

    # Optimize for best quality solution (quick passes)
    lambda <- optimize(
      f = hessian_optimize, interval = c(loadings_lambda_min, loadings_lambda_max),
      loadings_vector = best_result$par, zeros = zeros,
      R = R, loading_structure = loading_structure,
      rows = communities, n = data_dimensions[1],
      v = data_dimensions[2], constrained = FALSE,
      lower_triangle = lower_triangle,
      opt = opt, tol = 1e-03, iterations = 100
    )

    # Optimize over loadings
    result <- try(
      egm_optimize(
        loadings_vector = best_result$par, zeros = zeros,
        R = R, loading_structure = loading_structure,
        rows = communities, n = data_dimensions[1],
        v = data_dimensions[2], constrained = FALSE,
        lower_triangle = lower_triangle,
        opt = opt, lambda = lambda$minimum
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
# Updated 04.04.2025
srmr_network_cost <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  K <- solve(-network)
  I <- diag(sqrt(1 / diag(K)))

  # Return SRMR
  return(srmr(I %*% K %*% I, R))

}

#' @noRd
# SRMR network gradient
# Updated 24.04.2025
srmr_network_gradient <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  K <- solve(-network)
  I <- diag(sqrt(1 / diag(K)))

  # Compute error
  dError <- 2 * (I %*% K %*% I - R) / sum(lower_triangle)

  # Return gradient
  return((K %*% I %*% dError %*% I %*% K)[lower_triangle][zeros])

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Gaussian log-likelihood ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Log-likelihood network cost
# Updated 25.04.2025
logLik_network_cost <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  K <- solve(-network)
  I <- diag(sqrt(1 / diag(K)))

  # Return log-likelihood
  return(-log_likelihood(n, v, I %*% K %*% I, R, type = "zero"))

}

#' @noRd
# Log-likelihood network gradient
# Updated 24.04.2025
logLik_network_gradient <- function(network_vector, R, n, v, lower_triangle, zeros, ...)
{

  # Initialize network
  network <- matrix(0, nrow = v, ncol = v)
  network[lower_triangle][zeros] <- network_vector
  network <- network + t(network)

  # Convert to correlation matrix
  diag(network) <- -1
  K <- solve(-network)
  I <- diag(sqrt(1 / diag(K)))
  inverse_R <- solve(I %*% K %*% I)

  # Compute error
  dError <- ((n/2) * (inverse_R - inverse_R %*% R %*% inverse_R))

  # Return gradient
  return((K %*% I %*% dError %*% I %*% K)[lower_triangle][zeros])

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Function ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM network optimization
# Updated 13.06.2025
egm_network_optimize <- function(
    network_vector, network_length,
    R, n, v, lower_triangle,
    zeros, opt, ...
)
{

  return(
    silent_call(
      nlminb(
        start = network_vector,
        objective = switch(
          opt,
          "loglik" = logLik_network_cost,
          "srmr" = srmr_network_cost
        ),
        gradient = switch(
          opt,
          "loglik" = logLik_network_gradient,
          "srmr" = srmr_network_gradient
        ),
        R = R, n = n, v = v, lower_triangle = lower_triangle, zeros = zeros,
        lower = rep(-1, network_length), upper = rep(1, network_length),
        control = list(
          eval.max = 1000, iter.max = 1000,
          step.min = .Machine$double.eps, step.max = 1,
          abs.tol = .Machine$double.eps, rel.tol = .Machine$double.eps,
          x.tol = .Machine$double.eps, xf.tol = .Machine$double.eps,
          sing.tol = .Machine$double.eps
        )
      )
    )
  )

}
