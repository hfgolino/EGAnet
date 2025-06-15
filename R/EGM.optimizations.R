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
# Updated 09.06.2025
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
# Updated 09.06.2025
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
# Updated 02.06.2025
logLik_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, lambda, ...
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
    return(
      -log_likelihood(n, v, implied_R, R, type = "zero") + sum(differences > 0) + l2_cost
    )

  }else{

    # Without constraints, send it
    return(
      -log_likelihood(n, v, obtain_implied(loadings_vector, rows), R, type = "zero") +
        l2_cost
    )

  }

}

#' @noRd
# Log-likelihood gradient
# Updated 02.06.2025
logLik_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, lambda, ...
)
{

  # Set l2 gradient
  l2_gradient <- 2 * lambda * loadings_vector

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
        l2_gradient
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
      as.vector(t(crossprod(attributes(implied_R)$calculations$loadings, I %*% dError %*% I))) +
        l2_gradient
    )

  }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM optimization ----
# Updated 10.06.2025
egm_optimize <- function(
    loadings_vector, zeros,
    R, loading_structure, rows, n, v,
    constrained, lower_triangle, lambda, opt, ...
)
{

  # Get optimization functions
  cost <- switch(opt, "loglik" = logLik_cost, "srmr" = srmr_cost)
  gradient <- switch(opt, "loglik" = logLik_gradient, "srmr" = srmr_gradient)

  # Set bounds (based on revised network loadings paper)
  bounds <- zeros * 0.70

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
          eval.max = 1000, iter.max = 1000,
          step.min = 1e-12, step.max = 0.01
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
# Updated 13.06.2025
hessian_optimize <- function(
    lambda, loadings_vector, zeros,
    R, loading_structure, rows, n, v,
    constrained, lower_triangle, opt, ...
)
{

  # Optimize over loadings
  result <- try(
    egm_optimize(
      loadings_vector = loadings_vector, zeros = zeros,
      R = R, loading_structure = loading_structure,
      rows = rows, n = n, v = v, constrained = constrained,
      lower_triangle = lower_triangle,
      lambda = exp(lambda), opt = opt
    ), silent = TRUE
  )

  # Get error flag
  error_flag <- is(result, "try-error")

  # Check for error
  return(
    swiftelse(
      is(result, "try-error"), -1e10,
      round(min(matrix_eigenvalues(result$hessian)), 3)
    )
  )

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
          step.min = 1e-10, step.max = 1e-06
        )
      )
    )
  )

}