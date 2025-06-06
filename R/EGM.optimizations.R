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
# Updated 29.05.2025
srmr_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    ...
)
{

  # Set l2 cost (lambda = 0.01)
  l2_cost <- 0.01 * sum(loadings_vector^2)

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
# Updated 29.05.2025
srmr_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle,
    ...
)
{

  # Set l2 gradient (lambda = 0.01)
  l2_gradient <- 0.02 * loadings_vector

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
# Updated 02.06.2025
egm_optimize <- function(
    loadings_vector, loadings_length,
    zeros, R, loading_structure, rows, n, v,
    constrained, lower_triangle, lambda, opt, ...
)
{

  return(
    silent_call(
      nlminb(
        start = loadings_vector,
        objective = switch(
          opt,
          "loglik" = logLik_cost,
          "srmr" = srmr_cost
        ),
        gradient = switch(
          opt,
          "loglik" = logLik_gradient,
          "srmr" = srmr_gradient
        ),
        R = R, loading_structure = loading_structure,
        rows = rows, n = n, v = v, constrained = constrained,
        lower_triangle = lower_triangle, lambda = lambda,
        lower = rep(-1, loadings_length) * zeros, upper = zeros,
        control = list(
          eval.max = 10000, iter.max = 10000,
          step.min = 1e-12, step.max = 0.01,
          rel.tol = 1e-12, abs.tol = 1e-12
        )
      )
    )
  )

}

# # @noRd
# # lambda search ----
# # Updated 02.06.2025
# random_start <- function(
#     loadings, communities, data_dimensions, empirical_R, opt
# )
# {
#
#   # Get loading dimensions
#   dimensions <- dim(loadings)
#   dimension_names <- dimnames(loadings)
#
#   # Obtain loadings vector
#   loadings_vector <- as.vector(loadings)
#
#   # Get length and set zeros
#   loadings_length <- length(loadings_vector)
#
#   # Update loadings vector
#   loadings_vector <- loadings_vector * runif_xoshiro(
#     loadings_length, min = 1e-04, max = 1e-02
#   )
#
#   # Allow zeros to be estimated
#   zeros <- rep(1, loadings_length)
#
#   # Set up loading structure
#   loading_structure <- matrix(
#     TRUE, nrow = dimensions[1], ncol = dimensions[2],
#     dimnames = list(dimension_names[[1]], dimension_names[[2]])
#   )
#
#   # Optimize over loadings
#   result <- try(
#     egm_optimize(
#       loadings_vector = loadings_vector,
#       loadings_length = loadings_length,
#       zeros = zeros, R = empirical_R,
#       loading_structure = loading_structure,
#       rows = communities, n = data_dimensions[1],
#       v = data_dimensions[2], constrained = FALSE,
#       lower_triangle = lower.tri(empirical_R),
#       lambda = lambda, opt = opt
#     ), silent = TRUE
#   )
#
#   # Return values
#   if(is(result, "try-error")){
#     return(list(loadings = NULL, fit = NA, convergence = 1))
#   }else{
#
#     # Check Hessian
#     hessian <- try(
#       optimHess(
#         par = result$par,
#         fn = switch(
#           opt,
#           "loglik" = logLik_cost,
#           "srmr" = srmr_cost
#         ),
#         gr = switch(
#           opt,
#           "loglik" = logLik_gradient,
#           "srmr" = srmr_gradient
#         ),
#         loadings_length = loadings_length,
#         zeros = zeros, R = empirical_R,
#         loading_structure = loading_structure,
#         rows = communities, n = data_dimensions[1],
#         v = data_dimensions[2], constrained = FALSE,
#         lower_triangle = lower.tri(empirical_R)
#       ), silent = TRUE
#     )
#
#     # Get minimum eigenvalue
#     min_eigenvalue <- try(min(matrix_eigenvalues(hessian)), silent = TRUE)
#     min_eigenvalue <- swiftelse(is(min_eigenvalue, "try-error"), Inf, min_eigenvalue)
#
#     # Get condition number
#     condition_number <- try(kappa(hessian), silent = TRUE)
#     condition_number <- swiftelse(
#       is(condition_number, "try-error"), Inf, condition_number
#     )
#
#     # Format loadings
#     loadings <- matrix(
#       result$par,
#       nrow = data_dimensions[2], ncol = communities,
#       dimnames = dimnames(loadings)
#     )
#
#     # Return value
#     return(
#       list(
#         loadings = loadings, fit = result$objective,
#         convergence = 0, # accept false convergences (due to ridge penalty)
#         min_eigenvalue_sign = sign(min_eigenvalue),
#         min_eigenvalue = min_eigenvalue,
#         condition_number = condition_number
#       )
#     )
#
#   }
#
# }

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
# Updated 02.06.2025
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
          eval.max = 10000, iter.max = 10000,
          step.min = 1e-10, step.max = 1e-06
        )
      )
    )
  )

}