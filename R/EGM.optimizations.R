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
# Updated 04.04.2025
srmr_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Check for positive definite
    return(srmr(R, implied_R) + sum(differences > 0))

  }else{ # Without constraints, send it
    return(srmr(R, obtain_implied(loadings_vector, rows)))
  }

}

#' @noRd
# SRMR gradient
# Updated 04.04.2025
srmr_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Compute error
    error <- (implied_R - R)[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- 2 * error / length(error)
    dError <- dError + t(dError)
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector( # (2x leads to fewer iterations)
        t(crossprod(2 * loadings_matrix, I %*% tcrossprod(dError, I))) + (differences > 0)
      )
    )

  }else{ # Without constraints, send it

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Compute error
    error <- (implied_R - R)[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- 2 * error / length(error)
    dError <- dError + t(dError)
    I <- attributes(implied_R)$calculations$I

    # Return gradient (2x leads to fewer iterations)
    return(
      as.vector(
        t(crossprod(2 * attributes(implied_R)$calculations$loadings, I %*% tcrossprod(dError, I)))
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
# Updated 04.04.2025
logLik_cost <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Get implied R
    implied_R <- obtain_implied(loadings_vector, rows)

    # Assemble loading matrix
    loadings_matrix <- attributes(implied_R)$calculations$loadings

    # Obtain differences
    differences <- abs(loadings_matrix) - loadings_matrix[loading_structure]

    # Check for positive definite
    return(-log_likelihood(n, v, implied_R, R, type = "zero") + sum(differences > 0))

  }else{ # Without constraints, send it
    return(-log_likelihood(n, v, obtain_implied(loadings_vector, rows), R, type = "zero"))
  }

}

#' @noRd
# Log-likelihood gradient
# Updated 04.04.2025
logLik_gradient <- function(
    loadings_vector, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

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
    error <- ((n/2) * (inverse_R - inverse_R %*% R %*% inverse_R))[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- error
    dError <- dError + t(dError)
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector(t(crossprod(loadings_matrix, I %*% tcrossprod(dError, I))) + (differences > 0))
    )

  }else{ # Without constraints, send it

    # Get implied and inverse R
    implied_R <- obtain_implied(loadings_vector, rows)
    inverse_R <- solve(implied_R)

    # Compute error
    error <- ((n/2) * (inverse_R - inverse_R %*% R %*% inverse_R))[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- error
    dError <- dError + t(dError)
    I <- attributes(implied_R)$calculations$I

    # Return gradient
    return(
      as.vector(t(crossprod(attributes(implied_R)$calculations$loadings, I %*% tcrossprod(dError, I))))
    )

  }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Optimization Function ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# EGM optimization
# Updated 11.04.2025
egm_optimize <- function(
    loadings_vector, loadings_length,
    zeros, R, loading_structure, rows, n, v,
    constrained, lower_triangle, opt, ...
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
        rows = rows, n = n, v = v,
        constrained = constrained, lower_triangle = lower_triangle,
        lower = rep(-1, loadings_length) * zeros, upper = zeros,
        control = list(eval.max = 10000, iter.max = 10000)
      )
    )
  )

}
