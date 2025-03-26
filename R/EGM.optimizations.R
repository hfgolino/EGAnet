#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Optimization functions for EGM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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
  return(
    sqrt(mean((base[lower_triangle] - comparison[lower_triangle])^2))
  )

}

#' @noRd
# SRMR cost
# Updated 20.03.2025
srmr_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

    # Obtain assign loadings
    assign_loadings <- loadings_matrix[loading_structure]

    # Transpose loadings matrix
    loadings_matrix <- t(loadings_matrix)

    # Obtain differences
    differences <- abs(loadings_matrix) - assign_loadings

    # Obtain difference values
    difference_values <- (differences * (differences > 0))^2

    # Convert loadings to implied correlations following `nload2cor`
    # Uses raw code to avoid overhead of additional function calls

    # Compute partial correlation (decrease loadings based on constraints)
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- interdependence <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / interdependence))

    # Get implied R
    implied_R <- I %*% implied_P %*% I

    # Check for positive definite
    return(
      swiftelse(
        !anyNA(implied_R) && is_positive_definite(implied_R), # return SRMR
        srmr(R, implied_R) + sqrt(mean(difference_values)), Inf
      )
    )

  }else{ # Without constraints, send it

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, ncol = rows)

    # Compute partial correlation
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / diag(implied_P)))

    # Get implied R
    implied_R <- I %*% implied_P %*% I

    # Check for positive definite
    return(
      swiftelse(
        !anyNA(implied_R) && is_positive_definite(implied_R), # return SRMR
        srmr(R, implied_R), Inf
      )
    )

  }

}

#' @noRd
# SRMR gradient
# Updated 20.03.2025
srmr_gradient <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

    # Obtain assign loadings
    assign_loadings <- loadings_matrix[loading_structure]

    # Transpose loadings matrix
    loadings_matrix <- t(loadings_matrix)

    # Obtain differences
    differences <- abs(loadings_matrix) - assign_loadings

    # Obtain difference values
    difference_values <- (differences * (differences > 0))

    # Convert loadings to implied correlations following `nload2cor`
    # Uses raw code to avoid overhead of additional function calls

    # Compute partial correlation
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- interdependence <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / interdependence))

    # Compute error
    error <- (I %*% implied_P %*% I - R)[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- 2 * error / length(error)
    dError <- dError + t(dError)

    # Return gradient
    return(
      as.vector( # (2x leads to fewer iterations)
        t(crossprod(2 * loadings_matrix, I %*% tcrossprod(dError, I))) +
          difference_values
      ) * zeros
    )

  }else{ # Without constraints, send it

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, ncol = rows)

    # Compute partial correlation
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- interdependence <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / interdependence))

    # Compute error
    error <- (I %*% implied_P %*% I - R)[lower_triangle]

    # Derivative of error with respect to P (covariance)
    dError <- matrix(0, nrow = v, ncol = v)
    dError[lower_triangle] <- 2 * error / length(error)
    dError <- dError + t(dError)

    # Return gradient (2x leads to fewer iterations)
    return(as.vector(t(crossprod(2 * loadings_matrix, I %*% tcrossprod(dError, I)))) * zeros)

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
# Updated 20.03.2025
logLik_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, n, v,
    constrained, lower_triangle, ...
)
{

  # Check for constraint
  if(constrained){

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

    # Obtain assign loadings
    assign_loadings <- loadings_matrix[loading_structure]

    # Transpose loadings matrix
    loadings_matrix <- t(loadings_matrix)

    # Obtain differences
    differences <- abs(loadings_matrix) - assign_loadings

    # Obtain difference values
    difference_values <- (differences * (differences > 0))^2

    # Convert loadings to implied correlations following `nload2cor`
    # Uses raw code to avoid overhead of additional function calls

    # Compute partial correlation (decrease loadings based on constraints)
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- interdependence <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / interdependence))

    # Get implied R
    implied_R <- I %*% implied_P %*% I

    # Check for positive definite
    return(
      swiftelse(
        !anyNA(implied_R) && is_positive_definite(implied_R), # return log-likelihood
        -log_likelihood(n, v, implied_R, R, type = "zero") + sum(difference_values),
        Inf
      )
    )

  }else{ # Without constraints, send it

    # Assemble loading matrix
    loadings_matrix <- matrix(loadings_vector * zeros, ncol = rows)

    # Compute partial correlation
    implied_P <- tcrossprod(loadings_matrix)

    # Compute interdependence
    diag(implied_P) <- interdependence <- sqrt(diag(implied_P))

    # Compute matrix I
    I <- diag(sqrt(1 / interdependence))

    # Get implied R
    implied_R <- I %*% implied_P %*% I

    # Check for positive definite
    return(
      swiftelse(
        !anyNA(implied_P) && is_positive_definite(implied_P), # return log-likelihood
        -log_likelihood(n, v, implied_R, R, type = "zero"),
        Inf
      )
    )

  }

}


#' @noRd
# Compute log-likelihood metrics ----
# Updated 06.10.2024
likelihood <- function(n, p, R, S, loadings, type)
{

  # Get number of communities
  m <- dim(loadings)[2]

  # Log-likelihood
  loglik <- log_likelihood(n, p, R, S, type)

  # Total number of parameters
  parameters <- (p * m) + p + ((m * (m - 1)) / 2)

  # Model parameters
  model_parameters <- parameters - sum(loadings == 0)

  # Return log-likelihood
  return(
    c(
      logLik = loglik,
      AIC = -2 * loglik + 2 * model_parameters, # -2L + 2k
      BIC = -2 * loglik + model_parameters * log(n) # -2L + klog(n)
      # EBIC = -2 * loglik + model_parameters * log(n) + 2 * gamma * log(
      #   choose(parameters, model_parameters)
      # ), # -2L + klog(n) + 2 gamma log(binom(pk))
      # GFI = 1 - sum((R - S)^2) / sum(S^2)
    )
  )

}