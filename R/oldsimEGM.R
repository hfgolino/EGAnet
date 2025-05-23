#' Simulate data following a Exploratory Graph Model (\code{\link[EGAnet]{EGM}})
#'
#' @description Function to simulate data based on \code{\link[EGAnet]{EGM}}
#'
#' @param communities Numeric (length = 1).
#' Number of communities to generate
#'
#' @param variables Numeric vector (length = 1 or \code{communities}).
#' Number of variables per community
#'
#' @param loadings Numeric (length = 1).
#' Magnitude of the assigned network loadings.
#' Uses the same magnitude as factors loadings
#'
#' Uses \code{runif(n, min = value - 0.025, max = value + 0.025)} for some jitter in the loadings
#'
#' @param cross.loadings Numeric (length = 1).
#' Standard deviation of a normal distribution with a mean of zero (\code{n, mean = 0, sd = value}).
#' Defaults to \code{0.01}
#'
#' @param correlations Numeric (length = 1).
#' Magnitude of the community correlations
#'
#' Uses \code{runif(n, min = value - 0.015, max = value + 0.015)}
#' for some jitter in the correlations
#'
#' @param sample.size Numeric (length = 1).
#' Number of observations to generate
#'
#' @param p.in Numeric (length = 1).
#' Sets the probability of retaining an edge \emph{within} communities.
#' Single values are applied to all communities.
#' Defaults to \code{0.95}
#'
#' @param p.out Numeric (length = 1 or \code{communities}).
#' Sets the probability of retaining an edge \emph{between} communities.
#' Single values are applied to all communities.
#' Defaults to \code{0.80}
#'
#' @param max.iterations Numeric (length = 1).
#' Number of iterations to attempt to get convergence before erroring out.
#' Defaults to \code{1000}
#'
#' @examples
#' simulated <- simEGM(
#'   communities = 2, variables = 6,
#'   loadings = 0.55, # use standard factor loading sizes
#'   correlations = 0.30,
#'   sample.size = 1000
#' )
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Simulate EGM ----
# Updated 07.11.2024
simEGM <- function(
    communities, variables,
    loadings, cross.loadings = 0.01, correlations,
    sample.size, p.in = 0.95, p.out = 0.80, max.iterations = 1000
)
{

  # Check for missing arguments (argument, default, function)
  # loadings <- set_default(loadings, "moderate", simEGM)
  # correlations <- set_default(correlations, "moderate", simEGM)

  # Argument errors
  simEGM_errors(
    communities, variables, loadings, cross.loadings,
    correlations, sample.size, p.in, p.out, max.iterations
  )

  # Set up membership based on communities and variables
  membership <- swiftelse(
    length(variables) == 1,
    rep(seq_len(communities), each = variables),
    rep(seq_len(communities), times = variables)
  )

  # Set community sequence
  community_sequence <- seq_len(communities)

  # Table variables
  variables <- as.numeric(table(membership))

  # Get total variables
  total_variables <- sum(variables)

  # Get the start and end of the variables
  end <- cumsum(variables)
  start <- (end + 1) - variables

  # Derived using:
  # x <- c(0.00, 0.40, 0.55, 0.70, 1.00)
  # y <- c(0.00, 0.20, 0.35, 0.55, 1.00)
  # fit <- nls(
  #   y ~ a * x^2 + b * x + c, start = list(a = 0.1, b = 0.1, c = 0),
  #   algorithm = "port",
  #   control = nls.control(maxiter = 100000, tol = 1e-20, minFactor = 1e-20)
  # ); coef(fit)

  # Obtain loadings signs
  loadings_signs <- sign(loadings)
  loadings <- abs(loadings)

  # Determine loading ranges
  loading_range <- (
    0.797267042 * loadings^2 + 0.210094592 * loadings - 0.002682913
  ) * loadings_signs

  # Categorize loadings into small, moderate, and large
  loadings <- c("small", "moderate", "large")[
    as.numeric(cut(
      loadings,
      breaks = c(0.00, 0.475, 0.625, 1.00)
    ))
  ]

  # Derived using:
  # x <- c(0.00, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
  # y <- c(0.00, 0.1125, 0.15, 0.18, 0.24, 0.29, 0.36)
  # fit <- nls(
  #   y ~ a * x^2 + b * x + c, start = list(a = 0.1, b = 0.1, c = 0),
  #   algorithm = "port",
  #   control = nls.control(maxiter = 10000, tol = 1e-20, minFactor = 1e-20)
  # ); coef(fit)

  # Determine correlation ranges
  correlation_range <- (
    (0.091376419 * correlations^2 + 0.427398177 * correlations + 0.007358814) +
      switch(
        loadings,
        "small" = -0.05,
        "moderate" = 0.000,
        "large" = 0.05
      )
  ) / sqrt(log(variables^2))

  # Ensure zero is minimum
  correlation_range <- swiftelse(correlation_range < 0, 0, correlation_range)

  # Set up for sparsity adjustment
  sparsity_adjustment <- c(
    log10(sample.size), # none
    1, # small
    1/3, # moderate
    0.25, # large
    0 # very large
  )[
    as.numeric(cut(
      correlations,
      breaks = c(0.00, 0.10, 0.30, 0.50, 0.70, 1.00)
    ))
  ]

  # Set cross-loading sparsity parameter
  sparsity <- 1 / log10(sample.size) * sparsity_adjustment

  # Initialize R
  R <- NA

  # Count iterations
  count <- 0

  # Iterate until positive definite
  while(anyNA(R) || any(matrix_eigenvalues(R) < 0)){

    # Generate loadings
    loadings_matrix <- matrix(
      0, nrow = total_variables, ncol = communities,
      dimnames = list(
        paste0("V", format_integer(seq_len(total_variables), digits(total_variables) - 1)),
        paste0("C", format_integer(community_sequence, digits(communities) - 1))
      )
    )

    # Increase count
    count <- count + 1

    # Populate loadings
    for(i in community_sequence){

      # Populate assigned loadings
      loadings_matrix[start[i]:end[i], i] <- runif_xoshiro(
        variables[i], min = loading_range - 0.025, max = loading_range + 0.025
      )

      # Get indices
      indices <- loadings_matrix[start[i]:end[i], -i]

      # Number of variables
      index_length <- length(indices)

      # Add correlations on cross-loadings
      loadings_matrix[start[i]:end[i], -i] <- runif_xoshiro(
        variables[i], min = correlation_range[i] - 0.015, max = correlation_range[i] + 0.015
      )

      # Populate cross-loading
      loadings_matrix[start[i]:end[i], -i] <- loadings_matrix[start[i]:end[i], -i] +
        (0.00 + rnorm_ziggurat(index_length) * cross.loadings)
      # rnorm(index_length, mean = 0.00, sd = cross.loadings)

      # Set sparsity in cross-loadings
      loadings_matrix[start[i]:end[i], -i] <- loadings_matrix[start[i]:end[i], -i] *
        sample( # Probability of cross-loading being included is set by `log10(sample.size)`
          c(0, 1), size = index_length, replace = TRUE, prob = c(sparsity, 1 - sparsity)
        )

    }

    # Obtain correlation matrices
    loadings_matrix <- try(
      update_loadings(
        total_variables, communities, membership,
        p.in, p.out, loadings_matrix
      ), silent = TRUE
    )

    # Set correlations
    if(!is(loadings_matrix, "try-error")){
      R <- nload2cor(loadings_matrix)
    }

    # Check for max iterations
    if(count >= max.iterations){

      # Stop and send error
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Maximum iterations reached with no convergence. \n\n",
          "Try:\n",
          "  +  increasing 'loadings'\n",
          "  +  increasing 'sample.size'\n",
          "  +  decreasing 'correlations'"
        ),
        call = "simEGM"
      )

    }

  }

  # Perform Cholesky decomposition
  cholesky <- chol(R)

  # Generate data
  data <- MASS_mvrnorm_quick(
    p = total_variables, np = total_variables * sample.size,
    coV = diag(total_variables)
  ) %*% cholesky

  # Get population correlation
  P <- cor2pcor(R)

  # Set names
  row.names(P) <- colnames(P) <-
  row.names(R) <- colnames(R) <-
  colnames(data) <- names(membership) <-
  row.names(loadings_matrix)

  # Return results
  return(
    list(
      data = data,
      population_correlation = R,
      population_partial_correlation = P,
      parameters = list(
        loadings = loadings_matrix,
        correlations = correlations,
        p.in = p.in, p.out = p.out,
        membership = membership,
        iterations = count
      )
    )
  )

}

#' @noRd
# Errors ----
# Updated 07.11.2024
simEGM_errors <- function(
    communities, variables, loadings, cross.loadings,
    correlations, sample.size, p.in, p.out, max.iterations
)
{

  # 'communities'
  length_error(communities, 1, "simEGM")
  typeof_error(communities, "numeric", "simEGM")
  range_error(communities, c(1, Inf), "simEGM")

  # 'variables'
  length_error(variables, c(1, communities), "simEGM")
  typeof_error(variables, "numeric", "simEGM")
  range_error(variables, c(1, Inf), "simEGM")

  # 'loadings'
  length_error(loadings, c(1, communities), "simEGM")
  typeof_error(loadings, "numeric", "simEGM")
  range_error(loadings, c(-1, 1), "simEGM")

  # 'cross.loadings'
  length_error(cross.loadings, 1, "simEGM")
  typeof_error(cross.loadings, "numeric", "simEGM")
  range_error(cross.loadings, c(0, 1), "simEGM")

  # 'correlations'
  length_error(correlations, c(1, communities), "simEGM")
  typeof_error(correlations, "numeric", "simEGM")
  range_error(correlations, c(0, 1), "simEGM")

  # 'sample.size'
  length_error(sample.size, 1, "simEGM")
  typeof_error(sample.size, "numeric", "simEGM")
  range_error(sample.size, c(1, Inf), "simEGM")

  # 'p.in'
  length_error(p.in, 1, "simEGM")
  typeof_error(p.in, "numeric", "simEGM")
  range_error(p.in, c(0, 1), "simEGM")

  # 'p.out'
  length_error(p.out, 1, "simEGM")
  typeof_error(p.out, "numeric", "simEGM")
  range_error(p.out, c(0, 1), "simEGM")

  # 'max.iterations'
  length_error(max.iterations, 1, "simEGM")
  typeof_error(max.iterations, "numeric", "simEGM")
  range_error(max.iterations, c(1, Inf), "simEGM")

}

#' @noRd
# Update loadings to align with network ----
# Updated 04.11.2024
update_loadings <- function(
    total_variables, communities, membership,
    p.in, p.out, loadings_matrix
)
{

  # Set number of communities
  communities <- unique_length(membership)

  # Set up community variables
  community_variables <- lapply(
    seq_len(communities), function(community){membership == community}
  )

  # Obtain partial correlation matrix
  ## Function is in `EGM.R`
  P <- create_community_structure(
    P = nload2pcor(loadings_matrix),
    total_variables = total_variables,
    communities = communities,
    community_variables = community_variables,
    p.in = p.in, p.out = p.out
  )

  # Set lower triangle
  lower_triangle <- lower.tri(P)

  # Obtain the lower triangle
  P_lower <- P[lower_triangle]

  # Get zeros
  zeros <- P_lower != 0
  P_length <- sum(zeros)

  # Use optimize to minimize the RMSE
  result <- silent_call(
    nlminb(
      start = P_lower[zeros], objective = P_cost,
      gradient = P_gradient,
      P_lower = P_lower, zeros = zeros,
      R = nload2cor(loadings_matrix),
      lower_triangle = lower_triangle,
      total_variables = total_variables,
      lower = rep(-1, P_length),
      upper = rep(1, P_length),
      control = list(eval.max = 10000, iter.max = 10000)
    )
  )

  # Update P vector
  P_lower[zeros] <- result$par

  # Fill out matrix
  P <- matrix(0, nrow = total_variables, ncol = total_variables)
  P[lower_triangle] <- P_lower
  P <- t(P) + P

  # Set bounds
  loadings_vector <- as.vector(loadings_matrix)
  loadings_length <- length(loadings_vector)
  zeros <- loadings_vector != 0

  # Use optimize to minimize the RMSE
  result <- silent_call(
    nlminb(
      start = loadings_vector,
      objective = N_cost, gradient = N_gradient,
      P = P, zeros = zeros, total_variables = total_variables,
      lower = rep(-1, loadings_length),
      upper = rep(1, loadings_length),
      control = list(eval.max = 10000, iter.max = 10000)
    )
  )

  # Extract optimized loadings
  loadings_matrix <- matrix(
    result$par, nrow = total_variables,
    dimnames = dimnames(loadings_matrix)
  )

  # Return results
  return(loadings_matrix)

}

#' @noRd
# Partial correlation cost ----
# Updated 15.10.2024
P_cost <- function(P_nonzero, P_lower, zeros, R, total_variables, lower_triangle)
{

  # Set up P vector
  P_lower[zeros] <- P_nonzero

  # Get partial correlations
  P_matrix <- matrix(0, nrow = total_variables, ncol = total_variables)

  # Set lower triangle
  P_matrix[lower_triangle] <- P_lower

  # Transpose
  P_matrix <- t(P_matrix) + P_matrix

  # Set diagonal to negative 1
  diag(P_matrix) <- -1

  # Obtain inverse
  INV <- solve(-P_matrix)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Error
  error <- (D %*% INV %*% D - R)[lower_triangle]^2

  # Return RMSE
  return(sqrt(mean(error)))

}

#' @noRd
# Partial correlation gradient ----
# Updated 10.11.2024
P_gradient <- function(P_nonzero, P_lower, zeros, R, total_variables, lower_triangle)
{

  # Set up P vector
  P_lower[zeros] <- P_nonzero

  # Get partial correlations
  P_matrix <- matrix(0, nrow = total_variables, ncol = total_variables)

  # Set lower triangle
  P_matrix[lower_triangle] <- P_lower

  # Transpose
  P_matrix <- t(P_matrix) + P_matrix

  # Set diagonal to negative 1
  diag(P_matrix) <- -1

  # Obtain inverse
  INV <- solve(-P_matrix)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Compute error
  error <- (D %*% INV %*% D - R)

  # Return gradient
  return(2 * error[lower_triangle][zeros])

}

#' @noRd
# Loadings partial correlation cost ----
# Updated 04.11.2024
N_cost <- function(loadings_vector, P, zeros, lower_triangle, total_variables)
{

  # Set up loadings matrix
  loadings_matrix <- matrix(loadings_vector * zeros, nrow = total_variables)

  # Compute partial correlation
  implied_P <- tcrossprod(loadings_matrix)

  # Obtain interdependence
  interdependence <- sqrt(rowSums(loadings_matrix^2))

  # Set diagonal to interdependence
  diag(implied_P) <- interdependence

  # Compute matrix I
  I <- diag(sqrt(1 / interdependence))

  # Compute zero-order correlations
  R <- I %*% implied_P %*% I

  # Obtain inverse covariance matrix
  INV <- solve(R)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Compute partial correlations
  implied_P <- -D %*% INV %*% D

  # Set diagonal to zero
  diag(implied_P) <- 0

  # Error
  error <- (implied_P - P)[lower_triangle]^2

  # Return SRMR cost
  return(sqrt(mean(error)))

}

#' @noRd
# Loadings partial correlation gradient ----
# Updated 06.11.2024
N_gradient <- function(loadings_vector, P, zeros, lower_triangle, total_variables)
{

  # Set up loadings matrix
  loadings_matrix <- matrix(loadings_vector * zeros, nrow = total_variables)

  # Compute partial covariance
  implied_P <- tcrossprod(loadings_matrix)

  # Obtain interdependence
  interdependence <- sqrt(rowSums(loadings_matrix^2))

  # Set diagonal to interdependence
  diag(implied_P) <- interdependence

  # Compute matrix I
  I <- diag(sqrt(1 / interdependence))

  # Compute zero-order correlations
  R <- I %*% implied_P %*% I

  # Obtain inverse covariance matrix
  INV <- solve(R)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Compute partial correlations
  implied_P <- -D %*% INV %*% D

  # Set diagonal to zero
  diag(implied_P) <- 0

  # Compute error
  error <- (implied_P - P)[lower_triangle]

  # Derivative of error with respect to P (covariance)
  dError <- matrix(0, nrow = total_variables, ncol = total_variables)
  dError[lower_triangle] <- 2 * error / length(error)
  dError <- dError + t(dError)

  # Derivative of INV with respect to P
  # dINV <- -D %*% dError %*% D

  # Derivative with respect to R
  dR <- -INV %*% (-D %*% dError %*% D) %*% INV

  # Derivative with respect to P
  # dP <- I %*% tcrossprod(dR, I)

  # Return gradient (2x leads to fewer iterations)
  return(as.vector(t(crossprod(2 * loadings_matrix, I %*% tcrossprod(dR, I)))) * zeros)

}
