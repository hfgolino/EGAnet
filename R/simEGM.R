#' Simulate data following a Exploratory Graph Model
#'
#' @description Function to simulate data based on the Exploratory Graph Model
#'
#' @param communities Numeric (length = 1).
#' Number of communities to generate
#'
#' @param variables Numeric vector (length = 1 or \code{communities}).
#' Number of variables per community
#'
#' @param loadings Character (length = 1).
#' Magnitude of the assigned network loadings.
#' Available options (revised network loadings with \code{scaling = 2} in parentheses):
#'
#' \itemize{
#'
#' \item \code{"small"} --- 0.15 (approximates 0.20)
#'
#' \item \code{"moderate"} --- 0.20 (approximates 0.35)
#'
#' \item \code{"large"} --- 0.25 (approximates 0.50)
#'
#' }
#'
#' Values provided are for \code{scaling = 2/3} in \code{\link[EGAnet]{net.loads}}.
#' Uses \code{runif(n, min = value - 0.025, max = value + 0.025)} for some jitter in the loadings
#'
#' @param cross.loadings Numeric (length = 1).
#' Standard deviation of a normal distribution with a mean of zero.
#' Defaults to \code{0.01}
#'
#' @param correlations Character (length = 1).
#' Magnitude of the community correlations.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"none"} --- 0.00
#'
#' \item \code{"small"} --- 0.15 (adds about 0.015 to the cross-loadings)
#'
#' \item \code{"moderate"} --- 0.30 (adds about 0.03 to the cross-loadings)
#'
#' \item \code{"large"} --- 0.50 (adds about 0.06 to the cross-loadings)
#'
#' \item \code{"very large"} --- 0.70 (adds about 0.09 to the cross-loadings)
#'
#' }
#'
#' Uses \code{rnorm(n, mean = value, sd = 0.01)} for some jitter in the correlations
#'
#' @param sample.size Numeric (length = 1).
#' Number of observations to generate
#'
#' @param max.iterations Numeric (length = 1).
#' Number of iterations to attempt to get convergence before erroring out.
#' Defaults to \code{1000}
#'
#' @examples
#' \dontrun{
#' # Estimate EGM
#' EGM_data <- simEGM(
#'   communities = 2, variables = 6,
#'   loadings = "moderate", correlations = "moderate",
#'   sample.size = 1000
#' )}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#
# Simulate EGM ----
# Updated 29.09.2024
simEGM <- function(
    communities, variables,
    loadings = c("small", "moderate", "large"), cross.loadings = 0.01,
    correlations = c("none", "small", "moderate", "large", "very large"),
    sample.size, max.iterations = 1000
)
{

  # Check for missing arguments (argument, default, function)
  loadings <- set_default(loadings, "moderate", simEGM)
  correlations <- set_default(correlations, "moderate", simEGM)

  # Argument errors (return data in case of tibble)
  simEGM_errors(communities, variables, cross.loadings, sample.size, max.iterations)

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

  # Determine loading ranges
  loading_range <- switch(
    loadings,
    "small" = 0.40,
    "moderate" = 0.55,
    "large" = 0.70
  )

  # Determine correlation ranges
  correlation_range <- (switch(
    correlations,
    "none" = 0.00,
    "small" = 0.15,
    "moderate" = 0.30,
    "large" = 0.40,
    "very large" = 0.50
  ) + switch(
    loadings,
    "small" = -0.10,
    "moderate" = 0.00,
    "large" = 0.05
  )) / sqrt(log(variables^2))

  # Ensure zero is minimum
  correlation_range <- swiftelse(correlation_range < 0, 0, correlation_range)

  # Sparsity adjustment based on correlations
  sparsity_adjustment <- switch(
    correlations,
    "none" = log10(sample.size),
    "small" = 1,
    "moderate" = 1/3,
    "large" = 0.25,
    "very large" = 0
  )


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

    # Store community sums
    community_sums <- numeric(communities)

    # Populate loadings
    for(i in community_sequence){

      # Populate assigned loadings
      loadings_matrix[start[i]:end[i], i] <- runif_xoshiro(
        variables[i], min = loading_range - 0.05, max = loading_range + 0.05
      )
      # loading_range + rnorm_ziggurat(variables[i]) * 0.02
      # rnorm(variables[i], mean = loading_range, sd = 0.01)

      # Obtain the sum
      community_sums[i] <- sum(abs(loadings_matrix[start[i]:end[i], i]))

      # Get indices
      indices <- loadings_matrix[start[i]:end[i], -i]

      # Number of variables
      index_length <- length(indices)

      # Add correlations on cross-loadings
      loadings_matrix[start[i]:end[i], -i] <- runif_xoshiro(
        variables[i], min = correlation_range[i] - 0.05, max = correlation_range[i] + 0.05
      )
      # correlation_range[i] + rnorm_ziggurat(index_length) * 0.01
      # rnorm(index_length, mean = correlation_range, sd = 0.01)

      # Populate cross-loading
      loadings_matrix[start[i]:end[i], -i] <- loadings_matrix[start[i]:end[i], -i] +
        (0.00 + rnorm_ziggurat(index_length) * cross.loadings) *
        sample( # Probability of cross-loading being included is set by `log10(sample.size)`
          c(0, 1), size = index_length, replace = TRUE, prob = c(sparsity, 1 - sparsity)
        )
      # rnorm(index_length, mean = 0.00, sd = cross.loadings)

      # Set sparsity in cross-loadings
      loadings_matrix[start[i]:end[i], -i] <- loadings_matrix[start[i]:end[i], -i] *
        sample( # Probability of cross-loading being included is set by `log10(sample.size)`
          c(0, 1), size = index_length, replace = TRUE, prob = c(sparsity, 1 - sparsity)
        )

    }

    # Adjust loadings matrix for number of variables in each community
    loadings_matrix <- t(t(loadings_matrix) / (community_sums^(1 / log(2 * variables))))

    # Obtain partial correlations from loadings
    P <- silent_call(nload2pcor(loadings_matrix))

    # Obtain zero-order correlations from loadings
    R <- silent_call(nload2cor(loadings_matrix))

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
  )

  # Return results
  return(
    list(
      data = data %*% cholesky,
      population_partial_correlation = P,
      population_correlation = R,
      parameters = list(
        loadings = loadings_matrix,
        correlations = correlations,
        membership = membership,
        iterations = count
      )
    )
  )

}

#' @noRd
# Errors ----
# Updated 25.09.2024
simEGM_errors <- function(communities, variables, cross.loadings, sample.size, max.iterations)
{

  # 'communities'
  length_error(communities, 1, "simEGM")
  typeof_error(communities, "numeric", "simEGM")
  range_error(communities, c(1, Inf), "simEGM")

  # 'variables'
  length_error(variables, c(1, communities), "simEGM")
  typeof_error(variables, "numeric", "simEGM")
  range_error(variables, c(2, Inf), "simEGM")

  # 'cross.loadings'
  length_error(cross.loadings, 1, "simEGM")
  typeof_error(cross.loadings, "numeric", "simEGM")
  range_error(cross.loadings, c(0, 1), "simEGM")

  # 'sample.size'
  length_error(sample.size, 1, "simEGM")
  typeof_error(sample.size, "numeric", "simEGM")
  range_error(sample.size, c(1, Inf), "simEGM")

  # 'max.iterations'
  length_error(max.iterations, 1, "simEGM")
  typeof_error(max.iterations, "numeric", "simEGM")
  range_error(max.iterations, c(1, Inf), "simEGM")

}

#' @noRd
# Network loadings to partial correlations ----
# Updated 29.09.2024
nload2pcor <- function(loadings)
{

  # Obtain uniqueness
  uniqueness <- 1 - rowSums(loadings^2)

  # Compute partial correlation
  P <- tcrossprod(loadings) * tcrossprod(sqrt(uniqueness))
  diag(P) <- 0

  # # Obtain uniqueness
  # uniqueness <- 1 - rowSums(loadings^2)
  #
  # # Compute A
  # A <- -tcrossprod(loadings) * tcrossprod(sqrt(uniqueness))
  #
  # # Add uniqueness to diagonal
  # diag(A) <- uniqueness
  #
  # # Obtain anti-image
  # anti_image <- A %*% pcor2inv(A) %*% A
  #
  # # Obtain partial correlation matrix
  # P <- -anti_image / tcrossprod(sqrt(uniqueness))
  # diag(P) <- 0

  # Return
  return(P)

}

#' @noRd
# Network loadings to correlations ----
# Updated 29.09.2024
nload2cor <- function(loadings)
{

  # Obtain uniqueness
  uniqueness <- 1 - rowSums(loadings^2)

  # Compute partial correlation
  P <- tcrossprod(loadings) * tcrossprod(sqrt(uniqueness))
  diag(P) <- sqrt(1 - uniqueness)

  # Return correlation
  return(cov2cor(P))

  # # Obtain uniqueness
  # uniqueness <- 1 - rowSums(loadings^2)
  #
  # # Compute A
  # A <- -tcrossprod(loadings) * tcrossprod(sqrt(uniqueness))
  #
  # # Add uniqueness to diagonal
  # diag(A) <- uniqueness
  #
  # # Obtain anti-image
  # anti_image <- A %*% pcor2inv(A) %*% A
  #
  # # Obtain partial correlation matrix
  # P <- -anti_image / tcrossprod(sqrt(uniqueness))
  # diag(P) <- sqrt(1 - uniqueness)
  #
  # # Return correlation
  # return(cov2cor(solve(-P)))

}
