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
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"small"} --- about 0.20
#'
#' \item \code{"moderate"} --- about 0.35
#'
#' \item \code{"large"} --- about 0.50
#'
#' }
#'
#' @param cross.loadings Numeric (length = 1).
#' Standard deviation of a normal distribution with a mean of zero.
#' Defaults to \code{1}
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
#' @param sample.size Numeric (length = 1).
#' Number of observations to generate
#'
#' @param max.iterations Numeric (length = 1).
#' Number of iterations to attempt to get convergence before erroring out.
#' Defaults to \code{100}
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
simEGM <- function(
    communities, variables,
    loadings = c("small", "moderate", "large"), cross.loadings = 0.01,
    correlations = c("none", "small", "moderate", "large", "very large"),
    sample.size, max.iterations = 100
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
  variables <- table(membership)

  # Get total variables
  total_variables <- sum(variables)

  # Get the start and end of the variables
  end <- cumsum(variables)
  start <- (end + 1) - variables

  # Generate loadings
  loadings_matrix <- matrix(
    0, nrow = total_variables, ncol = communities,
    dimnames = list(
      paste0("V", format_integer(seq_len(total_variables), digits(total_variables) - 1)),
      paste0("C", format_integer(community_sequence, digits(communities) - 1))
    )
  )

  # Determine loading ranges
  loading_range <- switch(
    loadings,
    "small" = 0.15,
    "moderate" = 0.20,
    "large" = 0.25
  )

  # Determine correlation ranges
  correlation_range <- swiftelse(
    correlations == "none", 0.00,
    switch(
      correlations,
      "small" = 0.015,
      "moderate" = 0.03,
      "large" = 0.06,
      "very large" = 0.09
    )
  )

  # Initialize R
  R <- NA

  # Count iterations
  count <- 0

  # Set cross-loading sparsity parameter
  sparsity <- 1 / log10(sample.size)

  # Iterate until positive definite
  while(anyNA(R) || any(matrix_eigenvalues(R) < 0)){

    # Increase count
    count <- count + 1

    # Populate loadings
    for(i in community_sequence){

      # Populate assigned loadings
      loadings_matrix[start[i]:end[i], i] <- loading_range + rnorm_ziggurat(variables[i]) * 0.02
      # rnorm(variables[i], mean = loading_range, sd = 0.01)

      # Get indices
      indices <- loadings_matrix[start[i]:end[i], -i]

      # Number of variables
      index_length <- length(indices)

      # Add correlations on cross-loadings
      loadings_matrix[start[i]:end[i], -i] <- correlation_range + rnorm_ziggurat(index_length) * 0.01
      # rnorm(index_length, mean = correlation_range, sd = 0.01)

      # Generate cross-loadings
      cross_loading <- sample( # Probability of cross-loading being included is set by `log10(sample.size)`
        c(0, 1), size = index_length, replace = TRUE, prob = c(sparsity, 1 - sparsity)
      ) * (0.00 + rnorm_ziggurat(index_length) * cross.loadings)
      # rnorm(index_length, mean = 0.00, sd = cross.loadings)

      # Populate cross-loading
      loadings_matrix[start[i]:end[i], -i] <- loadings_matrix[start[i]:end[i], -i] + cross_loading

    }

    # Compute A
    A <- -tcrossprod(loadings_matrix)

    # Obtain uniqueness
    uniqueness <- 1 - rowSums(loadings_matrix^2)

    # Add uniqueness to diagonal
    diag(A) <- uniqueness

    # Obtain anti-image
    anti_image <- A %*% pcor2inv(A) %*% A

    # Obtain partial correlation matrix
    P <- -anti_image / tcrossprod(sqrt(uniqueness))
    diag(P) <- 0

    # Convert to correlation matrix
    R <- silent_call(pcor2cor(P))

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
      population_partial_correlatiion = P,
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
  length_error(variables, c(1, length(communities)), "simEGM")
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
