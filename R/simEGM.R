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
#' For reference, small (0.20), moderate (0.35), and large (0.50)
#'
#' Uses \code{runif(n, min = value - 0.075, max = value + 0.075)} for some jitter in the loadings
#'
#' @param cross.loadings Numeric (length = 1).
#' Standard deviation of a normal distribution with a mean of zero (\code{n, mean = 0, sd = value}).
#' Defaults to \code{0.01}.
#' Not recommended to change too drastically (small increments such as \code{0.01} work best)
#'
#' @param correlations Numeric (length = 1).
#' Magnitude of the community correlations
#'
#' @param sample.size Numeric (length = 1).
#' Number of observations to generate
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
#' @export
#'
# Simulate EGM ----
# Updated 13.04.2025
simEGM <- function(
    communities, variables,
    loadings, cross.loadings = 0.01,
    correlations, sample.size,
    max.iterations = 1000
){

  # Argument errors
  simEGM_errors(
    communities, variables, loadings, cross.loadings,
    correlations, sample.size, max.iterations
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

  # Set names
  names(membership) <- node_names <- paste0(
    "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
  )

  # Initialize structure matrix
  loading_structure <- matrix(
    0, nrow = total_variables, ncol = communities,
    dimnames = list(
      node_names,
      format_integer(community_sequence, digits(communities) - 1) # community names
    )
  )

  # Fill structure
  for(i in seq_len(total_variables)){
    loading_structure[i, membership[i]] <- 1
  }

  # Count iterations
  count <- 0

  # Initialize checks
  PD_check <- CC_check <- FALSE

  # Ensure proper matrix
  while(!PD_check){

    # Increase count
    count <- count + 1

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

    # Generate within-community loadings
    for(i in community_sequence){

      # Get block index
      block_index <- loading_structure[,i] == 1

      # Set assigned loadings
      loading_structure[block_index, i] <- runif_xoshiro(
        sum(block_index), min = loadings - 0.075, max = loadings + 0.075
      )

      # Get off-diagonal indices
      off_index <- loading_structure[,i] == 0

      # Set correlations
      loading_structure[off_index, i] <- 0.25 * correlations /
        ((1 - loading_structure[block_index, i]) * sqrt(log(total_variables)))

      # Add cross-loadings
      loading_structure[off_index, i] <- loading_structure[off_index, i] +
        rnorm_ziggurat(sum(off_index)) * cross.loadings

    }

    # Obtain population correlation matrix
    R <- nload2cor(loading_structure)

    # Check for positive definite
    PD_check <- is_positive_definite(R)

  }

  # Perform Cholesky decomposition
  cholesky <- chol(R)

  # Generate data
  data <- MASS_mvrnorm_quick(
    p = total_variables, np = total_variables * sample.size,
    coV = diag(total_variables)
  ) %*% cholesky

  # Set variable names
  colnames(data) <- node_names

  # Obtain precision, partial correlations, and adjacency matrices
  K <- solve(R)
  P <- -cov2cor(K); diag(P) <- 0
  adjacency <- beta_min(P, membership, K, total_variables, sample.size)

  # Return results
  return(
    list(
      data = data,
      population_correlation = R,
      population_precision = K,
      population_partial_correlation = P,
      parameters = list(
        adjacency = adjacency,
        network = P * adjacency,
        loadings = loading_structure,
        correlations = correlations,
        membership = membership,
        iterations = count
      )
    )
  )

}

#' @noRd
# Errors ----
# Updated 06.04.2025
simEGM_errors <- function(
    communities, variables, loadings, cross.loadings,
    correlations, sample.size, max.iterations
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

  # 'max.iterations'
  length_error(max.iterations, 1, "simEGM")
  typeof_error(max.iterations, "numeric", "simEGM")
  range_error(max.iterations, c(1, Inf), "simEGM")

}

#' @noRd
# beta-min criterion ----
# Updated 13.04.2025
beta_min <- function(P, membership = NULL, K, total_variables, sample_size)
{

  # Obtain modularity-adjusted beta-min criterion
  minimum <- swiftelse(
    is.null(membership), 1, modularity(P, membership)
  ) * sqrt(log(total_variables) / sample_size)

  # Obtain inverse variances
  inverse_variances <- diag(K)

  # Obtain betas
  beta <- P * sqrt(outer(inverse_variances, inverse_variances, FUN = "/"))

  # Set adjacency matrix
  adjacency <- abs(beta) > minimum

  # Attach minimum
  attr(adjacency, "beta.min") <- minimum

  # Return adjacency matrix
  return(adjacency)

}

