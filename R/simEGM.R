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
#' @param loadings Numeric (length = 1, \code{communities}, or
#' total variables \eqn{\times} \code{communities}).
#' Magnitude of the assigned network loadings.
#' For reference, small (0.20), moderate (0.35), and large (0.50).
#' Input can be a loading matrix but must have the dimensions:
#' total variables \eqn{\times} \code{communities}
#'
#' Uses \code{runif(n, min = value - 0.075, max = value + 0.075)} for some jitter in the loadings
#'
#' @param cross.loadings Numeric (length = 1).
#' Standard deviation of a normal distribution with a mean of zero (\code{n, mean = 0, sd = value}).
#' Defaults to \code{0.01}.
#' Not recommended to change too drastically (small increments such as \code{0.01} work best)
#'
#' @param correlations Numeric (length = 1 or
#' \code{communities} \eqn{\times} \code{communities} matrix).
#' Magnitude of the community correlations.
#' Input can be a correlations matrix but must have the dimensions:
#' \code{communities} \eqn{\times} \code{communities}
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
# Updated 28.05.2025
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

  # Initialize checks
  PD_check <- FALSE

  # Check that loadings matrix is not already supplied
  if(length(as.matrix(loadings)) == (total_variables * communities)){

    # Obtain population correlation matrix
    R <- nload2cor(loadings)

    # Check for positive definite
    PD_check <- is_positive_definite(R)

    # Set error if not positive definite
    if(!PD_check){

      # Stop and send error
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Matrix input into 'loadings' does not produce a positive definite ",
          "correlation matrix.\n\n",
          "Check your matrix using: `eigen(EGAnet:::nload2cor(loadings))$values`"
        ),
        call = "simEGM"
      )

    }

    # Initialize correlation matrix
    correlations <- matrix(0, nrow = communities, ncol = communities)

    # Loop over loadings to derive correlations
    for(i in community_sequence){

      # Get block index
      block_index <- membership == i

      # Get mean loadings for block
      mean_assigned <- 1 - mean(loadings[block_index, i])

      # Loop over other blocks
      for(j in community_sequence){

        # Except for itself
        if(i != j){

          # Estimate correlations
          correlations[j,i] <- correlations[i,j] <- 4 * mean(loadings[block_index,j]) *
            (mean_assigned * sqrt(log(total_variables)))

        }

      }

    }

    # Set diagonal to 1 (add count for return)
    diag(correlations) <- count <- 1

    # Set loading structure for return
    loading_structure <- loadings

  }else{ # Generate the matrix

    # Initialize structure matrix
    loading_structure <- matrix(
      0, nrow = total_variables, ncol = communities,
      dimnames = list(
        node_names,
        format_integer(community_sequence, digits(communities) - 1) # community names
      )
    )

    # Set up loadings
    if(!is.matrix(loadings)){

      # Ensure length of communities
      if(length(loadings) == 1){
        loadings <- rep(loadings, communities)
      }

    }

    # Set up correlations
    if(!is.matrix(correlations)){

      # Initialize correlation matrix
      correlations <- matrix(correlations, nrow = communities, ncol = communities)

      # Set diagonal to 1
      diag(correlations) <- 1

    }

    # Count iterations
    count <- 0

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
        block_index <- membership == i

        # Get number of variables
        block_variables <- sum(block_index)

        # Set assigned loadings
        loading_structure[block_index, i] <- runif_xoshiro(
          block_variables, min = loadings[i] - 0.05, max = loadings[i] + 0.05
        )

        # Loop over cross-loadings
        for(j in community_sequence){

          # Except for itself
          if(i != j){

            # Set correlations
            loading_structure[block_index, j] <- 0.25 * correlations[i,j] /
              ((1 - loading_structure[block_index, i]) * sqrt(log(total_variables)))

            # Add cross-loadings
            loading_structure[block_index, j] <- loading_structure[block_index, j] +
              rnorm_ziggurat(block_variables) * cross.loadings

          }

        }

      }

      # Obtain population correlation matrix
      R <- nload2cor(loading_structure)

      # Check for positive definite
      PD_check <- is_positive_definite(R)

    }

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
# Updated 14.04.2025
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

  # Check for variables
  total_variables <- swiftelse(
    length(variables) == 1,
    communities * communities * variables,
    communities * sum(variables)
  )

  # 'loadings'
  length_error(
    as.matrix(loadings), # set matrix for length detection
    c(1, communities, total_variables), # ensure one, communities, or matrix
    "simEGM"
  )
  typeof_error(loadings, "numeric", "simEGM")
  range_error(loadings, c(-1, 1), "simEGM")

  # 'cross.loadings'
  length_error(cross.loadings, 1, "simEGM")
  typeof_error(cross.loadings, "numeric", "simEGM")
  range_error(cross.loadings, c(0, 1), "simEGM")

  # 'correlations'
  length_error(
    as.matrix(correlations), # set matrix for length detection
    c(1, communities * communities), # ensure one value or matrix
    "simEGM"
  )
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
