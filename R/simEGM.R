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
# Updated 21.07.2025
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
  PD_check <- FALSE; cross_check <- TRUE

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

    # Obtain simple structure
    simple_structure <- loadings
    for(i in community_sequence){
      simple_structure[membership == i, -i] <- 0
    }

    # Compute correlations
    correlations <- community_correlations(simple_structure, loadings)

    # Set loading structure for return
    loading_structure <- loadings

  }else{ # Generate the matrix

    # Initialize structure matrix
    between_indices <- loading_structure <- matrix(
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
    while(!PD_check || cross_check){

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

            # Set cross-loading probability
            cross_probability <- abs(correlations[i,j])^(1/3)

            # Set zero cross-loading indices
            between_indices[block_index, j] <- c(
              correlations[i,j] != 0,
              sample(
                c(FALSE, TRUE), block_variables - 1,
                replace = TRUE, prob = c(1 - cross_probability, cross_probability)
              )
            )

          }

        }

      }

      # Obtain full loadings and correlations
      loadings_output <- egm_correlation_optimize(
        loading_structure = loading_structure, correlations = correlations,
        between_indices = as.logical(between_indices)
      )

      # Update parameters
      loading_structure <- loadings_output$loadings # updated with cross-loadings

      # Add cross-loadings that do not contribute to correlations
      for(i in community_sequence){

        # Get block index
        block_index <- membership == i

        # Get number of variables
        block_variables <- sum(block_index)

        # Loop over cross-loadings
        for(j in community_sequence){

          # Except for itself
          if(i != j){

            # Target cross-loadings
            target_cross <- loading_structure[block_index, j]

            # Add cross-loadings
            loading_structure[block_index, j] <-  target_cross +
              rnorm_ziggurat(block_variables) * cross.loadings *
              (target_cross != 0) # ensure sparsity (if sparse)

          }

        }

      }

      # Obtain population correlation matrix
      network <- set_network(loading_structure, membership, c(sample.size, total_variables))

      # Convert network to zero-order correlations
      R <- pcor2cor(network)

      # Check for positive definite
      PD_check <- is_positive_definite(R)

      # Check for cross-loadings that are larger than their assigned loadings
      cross_check <- any(max.col(abs(loading_structure)) != membership)

    }

  }

  # Obtain precision, partial correlations, and adjacency matrices
  K <- solve(R)

  # Perform Cholesky decomposition
  cholesky <- chol(R)

  # Generate data
  data <- MASS_mvrnorm_quick(
    p = total_variables, np = total_variables * sample.size,
    coV = diag(total_variables)
  ) %*% cholesky

  # Set variable names
  colnames(data) <- node_names

  # Return results
  return(
    list(
      data = data,
      population_correlation = R,
      population_precision = K,
      population_network = network,
      parameters = list(
        adjacency = network != 0,
        loadings = loading_structure,
        correlations = loadings_output$correlations, # precise correlations
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

#' @noRd
# beta-min criterion ----
# Updated 06.07.2025
beta_min <- function(P, Q, K, total_variables, sample_size)
{

  # Calculate beta-min
  minimum <- sqrt(Q * log(total_variables) / sample_size)

  # Obtain inverse variances
  inverse_variances <- diag(K)

  # Obtain betas
  inv_K <- outer(inverse_variances, inverse_variances, FUN = "/")
  beta <- P * ((inv_K + t(inv_K)) / 2)

  # Set adjacency matrix
  adjacency <- abs(beta) > minimum

  # Attach minimum
  attr(adjacency, "beta.min") <- minimum

  # Return adjacency matrix
  return(adjacency)

}
