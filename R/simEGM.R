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
#' @param density.power Numeric (length = 1).
#' Controls the sparsity of the network loading matrix.
#' The probability any given cross-loading will be set to zero
#' is based on the community correlation between two communities
#' raised to \eqn{\frac{1}{density.power}}. Accepts values between
#'  \code{1} and \code{Inf} but sparser loadings (i.e., values less than
#'  \code{3}) can lead to convergence issues.
#'  Defaults to \code{3}
#'
#' @param quality Character (length = 1).
#' Quality metrics related to the alignment of the correlations
#' implied by the loadings and network are computed with certain
#' standards in place to accept a solution.
#' These metrics include:
#'
#' \itemize{
#'
#' \item SRMR (or RMSE) --- standardized root mean residual where
#' acceptable equals 0.02 and robust equals 0.01
#'
#' \item MAE --- mean absolute error where
#' acceptable equals 0.02 and robust equals 0.01
#'
#' \item \code{\link[EGAnet]{frobenius}} --- Frobenius norm where
#' acceptable equals 0.90 and robust equals 0.95
#'
#' \item \code{\link[EGAnet]{jsd}} --- Jensen-Shannon Distance where
#' acceptable equals 0.05 and robust equals 0.025
#'
#' }
#'
#' Defaults to \code{"acceptable"}.
#' \code{"robust"} is available but most often needs \code{density.power}
#' to be increased to allow for more cross-loadings to converge
#'
#' @param max.iterations Numeric (length = 1).
#' Number of iterations to attempt to get convergence before erroring out.
#' Defaults to \code{100}
#'
#' @examples
#' simulated <- simEGM(
#'   communities = 2, variables = 6,
#'   loadings = 0.35, # use network loading sizes
#'   correlations = 0.30, sample.size = 1000
#' )
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Simulate EGM ----
# Updated 29.07.2025
simEGM <- function(
    communities, variables,
    loadings, cross.loadings = 0.01,
    correlations, sample.size,
    quality = c("acceptable", "robust"),
    max.iterations = 100
){

  # Set quality argument
  quality <- set_default(quality, "acceptable", simEGM)

  # Set quality comparisons
  quality_comp <- data.frame(
    acceptable = c(0.02, 0.02, 0.90, 0.05),
    robust = c(0.01, 0.01, 0.95, 0.025)
  )

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

  # Multidimensional flag
  dimensional_flag <- communities > 1

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
  PD_check <- FALSE; cross_check <- TRUE; quality_check <- FALSE

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
    while(!PD_check || cross_check || !quality_check){

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
            "  +  decreasing 'correlations'\n",
            "  +  setting quality to 'acceptable'"
          ),
          call = "simEGM"
        )

      }

      # Initialize structure matrix
      between_indices <- loading_structure <- matrix(
        0, nrow = total_variables, ncol = communities,
        dimnames = list(
          node_names,
          format_integer(community_sequence, digits(communities) - 1) # community names
        )
      )

      # Generate within-community loadings
      for(i in community_sequence){

        # Get block index
        block_index <- membership == i

        # Get number of variables
        block_variables <- sum(block_index)

        # Set assigned loadings
        loading_structure[block_index, i] <- runif_xoshiro(
          block_variables, min = loadings[i] - 0.025, max = loadings[i] + 0.025
        )

        # Check for multidimensional
        if(dimensional_flag){

          # Loop over cross-loadings
          for(j in community_sequence){

            # Except for itself
            if(i != j){

              # Set zero cross-loading indices based on cross-loading probability
              between_indices[block_index, j] <- shuffle( # ensure at least one cross-loading with correlations
                c(correlations[i,j] == 0, runif_xoshiro(block_variables - 1))
              ) < abs(correlations[i,j])^(1 / variables[i])

            }

          }

        }

      }

      # Check for multidimensional
      if(dimensional_flag){

        # Set within indices for simple structure
        within_indices <- loading_structure != 0

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

      }

      # Set correlations for loadings
      R <- nload2cor(loading_structure)
      P <- cor2pcor(R)

      # Obtain network matrix based on Chung-Lu expectation of simple structure
      network <- expected_network(loading_structure, membership, total_variables)
      network_R <- silent_call(try(pcor2cor(network), silent = TRUE))

      # Check for issues
      if(is(network_R, "try-error") || anyNA(network_R) || !is_positive_definite(network_R)){
        next
      }

      # Set lower triangle
      lower_triangle <- lower.tri(R)

      # Optimize network toward loadings
      network_vector <- as.vector(network[lower_triangle])
      zeros <- network_vector != 0

      # Update network edge parameters
      network_vector[zeros] <- egm_network_optimize(
        network_vector = network_vector[zeros],
        R = R, n = sample.size, v = total_variables,
        lower_triangle = lower_triangle,
        zeros = zeros, opt = "srmr"
      )$par

      # Update network
      network[lower_triangle] <- network_vector
      network <- t(network)
      network[lower_triangle] <- network_vector

      # Check for positive definite
      PD_check <- is_positive_definite(R)

      # Check for cross-loadings that are larger than their assigned loadings
      cross_check <- any(max.col(abs(loading_structure)) != membership)

      # Set quality metrics
      quality_metrics <- c(
        srmr(P, network), mean(abs(P - network)),
        frobenius(P, network), jsd(P, network)
      )

      # Quality metric check
      quality_df <- data.frame(
        Metric = c("SRMR", "MAE", "Frobenius", "JSD"),
        Value = quality_metrics,
        Acceptable = quality_comp$acceptable,
        Robust = quality_comp$robust,
        Quality = c( # Use robust
          quality_metrics[1] < quality_comp[[quality]][1],
          quality_metrics[2] < quality_comp[[quality]][2],
          quality_metrics[3] > quality_comp[[quality]][3],
          quality_metrics[4] < quality_comp[[quality]][4]
        )
      )

      # Check that quality metrics are all satisfied
      quality_check <- sum(quality_df$Quality) == 4

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
        correlations = correlations,
        membership = membership,
        iterations = count
      )
    )
  )

}

#' @noRd
# Errors ----
# Updated 29.07.2025
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
# Expected network ----
# Updated 30.07.2025
expected_network <- function(loading_structure, membership, total_variables)
{

  # Obtain partial correlations
  P <- nload2pcor(loading_structure)

  # Set edges to zero for zero loadings
  for(i in seq_len(total_variables)){

    # Check for zero loadings
    zero_index <- which(loading_structure[i,] == 0)

    # Identify whether zero loadings exist
    if(length(zero_index) != 0){

      # Loop over zero loading memberships
      for(community in zero_index){

        # Get index
        index <- membership == community

        # Set values to zero
        P[index, i] <- P[i, index] <- 0

      }

    }

  }

  # Set Chung-Lu configuration based on maximum cross-loading
  max_cross <- nvapply(seq_len(total_variables), function(i){
    max(abs(loading_structure[i, -membership[i]]))
  })

  # Update P and get implied
  updated_P <- P * (abs(P) > (tcrossprod(max_cross) / sum(max_cross)))
  implied_R <- silent_call(pcor2cor(updated_P))

  # Check for positive definite
  if(anyNA(implied_R) || !is_positive_definite(implied_R)){
    return(P) # cannot make more sparse
  }else{
    return(updated_P) # sparsity addition succeeded
  }

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
