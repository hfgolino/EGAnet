#' Simulate data following a Exploratory Graph Model (\code{EGM})
#'
#' @description Function to simulate data based on \code{EGM}
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
#' Uses \code{runif(n, min = value - 0.025, max = value + 0.025)} for some jitter in the loadings
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
#' acceptable equals 0.10 and robust equals 0.05
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
# Updated 25.02.2026
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
    acceptable = c(0.02, 0.02, 0.90, 0.10),
    robust = c(0.01, 0.01, 0.95, 0.05)
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
  PD_check <- FALSE; cross_check <- TRUE; quality_check <- FALSE; stable <- FALSE

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
    while(!PD_check || cross_check || !quality_check || !stable){

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

      # Initialization flag
      initalize_flag <- count == 1

      # Initialization
      if(initalize_flag){

        # Initialize structure matrix
        loading_structure <- matrix(
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

        }

      }

      # Check for between-community loadings
      if(dimensional_flag){

        # Check if they need to be initialized or re-generated
        if(initalize_flag || stable){

          # Initialize between indices matrix
          between_indices <- matrix(
            0, nrow = total_variables, ncol = communities,
            dimnames = list(
              node_names,
              format_integer(community_sequence, digits(communities) - 1) # community names
            )
          )

          # Re-generate cross-loadings
          for(i in community_sequence){

            # Get block index
            block_index <- membership == i

            # Get number of variables
            block_variables <- sum(block_index)

            # Loop over cross-loadings
            for(j in community_sequence){

              # Except for itself
              if(i != j){

                # Set zero cross-loading indices based on cross-loading probability
                between_indices[block_index, j] <- shuffle( # ensure at least one cross-loading with correlations
                  c(correlations[i,j] == 0, runif_xoshiro(block_variables - 1))
                ) < abs(correlations[i,j])^loading_structure[block_index, i]

              }

            }

          }

          # Obtain full loadings and correlations
          loadings_output <- egm_correlation_optimize(
            loading_structure = loading_structure, correlations = correlations,
            between_indices = as.logical(between_indices),
            membership = membership, total_variables = total_variables
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

      }

      # E-STEP: given loadings, derive network structure
      network <- expected_network(loading_structure, membership, total_variables)
      network_R <- silent_call(try(pcor2cor(network), silent = TRUE))

      # Check for issues
      if(is(network_R, "try-error") || anyNA(network_R) || !is_positive_definite(network_R)){
        next
      }

      # M-STEP: given fixed sparsity and network structure, optimize toward network
      loading_structure[] <- loadings_optimization(
        iter = 10, loadings_vector = as.vector(loading_structure),
        zeros = loading_structure != 0, R = network_R,
        communities = communities,
        data_dimensions = c(sample.size, total_variables),
        opt = "loglik"
      )$par

      # Compute correlations from loadings
      R <- nload2cor(loading_structure)

      # Stable structure
      stable <- all(
        (expected_network(loading_structure, membership, total_variables) != 0) == (network != 0)
      )

      # Check for positive definite
      PD_check <- is_positive_definite(R)

      # Check for cross-loadings that are larger than their assigned loadings
      cross_check <- any(max.col(abs(loading_structure)) != membership)

      # Set quality metrics
      quality_metrics <- c(
        srmr(R, network_R), mean(abs(R - network_R)),
        frobenius(R, network_R), jsd(R, network_R)
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
      ),
      quality = quality_df
    )
  )

}

# Bug checking ----
# communities = 3; variables = 6
# loadings = 0.35; cross.loadings = 0.01
# correlations = 0.30; sample.size = 1000
# quality = "acceptable"; max.iterations = 100
# source("/home/alextops/R/R-packages/EGAnet/R/utils-EGAnet.R")
# source("/home/alextops/R/R-packages/EGAnet/R/helpers.R")
# source("/home/alextops/R/R-packages/EGAnet/R/EGM.optimizations.R")

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
# Updated 25.02.2026
expected_network <- function(loading_structure, membership, total_variables)
{

  # Obtain partial correlations
  P <- nload2pcor(loading_structure)

  # Set Chung-Lu configuration based on interdependence
  assigned_loading <- sqrt(rowSums(loading_structure^2))

  # # Set Chung-Lu configuration based on maximum loading
  # assigned_loading <- nvapply(seq_len(total_variables), function(i){
  #   abs(loading_structure[i, membership[i]])
  # })

  # Return partial correlations
  return(P * (abs(P) > (tcrossprod(assigned_loading) / sum(assigned_loading))))

}

#' @noRd
# Community correlations from loadings ----
# Updated 21.07.2025
community_correlations <- function(simple_structure, loading_structure)
{
  return(cov2cor(crossprod(simple_structure, nload2cor(loading_structure)) %*% simple_structure))
}

#' @noRd
# Network loadings to partial correlations ----
# Updated 03.08.2025
nload2pcor <- function(loadings)
{

  # Compute covariance matrix
  P <- tcrossprod(loadings)

  # Set diagonal to interdependence
  diag(P) <- sqrt(diag(P))

  # Obtain inverse covariance matrix
  INV <- solve(P)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Compute partial correlations
  P <- -D %*% INV %*% D

  # Set diagonal to zero
  diag(P) <- 0

  # Return partial correlation
  return(P)

}

#' @noRd
# Network loadings to correlations ----
# Updated 03.08.2025
nload2cor <- function(loadings)
{

  # Compute covariance matrix
  P <- tcrossprod(loadings)

  # Compute interdependence
  diag(P) <- interdependence <- sqrt(diag(P))

  # Compute matrix I
  I <- diag(sqrt(1 / interdependence))

  # Return correlation
  return(I %*% P %*% I)

}

#' @noRd
# Expected edge values ----
# Updated 26.06.2025
expected_edges <- function(network, data_dimensions = NULL)
{

  # Compute node strength
  strength <- colSums(network)

  # Obtain the normalized cross-product
  EE <- tcrossprod(strength) / sum(strength)

  # Obtain differences as SE
  if(!is.null(data_dimensions)){

    # Delta method for the standard error

    # Ensure absolute
    network <- abs(network)

    # Convert network to Fisher's z
    z_network <- r2z(network)

    # Compute node strength
    strength <- colSums(z_network)
    total_strength <- sum(strength)
    standard_strength <- strength / total_strength
    total_squared <- total_strength^2
    cross_strength <- tcrossprod(strength)
    dEE <- cross_strength / total_squared

    # Get variables minus one
    p_minus_one <- data_dimensions[2] - 1

    # Because of Chung-Lu configuration model, all nodes are
    # assumed to have similar variability
    variance <- p_minus_one * (1 / (data_dimensions[1] - data_dimensions[2] - 1))

    # Compute gradient
    gradient_j <- (matrix(standard_strength, nrow = data_dimensions[2], ncol = data_dimensions[2]) - dEE)^2

    # Get SE (still in Fisher's z)
    # t(gradient_j) = shorthand for gradient i
    # dEE^2 = shorthand for gradient k
    SE <- sqrt(variance * (t(gradient_j) + gradient_j + (p_minus_one - 1) * dEE^2))

    # Absolute partial correlations (assumes similar transformation deviations)
    lower_triangle <- lower.tri(EE)
    average <- sum(EE[lower_triangle] * network[lower_triangle]) / sum(EE[lower_triangle])

    # Calculate the Jacobian
    attr(EE, "SE") <- SE * sqrt((1 - average^2)^2)

  }

  # Return result
  return(EE)

}

#' @noRd
# Creates community structure ----
# Updated 08.10.2024
create_community_structure <- function(
    P, total_variables, communities, community_variables, p.in, p.out
)
{

  # Set vectors of 'p.in' and 'p.out'
  p.in <- swiftelse(length(p.in) == 1, rep(p.in, communities), p.in)
  p.out <- swiftelse(length(p.out) == 1, rep(p.out, communities), p.out)

  # Set diagonal to missing
  diag(P) <- NA

  # Set community blocks
  for(i in seq_len(communities)){

    # Randomly set zero in block
    indices <- P[community_variables[[i]], community_variables[[i]]]

    # Check for current sparsity
    if(compute_density(indices) > p.in[i]){

      # Get lower triangle
      lower_triangle <- lower.tri(indices)

      # Sample to set to zero
      indices[
        abs(indices) < quantile(
          abs(indices[lower_triangle]), probs = 1 - p.in[i], na.rm = TRUE
        )
      ] <- 0

      # Set back into block
      P[community_variables[[i]], community_variables[[i]]] <- indices

    }

    # Set between-community indices
    indices <- P[community_variables[[i]], -unlist(community_variables[-i])]

    # Check for current sparsity
    if(compute_density(indices) > p.out[i]){

      # Get threshold value
      threshold_value <- quantile(abs(indices), probs = 1 - p.out[i], na.rm = TRUE)

      # Set below to zero
      indices[abs(indices) < threshold_value] <- 0

      # Set back into block
      P[community_variables[[i]], -unlist(community_variables[-i])] <- indices

      # Do other side
      indices <- P[-unlist(community_variables[-i]), community_variables[[i]]]

      # Set below to zero
      indices[abs(indices) < threshold_value] <- 0

      # Set back into block
      P[-unlist(community_variables[-i]), community_variables[[i]]] <- indices

    }

  }

  # Set diagonal to zero
  diag(P) <- 0

  # Return community network
  return(P)

}

#' @noRd
# Computes density ----
# Updated 08.10.2024
compute_density <- function(network)
{

  # Get dimensions
  dimensions <- dim(network)

  # Obtain total number of edges
  edges <- prod(dimensions)

  # Check for on-diagonal
  if(dimensions[1] == dimensions[2]){

    # Total possible edges
    total_possible <- (edges - dimensions[2]) / 2

    # Compute sparsity
    return((total_possible - (sum(network == 0, na.rm = TRUE) / 2)) / total_possible)

  }else{ # Otherwise, treat as off-diagonal

    # Compute sparsity
    return((edges - sum(network == 0, na.rm = TRUE)) / edges)

  }

}