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
#' Uses \code{runif(n, min = value - 0.025, max = value + 0.025)}
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
#' @export
#'
# Simulate EGM ----
# Updated 26.03.2025
simEGM <- function(
    communities, variables,
    loadings, cross.loadings = 0.01,
    correlations, sample.size,
    p.in = 0.95, p.out = 0.80,
    max.iterations = 1000
){

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

  # Set names
  names(membership) <- node_names <- paste0(
    "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
  )
  community_names <- format_integer(community_sequence, digits(communities) - 1)

  # Set up network structure
  network_structure <- matrix(
    0, nrow = total_variables, ncol = total_variables,
    dimnames = list(node_names, node_names)
  )

  # Initialize structure matrix
  structure_matrix <- matrix(
    0, nrow = total_variables, ncol = communities,
    dimnames = list(node_names, community_names)
  )

  # Loop over to create communities (block-diagonal)
  for(i in community_sequence){

    # Structure index
    structure_index <- membership == i

    # Set structure matrix
    structure_matrix[structure_index, i] <- 1

    # Identity community
    within_network <- network_structure[structure_index, structure_index]

    # Obtain lower triangle
    lower_triangle <- lower.tri(within_network)

    # Get number of indices
    n_lower <- sum(lower_triangle)

    # Initialize indices
    index <- numeric(n_lower)

    # Set density
    index[seq_len(round(n_lower * p.in))] <- 1

    # Shuffle indices back into network
    network_structure[structure_index, structure_index][lower_triangle] <- sample(index)

  }

  # Loop over to create non-communities (off-diagonal)
  for(i in community_sequence){
    for(j in community_sequence){

      # Only do greater current community (lower triangle)
      if(j > i){

        # Structure index
        block_index <- membership == i
        off_index <- membership == j

        # Non-community
        between_network <- network_structure[off_index, block_index]

        # Get number of indices
        n_between <- length(between_network)

        # Initialize indices
        index <- numeric(n_between)

        # Set density
        index[seq_len(round(n_between * p.out))] <- 1

        # Shuffle indices back into network
        network_structure[off_index, block_index] <- sample(index)

      }
    }
  }

  # Make adjacency matrix symmetric
  network_structure <- network_structure + t(network_structure)

  # Get loading structure
  loading_structure <- network_structure %*% structure_matrix

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
      block_index <- structure_matrix[,i] == 1

      # Obtain degree
      block_degree <- loading_structure[block_index, i]

      # Order by degree
      block_order <- order(block_degree, decreasing = TRUE)

      # Set assigned loadings
      loading_structure[block_index, i][block_order] <- sort(
        runif_xoshiro(sum(block_index), min = loadings - 0.075, max = loadings + 0.075),
        decreasing = TRUE
      ) * swiftelse(block_degree[block_order] == 0, 0, 1)

      # Get off-diagonal indices
      off_index <- structure_matrix[,i] == 0

      # Obtain degree
      off_degree <- loading_structure[off_index, i]

      # Order by degree
      off_order <- order(off_degree, decreasing = TRUE)

      # Off length
      off_length <- sum(off_index)

      # Set range
      correlation_range <- 0.25 * correlations /
        ((1 - range(loading_structure[block_index, i])) * sqrt(log(total_variables)))

      # Set correlations
      loading_structure[off_index, i] <- runif_xoshiro(
        off_length, min = min(correlation_range), max = max(correlation_range)
      )

      # Add cross-loadings
      loading_structure[off_index, i][off_order] <- (
        loading_structure[off_index, i][off_order] +
          sort(rnorm_ziggurat(off_length) * cross.loadings, decreasing = TRUE)
      ) * swiftelse(off_degree[off_order] == 0, 0, 1)

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

  # Return results
  return(
    list(
      data = data,
      population_correlation = R,
      population_partial_correlation = cor2pcor(R),
      parameters = list(
        network = network_structure,
        loadings = loading_structure,
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