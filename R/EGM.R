#' Exploratory Graph Model
#'
#' @description Function to fit the Exploratory Graph Model
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param EGM.type Character vector (length = 1).
#' Sets the procedure to conduct \code{EGM}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"search"} --- Searches over \code{p.in} and \code{p.out}
#' parameters for best fit based on Bayesian information criterion (BIC).
#' Uses \code{EGM.type = "standard"} under the hood. Only the argument
#' \code{p.in} is used such that \code{p.in} searches over
#' \code{seq(p.in, 1, 0.05)} and \code{p.out} searches over
#' \code{seq(0.00, p.in, 0.05)}
#'
#' \item \code{"standard"} --- Applies the standard EGM model which
#' estimates communities based on the non-regularized empirical partial
#' correlation matrix and sparsity is set using \code{p.in} and \code{p.out}
#'
#' \item \code{"EGA"} --- Applies \code{\link[EGAnet]{EGA}} to obtain the
#' (sparse) regularized network structure, communities, and memberships
#'
#' }
#'
#' @param communities Numeric vector (length = 1).
#' Number of communities to use for the \code{"standard"} type of EGM.
#' Defaults to \code{NULL}.
#' Providing no input will use the communities and memberships output
#' from the Walktrap algorithm (\code{\link[igraph]{cluster_walktrap}}) based
#' on the empirical non-regularized partial correlation matrix
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}.
#' Defaults to \code{NULL}
#'
#' @param p.in Numeric vector (length = 1).
#' Probability that a node is randomly linked to other nodes in the same community.
#' Within community edges are set to zero based on \code{quantile(x, prob = 1 - p.in)}
#' ensuring the lowest edge values are set to zero (i.e., most probable to \emph{not}
#' be randomly connected).
#' Only used for \code{EGM.type = "standard"}.
#' Defaults to \code{NULL} but must be set
#'
#' @param p.out Numeric vector (length = 1).
#' Probability that a node is randomly linked to other nodes \emph{not} in the same community.
#' Between community edges are set to zero based on \code{quantile(x, prob = 1 - p.out)}
#' ensuring the lowest edge values are set to zero (i.e., most probable to \emph{not}
#' be randomly connected).
#' Only used for \code{EGM.type = "standard"}.
#' Defaults to \code{NULL} but must be set
#'
#' @param ncores Numeric (length = 1).
#' Number of cores to use in computing results.
#' Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then type: \code{parallel::detectCores()}
#'
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{community.unidimensional}},
#' \code{\link[EGAnet]{EGA}}, and
#' \code{\link[EGAnet]{net.loads}}
#'
#' @examples
#' # Estimate EGM
#' wmt_egm <- EGM(wmt2[,7:24])
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Estimate EGM ----
# Updated 07.10.2024
EGM <- function(
    data, EGM.type = c("search", "standard", "EGA"),
    communities = NULL, structure = NULL,
    p.in = NULL, p.out = NULL, verbose = TRUE, ...
)
{

  # Set default
  EGM.type <- set_default(EGM.type, "standard", EGM)

  # Check data and structure
  data <- EGM_errors(
    data, EGM.type, communities, structure,
    p.in, p.out, verbose, ...
  )

  # Switch and return results based on type
  return(
    switch(
      EGM.type,
      "search" = EGM.search(data, communities, structure, p.in, verbose, ...),
      "standard" = EGM.standard(data, communities, structure, p.in, p.out, ...),
      "ega" = EGM.EGA(data, structure, ...)
    )
  )

}

#' @noRd
# EGM Errors ----
# Updated 07.10.2023
EGM_errors <- function(
    data, EGM.type, communities, structure,
    p.in, p.out, verbose, ...
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "EGM")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Check for NULL communities
  if(!is.null(communities)){

    # If not NULL, check 'communities' errors
    typeof_error(communities, "numeric", "EGM")
    object_error(communities, "vector", "EGM")
    length_error(communities, 1, "EGM")

  }

  # Check for NULL structure
  if(!is.null(structure)){

    # If not NULL, check 'structure' errors
    typeof_error(structure, c("character", "numeric"), "EGM")
    object_error(structure, "vector", "EGM")
    length_error(structure, dim(data)[2], "EGM")

  }

  # Check first for parameters involved in both search and standard
  if(EGM.type != "ega"){

    # Check for NULL in 'p.in'
    if(is.null(p.in)){
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Input for 'p.in' is `NULL`. A value between 0 and 1 for 'p.in' must be provided."
        ),
        call = "EGM"
      )
    }

    # Check 'p.in' errors
    typeof_error(p.in, "numeric", "EGM")
    range_error(p.in, c(0, 1), "EGM")
    length_error(p.in, c(1, communities), "EGM")

  }

  # Check for EGM type
  if(EGM.type == "search"){

    # 'verbose' errors
    length_error(verbose, 1, "EGM")
    typeof_error(verbose, "logical", "EGM")

  }else if(EGM.type == "standard"){

    # Check for NULL in 'p.out'
    if(is.null(p.out)){
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Input for 'p.out' is `NULL`. A value between 0 and 1 for 'p.out' must be provided."
        ),
        call = "EGM"
      )
    }

    # Check 'p.out' errors
    typeof_error(p.out, "numeric", "EGM")
    range_error(p.out, c(0, 1), "EGM")
    length_error(p.out, c(1, communities), "EGM")

  }

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, FALSE)
  }

  # Return data in case of tibble
  return(data)

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
  diag(P) <- sqrt(1 - uniqueness)

  # Return partial correlation
  return(cor2pcor(cov2cor(P)))

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

}

# Estimated loadings cost (based on SRMR) ----
# Updated 06.10.2024
estimated_N_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, ...
)
{

  # Assemble loading matrix
  loading_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

  # Obtain assign loadings
  assign_loadings <- loading_matrix[loading_structure]

  # Obtain differences
  differences <- abs(t(loading_matrix)) - assign_loadings

  # Obtain difference values
  difference_values <- differences * (differences > 0)

  # Set penalties
  penalty <- sqrt(mean((difference_values)^2))

  # Try result
  implied_R <- try(nload2cor(t(loading_matrix)), silent = TRUE)

  # Check for error
  return(
    swiftelse( # Return infinite on error
      is(implied_R, "try-error"), Inf,
      srmr(R, implied_R) + penalty
      # Return SRMR otherwise
    )
  )

}

#' @noRd
# Log-likelihood only ----
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
# Compute log-likelihood metrics ----
# Updated 06.10.2024
likelihood <- function(n, p, R, S, loadings)
{

  # Get number of communities
  m <- dim(loadings)[2]

  # Log-likelihood
  loglik <- log_likelihood(n, p, R, S)

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
      # ),
      # # -2L + klog(n) + 2 gamma log(binom(pk))
      # GFI = 1 - sum((R - S)^2) / sum(S^2)
    )
  )

}

#' @noRd
# EGM | Standard ----
# Updated 07.10.2024
EGM.standard <- function(data, communities, structure, p.in, p.out, ...)
{

  # Get dimensions
  dimensions <- dim(data)

  # Estimate zero-order and partial correlations
  empirical_R <- auto.correlate(data, ...)
  empirical_P <- cor2pcor(empirical_R)

  # Check for whether structure is provided
  if(is.null(structure)){

    # Obtain structure if not provided
    walktrap <- igraph::cluster_walktrap(convert2igraph(abs(empirical_P)))

    # Obtain memberships based on number of communities
    structure <- swiftelse(
      is.null(communities) || unique_length(walktrap$membership) == communities,
      walktrap$membership, cutree(as.hclust(walktrap), k = communities)
    )

  }

  # Set number of communities
  communities <- unique_length(structure)

  # Set up community variables
  community_variables <- lapply(
    seq_len(communities), function(community){structure == community}
  )

  # Obtain community structure
  community_P <- create_community_structure(
    P = empirical_P, total_variables = dimensions[2],
    communities = communities,
    community_variables = community_variables,
    p.in = p.in, p.out = p.out
  )

  # Update loadings
  output <- silent_call(net.loads(community_P, structure, ...))
  output$std <- output$std[colnames(data),]

  # Compute network scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  standard_correlations <- cor(standard_scores, use = "pairwise")

  # Compute model-implied correlations
  standard_R <- nload2cor(output$std)
  standard_P <- cor2pcor(standard_R)

  # Obtain loadings vector and get bounds
  loadings_vector <- as.vector(output$std)
  zeros <- loadings_vector != 0

  # Set up loading structure
  # Uses transpose for 2x speed up in optimization
  loading_structure <- matrix(FALSE, nrow = communities, ncol = dimensions[2])

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[structure[i], i] <- TRUE
  }

  # Use optimize to minimize the SRMR
  result <- silent_call(
    nlm(
      p = loadings_vector, f = estimated_N_cost,
      zeros = zeros, R = empirical_R,
      loading_structure = loading_structure,
      rows = communities, iterlim = 1000
    )
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$estimate, nrow = dimensions[2],
    dimnames = dimnames(output$std)
  )

  # Replace loadings output with optimized loadings
  output$std <- optimized_loadings

  # Compute network scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Compute model-implied correlations
  optimized_R <- nload2cor(optimized_loadings)
  optimized_P <- cor2pcor(optimized_R)

  # Set up EGA
  ega_list <- list(
    dim.variables = data.frame(
      items = dimnames(data)[[2]],
      dimension = structure
    ),
    network = community_P, wc = structure,
    n.dim = unique_length(structure), correlation = empirical_R,
    n = dimensions[2], TEFI = tefi(empirical_R, structure)
  ); class(ega_list) <- "EGA"

  # Set up results
  results <- list(
    EGA = ega_list, structure = structure,
    model = list( # general list
      standard = list( # using standard parameters
        loadings = output$std,
        scores = standard_scores,
        correlations = standard_correlations,
        fit = c(
          R.srmr = srmr(empirical_R, standard_R),
          P.srmr = srmr(cor2pcor(empirical_R), standard_P),
          likelihood(
            n = dimensions[1], p = dimensions[2],
            R = standard_R, S = empirical_R, loadings = output$std
          ),
          TEFI = tefi(standard_R, structure = structure)$VN.Entropy.Fit
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = c(
          R.srmr = srmr(empirical_R, optimized_R),
          P.srmr = srmr(cor2pcor(empirical_R), optimized_P),
          likelihood(
            n = dimensions[1], p = dimensions[2],
            R = optimized_R, S = empirical_R, loadings = optimized_loadings
          ),
          TEFI = tefi(optimized_R, structure = structure)$VN.Entropy.Fit
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- c("EGM", "standard")

  # Return results
  return(results)

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

#' @noRd
# EGM | EGA ----
# Updated 07.10.2024
EGM.EGA <- function(data, structure, ...)
{

  # Obtain data dimensions
  data_dimensions <- dim(data)

  # Estimate EGA
  ega <- EGA(data, plot.EGA = FALSE, ...)

  # Obtain variable names from the network
  variable_names <- dimnames(ega$network)[[2]]

  # Set memberships based on structure
  if(!is.null(structure)){
    ega$wc[] <- structure
  }else{
    structure <- ega$wc
  }

  # Obtain standard network loadings
  output <- silent_call(net.loads(A = ega$network, wc = ega$wc, ...))

  # Obtain standard loadings
  standard_loadings <- output$std[variable_names,]

  # Compute network scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  standard_correlations <- cor(standard_scores, use = "pairwise")

  # Obtain model-implied correlations
  standard_R <- nload2cor(standard_loadings)
  standard_P <- cor2pcor(standard_R)

  # Get loading dimensions
  dimensions <- dim(standard_loadings)
  dimension_names <- dimnames(standard_loadings)

  # Obtain loadings vector and get bounds
  loadings_vector <- as.vector(standard_loadings)
  zeros <- loadings_vector != 0

  # Set up loading structure
  # Uses transpose for 2x speed up in optimization
  loading_structure <- matrix(
    FALSE, nrow = dimensions[2],
    ncol = dimensions[1],
    dimnames = list(dimension_names[[2]], dimension_names[[1]])
  )

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[structure[i], i] <- TRUE
  }

  # Use optimize to minimize the SRMR
  result <- silent_call(
    nlm(
      p = loadings_vector, f = estimated_N_cost,
      zeros = zeros, R = ega$correlation,
      loading_structure = loading_structure,
      rows = dimensions[2],
      iterlim = 1000
    )
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$estimate, nrow = dimensions[1],
    dimnames = dimension_names
  )

  # Replace loadings output with optimized loadings
  output$loadings$std <- optimized_loadings

  # Compute network scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Obtain model-implied correlations
  optimized_R <- nload2cor(optimized_loadings)
  optimized_P <- cor2pcor(optimized_R)

  # Set up results
  results <- list(
    EGA = ega, structure = structure,
    model = list( # general list
      standard = list( # using standard parameters
        loadings = standard_loadings,
        scores = standard_scores,
        correlations = standard_correlations,
        fit = c(
          R.srmr = srmr(ega$correlation, standard_R),
          P.srmr = srmr(cor2pcor(ega$correlation), standard_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = standard_R, S = ega$correlation, loadings = standard_loadings
          ),
          TEFI = tefi(standard_R, structure = ega$wc)$VN.Entropy.Fit
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = c(
          R.srmr = srmr(ega$correlation, optimized_R),
          P.srmr = srmr(cor2pcor(ega$correlation), optimized_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = optimized_R, S = ega$correlation, loadings = optimized_loadings
          ),
          TEFI = tefi(optimized_R, structure = ega$wc)$VN.Entropy.Fit
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- c("EGM", "EGA")

  # Return results
  return(results)

}

#' @noRd
# EGM | Search ----
# Updated 07.10.2024
EGM.search <- function(data, communities, structure, p.in, verbose, ...)
{

  # Perform search based on 'p.in'
  p_grid <- expand.grid(
    p_in = seq(p.in, 1, 0.05), p_out = seq(0.00, p.in, 0.05)
  )

  # Get number of search
  p_number <- dim(p_grid)[1]

  # Loop over grid search
  grid_search <- parallel_process(
    iterations = p_number, datalist = seq_len(p_number),
    FUN = function(i, data, p_grid, communities){

      # First, try
      output <- silent_call(
        try(
          EGM(
            data = data, EGM.type = "standard",
            communities = communities, p.in = p_grid$p_in[i],
            p.out = p_grid$p_out[i]
          ), silent = TRUE
        )
      )

      # Return result
      return(swiftelse(is(output, "try-error"), NULL, output))

    }, data, p_grid, communities, ncores = 1, progress = verbose
  )

  # Obtain fits
  optimized_fits <- nvapply(
    grid_search, function(x){
      swiftelse(is.null(x), NA, x$model$optimized$fit[["BIC"]])
    }
  )

  # Obtain index
  min_index <- which.min(optimized_fits)

  # Set up final model
  results <- grid_search[[min_index]]

  # Add 'p.in' and 'p.out' parameters
  results$search <- c(
    p.in = p_grid[min_index, "p_in"],
    p.out = p_grid[min_index, "p_out"]
  )

  # Overwrite class
  class(results) <- c("EGM", "search")

  # Return results
  return(results)

}