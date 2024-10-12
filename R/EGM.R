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
#' \item \code{"EGA"} (default) --- Applies \code{\link[EGAnet]{EGA}} to obtain the
#' (sparse) regularized network structure, communities, and memberships
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
#' @param opt Character vector (length = 1).
#' Fit index to use for optimization of network loadings to the
#' zero-order correlation matrix.
#' Available options include: \code{"AIC"}, \code{"BIC"}, and \code{"SRMR"}.
#' Defaults to \code{"SRMR"}
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
#' # Get depression data
#' data <- na.omit(depression[,24:44])
#'
#' # Estimate EGM (using EGA)
#' egm_ega <- EGM(data)
#'
#' # Estimate EGM (using standard)
#' egm_standard <- EGM(
#'   data, EGM.type = "standard",
#'   communities = 3, # specify number of communities
#'   p.in = 0.95, # probability of edges *in* each community
#'   p.out = 0.80 # probability of edges *between* each community
#' )
#'
#' \dontrun{
#' # Estimate EGM (using search)
#' egm_search <- EGM(
#'   data, EGM.type = "search", communities = 3,
#'   p.in = 0.95 # only need 'p.in'
#' )}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Estimate EGM ----
# Updated 10.10.2024
EGM <- function(
    data, EGM.type = c("search", "standard", "EGA"),
    communities = NULL, structure = NULL,
    p.in = NULL, p.out = NULL,
    opt = c("AIC", "BIC", "SRMR"),
    verbose = TRUE, ...
)
{

  # Set default
  EGM.type <- set_default(EGM.type, "ega", EGM)
  opt <- set_default(opt, "srmr", EGM)

  # Check data and structure
  data <- EGM_errors(
    data, EGM.type, communities, structure,
    p.in, p.out, verbose, ...
  )

  # Switch and return results based on type
  return(
    switch(
      EGM.type,
      "search" = EGM.search(data, communities, structure, p.in, opt, verbose, ...),
      "standard" = EGM.standard(data, communities, structure, p.in, p.out, opt, ...),
      "ega" = EGM.EGA(data, structure, opt, ...)
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
# Updated 09.10.2024
nload2pcor <- function(loadings)
{

  # Compute partial correlation
  P <- tcrossprod(loadings)
  diag(P) <- sqrt(rowSums(loadings^2))

  # Return partial correlation
  return(cor2pcor(cov2cor(P)))

}

#' @noRd
# Network loadings to correlations ----
# Updated 09.10.2024
nload2cor <- function(loadings)
{

  # Compute partial correlation
  P <- tcrossprod(loadings)
  diag(P) <- sqrt(rowSums(loadings^2))

  # Return correlation
  return(cov2cor(P))

}

#' @noRd
# Estimated loadings cost (based on SRMR) ----
# Updated 12.10.2024
srmr_N_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, ...
)
{

  # Assemble loading matrix
  loading_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

  # Obtain assign loadings
  assign_loadings <- loading_matrix[loading_structure]

  # Transpose loadings matrix
  loading_matrix <- t(loading_matrix)

  # Obtain differences
  differences <- abs(loading_matrix) - assign_loadings

  # Obtain difference values
  difference_values <- (differences * (differences > 0))^2

  # Return SRMR
  return(
    srmr(R, nload2cor(loading_matrix)) + # SRMR term
    sqrt(mean(difference_values)) # penalty term
  )

}

# @noRd
# Estimated loadings cost (based on SRMR)
# Updated 12.10.2024
# srmr_N_gradient <- function(
#     loadings_vector, zeros, R,
#     loading_structure, rows, ...
# )
# {
#
#   # NOT USED!!
#   # CLOSEST SO FAR TO ACTUAL GRADIENT BUT NOT CORRECT!!
#
#   # Assemble loading matrix
#   loading_matrix <- t(matrix(loadings_vector, nrow = rows, byrow = TRUE))
#
#   # Compute loading differences
#   differences <- as.vector(
#     net.loads((nload2cor(loading_matrix) - R)^2, ega$wc)$std[colnames(ega$network),]
#   )
#
#   # Return gradient
#   return(-differences / length(differences) * sqrt(mean(differences^2)) * zeros)
#
# }

#' @noRd
# Estimated loadings cost (based on log-likelihood) ----
# Updated 12.10.2024
likelihood_N_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, rows, n, opt, ...
)
{

  # Assemble loading matrix
  loading_matrix <- matrix(loadings_vector * zeros, nrow = rows, byrow = TRUE)

  # Obtain assign loadings
  assign_loadings <- loading_matrix[loading_structure]

  # Transpose loadings matrix
  loading_matrix <- t(loading_matrix)

  # Obtain differences
  differences <- abs(loading_matrix) - assign_loadings

  # Obtain difference values
  difference_values <- differences * (differences > 0)

  # Check for error
  return(
    likelihood(n, rows, nload2cor(loading_matrix), R, loading_matrix)[[opt]] + # likelihood term
    sqrt(mean((difference_values)^2)) # penalty term
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
likelihood <- function(n, p, R, S, loadings, type)
{

  # Get number of communities
  m <- dim(loadings)[2]

  # Log-likelihood
  loglik <- log_likelihood(n, p, R, S, type)

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
# Compute model parameters for TEFI adjustment ----
# Updated 11.10.2024
compute_tefi_adjustment <- function(loadings, correlations)
{

  # Obtain dimensions
  dimensions <- dim(loadings)

  # Total number of parameters
  parameters <- (dimensions[1] * dimensions[2]) + dimensions[1] + ((dimensions[2] * (dimensions[2] - 1)) / 2)

  # Obtain model parameters
  model_parameters <- parameters - sum(loadings == 0)

  # Return TEFI adjustment
  # return(model_parameters * log(abs(mean(correlations[lower.tri(correlations)]))) - log(mean(loadings)))
  # Maybe...
  return(model_parameters * log(mean(abs(correlations[lower.tri(correlations)]))) - log(mean(abs(loadings))))

}

#' @noRd
# EGM | Standard ----
# Updated 12.10.2024
EGM.standard <- function(data, communities, structure, p.in, p.out, opt, ...)
{

  # Get ellipse
  ellipse <- list(...)

  # Get dimensions
  data_dimensions <- dim(data)
  dimension_names <- dimnames(data)

  # Determine if sample size was provided
  if(data_dimensions[1] == data_dimensions[2]){

    # Check if sample size was provided
    if("n" %in% names(ellipse)){

      # Set sample size
      data_dimensions[1] <- ellipse$n

      # Set empirical zero-order correlations
      empirical_R <- data

    }else{

      # Stop and send error
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "A symmetric (", data_dimensions[1], " x ",
          data_dimensions[2], ") matrix was input but sample size was not provided.\n\n",
          "If you'd like to use a correlation matrix, please set the argument 'n' to your sample size"
        ),
        call = "EGM"
      )

    }

  }else{

    # Estimate zero-order correlations
    empirical_R <- auto.correlate(data, ...)

  }

  # Obtain partial correlations
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
    P = empirical_P, total_variables = data_dimensions[2],
    communities = communities,
    community_variables = community_variables,
    p.in = p.in, p.out = p.out
  )

  # Update loadings
  output <- silent_call(net.loads(community_P, structure, ...))
  output$std <- output$std[dimension_names[[2]],]

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
  loading_structure <- matrix(FALSE, nrow = communities, ncol = data_dimensions[2])

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[structure[i], i] <- TRUE
  }

  # Switch out optimization criterion
  cost_FUN <- switch(
    opt,
    "aic" = likelihood_N_cost,
    "bic" = likelihood_N_cost,
    "srmr" = srmr_N_cost
  )

  # Optimize network loadings
  result <- silent_call(
    nlm(
      p = loadings_vector, f = cost_FUN,
      zeros = zeros, R = empirical_R,
      loading_structure = loading_structure,
      rows = communities, n = data_dimensions[1],
      opt = toupper(opt), iterlim = 1000
    )
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$estimate, nrow = data_dimensions[2],
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

  # Compute adjusted TEFI (for optimized)
  optimized_tefi_adjusted <- compute_tefi_adjustment(optimized_loadings, optimized_correlations)

  # Set up EGA
  ega_list <- list(
    dim.variables = data.frame(
      items = dimension_names[[2]],
      dimension = structure
    ),
    network = community_P, wc = structure,
    n.dim = unique_length(structure), correlation = empirical_R,
    n = data_dimensions[1], TEFI_adj = tefi(empirical_R, structure)$VN.Entropy.Fit + optimized_tefi_adjusted
  ); class(ega_list) <- "EGA"

  # Attach methods to network
  attr(ega_list$network, which = "methods") <- list(
    model = "egm", communities = communities, p.in = p.in, p.out = p.out
  )

  # Attach class to memberships
  class(ega_list$wc) <- "EGA.community"

  # Attach methods to memberships
  names(ega_list$wc) <- dimension_names[[2]]
  attr(ega_list$wc, which = "methods") <- list(
    algorithm = swiftelse(is.null(structure), NULL, "walktrap"),
    objective_function = NULL
  )

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
          P.srmr = srmr(empirical_P, standard_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = standard_R, S = empirical_R, loadings = output$std
          ),
          TEFI_adj = tefi(standard_R, structure = structure)$VN.Entropy.Fit + compute_tefi_adjustment(
            output$std, standard_correlations
          )
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = c(
          R.srmr = srmr(empirical_R, optimized_R),
          P.srmr = srmr(empirical_P, optimized_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = optimized_R, S = empirical_R, loadings = optimized_loadings
          ),
          TEFI_adj = tefi(optimized_R, structure = structure)$VN.Entropy.Fit + optimized_tefi_adjusted
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- "EGM"

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
# Updated 11.10.2024
EGM.EGA <- function(data, structure, opt, ...)
{

  # Obtain data dimensions
  data_dimensions <- dim(data)

  # Estimate EGA
  ega <- EGA(data, plot.EGA = FALSE, ...)
  empirical_P <- cor2pcor(ega$correlation)

  # Update number of rows based on EGA
  data_dimensions[1] <- ega$n

  # Obtain variable names from the network
  variable_names <- dimnames(ega$network)[[2]]

  # Set memberships based on structure
  if(!is.null(structure)){

    # Update EGA features
    ega$wc[] <- structure
    ega$dim.variables$dimension <- structure

  }else{
    structure <- ega$wc
  }

  # Set communities
  ega$n.dim <- communities <- unique_length(structure)

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
  loadings_length <- length(loadings_vector)

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

  # Switch out optimization criterion
  cost_FUN <- switch(
    opt,
    "aic" = likelihood_N_cost,
    "bic" = likelihood_N_cost,
    "srmr" = srmr_N_cost
  )

  # Optimize network loadings
  result <- silent_call(
    nlm(
      p = loadings_vector, f = cost_FUN,
      zeros = zeros, R = ega$correlation,
      loading_structure = loading_structure,
      rows = communities, n = data_dimensions[1],
      opt = toupper(opt), iterlim = 1000
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
          P.srmr = srmr(empirical_P, standard_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = standard_R, S = ega$correlation, loadings = standard_loadings
          ),
          TEFI_adj = tefi(standard_R, structure = ega$wc)$VN.Entropy.Fit + compute_tefi_adjustment(
            standard_loadings, standard_correlations
          )
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = c(
          R.srmr = srmr(ega$correlation, optimized_R),
          P.srmr = srmr(empirical_P, optimized_P),
          likelihood(
            n = data_dimensions[1], p = data_dimensions[2],
            R = optimized_R, S = ega$correlation, loadings = optimized_loadings
          ),
          TEFI_adj = tefi(optimized_R, structure = ega$wc)$VN.Entropy.Fit + compute_tefi_adjustment(
            optimized_loadings, optimized_correlations
          )
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- "EGM"

  # Return results
  return(results)

}

#' @noRd
# EGM | Search ----
# Updated 10.10.2024
EGM.search <- function(data, communities, structure, p.in, opt, verbose, ...)
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
            p.out = p_grid$p_out[i], opt = opt
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
      swiftelse(is.null(x), NA, x$model$optimized$fit[[toupper(opt)]])
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
  class(results) <- "EGM"

  # Return results
  return(results)

}