#' Exploratory Graph Model
#'
#' @description Function to fit the Exploratory Graph Model
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{community.unidimensional}}, and
#' \code{\link[EGAnet]{EGA}}
#'
#' @examples
#' # Estimate EGM
#' wmt_egm <- EGM(wmt2[,7:24])
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#
# Estimate EGM ----
# Updated 04.10.2024
EGM <- function(data, structure = NULL, ...)
{

  # Check data and structure
  data <- EGM_errors(data, structure, ...)

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
  loading_structure <- matrix(
    FALSE, nrow = dimensions[1],
    ncol = dimensions[2],
    dimnames = dimension_names
  )

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[i, structure[i]] <- TRUE
  }

  # Use optimize to minimize the SRMR
  result <- silent_call(
    nlm(
      p = loadings_vector, f = estimated_N_cost,
      zeros = zeros, R = ega$correlation,
      loading_structure = loading_structure,
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
          likelihood(data, standard_R, ega$correlation, standard_loadings),
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
          likelihood(data, optimized_R, ega$correlation, optimized_loadings),
          TEFI = tefi(optimized_R, structure = ega$wc)$VN.Entropy.Fit
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
# EGM Errors ----
# Updated 04.10.2023
EGM_errors <- function(data, structure, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "EGA")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Check for NULL structure
  if(!is.null(structure)){

    # If not NULL, check 'structure' errors
    typeof_error(structure, c("character", "numeric"), "EGM")
    object_error(structure, "vector", "EGM")
    length_error(structure, dim(data)[2], "EGM")

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

#' @noRd
# Estimated loadings cost (based on SRMR) ----
# Updated 04.10.2024
estimated_N_cost <- function(
    loadings_vector, zeros, R,
    loading_structure, ...
)
{

  # Assemble loading matrix
  loading_matrix <- matrix(loadings_vector * zeros, nrow = nrow(R))

  # Obtain assign loadings
  assign_loadings <- apply(
    loading_matrix * loading_structure, 1, function(x){
    x[x != 0]
  })

  # Obtain differences
  differences <- abs(loading_matrix) - abs(assign_loadings)

  # Obtain difference values
  difference_values <- differences * sweep(
    x = differences, MARGIN = 2, STATS = 0, FUN = ">"
  )

  # Set penalties
  penalty <- sqrt(mean((difference_values)^2))

  # Try result
  implied_R <- try(nload2cor(loading_matrix), silent = TRUE)

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
# Compute log-likelihood metrics ----
# Updated 04.10.2024
likelihood <- function(data, R, S, loadings)
{

  # Obtain dimensions
  dimensions <- dim(data)
  n <- dimensions[1] # sample size
  p <- dimensions[2] # number of variables
  m <- dim(loadings)[2] # number of dimensions

  # Log-likelihood
  loglik <- -(n / 2) * (
    p * log(2 * pi) + log(det(R)) + sum(diag(
      S %*% solve(R)
    ))
  )

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
