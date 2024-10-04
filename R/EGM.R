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
  }

  # Obtain standard network loadings
  output <- silent_call(net.loads(A = ega$network, wc = ega$wc, ...))

  # Obtain standard loadings
  standard_loadings <- output$std[variable_names,]

  # Compute network scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  standard_correlations <- cor(standard_scores)

  # Obtain model-implied correlations
  standard_R <- nload2cor(standard_loadings)
  standard_P <- cor2pcor(standard_R)

  # Obtain loadings vector and get bounds
  loading_vector <- as.vector(standard_loadings)
  zeros <- loading_vector != 0

  # Use optimize to minimize the SRMR
  result <- silent_call(
    nlm(
      p = loading_vector, f = estimated_N_cost,
      zeros = zeros, R = ega$correlation,
      iterlim = 1000
    )
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$estimate, nrow = nrow(standard_loadings),
    dimnames = dimnames(standard_loadings)
  )

  # Replace loadings output with optimized loadings
  output$loadings$std <- optimized_loadings

  # Compute network scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores)

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

#' Estimated loadings cost (based on SRMR) ----
#' @noRd
# Updated 04.10.2024
estimated_N_cost <- function(loadings_vector, zeros, R, ...)
{
  return(srmr(R, nload2cor(matrix(loadings_vector * zeros, nrow = nrow(R)))))
}

#' Function to compute log-likelihood metrics ----
#' @noRd
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
