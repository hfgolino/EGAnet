#' @title Traditional Fit Metrics for Networks
#'
#' @description Computes several traditional fit metrics for networks including
#'
#' \itemize{
#'
#' \item chi-square (\eqn{\chi^2})
#'
#' \item root mean square error of approximation (RMSEA) with confidence intervals
#'
#' \item confirmatory fit index (CFI)
#'
#' \item Tucker-Lewis inded (TLI)
#'
#' \item standardized root mean residual (SRMR)
#'
#' \item log-likelihood
#'
#' \item Akaike's information criterion (AIC)
#'
#' \item Bayesian information criterion (BIC)
#'
#' }
#'
#' @param network Matrix or data frame.
#' A p by p sqaure network matrix
#'
#' @param n Numeric (length = 1).
#' Sample size
#'
#' @param S Matrix or data frame.
#' A p by p sqaure zero-order correlation matrix corresponding
#' with the input \code{network}
#'
#' @param ci Numeric (length = 1).
#' Confidence interval for RMSEA
#'
#' @return Returns a named vector of fit statistics
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain correlation matrix
#' S <- auto.correlate(wmt)
#'
#' # EBICglasso (default for EGA functions)
#' glasso_network <- network.estimation(
#'   data = wmt, model = "glasso"
#' )
#'
#' # Obtain fit (expects continuous variables!)
#' network.fit(network = glasso_network, n = nrow(wmt), S = S)
#' # Scaled metrics are not yet available for
#' # dichotomous or polytomous data!
#'
#' @references
#' Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017).
#' Generalized network psychometrics: Combining network and latent variable models.
#' \emph{Psychometrika}, \emph{82}(4), 904–927.
#'
#' @export
#'
# Compute network fit statistics ----
# Updated 19.11.2025
network.fit <- function(network, n, S, ci = 0.95)
{

  # Obtain lower triangle indices
  lower_triangle <- lower.tri(network)

  # Obtain model-implied zero-order correlation matrix
  R <- pcor2cor(network)

  # Compute number of parameters
  p <- dim(network)[2]
  zero_parameters <- p * (p - 1) / 2
  model_parameters <- sum(network[lower_triangle] != 0)

  # Compute baseline
  baseline <- diag(1, nrow = p, ncol = p)
  baseline_ML <- log(det(baseline)) + sum(diag(S %*% solve(baseline))) - log(det(S)) - p
  baseline_chi_square <- n * baseline_ML
  baseline_tli <- baseline_chi_square / zero_parameters

  # Compute traditional SEM measures
  loglik_ML <- log(det(R)) + sum(diag(S %*% solve(R))) - log(det(S)) - p
  chi_square <- n * loglik_ML
  df <- zero_parameters - model_parameters
  chi_max <- max(chi_square - df, 0)
  nDF <- n * df
  rmsea_null <- nDF * 0.0025 # 0.05^2

  # log-likelihood
  loglik <- -(n / 2) * (p * log(2 * pi) + log(det(R)) + sum(diag(S %*% solve(R))))
  # Assumes no mean structure (or that all means are equal to zero)

  # Obtain RMSEA confidence intervals
  rmsea_cis <- rmsea_ci(chi_square, df, n, nDF, ci)

  # Get fit indices
  fit_indices <- c(
    # Traditional fit measures
    chisq = chi_square, df = df, chisq.p.value = 1 - pchisq(chi_square, df = df),
    RMSEA = sqrt(chi_max / nDF),
    rmsea_cis,
    RMSEA.p.value = 1 - pchisq(chi_max, df = df, ncp = rmsea_null),
    CFI = 1 - (chi_max / max(baseline_chi_square - zero_parameters, 0)),
    TLI = (baseline_tli - (chi_square / df)) / (baseline_tli - 1),
    SRMR = sqrt(mean((R[lower_triangle] - S[lower_triangle])^2)),
    # Gaussian log-likelihood measures
    logLik = loglik,
    AIC = -2 * loglik + 2 * model_parameters,
    BIC = -2 * loglik + model_parameters * log(n)
  )

  # Rename confidence intervals
  names(fit_indices)[names(fit_indices) %in% c("lower", "upper")] <-
    paste("RMSEA", format_integer(ci * 100, 1), c("lower", "upper"), sep = ".")

  # Return log-likelihood
  return(fit_indices)

}

#' @noRd
# RMSEA Confidence Intervals ----
# Follows {lavaan} version 0.6.19
# Updated 01.11.2024
rmsea_ci <- function(chi_square, df, n, nDF, ci)
{

  # Set up CI
  lower_ci <- 1 - (1 - ci) / 2
  upper_ci <- 1 - lower_ci

  # Internal function for finding RMSEA confidence intervals
  # (same as {lavann} version 0.6.19)
  find_lambda <- function(lambda, ci){
    pchisq(chi_square, df = df, ncp = lambda) - ci
  }

  # Find lower bound
  if(df < 1 || find_lambda(0, lower_ci) < 0){
    rmsea_lower <- 0
  }else{

    # Try uniroot
    lambda <- try(
      uniroot(f = find_lambda, lower = 0, upper = chi_square, ci = lower_ci)$root,
      silent = TRUE
    )

    # Determine if there was an error
    rmsea_lower <- swiftelse(is(lambda, "try-error"), NA, sqrt(lambda / nDF))

  }

  # Find upper bound
  N_RMSEA <- max(n, chi_square * 4)
  if(df < 1 || find_lambda(N_RMSEA, upper_ci) > 0 || find_lambda(0, upper_ci) < 0){
    rmsea_upper <- 0
  }else{

    # Try uniroot
    lambda <- try(
      uniroot(f = find_lambda, lower = 0, upper = N_RMSEA, ci = upper_ci)$root,
      silent = TRUE
    )

    # Determine if there was an error
    rmsea_upper <- swiftelse(is(lambda, "try-error"), NA, sqrt(lambda / nDF))

  }

  # Return confidence interval
  return(c(lower = rmsea_lower, upper = rmsea_upper))

}
