#' Compare \code{\link[EGAnet]{EGM}} with EFA
#'
#' @description Estimates an \code{\link[EGAnet]{EGM}} based on \code{\link[EGAnet]{EGA}} and
#' uses the number of communities as the number of dimensions in exploratory factor analysis
#' (EFA) using \code{\link[psych]{fa}}
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param constrained Boolean (length = 1).
#' Whether memberships of the communities should
#' be added as a constraint when optimizing the network loadings.
#' Defaults to \code{FALSE} to freely estimate each loading similar to
#' exploratory factor analysis.
#'
#' \emph{Note: This default differs from} \code{\link[EGAnet]{EGM}}\emph{.
#' Constraining loadings puts EGM at a deficit relative to EFA and therefore
#' biases the comparability between the methods. It's best to leave the
#' default of unconstrained when using this function.}
#'
#' @param rotation Character.
#' A rotation to use to obtain a simpler structure for EFA.
#' For a list of rotations, see \code{\link[GPArotation]{rotations}} for options.
#' Defaults to \code{"geominQ"}
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{community.unidimensional}},
#' \code{\link[EGAnet]{EGA}},
#' \code{\link[EGAnet]{EGM}},
#' \code{\link[EGAnet]{net.loads}}, and
#' \code{\link[lavaan]{efa}}
#'
#' @examples
#' # Get depression data
#' data <- na.omit(depression[,24:44])
#'
#' # Compare EGM (using EGA) with EFA
#' \dontrun{
#' results <- EGM.compare(data)
#'
#' # Print summary
#' summary(results)}
#'
#'
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#
# Compare EGM to EFA ----
# Updated 02.11.2024
EGM.compare <- function(data, constrained = FALSE, rotation = "geomin", ...)
{

  # Obtain ellipse
  ellipse <- list(...)

  # Check data and structure
  data <- EGM.compare_errors(data, ...)

  # Obtain dimensions of the data
  dimensions <- dim(data)
  dimension_names <- dimnames(data)

  # Estimate EGM
  egm <- EGM(data, constrained = constrained, ...)

  # Set communities
  communities <- unique_length(egm$EGA$wc)

  # Obtain the number of categories for each variables
  categories <- data_categories(data)

  # Determine categorical variables
  categorical_variables <- categories <= swiftelse(
    "ordinal.categories" %in% names(ellipse), ellipse$ordinal.categories, 7
  )

  # Set up arguments
  efa_ARGS <- obtain_arguments(
    FUN = lavaan::efa,
    FUN.args = c(
      list(
        data = data, nfactors = communities,
        ov.names = dimension_names[[2]],
        rotation = rotation,
        rotation.args = list(
          geomin.epsilon = switch(
            as.character(communities), `2` = 0.0001, `3` = 0.001, 0.01
          )
        ),
        estimator = swiftelse(any(categorical_variables), "WLSMV", "ML"),
        ordered = dimension_names[[2]][categorical_variables]
      ), ellipse
    )
  )

  # Obtain EFA
  efa <- silent_call(
    get_factor_results(
      output = do.call(what = lavaan::efa, args = efa_ARGS),
      egm = egm, dimensions = dimensions, ...
    )
  )

  # Set up results
  results <- list(
    EGM = egm,
    EFA = efa,
    fit = data.frame(
      EGM = egm$model$optimized$fit,
      EFA = efa$fit
    )
  )

  # Set class
  class(results) <- "EGM.compare"

  # Send warning about categorical data
  if(any(categorical_variables)){

    warning(
      paste0(
        "Categorical variables were detected in your data.\n\n",
        "`EGM.compare` does not currently support scaled fit measures ",
        "for non-continuous data. ",
        "Fit statistics are reported based on maximum likelihood and are ",
        "likely to ", styletext("underestimate", "italics"), " fit.\n\n",
        "Use caution when interpreting the fit statistics."
      ), call. = FALSE
    )

  }

  # Return results
  return(results)

}

#' @noRd
# EGM Compare Errors ----
# Updated 04.10.2023
EGM.compare_errors <- function(data, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "EGA")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, FALSE)
  }

  # Return data in case of tibble
  return(data)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 01.11.2024
print.EGM.compare <- function(x, ...)
{

  # Control exponential notation of numbers
  user_scipen <- options("scipen")$scipen

  # Set option
  options(scipen = 999)

  # Find smallest difference
  absolute <- abs(x$fit)
  avoid_ps <- which(
    row.names(x$fit) %in%
    c("df", "chisq.p.value", "RMSEA.95.lower", "RMSEA.p.value")
  )
  smallest_difference <- min(
    abs(absolute[-avoid_ps, "EGM"] - absolute[-avoid_ps, "EFA"])
  )

  # Get decimal split
  decimal_split <- strsplit(
    as.character(smallest_difference), split = "\\."
  )[[1]][[2]]

  # Get number of leading zeros
  zeros <- strsplit(decimal_split, split = "")[[1]] == "0"

  # Find smallest digit after decimal
  smallest_decimal <- swiftelse(
    smallest_difference > 1, 1,
    setdiff(seq_along(zeros), which(zeros))[1]
  )

  # Round likelihood to smallest decimal
  rounded <- round(x$fit, smallest_decimal)

  # Get minimum
  minimums <- apply(rounded[-2,], 1, which.min)

  # Replace maximums
  maximums <- c("CFI", "TLI", "logLik")

  # Replace log-likelihood with maximum
  minimums[maximums] <- swiftelse(minimums[maximums] == 1, 2, 1)

  # Add lowest of each column to bottom row
  rounded[-2, "best"] <- c("EGM", "EFA")[minimums]

  # Print to smallest decimal
  print(rounded)

  # Set back user's options
  options(scipen = user_scipen)

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 09.10.2024
summary.EGM.compare <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# Get factor results ----
# Updated 01.11.2024
get_factor_results <- function(output, rotation, egm, dimensions, ...)
{

  # Get standardized parameters
  parameters <- lavaan::inspect(output[[1]], what = "std")

  # Align rotated loadings
  aligned_output <- fungible::faAlign(
    F1 = egm$model$optimized$loadings,
    F2 = parameters$lambda,
    Phi2 = parameters$psi
  )

  # Obtain model-implied correlations
  implied_R <- aligned_output$F2 %*% tcrossprod(
    aligned_output$Phi, aligned_output$F2
  ); diag(implied_R) <- 1

  # Add parameters to EFA
  output$loadings <- aligned_output$F2
  output$factor_correlations <- aligned_output$Phi2
  output$implied <- list(R = implied_R, P = cor2pcor(implied_R))

  # Get {lavaan} fit measures
  lavaan_fit <- lavaan::fitmeasures(output[[1]])

  # Compute likelihood for EFA
  output$fit <- fit(
    n = dimensions[1], p = dimensions[2], R = output$implied$R,
    S = egm$EGA$correlation, loadings = output$loadings,
    correlations = output$factor_correlations,
    structure = egm$EGA$wc, ci = 0.95,
    scaling_factor = swiftelse(
      "chisq.scaling.factor" %in% names(lavaan_fit),
      lavaan_fit[["chisq.scaling.factor"]], NULL
    )
  )

  # Return updated output
  return(output)

}

#' @noRd
# Compute fit metrics ----
# Updated 01.11.2024
fit <- function(n, p, R, S, loadings, correlations, structure, ci, scaling_factor = NULL)
{

  # Get number of communities
  m <- dim(loadings)[2]

  # Total number of parameters
  zero_parameters <- p * (p - 1) / 2
  loading_parameters <- p * m - sum(loadings == 0)
  correlation_parameters <- ((m * (m - 1)) / 2)
  model_parameters <- loading_parameters + correlation_parameters

  # Baseline
  baseline <- diag(1, nrow = p, ncol = p)
  baseline_ML <- log(det(baseline)) + sum(diag(S %*% solve(baseline))) - log(det(S)) - p
  baseline_chi_square <- n * baseline_ML
  baseline_tli <- baseline_chi_square / zero_parameters

  # Compute traditional SEM measures
  loglik_ML <- log(det(R)) + sum(diag(S %*% solve(R))) - log(det(S)) - p
  chi_square <- n * loglik_ML
  df <- zero_parameters - loading_parameters + correlation_parameters
  chi_max <- max(chi_square - df, 0)
  nDF <- n * df
  rmsea_null <- nDF * 0.0025 # 0.05^2

  # log-likelihood
  loglik <- -(n / 2) * (p * log(2 * pi) + log(det(R)) + sum(diag(S %*% solve(R))))
  # Assumes no mean structure (or that all means are equal to zero)

  # Compute TEFI
  TEFI <- tefi(R, structure = structure)$VN.Entropy.Fit

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
    SRMR = srmr(S, R),
    # Gaussian log-likelihood measures
    logLik = loglik,
    AIC = -2 * loglik + 2 * model_parameters,
    BIC = -2 * loglik + model_parameters * log(n),
    # TEFI measures
    TEFI = TEFI,
    TEFI.adj = TEFI - (-2 * log(model_parameters) + mean(abs(correlations)))
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
