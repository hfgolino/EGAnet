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
#' \code{\link[psych]{fa}}
#'
#' @examples
#' # Get depression data
#' data <- na.omit(depression[,24:44])
#'
#' # Estimate EGM (using EGA)
#' EGM.compare(data)
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#
# Compare EGM to EFA ----
# Updated 23.10.2024
EGM.compare <- function(data, constrained = FALSE, rotation = "geominQ", ...)
{

  # Check data and structure
  data <- EGM.compare_errors(data, ...)

  # Obtain dimensions of the data
  dimensions <- dim(data)

  # Estimate EGM
  egm <- EGM(data, constrained = constrained, ...)

  # Set communities
  communities <- unique_length(egm$EGA$wc)

  # Obtain EFA
  efa <- get_factor_results(
    output = psych::fa(
      egm$EGA$correlation, n.obs = dimensions[2],
      nfactors = communities, rotate = "none", ...
    ), rotation = rotation, egm = egm,
    dimensions = dimensions, ...
  )

  # Set up results
  results <- list(
    EGM = egm,
    EFA = efa,
    likelihood = as.data.frame(
      rbind(EGM = egm$model$optimized$fit, EFA = efa$fit)
    )
  )

  # Set class
  class(results) <- "EGM.compare"

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
# Updated 10.10.2024
print.EGM.compare <- function(x, ...)
{

  # Control exponential notation of numbers
  user_scipen <- options("scipen")$scipen

  # Set option
  options(scipen = 999)

  # Find smallest difference
  absolute <- abs(x$likelihood)
  smallest_difference <- min(abs(absolute[1,] - absolute[2,]))

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
  rounded <- round(x$likelihood, smallest_decimal)

  # Get minimum
  minimums <- nvapply(rounded, which.min)

  # Replace log-likelihood with maximum
  minimums["logLik"] <- swiftelse(minimums["logLik"] == 1, 2, 1)

  # Add lowest of each column to bottom row
  rounded["best",] <- c("EGM", "EFA")[minimums]

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
# Updated 31.10.2024
get_factor_results <- function(output, rotation, egm, dimensions, ...)
{

  # Check for errors in rotation
  # (this code and the following rotation code is in `net.loads`)
  rotation <- rotation_errors(rotation)

  # If rotation exists, then obtain it
  rotation_FUN <- get(rotation, envir = asNamespace("GPArotation"))

  # Get ellipse arguments
  ellipse <- list(...)

  # Get arguments for function
  rotation_ARGS <- obtain_arguments(rotation_FUN, ellipse)

  # Supply loadings
  rotation_ARGS$A <- output$loadings[]

  # Set default arguments for rotations
  rotation_ARGS <- rotation_defaults(rotation, rotation_ARGS, ellipse)

  # Perform rotations
  rotation_OUTPUT <- do.call(rotation_FUN, rotation_ARGS)

  # Align rotated loadings
  aligned_output <- fungible::faAlign(
    F1 = egm$model$optimized$loadings,
    F2 = rotation_OUTPUT$loadings,
    Phi2 = rotation_OUTPUT$Phi
  )

  # Obtain model-implied correlations
  implied_R <- aligned_output$F2 %*% tcrossprod(
    aligned_output$Phi, aligned_output$F2
  ); diag(implied_R) <- 1

  # Add parameters to EFA
  output$loadings <- aligned_output$F2
  output$factor_correlations <- aligned_output$Phi2
  output$implied <- list(R = implied_R, P = cor2pcor(implied_R))

  # Compute likelihood for EFA
  output$fit <- c(
    R.srmr = srmr(egm$EGA$correlation, implied_R),
    P.srmr = srmr(cor2pcor(egm$EGA$correlation), output$implied$P),
    likelihood(
      n = dimensions[1], p = dimensions[2], R = output$implied$R,
      S = egm$EGA$correlation, loadings = output$factor_correlations
    ),
    TEFI_adj = tefi(output$implied$R, structure = egm$EGA$wc)$VN.Entropy.Fit - compute_tefi_adjustment(
      output$loadings, output$factor_correlations
    )
  )

  # Return updated output
  return(output)

}
