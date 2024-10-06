#' Compare \code{\link[EGAnet]{EGM}} with Exploratory Factor Analysis
#'
#' @description Estimates an \code{\link[EGAnet]{EGM}} based on \code{\link[EGAnet]{EGA}} and
#' uses the number of communities as the number of dimensions in exploratory factor analysis
#' (EFA) using \code{\link[psych]{fa}}
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{community.unidimensional}},
#' \code{\link[EGAnet]{EGA}}, and
#' \code{\link[psych]{fa}}
#'
#' @examples
#' # Compare EGM with EFA
#' EGM.compare(wmt2[,7:24])
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#
# Compare EGM to EFA ----
# Updated 06.10.2024
EGM.compare <- function(data, rotation = "geominQ", ...)
{

  # Check data and structure
  data <- EGM.compare_errors(data, ...)

  # Obtain dimensions of the data
  dimensions <- dim(data)

  # Estimate EGM
  egm <- EGM(data, ...)

  # Set up EFA
  efa <- psych::fa(
    egm$EGA$correlation, n.obs = dimensions[2],
    nfactors = egm$EGA$n.dim, rotate = "none"
  )

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
  rotation_ARGS$A <- efa$loadings[,seq_len(egm$EGA$n.dim)]

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

  # Compute likelihood for EFA
  efa_fit <- c(
    R.srmr = srmr(egm$EGA$correlation, implied_R),
    P.srmr = srmr(cor2pcor(egm$EGA$correlation), cor2pcor(implied_R)),
    likelihood(
      n = dimensions[1], p = dimensions[2], R = implied_R,
      S = egm$EGA$correlation, loadings = aligned_output$F2
    ),
    TEFI = tefi(implied_R, structure = egm$EGA$wc)$VN.Entropy.Fit
  )

  # Set up results
  results <- list(
    EGM = egm,
    EFA = efa,
    likelihood = as.data.frame(rbind(EGM = egm$model$optimized$fit, EFA = efa_fit))
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
