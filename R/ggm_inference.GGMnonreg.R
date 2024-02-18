#' @title \code{ggm_inference} from \code{GGMnonreg} 1.0.0
#'
#' @description This minimal implementation of the \code{ggm_inference} function in
#' \code{GGMnonreg} is modified to operate with \code{\link{EGAnet}} with maximum efficiency
#' (i.e., significantly faster computation for categorical correlations). The
#' goal is to maintain identical functionality to the original function
#' but faster (for more details, see Williams et al., 2019)
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Raw data must be provided (i.e., no correlation matrices)
#'
#' @param n \code{NULL}.
#' Not actually used but makes it easier for general functionality
#' in the package
#'
#' @param p.value Numeric (length = 1).
#' Significance level for edge inclusion.
#' Defaults to \code{0.10} based on Williams and Rodriguez (2022)
#'
#' @param bootstrap Boolean (length = 1).
#' Defaults to \code{TRUE}.
#' Can only be \code{FALSE} for Pearson's
#' (\code{corr = "pearson"}) and Spearman's
#' (\code{corr = "spearman"}) correlations
#'
#' @param iter Numeric (length = 1).
#' Number of replica samples to generate from the bootstrap analysis.
#' Defaults to \code{100} (minimum allowed).
#' Default in \code{GGMnonreg} is \code{1000}
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)
#'
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to
#' compute correlations. Arguments can be passed along to the function
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for
#' all variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is
#' computed for all variables regardless of categories
#'
#' }
#'
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available
#' cases between two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete
#' cases in the dataset
#'
#' }
#'
#' @param returnAllResults Boolean (length = 1).
#' Not actually used but makes it easier for general functionality
#' in the package
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ...
#' Not actually used but makes it easier for general functionality
#' in the package
#'
#' @return A partial correlation matrix
#'
#' @references
#' \strong{Seminal paper on Maximum Likelihood Estimation} \cr
#' Williams, D. R., Rhemtulla, M., Wysocki, A. C., & Rast, P. (2019).
#' On nonregularized estimation of psychological networks.
#' \emph{Multivariate Behavioral Research}, \emph{54}(5), 719-750.
#'
#' \strong{Original implementation (GGMnonreg package)} \cr
#' Williams, D. R. (2021).
#' GGMnonreg: Non-regularized Gaussian Graphical Models in R.
#' \emph{Journal of Open Source Software}, \emph{6}(67), 3308.
#'
#' \strong{(Non-)regularization network generalizability} \cr
#' Williams, D. R., & Rodriguez, J. E. (2022).
#' Why overfitting is not (usually) a problem in partial correlation networks.
#' \emph{Psychological Methods}, \emph{27}(5), 822-840.
#'
#' @author Donald R. Williams; for maintanence,
#' Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen at gmail.com>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Compute non-regularized network
#' nonreg <- ggm_inference.GGMnonreg(data = wmt)}
#'
#' @export
#'
# Computes non-regularized network using Maximum Likelihood ----
# Updated 18.02.2024
ggm_inference.GGMnonreg <- function(
    data, n = NULL, p.value = 0.10, bootstrap = TRUE, iter = 100,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    returnAllResults = FALSE,
    verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", ggm_inference.GGMnonreg)
  na.data <- set_default(na.data, "pairwise", auto.correlate)

  # Argument errors (return data in case of tibble)
  data <- ggm_inference.GGMnonreg_errors(data, p.value, bootstrap, iter, verbose, ...)

  # Check for correlation conflicts with bootstrap
  if(!bootstrap && !corr %in% c("pearson", "spearman")){

    # Send error
    stop(
      "With `bootstrap = FALSE`, argument 'corr' must either be \"pearson\" or \"spearman\"",
      call. = FALSE
    )

  }

  # Get dimensions
  dimensions <- dim(data)
  dimension_names <- dimnames(data)
  row_sequence <- seq_len(dimensions[1])

  # Check for bootstrap
  if(bootstrap){

    # Get bootstrapped partial correlations
    bootstrap_partial <- lapply(
      seq_len(iter), function(iteration){

        # Return partial correlations
        return(
          cor2pcor(
            obtain_sample_correlations(
              data = data[shuffle_replace(row_sequence, size = dimensions[1]),],
              n = dimensions[1], corr = corr, na.data = na.data, verbose = verbose
            )$correlation_matrix
          )
        )

      }
    )

    # Set confidence interval
    ci <- p.value / 2

    # Get means of network values
    boot_means <- symmetric_matrix_lapply(bootstrap_partial, mean)
    boot_low_ci <- symmetric_matrix_lapply(bootstrap_partial, quantile, probs = ci)
    boot_high_ci <- symmetric_matrix_lapply(bootstrap_partial, quantile, probs = 1 - ci)

    # Get significant indices
    boot_indices <- !(boot_low_ci < 0 & boot_high_ci > 0)

    # Set up network
    network <- matrix(
      0, nrow = dimensions[2], ncol = dimensions[2],
      dimnames = list(dimension_names[[2]], dimension_names[[2]])
    )

    # Get network
    network[boot_indices] <- boot_means[boot_indices]

  }else{

    # Compute partial correlations
    partial <- cor2pcor(cor(data, use = na.data, method = corr))

    # Convert to Fisher's z
    partial_z <- abs(0.5 * log((1 + partial) / (1 - partial)))

    # Get test statistic
    statistic <- partial_z / sqrt(
      swiftelse(corr == "pearson", 1, 1.06)  /
      (dimensions[1] - dimensions[2] - 5)
    )

    # Get network
    network <- partial * ((pnorm(statistic, lower.tail = FALSE) * 2) < p.value)

  }

  # Set methods in attributes
  attr(network, "methods") <- list(
    corr = corr, na.data = na.data,
    p.value = p.value, bootstrap = bootstrap,
    iter = iter
  )

  # Return network
  return(network)

}

# Bug Checking ----
# ## Basic input
# data = wmt2[,7:24]; alpha = 0.10; bootstrap = TRUE
# iter = 100; corr = "auto"; na.data = "pairwise"
# verbose = FALSE

#' @noRd
# Errors ----
# Updated 18.02.2024
ggm_inference.GGMnonreg_errors <- function(data, p.value, bootstrap, iter, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "ggm_inference.GGMnonreg")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'p.value' errors
  length_error(p.value, 1, "ggm_inference.GGMnonreg")
  typeof_error(p.value, "numeric", "ggm_inference.GGMnonreg")
  range_error(p.value, c(0.01, 0.99), "ggm_inference.GGMnonreg")

  # 'bootstrap' errors
  length_error(bootstrap, 1, "ggm_inference.GGMnonreg")
  typeof_error(bootstrap, "logical", "ggm_inference.GGMnonreg")

  # 'iter' errors
  length_error(iter, 1, "ggm_inference.GGMnonreg")
  typeof_error(iter, "numeric", "ggm_inference.GGMnonreg")
  range_error(iter, c(100, Inf), "ggm_inference.GGMnonreg")

  # 'verbose' errors
  length_error(verbose, 1, "ggm_inference.GGMnonreg")
  typeof_error(verbose, "logical", "ggm_inference.GGMnonreg")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return usable data in case of tibble
  return(data)

}

