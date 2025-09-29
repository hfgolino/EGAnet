#' @title Intra- and Inter-individual \code{\link[EGAnet]{dynEGA}}
#'
#' @description A wrapper function to estimate both intra-individual
#' (\code{level = "individual"}) and inter-individual (\code{level = "population"})
#' structures using \code{\link[EGAnet]{dynEGA}}. This wrapper enforces estimation
#' at both levels simultaneously and simplifies the interface for users who want
#' both outputs in a single call.
#'
#' @param data Matrix or data frame in long format.
#' Each row \emph{t} should represent observations for all variables at time point
#' \emph{t} for a given participant. The next row, \emph{t + 1}, represents the next
#' measurement occasion for the same participant, and so forth. After one participant’s
#' series ends, the next participant’s data should immediately follow in the same format.
#'
#' \code{data} must contain an ID column (named \code{"ID"} by default or specified via \code{id});
#' otherwise, it is assumed the entire dataset represents the population.
#'
#' If group-level analyses are desired, \code{data} should contain a Group column
#' (named \code{"Group"} by default or specified via \code{group}).
#'
#' A measurement occasion variable is not required and should be removed prior to analysis.
#'
#' @param id Numeric or character (length = 1).
#' Number or name of the column identifying each individual.
#' Defaults to \code{NULL}.
#'
#' @param n.embed Numeric (length = 1).
#' Defaults to \code{5}.
#' Number of embedded dimensions (the number of consecutive observations used in the
#' \code{\link[EGAnet]{Embed}} function). For example, \code{n.embed = 5} uses five
#' consecutive observations to estimate a single derivative.
#'
#' @param tau Numeric (length = 1). Defaults to \code{1}.
#' Number of observations to offset successive embeddings in
#' the \code{\link[EGAnet]{Embed}} function.
#'
#' @param delta Numeric (length = 1). Defaults to \code{1}.
#' The time between successive observations in the time series (i.e., lag).
#'
#' @param use.derivatives Numeric (length = 1). Defaults to \code{1}.
#' Order of the derivative to be used. Options:
#' \itemize{
#'   \item \code{0} — No derivatives; moving average
#'   \item \code{1} — First-order derivatives; velocity (rate of change)
#'   \item \code{2} — Second-order derivatives; acceleration (rate of change of rate of change)
#' }
#'
#' @param corr Character. Correlation method (default = \code{"auto"}).
#' Options: \code{"auto"}, \code{"cor_auto"}, \code{"pearson"}, \code{"spearman"}.
#' See \code{\link[EGAnet]{dynEGA}} for details.
#'
#' @param na.data Character. Missing data handling (default = \code{"pairwise"}).
#' Options: \code{"pairwise"}, \code{"listwise"}.
#'
#' @param model Character. Network estimation method (default = \code{"glasso"}).
#' Options: \code{"BGGM"}, \code{"glasso"}, \code{"TMFG"}.
#'
#' @param algorithm Character. Community detection method (default = \code{"walktrap"}).
#' Options: \code{"leiden"}, \code{"louvain"}, \code{"walktrap"}.
#'
#' @param uni.method Character. Method for testing unidimensionality (default = \code{"louvain"}).
#' Options: \code{"expand"}, \code{"LE"}, \code{"louvain"}.
#'
#' @param ncores Numeric. Number of cores for parallel processing.
#' Defaults to half of available cores (\code{ceiling(parallel::detectCores()/2)}).
#' Set \code{ncores = 1} to disable parallelization.
#'
#' @param verbose Logical. Should progress messages be displayed? Default = \code{TRUE}.
#'
#' @param optimization Logical. Default = \code{FALSE}.
#' If \code{TRUE}, the function performs automated optimization of the embedding
#' dimension (\code{n.embed}) for each individual by minimizing TEFI values.
#' Population-level structure is then estimated from the optimized individual embeddings.
#'
#' @param n.embed.all.ind Vector of candidate \code{n.embed} values for optimization.
#' Default = \code{3:25}. For each individual with \eqn{T_i} time points, the optimizer
#' restricts to \code{n.embed <= T_i - 3}. If none are valid, that individual is skipped.
#'
#' @param jitter.eps Numeric scalar. Default = \code{1e-3}.
#' Small Gaussian noise added to any zero-variance columns to prevent estimation
#' failures. This preserves the overall structure of the data but avoids singular
#' covariance matrices.
#'
#' @param ... Additional arguments passed to \code{\link[EGAnet]{dynEGA}} and
#' related functions.
#'
#' @return A list with \code{dynEGA} results at the population and individual levels.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' # Simulated example
#' dyn.ega1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA,
#'   n.embed = 5, tau = 1, delta = 1,
#'   id = 25, use.derivatives = 1,
#'   ncores = 2, corr = "pearson",
#'   optimization = TRUE
#' )
#' }
#'
#' @seealso \code{\link[EGAnet]{dynEGA}}, \code{\link[EGAnet]{plot.EGAnet}}
#' @export
#'
# Intra- and Interindividual dynEGA
# Updated 24.9.2024
dynEGA.ind.pop <- function(
  # `dynEGA` arguments
  data,  id = NULL,
  n.embed = 5, tau = 1, delta = 1, use.derivatives = 1,
  # `EGA` arguments
  corr = c("auto", "cor_auto", "pearson", "spearman"),
  na.data = c("pairwise", "listwise"),
  model = c("BGGM", "glasso", "TMFG"),
  algorithm = c("leiden", "louvain", "walktrap"),
  uni.method = c("expand", "LE", "louvain"),
  ncores, verbose = TRUE, optimization = FALSE, ...
){
  # Use `dynEGA` (input handling occurs inside `dynEGA`)
  return(
    dynEGA(
      data = data, id = id, n.embed = n.embed, tau = tau,
      delta = delta, use.derivatives = use.derivatives,
      level = c("population", "individual"),
      corr = corr, na.data = na.data,
      model = model, algorithm = algorithm, uni.method = uni.method,
      ncores = ncores, verbose = verbose, optimization = optimization, ...
    )
  )

}
