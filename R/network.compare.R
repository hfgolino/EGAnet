#' @title Compares Network Structures Using Permutation
#'
#' @description A permutation implementation to determine statistical
#' significance of whether the network structures are different from one another
#'
#' @param base Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' First dataset
#'
#' @param comparison Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Second dataset
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
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
#' Arguments can be passed along to the function
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
#'
#' }
#'
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#'
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"BGGM"} --- Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{?BGGM::estimate} for more details
#'
#' \item \code{"glasso"} --- Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details
#'
#' \item \code{"TMFG"} --- Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details
#'
#' }
#'
#' @param iter Numeric (length = 1).
#' Number of permutations to perform.
#' Defaults to \code{1000} (recommended)
#'
#' @param ncores Numeric (length = 1).
#' Number of cores to use in computing results.
#' Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results.
#' Set for reproducible results.
#' See \href{https://r-ega.net/articles/reproducibility-prng.html}{Reproducibility and PRNG}
#' for more details on random number generation in \code{\link{EGAnet}}
#'
#' @param ... Additional arguments that can be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{EGA}}
#'
#' @return Returns data frame with row names of each measure, empirical value (\code{statisticl}), and p-value based on the permutation test (\code{p.value})
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Set groups (if necessary)
#' groups <- rep(1:2, each = nrow(wmt) / 2)
#'
#' # Groups
#' group1 <- wmt[groups == 1,]
#' group2 <- wmt[groups == 2,]
#'
#' \dontrun{
#' results <- network.compare(group1, group2)}
#'
#' @references
#' \strong{Frobenius Norm} \cr
#' Ulitzsch, E., Khanna, S., Rhemtulla, M., & Domingue, B. W. (2023).
#' A graph theory based similarity metric enables comparison of subpopulation psychometric networks.
#' \emph{Psychological Methods}.
#'
#' \strong{Jensen-Shannon Similarity (1 - Distance)} \cr
#' De Domenico, M., Nicosia, V., Arenas, A., & Latora, V. (2015).
#' Structural reducibility of multilayer networks.
#' \emph{Nature Communications}, \emph{6}(1), 1–9.
#'
#' \strong{Total Network Strength} \cr
#' van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023).
#' Comparing network structures on three aspects: A permutation test.
#' \emph{Psychological Methods}, \emph{28}(6), 1273–1285.
#'
#'
#' @export
#'
# Perform permutations for network structures ----
# Updated 11.07.2024
network.compare <- function(
    base, comparison,
    # EGA arguments
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    # Permutation arguments
    iter = 1000, ncores, verbose = TRUE, seed = NULL,
    ...
)
{

  # Experimental warning
  experimental("network.compare")

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", EGA.estimate)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)

  # Set cores
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}

  # Check for input errors
  error_return <- network.compare_errors(base, comparison, iter, verbose, seed, ...)

  # Get ellipse
  ellipse <- list(...)

  # Check for seed
  if(!is.null(seed)){
    seeds <- reproducible_seeds(iter, seed)
  }else{

    # Set all seeds to zero (or random)
    seeds <- rep(0, iter)

    # Check for external suppression (from `invariance`)
    if(!"suppress" %in% names(ellipse) || !ellipse$suppress){
      message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results")
    }

  }

  # Get empirical networks
  base_empirical_network <- EGA(
    base, corr = corr, na.data = na.data,
    model = model, plot.EGA = FALSE, ...
  )$network
  comparison_empirical_network <- EGA(
    comparison, corr = corr, na.data = na.data,
    model = model, plot.EGA = FALSE, ...
  )$network

  # Get empirical estimates
  empirical_values = abs(
    c(
      "Frobenius" = frobenius(base_empirical_network, comparison_empirical_network),
      "JSS" = 1 - jsd(base_empirical_network, comparison_empirical_network),
      "Total Strength" = sum(colSums(abs(base_empirical_network), na.rm = TRUE), na.rm = TRUE) -
        sum(colSums(abs(comparison_empirical_network), na.rm = TRUE), na.rm = TRUE)
    )
  )

  # Create combined dataset
  combined <- rbind(base, comparison)

  # Set up indices
  combined_index <- nrow_sequence(combined)
  base_length <- dim(base)[1]

  # Perform permutations
  permutated_values <- do.call(
    rbind, parallel_process(
      iterations = iter, datalist = seeds, FUN = function(seed, ...){

        # Get shuffled indices
        base_shuffled <- shuffle(combined_index, size = base_length, seed = seed)

        # Get permutated networks
        base_network <- EGA(
          combined[base_shuffled,], corr = corr, na.data = na.data,
          model = model, plot.EGA = FALSE, ...
        )$network
        comparison_network <- EGA(
          combined[-base_shuffled,], corr = corr, na.data = na.data,
          model = model, plot.EGA = FALSE, ...
        )$network

        # Return permutated estimates
        return(
          c(
            "Frobenius" = frobenius(base_network, comparison_network),
            "JSS" = 1 - jsd(base_network, comparison_network),
            "Total Strength" = sum(colSums(abs(base_network), na.rm = TRUE), na.rm = TRUE) -
              sum(colSums(abs(comparison_network), na.rm = TRUE), na.rm = TRUE)
          )
        )

      }, ncores = ncores, progress = verbose, ...
    )
  )

  # Return statistics
  return(
    t(data.frame(
      "statistic" = empirical_values,
      "p.value" = c(
        mean(permutated_values[,1] <= empirical_values[1]),
        mean(permutated_values[,2] <= empirical_values[2]),
        mean(abs(permutated_values[,3]) >= empirical_values[3])
      ),
      "M_permutated" = colMeans(permutated_values),
      "SD_permutated" = apply(permutated_values, 2, sd)
    ))
  )

}

# Bug Checking ----
# wmt <- wmt2[-1,7:24]
# groups <- rep(1:2, each = nrow(wmt) / 2)
# base <- wmt[groups == 1,]; comparison <- wmt[groups == 2,]
# corr = "auto"; na.data = "pairwise"; model = "glasso"
# iter = 1000; ncores = 8; verbose = TRUE; seed = NULL

#' @noRd
# Errors ----
# Updated 10.07.2024
network.compare_errors <- function(base, comparison, iter, verbose, seed, ...)
{

  # 'base' errors
  object_error(base, c("matrix", "data.frame", "tibble"), "network.compare")

  # Check for tibble
  if(get_object_type(base) == "tibble"){
    base <- as.data.frame(base)
  }

  # 'comparison' errors
  object_error(comparison, c("matrix", "data.frame", "tibble"), "network.compare")

  # Check for tibble
  if(get_object_type(comparison) == "tibble"){
    comparison <- as.data.frame(comparison)
  }

  # 'iter' errors
  length_error(iter, 1, "network.compare")
  typeof_error(iter, "numeric", "network.compare")

  # 'verbose' errors
  length_error(verbose, 1, "network.compare")
  typeof_error(verbose, "logical", "network.compare")

  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "network.compare")
    typeof_error(seed, "numeric", "network.compare")
    range_error(seed,  c(0, Inf), "network.compare")
  }

  # Check for usable data
  if(needs_usable(list(...))){
    base <- usable_data(base, verbose)
  }

  # Check for usable data
  if(needs_usable(list(...))){
    comparison <- usable_data(comparison, verbose)
  }

  # Error checks complete

  # Check for length of input in expected length
  if(dim(base)[2] != dim(comparison)[2]){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Number of columns for '", deparse(substitute(base)),
        "' (", dim(base)[2], ") did not match number of columns for '",
        deparse(substitute(comparison)), "' (", dim(comparison)[2], ").\n",
        "Number of columns in these datasets must match"
      ),
      call = "network.compare"
    )
  }

  # Return 'base' and 'comparison'
  return(list(base = base, comparison = comparison))

}
