#' @title Compares Dynamic Network Structures Using Permutation
#'
#' @description A permutation implementation to determine statistical
#' significance of whether the dynamic network structures are different from one another
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis as well as
#' an ID column
#'
#' @param groups Numeric or character vector (length = \code{nrow(data)}).
#' Group membership corresponding to each case in data
#'
#' @param paired Boolean (length = 1).
#' Whether groups are repeated measures representing
#' paired samples.
#' Defaults to \code{FALSE}
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
#' \item \code{"cosine"} --- Uses \code{\link[EGAnet]{cosine}} to compute cosine similarity
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
#' @param id Numeric or character (length = 1).
#' Number or name of the column identifying each individual.
#' Defaults to \code{NULL}
#'
#' @param n.embed Numeric (length = 1).
#' Defaults to \code{5}.
#' Number of embedded dimensions (the number of observations to
#' be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"n.embed = 5"} will use five consecutive observations
#' to estimate a single derivative
#'
#' @param n.embed.optimize Boolean (length = 1).
#' If \code{TRUE}, performs optimization of \code{n.embed} for each individual,
#' then constructs the population based on optimized derivatives. When \code{TRUE},
#' individual networks are considered of interest and will always be output.
#' Defaults to \code{FALSE}
#'
#' @param tau Numeric (length = 1).
#' Defaults to \code{1}.
#' Number of observations to offset successive embeddings in
#' the \code{\link[EGAnet]{Embed}} function.
#' Generally recommended to leave "as is"
#'
#' @param delta Numeric (length = 1).
#' Defaults to \code{1}.
#' The time between successive observations in the time series (i.e, lag).
#' Generally recommended to leave "as is"
#'
#' @param use.derivatives Numeric (length = 1).
#' Defaults to \code{1}.
#' The order of the derivative to be used in the analysis.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{0} --- No derivatives; consistent with moving average
#'
#' \item \code{1} --- First-order derivatives; interpreted as "velocity" or
#' rate of change over time
#'
#' \item \code{2} --- Second-order derivatives; interpreted as "acceleration" or
#' rate of the rate of change over time
#'
#' }
#'
#' Generally recommended to leave "as is"
#'
#' @param na.derivative Character (length = 1).
#' How should missing data in the embeddings be handled?
#' Available options (see Boker et al. (2018) in \code{\link[EGAnet]{glla}} references for more details):
#'
#' \itemize{
#'
#' \item \code{"none"} (default) --- does nothing and leaves \code{NA}s in data
#'
#' \item \code{"kalman"} --- uses Kalman smoothing (\code{\link[stats]{KalmanSmooth}}) with
#' structural time series models (\code{\link[stats]{StructTS}}) to impute missing values.
#' This approach models the underlying temporal dependencies (trend, seasonality, autocorrelation)
#' to generate estimates for missing observations while preserving the original time scale.
#' More computationally intensive than the other methods but typically provides the
#' most accurate imputation by respecting the stochastic properties of the time series
#'
#' \item \code{"rowwise"} --- adjusts time interval with respect to each embedding ensuring
#' time intervals are adaptive to the missing data (tends to be more accurate than \code{"none"})
#'
#' \item \code{"skipover"} --- "skips over" missing data and treats the non-missing points
#' as continuous points in time (note that the time scale shifts to the "per mean time interval,"
#' which is different and \emph{larger} than the original scale)
#'
#' }
#'
#' @param zero.jitter Numeric (length = 1).
#' Small amount of Gaussian noise added to zero variance derivatives to prevent
#' estimation failures. For more than one variable, noise is generated
#' multivariate normal distribution to ensure orthogonal noise is added.
#' The jitter preserves the overall structure but avoids singular
#' covariance matrices during network estimation.
#' Defaults to \code{0.001}
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
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results.
#' Set for reproducible results.
#' See \href{https://r-ega.net/articles/reproducibility-prng.html}{Reproducibility and PRNG}
#' for more details on random number generation in \code{\link{EGAnet}}
#'
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments that can be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{EGA}}, and
#' \code{\link[EGAnet]{jsd}}
#'
#' @return Returns a list:
#'
#' \item{network}{Data frame with row names of each measure, empirical value (\code{statistic}), and
#' \emph{p}-value based on the permutation test (\code{p.value})}
#'
#' \item{edges}{List containing matrices of values for empirical values (\code{statistic}),
#' \emph{p}-values (\code{p.value}), and Benjamini-Hochberg corrected \emph{p}-values (\code{p.adjusted})}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Three similar groups
#'
#' # Set seed
#' set.seed(42)
#'
#' # Simulate dynamic data
#' participants <- lapply(
#'   seq_len(50), function(i){
#'
#'     # Get output
#'     output <- simDFM(
#'       variab = 6, timep = 15,
#'       nfact = 2, error = 0.100,
#'       dfm = "DAFS", loadings = 0.60,
#'       autoreg = 0.80, crossreg = 0.10,
#'       var.shock = 0.36, cov.shock = 0.18,
#'       burnin = 2000
#'     )
#'
#'     #  Add ID
#'     df <- data.frame(
#'       ID = i,
#'       Group = rep(1:3, each = 5),
#'       output$data
#'     )
#'
#'     # Return data
#'     return(df)
#'
#'   }
#' )
#'
#' # Put participants into a data frame
#' df <- do.call(rbind.data.frame, participants)
#'
#' \dontrun{
#' # Perform comparison
#' dynamic.network.compare(
#'   data = df[,-2], groups = df[,2], paired = TRUE,
#'   # EGA arguments
#'   corr = "auto", na.data = "pairwise", model = "glasso",
#'   # dynEGA arguments
#'   id = NULL, n.embed = 3,
#'   tau = 1, delta = 1, use.derivatives = 1,
#'   # Permutation arguments
#'   iter = 1000, ncores = 2, verbose = TRUE, seed = 42
#' )}
#'
#' # Two similar groups and one different
#'
#' # Simulate dynamic data
#' participants <- lapply(
#'   seq_len(50), function(i){
#'
#'     # Get output
#'     output <- simDFM(
#'       variab = 4, timep = 5,
#'       nfact = 3, error = 0.100,
#'       dfm = "DAFS", loadings = 0.60,
#'       autoreg = 0.80, crossreg = 0.10,
#'       var.shock = 0.36, cov.shock = 0.18,
#'       burnin = 2000
#'     )
#'
#'     #  Add ID
#'     df <- data.frame(
#'       ID = i,
#'       Group = rep(3, each = 5),
#'       output$data
#'     )
#'
#'     # Return data
#'     return(df)
#'
#'   }
#' )
#'
#' # Replace group 3
#' new_group <- do.call(rbind.data.frame, participants)
#' df[df$Group == 3,] <- new_group
#'
#' \dontrun{
#' # Perform comparison
#' dynamic.network.compare(
#'   data = df[,-2], groups = df[,2], paired = TRUE,
#'   # EGA arguments
#'   corr = "auto", na.data = "pairwise", model = "glasso",
#'   # dynEGA arguments
#'   id = NULL, n.embed = 3,
#'   tau = 1, delta = 1, use.derivatives = 1,
#'   # Permutation arguments
#'   iter = 1000, ncores = 2, verbose = TRUE, seed = 42
#' )}
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
# Updated 19.11.2025
dynamic.network.compare <- function(
    data, groups, paired = FALSE,
    # EGA arguments
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    # dynEGA arguments
    id = NULL, n.embed = 5, n.embed.optimize = FALSE,
    tau = 1, delta = 1, use.derivatives = 1,
    na.derivative = c("none", "kalman", "rowwise", "skipover"),
    zero.jitter = 0.001,
    # Permutation arguments
    iter = 1000, ncores, seed = NULL, verbose = TRUE,
    ...
)
{

  # Experimental warning
  experimental("dynamic.network.compare")

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", dynamic.network.compare)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}

  # Argument errors (return data in case of tibble)
  data <- dynEGA_errors(
    data, id, NULL, n.embed, tau, delta,
    use.derivatives, zero.jitter, n.embed.optimize,
    ncores, verbose
  )

  # Update 'n.embed.optimize'
  n.embed.optimize <- attributes(data)$n.embed.optimize

  # Check for input errors
  groups <- dynamic.network.compare_errors(data, groups, paired, iter, seed, ...)

  # Set ellipse
  ellipse <- list(...)

  # Check for seed
  if(!is.null(seed)){
    seeds <- reproducible_seeds(iter - 1, seed)
  }else{

    # Set all seeds to zero (or random)
    seeds <- rep(0, iter)

    # Check for external suppression (from `invariance`)
    if(!"suppress" %in% names(ellipse) || !ellipse$suppress){
      message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results")
    }

  }

  # Ensure data has variable names and get dimensions
  data <- ensure_dimension_names(data)
  dimensions <- dim(data)
  dimension_names <- dimnames(data)

  # Check for ID
  if(is.null(id)){

    # Check for column with name ID
    index <- "id" %in% tolower(dimension_names[[2]])

    # Check for existing column
    if(any(index)){
      id <- which(index)
    }else{

      # Send error
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Argument 'id' is `NULL` and no column in 'data' could be ",
          "identified with the label 'id'.\n\n",
          "Please input a variable name or index into 'id' or include a column in 'data' ",
          "named 'ID'"
        ),
        call = "dynamic.network.compare"
      )

    }

  }


  # Get unique groups (factored)
  unique_factors <- na.omit(unique(groups))

  # Set groups as factors
  groups <- as.numeric(factor(groups, levels = unique_factors))

  # Get unique groups (numeric)
  unique_groups <- na.omit(unique(groups))

  # Generate all possible group pairings
  group_pairs <- combn(unique_groups, 2, simplify = FALSE)
  pairs_length <- length(group_pairs)
  pairs_sequence <- seq_len(pairs_length)

  # Get dynamic EGA lists
  dynamic_ega <- lapply(
    unique_groups, function(group){

      # Return output
      return(
        dynEGA(
          data = data[groups == group,], corr = corr,
          na.data = na.data, model = model,
          id = id, n.embed = n.embed,
          n.embed.optimize = n.embed.optimize,
          tau = tau, delta = delta,
          use.derivatives = use.derivatives,
          na.derivative = na.derivative,
          zero.jitter = zero.jitter,
          level = "population",
          seed = seed,
          verbose = FALSE,
          ...
        )
      )

    }
  )

  # Add names
  names(dynamic_ega) <- unique_groups

  # Get statistics
  empirical <- lapply(group_pairs, function(pair){

    # Simplify grabbing networks
    network1 <- dynamic_ega[[pair[1]]]$dynEGA$population$network
    network2 <- dynamic_ega[[pair[2]]]$dynEGA$population$network

    # Return statistics
    return(
      list(
        empirical_values = c(
          "Frobenius" = frobenius(network1, network2),
          "JSS" = 1 - jsd(network1, network2, ...),
          "Total Strength" = sum(colSums(abs(network1), na.rm = TRUE), na.rm = TRUE) -
            sum(colSums(abs(network2), na.rm = TRUE), na.rm = TRUE)
        ),
        empirical_matrix = network1 - network2
      )
    )

  })

  # Extract derivatives
  derivatives <- lapply(
    unique_groups, function(group){

      # Extract index
      index <- dynamic_ega[[group]]$Derivatives$EstimatesDF

      # Column names
      column_names <- colnames(index)

      # Target data frame
      df <- index[
        , c(
          grep(paste0("Ord", use.derivatives), column_names),
          which(column_names == "id")
        )
      ]

      # Create new 'id', 'group', and 'id_group'
      df$id <- as.numeric(factor(df$id))
      df$group <- group
      df$id_group <- paste0(df$id, "_", df$group)

      # Return derivatives
      return(df)

    }
  )

  # Create individual sequences
  individual_sequence <- lapply(derivatives, function(x){seq_along(unique(x$id))})
  group_sequence <- rep(unique_groups, times = nvapply(individual_sequence, length))

  # Set up derivatives data frame
  derivatives <- do.call(rbind.data.frame, derivatives)

  # Remove row names (creates issues in `EGA`)
  row.names(derivatives) <- NULL

  # Get indices for variables
  variable_indices <- !colnames(derivatives) %in% c("id", "group", "id_group")

  # Loop some number of iterations
  permutation_list <- parallel_process(
    iterations = (iter - 1), datalist = seeds, FUN = function(seed, ...){

      # Create new group derivatives
      new_derivatives <- lapply(group_pairs, function(pair){

        # Obtain new group membership
        new_membership <- shuffle(
          group_sequence[group_sequence %in% pair], seed = seed
        )

        # Target sequence
        target_sequence <- individual_sequence[[pair[1]]]

        # If paired groups
        if(paired){

          # Re-factor the memberships for easier assignments
          numeric_membership <- as.numeric(factor(new_membership, levels = pair))

          # Get new groups for first group
          new_group <- numeric_membership[target_sequence]

          # Set up for the second group
          new_membership <- c(new_group, 3 - new_group)

          # Replace new membership with levels
          new_membership <- pair[new_membership]

        }

        # Create new ID memberships
        index <- paste0(unlist(individual_sequence[pair]), "_", new_membership)

        # Find IDs in derivatives
        return(list(
          derivatives[
            derivatives$id_group %in% index[target_sequence],
            variable_indices
          ],
          derivatives[
            derivatives$id_group %in% index[-target_sequence],
            variable_indices
          ]
        ))

      })

      # Get statistics
      permutated <- lapply(new_derivatives, function(pair){

        # Simplify grabbing networks
        network1 <- EGA(pair[[1]], plot.EGA = FALSE, ...)$network
        network2 <- EGA(pair[[2]], plot.EGA = FALSE, ...)$network

        # Return statistics
        return(
          list(
            empirical_values = c(
              "Frobenius" = frobenius(network1, network2),
              "JSS" = 1 - jsd(network1, network2, ...),
              "Total Strength" = sum(colSums(abs(network1), na.rm = TRUE), na.rm = TRUE) -
                sum(colSums(abs(network2), na.rm = TRUE), na.rm = TRUE)
            ),
            empirical_matrix = network1 - network2
          )
        )

      })

      # Return permutated results
      return(permutated)

    }, ncores = ncores, progress = verbose, ...
  )

  # Separate into pairwise comparisons
  ## Network
  permutated_values <- lapply(pairs_sequence, function(i){

    # Nest loop over iterations
    return(
      do.call(rbind, lapply(permutation_list, function(x){x[[i]]$empirical_values}))
    )

  })
  ## Edges
  permutated_matrices <- lapply(pairs_sequence, function(i){

    # Nest loop over iterations
    return(lapply(permutation_list, function(x){x[[i]]$empirical_matrix}))

  })

  # Compute results
  ## Network
  result_values <- lapply(pairs_sequence, function(i){

    # Return data frame
    return(
      t(
        data.frame(
          "statistic" = empirical[[i]]$empirical_values,
          "p.value" = c(
            mean( # Frobenius
              c(TRUE, permutated_values[[i]][,1] <= empirical[[i]]$empirical_values[1]),
              na.rm = TRUE
            ),
            mean( # JSS
              c(TRUE, permutated_values[[i]][,2] <= empirical[[i]]$empirical_values[2]),
              na.rm = TRUE
            ),
            mean( # Total strength
              c(TRUE, abs(permutated_values[[i]][,3]) >= abs(empirical[[i]]$empirical_values[3])),
              na.rm = TRUE
            )
          ),
          "M_permutated" = colMeans(permutated_values[[i]], na.rm = TRUE),
          "SD_permutated" = apply(permutated_values[[i]], 2, sd, na.rm = TRUE)
        )
      )
    )

  })
  ## Edges
  ## Get lower triangle
  lower_triangle <- lower.tri(dynamic_ega[[1]]$dynEGA$population$network)

  ## Get the p-values for edges
  result_edges <- lapply(
    pairs_sequence, function(i){

      # Get p-values
      edge_p_adjusted <- edge_p <- apply(
        simplify2array(
          lapply(permutated_matrices[[i]], function(x){
            abs(x) >= abs(empirical[[i]]$empirical_matrix)
          })
        ), 1:2, mean, na.rm = TRUE
      )

      # Compute adjusted p-values
      edge_p_adjusted_lower <- p.adjust(edge_p[lower_triangle], method = "BH")
      edge_p_adjusted[lower_triangle] <- edge_p_adjusted_lower
      edge_p_adjusted <- t(edge_p_adjusted)
      edge_p_adjusted[lower_triangle] <- edge_p_adjusted_lower

      # Return results
      return(
        list(
          statistic = empirical[[i]]$empirical_matrix,
          p.value = edge_p,
          p.adjusted = edge_p_adjusted,
          M_permutated = apply(simplify2array(permutated_matrices[[i]]), 1:2, mean, na.rm = TRUE),
          SD_permutated = apply(simplify2array(permutated_matrices[[i]]), 1:2, sd, na.rm = TRUE)
        )
      )

    }
  )

  # Set up results
  results <- lapply(pairs_sequence, function(i){

    # Create as a 'network.compare' class
    result <- list(network = result_values[[i]], edges = result_edges[[i]])

    # Set class
    class(result) <- "network.compare"

    # Return result
    return(result)

  })

  # Add names
  names(results) <- ulapply(group_pairs, function(x){
    paste0(unique_factors[x][1], " - ", unique_factors[x][2])
  })

  # Add dynamic EGA results to results
  results$dynEGA <- dynamic_ega

  # Set class
  class(results) <- "dynamic.network.compare"

  # Set up results
  return(results)

}

#' @noRd
# Errors ----
# Updated 03.06.2025
dynamic.network.compare_errors <- function(data, groups, paired, iter, seed, ...)
{

  # 'groups' errors
  object_error(groups, c("vector", "matrix", "data.frame", "factor"), "dynamic.network.compare")
  groups <- force_vector(groups)
  length_error(groups, dim(data)[1], "dynamic.network.compare")

  # 'paired' errors
  length_error(paired, 1, "dynamic.network.compare")
  typeof_error(paired, "logical", "dynamic.network.compare")

  # 'iter' errors
  length_error(iter, 1, "dynamic.network.compare")
  typeof_error(iter, "numeric", "dynamic.network.compare")

  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "dynamic.network.compare")
    typeof_error(seed, "numeric", "dynamic.network.compare")
    range_error(seed,  c(0, Inf), "dynamic.network.compare")
  }

  # Return groups
  return(groups)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 14.05.2025
print.dynamic.network.compare <- function(x, ...)
{

  # Print results (relies on `network.compare`)
  print(x[names(x) != "dynEGA"])

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 14.05.2025
summary.dynamic.network.compare <- function(object, ...)
{

  # Same as print
  print(object, ...)

}