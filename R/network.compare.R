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
#' \code{\link[EGAnet]{community.consensus}},
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
#' \dontrun{# Perform comparison
#' results <- network.compare(group1, group2)
#'
#' # Print results
#' print(results)
#'
#' # Plot edge differences
#' plot(results)}
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
# Updated 13.05.2025
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
  empirical_values <- abs(
    c(
      "Frobenius" = frobenius(base_empirical_network, comparison_empirical_network),
      "JSS" = 1 - jsd(base_empirical_network, comparison_empirical_network, ...),
      "Total Strength" = sum(colSums(abs(base_empirical_network), na.rm = TRUE), na.rm = TRUE) -
        sum(colSums(abs(comparison_empirical_network), na.rm = TRUE), na.rm = TRUE)
    )
  )

  # Empirical differences
  empirical_matrix <- base_empirical_network - comparison_empirical_network

  # Create combined dataset
  combined <- rbind(base, comparison)

  # Set up indices
  combined_index <- nrow_sequence(combined)
  base_length <- dim(base)[1]

  # Perform permutations
  permutated <- parallel_process(
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
        list(
          empirical_values = c(
            "Frobenius" = frobenius(base_network, comparison_network),
            "JSS" = 1 - jsd(base_network, comparison_network, ...),
            "Total Strength" = sum(colSums(abs(base_network), na.rm = TRUE), na.rm = TRUE) -
              sum(colSums(abs(comparison_network), na.rm = TRUE), na.rm = TRUE)
          ),
          empirical_matrix = base_network - comparison_network
        )
      )

    }, ncores = ncores, progress = verbose, ...
  )

  # Separate values from matrices
  permutated_values <- do.call(
    rbind, lapply(permutated, function(x){x$empirical_values})
  )
  permutated_matrices <- lapply(permutated, function(x){x$empirical_matrix})

  # Get the p-values for edges
  edge_p <- apply(
    simplify2array(
      lapply(permutated_matrices, function(x){
        abs(x) >= abs(empirical_matrix)
      })
    ), 1:2, mean, na.rm = TRUE
  )

  # Get lower triangle
  lower_triangle <- lower.tri(edge_p)
  edge_p_adjusted <- edge_p
  edge_p_adjusted_lower <- p.adjust(edge_p[lower_triangle], method = "BH")
  edge_p_adjusted[lower_triangle] <- edge_p_adjusted_lower
  edge_p_adjusted <- t(edge_p_adjusted)
  edge_p_adjusted[lower_triangle] <- edge_p_adjusted_lower

  # Set up results
  results <- list(
    network = t(data.frame(
      "statistic" = empirical_values,
      "p.value" = c(
        mean(permutated_values[,1] <= empirical_values[1]),
        mean(permutated_values[,2] <= empirical_values[2]),
        mean(abs(permutated_values[,3]) >= abs(empirical_values[3]))
      ),
      "M_permutated" = colMeans(permutated_values),
      "SD_permutated" = apply(permutated_values, 2, sd)
    )),
    edges = list(
      statistic = empirical_matrix,
      p.value = edge_p,
      p.adjusted = edge_p_adjusted,
      M_permutated = apply(simplify2array(permutated_matrices), 1:2, mean, na.rm = TRUE),
      SD_permutated = apply(simplify2array(permutated_matrices), 1:2, sd, na.rm = TRUE)
    )
  )

  # Set class
  class(results) <- "network.compare"

  # Return statistics
  return(results)

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

#' @exportS3Method
# S3 Print Method ----
# Updated 02.09.2024
print.network.compare <- function(x, ...)
{

  # Print network results
  print(x$network, digits = 4)

  # Print edges
  cat(
    paste0(
      "\nNumber of significant edges (p <= 0.05): ",
      length(which(x$edges$p.value <= 0.05)) / 2,
      "\nNumber of significant edges (p_BH <= 0.10): ",
      length(which(x$edges$p.adjusted <= 0.10)) / 2
    )
  )

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 02.09.2024
summary.network.compare <- function(object, ...)
{

  # Same as print
  print(object, ...)

}

#' @exportS3Method
# S3 Plot Method ----
# Updated 02.09.2024
plot.network.compare <- function(x, p_type = c("p", "p_BH"), p_value = 0.05, ...)
{

  # Get p errors
  p_type <- swiftelse(missing(p_type), "p", match.arg(p_type))
  range_error(p_value, c(0, 1), "plot.network.compare")

  # Get number of nodes
  nodes <- dim(x$edges$statistic)[2]
  node_sequence <- seq_len(nodes)

  # Set up data frame
  plot_df <- data.frame(
    Rows = rep(node_sequence, each = nodes),
    Columns = rep(node_sequence, times = nodes),
    statistic = swiftelse(
      x$edges$statistic == 0, "",
      format_decimal(as.numeric(x$edges$statistic), 2)
    ),
    p.value = swiftelse(
      as.numeric(
        x$edges[[swiftelse(p_type == "p", "p.value", "p.adjusted")]]
      ) <= p_value, 1, 0
    )
  )

  # Set names (if possible)
  node_names <- dimnames(x$edges$statistic)[[2]]
  if(!is.null(node_names)){

    # Replace in plot data frame
    plot_df$Rows <- factor(node_names[plot_df$Rows], levels = node_names)
    plot_df$Columns <- factor(node_names[plot_df$Columns], levels = rev(node_names))

  }

  # Plot
  ggplot2::ggplot(
    data = plot_df,
    ggplot2::aes(x = Rows, y = Columns, fill = p.value, label = statistic)
  ) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::geom_text() +
    ggplot2::scale_fill_gradient(low = "white", high = "grey", limits = c(0, 1)) +
    ggplot2::labs(
      title = "Significant Edge Differences",
      subtitle = swiftelse(
        p_type == "p",
        bquote(paste(italic(p), " < ", .(p_value))),
        bquote(paste(italic(p)[adj.], " < ", .(p_value)))
      )
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none"
    )

}

#' @noRd
# Global variables needed for CRAN checks ----
# Updated 09.02.2024
utils::globalVariables(c("p.value", "statistic"))
