#' @title Hierarchical \code{\link[EGAnet]{EGA}}
#'
#' @description Estimates EGA using the lower-order solution of the Louvain
#' algorithm (\code{\link[igraph]{cluster_louvain}})to identify the lower-order
#' dimensions and then uses factor or network loadings to estimate factor
#' or network scores, which are used to estimate the higher-order dimensions
#' (for more details, see Jiménez et al., 2023)
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#' (does not accept correlation matrices)
#'
#' @param loading.method Character (length = 1).
#' Sets network loading calculation based on implementation
#' described in \code{"BRM"} (Christensen & Golino, 2021) or
#' an \code{"experimental"} implementation.
#' Defaults to \code{"BRM"}
#'
#' @param rotation Character.
#' A rotation to use to obtain a simpler structure.
#' For a list of rotations, see \code{\link[GPArotation]{rotations}} for options.
#' Defaults to \code{NULL} or no rotation.
#' By setting a rotation, \code{scores} estimation will be
#' based on the rotated loadings rather than unrotated loadings
#'
#' @param scores Character (length = 1).
#' How should scores for the higher-order structure be estimated?
#' Defaults to \code{"network"} for network scores computed using
#' the \code{\link[EGAnet]{net.scores}} function.
#' Set to \code{"factor"} for factor scores computed using
#' \code{\link[psych]{fa}}. Factors scores will be based on
#' \strong{EFA} (as in Jiménez et al., 2023)
#'
#' \emph{Factor scores use the number of communities from
#' \code{\link[EGAnet]{EGA}}. Estimated factor loadings may
#' not align with these communities. The plots using factor scores
#' will have higher order factors that may not completely map on to
#' the lower order communities. Look at
#' \code{$hierarchical$higher_order$lower_loadings} to determine the
#' composition of the lower order factors.}
#'
#' @param loading.structure Character (length = 1).
#' Whether simple structure or the saturated loading matrix
#' should be used when computing scores (\code{scores = "network"} only).
#' Defaults to \code{"simple"}
#'
#' \code{"simple"} structure more closely mirrors traditional
#' hierarchical factor analytic methods such as CFA; \code{"full"}
#' structure more closely mirrors EFA methods
#'
#' Simple structure is the more conservative (established) approach
#' and is therefore the default. Treat \code{"full"} as experimental
#' as proper vetting and validation has not been established
#'
#' @param impute Character (length = 1).
#' If there are any missing data, then imputation can be implemented.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"none"} --- Default. No imputation is performed
#'
#' \item \code{"mean"} --- The mean value of each variable is used to replace missing data
#' for that variable
#'
#' \item \code{"median"} --- The median value of each variable is used to replace missing data
#' for that variable
#'
#' }
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
#' @param lower.algorithm Character or
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to the lower order \code{"louvain"} with most common
#' consensus clustering (1000 iterations; see
#' \code{\link[EGAnet]{community.consensus}} for more details)
#'
#' Louvain with consensus clustering is \emph{strongly}
#' recommended. Using any other algorithm is considered
#' \emph{experimental} as they have not been designed to
#' capture lower order communities
#'
#' @param higher.algorithm Character or
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#'
#' \itemize{
#'
#' \item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details
#'
#' \item \code{"louvain"} --- By default, \code{"louvain"} will implement the Louvain algorithm using
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}}
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000}
#' unless specified otherwise
#'
#' \item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details
#'
#' }
#'
#' Using \code{algorithm} will set only \code{higher.algorithm} and
#' \code{lower.algorithm} will default to Louvain with most common
#' consensus clustering (1000 iterations)
#'
#' @param uni.method Character (length = 1).
#' What unidimensionality method should be used?
#' Defaults to \code{"louvain"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"expand"} --- Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation
#'
#' \item \code{"LE"} --- Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023)
#' \emph{Behavior Research Methods} simulation
#'
#' \item \code{"louvain"} --- Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Louvain solution is used; otherwise, regular EGA is used.
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}
#'
#' }
#'
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{EGA}}, and
#' \code{\link[GPArotation]{rotations}}
#'
#' @return Returns a list of lists containing:
#'
#' \item{lower_order}{\code{\link[EGAnet]{EGA}} results for the lower order structure}
#'
#' \item{higher_order}{\code{\link[EGAnet]{EGA}} results for the higher order structure}
#'
#' \item{parameters}{A list containing \code{lower_loadings} and \code{lower_scores}
#' that were used to estimate scores and the higher order \code{\link[EGAnet]{EGA}}
#' results, respectively}
#'
#' \item{dim.variables}{A data frame with variable names and their lower and higher
#' order assignments}
#'
#' \item{TEFI}{Generalized TEFI using \code{\link[EGAnet]{tefi}}}
#'
#' \item{plot.hierEGA}{Plot output if \code{plot.EGA = TRUE}}
#'
#' @references
#' \strong{Hierarchical EGA simulation} \cr
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A. P., & Garrido, L. E. (2023).
#' Dimensionality assessment in bifactor structures with multiple general factors: A network psychometrics approach.
#' \emph{Psychological Methods}.
#'
#' \strong{Conceptual implementation} \cr
#' Golino, H., Thiyagarajan, J. A., Sadana, R., Teles, M., Christensen, A. P., & Boker, S. M. (2020).
#' Investigating the broad domains of intrinsic capacity, functional ability and
#' environment: An exploratory graph analysis approach for improving analytical
#' methodologies for measuring healthy aging.
#' \emph{PsyArXiv}.
#'
#' @author
#' Marcos Jiménez <marcosjnezhquez@gmailcom>,
#' Francisco J. Abad <fjose.abad@uam.es>,
#' Eduardo Garcia-Garzon <egarcia@ucjc.edu>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>, and
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu.do>
#'
#' @examples
#' # Example using network scores
#' opt.hier <- hierEGA(
#'   data = optimism, scores = "network",
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' \donttest{
#' # Plot multilevel plot
#' plot(opt.hier, plot.type = "multilevel")
#'
#' # Plot multilevel plot with higher order
#' # border color matching the corresponding
#' # lower order color
#' plot(opt.hier, color.match = TRUE)
#'
#' # Plot levels separately
#' plot(opt.hier, plot.type = "separate")}
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
#'
# Hierarchical EGA ----
# Updated 24.10.2023
hierEGA <- function(
    data,
    # `net.scores` arguments
    loading.method = c("BRM", "experimental"),
    rotation = NULL, scores = c("factor", "network"),
    loading.structure = c("simple", "full"),
    impute = c("mean", "median", "none"),
    # `EGA` arguments
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    lower.algorithm = "louvain",
    higher.algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    plot.EGA = TRUE, verbose = FALSE,
    ...
)
{

  # Send experimental message (for now)
  experimental("hierEGA")

  # Argument errors (return data in case of tibble)
  data <- hierEGA_errors(data, plot.EGA, verbose, ...)

  # Get ellipse arguments
  ellipse <- list(needs_usable = FALSE, ...)

  # Check for missing arguments (argument, default, function)
  ## `net.scores`
  loading.method <- set_default(loading.method, "brm", net.loads)
  scores <- set_default(scores, "network", hierEGA)
  loading.structure <- set_default(loading.structure, "simple", hierEGA)
  impute <- set_default(impute, "none", net.scores)
  ## `EGA`
  corr <- set_default(corr, "auto", hierEGA)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)

  ## Determine if single algorithm is entered
  ## Overwrite 'lower.algorithm' and 'higher.algorithm'
  if("algorithm" %in% names(ellipse)){

    # Set lower and higher order algorithm
    lower.algorithm <- "louvain"
    ellipse$consensus.method <- "most_common"
    ellipse$consensus.iter <- 1000
    higher.algorithm <- ellipse$algorithm

    # Remove 'algorithm' from ellipse to avoid conflicts
    ellipse <- ellipse[names(ellipse) != "algorithm"]

  }

  ## Handle 'lower.algorithm' and 'higher.algorithm' arguments (workaround)
  algorithm <- lower.algorithm
  lower.algorithm <- set_default(algorithm, "louvain", community.detection)
  algorithm <- higher.algorithm
  higher.algorithm <- set_default(algorithm, "walktrap", community.detection)

  # Get EGA
  lower_order_result <- do.call(
    what = EGA,
    args = c(
      # Standard `EGA` arguments
      list(
        data = data, corr = corr, na.data = na.data,
        model = model, algorithm = lower.algorithm, uni.method = uni.method,
        plot.EGA = FALSE, verbose = verbose, order = "lower"
      ),
      # Send ellipse arguments
      ellipse
    )
  )

  # Ensure data has names (needed in `net.scores`)
  if(is.null(dimnames(data)[[2]])){
    dimnames(data)[[2]] <- dimnames(lower_order_result$network)[[2]]
  }

  # Check for no dimensions in lower order
  if(lower_order_result$n.dim == 0){

    # Return empty results except lower order
    return(
      list(
        lower_order = lower_order_result,
        higher_order = NULL,
        parameters = list(
          lower_loadings = NULL,
          lower_scores = NULL
        )
      )
    )

  }else if(scores == "factor"){

    # Send warning
    warning(
      paste0(
        "Factor dimensions may not align with `EGA` detected communities. ",
        "Interpret results with caution. ",
        "\n\nUse `$higher_order$lower_loadings` to guide interpretations"
      ),
      call. = FALSE
    )

    # Get arguments for EFA
    efa_ARGS <- obtain_arguments(efa_scores, list(...))

    # Set data, correlation matrix, and number of factors
    efa_ARGS[
      c("data", "correlation_matrix", "nfactors", "impute")
    ] <- list(
      data, lower_order_result$correlation,
      lower_order_result$n.dim, impute
    )

    # Estimate EFA
    efa_output <- do.call(efa_scores, efa_ARGS)

    # Get scores
    score_estimates <- efa_output$scores

    # Get loading names
    loading_names <- dimnames(efa_output$loadings)

    # Obtain highest loadings for each variable
    lower_wc <- max.col(abs(efa_output$loadings), "first")

    # Get max `wc`
    max_wc <- max(lower_wc, na.rm = TRUE)

    # Format `wc` to be characters
    lower_wc <- format_integer(lower_wc, digits(max_wc) - 1)
    names(lower_wc) <- loading_names[[1]]

    # Set up dimension names to mirror network output
    dimnames(efa_output$loadings)[[2]] <- format_integer(
      as.numeric(gsub("MR", "", loading_names[[2]])),
      digits(max_wc) - 1
    )
    dimnames(score_estimates)[[2]] <- dimnames(efa_output$loadings)[[2]]

    # Put loadings in descending order
    lower_loadings <- descending_order(
      standardized = efa_output$loadings, wc = lower_wc,
      unique_communities = sort(dimnames(efa_output$loadings)[[2]]),
      node_names = loading_names[[1]]
    )

  }else if(scores == "network"){

    # Compute network scores
    network_output <- net.scores(
      data = data, A = lower_order_result,
      rotation = rotation, loading.method = loading.method,
      scoring.method = "network",
      loading.structure = loading.structure,
      impute = impute,
      ...
    )

    # Score estimates
    if(is.null(rotation)){
      score_estimates <- network_output$scores$std.scores
      lower_loadings <- network_output$loadings$std
    }else{
      score_estimates <- network_output$scores$rot.scores
      lower_loadings <- network_output$loadings$rotated
    }

  }

  # Check for unidimensional lower order result
  # (rare but can happen)
  if(lower_order_result$n.dim == 1){

    # Set up results
    results <- list(
      lower_order = lower_order_result,
      higher_order = lower_order_result,
      parameters = list(
        lower_loadings = lower_loadings,
        lower_scores = score_estimates
      )
    )

    # Set higher order flag
    higher_order <- FALSE

  }else{

    # Set up results
    results <- list(
      lower_order = lower_order_result,
      higher_order = do.call(
        what = EGA,
        args = c(
          # Standard `EGA` arguments
          list(
            data = score_estimates, corr = corr, na.data = na.data,
            model = model, algorithm = higher.algorithm, uni.method = uni.method,
            plot.EGA = FALSE, verbose = FALSE
            # issues at the this level may be inconsistent
            # with the output (e.g., singleton communities)
            # for this reason, 'verbose = FALSE'
          ),
          # Send ellipse arguments
          ellipse
        )
      ),
      parameters = list(
        lower_loadings = lower_loadings,
        lower_scores = score_estimates
      )
    )

    # For higher order, allow singleton communities
    if(anyNA(results$higher_order$wc)){

      # Get missing indices
      missing_index <- is.na(results$higher_order$wc)

      # Update missing indices with singleton values
      results$higher_order$wc[missing_index] <-
        seq_len(sum(missing_index)) + results$higher_order$n.dim
    }

    # Set higher order flag
    higher_order <- TRUE

  }

  # Add dimension variables like `EGA`
  results$dim.variables <- fast.data.frame(
    c(
      names(results$lower_order$wc),
      as.vector(results$lower_order$wc),
      as.vector(
        single_revalue_memberships( # function in `bootEGA` internals
          results$lower_order$wc, results$higher_order$wc
        )
      )
    ),
    nrow = length(results$lower_order$wc), ncol = 3,
    colnames = c("items", "lower", "higher")
  )

  # Set "methods" attributes
  attr(results, "methods") <- list(
    loading.method = loading.method,
    rotation = rotation,
    scores = scores
  )

  # Set class
  class(results) <- "hierEGA"

  # Obtain generalized TEFI
  gTEFI <- tefi(results)

  # Set up TEFI results
  results$lower_order$TEFI <- gTEFI$Lower.Order.VN
  results$TEFI <- gTEFI$VN.Entropy.Fit

  # Set up check for higher order results
  if(higher_order){

    # Set higher order TEFI
    results$higher_order$TEFI <- gTEFI$Higher.Order.VN

    # Message for correlated factor vs. bifactor

    # Set up messages
    ## General start
    general_start <- paste0(
      "Based on lower (", round(results$lower_order$TEFI, 3),
      ") and higher (", round(results$higher_order$TEFI, 3),
      ") TEFI, there is better fit for a "
    )

    ## Alternates
    first_order <- paste0("correlated lower order factor structure than bifactor structure")
    bifactor <- paste0("bifactor structure than correlated lower order factor structure")

    # Compare
    attr(results, "methods")$interpretation <- swiftelse(
      results$lower_order$TEFI < results$higher_order$TEFI,
      paste0(general_start, first_order), paste0(general_start, bifactor)
    )

  }

  # Check for plot
  if(higher_order && plot.EGA){

    # Get plot
    results$plot.hierEGA <- plot(results, ...)

    # Actually send plot
    plot(results$plot.hierEGA)

  }

  # Return results
  return(results)

}

# Bug checking ----
# data = NetworkToolbox::neoOpen; loading.method = "BRM"
# loading.structure = "simple"
# rotation = "geominQ"; scores = "network"
# scoring.method = "network"; impute = "none"
# corr = "auto"; na.data = "pairwise"; model = "glasso"
# ellipse = list(algorithm = "walktrap"); uni.method = "louvain"
# plot.EGA = FALSE; verbose = FALSE

#' @noRd
# Argument errors ----
# Updated 07.09.2023
hierEGA_errors <- function(data, plot.EGA, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "hierEGA")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'plot.EGA' errors
  length_error(plot.EGA, 1, "hierEGA")
  typeof_error(plot.EGA, "logical", "hierEGA")

  # 'verbose' errors
  length_error(verbose, 1, "hierEGA")
  typeof_error(verbose, "logical", "hierEGA")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return usable data in case of tibble
  return(data)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 17.11.2023
print.hierEGA <- function(x, ...)
{

  # Print lower order
  cat(
    styletext(
      text = styletext(
        text =  "Lower Order\n\n",
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )

  # Print lower network
  print(x$lower_order)

  # Add breakspace
  cat("\n\n------------\n\n")

  # Print higher order
  cat(
    styletext(
      text = styletext(
        text =  "Higher Order\n\n",
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )

  # Print higher network
  print(x$higher_order)

  # Add break space
  cat("\n\n----\n\n")

  # Print TEFI
  cat(paste0("Generalized TEFI: ", round(x$TEFI, 3)))

  # Check for interpretation
  if("interpretation" %in% names(attr(x, "methods"))){

    # Add break space
    cat("\n\n----\n\n")

    # Print interpretation
    cat(paste0("Interpretation: ", attr(x, "methods")$interpretation))

  }

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 17.07.2023
summary.hierEGA <- function(object, ...)
{
  print(object, ...) # same as print
}


#' @exportS3Method
# S3 Plot Method ----
# Updated 01.03.2024
plot.hierEGA <- function(
    x, plot.type = c("multilevel", "separate"),
    color.match = FALSE, ...
)
{

  # Get ellipse arguments for defaults
  ellipse <- list(...)

  # Get attributes
  methods_attributes <- attr(x, "methods")

  # Check for factor scores
  plot.type <- swiftelse(
    methods_attributes$scores == "factor",
    "separate", # must be separate
    set_default(plot.type, "multilevel", plot.hierEGA)
  )

  # Multilevel plot
  if(plot.type == "multilevel"){

    # Set edge size
    if(!"edge.size" %in% ellipse){
      ellipse$edge.size <- 8 # default in `basic_plot_setup`
    }

    # Get names for levels
    lower_names <- names(x$lower_order$wc)
    higher_names <- names(x$higher_order$wc)
    all_names <- c(lower_names, higher_names)

    # Number of total nodes
    total_nodes <- length(all_names)

    # Initialize network
    hierarchical_network <- matrix(
      0, nrow = total_nodes, ncol = total_nodes,
      dimnames = list(all_names, all_names)
    )

    # Population hybrid network with lower and higher order networks
    hierarchical_network[lower_names, lower_names] <- x$lower_order$network
    hierarchical_network[higher_names, higher_names] <- x$higher_order$network

    # Add assignment loading to hierarchical network
    for(community in dimnames(x$parameters$lower_loadings)[[2]]){

      # Get assignment loadings
      assignment_loadings <- x$parameters$lower_loadings[
        names(x$lower_order$wc)[x$lower_order$wc == as.numeric(community)],
        community
      ]

      # Add to hierarchical network
      hierarchical_network[names(assignment_loadings), community] <-
        hierarchical_network[community, names(assignment_loadings)] <-
        assignment_loadings

    }

    # Get digits based on lower order communities (includes `NA`)
    places <- digits(length(unique(x$lower_order$wc))) - 1

    # Set up names for higher order memberships
    higher_order_names <- paste0(
      "Higher_", format_integer(x$higher_order$wc, places)
    )

    # Set up plot list as standard `EGA`
    plot_list <- list(
      # Round hierarchical network to 4
      # (same as what's used in `basic_plot_setup`)
      network = round(hierarchical_network, 4),
      wc = c(
        paste0("Lower_", format_integer(x$lower_order$wc, places)),
        higher_order_names
      )
    )

    # Set class
    class(plot_list) <- "EGA"

    # Create the initial plot
    initial_plot <- plot(plot_list, ..., arguments = TRUE)

    # Make a copy of the hierarchical network
    hierarchical_copy <- plot_list$network

    # Add assignment loading to hierarchical network
    for(community in dimnames(x$parameters$lower_loadings)[[2]]){

      # Get assignment loadings
      assignment_loadings <- x$parameters$lower_loadings[
        names(x$lower_order$wc)[x$lower_order$wc == as.numeric(community)],
        community
      ]

      # Add to hierarchical network
      hierarchical_copy[names(assignment_loadings), community] <-
        hierarchical_copy[community, names(assignment_loadings)] <-
        1 # Set all values to 1

    }

    # Update plot list
    plot_list$network <- hierarchical_copy

    # Create the second plot
    second_plot <- plot(plot_list, arguments = TRUE)

    # Update multilevel edge appearances

    # Get edge indices
    edge_index <- second_plot$ARGS$edge.size == ellipse$edge.size

    # Set edge color
    edge_color <- initial_plot$ARGS$edge.color
    edge_color[edge_index] <- "grey"

    ## Edge line type
    line_type <- initial_plot$ARGS$edge.lty
    line_type[edge_index] <- "dashed"

    ## Edge size
    edge_size <- initial_plot$ARGS$edge.size
    edge_size[edge_index] <- edge_size[edge_index] * 0.50

    ## Edge alpha
    edge_alpha <- initial_plot$ARGS$edge.alpha
    edge_alpha[edge_index] <- 0.50

    # Get mode (layout)
    mode <- scale(initial_plot$ARGS$mode, scale = FALSE)

    # Get maximum y-position
    y_position <- max(mode[,2]) + 15 # for some space

    # Get higher order indices
    higher_index <- seq_along(higher_names)

    # Shift higher order nodes up
    mode[higher_index, 2] <- mode[higher_index, 2] + y_position +
      abs(min(mode[higher_index, 2]))

    # Re-scale x-position
    for(index in higher_index){

      # Get lower mode position
      lower_mode_position <- match(
        lower_names[x$lower_order$wc == index],
        initial_plot$ARGS$node.label
      )

      # Get higher mode position
      higher_mode_position <- match(
        higher_names[index],
        initial_plot$ARGS$node.label
      )

      # Update x-position for higher order node
      mode[higher_mode_position, 1] <-
        median(mode[lower_mode_position, 1])

    }

    # Check for community matching between layers
    if(color.match){

      # Get first plot list
      first_plot_list <- silent_plot(
        plot_list,
        mode = mode,
        edge.size = edge_size,
        edge.color = edge_color,
        edge.alpha = edge_alpha,
        edge.lty = line_type,
        node.size = -1,
        ...,
        arguments = TRUE
      )

      # Separate plot and arguments
      first_layer <- first_plot_list$network_plot
      plot_ARGS <- first_plot_list$ARGS

      # Create copy of node colors for border colors
      border_color <- plot_ARGS$node.color

      # Get node names
      node_names <- dimnames(plot_ARGS$net)[[2]]

      # Get targets
      higher_order_targets <- match(names(x$higher_order$wc), node_names)

      # Target higher order colors
      border_color[higher_order_targets] <- unique(
        border_color[match(names(x$lower_order$wc), node_names)]
      )

      # Set stroke size to be thicker for higher order only
      ## Get stroke size
      stroke_size <- rep(
        swiftelse(
          is.null(first_layer$guides$colour$override.aes$stroke),
          1.5, first_layer$guides$colour$override.aes$stroke
        ), total_nodes
      )

      ## Adjust higher order only
      stroke_size[higher_order_targets] <- 2.5

      # Custom nodes: transparent insides and dark borders
      second_layer <- first_layer +
        ggplot2::geom_point( # transparent insides
          size = second_plot$ARGS$node.size + 0.50, shape = 19,
          color = plot_ARGS$node.color,
          alpha = plot_ARGS$node.alpha,
          show.legend = FALSE
        ) +
        ggplot2::geom_point( # dark borders
          size = second_plot$ARGS$node.size,
          color = border_color,
          shape = 1, stroke = stroke_size, alpha = 1
        ) +
        ggplot2::geom_text( # put text back on top
          ggplot2::aes(label = first_layer$data$label),
          color = "black",
          size = plot_ARGS$label.size
        ) +
        ggplot2::guides( # create legend with these settings
          color = ggplot2::guide_legend(
            override.aes = list(
              shape = 21,
              fill = unique(plot_ARGS$node.color),
              size = median(second_plot$ARGS$node.size, na.rm = TRUE),
              alpha = median(plot_ARGS$node.alpha, na.rm = TRUE),
              stroke = 1.5
            ),
            title = swiftelse(
              "legend.title" %in% names(ellipse),
              ellipse$legend.title, ""
            )
          )
        )

      # Send it
      return(silent_plot(second_layer))

    }else{
      return(
        silent_plot(
          plot_list,
          mode = mode,
          edge.size = edge_size,
          edge.color = edge_color,
          edge.alpha = edge_alpha,
          edge.lty = line_type,
          ...
        )
      )
    }

  }else if(plot.type == "separate"){ # Separate plot

    # Set labels
    if(!"labels" %in% ellipse){
      ellipse$labels <- c("Lower", "Higher")
    }

    # Plot lower and higher order side-by-side
    return(
      ggpubr::ggarrange(
        silent_plot(x$lower_order, ...),
        silent_plot(x$higher_order, ...),
        labels = ellipse$labels,
        ...
      )
    )

  }

}

#' @noRd
# Oblimin fit ----
# Updated 15.07.2023
oblimin <- function(L, gamma = 0)
{

  # Get dimensions of loadings
  dimensions <- dim(L)

  # Check for special case of gamma = 0
  if(gamma == 0){ # Set IgC matrix as identity
    IgC <- diag(dimensions[1])
  }else{

    # Get gC matrix
    gC <- matrix(
      gamma / dimensions[1],
      nrow = dimensions[1], ncol = dimensions[1]
    )

    # Subtract gC from identity matrix
    IgC <- diag(dimensions[1]) - gC

  }

  # Initialize N matrix
  N <- matrix(1, nrow = dimensions[2], ncol = dimensions[2])

  # Set diagonal to zero
  diag(N) <- 0

  # Square loadings matrix
  L2 <- L^2 # nanoseconds faster than L * L

  # Return f
  return(sum(diag(crossprod(L2, IgC %*% L2 %*% N))) * 0.25)

}

#' @noRd
# Oblimin rotation ----
# Updated 15.07.2023
oblimin_rotate <- function(
    loadings, nfactors,
    n.rotations = 10, maxit = 10000,
    eps = 1e-5, gamma = 0,
    rotate = "oblimin"
)
{

  # Check for unidimensional structure (28.11.2022)
  if(nfactors == 1){
    return(list(loadings = loadings, Phi = 1))
  }

  # Get factor sequence
  factor_sequence <- seq_len(nfactors)

  # Initialize minimum fit value
  fit_value <- Inf

  # Perform loop
  for(i in seq_len(n.rotations)){

    # Random values
    X <- nvapply(factor_sequence, function(x){rnorm(nfactors)}, LENGTH = nfactors)

    # Loading results
    result <- GPArotation::GPFoblq(
      A = loadings, method = rotate, Tmat = qr.Q(qr(X)),
      maxit = maxit, eps = eps
    )

    # Check for improvement in fit
    if(oblimin(result$loadings, gamma = gamma) < fit_value){
      return_result <- result[c("loadings", "Phi")]
    }

  }

  # Return results
  return(return_result)

}

#' @noRd
# Exploratory Factor Analysis Scores ----
# Updated 15.07.2023
efa_scores <- function(
    data, correlation_matrix, nfactors,
    fm = "minres", rotate = "oblimin",
    n.rotations = 10, maxit = 1e4,
    factor.scores = "Thurstone", gamma = 0,
    eps = 1e-5, impute
)
{

  # Start with {psych}'s factor analysis
  fa <- psych::fa(
    correlation_matrix,
    nfactors = nfactors,
    fm = fm, rotate = "none"
  )

  # Perform oblimin rotation
  rotation <- oblimin_rotate(
    loadings = fa$loadings, nfactors = nfactors,
    n.rotations = n.rotations, maxit = maxit,
    eps = eps, gamma = gamma, rotate = rotate
  )

  # Return result
  return(
    list(
      loadings = rotation$loadings,
      scores = psych::factor.scores(
        x = data, f = rotation$loadings,
        Phi = rotation$Phi, method = factor.scores,
        rho = correlation_matrix, missing = TRUE,
        impute = impute
      )$scores
    )
  )

  # At this point, only scores are used and reported
  # return(
  #   list(
  #     loadings = rotation$loadings, Phi = rotation$Phi,
  #     Shat = fa$model, uniquenesses = fa$uniquenesses,
  #     factor.scores = psych::factor.scores(
  #       x = data, f = rotation$loadings,
  #       Phi = rotation$Phi, method = factor.scores,
  #       rho = correlation_matrix, missing = TRUE,
  #       impute = impute
  #     )
  #   )
  # )

}



