#' @title Algorithmic Complexity for Emergent Subgraphs (ACES)
#'
#' @description Iteratively tests whether each community is internally
#' cohesive and, whenever a community's cohesion falls below
#' \code{cohesion.cut}, re-estimates communities within that community's
#' subgraph. Iteration continues until every community is cohesive (or can
#' no longer be split), at which point the partition is returned. A
#' network that is already cohesive as a whole is returned as a single
#' community. Cohesion is the algebraic connectivity of a (sub)graph --
#' the second-smallest eigenvalue of its normalized Laplacian
#' (Fiedler, 1973; Chung, 1997) -- rescaled to range from 0 (disconnected)
#' to 1 (a complete graph with equal weights)
#'
#' @param network Matrix or \code{igraph} network object
#'
#' @param algorithm Character or \code{igraph} \code{cluster_*} function
#' (length = 1).
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"edge_betweenness"} --- See \code{\link[igraph]{cluster_edge_betweenness}} for more details
#'
#' \item \code{"fast_greedy"} --- See \code{\link[igraph]{cluster_fast_greedy}} for more details
#'
#' \item \code{"fluid"} --- See \code{\link[igraph]{cluster_fluid_communities}} for more details
#'
#' \item \code{"infomap"} --- See \code{\link[igraph]{cluster_infomap}} for more details
#'
#' \item \code{"label_prop"} --- See \code{\link[igraph]{cluster_label_prop}} for more details
#'
#' \item \code{"leading_eigen"} --- See \code{\link[igraph]{cluster_leading_eigen}} for more details
#'
#' \item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details.
#' \emph{Note}: The Leiden algorithm will default to the
#' modularity objective function (\code{objective_function = "modularity"}).
#' Set \code{objective_function = "CPM"} to use the
#' Constant Potts Model instead (see examples)
#'
#' \item \code{"louvain"} --- Applies \code{EGAnet}'s own C implementation of the
#' Louvain algorithm (Blondel et al., 2008), which replaces
#' \code{\link[igraph]{cluster_louvain}}, via \code{\link{community.consensus}}.
#' Accepts \code{resolution} (defaults
#' to \code{1}), \code{seed} (defaults to \code{NULL}, i.e., not reproducible),
#' and \code{order} (\code{"higher"} or \code{"lower"}; defaults to
#' \code{"higher"}). Unlike \code{igraph::cluster_louvain}, which draws on R's
#' global random number generator, this implementation's randomness depends
#' only on \code{seed} -- making it reproducible regardless of R's RNG state,
#' call order, or parallelization
#'
#' \item \code{"optimal"} --- See \code{\link[igraph]{cluster_optimal}} for more details
#'
#' \item \code{"spinglass"} --- See \code{\link[igraph]{cluster_spinglass}} for more details
#'
#' \item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details
#'
#' }
#'
#' @param cohesion.cut Numeric (length = 1).
#' Cohesion threshold between \code{0} and \code{1}.
#' Communities with cohesion below this value are re-estimated
#' within their own subgraph; communities at or above this value
#' are retained as is. Higher values split communities more
#' aggressively (and yield more communities). Defaults to \code{0.50}
#'
#' @param allow.singleton Boolean (length = 1).
#' Whether singleton or single node communities should be allowed.
#' When \code{FALSE}, singleton communities will be set to
#' missing (\code{NA}); otherwise, when \code{TRUE}, singleton
#' communities will be allowed.
#' Singletons are always allowed \emph{during} the iterative procedure
#' and are only set to missing at the end. Defaults to \code{FALSE}
#'
#' @param seed Numeric (length = 1).
#' Sets seed for reproducible results.
#' Defaults to \code{NULL} or random results.
#' Passed on to the function selected via \code{algorithm} (relevant for
#' \code{algorithm = "louvain"}) for \emph{every} community estimation
#' performed during the iterative procedure, making the
#' entire ACES procedure reproducible
#'
#' @param ... Additional arguments to be passed on to the community
#' detection function selected via \code{algorithm}
#' (see \code{\link{community.consensus}} or \code{\link{community.detection}}
#' for their respective arguments)
#'
#' @return Returns a named vector of final community memberships obtained
#' by iteratively re-estimating communities within subgraphs whose cohesion
#' is below \code{cohesion.cut} until convergence (i.e., until all communities
#' are cohesive or can no longer be split). Nodes without any connections
#' are set to missing (\code{NA})
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Estimate network
#' network <- EBICglasso.qgraph(data = wmt)
#'
#' # Compute ACES memberships (Louvain-based consensus)
#' community.aces(network)
#'
#' # Compute ACES memberships (single-algorithm community detection)
#' community.aces(network, algorithm = "walktrap")
#'
#' @references
#' \strong{Algebraic connectivity} \cr
#' Fiedler, M. (1973).
#' Algebraic connectivity of graphs.
#' \emph{Czechoslovak Mathematical Journal}, \emph{23}(2), 298-305.
#'
#' \strong{Normalized Laplacian} \cr
#' Chung, F. R. K. (1997).
#' Spectral graph theory.
#' \emph{American Mathematical Society}.
#'
#' \strong{Louvain algorithm} \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#'
#' @export
#'
# Algorithmic Complexity for Emergent Subgraphs ----
# Updated 30.10.2026
community.aces <- function(
    network, algorithm = c(
      "edge_betweenness", "fast_greedy",
      "fluid", "infomap", "label_prop",
      "leading_eigen", "leiden", "louvain", "optimal",
      "spinglass", "walktrap"
    ), cohesion.cut = 0.50, allow.singleton = FALSE, seed = NULL, ...
)
{

  # Obtain ellipse
  ellipse <- list(...)

  # Check for missing arguments (argument, default, function)
  algorithm <- set_default(algorithm, "louvain", community.detection)

  # Argument errors (return network in case of {igraph} or tibble)
  network <- community.aces_errors(network, cohesion.cut, allow.singleton, seed)

  # Check for names
  network <- ensure_dimension_names(network)

  # Set original node names
  original_node_names <- dimnames(network)[[2]]

  # Remove self-loops (cohesion assumes no self-loops)
  diag(network) <- 0

  # Check for disconnected nodes
  include <- colSums(network) != 0

  # Create final memberships
  final_wc <- structure(
    rep(NA, length(include)),
    names = original_node_names
  )

  # Check for any connected nodes
  if(any(include)){

    # Update network
    network <- network[include, include, drop = FALSE]

    # Set algorithm
    algorithm_FUN <- switch(
      algorithm,
      "louvain" = community.consensus,
      community.detection
    )

    # Check for whole-network cohesion
    if(cohesion(network) >= cohesion.cut){

      # Return single community
      wc <- structure(
        rep(1, dim(network)[2]),
        names = dimnames(network)[[2]]
      )

    }else{ # Estimate communities

      # Initial communities
      wc <- reindex_memberships(
        algorithm_FUN(network, allow.singleton = TRUE, seed = seed, ...)
      )

      # Set while loop
      while(TRUE){

        # Set previous memberships (for this pass)
        previous_wc <- wc

        # Set max membership
        max_wc <- unique_length(wc)

        # Loop over communities
        for(i in seq_len(max_wc)){

          # Set index
          index <- which(previous_wc == i)

          # Single nodes cannot be split
          if(length(index) < 2){
            next
          }

          # Set target network
          target_network <- network[index, index, drop = FALSE]

          # Check cohesion
          if(cohesion(target_network) < cohesion.cut){

            # Initialize target memberships
            target_wc <- rep(NA_real_, length(index))

            # Determine connected nodes within subgraph
            connected <- colSums(target_network) != 0

            # Estimate communities among connected nodes
            if(sum(connected) > 1){
              target_wc[connected] <- algorithm_FUN(
                target_network[connected, connected, drop = FALSE],
                allow.singleton = TRUE, seed = seed, ...
              )
            }

            # Nodes not connected within the subgraph become singletons
            missing <- is.na(target_wc)
            target_wc[missing] <- max(target_wc, 0, na.rm = TRUE) + seq_len(sum(missing))

            # Update communities
            wc[index] <- target_wc + max_wc

            # Reindex communities
            wc <- reindex_memberships(wc)

            # Update max membership
            max_wc <- unique_length(wc)

          }

        }

        # Check that there was an update
        if(all(wc == previous_wc)){
          break
        }

      }

    }

    # Check singleton behavior
    if(!allow.singleton){

      # Determine whether there are any singleton communities
      membership_frequency <- fast_table(wc)

      # Singletons
      singletons <- membership_frequency == 1

      # Check for frequencies equal to one
      if(any(singletons)){

        # Identify communities
        singleton_communities <- as.numeric(
          names(membership_frequency)[singletons]
        )

        # Set values to NA
        wc[wc %in% singleton_communities] <- NA

      }

      # Reindex memberships
      wc <- reindex_memberships(wc)

    }

    # Set memberships
    final_wc[include] <- wc

  }

  # Add methods to membership attributes
  attr(final_wc, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm),
    objective_function = ellipse$objective_function
    # `objective_function` will be NULL unless it's there!
  )

  # Make memberships have S3 class
  class(final_wc) <- "EGA.community"

  # Return memberships
  return(final_wc)

}

# Bug checking ----
## Basic input
# network = EBICglasso.qgraph(wmt2[,7:24])
# algorithm = "louvain"; cohesion.cut = 0.85
# allow.singleton = FALSE; seed = NULL
# ellipse = list()

#' @noRd
# Errors ----
# Updated 30.10.2026
community.aces_errors <- function(network, cohesion.cut, allow.singleton, seed)
{

  # 'network' errors
  if(is(network, "igraph")){ # Check for {igraph} network
    network <- igraph2matrix(network)
  }else{
    object_error(network, c("matrix", "data.frame", "tibble"), "community.aces")
  }

  # Check for tibble
  if(get_object_type(network) == "tibble"){
    network <- as.data.frame(network)
  }

  # Use network matrix
  network <- abs(as.matrix(network))

  # Check for square network
  if(dim(network)[1] != dim(network)[2]){
    stop(
      "'network' must be a square matrix. Input has dimensions: ",
      paste0(dim(network), collapse = " x "), call. = FALSE
    )
  }

  # Check for missing values
  if(anyNA(network)){
    stop("'network' contains missing values (NA). Missing values are not allowed in the network.", call. = FALSE)
  }

  # Check for symmetric network
  if(!is_symmetric(network)){
    stop("'network' must be symmetric.", call. = FALSE)
  }

  # 'cohesion.cut' errors
  length_error(cohesion.cut, 1, "community.aces")
  typeof_error(cohesion.cut, "numeric", "community.aces")
  range_error(cohesion.cut, c(0, 1), "community.aces")

  # 'allow.singleton' errors
  length_error(allow.singleton, 1, "community.aces")
  typeof_error(allow.singleton, "logical", "community.aces")

  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "community.aces")
    typeof_error(seed, "numeric", "community.aces")
    range_error(seed, c(0, Inf), "community.aces")
  }

  # Return network in case of {igraph} or tibble
  return(network)

}

#' @noRd
# Cohesion function ----
# Updated 30.10.2026
cohesion <- function(network)
{

  # Number of nodes
  nodes <- dim(network)[2]

  # A single node cannot be split any further
  if(nodes < 2){
    return(1)
  }

  # Node strength
  strength <- colSums(network)

  # A zero-strength node is disconnected from the community (Laplacian is undefined, 0 / 0)
  if(any(strength == 0)){
    return(0)
  }

  # Normalized Laplacian
  laplacian <- diag(nodes) - network / sqrt(outer(strength, strength))

  # Eigenvalues are returned in decreasing order
  values <- eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values

  # Second smallest eigenvalue never exceeds nodes / (nodes - 1) for any weighted graph
  # (equality only for an equal-weight complete graph), so rescaling by that bound puts
  # cohesion on [0, 1] with a perfect community equal to 1 regardless of the number of nodes
  # `min` and `max` only guard against floating-point error; the bound is exact
  return(min(1, max(values[nodes - 1] * (nodes - 1) / nodes, 0)))

}
