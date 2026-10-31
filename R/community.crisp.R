#' @title Cohesion and Resistance In Subgraph Partitioning (CRISP)
#'
#' @description Iteratively tests whether re-estimating communities within
#' each community's subgraph yields a stronger structure. Each community is
#' scored by the geometric mean of its \emph{cohesion} and its
#' \emph{separation} (the \code{strong_structure} score). A split is accepted
#' when the size-weighted mean score of the new communities exceeds the score
#' of the community being split. Iteration continues until no community can be
#' improved by splitting, at which point the partition is returned. A network
#' that cannot be improved by splitting is returned as a single community.
#'
#' \strong{Cohesion} is the algebraic connectivity of a (sub)graph -- the
#' second-smallest eigenvalue of its normalized Laplacian
#' (Fiedler, 1973; Chung, 1997) -- rescaled to range from 0 (disconnected)
#' to 1 (a complete graph with equal weights).
#' \strong{Separation} is one minus the largest \emph{merge propensity} with
#' any other community: the observed weight between two communities relative
#' to the weight expected under a configuration-model null within their union
#' (capped at 1). Singletons attached only to a sibling therefore have
#' separation 0 and are not rewarded for being split off. Edge weights are
#' treated as absolute values
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
#' entire CRISP procedure reproducible
#'
#' @param ... Additional arguments to be passed on to the community
#' detection function selected via \code{algorithm}
#' (see \code{\link{community.consensus}} or \code{\link{community.detection}}
#' for their respective arguments)
#'
#' @return Returns a named vector of final community memberships obtained
#' by iteratively re-estimating communities within subgraphs until no split
#' increases the size-weighted mean structure score (geometric mean of
#' cohesion and separation). Nodes without any connections
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
#' # Compute CRISP memberships (Louvain-based consensus)
#' community.crisp(network)
#'
#' # Compute CRISP memberships (single-algorithm community detection)
#' community.crisp(network, algorithm = "walktrap")
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
# Cohesion and Resistance In Subgraph Partitioning ----
# Updated 31.10.2026
community.crisp <- function(
    network, algorithm = c(
      "edge_betweenness", "fast_greedy",
      "fluid", "infomap", "label_prop",
      "leading_eigen", "leiden", "louvain", "optimal",
      "spinglass", "walktrap"
    ), allow.singleton = FALSE, seed = NULL, ...
)
{

  # Obtain ellipse
  ellipse <- list(...)

  # Check for missing arguments (argument, default, function)
  algorithm <- set_default(algorithm, "louvain", community.detection)

  # Argument errors (return network in case of {igraph} or tibble)
  network <- community.crisp_errors(network, allow.singleton, seed)

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

    # Initialize memberships
    wc <- structure(
      rep(1, dim(network)[2]),
      names = dimnames(network)[[2]]
    )

    # Set while loop
    while(TRUE){

      # Set previous memberships (snapshot for this pass)
      previous_wc <- wc

      # Set number of clusters
      n_clusters <- unique_length(previous_wc)

      # Compute structure score (cohesion and separation) of each community
      scores <- strong_structure(network, wc)$overall

      # Loop over communities
      for(i in seq_len(n_clusters)){

        # Obtain index
        index <- previous_wc == i

        # Single nodes cannot be split any further
        if(sum(index) < 2){
          next
        }

        # Set target network
        target_network <- network[index, index, drop = FALSE]

        # Determine connected nodes within subgraph
        connected <- colSums(target_network) != 0

        # Initialize new memberships
        new_wc <- rep(NA_real_, length(connected))

        # Estimate communities among connected nodes
        # (isolated nodes are not passed to the algorithm)
        if(any(connected)){
          new_wc[connected] <- reindex_memberships(
            algorithm_FUN(
              target_network[connected, connected, drop = FALSE],
              allow.singleton = TRUE, seed = seed, ...
            )
          )
        }

        # Nodes isolated within the subgraph become singletons
        isolated <- is.na(new_wc)
        new_wc[isolated] <- max(new_wc, 0, na.rm = TRUE) + seq_len(sum(isolated))

        # Reindex new memberships
        new_wc <- reindex_memberships(new_wc)

        # Set number of new communities
        n_new <- unique_length(new_wc)

        # Compute structure score of each new community
        new_scores <- strong_structure(target_network, new_wc)$overall

        # Determine whether to update (size-weighted mean structure score)
        if(weighted.mean(new_scores, tabulate(new_wc, n_new)) > scores[i]){

          # Update memberships (offset avoids collisions with existing labels)
          wc[index] <- new_wc + max(wc)

        }

      }

      # Reindex memberships
      wc <- reindex_memberships(wc)

      # Check that there was an update
      if(all(previous_wc == wc)){
        break
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
# algorithm = "louvain"
# allow.singleton = FALSE; seed = NULL
# ellipse = list()

#' @noRd
# Errors ----
# Updated 31.10.2026
community.crisp_errors <- function(network, allow.singleton, seed)
{

  # 'network' errors
  if(is(network, "igraph")){ # Check for {igraph} network
    network <- igraph2matrix(network)
  }else{
    object_error(network, c("matrix", "data.frame", "tibble"), "community.crisp")
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

  # 'allow.singleton' errors
  length_error(allow.singleton, 1, "community.crisp")
  typeof_error(allow.singleton, "logical", "community.crisp")

  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "community.crisp")
    typeof_error(seed, "numeric", "community.crisp")
    range_error(seed, c(0, Inf), "community.crisp")
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

#' @noRd
# Strong structure function ----
# Updated 31.10.2026
strong_structure <- function(network, membership)
{

  # Absolute weights, no self-loops
  network <- abs(network)
  diag(network) <- 0

  # Communities
  communities <- sort(unique(membership))
  K <- length(communities)

  # Initialize
  community_cohesion <- rep(NA_real_, K)
  merge <- matrix(0, nrow = K, ncol = K)

  # Loop over communities
  for(i in seq_len(K)){

    # Target community
    target <- membership == communities[i]

    # Cohesion (needs enough nodes for the spectral cut to be meaningful)
    community_cohesion[i] <- cohesion(network[target, target, drop = FALSE])

    # Merge propensity with every other community (local null in the union)
    for(j in seq_len(K)[-i]){

      other <- membership == communities[j]
      union <- target | other

      # Strengths within the union
      k_A <- colSums(network[union, target, drop = FALSE])
      k_B <- colSums(network[union, other, drop = FALSE])

      # Observed and expected cross-weight
      W_AB <- sum(network[target, other, drop = FALSE])
      E_AB <- sum(k_A) * sum(k_B) / (sum(k_A) + sum(k_B))

      # Observed / expected (capped at 1)
      # E_AB is NaN (0 / 0) when both communities have no weight inside their union
      merge[i, j] <- if(!is.finite(E_AB) || E_AB <= 0) 1 else min(W_AB / E_AB, 1)

    }

  }

  # Separation: resistance to the most mergeable partner
  separation <- if(K > 1) 1 - apply(merge, 1, max) else 1

  # Overall: geometric mean of cohesion and separation (per community)
  overall <- sqrt(community_cohesion * separation)

  # Return results
  return(list(cohesion = community_cohesion, separation = separation, overall = overall))

}