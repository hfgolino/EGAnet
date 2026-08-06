#' @title Pairwise Resolution Iteration via Subgraph Modularity (PRISM)
#'
#' @description Iteratively tests pairwise unions of communities for a
#' modularity gain at the subgraph level and merges communities whenever
#' the union increases modularity relative to the communities' separate
#' memberships. Iteration continues until no pair of communities yields a
#' consistent, positive gain, at which point the partition is returned.
#' The approach operates on single subgraphs rather than the full network
#' at each comparison, which is intended to overcome the modularity
#' resolution limit described by Fortunato and Barthelemy (2007)
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
#' \code{\link[igraph]{cluster_louvain}}. Accepts \code{resolution} (defaults
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
#' Passed on to the function selected via \code{algorithm}.
#' When \code{FALSE}, singleton communities will be set to
#' missing (\code{NA}); otherwise, when \code{TRUE}, singleton
#' communities will be allowed
#'
#' @param ... Additional arguments to be passed on to the community
#' detection function selected via \code{algorithm}
#' (see \code{\link{community.consensus}} or \code{\link{community.detection}}
#' for their respective arguments)
#'
#' @return Returns a named vector of final community memberships obtained
#' by iterating pairwise resolution over subgraph modularity until
#' convergence (i.e., until no further consistent modularity gain is found)
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
#' # Compute PRISM memberships (Louvain-based consensus)
#' community.prism(network)
#'
#' # Compute PRISM memberships (single-algorithm community detection)
#' community.prism(network, algorithm = "walktrap")
#'
#' @references
#' Granell, C., Gomez, S., & Arenas, A. (2012).
#' Hierarchical multiresolution method to overcome the resolution limit in complex networks.
#' \emph{International Journal of Bifurcation and Chaos}, \emph{22}(07), 1250171.
#'
#' Fortunato, S., & Barthelemy, M. (2007).
#' Resolution limit in community detection.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{104}(1), 36-41.
#'
#' @export
#'
# Pairwise Resolution Iteration via Subgraph Modularity ----
# Updated 05.08.2026
community.prism <- function(
    network, algorithm = c(
      "edge_betweenness", "fast_greedy",
      "fluid", "infomap", "label_prop",
      "leading_eigen", "leiden", "louvain", "optimal",
      "spinglass", "walktrap"
    ),
    allow.singleton = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  algorithm <- set_default(algorithm, "louvain", community.detection)

  # Check for {igraph} network
  if(is(network, "igraph")){
    network <- igraph2matrix(network)
  }

  # Use network matrix
  network <- abs(as.matrix(network))

  # Check for names
  network <- ensure_dimension_names(network)

  # Set algorithm
  algorithm_FUN <- switch(
    algorithm,
    "louvain" = community.consensus,
    community.detection
  )

  # Estimate initial memberships
  wc <- algorithm_FUN(network, allow.singleton = allow.singleton, ...)

  # Set node names
  node_names <- names(wc)

  # Set while loop
  while(TRUE){

    # Set updated and previous memberships
    updated_wc <- previous_wc <- wc

    # Number of clusters
    n_clusters <- unique_length(wc)
    seq_clusters <- seq_len(n_clusters)

    # Initialize gain
    gain <- -Inf

    # Store results
    gain_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters)

    # Store results
    store <- lapply(seq_clusters, function(x){
      vector("list", n_clusters)
    })

    # Set indices
    indices <- lapply(seq_clusters, function(i){
      which(wc == seq_clusters[i])
    })

    # Loop over target
    for(i in seq_clusters){

      # Loop over comparator
      for(j in seq_clusters){

        # Set index of pairwise clusters
        index <- swiftelse(i != j, c(indices[[i]], indices[[j]]), indices[[i]])

        # Set subgraph and membership
        subgraph <- network[index, index]
        subgraph_wc <- algorithm_FUN(subgraph, allow.singleton = allow.singleton, ...)

        # Compute modularity gain
        gain <- modularity(subgraph, subgraph_wc) - modularity(subgraph, wc[index])

        # Check for any gain
        if(gain < 1e-08){
          break
        }

        # Update gain matrix
        gain_matrix[i,j] <- gain

        # Store memberships
        store[[i]][[j]] <- reindex_memberships(subgraph_wc[node_names[indices[[i]]]])

      }

    }

    # Check for consistent gains
    expand_index <- which(rowSums(gain_matrix > 0) == n_clusters)

    # Check for no expansion
    if(length(expand_index) == 0){
      break
    }

    # Loop over to expand
    for(index in expand_index){

      # Obtain clusters
      expansion <- count_table(do.call(rbind, store[[index]]))
      ordered <- expansion[order(expansion$Value, decreasing = TRUE),]

      # Maximum counts
      max_counts <- ordered$Value == max(ordered$Value)

      # Check for ties
      if(sum(max_counts) > 1){

        # Obtain solutions
        solutions <- ordered[max_counts, -dim(ordered)[2]]

        # Obtain modularity from solutions
        wc_index <- wc == index
        Qs <- apply(solutions, 1, modularity, network = network[wc_index, wc_index])

        # Use highest modularity gain
        new_wc <- solutions[which.max(Qs),]

      }else{

        # Use most common
        new_wc <- ordered[1, -dim(expansion)[2]]

      }

      # Update memberships
      updated_wc[names(new_wc)] <- as.numeric(n_clusters + new_wc)
      updated_wc <- reindex_memberships(updated_wc)
      n_clusters <- unique_length(updated_wc)

    }

    # Check that there was an update
    if(all(updated_wc == previous_wc)){
      break
    }

    # Update memberships
    wc <- updated_wc

  }

  # Return memberships
  return(wc)

}