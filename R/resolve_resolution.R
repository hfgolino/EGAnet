#' @noRd
# Updated 26.07.2026
resolve_resolution <- function(network, algorithm, ...)
{

  # Based on single subgraph:
  # Granell, C., Gomez, S., & Arenas, A. (2012).
  # Hierarchical multiresolution method to overcome the resolution limit in complex networks.
  # International Journal of Bifurcation and Chaos, 22(07), 1250171.

  # Mentioned single subgraph:
  # Fortunato, S., & Barthelemy, M. (2007).
  # Resolution limit in community detection.
  # Proceedings of the National Academy of Sciences, 104(1), 36–41.

  # Set algorithm
  algorithm_FUN <- switch(
    algorithm,
    "louvain" = community.consensus,
    community.detection
  )

  # Estimate initial memberships
  wc <- algorithm_FUN(network, algorithm, ...)

  # Set while loop
  while(TRUE){

    # Set previous memberships
    previous_wc <- wc

    # Get number of nodes
    nodes <- length(wc)

    # Check for names on memberships
    if(is.null(names(wc))){
      names(wc) <- paste0("V", format_integer(nodes, digits(nodes) - 1))
    }

    # Set node names
    node_names <- names(wc)

    # Number of clusters
    n_clusters <- unique_length(wc)
    seq_clusters <- seq_len(n_clusters)

    # Initialize gain
    gain <- -Inf

    # Store results
    gain_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters)

    # Store results
    store <- lapply(seq_len(n_clusters), function(x){
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
        index <- c(indices[[i]], indices[[j]])

        # Set subgraph and membership
        subgraph <- network[index, index]
        subgraph_wc <- algorithm_FUN(subgraph, algorithm, allow.singleton = TRUE, ...)

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
      new_wc <- expansion[order(expansion$Value, decreasing = TRUE),][1, -dim(expansion)[2]]

      # Update memberships
      wc[names(new_wc)] <- as.numeric(n_clusters + new_wc)
      wc <- reindex_memberships(wc)
      n_clusters <- unique_length(wc)

    }

    # Check that there was an update
    if(all(wc == previous_wc)){
      break
    }

  }

  # Return memberships
  return(wc)

}
