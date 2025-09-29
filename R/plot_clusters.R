#' @title Plot Clustered Individual Networks
#'
#' @description Visualize clusters of individual networks identified by
#' \code{\link[EGAnet]{infoCluster}} applied to a \code{\link[EGAnet]{dynEGA}} object.
#'
#' Provides three visualization modes:
#' \itemize{
#'   \item \code{"population"}: Cluster-level "population" networks, obtained by stacking optimized derivatives across individuals and estimating a new EGA
#'   \item \code{"average"}: Cluster-average networks, obtained by averaging adjacency matrices across individuals in the cluster
#' }
#'
#' @param dynEGA.object A \code{dynEGA} object containing optimized individual results (see \code{\link[EGAnet]{dynEGA}}).
#'
#' @param clustering An \code{infoCluster} object containing cluster assignments for individuals
#' (see \code{\link[EGAnet]{infoCluster}}).
#' @param include Numeric vector. Specifies which cluster memberships should be
#' explicitly included in the plot. By default, all clusters are shown. Use this
#' argument to restrict the visualization to a subset of clusters (e.g.,
#' `include = c(1, 3, 5)` will only display clusters 1, 3, and 5). This option
#' applies to all plot types (`"population"` and ``"average"`).
#' @param type Character. Type of visualization to produce:
#' \itemize{
#'   \item \code{"population"} (default) -  Cluster-level population networks by stacking derivatives across individuals and re-estimating with \code{\link[EGAnet]{EGA}}
#'   \item \code{"average"} - Average adjacency networks with prototype community assignment (majority vote) per cluster
#' }
#'
#' @param node.size Numeric. Node size in network visualizations (\code{type = "average"} or \code{"cluster_population"}). Default = \code{3}.
#'
#' @param label.size Numeric. Label size in network visualizations (\code{type = "average"} or \code{"cluster_population"}). Default = \code{1.2}.
#'
#' @param ... Additional arguments passed to \code{\link[EGAnet]{EGA}}
#' (used only when \code{type = "cluster_population"}).
#'
#' @return
#' Depending on the chosen \code{type}:
#' \itemize{
#'   \item \code{"population"} - A \code{patchwork} object combining cluster-level population network plots
#'   \item \code{"average"} - A \code{patchwork} object combining cluster-average network plots
#' }
#'
#' @details
#' This function provides flexible visualization of individual clusters obtained from
#' \code{\link[EGAnet]{infoCluster}}:
#' \itemize{
#'   \item Heatmaps are useful for inspecting how community memberships differ across clusters
#'   \item Average networks show the mean connectivity structure for each cluster, with prototype memberships by majority vote
#'   \item Cluster-level population networks treat each cluster as a "mini-population," stacking derivatives across individuals to re-estimate an EGA
#' }
#'
#' @examples
#' \dontrun{
#' # Run dynEGA with optimization
#' optimized_all <- dynEGA(data = sim.dynEGA,
#'                         level = c("individual", "population"),
#'                         optimization = TRUE)
#'
#' # Cluster individuals
#' clust <- infoCluster(dynEGA.object = optimized_all)
#'
#'
#' # Average networks per cluster
#' plot_clusters(optimized_all, clust, type = "average")
#'
#' # Cluster-level population networks
#' plot_clusters(optimized_all, clust, type = "population")
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @seealso \code{\link[EGAnet]{dynEGA}}, \code{\link[EGAnet]{infoCluster}}, \code{\link[EGAnet]{EGA}}
#'
#' @export
plot_clusters <- function(dynEGA.object, clustering, include,
                          type = c("population", "average"),
                          node.size = 3,
                          label.size = 1.2,
                          ...) {
  type <- match.arg(type)

  if(missing(include) == TRUE){
    include <- unique(clustering$clusters)
  }


  # Select a subset of clusters to include:
  clustering_subset <- clustering
  clustering_subset$clusters <- clustering$clusters[which(clustering$clusters==include)]

  # ID --> cluster
  cluster_df <- data.frame(
    ID = names(clustering_subset$clusters),
    cluster = clustering_subset$clusters,
    stringsAsFactors = FALSE
  )

 if (type == "average") {
    average_network <- function(ids) {
      nets <- lapply(ids, function(id) {
        dynEGA.object$dynEGA$individual[[id]]$network
      })
      Reduce("+", nets) / length(nets)
    }

    plots <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
      ids_in_cluster <- cluster_df$ID[cluster_df$cluster == cl]
      net_mean <- average_network(ids_in_cluster)

      # Get community memberships for all individuals in cluster
      wc_list <- lapply(ids_in_cluster, function(id) {
        dynEGA.object$dynEGA$individual[[id]]$wc
      })

      # Create membership matrix (individuals across columns, nodes down rows)
      wc_mat <- do.call(cbind, lapply(wc_list, as.integer))

      # Initialize consensus matrix
      n_nodes <- nrow(wc_mat)
      consensus_matrix <- matrix(
        nrow = n_nodes, ncol = n_nodes
      )

      # Loop over to get proportions
      for(i in seq_len(n_nodes)){
        for(j in i:n_nodes){
          # Fill consensus matrix
          consensus_matrix[i,j] <- consensus_matrix[j,i] <- sum(
            wc_mat[i,] == wc_mat[j,],
            na.rm = TRUE
          )
        }
      }

      # Divide by number of individuals
      consensus_matrix <- consensus_matrix / length(ids_in_cluster)

      # Derive consensus communities from consensus matrix
      wc_proto <- apply(wc_mat, 1, function(x) {
        as.integer(names(sort(table(x), decreasing = TRUE))[1])
      })

      EGAnet:::basic_plot_setup(
        network = net_mean,
        wc = wc_proto,
        title = paste("Cluster", cl),
        node.size = node.size, label.size = label.size, ...
      )
    })

    return(do.call(ggpubr::ggarrange, plots))

  } else if (type == "population") {
    plots <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
      ids_in_cluster <- cluster_df$ID[cluster_df$cluster == cl]

      # Stack derivatives
      pooled_data <- do.call(
        rbind,
        lapply(ids_in_cluster, function(id) {
          dynEGA.object$dynEGA$individual[[id]]$data_used
        })
      )

      # Run EGA on pooled derivatives
      pop_res <- EGA(pooled_data, plot.EGA = FALSE, verbose = FALSE, ...)

      # Plot
      EGAnet:::basic_plot_setup(
        network = pop_res$network,
        wc = pop_res$wc,
        title = paste("Cluster", cl),
        node.size = node.size, label.size = label.size, ...
      )
    })

    return(do.call(ggpubr::ggarrange, plots))
  }
}
