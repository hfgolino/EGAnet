#' @title Plot Clustered Individual Networks
#'
#' @description Visualize clusters of individual networks identified by
#' \code{\link[EGAnet]{infoCluster}} applied to a \code{\link[EGAnet]{dynEGA}} object.
#'
#' Provides three visualization modes:
#' \itemize{
#'   \item \code{"heatmap"}: Heatmap of variables \code{x} individuals, colored by community assignment within each cluster
#'   \item \code{"average"}: Cluster-average networks, obtained by averaging adjacency matrices across individuals in the cluster
#'   \item \code{"cluster_population"}: Cluster-level "population" networks, obtained by stacking optimized derivatives across individuals and estimating a new EGA
#' }
#'
#' @param dynEGA.object A \code{dynEGA} object containing optimized individual results (see \code{\link[EGAnet]{dynEGA}}).
#'
#' @param clustering An \code{infoCluster} object containing cluster assignments for individuals
#' (see \code{\link[EGAnet]{infoCluster}}).
#'
#' @param type Character. Type of visualization to produce:
#' \itemize{
#'   \item \code{"cluster_population"} (default) -  Cluster-level population networks by stacking derivatives across individuals and re-estimating with \code{\link[EGAnet]{EGA}}
#'   \item \code{"average"} - Average adjacency networks with prototype community assignment (majority vote) per cluster
#'   \item \code{"heatmap"} - Heatmap of variable–community assignments across individuals, grouped by cluster
#' }
#'
#' @param palette Character. Palette name passed to \code{\link[ggplot2]{scale_fill_brewer}}
#' when \code{type = "heatmap"}. Default = \code{"Set3"}.
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
#'   \item \code{"cluster_population"} - A \code{patchwork} object combining cluster-level population network plots
#'   \item \code{"average"} - A \code{patchwork} object combining cluster-average network plots
#'   \item \code{"heatmap"} - A \code{ggplot2} object showing variables × individuals, colored by community and faceted by cluster
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
#' # Heatmap of variable–community assignments
#' plot_clusters(optimized_all, clust, type = "heatmap")
#'
#' # Average networks per cluster
#' plot_clusters(optimized_all, clust, type = "average")
#'
#' # Cluster-level population networks
#' plot_clusters(optimized_all, clust, type = "cluster_population")
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @seealso \code{\link[EGAnet]{dynEGA}}, \code{\link[EGAnet]{infoCluster}}, \code{\link[EGAnet]{EGA}}
#'
#' @export

plot.clusters <- function(dynEGA.object, clustering,
                          type = c("cluster_population", "average", "heatmap"),
                          palette = "Set3",
                          node.size = 3,
                          label.size = 1.2,
                          ...) {
  type <- match.arg(type)

  # ID --> cluster
  cluster_df <- data.frame(
    ID = names(clustering$clusters),
    cluster = clustering$clusters,
    stringsAsFactors = FALSE
  )

  if (type == "heatmap") {
    wc_df <- lapply(names(dynEGA.object$dynEGA$individual), function(id) {
      res <- dynEGA.object$dynEGA$individual[[id]]
      if (!is.null(res$wc)) {
        data.frame(
          ID = id,
          variable = names(res$wc),
          community = as.integer(res$wc),
          stringsAsFactors = FALSE
        )
      } else NULL
    }) %>% dplyr::bind_rows()

    wc_df <- wc_df %>%
      dplyr::left_join(cluster_df, by = "ID")

    return(
      ggplot2::ggplot(wc_df, aes(x = variable, y = ID, fill = factor(community))) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(~ cluster, scales = "free_y") +
        ggplot2::scale_fill_brewer(palette = palette) +
        ggplot2::theme_minimal() +
        ggplot2::labs(fill = "Community", x = "Variable", y = "Individual")
    )

  } else if (type == "average") {
    average_network <- function(ids) {
      nets <- lapply(ids, function(id) {
        dynEGA.object$dynEGA$individual[[id]]$network
      })
      Reduce("+", nets) / length(nets)
    }

    plots <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
      ids_in_cluster <- cluster_df$ID[cluster_df$cluster == cl]
      net_mean <- average_network(ids_in_cluster)

      wc_list <- lapply(ids_in_cluster, function(id) {
        dynEGA.object$dynEGA$individual[[id]]$wc
      })
      wc_mat <- do.call(cbind, lapply(wc_list, as.integer))
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

    return(patchwork::wrap_plots(plots))

  } else if (type == "cluster_population") {
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

    return(patchwork::wrap_plots(plots))
  }
}
