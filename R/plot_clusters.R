#' @title Plot Clustered Individual Networks
#'
#' @description Visualize clusters of individual networks identified by
#' \code{\link[EGAnet]{infoCluster}} applied to a \code{\link[EGAnet]{dynEGA}} object.
#'
#' Provides visualization modes:
#' \itemize{
#'   \item \code{"population"}: Cluster-level "population" networks, obtained by stacking derivatives across individuals and estimating a new EGA
#'   \item \code{"average"}: Cluster-average networks, obtained by averaging adjacency matrices across individuals in the cluster
#' }
#'
#' @param dynEGA.object A \code{dynEGA} object containing individual results (see \code{\link[EGAnet]{dynEGA}}).
#' Can be from either optimized or non-optimized runs.
#'
#' @param clustering An \code{infoCluster} object containing cluster assignments for individuals
#' (see \code{\link[EGAnet]{infoCluster}}).
#'
#' @param include Numeric vector. Specifies which cluster memberships should be
#' explicitly included in the plot. By default, all clusters are shown. Use this
#' argument to restrict the visualization to a subset of clusters (e.g.,
#' `include = c(1, 3, 5)` will only display clusters 1, 3, and 5). This option
#' applies to all plot types (`"population"` and `"average"`).
#'
#' @param type Character. Type of visualization to produce:
#' \itemize{
#'   \item \code{"population"} (default) - Cluster-level population networks by stacking derivatives across individuals and re-estimating with \code{\link[EGAnet]{EGA}}
#'   \item \code{"average"} - Average adjacency networks with prototype community assignment (majority vote) per cluster
#' }
#'
#' @param node.size Numeric. Node size in network visualizations. Default = \code{3}.
#'
#' @param label.size Numeric. Label size in network visualizations. Default = \code{1.2}.
#'
#' @param ... Additional arguments passed to \code{\link[EGAnet]{EGA}}
#' (used only when \code{type = "population"}).
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
#'   \item Average networks show the mean connectivity structure for each cluster, with prototype memberships by majority vote
#'   \item Cluster-level population networks treat each cluster as a "mini-population," stacking derivatives across individuals to re-estimate an EGA
#' }
#'
#' When working with non-optimized dynEGA results, the function will extract the appropriate
#' derivative data from the \code{Derivatives$Estimates} component. When working with optimized
#' results, it will use the \code{data_used} component from each individual's results.
#'
#' @examples
#' \dontrun{
#' # Run dynEGA with optimization
#' optimized_all <- dynEGA(data = sim.dynEGA,
#'                         level = c("individual", "population"),
#'                         optimization = TRUE,
#'                         n.embed.all.ind = 3:15)
#'
#' # Cluster individuals
#' clust <- infoCluster(dynEGA.object = optimized_all)
#'
#' # Average networks per cluster
#' plot_clusters(optimized_all, clust, type = "average")
#'
#' # Cluster-level population networks
#' plot_clusters(optimized_all, clust, type = "population")
#'
#' # Cluster-level population networks, including only Cluster 2:
#' plot_clusters(optimized_all, clust, include = c(2), type = "population")
#'
#' # Also works with non-optimized dynEGA
#' standard_all <- dynEGA(data = sim.dynEGA,
#'                        level = c("individual", "population"))
#'
#' clust_standard <- infoCluster(dynEGA.object = standard_all)
#' plot_clusters(standard_all, clust_standard, type = "population")
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

  # Match argument
  type <- match.arg(type)

  # Set include to all clusters if not specified
  if(missing(include)){
    include <- unique(clustering$clusters)
  }

  # Filter clusters based on include
  selected_clusters <- clustering$clusters[clustering$clusters %in% include]

  # If no clusters match, stop
  if(length(selected_clusters) == 0){
    stop("No clusters found matching the 'include' specification.")
  }

  # Create cluster dataframe
  cluster_df <- data.frame(
    ID = names(selected_clusters),
    cluster = selected_clusters,
    stringsAsFactors = FALSE
  )

  # Check if we have optimization results (data_used exists)
  has_optimization <- !is.null(dynEGA.object$dynEGA$individual[[1]]$data_used)

  # Helper function to get derivative data for an individual
  get_individual_derivatives <- function(id){
    if(has_optimization){
      # Use optimized derivatives if available
      return(dynEGA.object$dynEGA$individual[[id]]$data_used)
    } else {
      # Extract from Derivatives$Estimates for non-optimized results
      # Find the derivatives for this ID
      derivatives_list <- dynEGA.object$Derivatives$Estimates

      # Check if we can match by names
      if(!is.null(names(derivatives_list)) && id %in% names(derivatives_list)){
        deriv_data <- derivatives_list[[id]]
      } else {
        # Try to match by ID attribute
        id_match <- sapply(derivatives_list, function(x) {
          id_attr <- attr(x, "ID")
          if(is.null(id_attr)) return(FALSE)
          return(as.character(id_attr) == as.character(id))
        })

        if(any(id_match)){
          deriv_data <- derivatives_list[[which(id_match)[1]]]
        } else {
          stop(paste("Cannot find derivative data for individual:", id))
        }
      }

      # Extract the appropriate derivatives based on use.derivatives
      use_deriv <- attr(dynEGA.object, "glla")$use.derivatives

      # Get derivative columns
      deriv_cols <- grep(
        paste0(".Ord", use_deriv, "$"),
        colnames(deriv_data)
      )

      if(length(deriv_cols) == 0){
        # Fallback: try to match any Ord columns
        deriv_cols <- grep("Ord", colnames(deriv_data))
      }

      return(deriv_data[, deriv_cols, drop = FALSE])
    }
  }

  if (type == "average") {
    # Function to compute average network
    average_network <- function(ids) {
      nets <- lapply(ids, function(id) {
        ind_result <- dynEGA.object$dynEGA$individual[[id]]
        if(is.null(ind_result$network)){
          warning(paste("No network found for individual:", id))
          return(NULL)
        }
        return(ind_result$network)
      })

      # Remove NULL networks
      nets <- nets[!sapply(nets, is.null)]

      if(length(nets) == 0){
        stop("No valid networks found for averaging")
      }

      Reduce("+", nets) / length(nets)
    }

    # Create plots for each cluster
    plots <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
      ids_in_cluster <- cluster_df$ID[cluster_df$cluster == cl]

      # Get average network
      net_mean <- average_network(ids_in_cluster)

      # Get community memberships for all individuals in cluster
      wc_list <- lapply(ids_in_cluster, function(id) {
        wc <- dynEGA.object$dynEGA$individual[[id]]$wc
        if(is.null(wc)){
          warning(paste("No community membership found for individual:", id))
          return(NULL)
        }
        return(wc)
      })

      # Remove NULL memberships
      wc_list <- wc_list[!sapply(wc_list, is.null)]

      if(length(wc_list) == 0){
        stop(paste("No valid community memberships found for cluster:", cl))
      }

      # Create membership matrix
      wc_mat <- do.call(cbind, lapply(wc_list, as.integer))

      # Compute prototype membership (majority vote)
      wc_proto <- apply(wc_mat, 1, function(x) {
        tab <- table(x)
        as.integer(names(tab)[which.max(tab)])
      })

      # Create plot
      EGAnet:::basic_plot_setup(
        network = net_mean,
        wc = wc_proto,
        title = paste("Cluster", cl),
        node.size = node.size,
        label.size = label.size
      )
    })

    return(do.call(ggpubr::ggarrange, plots))

  } else if (type == "population") {

    # Create plots for each cluster
    plots <- lapply(sort(unique(cluster_df$cluster)), function(cl) {
      ids_in_cluster <- cluster_df$ID[cluster_df$cluster == cl]

      # Stack derivatives
      pooled_data <- tryCatch({
        do.call(rbind, lapply(ids_in_cluster, get_individual_derivatives))
      }, error = function(e){
        stop(paste("Error pooling derivatives for cluster", cl, ":", e$message))
      })

      # Run EGA on pooled derivatives
      pop_res <- tryCatch({
        EGA(pooled_data, plot.EGA = FALSE, verbose = FALSE, ...)
      }, error = function(e){
        stop(paste("Error running EGA for cluster", cl, ":", e$message))
      })

      # Create plot
      EGAnet:::basic_plot_setup(
        network = pop_res$network,
        wc = pop_res$wc,
        title = paste("Cluster", cl),
        node.size = node.size,
        label.size = label.size
      )
    })

    return(do.call(ggpubr::ggarrange, plots))
  }
}
