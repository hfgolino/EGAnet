#' @title Information Theoretic Mixture Clustering for \code{\link[EGAnet]{dynEGA}}
#'
#' @description Performs hierarchical clustering using Jensen-Shannon distance
#' followed by the Louvain algorithm with consensus clustering. The method
#' iteratively identifies smaller and smaller clusters until there is no
#' change in the clusters identified
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object that is used to match 
#' the arguments of the EII object
#' 
#' @param plot.cluster Boolean (length = 1).
#' Should plot of optimal and hierarchical clusters be output?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not plot
#'
#' @examples
#'# Obtain data
#' sim.dynEGA <- sim.dynEGA # bypasses CRAN checks
#' 
#' \dontrun{
#' # Dynamic EGA individual and population structure
#' dyn.ega1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 25, use.derivatives = 1, 
#'   ncores = 2, corr = "pearson"
#' )
#' 
#' # Perform information-theoretic clustering
#' clust1 <- infoCluster(dynEGA.object = dyn.ega1)}
#'
#' @return Returns a list containing:
#' 
#' \item{clusters}{A vector corresponding to cluster each participant belongs to}
#' 
#' \item{clusterTree}{The dendogram from \code{\link[stats]{hclust}} the hierarhical clustering}
#'
#' \item{clusterPlot}{Plot output from results}
#' 
#' \item{JSD}{Jensen-Shannon Distance}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#' 
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#' 
#' @export
#' 
# Information Theoretic Clustering for dynEGA
# Updated 12.11.2023
infoCluster <- function(dynEGA.object, plot.cluster = TRUE)
{
  
  # Send experimental message (for now)
  experimental("infoCluster")
  
  # Check for appropriate class ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "infoCluster")
  }
  
  # Get proper objects (if not, send an error)
  dynega_objects <- get_dynEGA_object(dynEGA.object)
  
  # Get individual networks 
  individual_networks <- lapply(
    dynega_objects$individual, function(x){x$network}
  )
  
  # Get pairwise JSD
  jsd_matrix <- pairwise_spectral_JSD(individual_networks)

  # Make diagonal NA
  diag(jsd_matrix) <- NA
  
  # Remove all NAs
  rm_cols <- !lvapply(
    as.data.frame(jsd_matrix), function(x){all(is.na(x))}
  )
    
  # Remove missing data points and convert to distance
  jsd_matrix <- jsd_matrix[rm_cols, rm_cols]
  
  # Get similarity matrix
  jss_matrix <- 1 - jsd_matrix
  
  # Make diagonal 0 again
  diag(jss_matrix) <- diag(jsd_matrix) <- 0

  # Use agglomerative Ward's D
  hier_clust <- hclust(
    d = as.dist(jsd_matrix),
    method = "ward.D2"
  )

  # Get number of cases
  cases <- dim(jsd_matrix)[2]

  # Get cut sequence
  cut_sequence <- seq_len(cases)[-1]

  # Obtain cuts
  hier_cuts <- lapply(cut_sequence, function(i){
    cutree(hier_clust, i)
  })

  # Compute modularity over solutions
  Qs <- nvapply(
    hier_cuts, function(cut){
      modularity(jss_matrix, cut, resolution = 1.01)
    }
  )
  
  # Switch based on positive modularity
  Q_index <- which.max(Qs)
  
  # Obtain clusters
  clusters <- hier_cuts[[Q_index]]
  
  # Set up results
  results <- list(
    clusters = clusters,
    modularity = Qs[Q_index],
    clusterTree = hier_clust,
    JSD = jsd_matrix
  )

  # Set class
  class(results) <- "infoCluster"
  
  # Check for plot
  if(isTRUE(plot.cluster)){
    
    # Get plot
    results$plot_cluster <- plot(results)
    
    # Actually send plot
    silent_plot(results$plot_cluster)
    
  }
  
  # Return results
  return(results)
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 14.07.2023
print.infoCluster <- function(x, ...)
{
 
  # Print clusters
  cat("Number of cases: ", length(x$clusters), "\n")
  cat("Number of clusters: ", unique_length(x$clusters))
  
  # Add breakspace
  cat("\n\n")
  
  # Print cluster assignments
  print(x$clusters)

}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 14.07.2023
summary.infoCluster <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method ----
# Works fast enough, so leaving as original code
# Updated 16.11.2023
plot.infoCluster <- function(x, label_size = 3, ...)
{
  
  # Get maximum clusters
  max_clusters <- max(x$clusters)
  
  # Get cluster sequence
  cluster_sequence <- seq_len(max_clusters)

  # Get color palette
  color_palette <- color_palette_EGA("polychrome", wc = rev(cluster_sequence))
  
  # Get dendrogram and color branches
  cluster_data <- dendextend::as.ggdend(
    dendextend::color_branches(
      as.dendrogram(x$clusterTree),
      clusters =  x$clusters[x$clusterTree$order],
      col = color_palette
    )
  )
  
  # Change NA to grey
  cluster_data$segments$col[is.na(cluster_data$segments$col)] <- "grey"
  
  # Update color palette
  color_palette <- c("grey", rev(color_palette))
  
  # Factor segments by color
  cluster_data$segments$col <- factor(
    cluster_data$segments$col, levels = color_palette
  )

  # Set up plot
  cluster_plot <- ggplot2::ggplot() + 
    ggplot2::geom_segment(
      data = cluster_data$segment,
      ggplot2::aes(
        x = x, y = y, 
        xend = xend, yend = yend, 
        color = col
      )
    ) +
    ggplot2::geom_text(
      data = cluster_data$labels,
      ggplot2::aes(x, y, label = label, hjust = 0),
      size = label_size
    ) +
    ggplot2::scale_color_manual(
      labels = c("", cluster_sequence),
      values = color_palette
    ) +
    ggplot2::coord_flip() + 
    ggplot2::scale_y_reverse(expand = c(0.2, 0)) + 
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = "Cluster")
    )
 
  # Remove clusters if none
  if(all(x$clusters == ncol_sequence(x$JSD))){
    cluster_plot <- cluster_plot +
      ggplot2::theme(legend.position = "none")
  }
  
  # Return plot
  return(cluster_plot)
  
}

#' @noRd
# Global variables needed for CRAN checks ----
# Updated 17.11.2023
utils::globalVariables(c("x", "y", "xend", "yend", "cluster", "label")) 



