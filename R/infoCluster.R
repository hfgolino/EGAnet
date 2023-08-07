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
# Updated 03.08.2023
infoCluster <- function(dynEGA.object, plot.cluster = TRUE)
{
  
  # Send experimental message (for now)
  experimental("infoCluster")
  
  # Check for appropriate class ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA")
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
  rm_cols <- lvapply(
    as.data.frame(jsd_matrix), function(x){all(is.na(x))}
  )
    
  # Remove missing data points
  jsd_matrix <- jsd_matrix[!rm_cols, !rm_cols]
  
  # Get similarity matrix
  jss_matrix <- 1 - jsd_matrix

  # Make diagonal 0 again
  diag(jss_matrix) <- diag(jsd_matrix) <- 0
  
  # Perform hierarchical clustering (follows Walktrap)
  hier_clust <- hclust(
    d = as.dist(jsd_matrix),
    method = "complete"
  )
  
  # Get number of cases
  cases <- dim(jsd_matrix)[2]

  # Get cut sequence
  cut_sequence <- seq_len(cases)

  # Obtain cuts
  hier_cuts <- lapply(cut_sequence, function(i){
    cutree(hier_clust, i)
  })

  # Name cuts
  names(hier_cuts) <- cut_sequence

  # Make any cuts with singleton clusters NULL
  remaining_cuts <- !lvapply(
    hier_cuts, function(x){any(table(x) == 1)}
  )

  # Obtain cuts
  cuts <- as.numeric(names(remaining_cuts)[remaining_cuts])

  # Get modularities
  Qs <- nvapply(
    cuts, function(cut_value){
      modularity(jss_matrix, cutree(hier_clust, cut_value))
    }
  )
    
  # Maybe just use Walktrap?
  # clusters <- remove_attributes(
  #   community.detection(jss_matrix, algorithm = "walktrap")
  # )
  
  # Possibility for Louvain with consensus?
  # clusters <- remove_attributes(
  #   community.consensus(jss_matrix)
  # )
  
  # Obtain clusters
  clusters <- cutree(hier_clust, cuts[which.max(Qs)])
  
  # Check if single cluster
  if(unique_length(clusters) == 1){

    # Get number of nodes to initialize matrices
    nodes <- dim(individual_networks[[1]])[2]
    
    # Get indices of upper triangle
    upper_indices <- which(upper.tri(diag(nodes)))
    
    # Generate random networks
    random_networks <- lapply(individual_networks, function(network){

      # Initialize new matrix
      new_network <- matrix(0, nrow = nodes, ncol = nodes)
      
      # Set shuffled indices up to edges to 1
      new_network[
        shuffle(upper_indices, size = edge_count(network))
      ] <- 1
      
      # Make network symmetric
      return(new_network + t(new_network))
    
    })
    
    # Get the random JSD matrix
    jsd_random_matrix <- pairwise_spectral_JSD(random_networks)
    
    # Make diagonal NA
    diag(jsd_random_matrix) <- NA
    
    # Remove all NAs
    rm_cols <- lvapply(
      as.data.frame(jsd_random_matrix), function(x){all(is.na(x))}
    )
    
    # Remove missing data points
    jsd_random_matrix <- jsd_random_matrix[!rm_cols, !rm_cols]

    # Make diagonal 0 again
    diag(jsd_random_matrix) <- 0

    # Compare to empirical
    comparison <- t.test(
      jsd_matrix[upper_indices],
      jsd_random_matrix[upper_indices],
      paired = TRUE,
      var.equal = FALSE
    )
    
    # Obtain sign of statistic
    comparison_sign <- sign(comparison$statistic)
    
    # Compute adaptive alpha
    adaptive_p <- adapt.a(
      test = "paired",
      n = length(upper_indices),
      alpha = .001,
      power = 0.80,
      efxize = "large"
    )
    
    # Check for empirical JSD > random JSD OR non-significant t-test
    if(comparison_sign == 1 | comparison$p.value > adaptive_p$adapt.a){
      
      # Set clusters to all individuals
      clusters <- cut_sequence
      names(clusters) <- dimnames(jsd_matrix)[[2]]
      
    }
    
    # Compile results
    single_cluster <- list(
      JSD_random = jsd_random_matrix,
      t.test = comparison,
      adaptive.p.value = adaptive_p,
      d = d(
        jsd_matrix[upper_indices],
        jsd_random_matrix[upper_indices],
        paired = TRUE
      )
    )

  }
  
  # Set up results
  results <- list(
    clusters = clusters,
    modularity = Qs[which.max(Qs)],
    clusterTree = hier_clust,
    JSD = jsd_matrix
  )
  
  # Check for single cluster test
  if(exists("single_cluster")){
    results$single.cluster.test <- single_cluster
  }

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
# Updated 13.07.2023
plot.infoCluster <- function(x, ...)
{
  
  # Prepare data for {ggplot2}
  cluster_data <- ggdendro::dendro_data(x$clusterTree)
  
  # Get clusters
  clusters <- x$clusters
  
  # Create data frame
  cluster_df <- data.frame(
    label = names(clusters),
    cluster = factor(clusters)
  )
  
  # Merge data
  cluster_data$labels <- merge(
    cluster_data$labels, cluster_df,
    by = "label"
  )
  
  # Split dendrogram into upper grey section and lower coloured section
  cut <- unique_length(clusters)
  height <- unique(cluster_data$segments$y)[order(unique(cluster_data$segments$y), decreasing = TRUE)]
  cut.height <- mean(c(height[cut], height[cut-1]))
  cluster_data$segments$line <- swiftelse(cluster_data$segments$y == cluster_data$segments$yend &
                                         cluster_data$segments$y > cut.height, 1, 2)
  cluster_data$segments$line <- swiftelse(cluster_data$segments$yend  > cut.height, 1, cluster_data$segments$line)
  
  # Number the clusters
  cluster_data$segments$cluster <- c(-1, diff(cluster_data$segments$line))
  change <- which(cluster_data$segments$cluster == 1)
  for (i in 1:cut) cluster_data$segments$cluster[change[i]] = i + 1
  cluster_data$segments$cluster <-  swiftelse(cluster_data$segments$line == 1, 1, 
                                           swiftelse(cluster_data$segments$cluster == 0, NA, cluster_data$segments$cluster))
  
  
  # Replace NA values in cluster
  if(!all(is.na(cluster_data$segments$cluster))){
    for(i in seq_along(cluster_data$segments$cluster)){
      
      if(is.na(cluster_data$segments$cluster[i])){
        cluster_data$segments$cluster[i] <- cluster_data$segments$cluster[i-1]
      }
      
    }
  }else{
    cluster_data$segments$cluster <- -1
  }
  
  
  # Consistent numbering between segment$cluster and label$cluster
  cluster_df$label <- factor(cluster_df$label, levels = cluster_data$labels$label)
  cluster_df$cluster <- factor((cluster_df$cluster), levels = unique(cluster_df$cluster), labels = (1:cut) + 1)
  cluster_data[["labels"]] <- merge(cluster_data[["labels"]], cluster_df, by = "label")
  
  # Positions for cluster labels
  n.rle <- rle(cluster_data$segments$cluster)
  N <- cumsum(n.rle$lengths)
  N <- N[seq(1, length(N), 2)] + 1
  N.df <- cluster_data$segments[N, ]
  N.df$cluster <- N.df$cluster - 1
  
  # Check for all the same cluster
  if(all(cluster_data$segments$cluster == -1)){
    cluster_data$segments$cluster <- rep(
      1, length(cluster_data$segments$cluster)
    )
  }
  
  # Ensure clusters are factors
  cluster_data$segments$cluster <- as.factor(
    cluster_data$segments$cluster
  )
  
  # Make labels
  if(max(clusters) == 1){
    label <- "1"
  }else{
    label <- c("", 1:max(clusters))
  }
  
  # Set up plot
  cluster_plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = cluster_data$segment,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = cluster)
    ) +
    ggplot2::geom_text(
      data = cluster_data$label,
      ggplot2::aes(x, y, label = label, hjust = 0),
      size = 3
    ) +
    ggplot2::geom_text(
      data = N.df,
      ggplot2::aes(
        x = x, y = y, label = factor(cluster),
        colour = factor(cluster + 1)
      ),
      hjust = 1.5, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      labels = label,
      values = c(
        "grey", color_palette_EGA(
              "polychrome", wc = 1:max(clusters)
        )
      )
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
  if(all(clusters == ncol_sequence(x$JSD))){
    cluster_plot <- cluster_plot +
      ggplot2::theme(
        legend.position = "none"
      )
  }
  
  # Return plot
  return(cluster_plot)
  
}

#' @noRd
# Global variables needed for CRAN checks ----
# Updated 04.08.2023
utils::globalVariables(c("x", "y", "xend", "yend", "cluster")) 



