#' Information Theoretic Mixture Clustering for \code{\link[EGAnet]{dynEGA}}
#'
#' @description Performs hierarchical clustering using Jensen-Shannon distance
#' followed by the Louvain algorithm with consensus clustering. The method
#' iteratively identifies smaller and smaller clusters until there is no
#' change in the clusters identified
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object that is used to match the arguments of the EII object.
#' 
#' @param plot.cluster Boolean.
#' Should plot of optimal and hierarchical clusters be output?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not plot
#'
#' @examples
#'
#' \dontrun{
#' \donttest{
#' # Perform dynEGA
#' dyn1 <- dynEGA.ind.pop(data = sim.dynEGA[,-c(22)], n.embed = 5, tau = 1,
#'                       delta = 1, id = 21, use.derivatives = 1,
#'                     model = "glasso", ncores = 2, corr = "pearson")
#'
#' # Perform hierarchical clustering
#' clust1 <- infoCluster(
#'   dynEGA.object = dyn1,
#'   plot.cluster = TRUE
#' )
#' }}
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
#' @importFrom stats hclust
#' 
#' @export
# Information Theoretic Clustering for dynEGA
# Updated 16.07.2022
infoCluster <- function(dynEGA.object, plot.cluster = TRUE)
{
  
  # Check for class
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    stop(
      paste(
        "Input into the `dynEGA.object` argument's class is not `dynEGA` or `dynEGA.ind.pop`.\n\n",
        "Class of dynEGA.object = ", paste(
          class(dynEGA.object), sep = "", collapse = ", "
        ),
        sep = ""
      )
    )
  }else if(is(dynEGA.object, "dynEGA.ind.pop")){
    dynEGA.pop <- dynEGA.object$dynEGA.pop
  }else if(is(dynEGA.object, "dynEGA")){
    dynEGA.pop <- dynEGA.object
  }
  
  # Obtain individual dynEGA objects only
  dynEGA.ind <- dynEGA.object$dynEGA.ind$dynEGA
  
  # Remove methods from dynEGA.ind
  if("methods" %in% tolower(names(dynEGA.ind))){
    dynEGA.ind <- dynEGA.ind[-which(tolower(names(dynEGA.ind)) == "methods")]
  }
  
  # Obtain IDs
  IDs <- names(dynEGA.ind)
  
  # Obtain networks
  networks <- lapply(dynEGA.ind, function(x){
    x$network
  })
  
  # Message user
  message("Computing Jensen-Shannon Distance...\n", appendLF = FALSE)

  # Initialize JSD matrix
  jsd_matrix <- matrix(
    0,
    nrow = length(networks),
    ncol = length(networks)
  )
  
  # Set up progess bar
  pb <- txtProgressBar(
    max = length(networks),
    style = 3
  )
  
  # Populate JSD matrix
  for(i in 2:length(networks)){
    
    for(j in 1:(i-1)){
      
      # Try
      jsd_value <- try(
        jsd(
          networks[[i]], networks[[j]]
        ),
        silent = TRUE
      )
      
      # Check if value is OK
      if(!is(jsd_value, "try-error")){
        jsd_matrix[i,j] <- jsd(
          networks[[i]], networks[[j]]
        )
      }else{
        jsd_matrix[i,j] <- NA
      }
      
    }
    
    # Update progress bar
    setTxtProgressBar(pb, i)
    
  }
  
  # Close progress bar
  close(pb)
  
  # Make symmetric
  jsd_sym <- jsd_matrix + t(jsd_matrix)
  
  # Add names
  colnames(jsd_sym) <- names(networks)
  row.names(jsd_sym) <- names(networks)
  
  # Make jsdist
  jsdist <- jsd_sym
  
  # Make diagonal NA
  diag(jsdist) <- NA
  
  # Remove all NAs
  rm_cols <- apply(jsdist, 2, function(x){all(is.na(x))})
  
  # Remove missing data points
  jsdist <- jsdist[!rm_cols, !rm_cols]
  
  # Make diagonal 0 again
  diag(jsdist) <- 0
  
  # Compute Louvain
  consensus <- most_common_consensus(
    1 - jsdist,
    order = "lower",
    consensus.iter = 1000
  )$most_common
  
  # Unique consensus
  unique_consensus <- length(na.omit(unique(consensus)))
  
  # Perform hierarchical clustering
  hier_clust <- hclust(as.dist(jsdist))
  
  # Check for single cluster
  if(
    unique_consensus == 1 | # consensus = 1 OR
    unique_consensus == length(consensus) # consensus all individuals
  ){
    
    # Obtain clusters
    clusters <- walktrap$membership
    names(clusters) <- colnames(jsdist)
    
  }else{
    
    # Initialize silhouette vector
    silhouette_vec <- numeric(length = ncol(jsdist) - 1)
    
    # Make names the number of clusters
    names(silhouette_vec) <- 2:ncol(jsdist)
    
    # Loop through cuts
    for(i in 2:length(silhouette_vec)){
      
      # Compute silhouette
      hier_silho <- cluster::silhouette(
        x = cutree(hier_clust, i),
        dist = as.dist(jsdist)
      )
      
      # Obtain summary
      silho_summ <- summary(hier_silho)
      
      # Obtain average silhouette
      silhouette_vec[i-1] <- mean(silho_summ$clus.avg.widths)
      
    }
    
    # Obtain maximum average silhouette
    optimal_cut <- as.numeric(names(which.max(silhouette_vec)))
    
    # Obtain clusters
    clusters <- cutree(hier_clust, optimal_cut)
    
  }
  
  # Obtain optimal silhouette
  optimal_silho <- cluster::silhouette(
    x = clusters,
    dist = as.dist(jsdist)
  )
  
  # Convert for ggplot2
  cluster_data <- ggdendro::dendro_data(
    hier_clust
  )
  
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
  
  # Set colors
  cluster_data$segments$cluster <- cluster_data$labels$cluster[
    match(
      floor(cluster_data$segments$x),
      cluster_data$labels$x
    )
  ]
  
  # Set up plot
  cluster_plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = ggdendro::segment(cluster_data),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = cluster)
    ) +
    ggplot2::geom_text(
      data = ggdendro::label(cluster_data),
      ggplot2::aes(x, y, label = label, hjust = 0), 
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = color_palette_EGA(
        "polychrome", wc = 1:max(clusters)
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
      legend.key = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = "Cluster")
    )
  
  # Check if plot should be plotted
  if(isTRUE(plot.cluster)){
    plot(cluster_plot)
  }
  
  ## Return data
  results <- list()
  results$clusters <- clusters
  results$clusterTree <- hier_clust
  results$clusterPlot <- cluster_plot
  results$JSD <- jsdist
  
  ## Set class
  class(results) <- "infoCluster"
  
  return(results)
  
}
