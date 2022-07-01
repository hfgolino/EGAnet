#' Variation of Information Hierarchical Clustering
#'
#' @description Performs hierarchical k-means clustering using variation of information
#' of the community memberships derived from individuals in \code{\link[EGAnet]{dynEGA}}
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.pop.ind}} object that is used to match the arguments of the EII object.
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
#' \item{id_clusters}{A list containing the IDs for each cluster. 
#' \code{missing} contains IDs that were not able to be included due
#' to missing too many memberships (NA >= 0.50)}
#'
#' \item{optimal}{Optimal clusters determine by *k*-means clustering from
#' \code{\link[factoextra]{fviz_nbclust}}}
#' 
#' \item{hierarchical}{Results of the hierarchical clustering based on the optimal
#' number of clusts. Reports results from \code{\link[factoextra]{hcut}}}
#'
#' \item{plot}{Plot output from results}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Information Theoretic Clustering
# Updated 01.07.2022
infoCluster <- function(
    dynEGA.object,
    plot.cluster = TRUE
)
{
  
  # Check for class
  if(!is(dynEGA.object, "dynEGA.Individuals") & !is(dynEGA.object, "dynEGA.ind.pop")){
    stop(
      paste(
        "Input into the `dynEGA.object` argument's class is not `dynEGA.Individuals` or `dynEGA.ind.pop`.\n\n",
        "Class of dynEGA.object = ", paste(
          class(dynEGA.object), sep = "", collapse = ", "
        ),
        sep = ""
      )
    )
  }else if(is(dynEGA.object, "dynEGA.ind.pop")){
    ind_ega <- dynEGA.object$dynEGA.ind$dynEGA
  }else if(is(dynEGA.object, "dynEGA")){
    ind_ega <- dynEGA.object
  }
  
  # Remove methods
  if("methods" %in% tolower(names(ind_ega))){
    ind_ega <- ind_ega[-which(tolower(names(ind_ega)) == "methods")]
  }
  
  # Obtain memberships
  membership_list <- lapply(
    ind_ega, function(x){
      x$wc
    }
  )
  
  # Obtain all IDs
  all_ids <- names(membership_list)
  
  # Remove memberships where NAs are >= .50 of memberships
  membership_NA <- unlist(
    lapply(
      membership_list, function(x){
        mean(is.na(x)) >= 0.50
      }
    )
  )
  
  # Remove missing memberships
  membership_list <- membership_list[!membership_NA]
  
  # Obtain usable IDs
  ids <- names(membership_list)
  
  # Make grid with IDs
  grid_ids <- expand.grid(x = ids, y = ids)
  
  # Message user
  message("Computing variation of information...", appendLF = FALSE)
  
  # Obtain variation of information
  vi_list <- lapply(1:nrow(grid_ids), function(i){
    
    # Remove NAs
    wc_1 <- na.omit(membership_list[[grid_ids[i,1]]])
    wc_2 <- na.omit(membership_list[[grid_ids[i,2]]])
    
    # Obtain common variables
    common_wc <- intersect(names(wc_1), names(wc_2))

    # Compute variation of information
    igraph::compare(
      comm1 = wc_1[common_wc],
      comm2 = wc_2[common_wc],
      method = "vi"
    )
    
  })
  
  # Message user
  message("done", appendLF = TRUE)
  
  # Convert to matrix
  vi_matrix <- matrix(
    unlist(vi_list),
    nrow = length(ids),
    ncol = length(ids),
    byrow = TRUE
  )
  
  # Name rows and columns
  row.names(vi_matrix) <- ids
  colnames(vi_matrix) <- ids
  
  # # Scale variation of information matrix to be zero to one
  # scaled_vi <- custom.min.max(
  #   vi_matrix, c(0, 1)
  # )
  # 
  # # Reverse so 1 = greater similarity, 0 = no similarity
  # scaled_vi <- 1 - scaled_vi
  # 
  # # Convert to igraph
  # g <- convert2igraph(scaled_vi)
  # 
  # # Apply Louvain algorithm
  # clusters <- cluster_louvain(g)$membership
  
  # Initialize optimal clusters
  optimal_clusters <- list()
  class(optimal_clusters) <- "try-error"
  
  # Initialize maximum clusters
  cluster.max <- length(ids)
  
  # Message user
  message("Searching for optimal clusters...", appendLF = FALSE)
  
  # Loop through maximum clusters
  while(is(optimal_clusters, "try-error")){
    
    # Optimal clusters using k-means
    optimal_clusters <- try(
      factoextra::fviz_nbclust(
        scaled_vi,
        FUNcluster = kmeans,
        method = "silhouette",
        k.max = cluster.max,
        nboot = 500
      ),
      silent = TRUE
    )
    
    # Subtract 1 from max clusters
    if(is(optimal_clusters, "try-error")){
      cluster.max <- cluster.max - 1 
    }
    
  }
  
  # Message user
  message("done", appendLF = TRUE)
  
  # Obtain optimal clusters
  optimal_number <- as.numeric(
    as.character(
      optimal_clusters$data$clusters[
        which.max(optimal_clusters$data$y)
      ]
    )
  )
  
  # Compute hierarchical clustering and cut into maximum clusters
  cluster_result <- factoextra::hcut(
    scaled_vi, k = optimal_number, stand = TRUE
  )
  
  # Visualize
  cluster_plot <- suppressWarnings(
    factoextra::fviz_dend(
      cluster_result, rect = TRUE, cex = 0.5,
      k_colors = color_palette_EGA(
        "polychrome", wc = 1:cluster.max
      )
    )
  )
  
  # Set up optimal cluster plot
  optimal_clusters <- optimal_clusters +
    ggplot2::scale_x_discrete(
      breaks = as.character(
        sort(
          c(
            seq(optimal_number, 1, -floor(length(ids) / 10)),
            optimal_number,
            seq(optimal_number, length(ids), floor(length(ids) / 10)) 
          )
        )
      )
    )
  
  # Organize plots
  cluster_plot_arrange <- ggpubr::ggarrange(
    optimal_clusters,
    cluster_plot
  )
  
  # Plot
  if(isTRUE(plot.cluster)){
    cluster_plot_arrange
  }
  
  # Obtain clusters of IDs
  id_cluster_list <- list()
  
  # Loop through to add IDs
  for(i in 1:optimal_number){
    
    id_cluster_list[[as.character(i)]] <-
      names(cluster_result$cluster)[
        cluster_result$cluster == i
      ]
    
  }
  
  # Add missing IDs to list
  id_cluster_list$missing <- setdiff(all_ids, ids)
  
  # Return results
  results <- list()
  results$id_clusters <- id_cluster_list
  results$optimal <- optimal_clusters
  results$hierarhical <- cluster_result
  results$plot <- cluster_plot_arrange
  
  # Add class
  class(results) <- "infoCluster"
  
  return(results)
  
  
}