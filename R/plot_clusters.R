#' @title Plot Clustered Individual Networks
#'
#' @description Visualize clusters of individual networks identified by
#' \code{\link[EGAnet]{infoCluster}} (or another clustering method) applied to a
#' \code{\link[EGAnet]{dynEGA}} object
#'
#' Provides visualization modes:
#' \itemize{
#'   \item \code{"population"}: Cluster-level "population" networks, obtained by stacking derivatives across individuals and estimating a new EGA
#'   \item \code{"average"}: Cluster-average networks, obtained by averaging adjacency matrices across individuals in the cluster
#' }
#'
#' @param dynEGA.object A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object
#'
#' @param clustering Vector (length of individuals).
#' A vector of cluster membership for each individual in the \code{dynEGA.object}.
#' A common and easy option is to use an \code{infoCluster} object containing cluster
#' assignments for individuals (see \code{\link[EGAnet]{infoCluster}}).
#' Accepts any vector of memberships
#'
#' @param include Numeric vector. Specifies which cluster memberships should be
#' explicitly included in the plot. By default, all clusters are shown. Use this
#' argument to restrict the visualization to a subset of clusters (e.g.,
#' `include = c(1, 3, 5)` will only display clusters 1, 3, and 5)
#'
#' @param type Character (length = 1).
#' Type of visualization to produce:
#'
#' \itemize{
#'
#' \item \code{"population"} (default) --- Cluster-level population networks by stacking
#' derivatives across individuals in the cluster and re-estimating with \code{\link[EGAnet]{EGA}}
#'
#' \item \code{"average"} --- Averaged networks with a consensus clustering matrix (pairwise)
#' membership similarity that is supplied to \code{\link[EGAnet]{community.detection}}
#'
#' }
#'
#' In either type, the argument setting applied in your `dynEGA.object` will be carried
#' forward into the network estimation (`"population"`) and community detection (both)
#'
#' @param node.size Numeric (length = 1 or number of variables).
#' Node size in network visualizations.
#' Defaults to \code{3}
#'
#' @param label.size Numeric (length = 1 or number of nodes.
#' Label size in network visualizations.
#' Defaults to \code{1.2}
#'
#' @param ... Additional arguments passed to
#' \code{\link[EGAnet]{EGA}} and
#' \code{\link[EGAnet]{compare.EGA.plots}}
#'
#' @return
#' Plots created using \code{\link[EGAnet]{compare.EGA.plots}}:
#'
#' \item{"population"}{A plot combining cluster-level population network plots}
#'
#' \item{"average"}{A plot combining cluster-average network plots}
#'
#' @details
#' This function provides flexible visualization of individual clusters obtained from
#' \code{\link[EGAnet]{infoCluster}} (or other clustering input; see examples):
#'
#' \itemize{
#'
#' \item Cluster-level population networks treat each cluster as a "mini-population," stacking
#' derivatives across individuals to re-estimate an EGA
#'
#' \item Average networks show the mean connectivity structure for each cluster,
#' with consensus clustering memberships
#'
#' }
#'
#' Function automatically extracts the appropriate derivatives (\code{"population"}) and
#' networks (\code{"average"}) from the \code{dynEGA.object} results.
#'
#' @examples
#' # Load data
#' data <- sim.dynEGA
#'
#' \dontrun{
#' # Run dynEGA with optimization
#' optimized_all <- dynEGA(
#'   data = data, level = c("individual", "population"),
#'   n.embed = 3:25, n.embed.optimize = TRUE
#' )
#'
#' # Cluster individuals
#' clust <- infoCluster(dynEGA.object = optimized_all)
#'
#' # Cluster-level population networks
#' plot_clusters(
#'   dynEGA.object = optimized_all,
#'   clustering = clust,
#'   type = "population"
#' )
#'
#' # Average networks per cluster
#' plot_clusters(
#'   dynEGA.object = optimized_all,
#'   clustering = clust,
#'   type = "average"
#' )
#'
#' # Cluster-level population networks, including only Cluster 2:
#' plot_clusters(
#'   dynEGA.object = optimized_all,
#'   clustering = clust, include = 2,
#'   type = "population"
#' )
#'
#' # Using alternative clusters
#' plot_clusters(
#'   dynEGA.object = optimized_all,
#'   clustering = rep(1:2, each = 50), # vector of memberships
#'   type = "population"
#')
#'
#' # Run with non-optimized dynEGA
#' standard_all <- dynEGA(
#'   data = data,
#'   level = c("individual", "population")
#' )
#'
#' # Obtain clusters
#' clust_standard <- infoCluster(dynEGA.object = standard_all)
#'
#' # Plot clusters with population
#' plot_clusters(
#'   dynEGA.object = standard_all,
#'   clustering = clust_standard,
#'   type = "population"
#' )}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @seealso \code{\link[EGAnet]{dynEGA}}, \code{\link[EGAnet]{infoCluster}}, \code{\link[EGAnet]{EGA}},
#' \code{\link[EGAnet]{compare.EGA.plots}}
#'
#' @export
#'
# Plots EGA Sub-Groups ----
# Updated 18.11.2025
plot_clusters <- function(
    dynEGA.object, clustering, include, type = c("population", "average"),
    node.size = 3, label.size = 1.2, ...
)
{

  # Check for missing arguments (argument, default, function)
  type <- set_default(type, "population", plot_clusters)

  # Collect inputs
  inputs <- plot_clusters_errors(
    dynEGA.object, clustering, include, node.size, label.size
  )

  # Obtain EGA settings
  uni.method <- attributes(inputs$dynega_objects$dynEGA$population)$unidimensional$uni.method
  network_methods <- attributes(inputs$dynega_objects$dynEGA$population$network)$methods
  community_methods <- tolower(
    attributes(inputs$dynega_objects$dynEGA$population$wc)$methods$algorithm
  )

  # Switch based on type
  if(type == "population"){

    # Split derivatives by clusters
    clustered_data <- lapply(
      inputs$include, function(group){
        do.call(
          rbind, inputs$dynega_objects$Derivatives$Estimates[inputs$clustering == group]
        )
      }
    )

    # Apply EGA to clustered data
    clustered_ega <- lapply(
      clustered_data, EGA, corr = network_methods$corr,
      na.data = network_methods$na.data,
      model = network_methods$model, algorithm = community_methods,
      uni.method = uni.method, plot.EGA = FALSE, ...
    )

  }else{

    # Separate networks into lists of groups
    group_list <- lapply(
      inputs$include, function(group){
        lapply(
          inputs$dynega_objects$dynEGA$individual[inputs$clustering == group],
          function(x){x$network}
        )
      }
    )

    # Obtain average networks
    average_networks <- lapply(
      group_list, symmetric_matrix_lapply, function(x){
        z2r(mean(r2z(x))) # Convert to Fisher's z, average, convert back to correlations
      }
    )

    # Obtain memberships
    memberships <- lapply(
      inputs$include, function(group){
        do.call(
          rbind, lapply(
            inputs$dynega_objects$dynEGA$individual[inputs$clustering == group],
            function(x){x$wc}
          )
        )
      }
    )

    # Obtain consensus matrices
    consensus <- lapply(memberships, create_consensus)

    # Set up EGA objects for comparison
    clustered_ega <- lapply(inputs$include, function(i){

      # EGA object
      ega_object <- list(
        network = average_networks[[i]],
        wc = community.detection(
          consensus[[i]], algorithm = community_methods
        )

      )

      # Set class
      class(ega_object) <- "EGA"

      # Return object
      return(ega_object)

    })

  }

  # Compare plots
  clustered_plot <- compare.EGA.plots(
    input.list = clustered_ega, labels = paste("Cluster", inputs$include), ...
  )

  # Return output
  return(clustered_plot)

}

#' @noRd
# Errors ----
# Updated 17.11.2025
plot_clusters_errors <- function(
    dynEGA.object, clustering, include, type = c("population", "average"),
    node.size = 3, label.size = 1.2, ...
)
{

  # Check for "dynEGA" class
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "plot_clusters")
  }

  # Ensure population and individual
  if(!all(c("population", "individual") %in% names(dynEGA.object$dynEGA))){
    stop(
      "'dynEGA.object' must include both \"population\" and \"individual\" levels.",
      call. = FALSE
    )
  }

  # Check for clustering
  if(missing(clustering)){
    stop("Input for 'clustering' must be provided", call. = FALSE)
  }

  # Check for `infoCluster` object
  if(is(clustering, "infoCluster")){
    clustering <- clustering$clusters
  }

  # Ensure clustering can be vector
  clustering <- as.vector(clustering)

  # Obtain lengths
  n_individuals <- length(dynEGA.object$dynEGA$individual)
  n_clustering <- length(clustering)

  # Check for same length
  if(n_individuals != n_clustering){
    stop(
      paste0(
        "Number of individuals (", n_individuals, "), is not the same ",
        "length as the clustering elements (", n_clustering, ")"
      ),
      call. = FALSE
    )
  }

  # Check for missing include
  if(missing(include)){
    include <- unique(clustering)
  }

  # Ensure that all are available in 'include'
  if(!all(include %in% clustering)){
    stop(
      paste0(
        "Some clusters specified in 'include' are not available in 'clustering': ",
        paste0(setdiff(include, clustering), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # 'node.size' errors
  typeof_error(node.size, "numeric", "plot_clusters")

  # 'label.size' errors
  typeof_error(label.size, "numeric", "plot_clusters")

  # Return all output
  return(
    list(
      dynega_objects = dynEGA.object,
      clustering = clustering,
      include = include,
      node.size = node.size,
      label.size = label.size
    )
  )

}

#' @noRd
# Create consensus matrix ----
# Updated 17.11.2025
create_consensus <- function(membership_matrix)
{

  # Obtain dimensions
  dimensions <- dim(membership_matrix)

  # Initialize consensus matrix
  consensus_matrix <- matrix(
    nrow = dimensions[2], ncol = dimensions[2]
  )

  # Loop over to get proportions
  for(i in seq_len(dimensions[2])){
    for(j in i:dimensions[2]){

      # Fill consensus matrix
      consensus_matrix[i,j] <- consensus_matrix[j,i] <- sum(
        membership_matrix[,i] == membership_matrix[,j],
        na.rm = TRUE
      )

    }
  }

  # Divide by number of people
  return(consensus_matrix / dimensions[1])

}
