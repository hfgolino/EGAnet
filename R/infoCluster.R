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
#' \item{clusterTree}{The dendogram data frame showing the hierarhical clustering}
#'
#' \item{clusterPlot}{Plot output from results}
#' 
#' \item{JSS}{Jensen-Shannon Similarity based on 1 - Jensen-Shannon Distance}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Information Theoretic Clustering for dynEGA
# Updated 10.07.2022
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
  
  # Missing parallelization
  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
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
  
  # Get similarity
  jss <- 1 - jsd_sym
  
  # Add names
  colnames(jss) <- names(networks)
  row.names(jss) <- names(networks)
  
  # Make diagonal NA
  diag(jss) <- NA
  
  # Remove all NAs
  rm_cols <- apply(jss, 2, function(x){all(is.na(x))})
  
  # Remove missing data points
  jss <- jss[!rm_cols, !rm_cols]
  
  # Make diagonal 1 again
  diag(jss) <- 1
  
  # Remove values -1 > x > 1
  rm_cols <- apply(jss, 2, function(x){any(abs(x) > 1)})
  
  # Remove missing data points
  jss <- jss[!rm_cols, !rm_cols]
  
  # Message user
  message("Obtaining clusters...", appendLF = FALSE)
  
  # Obtain cluster list
  cluster_list <- list()
  
  # Obtain clusters
  cluster_list[[1]] <- most_common_consensus(
    jss,
    order = "lower",
    consensus.iter = 1000
  )$most_common
  
  # Set count
  counter <- 2
  
  # Obtain finer clusters
  while(TRUE){
    
    # Obtain unique clusters
    unique_clusters <- unique(cluster_list[[counter - 1]])
    
    # Initialize next clusters
    cluster_list[[counter]] <- cluster_list[[counter - 1]]
    
    # Initialize cluster number to add
    cluster_add <- 0
    
    # Apply clustering to clusters
    for(i in unique_clusters){
      
      # Obtain index
      index <- which(
        cluster_list[[counter - 1]] == i
      )
      
      # Skip over index if not matrix
      if(is.matrix(jss[index, index])){
        # Obtain lower clusters
        cluster_list[[counter]][index] <- most_common_consensus(
          jss[index, index],
          order = "lower",
          consensus.iter = 1000
        )$most_common + cluster_add
      }
      
      # Increase cluster number to add
      cluster_add <- max(cluster_list[[counter]][index])
    
    }
    
    # Break when all are equal
    if(all(cluster_list[[counter]] == cluster_list[[counter - 1]])){
      cluster_list <- cluster_list[-counter]
      break
    }
    
    # Increase count
    counter <- counter + 1
    
  }
  
  # Message user
  message("done", appendLF = TRUE)
  
  # Loop through cluster list
  # Replace all unique with one cluster
  cluster_list <- lapply(cluster_list, function(x){
    if(length(x) == length(unique(x))){
      one <- rep(1, length(x))
      names(one) <- names(x)
      x <- one
    }
    return(x)
  })

  ## Initialize tree matrix
  cluster_tree <- data.frame(
    cluster0 = 0,
    cluster1 = cluster_list[[1]]
  )
  
  ## Populate tree matrix
  if(length(cluster_list) > 1){
    for(i in 2:length(cluster_list)){
      
      # Make another tree
      cluster_tree[[paste("cluster", i, sep = "")]] <- cluster_list[[i]]
      
    }
  }
  
  ## Add ID
  cluster_tree$id = names(cluster_list[[1]])
  row.names(cluster_tree) <- NULL
  
  ## Initialize path string
  cluster_tree$pathString <- 0
  
  ## Add path string
  for(i in 1:nrow(cluster_tree)){
    
    cluster_tree$pathString[i] <- paste(
      cluster_tree[i,-ncol(cluster_tree)], sep = "", collapse = "/"
    )
    
  }

  ## Set as node
  cluster_node <- data.tree::as.Node(cluster_tree)
  
  ## Convert to phylo tree
  phylo_tree <- ape::as.phylo(
    cluster_node
  )
  
  ## Check for plot
  if(isTRUE(plot.cluster)){
    plot(phylo_tree)
  }
  
  
  # Prepare data for results
  ## Remove first and last column of cluster tree
  cluster_tree <- cluster_tree[,-c(1, ncol(cluster_tree))]
  
  ## Move ID to front
  cluster_tree <- cluster_tree[,c("id", colnames(cluster_tree)[-ncol(cluster_tree)])]
  
  ## Obtain best modularity
  jss_modularity <- apply(cluster_tree[,-1], 2, function(x){
    modularity(
      x, jss, resolution = 1
    )
  })
  return_cluster <- which.max(jss_modularity)
  
  
  ## Return data
  results <- list()
  results$clusters <- cluster_tree[,return_cluster + 1]
  names(results$clusters) <- cluster_tree$id
  results$clusterTree <- cluster_tree
  results$clusterPlot <- phylo_tree
  results$JSS <- jss
  
  ## Set class
  class(results) <- "infoCluster"
  
  return(results)
  
}
