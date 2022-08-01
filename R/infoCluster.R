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
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#' 
#' @param plot.cluster Boolean.
#' Should plot of optimal and hierarchical clusters be output?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not plot
#'
#' @examples
#'# Obtain data
#' sim.dynEGA <- sim.dynEGA # bypasses CRAN checks
#' 
#' \donttest{# Dynamic EGA individual and population structure
#' dyn.ega1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 21, use.derivatives = 1, 
#'   ncores = 2, corr = "pearson"
#' )
#' 
#' # Perform information-theoretic clustering
#' clust1 <- infoCluster(
#'   dynEGA.object = dyn.ega1,
#'   ncores = 2, # Two cores for CRAN checks
#'   plot.cluster = FALSE # No plot for CRAN checks
#' )}
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
#' @importFrom stats hclust as.dist cutree
#' @importFrom utils globalVariables
#' 
#' @export
# Information Theoretic Clustering for dynEGA
# Updated 29.07.2022
infoCluster <- function(
    dynEGA.object,
    ncores,
    plot.cluster = TRUE
)
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
  
  # # Obtain memberships
  # wcs <- lapply(dynEGA.ind, function(x){
  #   x$wc
  # })
  # 
  # # Message user
  # message("Computing Variation of Information...\n", appendLF = FALSE)
  # 
  # # Obtain lists
  # vi_lists <- lapply(
  #   X = 2:length(wcs),
  #   FUN = function(i){
  #     
  #     # Compute VI values
  #     vi_values <- lapply(1:(i-1), function(j){
  #       
  #       # Try
  #       vi_value <- try(
  #         vi(
  #           wc1 = wcs[[i]],
  #           wc2 = wcs[[j]]
  #         ),
  #         silent = TRUE
  #       )
  #       
  #       # Check if value is OK
  #       if(!is(vi_value, "try-error")){
  #         return(vi_value)
  #       }else{
  #         return(NA)
  #       }
  #       
  #     })
  #     
  #     # Return
  #     return(vi_values)
  #     
  #   }
  # )
  # 
  # # Organize data
  # vi_i <- lapply(vi_lists, unlist)
  # 
  # # Initialize JSD matrix
  # vi_matrix <- matrix(
  #   0,
  #   nrow = length(networks),
  #   ncol = length(networks)
  # )
  # 
  # # Loop through
  # for(i in 1:length(vi_i)){
  #   vi_matrix[i+1,1:(length(vi_i[[i]]))] <- vi_i[[i]]
  # }
  # 
  # # Make symmetric
  # vi_sym <- vi_matrix + t(vi_matrix)
  # 
  # # Add names
  # colnames(vi_sym) <- names(networks)
  # row.names(vi_sym) <- names(networks)
  # 
  # # Variation of Information distance
  # vidist <- vi_sym
  # 
  # # Message user
  # message("done.")
  
  # Message user
  message("Computing Jensen-Shannon Distance...\n", appendLF = FALSE)
  
  # Set cores (if missing)
  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }
  
  # Make cluster
  cl <- parallel::makeCluster(ncores)
  
  # Export
  # parallel::clusterExport(
  #   cl = cl,
  #   varlist = c(
  #     # "rescaled_laplacian",
  #     # "vn_entropy",
  #     "jsd",
  #     "networks"
  #   ),
  #   envir = environment()
  # )
  
  # Obtain lists
  jsd_lists <- pbapply::pblapply(
    cl = cl,
    X = 2:length(networks),
    FUN = function(i){
      
      # Compute JSD values
      jsd_values <- lapply(1:(i-1), function(j){
        
        # Try
        jsd_value <- try(
          jsd(
            network1 = networks[[i]],
            network2 = networks[[j]],
            method = "spectral"
          ),
          silent = TRUE
        )
        
        # Check if value is OK
        if(!is(jsd_value, "try-error")){
          return(jsd_value)
        }else{
          return(NA)
        }
        
      })
      
      # Return
      return(jsd_values)
      
    }
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  # Organize data
  jsd_i <- lapply(jsd_lists, unlist)
  
  # Initialize JSD matrix
  jsd_matrix <- matrix(
    0,
    nrow = length(networks),
    ncol = length(networks)
  )
  
  # Loop through
  for(i in 1:length(jsd_i)){
    jsd_matrix[i+1,1:(length(jsd_i[[i]]))] <- jsd_i[[i]]
  }
  
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
  
  # Combine VI and JSD
  # infodist <- vidist + jsdist
  # diag(infodist) <- 0
  
  # Perform hierarchical clustering
  hier_clust <- hclust(
    d = as.dist(jsdist),
    method = "complete"
  )
  
  # Jensen-Shannon Similarity
  jss <- 1 - jsdist

  # Make diagonal of Jensen-Shannon Similarity = 0
  diag(jss) <- 0
  
  # Compute modularity matrix
  Q_matrix <- modularity_matrix(
    A = jss,
    resolution = 1
  )
  
  # Maximize modularity
  Qs <- unlist(
    lapply(
      X = 1:ncol(jss),
      FUN = function(i){
        quick_modularity(
          communities = cutree(hier_clust, i),
          A = jss,
          Q_matrix = Q_matrix
        )
      }
    )
  )
  
  # Obtain clusters
  clusters <- cutree(hier_clust, which.max(Qs))
  
  # Check if single cluster
  if(length(unique(na.omit(clusters))) == 1){

    # Message user
    message("One cluster detected. Testing random networks to determine the validity...")

    # Generate random networks
    random_networks <- lapply(networks, function(x){

      # Obtain random network
      random_network <- randnet(
        nodes = ncol(x),
        edges = sum(ifelse(x != 0, 1, 0)) / 2
      )

      # Return random network
      return(random_network)

    })

    # Message user
    message("Computing Jensen-Shannon Distance for Random Networks...\n", appendLF = FALSE)

    # Make cluster
    cl <- parallel::makeCluster(ncores)

    # Export
    # parallel::clusterExport(
    #   cl = cl,
    #   varlist = c(
    #     # "rescaled_laplacian",
    #     # "vn_entropy",
    #     "jsd",
    #     "random_networks"
    #   ),
    #   envir = environment()
    # )

    # Obtain lists
    jsd_random_lists <- pbapply::pblapply(
      cl = cl,
      X = 2:length(random_networks),
      FUN = function(i){

        # Compute JSD values
        jsd_values <- lapply(1:(i-1), function(j){

          # Try
          jsd_value <- try(
            jsd(
              network1 = random_networks[[i]],
              network2 = random_networks[[j]],
              method = "spectral"
            ),
            silent = TRUE
          )

          # Check if value is OK
          if(!is(jsd_value, "try-error")){
            return(jsd_value)
          }else{
            return(NA)
          }

        })

        # Return
        return(jsd_values)

      }
    )

    # Stop cluster
    parallel::stopCluster(cl)

    # Organize data
    jsd_random_i <- lapply(jsd_random_lists, unlist)

    # Initialize JSD matrix
    jsd_random_matrix <- matrix(
      0,
      nrow = length(random_networks),
      ncol = length(random_networks)
    )

    # Loop through
    for(i in 1:length(jsd_random_i)){
      jsd_random_matrix[i+1,1:(length(jsd_random_i[[i]]))] <- jsd_random_i[[i]]
    }

    # Make symmetric
    jsd_random_sym <- jsd_random_matrix + t(jsd_random_matrix)

    # Add names
    colnames(jsd_random_sym) <- names(random_networks)
    row.names(jsd_random_sym) <- names(random_networks)

    # Make jsdist
    jsdist_random <- jsd_random_sym

    # Make diagonal NA
    diag(jsdist_random) <- NA

    # Remove all NAs
    rm_cols <- apply(jsdist_random, 2, function(x){all(is.na(x))})

    # Remove missing data points
    jsdist_random <- jsdist_random[!rm_cols, !rm_cols]

    # Make diagonal 0 again
    diag(jsdist_random) <- 0

    # Compare to empirical
    comparison <- t.test(
      jsdist[lower.tri(jsdist)],
      jsdist_random[lower.tri(jsdist_random)],
      paired = TRUE,
      var.equal = FALSE
    )
    
    # Obtain sign of statistic
    comparison_sign <- sign(comparison$statistic)
    
    # Compute adaptive alpha
    adaptive_p <- adapt.a(
      test = "paired",
      n = length(jsdist[lower.tri(jsdist)]),
      alpha = .001,
      power = 0.80,
      efxize = "large"
    )
    
    # Check for empirical JSD > random JSD OR
    # non-significant t-test
    if(
      comparison_sign == 1 |
      comparison$p.value > adaptive_p$adapt.a
    ){
      
      # Set clusters to all individuals
      clusters <- 1:ncol(jsdist)
      names(clusters) <- colnames(jsdist)
      
      # Let user know
      message(
        "Empirical Jensen-Shannon distance was no different or greater than random Jensen-Shannon Distance: No clusters exist."
      )
      
    }

  }
  
  # No plot if no clusters
  if(all(clusters != 1:ncol(jsdist))){
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
    
    # Split dendrogram into upper grey section and lower coloured section
    cut <- max(clusters)
    height <- unique(cluster_data$segments$y)[order(unique(cluster_data$segments$y), decreasing = TRUE)]
    cut.height <- mean(c(height[cut], height[cut-1]))
    cluster_data$segments$line <- ifelse(cluster_data$segments$y == cluster_data$segments$yend &
                                           cluster_data$segments$y > cut.height, 1, 2)
    cluster_data$segments$line <- ifelse(cluster_data$segments$yend  > cut.height, 1, cluster_data$segments$line)
    
    # Number the clusters
    cluster_data$segments$cluster <- c(-1, diff(cluster_data$segments$line))
    change <- which(cluster_data$segments$cluster == 1)
    for (i in 1:cut) cluster_data$segments$cluster[change[i]] = i + 1
    cluster_data$segments$cluster <-  ifelse(cluster_data$segments$line == 1, 1, 
                                             ifelse(cluster_data$segments$cluster == 0, NA, cluster_data$segments$cluster))
    
    
    # Replace NA values in cluster
    for(i in seq_along(cluster_data$segments$cluster)){
      
      if(is.na(cluster_data$segments$cluster[i])){
        cluster_data$segments$cluster[i] <- cluster_data$segments$cluster[i-1]
      }
      
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
    
    # Check if plot should be plotted
    if(isTRUE(plot.cluster)){
      suppressWarnings(
        plot(cluster_plot)
      )
    }
  }else{
    cluster_plot <- NULL
  }
  
  ## Return data
  results <- list(
    clusters = clusters,
    modularity = Qs[which.max(Qs)],
    clusterTree = hier_cluster,
    clusterPlot = cluster_plot,
    JSD = jsdist
  )

  ## Set class
  class(results) <- "infoCluster"
  
  return(results)
  
}
