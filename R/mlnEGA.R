# Multilayer behavioral and brain networks
#' Multilayer Behavioral and Brain Networks
#' 
#' @description Computes a multilayered network between behavioral and brain networks
#' following Brooks et al. (2020; see references)
#' 
#' @param behav.data Matrix or data frame.
#' Behavioral data where rows are participants and
#' columns are variables
#' 
#' @param behav.time Boolean.
#' Is the behavioral data time series?
#' Defaults to \code{FALSE}
#' 
#' @param brain.data Array or list.
#' Entries should be correlation matrices that correspond
#' to the participants in \code{behav.data}
#' 
#' @param brain.time Boolean.
#' Is the brain data time series?
#' Defaults to \code{FALSE}
#' 
#' @param alpha Numeric.
#' Significance level for interlayer links.
#' Defaults to \code{.05}
#' 
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#' 
#' @param EGA.scores Boolean.
#' Should \code{\link[EGAnet]{EGA}} be used to determine
#' dimensions and compute network scores (\code{\link[EGAnet]{net.scores}})?
#' Defaults to \code{TRUE}.
#' Uses \code{\link[EGAnet]{dynEGA}} for time series data.
#' Set to \code{FALSE} if \code{behav.data} is already scored
#' 
#' @param EGA.args List.
#' A list of additional arguments for computing EGA.
#' See \code{\link[EGAnet]{EGA}} for cross-sectional data and
#' \code{\link[EGAnet]{dynEGA}} for time series data
#' 
#' @param plot.mln Boolean.
#' If \code{TRUE}, returns a plot of the multilayer network.
#' Defaults to \code{TRUE}
#' 
#' @param plot.args List.
#' A list of additional arguments for the network plot:
#' 
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.4.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.7.}
#' 
#' }
#' 
#' @return Returns a list containing:
#' 
#' \item{behavioral.scores}{Either scores input as \code{data} or network scores estimated from
#' \code{\link[EGAnet]{EGA}}}
#' 
#' \item{brain.networks}{A list containing:
#' 
#' \itemize{
#' 
#' \item{\code{individuals}}
#' {A list containing the correlation matrix of each individual participant}
#' 
#' \item{\code{group}}
#' {A matrix of the group correlation matrix (mean of each ROI--ROI connection across participants)}
#' 
#' \item{\code{group.threshold}}
#' {The threshold group brain matrix. The threshold is the critical correlation value
#' based on the \code{alpha} that was input}
#' 
#' }
#' 
#' }
#' 
#' \item{behavioral.network}{A network estimation using the \code{behavioral.scores} and \code{\link[qgraph]{EBICglasso}}}
#' 
#' \item{interlayer.links}{Correlations between the \code{behavioral.scores} and node \code{\link[NetworkToolbox]{strength}}
#' of each node in the \code{brain.networks$individuals} network}
#' 
#' \item{multilayer.network}{The adjacency matrix of the multilayer network that combines the \code{behavioral.network},
#' \code{interlayer.links}, and \code{brain.networks$group.threshold} networks}
#' 
#' @examples
#' 
#' @references 
#' Brooks, D., Hulst, H. E., de Bruin, L., Glas, G., Geurts, J. J. G., & Douw, L. (2020).
#' The multilayer network approach in the study of personality neuroscience.
#' \emph{Brain Sciences}, \emph{10}, 915.
#' doi:\href{https://doi.org/10.3390/brainsci10120915}{10.3390/brainsci10120915}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
# Multilayer behavioral and brain networks ----
# Updated 07.12.2020
mlnEGA <- function (behav.data, behav.time = FALSE,
                    brain.data, brain.time = FALSE,
                    alpha = .05, ncores,
                    EGA.scores = TRUE, EGA.args = list(),
                    plot.mln = TRUE, plot.args = list())
{
  # Stop on behav.time = TRUE
  # NEED TO ADD DYNEGA CAPABILITY
  if(behav.time){
    stop("Functionality for time series behavioral data has not been implemented yet. Coming soon. . .")
  }
  
  # Missing arguments
  if(missing(ncores)){
    ncores <- parallel::detectCores() / 2
  }
  
  
  # Check for whether brain data is in an array
  if(is.array(brain.data)){
    
    # Initialize brain list
    brain.list <- vector("list", length = dim(brain.data)[3])
    
    # Input brains into list
    for(i in 1:dim(brain.data)[3]){
      brain.list[[i]] <- brain.data[,,i]
    }
    
  }else{brain.list <- brain.data}
  
  # Check for time series brain data
  if(brain.time){## Compute distance correlation
    brain.list <- NetworkToolbox::dCor.parallel(brain.list, ncores = ncores)
  }
  
  # Check if number of rows of behavioral data matches
  # the number of participants in brain data
  if(nrow(behav.data) != length(brain.list)){
    stop("Number of rows in 'behav.data' does not match the length of 'brain.data'")
  }
  
  # Compute group average brain network
  brain.network <- apply(simplify2array(brain.list), 1:2, mean)
  
  # Critical r value
  critical.r <- function (n, a)
  {
    df <- n - 2
    critical.t <- qt( a/2, df, lower.tail = FALSE )
    cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
    return(cvr)
  }
  
  # Critical value
  crit.val <- critical.r(nrow(behav.data), alpha)
  
  # Thresholded group network
  brain.threshold <- brain.network
  brain.threshold[abs(brain.threshold) < crit.val] <- 0
  
  # Compute strength for brain data nodes
  str.list <- lapply(brain.list, NetworkToolbox::strength)
  
  # Convert to matrix
  str.mat <- t(simplify2array(str.list))
  
  # Use EGA scores?
  if(EGA.scores){
    
    if(behav.time){
      
      # ADD DYNEGA HERE
      
    }else{
      
      # Compute EGA of behavioral data
      ## Get default arguments
      ega.defaults <- formals(EGA)[-which(names(formals(EGA)) == "...")]
      ## Set defaults
      ega.defaults$data <- behav.data
      ega.defaults$model <- "glasso"
      ega.defaults$algorithm <- "walktrap"
      ega.defaults$plot.EGA <- FALSE # Set to FALSE to speed up computation
      ega.defaults$plot.type <- "GGally"
      ega.defaults$verbose <- FALSE # Ignore messages for now
      
      ## Input arguments
      if(any(names(EGA.args) == names(ega.defaults))){
        
        ### Target arguments
        target.args <- names(ega.defaults[which(names(EGA.args) == names(ega.defaults))])
        
        ega.defaults[target.args] <- EGA.args[target.args]
      }
      
      ## EGA
      ega <- do.call(EGA, ega.defaults)
      
    }
    
    # Compute network scores
    network.scores <- net.scores(data = behav.data, A = ega, impute = "none")$std.scores
    
  }else{network.scores <- data}
  
  ## Get network of behavioral scores
  ega.defaults$data <- network.scores
  behav.network <- suppressWarnings(suppressMessages(do.call(EGA, ega.defaults)$network))
  
  # Compute interlayer links
  interlinks <- cor(str.mat, network.scores)
  
  # Set links that are not greater than critical value to zero
  interlinks[abs(interlinks) < crit.val] <- 0
  
  # Construct multilayer network matrix
  multilayer <- matrix(0, nrow = ncol(behav.network) + ncol(brain.network),
                       ncol = ncol(behav.network) + ncol(brain.network))
  colnames(multilayer) <- c(colnames(behav.network), colnames(brain.network))
  rownames(multilayer) <- colnames(multilayer)
  
  ## Input values
  ### Input interlayer
  multilayer[colnames(interlinks), rownames(interlinks)] <- interlinks
  multilayer <- t(multilayer)
  multilayer[colnames(interlinks), rownames(interlinks)] <- interlinks
  ### Input behavioral network
  multilayer[colnames(behav.network), colnames(behav.network)] <- behav.network
  ### Input brain network
  multilayer[colnames(brain.network), colnames(brain.network)] <- brain.threshold
  
  # Return results
  res <- list()
  res$behavioral.scores <- network.scores
  res$brain.networks$individuals <- brain.list
  res$brain.networks$group <- brain.network
  res$brain.networks$group.threshold <- brain.threshold
  res$behavioral.network <- behav.network
  res$interlayer.links <- interlinks
  res$multilayer.network <- multilayer
  
  # Plot
  if(plot.mln){
    
    ## Check for input plot arguments
    if(length(plot.args) == 0){
      plot.args <-list(vsize = 6, alpha = 0.4, label.size = 5, edge.alpha = 0.7)}
    else{
      plot.args <- plot.args
      plots.arg1 <- list(vsize = 6, label.size = 5, alpha = 0.4, edge.alpha = 0.7)
      plot.args.use <- plot.args
      
      if(any(names(plots.arg1) %in% names(plot.args.use))){
        
        plot.replace.args <- plots.arg1[na.omit(match(names(plot.args.use), names(plots.arg1)))]
        
        plot.args <- c(plot.args.use,plots.arg1[names(plots.arg1) %in% names(plot.args.use)==FALSE])}
    }
    
    # Get network
    network1 <- network::network(res$multilayer.network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    # Set attributes
    network::set.vertex.attribute(network1, attrname= "Layers", value = c(rep("Behavioral", ncol(res$behavioral.scores)),
                                                                               rep("Neural", ncol(res$brain.networks$group.threshold))))
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(res$multilayer.network))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(scales::rescale(as.vector(res$multilayer.network),
                                                         to = c(.001, 1.75)),
                                         nrow = nrow(res$multilayer.network),
                                         ncol = ncol(res$multilayer.network)))
    
    # Layout "Spring"
    graph1 <- NetworkToolbox::convert2igraph(res$multilayer.network)
    edge.list <- igraph::as_edgelist(graph1)
    layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                               weights =
                                                                 abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                               vcount = ncol(res$multilayer.network))
    
    
    set.seed(1234)
    res$mln.plot <- GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                   color = "Layers", edge.color = c("color"),
                   alpha = plot.args$alpha, size = plot.args$vsize,
                   edge.alpha = plot.args$edge.alpha,
                   label.size = plot.args$label.size,
                   mode =  layout.spring,
                   label = colnames(res$multilayer.network)) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)
      )
    
    plot(res$mln.plot)
  }
  
  class(res) <- "MLN"
  
  return(res)
  
}
#----