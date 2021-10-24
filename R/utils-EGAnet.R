#%%%%%%%%%%%%%%%%%
# DEVELOPMENT ----
#%%%%%%%%%%%%%%%%%

#' @noRd
#' 
# Polytomous IRT parameters
# Updated 30.06.2021
poly.irt <- function(loadings, data)
{
  # Check for standardized loadings
  if(any(class(loadings) == "NetLoads")){
    loadings <- loadings$std
  }
  
  # Convert loadings to matrix
  loadings <- as.matrix(loadings)
  
  # Unique variance
  s <- sqrt(1 - rowSums(loadings^2))
  
  # Estimated discrimination parameters
  est_a <- sweep(loadings, 1, s, "/")
  
  # Estimated threshold parameters
  thresholds <- lavaan::lavCor(
    data,
    ordered = colnames(data),
    output = "thresholds"
  )
  
  # Separate thresholds
  threshs <- list()
  
  for(i in colnames(data)){
    threshs[[i]] <- thresholds[grep(i, names(thresholds))]
  }
  
  # Estimate location parameters
  ## All variables have the same number of categories
  if(all(unlist(lapply(threshs, length)))){
    
    threshs <- t(simplify2array(threshs))
    
    est_d <- sweep(-threshs, 1, s, "/")
  }else{## Some variables have different numbers of categories
    
    est_d <- lapply(threshs, function(x, s){
      -x/s
    }, s = s)
    
  }
  
  # Results list
  results <- list()
  results$discrimination <- est_a
  results$location <- est_d
  
  return(results)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# network.descriptives ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# From WGCNA version 1.70-3
#' @noRd
#' @importFrom graphics hist
#' @importFrom stats lm
#'
#Scale-free fit index
#Updated 12.05.2021
scaleFreeFitIndex=function(k,nBreaks=10, removeFirst = FALSE)
{
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), 
                length = nBreaks + 1)
  hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk= as.numeric(log10(p.dk + 1e-09))
  lm1 = try(lm(log.p.dk ~ log.dk));
  if (inherits(lm1, "try-error")) browser();
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  datout=data.frame(Rsquared.SFT=summary(lm1)$r.squared,
                    slope.SFT=summary(lm1)$coefficients[2, 1], 
                    truncatedExponentialAdjRsquared= summary(lm2)$adj.r.squared)
  datout
}

#%%%%%%%%%%%%%%%
# net.loads ----
#%%%%%%%%%%%%%%%

#' Communicating Nodes
#' @description Computes the between-community strength for each node in the network
#'
#' @param A An adjacency matrix of network data
#'
#' @param comm Can be a vector of community assignments or community detection algorithms
#' (\code{"walktrap"} or \code{"louvain"}) can be used to determine the number of factors.
#' Defaults to \code{"walktrap"}.
#' Set to \code{"louvain"} for \code{\link[NetworkToolbox]{louvain}} community detection
#'
#' @param cent Centrality measure to be used.
#' Defaults to \code{"strength"}.
#'
#' @param absolute Should network use absolute weights?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for signed weights
#'
#' @param metric Whether the metric should be compute for across all of the communities
#' (a single value) or for each community (a value for each community).
#' Defaults to \code{"across"}.
#' Set to \code{"each"} for values for each community
#'
#' @param diagonal Sets the diagonal values of the \code{A} input.
#' Defaults to \code{0}
#'
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#'
#' @return A vector containing the between-community strength value for each node
#'
#' @references
#' Blanken, T. F., Deserno, M. K., Dalege, J., Borsboom, D., Blanken, P., Kerkhof, G. A., & Cramer, A. O. (2018).
#' The role of stabilizing and communicating symptoms given overlapping communities in psychopathology networks.
#' \emph{Scientific Reports}, \emph{8}, 5854.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
#Communicating
#Updated 18.03.2020
comcat <- function (A, comm = c("walktrap","louvain"),
                    cent = c("strength","degree"),
                    absolute = TRUE,
                    metric = c("across","each"),
                    diagonal = 0, ...)
{
  # MISSING ARGUMENTS
  
  if(missing(comm))
  {comm <- "walktrap"
  }else{comm <- comm}
  
  if(missing(cent))
  {cent <- "strength"
  }else{cent <- match.arg(cent)}
  
  
  if(missing(metric))
  {metric <- "each"
  }else{metric <- match.arg(metric)}
  
  if(missing(diagonal))
  {diagonal <- 0
  }else{diagonal <- diagonal}
  
  # MAIN FUNCTION
  
  # Set diagonal
  diag(A) <- diagonal
  
  # Change edges to absolute
  if(absolute)
  {A <- abs(A)}
  
  # Convert to communities
  if(any(eval(formals(NetworkToolbox::stable)$comm) %in% comm))
  {
    facts <- switch(comm,
                    walktrap = igraph::cluster_walktrap(NetworkToolbox::convert2igraph(A), ...)$membership,
                    louvain = igraph::cluster_louvain(NetworkToolbox::convert2igraph(A), ...)$membership
    )
  }else{facts <- comm}
  
  # Convert facts to characters
  facts <- paste(facts)
  
  # Check for names of nodes
  if(is.null(colnames(A)))
  {colnames(A) <- paste("V", 1:ncol(A), sep = "")}
  
  names(facts) <- colnames(A)
  
  # Unique communities
  uniq <- unique(facts)
  
  # Initialize list
  fact <- list()
  
  # Check for unidimensionality
  if(length(na.omit(uniq)) != 1)
  {
    if(metric=="across") # All communities
    {
      
      # Loop over all communities not in target community
      for(i in 1:length(uniq))
      {
        # Connections outside of target community
        Ah <- A[which(facts != uniq[i]), which(facts == uniq[i])]
        
        # Centrality
        com <- switch(cent,
                      degree = colSums(NetworkToolbox::binarize(Ah)),
                      strength = colSums(Ah)
        )
        
        # Input into list
        fact[[i]] <- com
      }
      
      # Unlist to vector
      commn <- unlist(fact)
      
      # Reorder to be consist with labels
      commat <- commn[names(facts)]
      
      # Label vector
      names(commat) <- colnames(A)
      
    }else if(metric=="each") # Individual communities
    {
      # Initialize item list
      item <- list()
      
      # Initialize matrix
      commat <- matrix(NA, nrow = nrow(A), ncol = length(uniq))
      
      # Add column names
      colnames(commat) <- paste(uniq)
      
      # Loop through each node
      for(i in 1:ncol(A))
      {
        # Connections for node 'i'
        Ah <- A[,i]
        
        # Communities outside of target community
        uniq.no <- uniq[which(uniq!=facts[i])]
        
        # Loop through each community
        for(j in 1:length(uniq.no))
        {
          # Edges in target outside community
          Aha <- Ah[which(facts==uniq.no[j])]
          
          # Centrality
          com <- switch(cent,
                        degree = sum(ifelse(Aha!=0,1,0)),
                        strength = sum(Aha)
          )
          
          # Input into matrix
          commat[i,paste(uniq.no[j])] <- com
        }
      }
      
      # Order and label matrix
      colnames(commat) <- uniq
      row.names(commat) <- colnames(A)
    }
  }else{
    # Initialize matrix
    commat <- as.matrix(rep(0,ncol(A)), ncol = 1)
    # Label matrix
    colnames(commat) <- paste(unique(comm))
    row.names(commat) <- colnames(A)
  }
  
  return(commat)
}

#' Stabilizing Nodes
#' @description Computes the within-community centrality for each node in the network
#'
#' @param A An adjacency matrix of network data
#'
#' @param comm Can be a vector of community assignments or community detection algorithms
#' (\code{"walktrap"} or \code{"louvain"}) can be used to determine the number of factors.
#' Defaults to \code{"walktrap"}.
#' Set to \code{"louvain"} for \code{\link[NetworkToolbox]{louvain}} community detection
#'
#' @param cent Centrality measure to be used.
#' Defaults to \code{"strength"}.
#'
#' @param absolute Should network use absolute weights?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for signed weights
#'
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#'
#' @param diagonal Sets the diagonal values of the \code{A} input.
#' Defaults to \code{0}
#'
#' @return A matrix containing the within-community centrality value for each node
#'
#' @references
#' Blanken, T. F., Deserno, M. K., Dalege, J., Borsboom, D., Blanken, P., Kerkhof, G. A., & Cramer, A. O. (2018).
#' The role of stabilizing and communicating symptoms given overlapping communities in psychopathology networks.
#' \emph{Scientific Reports}, \emph{8}, 5854.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Stabilizing
#Updated 18.03.2020
stable <- function (A, comm = c("walktrap","louvain"),
                    cent = c("betweenness","rspbc","closeness",
                             "strength","degree","hybrid"),
                    absolute = TRUE, diagonal = 0, ...)
{
  # MISSING ARGUMENTS
  
  if(missing(comm))
  {comm <- "walktrap"
  }else{comm <- comm}
  
  if(missing(diagonal))
  {diagonal <- 0
  }else{diagonal <- diagonal}
  
  if(missing(cent))
  {cent <- "strength"
  }else{cent <- match.arg(cent)}
  
  # MAIN FUNCTION
  
  # Set diagonal
  diag(A) <- diagonal
  
  # Make weights absolute
  if(absolute)
  {A <- abs(A)}
  
  # Convert to communities
  if(any(eval(formals(NetworkToolbox::stable)$comm) %in% comm))
  {
    facts <- switch(comm,
                    walktrap = igraph::cluster_walktrap(NetworkToolbox::convert2igraph(A), ...)$membership,
                    louvain = igraph::cluster_louvain(NetworkToolbox::convert2igraph(A), ...)$membership
    )
  }else{facts <- comm}
  
  # Convert facts to characters
  facts <- paste(facts)
  
  # Check for names of nodes
  if(is.null(colnames(A)))
  {colnames(A) <- paste("V", 1:ncol(A), sep = "")}
  
  names(facts) <- colnames(A)
  
  # Unique communities
  uniq <- unique(facts)
  
  # Initialize community list
  fact <- list()
  
  # Loop through computing within-community centrality
  for(i in 1:length(uniq))
  {
    # Nodes only in community 'i'
    Ah <- A[which(facts == uniq[i]), which(facts == uniq[i])]
    
    # Check for matrix size
    if(length(Ah) != 1)
    {
      # Centrality measure
      stab <- switch(cent,
                     betweenness = NetworkToolbox::betweenness(Ah),
                     rspbc = NetworkToolbox::rspbc(Ah),
                     closeness = NetworkToolbox::closeness(Ah),
                     strength = colSums(Ah),
                     degree = colSums(NetworkToolbox::binarize(Ah))
      )
    }else{
      # Input zero
      stab <- 0
      # Change name
      names(stab) <- names(which(facts == uniq[i]))
    }
    
    # Input into list
    fact[[i]] <- stab
  }
  
  # Unlist for vector
  stabil <- unlist(fact)
  
  # Reorder to be consist with labels
  stabil <- stabil[names(facts)]
  
  # Check for missing values (change to 0)
  stabil <- ifelse(is.na(stabil),0,stabil)
  
  return(stabil)
}

#' Add signs
#' @description Adds signs to network loading matrix
#'
#' @param comm.str Matrix.
#' Unstandardized network loading matrix (node strength split between dimensions)
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' An adjacency matrix of network data
#'
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param dims Numeric.
#' Vector corresponding to dimensions
#'
#' @param pos.manifold Boolean.
#' Should a positive manifold be applied (i.e., should
#' all dimensions be positively correlated)?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for a positive manifold
#'
#' @return A matrix with signs appropriately added
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Add signs
# Updated 18.03.2020
add.signs <- function(comm.str, A, wc, dims, pos.manifold)
{
  # Signs within dimension
  for(i in 1:length(dims))
  {
    # Target dimension
    target <- which(wc==dims[i])
    
    # Initialize signs
    signs <- numeric(length(target))
    
    # Target matrix
    target.mat <- A[target,target]
    
    # Sign matrix
    sign.mat <- sign(target.mat)
    
    for(j in 1:nrow(sign.mat))
    {
      # Save original sign.mat
      orig.sign <- sign.mat
      
      # Check for max sums (current max)
      curr.max <- sum(colSums(sign.mat,na.rm=TRUE),na.rm=TRUE)
      
      # New sign mat
      sign.mat[j,] <- -sign.mat[j,]
      
      # Check for new max sums (new max)
      new.max <- sum(colSums(sign.mat,na.rm=TRUE),na.rm=TRUE)
      
      if(new.max <= curr.max)
      {
        sign.mat <- orig.sign
        signs[j] <- 1
      }else{
        signs[j] <- -1
      }
    }
    
    comm.str[which(wc==dims[i]),i] <- comm.str[which(wc==dims[i]),i] * ifelse(signs==0,1,signs)
    A[,which(wc==dims[i])] <- sweep(A[,which(wc==dims[i])],2,ifelse(signs==0,1,signs),`*`)
  }
  
  # Signs between dimensions
  for(i in 1:length(dims))
    for(j in 1:length(dims))
    {
      if(i!=j)
      {
        # Target dimension
        target1 <- which(wc==dims[i])
        target2 <- which(wc==dims[j])
        
        # Initialize signs
        signs <- numeric(length(target1))
        
        # Target matrix
        target.mat <- A[target1,target2]
        
        # Sign matrix
        sign.mat <- sign(target.mat)
        
        for(k in 1:nrow(sign.mat))
        {
          # Save original sign.mat
          orig.sign <- sign.mat
          
          # Check for max sums (current max)
          curr.max <- sum(colSums(sign.mat,na.rm=TRUE),na.rm=TRUE)
          
          # New sign mat
          sign.mat[k,] <- -sign.mat[k,]
          
          # Check for new max sums (new max)
          new.max <- sum(colSums(sign.mat,na.rm=TRUE),na.rm=TRUE)
          
          if(new.max <= curr.max)
          {
            sign.mat <- orig.sign
            signs[k] <- 1
          }else{
            signs[k] <- -1
          }
        }
        
        comm.str[which(wc==dims[i]),j] <- comm.str[which(wc==dims[i]),j] * ifelse(signs==0,1,signs)
      }
    }
  
  
  # Flip dimensions (if necessary)
  if(!pos.manifold)
  {
    for(i in 1:length(dims))
    {
      wc.sign <- sign(sum(comm.str[which(wc==dims[i]),i]))
      
      if(wc.sign != 1)
      {comm.str[which(wc==dims[i]),] <- -comm.str[which(wc==dims[i]),]}
    }
  }
  
  res <- list()
  res$comm.str <- comm.str
  res$A <- A
  
  return(res)
}

#' Unstandardized network loading matrix
#' @description Computes the unstandardized network loading matrix
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' An adjacency matrix of network data
#'
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param metric Character.
#' Argument to be passed onto \code{comcat}
#'
#' @param absolute Boolean.
#' Should absolute values be computed?
#' Defaults to \code{TRUE}
#'
#' @param diagonal Numeric.
#' Values to be input on the diagonal of the network
#'
#' @return An unstandardized network loading matrix
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Unstandardized Network Loadings
# Updated 02.07.2020
mat.func <- function(A, wc, metric = "each", absolute, diagonal)
{
  # Compute within-community node strength
  stab <- stable(A = A, comm = wc, absolute = absolute, diagonal = diagonal)
  
  if(length(unique(wc)) == 1)
  {comm.str <- round(stab,3)
  }else{
    
    # Compute between-community node strength
    comc <- comcat(A = A, comm = wc, metric = metric, absolute = absolute, diagonal = diagonal)
    
    # Ensure variable ordering is correct
    stab <- stab[row.names(comc)]
    
    # Combine between- and within-community node strength
    indices <- which(is.na(comc), arr.ind = TRUE)
    indices <- indices[names(stab),]
    
    for(i in 1:nrow(indices)){
      comc[indices[i,1], indices[i,2]] <- stab[i]
    }
    
    # Round to 3 decimal places
    comm.str <- round(comc, 3)
    
  }
  
  return(comm.str)
}

#' @noRd
# Function to order loadings largest to smallest
# within their respective factors
descend.ord <- function(loads, wc){
  # Initialize ordering vector
  ord.names <- vector("character")
  
  # Loop through dimensions
  for(i in colnames(loads)){
    ord <- order(loads[names(which(wc == i)),i], decreasing = TRUE)
    ord.names <- c(ord.names, names(which(wc == i))[ord])
  }
  
  # Reorder
  reord <- loads[ord.names,]
  
  # Check for matrix
  if(!is.matrix(reord)){
    reord <- as.matrix(reord)
  }
  
  # Make sure names
  row.names(reord) <- ord.names
  colnames(reord) <- colnames(loads)
  
  return(reord)
}

#%%%%%%%%%%%%%%
# PLOTTING ----
#%%%%%%%%%%%%%%

#' @noRd
# Sub-routine for compare.EGA.plots
# Updated 05.06.2021
compare.EGA <- function(ega.object1, ega.object2)
{
  # Plots
  plot1 <- ega.object1
  plot2 <- ega.object2
  
  # Reorder node coordinates for plot2
  plot2$data <- plot2$data[row.names(plot1$data),]
  
  # Reorder edge coordinates for plot2
  for(i in 1:nrow(plot2$layers[[1]]$data)){
    
    plot2$layers[[1]]$data$X1[i] <- which(plot2$layers[[1]]$data$X1[i] == plot2$data$x)
    plot2$layers[[1]]$data$X2[i] <- which(plot2$layers[[1]]$data$X2[i] == plot2$data$x)
    plot2$layers[[1]]$data$Y1[i] <- which(plot2$layers[[1]]$data$Y1[i] == plot2$data$y)
    plot2$layers[[1]]$data$Y2[i] <- which(plot2$layers[[1]]$data$Y2[i] == plot2$data$y)
    
  }
  
  # Reassign edge coordinates based on plot1
  plot2$layers[[1]]$data$X1 <- plot1$data$x[plot2$layers[[1]]$data$X1] # X1
  plot2$layers[[1]]$data$X2 <- plot1$data$x[plot2$layers[[1]]$data$X2] # X2
  plot2$layers[[1]]$data$Y1 <- plot1$data$y[plot2$layers[[1]]$data$Y1] # Y1
  plot2$layers[[1]]$data$Y2 <- plot1$data$y[plot2$layers[[1]]$data$Y2] # Y2
  
  # Assign coordinates of plot1 to plot2
  plot2$data$x <- plot1$data$x
  plot2$data$y <- plot1$data$y
  
  # Make plot list
  plots.net <- list()
  plots.net[[1]] <- plot1
  plots.net[[2]] <- plot2
  
  # Return plot list
  return(plots.net)
  
}

#' @noRd
# Defaults for GGally plotting
# For plots and methods
# Updated 28.07.2021
GGally.args <- function(plot.args)
{
  default.args <- formals(GGally::ggnet2)
  ega.default.args <- list(node.size = 12, edge.size = 8,
                           alpha = 0.5, label.size = 5,
                           edge.alpha = 0.5, layout.exp = 0.2)
  default.args[names(ega.default.args)]  <- ega.default.args
  default.args <- default.args[-length(default.args)]
  
  
  if("node.alpha" %in% names(plot.args)){
    plot.args$alpha <- plot.args$node.alpha
    plot.args$node.alpha <- NULL
  }
  
  if("vsize" %in% names(plot.args)){
    plot.args$node.size <- plot.args$vsize
    plot.args$vsize <- NULL
  }
  
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
    plot.args$legend.names <- NULL
  }
  
  if(!"color.palette" %in% names(plot.args)){
    default.args$color.palette <- "polychrome"
  }
  
  if("color.palette" %in% names(plot.args)){
    
    if(tolower(plot.args$color.palette) == "greyscale" | tolower(plot.args$color.palette) == "grayscale" | tolower(plot.args$color.palette) == "colorblind"){
      plot.args$edge.color <- c("#293132", "grey25")
      plot.args$edge.lty <- c("solid", "dashed")
    }
    
  }
  
  if(!"edge.color" %in% names(plot.args)){
    plot.args$edge.color <- c("darkgreen", "red")
  }else if(length(plot.args$edge.color) != 2){
    stop("Two colors needed for 'edge.color'")
  }
  
  if(!"edge.lty" %in% names(plot.args)){
    plot.args$edge.lty <- c("solid", "solid")
  }else if(length(plot.args$edge.lty) != 2){
    stop("Two line types needed for 'edge.lty'")
  }
  
  
  if(any(names(plot.args) %in% names(default.args))){
    target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
    default.args[names(target.args)] <- target.args
  }
  
  plot.args <- default.args
  
  return(plot.args)
}

#' @importFrom graphics text
#' @noRd
# Color sorting for EGA palettes
# For EGA_color_palette
# Updated 17.12.2020
color.sort <- function (wc)
{
  unlist(lapply(sort(wc, na.last = TRUE), function(x, uniq){
    
    if(is.na(x)){
      NA
    }else{
      which(x == uniq)
    }
  }, uniq = sort(unique(wc), na.last = TRUE)))
}

#' @noRd
# Rescale edges for GGally
# For plots
# Updated 17.01.2021
rescale.edges <- function (network, size)
{
  # Set rescaling
  scaling <- seq(0, 1, .000001) * size
  names(scaling) <- seq(0, 1, .000001)
  
  # Vectorize edges
  edges <- round(as.vector(as.matrix(network)), 5)
  
  # Get absolute edges
  abs.edges <- abs(edges)
  
  # Get edge signs
  signs.edges <- sign(edges)
  
  # Rescale edges
  rescaled.edges <- unname(scaling[as.character(abs.edges)])
  
  return(rescaled.edges)
}

#' @noRd
# Compare plots fix
# Updated 09.10.2021
compare.plot.fix.EGA <- function(object.list,  plot.type = c("GGally","qgraph"),
                                 plot.args = list()){
  #### MISSING ARGUMENTS HANDLING
  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  ## Check for input plot arguments
  if(plot.type == "GGally"){
    if("legend.names" %in% names(plot.args)){
      legend.names <- plot.args$legend.names
    }
    plot.args <- GGally.args(plot.args)
    color.palette <- plot.args$color.palette
  }
  
  ## Original plot arguments
  original.plot.args <- plot.args
  
  ## Initialize plot list
  ega.plots <- list()
  
  # Loop through object list
  for(i in 1:length(object.list)){
    
    if(class(object.list[[i]]) == "EGA"){
      x <- object.list[[i]]
    }else if(class(object.list[[i]]) == "bootEGA"){
      x <- list(
        network = object.list[[i]]$typicalGraph$graph,
        wc = object.list[[i]]$typicalGraph$wc
      )
    }else if(class(object.list[[i]]) == "dynEGA"){
      x <- object.list[[i]]$dynEGA
    }
    
    ### Plot ###
    if(plot.type == "qgraph"){
      ega.plot <- qgraph::qgraph(x$network, layout = "spring", vsize = plot.args$vsize, groups = as.factor(x$wc))
    }else if(plot.type == "GGally"){
      
      # Insignificant values (keeps ggnet2 from erroring out)
      x$network <- ifelse(abs(as.matrix(x$network)) <= .00001, 0, as.matrix(x$network))
      
      if(exists("legend.names")){
        for(l in 1:length(unique(legend.names))){
          x$wc[x$wc == l] <- legend.names[l]
        }
      }
      
      # Reorder network and communities
      if(i == 1){
        x$network <- x$network[order(x$wc), order(x$wc)]
        x$wc <- x$wc[order(x$wc)]
        fix.order.wc <- names(x$wc)
      }else{
        x$network <- x$network[fix.order.wc, fix.order.wc]
        x$wc <- x$wc[fix.order.wc]
      }
      
      # weighted  network
      network1 <- network::network(x$network,
                                   ignore.eval = FALSE,
                                   names.eval = "weights",
                                   directed = FALSE)
      
      network::set.vertex.attribute(network1, attrname= "Communities", value = x$wc)
      network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
      network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
      network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
      network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$network))
      network::set.edge.value(network1,attrname="ScaledWeights",
                              value=matrix(rescale.edges(x$network, plot.args$edge.size),
                                           nrow = nrow(x$network),
                                           ncol = ncol(x$network)))
      
      # Layout "Spring"
      graph1 <- NetworkToolbox::convert2igraph(x$network)
      edge.list <- igraph::as_edgelist(graph1)
      layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                 weights =
                                                                   abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                                 vcount = ncol(x$network))
      
      
      set.seed(1234)
      plot.args$net <- network1
      plot.args$node.color <- "Communities"
      plot.args$node.alpha <- plot.args$alpha
      plot.args$node.shape <- plot.args$shape
      node.size <- plot.args$node.size
      plot.args$node.size <- 0
      plot.args$color.palette <- NULL
      plot.args$palette <- NULL
      plot.args$edge.color <- "color"
      plot.args$edge.lty <- "line"
      plot.args$edge.size <- "ScaledWeights"
      
      lower <- abs(x$network[lower.tri(x$network)])
      non.zero <- sqrt(lower[lower != 0])
      
      plot.args$edge.alpha <- non.zero
      plot.args$mode <- layout.spring
      plot.args$label <- colnames(x$network)
      plot.args$node.label <- rep("", ncol(x$network))
      if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
      if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
      
      palette <- color_palette_EGA(color.palette, as.numeric(factor(x$wc)))
      palette <- ifelse(is.na(palette), "white", palette)
      
      ega.plot <- suppressWarnings(
        suppressMessages(
          do.call(GGally::ggnet2, plot.args) + 
            ggplot2::theme(legend.title = ggplot2::element_blank())
        )
      )
      
      name <- colnames(x$network)
      
      name.split <- lapply(name, function(x){
        unlist(strsplit(x, split = " "))
      })
      
      name <- unlist(
        lapply(name.split, function(x){
          
          len <- length(x)
          
          if(len > 1){
            
            add.line <- round(len / 2)
            
            paste(
              paste(x[1:add.line], collapse = " "),
              paste(x[(add.line+1):length(x)], collapse = " "),
              sep = "\n"
            )
            
          }else{x}
          
        })
      )
      
      # Border color
      if(all(color.palette == "grayscale" |
             color.palette == "greyscale" |
             color.palette == "colorblind")){
        border.color <- ifelse(palette == "white", "white", "black")
      }else{border.color <- palette}
      
      # Custom nodes: transparent insides and dark borders
      ega.plot <- ega.plot + 
        ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                            color = border.color,
                            shape = 1, stroke = 1.5, alpha = .8) +
        ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
                            color = palette,
                            shape = 19, alpha = plot.args$alpha) +
        ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
        ggplot2::guides(
          color = ggplot2::guide_legend(override.aes = list(
            color = unique(palette),
            size = node.size,
            alpha = plot.args$alpha,
            stroke = 1.5
          ))
        )
    }
    
    set.seed(NULL)
    
    # Return to object.list
    ega.plots[[i]] <- ega.plot
    
    # Reset plot.args
    plot.args <- original.plot.args
    
  }
  
  return(ega.plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MULTI-FUNCTION SUB-ROUTINES ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{\link[qgraph]{EBICglasso}} from \code{\link{qgraph}} 1.4.4
#'
#' This function uses the \code{\link[glasso]{glasso}} package
#' (Friedman, Hastie and Tibshirani, 2011) to compute a
#' sparse gaussian graphical model with the graphical lasso
#' (Friedman, Hastie & Tibshirani, 2008).
#' The tuning parameter is chosen using the Extended Bayesian Information criterium
#' (EBIC) described by Foygel & Drton (2010).
#'
#' @param data Data matrix
#'
#' @param n Number of participants
#'
#' @param gamma EBIC tuning parameter. 0.5 is generally a good choice.
#' Setting to zero will cause regular BIC to be used.
#'
#' @param penalize.diagonal Should the diagonal be penalized?
#'
#' @param nlambda Number of lambda values to test.
#'
#' @param lambda.min.ratio Ratio of lowest lambda value compared to maximal lambda
#'
#' @param returnAllResults   If \code{TRUE} this function does not
#' return a network but the results of the entire glasso path.
#'
#' @param penalizeMatrix Optional logical matrix to indicate which elements are penalized
#'
#' @param countDiagonal     Should diagonal be counted in EBIC computation?
#' Defaults to \code{FALSE}. Set to \code{TRUE} to mimic qgraph < 1.3 behavior (not recommended!).
#'
#' @param refit Logical, should the optimal graph be refitted without LASSO regularization?
#' Defaults to \code{FALSE}.
#'
#' @param ... Arguments sent to \code{\link[glasso]{glasso}}
#'
#' @details The glasso is run for 100 values of the tuning parameter logarithmically
#' spaced between the maximal value of the tuning parameter at which all edges are zero,
#' lambda_max, and lambda_max/100. For each of these graphs the EBIC is computed and
#' the graph with the best EBIC is selected. The partial correlation matrix
#' is computed using \code{\link[qgraph]{wi2net}} and returned.
#'
#' @return A partial correlation matrix
#'
#' @references
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#' Sparse inverse covariance estimation with the graphical lasso.
#' \emph{Biostatistics}, \emph{9}, 432-441.
#' doi: \href{https://doi.org/10.1093/biostatistics/kxm045}{10.1093/biostatistics/kxm045}
#'
#' #glasso package
#' Jerome Friedman, Trevor Hastie and Rob Tibshirani (2011).
#' glasso: Graphical lasso-estimation of Gaussian graphical models.
#' R package version 1.7.
#' \url{https://CRAN.R-project.org/package=glasso}
#'
#' Foygel, R., & Drton, M. (2010).
#' Extended Bayesian information criteria for Gaussian graphical models.
#' In Advances in neural information processing systems (pp. 604-612).
#' \url{https://papers.nips.cc/paper/4087-extended-bayesian-information-criteria-for-gaussian-graphical-models}
#'
#' #psych package
#' Revelle, W. (2014) psych: Procedures for Personality and Psychological Research,
#' Northwestern University, Evanston, Illinois, USA.
#' R package version 1.4.4.
#' \url{https://CRAN.R-project.org/package=psych}
#'
#' #Matrix package
#' Douglas Bates and Martin Maechler (2014).
#' Matrix: Sparse and Dense Matrix Classes and Methods.
#' R package version 1.1-3.
#' \url{https://CRAN.R-project.org/package=Matrix}
#'
#' @author Sacha Epskamp <mail@sachaepskamp.com>
#'
#' @examples
#' ### Using wmt2 dataset from EGAnet ###
#' data(wmt2)
#'
#' \dontrun{
#' # Compute correlations:
#' CorMat <- cor_auto(wmt2[,7:24])
#'
#' # Compute graph with tuning = 0 (BIC):
#' BICgraph <- EBICglasso.qgraph(CorMat, nrow(wmt2), 0)
#'
#' # Compute graph with tuning = 0.5 (EBIC)
#' EBICgraph <- EBICglasso.qgraph(CorMat, nrow(wmt2), 0.5)
#'
#' }
#'
#' @noRd
#'
# Computes optimal glasso network based on EBIC:
# Updated 24.03.2020
EBICglasso.qgraph <- function(
  data, # Sample covariance matrix
  n = NULL,
  gamma = 0.5,
  penalize.diagonal = FALSE, # Penalize diagonal?
  nlambda = 100,
  lambda.min.ratio = 0.01,
  returnAllResults = FALSE, # If true, returns a list
  penalizeMatrix, # Optional logical matrix to indicate which elements are penalized
  countDiagonal = FALSE, # Set to TRUE to get old qgraph behavior: conting diagonal elements as parameters in EBIC computation. This is not correct, but is included to replicate older analyses
  refit = FALSE, # If TRUE, network structure is taken and non-penalized version is computed.
  ... # glasso arguments
) {
  
  # Codes originally implemented by Sacha Epskamp in his qgraph package version 1.4.4.
  # Selects optimal lamba based on EBIC for given covariance matrix.
  # EBIC is computed as in Foygel, R., & Drton, M. (2010, November). Extended Bayesian Information Criteria for Gaussian Graphical Models. In NIPS (pp. 604-612). Chicago
  
  # Simply computes the Gaussian log likelihood given sample covariance and estimate of precision:
  
  # Original:
  # logGaus <- function(S,K,n)
  # {
  #   SK = S %*% K
  #   tr = function(A) sum(diag(A))
  #   n/2 * (log(det(K)) - tr(SK))
  # }
  
  ## According to huge???
  logGaus <- function(S,K,n)
  {
    KS = K %*% S
    tr = function(A) sum(diag(A))
    return(n/2 * (log(det(K)) - tr(KS))  )
  }
  
  # Computes the EBIC:
  EBIC <- function(S,K,n,gamma = 0.5,E,countDiagonal=FALSE)
  {
    #   browser()
    L <- logGaus(S, K, n)
    if (missing(E)){
      E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
    }
    p <- nrow(K)
    
    # return EBIC:
    -2 * L + E * log(n) + 4 * E * gamma * log(p)
  }
  
  # Computes partial correlation matrix given precision matrix:
  wi2net <- function(x)
  {
    x <- -stats::cov2cor(x)
    diag(x) <- 0
    x <- Matrix::forceSymmetric(x)
    return(x)
  }
  
  if(is.null(n))
  {
    if(nrow(data)!=ncol(data))
    {n <- nrow(data)
    }else{stop("Number of participants 'n' need to be specified")}
  }
  
  # Compute correlations matrix
  if(nrow(data)!=ncol(data))
  {S <- qgraph::cor_auto(data)
  }else{
    S <- data
  }
  
  # Compute lambda sequence (code taken from huge package):
  lambda.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  lambda.min = lambda.min.ratio*lambda.max
  lambda = exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  
  # Run glasso path:
  if (missing(penalizeMatrix)){
    glas_path <- glasso::glassopath(S, lambda, trace = 0, penalize.diagonal=penalize.diagonal, ...)
  }else{
    glas_path <- list(
      w = array(0, c(ncol(S), ncol(S), length(lambda))),
      wi = array(0, c(ncol(S), ncol(S), length(lambda))),
      rholist = lambda
    )
    
    for (i in 1:nlambda){
      res <- glasso::glasso(S, penalizeMatrix * lambda[i], trace = 0, penalize.diagonal=penalize.diagonal, ...)
      glas_path$w[,,i] <- res$w
      glas_path$wi[,,i] <- res$wi
    }
  }
  
  
  # Compute EBICs:
  #     EBICs <- apply(glas_path$wi,3,function(C){
  #       EBIC(S, C, n, gamma)
  #     })
  
  lik <- sapply(seq_along(lambda),function(i){
    logGaus(S, glas_path$wi[,,i], n)
  })
  
  EBICs <- sapply(seq_along(lambda),function(i){
    EBIC(S, glas_path$wi[,,i], n, gamma, countDiagonal=countDiagonal)
  })
  
  # Smallest EBIC:
  opt <- which.min(EBICs)
  
  # Check if rho is smallest:
  #if (opt == 1){
  #  warning("Network with lowest lambda selected as best network. Try setting 'lambda.min.ratio' lower.")
  #}
  
  # Return network:
  net <- as.matrix(Matrix::forceSymmetric(wi2net(glas_path$wi[,,opt])))
  colnames(net) <- rownames(net) <- colnames(S)
  
  # Check empty network:
  if (all(net == 0)){
    message("An empty network was selected to be the best fitting network. Possibly set 'lambda.min.ratio' higher to search more sparse networks. You can also change the 'gamma' parameter to improve sensitivity (at the cost of specificity).")
  }
  
  # Refit network:
  # Refit:
  if (refit){
    message("Refitting network without LASSO regularization")
    glassoRes <- suppressWarnings(glasso::glasso(S, 0, zero = which(net == 0 & upper.tri(net), arr.ind=TRUE), trace = 0, penalize.diagonal=penalize.diagonal, ...))
    net <- as.matrix(Matrix::forceSymmetric(wi2net(glassoRes$wi)))
    colnames(net) <- rownames(net) <- colnames(S)
    optwi <- glassoRes$wi
  } else {
    optwi <- glas_path$wi[,,opt]
  }
  
  # Return
  if (returnAllResults){
    return(list(
      results = glas_path,
      ebic = EBICs,
      loglik = lik,
      optnet = net,
      lambda = lambda,
      optwi = optwi
    ))
  } else return(net)
}

#' Proportion table
#'
#' @param boot.mat Matrix.
#' A matrix of bootstrapped memberships
#'
#' @return Matrix of proportions based on dimensions
#'
#' @noRd
#'
# Proportion Table
# Updated 26.02.2021
proportion.table <- function (boot.mat)
{
  # Get maximum number of dimensions
  max.dim <- max(boot.mat, na.rm = TRUE)
  
  # Set up table
  tab <- matrix(0, nrow = nrow(boot.mat), ncol = max.dim)
  colnames(tab) <- 1:max.dim # assign column names
  
  # Check if there are row names
  if(!is.null(row.names(boot.mat))){
    row.names(tab) <- row.names(boot.mat)
  }
  
  # Loop through maximum dimensions
  for(i in 1:max.dim)
  {tab[,i] <- apply(boot.mat, 1, function(x){mean(x == i, na.rm = TRUE)})}
  
  return(tab)
}


#' @noRd
# Normalize DNN weights
# Updated 24.03.2021
min.max <- function(vec)
{
  exp.min <- exp(0) / exp(1)
  exp.max <- exp(1) / exp(0)
  
  return((vec - exp.min) / (exp.max - exp.min))
}

#' @noRd
# Custom range min-max
# Updated 26.05.2021
custom.min.max <- function(vec, ran)
{
  a <- ran[1]
  b <- ran[2]
  
  return((b - a) * ((vec - min(vec)) / (max(vec) - min (vec))) + a)
}

#' @noRd
#' New NA function
#' Updated 24.10.2021
is.NA <- function(x){
  
  # Regular check
  reg.na <- is.na(x)
  
  # Character check
  char.na <- x == "NA"
  
  # Get any NA
  return(as.logical(reg.na + char.na))
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%
# DATA GENERATION ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Expand correlation matrix (unidimensional EGA)
# For EGA
# Updated 04.02.2021
expand.corr <- function(corr)
{
  # Number of variables
  nvar <- 4
  
  # Create additional correlations
  add.corr <- matrix(.50, nrow = nvar, ncol = nvar)
  diag(add.corr) <- 1
  
  # Initialize new correlation matrix
  new.corr <- matrix(0, nrow = nvar + nrow(corr), ncol = nvar + ncol(corr))
  
  # Input correlations into new correlation matrix
  new.corr[1:nvar, 1:nvar] <- add.corr
  new.corr[(nvar+1):nrow(new.corr), (nvar+1):ncol(new.corr)] <- corr
  
  # Add names
  colnames(new.corr) <- c(paste("SIM",1:nvar, sep = ""), colnames(corr))
  row.names(new.corr) <- colnames(new.corr)
  
  return(new.corr)
}

#' A sub-routine to simulate data for \code{\link[EGAnet]{EGA}}
#'
#' @param data Data for number of cases
#'
#' @param nvar Number of variables
#'
#' @param nfact Number of factors
#'
#' @param load Magnitude of loadings
#'
#' @return Simulated data
#'
#' @importFrom utils data
#'
#' @author Loli Dolores Nieto and Luis Eduardo Garrido
#'
#' @noRd
#'
# Simulate data function
# Updated 04.12.2020
sim.func <- function(data, nvar, nfact, load, n = 500)
{
  
  # Check for unidimensional structure
  ## Set up data simulation
  if(missing(data)){
    data <- matrix(nrow = n, ncol = 0)
  }else{
    n <- nrow(data)
  }
  corf <- 0
  J <- nvar*nfact
  sdcross = 0
  
  ## GENERATE SAMPLE DATA MATRIX
  check.eig <- TRUE
  check.com <- TRUE
  
  while(check.eig == TRUE|check.com == TRUE)
  {
    SATF = matrix(0, J, nfact)
    
    for(j in 1:nfact)
    {
      SATF[(j*nvar-nvar+1):(j*nvar),j]<-runif(nvar, load-.10, load+.10)
      
      if(nfact>1)
      {
        CROSS.L <- apply(as.matrix(SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)]), 2, function(x) rnorm((nvar-2), 0, sdcross))
        
        SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)] <- CROSS.L
      }
    }
    
    #SATF # Population factor loading matrix with cross-loadings and marker items
    
    FCOR      = matrix(corf, nfact, nfact); diag(FCOR)<-1 ## Factor correlation matrix
    R         = SATF%*%FCOR%*%t(SATF)                          ## Rr
    check.com = any(diag(R) > .90)                                  ## Check communalities values
    diag(R)   = 1                                                                    ## Insert ones in the diagonal of Rr
    #R                                                                                       ## Rp
    check.eig = any(eigen(R)$values <= 0)                      ## Check eigenvalues
  }
  
  U = chol(R)                                                                       ## Cholesky decomposition of Rp
  Z = MASS::mvrnorm(n, mu = rep(0, J), Sigma = diag(J))                                  ## Obtain sample matrix of continuous variables
  X = Z%*%U
  colnames(X) <- paste0("X", 1:ncol(X))
  
  data.sim <- cbind(X, data)
  
  return(data.sim)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## DATA CATEGORIZATION FUNCTION
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Script: Function to categorize continuous data
## % DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## First, likely values for the skewness of the variables must be specified
## according to a set of values ranging from -2 to 2 in increments of 0.5.
## This is done with the skew.values argument. Second, for the selected
## number of response categories to simulate, an object with the thresholds
## to produce different skewness is generated. Then, the skewness for each
## variable is randomly assigned and each variable is categorized according
## to the thresholds previously defined. The thresholds used to categorize
## the data set are those defined by Garrido, Abad, & Ponsoda (2011, 2013).
## %% ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## data:        data set of continuous variables
## ncat:        number of response categories
## skew.values: a vector with several of the following values: -2, -1.5,
##              -1, -0.5, 0, 0.5, 1, 1.5, or 2. It can also be a positive
##              integer with a single value.
## %% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## data.cat:  categorized data
## SIM.SKEW:  simulated skewness for each variable
## TRUE.SKEW: real generated skewness for each variable
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Updated 04.12.2020
# Authors: Maria Dolores Nieto and Luis Eduardo Garrido
categorize<-function(data, ncat, skew.values){
  
  data.cat<-matrix(0, nrow = nrow(data), ncol = ncol(data))
  SIM.SKEW  = rep(0, ncol(data))
  TRUE.SKEW = rep(0, ncol(data))
  
  if(ncat==2){# SKEW = -2     -1.5     -1     -0.5  0   0.5     1     1.5     2
    thres<-matrix(c(-1.0518,-0.8416,-0.5936,-0.3088,0,0.3088,0.5936,0.8416,1.0518), ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values
    
    for(j in 1:ncol(data)){
      
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][data[,j] < thres[1,COLUMN]]<-1
      data.cat[,j][data[,j] >=thres[1,COLUMN]]<-2
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==3){
    thres.2  <-c(0.8518,  1.3754)  ## SKEW = 2
    thres.1.5<-c(0.6131,  1.1969)  ## SKEW = 1.5
    thres.1  <-c(0.3195,  0.9921)  ## SKEW = 1.0
    thres.0.5<-c(-0.0236, 0.7256)  ## SKEW = 0.5
    thres.0  <-c(-1.0000, 1.0000)  ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values
    
    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][ data[,j]>=thres[2,COLUMN]]<-3
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==4){
    thres.2  <-c(0.7515,  1.1341, 1.5980) ## SKEW = 2
    thres.1.5<-c(0.4945,  0.9299, 1.4359) ## SKEW = 1.5
    thres.1  <-c(0.1678,  0.6873, 1.2513) ## SKEW = 1.0
    thres.0.5<-c(-0.2057, 0.3706, 0.9809) ## SKEW = 0.5
    thres.0  <-c(-1.5000, 0.0000, 1.5000) ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values
    
    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][(data[,j]>=thres[2,COLUMN])&(data[,j]<thres[3,COLUMN])]<-3
      data.cat[,j][ data[,j]>=thres[3,COLUMN]]<-4
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==5){
    thres.2  <-c(0.6792,  1.0043,  1.3441, 1.7703) ## SKEW = 2
    thres.1.5<-c(0.4071,  0.7827,  1.1596, 1.6186) ## SKEW = 1.5
    thres.1  <-c(0.0502,  0.5117,  0.9432, 1.4462) ## SKEW = 1.0
    thres.0.5<-c(-0.3414, 0.1642,  0.6257, 1.1645) ## SKEW = 0.5
    thres.0  <-c(-1.8000, -0.6000, 0.6000, 1.8000) ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values
    
    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][(data[,j]>=thres[2,COLUMN])&(data[,j]<thres[3,COLUMN])]<-3
      data.cat[,j][(data[,j]>=thres[3,COLUMN])&(data[,j]<thres[4,COLUMN])]<-4
      data.cat[,j][ data[,j]>=thres[4,COLUMN]]<-5
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  
  RESULTS<-list(data.cat  = data.cat,
                SIM.SKEW  = SIM.SKEW,
                TRUE.SKEW = TRUE.SKEW)
  
  return(RESULTS)
}

#%%%%%%%%%
# LCT ----
#%%%%%%%%%

#' @noRd
# Dynamic organization
# Updated 26.05.2021
dyn.org <- function(data, gen)
{
  for(i in 1:ncol(data)){
    
    gen[,i] <- custom.min.max(
      sort(gen[,i])[rank(data[,i], ties.method = "first")],
      range(data[,i])
    )
    
  }
  
  return(gen)
}

#' A sub-routine to compute the deep learning neural network model
#' weights for \code{\link[EGAnet]{LCT}}
#'
#' @param loads Matrix of loadings
#'
#' @param weights Weights from specific model (see \code{\link[EGAnet]{dnn.weights}})
#'
#' @return A prediction or probability of the specified model
#'
#' @noRd
#'
# DNN weights function
# Updated 09.08.2020
dnn.model.weights <- function (loads, weights)
{
  wb <- seq(1, length(weights), 3)
  
  for(i in wb[-length(wb)])
  {
    if(i == 1)
    {
      input <- as.vector(t(weights[[i]]) %*% as.matrix(loads)) + weights[[(i+1)]]
      input <- ifelse(input > 0, input, input * weights[[i+2]])
      layer <- input
    }else{
      layer <- as.vector(t(weights[[i]]) %*% as.matrix(layer)) + weights[[(i+1)]]
      layer <- ifelse(layer > 0, layer, layer * weights[[i+2]])
    }
  }
  
  # Output
  output <- as.vector(t(weights[[wb[length(wb)]]]) %*% as.matrix(layer)) + weights[[length(weights)]]
  
  # Sigmoid activation function
  sigmoid <- function(x)
  {exp(x) / (exp(x) + 1)}
  
  # Prediction
  prediction <- sigmoid(output)
  
  return(prediction)
}

#' A sub-routine to predict the model for \code{\link[EGAnet]{LCT}}
#'
#' @param loads Matrix of loadings
#'
#' @return The model prediction
#'
#' @importFrom utils data
#'
#' @noRd
#'
# DNN prediction function
# Updated 30.03.2021
dnn.predict <- function (loads)
{
  # Load deep learning neural network weights
  dnn.weights <- get(data("dnn.weights", envir = environment()))
  
  # Compute ratios
  ## Small ratio
  small.ratio <- exp(loads[1]) / exp(loads[6])
  ## Moderate ratio
  moderate.ratio <- exp(loads[2]) / exp(loads[7])
  ## Large ratio
  large.ratio <- exp(loads[3]) / exp(loads[8])
  ## Dominant ratio
  dominant.ratio <- exp(loads[4]) / exp(loads[9])
  ## Cross ratio
  cross.ratio <- exp(loads[5]) / exp(loads[10])
  
  small.ratio <- min.max(small.ratio)
  moderate.ratio <- min.max(moderate.ratio)
  large.ratio <- min.max(large.ratio)
  dominant.ratio <- min.max(dominant.ratio)
  cross.ratio <- min.max(cross.ratio)
  
  # Factor versus network model
  f_n <- vector("numeric", length = 3)
  
  # Check for low correlation factor versus network model
  f_n[1] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$lf_n_weights)
  
  # Check for high correlation with variables greater than factors versus network model
  f_n[2] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$hvgf_n_weights)
  
  # Check for high correlation factor versus network model
  f_n[3] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$hvlf_n_weights)
  
  # Check for factor model
  ifelse(any(f_n >= .50), return(1), return(2))
  
}

#%%%%%%%%%%%%%
# bootEGA ----
#%%%%%%%%%%%%%

#' A sub-routine to generate typical network structure following \code{\link{EGAnet}{EGA}} approach
#'
#' @noRd
#'
# Typical network (bootEGA) function
# Updated 12.03.2020
typicalStructure.network <- function (A, corr, model, model.args, n = NULL, uni.method,
                                      algorithm, algorithm.args)
{
  
  # Convert to igraph
  graph <- suppressWarnings(NetworkToolbox::convert2igraph(abs(A)))
  
  # Check for unconnected nodes
  if(igraph::vcount(graph)!=ncol(A)){
    
    warning("Estimated network contains unconnected nodes:\n",
            paste(names(which(NetworkToolbox::degree(A)==0)), collapse = ", "))
    
    unconnected <- which(NetworkToolbox::degree(A)==0)
    
  }
  
  # Algorithm Arguments
  ## Check for algorithm
  if(!is.function(algorithm)){
    
    if(algorithm == "walktrap"){
      algorithm.formals <- formals(igraph::cluster_walktrap)
    }else if(algorithm == "louvain"){
      algorithm.formals <- formals(igraph::cluster_louvain)
    }
    
  }else{algorithm.formals <- formals(algorithm)}
  
  ## Check for input algorithm arguments
  if(length(algorithm.args) != 0){
    
    ### Check for matching arguments
    if(any(names(algorithm.args) %in% names(algorithm.formals))){
      
      algorithm.replace.args <- algorithm.args[na.omit(match(names(algorithm.formals), names(algorithm.args)))]
      
      algorithm.formals[names(algorithm.replace.args)] <- algorithm.replace.args
    }
    
  }
  
  ## Remove ellipses
  if("..." %in% names(algorithm.formals)){
    algorithm.formals[which(names(algorithm.formals) == "...")] <- NULL
  }
  
  ## Remove weights from igraph functions' arguments
  if("weights" %in% names(algorithm.formals)){
    algorithm.formals[which(names(algorithm.formals) == "weights")] <- NULL
  }
  
  # Multidimensional result
  ## Run community detection algorithm
  algorithm.formals$graph <- graph
  
  if(!is.function(algorithm)){
    
    multi.wc <- switch(algorithm,
                       "walktrap" = do.call(igraph::cluster_walktrap, as.list(algorithm.formals))$membership,
                       "louvain" = do.call(igraph::cluster_louvain, as.list(algorithm.formals))$membership
    )
    
  }else{multi.wc <- do.call(what = algorithm, args = as.list(algorithm.formals))$membership}
  
  # Get new data
  if(model == "glasso"){
    
    # Obtain inverse of network
    g <- -A
    diag(g) <- 1
    
  }else if(model == "TMFG"){
    
    # Generate data
    g.data <- MASS::mvrnorm(n, mu = rep(0, ncol(A)), Sigma = as.matrix(Matrix::nearPD(A, corr = TRUE, keepDiag = TRUE)$mat))
    g <- -suppressMessages(NetworkToolbox::LoGo(g.data, normal = TRUE, partial = TRUE))
    diag(g) <- 1
    
  }
  
  # New data
  data <- MASS::mvrnorm(n, mu = rep(0, ncol(g)), Sigma = corpcor::pseudoinverse(g))
  
  # Check for unidimensional structure
  if(uni.method == "expand"){
    
    # Check for Spinglass algorithm
    if(is.function(algorithm)){
      
      # spins argument is used to identify Spinglass algorithm
      if("spins" %in% methods::formalArgs(algorithm)){
        
        # Simulate data from unidimensional factor model
        sim.data <- sim.func(data = data, nvar = 4, nfact = 1, load = .70)
        
        ## Compute correlation matrix
        cor.data <- switch(corr,
                           "cor_auto" = qgraph::cor_auto(sim.data),
                           "pearson" = cor(sim.data, use = "pairwise.complete.obs"),
                           "spearman" = cor(sim.data, method = "spearman", use = "pairwise.complete.obs")
        )
        
      }else{
        
        ## Compute correlation matrix
        cor.data <- switch(corr,
                           "cor_auto" = qgraph::cor_auto(data),
                           "pearson" = cor(data, use = "pairwise.complete.obs"),
                           "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
        )
        
        ## Expand correlation matrix
        cor.data <- expand.corr(cor.data)
        
      }
      
    }else{# Do regular adjustment
      
      ## Compute correlation matrix
      cor.data <- switch(corr,
                         "cor_auto" = qgraph::cor_auto(data),
                         "pearson" = cor(data, use = "pairwise.complete.obs"),
                         "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
      )
      
      ## Expand correlation matrix
      cor.data <- expand.corr(cor.data)
      
    }
    
    # Unidimensional result
    uni.res <- EGA.estimate(data = cor.data, n = n,
                            model = model, model.args = model.args,
                            algorithm = algorithm, algorithm.args = algorithm.args)
    
    ## Remove simulated data for multidimensional result
    cor.data <- cor.data[-c(1:4),-c(1:4)]
    
    if(uni.res$n.dim <= 2){
      wc <- uni.res$wc[-c(1:4)]
    }else{
      wc <- multi.wc
    }
    
  }else if(uni.method == "LE"){
    
    ## Compute correlation matrix
    cor.data <- switch(corr,
                       "cor_auto" = qgraph::cor_auto(data),
                       "pearson" = cor(data, use = "pairwise.complete.obs"),
                       "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
    )
    
    # Leading eigenvalue approach for one and two dimensions
    wc <- igraph::cluster_leading_eigen(NetworkToolbox::convert2igraph(abs(cor.data)))$membership
    names(wc) <- colnames(cor.data)
    n.dim <- length(na.omit(unique(wc)))
    
    
    # Set up results
    if(n.dim != 1){
      wc <- multi.wc
    }
    
  }
  
  # Obtain community memberships
  init.wc <- as.vector(matrix(NA, nrow = 1, ncol = ncol(A)))
  names(init.wc) <- colnames(A)
  init.wc[1:length(wc)] <- wc
  wc <- init.wc
  
  # Replace unconnected nodes with NA communities
  if(exists("unconnected")){
    wc[unconnected] <- NA
  }
  
  names(wc) <- colnames(A)
  
  return(wc)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dynEGA and mctest.ergoInfo ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Dynamic EGA used in the mctest.ergoInfo function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Dynamic EGA used in the mctest.ergoInfo function
#' @description Dynamic EGA used in the mctest.ergoInfo function. DynEGA estimates dynamic factors in multivariate time series (i.e. longitudinal data, panel data, intensive longitudinal data) at multiple
#' time scales, in different levels of analysis: individuals (intraindividual structure) and population (structure of the population).
#' Exploratory graph analysis is applied in the derivatives estimated using generalized local linear approximation (\code{\link[EGAnet]{glla}}). Instead of estimating factors by modeling how variables are covarying, as in traditional
#' EGA, dynEGA is a dynamic model that estimates the factor structure by modeling how variables are changing together.
#' GLLA is a filtering method for estimating derivatives from data that uses time delay embedding and a variant of Savitzky-Golay filtering to accomplish the task.
#'
#' @param data A dataframe with the variables to be used in the analysis. The dataframe should be in a long format (i.e. observations for the same individual (for example, individual 1) are placed in order, from time 1 to time t, followed by the observations from individual 2, also ordered from time 1 to time t.)
#'
#' @param n.embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"n.embed = 5"} will use five consecutive observations to estimate a single derivative.
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings in the \code{\link[EGAnet]{Embed}} function. A tau of one uses adjacent observations.
#' Default is \code{"tau = 1"}.
#'
#' @param delta Integer.
#' The time between successive observations in the time series.
#' Default is \code{"delta = 1"}.
#'
#' @param id Numeric.
#' Number of the column identifying each individual.
#'
#'
#' @param use.derivatives Integer.
#' The order of the derivative to be used in the EGA procedure. Default to 1.
#'
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use. Defaults to \code{glasso}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#' \dontrun{
#' \donttest{# Population structure:
#' dyn.ega1 <- dynEGA.ind.pop(data = sim.dynEGA, n.embed = 5, tau = 1,
#' delta = 1, id = 21, use.derivatives = 1, model = "glasso", ncores = 2,
#' cor = "pearson")
#' }
#' }
#'
#' @importFrom stats cor rnorm runif na.omit
#' @export

dynEGA.ind.pop <- function(data, n.embed, tau = 1, delta = 1,
                           id = NULL,
                           use.derivatives = 1,
                           model = c("glasso", "TMFG"), model.args = list(),
                           algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                           corr = c("cor_auto", "pearson", "spearman"),
                           ncores, ...){
  
  # Get additional arguments
  add.args <- list(...)
  
  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      paste(
        "The 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")",
        sep = ""
      )
    )
    
    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }
  
  
  # MISSING ARGUMENTS HANDLING
  if(missing(id))
  {stop("The 'id' argument is missing! \n The number of the column identifying each individual must be provided!")
  }else{id <- id}
  
  if(missing(corr))
  {corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}
  
  # Setting the order:
  
  order = 2
  
  ### Spliting by ID:
  
  #number of cases
  cases <- unique(data[,id])
  
  #initialize data list
  datalist <- vector("list", length = length(cases))
  datalist <- split(data[,-c(id)],data[,id])
  
  ### Estimating the derivatives using GLLA:
  
  #let user know derivatives estimation has started
  message("\nComputing derivatives using GLLA...\n", appendLF = FALSE)
  
  
  #initialize derivatives list
  derivlist <- list()
  
  #Parallel processing
  cl <- parallel::makeCluster(ncores)
  
  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "derivlist", "cases"),
                          envir=environment())
  
  # GLLA Estimation:
  glla.multi <- function(data, n.embed = n.embed, tau = tau, delta = delta, order = order){
    order.deriv <- paste0("Ord",seq(from = 0, to = order))
    data.est <- vector("list")
    for(i in 1:ncol(data)){
      data.est[[i]] <- as.data.frame(EGAnet::glla(data[,i], n.embed = n.embed, tau = tau, delta = delta, order = order))
      data.est[[i]] <- as.data.frame(data.est[[i]])
    }
    data.est2 <- vector("list")
    for(i in 0:order+1){
      data.est2[[i]] <- sapply(data.est, "[[", i)
    }
    
    data.estimates <- data.frame(Reduce(cbind, data.est2))
    colnames(data.estimates) <- paste(colnames(data), rep(order.deriv, each = ncol(data)), sep = ".")
    return(data.estimates)
  }
  
  #Compute derivatives per ID
  derivlist <- pbapply::pblapply(X = datalist, cl = cl,
                                 FUN = glla.multi,
                                 n.embed = n.embed, tau = tau, delta = delta, order = order)
  
  ### Estimating the dimensionality structure using EGA:
  
  message("Estimating the dimensionality structure using EGA...\n", appendLF = FALSE)
  
  for(i in 1:length(cases)){
    derivlist[[i]]$ID <- data[which(data[,id]==cases[i]),id][-c(1:(n.embed-1))]
  }
  
  names(derivlist) <- paste0("ID", cases)
  # Population Level:
  message("Level: Population...\n", appendLF = FALSE)
  
  data.all <- data.frame(Reduce(rbind, derivlist))
  
  # EGA Part
  
  if(use.derivatives == 0){
    ega1 <- EGA.estimate(data = data.all[,1:ncol(data[,-c(id)])],
                         model = model, model.args = model.args,
                         algorithm = algorithm, algorithm.args = algorithm.args,
                         corr = corr)}
  if(use.derivatives == 1){
    ega1 <- EGA.estimate(data = data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)],
                         model = model, model.args = model.args,
                         algorithm = algorithm, algorithm.args = algorithm.args,
                         corr = corr)}
  if(use.derivatives==2){
    init <- (ncol(data[,-c(id)])*2)+1
    cols <- seq(from = init, to = init+ncol(data[,-c(id)])-1)
    ega1 <- EGA.estimate(data = data.all[,cols],
                         model = model, model.args = model.args,
                         algorithm = algorithm, algorithm.args = algorithm.args,
                         corr = corr)}
  
  parallel::stopCluster(cl)
  
  # Level: Individual (intraindividual structure):
  message("Level: Individual (Intraindividual Structure)...\n", appendLF = FALSE)
  
  data.all <- data.frame(Reduce(rbind, derivlist))
  
  # Which derivatives to use:
  if(use.derivatives == 0){
    colstouse <- colnames(data.all[,1:ncol(data[,-c(id)])])}
  if(use.derivatives == 1){
    colstouse <- colnames(data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)])}
  if(use.derivatives==2){
    init <- (ncol(data[,-c(id)])*2)+1
    cols <- seq(from = init, to = init+ncol(data[,-c(id)])-1)
    colstouse <- colnames(data.all[,cols])
  }
  
  #initialize data list
  data.individuals <- vector("list", length = length(cases))
  data.individuals <- split(data.all[,colstouse],data.all$ID)
  names(data.individuals) <- paste0("ID", cases)
  
  
  #Parallel processing
  cl <- parallel::makeCluster(ncores)
  
  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("data.individuals", "cases"),
                          envir=environment())
  
  # EGA estimates per individual:
  ega.list.individuals <- list()
  
  ega.list.individuals <- pbapply::pblapply(X = data.individuals, cl = cl,
                                            FUN = EGA.estimate,
                                            model = model, model.args = model.args,
                                            algorithm = algorithm, algorithm.args = algorithm.args,
                                            corr = corr)
  parallel::stopCluster(cl)
  
  #let user know results have been computed
  message("done", appendLF = TRUE)
  
  # Results:
  results <- vector("list")
  results$Derivatives <- vector("list")
  results$Derivatives$Estimates <- derivlist
  results$Derivatives$EstimatesDF <- data.all
  results$dynEGA.pop <- ega1
  results$dynEGA.ind <- ega.list.individuals
  if(use.derivatives == 0){
    results$data.all <- data.all[,1:ncol(data[,-c(id)])]}
  if(use.derivatives == 1){
    results$data.all <- data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)]}
  if(use.derivatives == 2){
    results$data.all <- data.all[,cols]}
  results$data.individuals <- data.individuals
  results$model <- model
  class(results) <- "dynEGA.ind.pop"
  return(results)
}

#%%%%%%%%%
# UVA ----
#%%%%%%%%%

#' @noRd
# Redundancy Processing
# Updated 14.04.2021
redundancy.process <- function(data, cormat, n, model, method, type, sig, plot.redundancy, plot.args)
{
  # Compute redundancy method
  if(method == "irt"){
    
    mod <- mirt::mirt(data,1)
    sink <- capture.output(tom <- mirt::residuals(mod,type="Q3"))
    
  }else{
    
    if(method == "wto"){
      
      if(model == "glasso"){
        
        for(i in c(0.50, 0.25, 0))
        {
          net <- EBICglasso.qgraph(data = cormat, n = n, gamma = i)
          
          if(all(colSums(net)!=0))
          {break}
        }
        
      }else if(model == "tmfg"){
        
        net <- NetworkToolbox::TMFG(cormat)$A
        
      }else{
        
        stop(paste(model, "does not exist as an option for the argument 'model'"))
        
      }
      
      tom <- wTO::wTO(net, sign = "sign")
      
    }else if(method == "pcor"){
      
      tom <- -cov2cor(solve(cormat))
      
    }else{tom <- cormat}
    
  }
  
  # Number of variables; diagonal zero; absolute values
  vars <- ncol(tom); diag(tom) <- 0; tom <- abs(tom)
  
  # Lower triangle
  lower <- tom[lower.tri(tom, diag = FALSE)]
  
  # Names lower triangle
  name1 <- colnames(tom); name2 <- name1
  
  # Initialize name matrix
  name.matrix <- tom
  
  # Generate name matrix
  for(i in 1:vars)
    for(j in 1:vars){
      name.matrix[i,j] <- paste(name1[j], name2[i], sep = "--")
    }
  
  # Get lower triangle of names
  names(lower) <- name.matrix[lower.tri(name.matrix, diag = FALSE)]
  
  # Obtain positive values only
  pos.vals <- na.omit(ifelse(lower == 0, NA, lower)); attr(pos.vals, "na.action") <- NULL
  
  # Get redundant pairings
  if(type == "threshold"){## Threshold values
    redund <- pos.vals[which(pos.vals >= sig)]
    aic <- NULL
    g.dist <- NULL
  }else{## Determine distribution
    
    # Distributions, initialize AIC vector
    distr <- c("norm", "gamma"); aic <- numeric(length(distr)); names(aic) <- c("normal", "gamma")
    
    ## Obtain distribution
    for(i in 1:length(distr)){
      capture.output(
        aic[i] <- fitdistrplus::fitdist(pos.vals, distr[i], method="mle")$aic
      )
    }
    
    ## Obtain parameters
    g.dist <- suppressWarnings(MASS::fitdistr(pos.vals, names(aic)[which.min(aic)]))
    
    # Estimate p-values
    pval <- switch(names(aic)[which.min(aic)],
                   
                   normal = 1 - unlist(lapply(pos.vals, # positive values
                                              pnorm, # probability in normal distribution
                                              mean = g.dist$estimate["mean"], #mean of normal
                                              sd = g.dist$estimate["sd"]) #standard deviation of normal
                   ),
                   
                   gamma = 1 - unlist(lapply(pos.vals, # positive values
                                             pgamma, # probability in gamma distribution
                                             shape = g.dist$estimate["shape"], # shape of gamma
                                             rate = g.dist$estimate["rate"]) # rate of gamma
                   ),
    )
    
    # Check if using adaptive alpha
    if(type == "adapt"){
      sig <- NetworkToolbox::adapt.a("cor", alpha = sig, n = length(pos.vals), efxize = "medium")$adapt.a
    }
    
    # Get redundant pairings
    redund <- pos.vals[which(pval <= sig)]
    
  }
  
  # Check for redundant pairings
  if(length(redund) == 0){
    message("No redundant variables identified.")
    res.list <- NA
  }else{
    
    # Create result matrix
    split.res <- unlist(strsplit(names(redund), split = "--"))
    res.mat <- t(simplify2array(sapply(names(redund), strsplit, split = "--")))
    
    # Initialize result list
    res.list <- list()
    
    # Initialize count
    count <- 0
    
    while(nrow(res.mat)!=0){
      # Increase count
      count <- count + 1
      
      # Get variable counts
      var.counts <- sort(table(split.res), decreasing = TRUE)
      
      if(!all(var.counts==1)){
        # Identify targets
        target <- which(res.mat == names(var.counts[1]), arr.ind = TRUE)[,"row"]
        
        # Insert values into list
        res.list[[names(var.counts[1])]] <- setdiff(unlist(strsplit(names(target),split="--")),names(var.counts[1]))
        
        # Remove rows from result matrix
        res.mat <- res.mat[-target,]
        
        # Remove variables from split result
        split.res <- as.vector(res.mat)
        
        # Force matrix
        if(is.vector(res.mat))
        {res.mat <- t(as.matrix(res.mat))}
        
      }else{
        for(i in 1:nrow(res.mat))
        {res.list[[res.mat[i,1]]] <- unname(res.mat[i,2])}
        
        res.mat <- res.mat[-c(1:nrow(res.mat)),]
      }
    }
    
  }
  
  # Revert tom to matrix
  tom <- as.matrix(tom)
  
  # Check for plot
  if(plot.redundancy){
    
    # Initialize plot matrix
    plot.mat <- matrix(0, nrow = nrow(tom), ncol = ncol(tom))
    colnames(plot.mat) <- colnames(tom)
    row.names(plot.mat) <- colnames(tom)
    
    for(i in 1:length(res.list)){
      plot.mat[names(res.list)[i],res.list[[i]]] <- tom[names(res.list)[i],res.list[[i]]]
      plot.mat[res.list[[i]],names(res.list)[i]] <- tom[res.list[[i]],names(res.list)[i]]
    }
    
    if(any(colSums(plot.mat) == 0)){
      
      rm.mat <- which(colSums(plot.mat) == 0)
      plot.mat <- plot.mat[-rm.mat, -rm.mat]
      
    }
    
    plot.args$title <- switch(method,
                              "wto" = "Weighted\nTopological\nOverlap",
                              "pcor" = "Partial\nCorrelation",
                              "cor" = "Zero-order\nCorrelation",
                              "irt" = "IRT\nCorrelated\nResiudals"
    )
    
    if(ncol(plot.mat) <= 2){
      warning("No plot was produced because there are only two redundant variables")
    }else{
      
      # Global suppress warnings (need a better workaround)
      warn <- options("warn")[[1]]
      options(warn = -1)
        suppressMessages(
          net.plot <- redund.plot(plot.mat, plot.args)
        )
      options(warn = warn)
    }
    
  }
  
  # Get redundancy descriptives
  desc <- redund.desc(pos.vals = pos.vals, method = method, type = type, sig = sig,
                      aic = aic, g.dist = g.dist)
  
  # Add p-values
  if(type != "threshold"){
    desc$basic <- cbind(round(sig, 5), desc$basic)
    colnames(desc$basic)[1] <- "Sig"
    desc$centralTendency <- cbind(round(pval[row.names(desc$centralTendency)], 5),
          desc$centralTendency)
    colnames(desc$centralTendency)[1] <- "p-value"
  }
  
  # Results list
  res <- list()
  res$redundant <- res.list
  res$data <- data
  if(method != "irt"){res$correlation <- cormat}
  res$weights <- tom
  if(exists("net")){res$network <- net}
  if(exists("net.plot")){res$plot <- net.plot}
  res$descriptives <- desc
  res$method <- method
  res$model <- model
  res$type <- type
  if(type != "threshold"){res$distribution <- names(aic)[which.min(aic)]}
  
  class(res) <- "node.redundant"
  
  return(res)
  
}

#' @noRd
# Redundancy Naming
# Updated 13.12.2020
redund.names <- function(node.redundant.obj, key)
{
  # Check for node.redundant object class
  if(class(node.redundant.obj) != "node.redundant")
  {stop("A 'node.redundant' object must be used as input")}
  
  # Obtain and remove data from node redundant object
  data <- node.redundant.obj$data
  
  # Check that columns match key
  if(ncol(data) != length(as.vector(key)))
  {stop("Number of columns in data does not match the length of 'key'")}
  
  # Names of node.redundant object
  nr.names <- names(node.redundant.obj$redundant)
  
  # Key names
  key.names <- colnames(data)
  
  # Key change
  key.chn <- key
  
  for(i in 1:length(nr.names))
  {
    # Target redundant node
    target.r <- match(names(node.redundant.obj$redundant)[i],key.names)
    
    # Replace item name with description
    names(node.redundant.obj$redundant)[i] <- as.character(key.chn[target.r])
    
    # Target other nodes
    target.o <- match(node.redundant.obj$redundant[[i]],key.names)
    
    # Replace item names with description
    node.redundant.obj$redundant[[i]] <- as.character(key.chn[target.o])
  }
  
  # Create key code
  names(key) <- colnames(data)
  node.redundant.obj$key <- key
  
  return(node.redundant.obj)
}

#' @noRd
# Redundancy Descriptives
# Updated 13.12.2020
redund.desc <- function(pos.vals, method, type, sig, aic, g.dist)
{
  # Initialize descriptives matrix
  desc <- matrix(0, nrow = 1, ncol = 9)
  
  # Row name
  row.names(desc) <- switch(method,
                            "wto" = "wTO",
                            "pcor"= "pcor",
                            "cor" = "cor",
                            "irt" = "IRT"
  )
  
  colnames(desc) <- c("Mean", "SD", "Median", "MAD", "3*MAD", "6*MAD", "Minimum", "Maximum", "Critical Value")
  
  desc[,"Mean"] <- mean(pos.vals, na.rm = TRUE)
  desc[,"SD"] <- sd(pos.vals, na.rm = TRUE)
  desc[,"Median"] <- median(pos.vals, na.rm = TRUE)
  desc[,"MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE)
  desc[,"3*MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE) * 3
  desc[,"6*MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE) * 6
  desc[,"Minimum"] <- range(pos.vals, na.rm = TRUE)[1]
  desc[,"Maximum"] <- range(pos.vals, na.rm = TRUE)[2]
  
  # Critical value
  if(type == "threshold"){
    desc[,"Critical Value"] <- sig
  }else{
    
    desc[,"Critical Value"] <- switch(names(aic)[which.min(aic)],
                                      
                                      normal = qnorm(sig, #significance
                                                     mean = g.dist$estimate["mean"], #mean of normal
                                                     sd = g.dist$estimate["sd"], #sd of normal
                                                     lower.tail = FALSE),
                                      
                                      gamma = qgamma(sig, #significance
                                                     shape = g.dist$estimate["shape"], #shape of gamma
                                                     rate = g.dist$estimate["rate"], #rate of gamma
                                                     lower.tail = FALSE),
                                      
    )
    
  }
  
  # Organize positive values output
  ordered.pos <- sort(pos.vals, decreasing = TRUE)
  sd.from.mean <- (ordered.pos - mean(pos.vals, na.rm = TRUE)) / sd(ordered.pos, na.rm = TRUE)
  mad.from.median <- (ordered.pos - median(pos.vals, na.rm = TRUE)) / mad(ordered.pos, constant = 1, na.rm = TRUE)
  pos.output <- round(cbind(ordered.pos, sd.from.mean, mad.from.median), 3)
  
  colnames(pos.output)[1] <- switch(method,
                                    "wto" = "wTO",
                                    "pcor"= "pcor",
                                    "cor" = "cor",
                                    "irt" = "IRT"
  )
  
  colnames(pos.output)[2:3] <- c("SD from Mean", "MAD from Median")
  
  res.desc <- list()
  res.desc$basic <- round(desc, 3)
  res.desc$centralTendency <- pos.output
  
  return(res.desc)
  
}

#' @noRd
# Redundancy Plot
# Updated 13.12.2020
redund.plot <- function(plot.matrix, plot.args, plot.reduce = FALSE)
{
  # Convert to plot.mat
  plot.mat <- plot.matrix
  
  # weighted  network
  network1 <- network::network(plot.mat,
                               ignore.eval = FALSE,
                               names.eval = "weights",
                               directed = FALSE)
  
  if(isTRUE(plot.reduce)){
    wc <- c("Target", rep("Possible", ncol(plot.mat)-1))
    network::set.vertex.attribute(network1, attrname= "Communities", value = wc)
  }else{
    wc <- rep(plot.args$title, ncol(plot.mat))
    network::set.vertex.attribute(network1, attrname= "Communities", value = wc)
  }
  
  network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
  network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(plot.mat))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(rescale.edges(plot.mat, 5),
                                       nrow = nrow(plot.mat),
                                       ncol = ncol(plot.mat)))
  
  # Layout "Spring"
  graph1 <- NetworkToolbox::convert2igraph(plot.mat)
  edge.list <- igraph::as_edgelist(graph1)
  layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                             weights =
                                                               abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                             vcount = ncol(plot.mat))
  
  
  lower <- abs(plot.mat[lower.tri(plot.mat)])
  non.zero <- sqrt(lower[lower != 0])
  
  set.seed(1234)
  plot.args$net <- network1
  plot.args$node.color <- "Communities"
  plot.args$node.alpha <- plot.args$alpha
  plot.args$node.shape <- plot.args$shape
  node.size <- plot.args$node.size
  plot.args$node.size <- 0
  color.palette <- plot.args$color.palette
  plot.args$color.palette <- NULL
  plot.args$palette <- NULL
  plot.args$edge.lty <- "line"
  plot.args$edge.color <- "color"
  plot.args$edge.size <- "ScaledWeights"
  
  lower <- abs(plot.mat[lower.tri(plot.mat)])
  non.zero <- sqrt(lower[lower != 0])
  
  plot.args$edge.alpha <- non.zero
  plot.args$mode <- layout.spring
  plot.args$label <- colnames(plot.mat)
  plot.args$node.label <- rep("", ncol(plot.mat))
  if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
  if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
  
  redund.net <- suppressMessages(
    do.call(GGally::ggnet2, plot.args) + 
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc)))),
                                  breaks = sort(as.numeric(factor(wc)))) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ))
      )
  )
  
  redund.net <- redund.net + ggplot2::annotate("text", x = -Inf, y = Inf,
                                               hjust = 0, vjust = 1,
                                               label = plot.args$title, size = 5.5)
  
  set.seed(NULL)
  
  name <- colnames(plot.mat)
  
  # Custom nodes: transparent insides and dark borders
  redund.net <- redund.net + 
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                        color = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc))), sorted = FALSE),
                        shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
                        color = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc))), sorted = FALSE),
                        shape = 19, alpha = plot.args$alpha) +
    ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size)
  
  return(redund.net)
}

#' @importFrom graphics text
#' @noRd
# Redundancy Reduction
# Updated 15.02.2021
redund.reduce <- function(node.redundant.obj, reduce.method, plot.args, lavaan.args, corr)
{
  # Check for node.redundant object class
  if(class(node.redundant.obj) != "node.redundant")
  {stop("A 'node.redundant' object must be used as input")}
  
  # Redundant list
  redund <- node.redundant.obj$redundant
  
  # Copied data
  new.data <- node.redundant.obj$data
  
  # Weights
  if("network" %in% names(node.redundant.obj)){
    weights <- as.matrix(node.redundant.obj$network)
  }else if("correlation" %in% names(node.redundant.obj)){
    weights <- as.matrix(node.redundant.obj$correlation)
  }else{weights <- as.matrix(node.redundant.obj$weights)}
  
  # Track merged items
  merged <- list()
  
  # Track changed names
  name.chn <- vector("character")
  
  # Initialize count
  count <- 0
  
  # Get key
  if("key" %in% names(node.redundant.obj))
  {
    key <- node.redundant.obj$key
    names(key) <- names(node.redundant.obj$key)
  }else{
    key <- colnames(node.redundant.obj$data)
    names(key) <- key
  }
  
  # Line break
  linebreak()
  
  # Previous state
  prev.state <- list(redund)
  prev.state.data <- list(new.data)
  
  # Loop through named node redundant list
  while(length(redund) != 0)
  {
    # Targeting redundancy
    target.item <- names(redund)[1]
    
    # Potential redundancies
    pot <- redund[[1]]
    
    # Construct potential redundancies and menu
    poss <- redundancy.menu(redund = redund, reduce.method = reduce.method,
                            pot = pot, target.item = target.item,
                            weights = weights, plot.args = plot.args,
                            key = key, node.redundant.obj = node.redundant.obj)
    
    # Get input
    input <- input.check(poss, type = "redund")
    
    # Check if going back is necessary
    if(any(input == "b")){
      
      # Let user know they can't go back
      if(is.null(unlist(prev.state[length(prev.state) - 1], recursive = FALSE))){
        message("\nCannot go back. This is the start of the redundancy list")
      }else{
        
        # Renew redund list
        redund <- unlist(prev.state[length(prev.state) - 1], recursive = FALSE)
        
        # Remove previous state
        prev.state <- prev.state[-length(prev.state)]
        
        # Renew new data
        new.data <- as.data.frame(prev.state.data[length(prev.state.data) - 1])
        
        # Remove previous state data
        prev.state.data <- prev.state.data[-length(prev.state.data)]
        
        # Reduce count
        count <- count - 1
        
        # Remove merge
        merged <- merged[-length(merged)]
        
      }
      
    }else if(all(input != "0")){
      
      # Convert to numeric
      re.items <- as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ",")))
      
      # Items to combine with target
      comb <- poss[re.items]
      
      # Index items
      idx <- names(key)[match(comb, key)]
      
      # Target index
      tar.idx <- names(key)[match(target.item, key)]
      
      # Update merged list
      count <- count + 1
      merged[[count]] <- c(key[tar.idx], key[idx])
      
      # Combine into target index
      if(reduce.method == "latent"){
        
        # Latent variable
        ## create model
        mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,c(tar.idx, idx)]), collapse = " + "))
        
        # Replace arguments
        lavaan.args$model <- mod
        lavaan.args$data <- new.data
        ## Get default estimator
        categories <- apply(new.data[,c(tar.idx, idx)], 2, function(x){
          length(unique(x))
        })
        
        ## Check categories
        if(any(categories < 6)){# Not all continuous
          lavaan.args$estimator <- "WLSMV"
          lavaan.args$missing <- "pairwise"
        }else{# All can be considered continuous
          lavaan.args$estimator <- "MLR"
          lavaan.args$missing <- "fiml"
        }
        
        ## get CFA function from lavaan
        FUN <- lavaan::cfa
        
        ## fit model
        fit <- suppressWarnings(
          suppressMessages(
            do.call(what = "FUN", args = as.list(lavaan.args))
          )
        )
        
        ## identify cases
        cases <- lavaan::inspect(fit, "case.idx")
        
        ## compute latent variable score
        latent <- as.numeric(lavaan::lavPredict(fit))
        
        ## check for missing cases and handle
        if(length(cases) != nrow(new.data))
        {
          new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
          new.vec[cases] <- latent
        }else{new.vec <- latent}
        
        ## check for reverse scoring/labelling
        corrs <- as.matrix(cor(cbind(latent,new.data[,c(tar.idx, idx)]), use = "complete.obs")[1,-1])
        row.names(corrs) <- c(key[tar.idx], key[idx])
        colnames(corrs) <- "latent"
        
        if(any(sign(corrs)==-1))
        {
          message("Some variables are reverse coded (negative correlations with latent variable were found). Correlations with latent variable:")
          print(corrs)
          
          input2 <- "o"
          
          while(input2 != "y" && input2 != "n")
          {input2 <- readline("Reverse code for positive labelling (y/n): ")}
          
          if(input2 == "y")
          {new.vec <- -new.vec}
        }
        
        # input new vector
        new.data[,tar.idx] <- new.vec
        
        # Ask for new label
        lab <- readline(prompt = "New label for latent variable (no quotations): ")
        name.chn[count] <- lab
        col.idx <- match(tar.idx, colnames(new.data))
        colnames(new.data)[col.idx] <- lab
        
        message(paste("\nNew LATENT variable called '", lab,"' was created. Redundant variables were REMOVED", sep = ""))
        
      }else if(reduce.method == "remove" | reduce.method == "sum"){
        
        target.key <- c(tar.idx, idx)
        target.data <- new.data[,target.key]
        
        means <- round(colMeans(target.data, na.rm = TRUE), 2)
        sds <- round(apply(target.data, 2, sd, na.rm = TRUE), 2)
        ranges <- round(apply(target.data, 2, range, na.rm = TRUE), 2)
        
        if(ncol(target.data) > 2){
          
          # Corrected item-total correlations
          cor.corr <- round(item.total(target.data, corr), 2)
          
          ## Use information utility?
          categories <- apply(target.data, 2, function(x){
            length(unique(x))
          })
          
          ## Check categories
          if(all(categories <= 7)){
            
            ## Use
            #util <- round(info.util(target.data), 2)
            
            tab <- cbind(#util,
              cor.corr, means, sds, t(ranges))
            colnames(tab) <- c(#"Utility Gain",
              "Item-Total r", "Mean", "SD", "Low", "High")
          }else{
            
            tab <- cbind(cor.corr, means, sds, t(ranges))
            colnames(tab) <- c("Item-Total r", "Mean", "SD", "Low", "High")
            
          }
          
        }else{
          tab <- cbind(means, sds, t(ranges))
          colnames(tab) <- c("Mean", "SD", "Low", "High")
        }
        row.names(tab) <- c("0 (Target)", 1:length(comb))
        
        tab[,1:(ncol(tab) - 2)] <- matrix(sprintf("%.2f", tab[,1:(ncol(tab) - 2)]), nrow = nrow(tab), ncol = ncol(tab) - 2)
        
        if(!isSymmetric(new.data)){
          gridExtra::grid.arrange(gridExtra::tableGrob(tab))
        }
        
        # Input check
        new.input <- input.check(poss = c(target.item, comb), type = "remove")
        
        # All variables
        ind <- names(key[match(c(target.item, comb), key)])
        
        # All except selected variable
        idx <- ind[-(as.numeric(new.input)+1)]
        
        # Input selected
        merged[[count]] <- key[idx]
        
        # Name merged input
        name.chn[count] <- key[setdiff(ind, idx)]
        
        # Message user
        message(paste("\nKEPT '", key[ind[as.numeric(new.input) + 1]],"' and REMOVED all others", sep = ""))
        
      }
      
      # Remove redundant variables from data
      rm.idx <- match(idx, colnames(new.data))
      
      if(isSymmetric(new.data)){
        new.data <- new.data[-rm.idx, -rm.idx]
      }else{
        new.data <- new.data[,-rm.idx]
      }
      
      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)
      
      # Remove target item
      redund[[1]] <- NULL
      
      # Make ind
      if(reduce.method == "latent"){ind <- idx}
      
      # Remove variables within future options
      for(i in 1:length(ind)){
        
        ## Get elements from other redundancy lists
        elements <- sapply(redund, function(x){key[ind[i]] %in% x})
        
        ## If there are any
        if(any(elements)){
          
          ### Target elements
          target.elem <- which(elements)
          
          for(j in 1:length(target.elem)){
            
            #### Target list
            list.target <- redund[[target.elem[j]]]
            
            #### Remove from target list
            redund[[target.elem[j]]] <- redund[[target.elem[j]]][-which(key[ind[i]] == list.target)]
            
          }
          
        }
        
      }
      
      # Remove empty list objects
      rm.list <- which(
        unlist(lapply(redund, function(x){
          length(x)
        })) == 0
      )
      
      if(length(rm.list) != 0){
        redund <- redund[-rm.list]
      }
      
      # Remove object names
      if(any(key[ind] %in% names(redund))){
        
        ## Get target names
        name.targets <- which(key[ind] %in% names(redund))
        
        ## Loop through
        for(i in 1:length(name.targets)){
          
          # If there is only one item left, then remove from list
          if(length(redund[[key[ind][name.targets[i]]]]) == 1){
            redund[key[ind][name.targets[i]]] <- NULL
          }else{# Otherwise, replace name with first element and remove first element from list
            names(redund[key[ind]][name.targets[i]]) <- redund[[key[ind][name.targets[i]]]][1]
            redund[[key[ind][name.targets[i]]]][1] <- NA
            redund[[key[ind][name.targets[i]]]] <- na.omit(redund[[key[ind][name.targets[i]]]])
            redund[key[ind][name.targets[i]]] <- NULL
          }
          
        }
        
      }
      
      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)
      
    }else{
      
      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)
      
      redund[[1]] <- NULL
      
      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)
      
    }
    
    if(!is.null(input)){
      linebreak()
      input <- NULL
    }
    
    # Artificial pause for smoothness of experience
    Sys.sleep(1)
    
  }
  
  # Transform merged list to matrix
  if(length(merged) != 0)
  {
    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))
    
    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))
    
    # Input into merged matrix
    for(i in 1:length(merged))
    {
      diff <- m.rows - length(merged[[i]])
      
      m.mat[,i] <- c(merged[[i]], rep("", diff))
    }
    
    colnames(m.mat) <- name.chn
  }
  
  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key)))
  {
    # Target names
    target.names <- which(colnames(new.data) %in% names(key))
    
    # new.data names
    new.data.names <- colnames(new.data)[target.names]
    
    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]
  }
  
  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)
    
    if(reduce.method == "latent"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else if(reduce.method == "remove" | reduce.method == "sum"){
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }
  
  # Check for "sum"
  if(reduce.method == "sum"){
    
    # Reinstate new.data
    new.data <- node.redundant.obj$data
    
    # Collapse across rows
    for(i in 1:nrow(m.mat)){
      
      # Collapse
      collapse <- row.names(m.mat)[i]
      
      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]
      
      # Collapse and insert into matrix
      new.data[,collapse] <- rowSums(new.data[,c(collapse, redunds)])
      
      # Remove redundant terms
      new.data <- new.data[,-match(redunds, colnames(new.data))]
    }
    
  }
  
  
  
  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat
  
  return(res)
  
}

#' @noRd
# Redundancy Reduction (Automated)
# Updated 20.06.2021
redund.reduce.auto <- function(node.redundant.obj,
                               reduce.method, lavaan.args, corr)
{
  # Check for node.redundant object class
  if(class(node.redundant.obj) != "node.redundant")
  {stop("A 'node.redundant' object must be used as input")}
  
  # Redundant list
  redund <- node.redundant.obj$redundant
  
  # Copied data
  new.data <- node.redundant.obj$data
  
  # Track merged items
  merged <- list()
  
  # Track changed names
  name.chn <- vector("character")
  
  # Initialize count
  count <- 0
  
  # Get key
  if("key" %in% names(node.redundant.obj))
  {
    key <- node.redundant.obj$key
    names(key) <- names(node.redundant.obj$key)
  }else{
    key <- colnames(node.redundant.obj$data)
    names(key) <- key
  }
  
  # Line break
  #linebreak()
  
  # Previous state
  prev.state <- list(redund)
  prev.state.data <- list(new.data)
  
  # Loop through named node redundant list
  while(length(redund) != 0)
  {
    # Targeting redundancy
    target.item <- names(redund)[1]
    
    # Potential redundancies
    pot <- redund[[1]]
    
    # Get input
    input <- as.character(1:length(pot))
    
    # Auto
    if(all(input != "0")){
      
      # Convert to numeric
      re.items <- as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ",")))
      
      # Items to combine with target
      comb <- pot
      
      # Index items
      idx <- names(key)[match(comb, key)]
      
      # Target index
      tar.idx <- names(key)[match(target.item, key)]
      
      # Update merged list
      count <- count + 1
      merged[[count]] <- c(key[tar.idx], key[idx])
      
      # Combine into target index
      if(reduce.method == "latent"){
        
        # Latent variable
        ## create model
        mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,c(tar.idx, idx)]), collapse = " + "))
        
        # Replace arguments
        lavaan.args$model <- mod
        lavaan.args$data <- new.data
        ## Get default estimator
        categories <- apply(new.data[,c(tar.idx, idx)], 2, function(x){
          length(unique(x))
        })
        
        ## Check categories
        if(any(categories < 6)){# Not all continuous
          lavaan.args$estimator <- "WLSMV"
          lavaan.args$missing <- "pairwise"
        }else{# All can be considered continuous
          lavaan.args$estimator <- "MLR"
          lavaan.args$missing <- "fiml"
        }
        
        ## get CFA function from lavaan
        FUN <- lavaan::cfa
        
        ## fit model
        fit <- suppressWarnings(
          suppressMessages(
            do.call(what = "FUN", args = as.list(lavaan.args))
          )
        )
        
        ## identify cases
        cases <- lavaan::inspect(fit, "case.idx")
        
        ## compute latent variable score
        latent <- as.numeric(lavaan::lavPredict(fit))
        
        ## check for missing cases and handle
        if(length(cases) != nrow(new.data))
        {
          new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
          new.vec[cases] <- latent
        }else{new.vec <- latent}
        
        ## check for reverse scoring/labelling
        corrs <- as.matrix(cor(cbind(latent,new.data[,c(tar.idx, idx)]), use = "complete.obs")[1,-1])
        row.names(corrs) <- c(key[tar.idx], key[idx])
        colnames(corrs) <- "latent"
        
        # input new vector
        new.data[,tar.idx] <- new.vec
        
        # Ask for new label
        lab <- paste("LV_", count, sep = "")
        name.chn[count] <- lab
        col.idx <- match(tar.idx, colnames(new.data))
        colnames(new.data)[col.idx] <- lab
        
      }
      
      # Remove redundant variables from data
      rm.idx <- match(idx, colnames(new.data))
      
      if(isSymmetric(new.data)){
        new.data <- new.data[-rm.idx, -rm.idx]
      }else{
        new.data <- new.data[,-rm.idx]
      }
      
      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)
      
      # Remove target item
      redund[[1]] <- NULL
      
      # Make ind
      if(reduce.method == "latent"){ind <- idx}
      
      # Remove variables within future options
      for(i in 1:length(ind)){
        
        ## Get elements from other redundancy lists
        elements <- sapply(redund, function(x){key[ind[i]] %in% x})
        
        ## If there are any
        if(any(elements)){
          
          ### Target elements
          target.elem <- which(elements)
          
          for(j in 1:length(target.elem)){
            
            #### Target list
            list.target <- redund[[target.elem[j]]]
            
            #### Remove from target list
            redund[[target.elem[j]]] <- redund[[target.elem[j]]][-which(key[ind[i]] == list.target)]
            
          }
          
        }
        
      }
      
      # Remove empty list objects
      rm.list <- which(
        unlist(lapply(redund, function(x){
          length(x)
        })) == 0
      )
      
      if(length(rm.list) != 0){
        redund <- redund[-rm.list]
      }
      
      # Remove object names
      if(any(key[ind] %in% names(redund))){
        
        ## Get target names
        name.targets <- which(key[ind] %in% names(redund))
        
        ## Loop through
        for(i in 1:length(name.targets)){
          
          # If there is only one item left, then remove from list
          if(length(redund[[key[ind][name.targets[i]]]]) == 1){
            redund[key[ind][name.targets[i]]] <- NULL
          }else{# Otherwise, replace name with first element and remove first element from list
            names(redund[key[ind]][name.targets[i]]) <- redund[[key[ind][name.targets[i]]]][1]
            redund[[key[ind][name.targets[i]]]][1] <- NA
            redund[[key[ind][name.targets[i]]]] <- na.omit(redund[[key[ind][name.targets[i]]]])
            redund[key[ind][name.targets[i]]] <- NULL
          }
          
        }
        
      }
      
      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)
      
    }else{
      
      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)
      
      redund[[1]] <- NULL
      
      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)
      
    }
    
    if(!is.null(input)){
      #linebreak()
      input <- NULL
    }
    
    # Artificial pause for smoothness of experience
    #Sys.sleep(1)
    
  }
  
  # Transform merged list to matrix
  if(length(merged) != 0)
  {
    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))
    
    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))
    
    # Input into merged matrix
    for(i in 1:length(merged))
    {
      diff <- m.rows - length(merged[[i]])
      
      m.mat[,i] <- c(merged[[i]], rep("", diff))
    }
    
    colnames(m.mat) <- name.chn
  }
  
  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key)))
  {
    # Target names
    target.names <- which(colnames(new.data) %in% names(key))
    
    # new.data names
    new.data.names <- colnames(new.data)[target.names]
    
    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]
  }
  
  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)
    
    if(reduce.method == "latent"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else if(reduce.method == "remove" | reduce.method == "sum"){
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }
  
  # Check for "sum"
  if(reduce.method == "sum"){
    
    # Reinstate new.data
    new.data <- node.redundant.obj$data
    
    # Collapse across rows
    for(i in 1:nrow(m.mat)){
      
      # Collapse
      collapse <- row.names(m.mat)[i]
      
      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]
      
      # Collapse and insert into matrix
      new.data[,collapse] <- rowSums(new.data[,c(collapse, redunds)])
      
      # Remove redundant terms
      new.data <- new.data[,-match(redunds, colnames(new.data))]
    }
    
  }
  
  
  
  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat
  
  return(res)
  
}

#' @noRd
# Redundancy Adhoc Reduction (Automated)
# Updated 22.07.2021
redund.adhoc.auto <- function(node.redundant.obj,
                              node.redundant.reduced,
                              node.redundant.original,
                              reduce.method, lavaan.args, corr)
{
  
  # Redundant list
  redund <- node.redundant.obj$redundant
  
  # Check for overlaps
  overlap <- unlist(
    lapply(redund, function(x, name){
      any(name %in% x)
    }, name = names(redund))
  )
  
  # Merge overlaps
  if(length(which(overlap)) != 0){
    
    for(i in which(overlap)){
      redund[[i]] <- unique(unlist(c(redund[[i]], redund[names(redund) %in% redund[[i]]])))
      
      redund <- redund[-which(names(redund) %in% redund[[i]])]
    }
    
  }

  # Copied data
  new.data <- node.redundant.original$data
  
  # Get reduced
  reduced.merged <- node.redundant.reduced$merged
  
  # Make merged
  merged <- lapply(apply(reduced.merged, 1, as.list), function(x){unname(unlist(x))})
  
  # Set counter
  count <- max(as.numeric(gsub("LV_", "", names(merged))))
  
  # Update merged
  for(i in 1:length(redund)){
    
    # Target variables
    target <- c(names(redund)[i], redund[[i]])
    
    # Expand latent variables
    lv <- target[grep("LV_", target)]
    
    # If any latent variables
    if(length(lv) != 0){
      
      # Create new latent variable
      new_LV <- unname(c(unlist(merged[lv]), target[-grep("LV_", target)]))
      
      # Insert into merged
      count <- count + 1
      merged[[paste("LV_", count, sep = "")]] <- new_LV
      
      # Remove old latent variables
      merged[lv] <- NULL
      
    }else{
      
      # Create new latent variable
      new_LV <- target
      
      # Insert into merged
      count <- count + 1
      merged[[paste("LV_", count, sep = "")]] <- new_LV
      
      # Remove old variables
      merged[names(redund)[i]] <- NULL
      
    }
  }
  
  # Remove all ""
  merged <- lapply(merged, function(x){
    x <- na.omit(ifelse(x == "", NA, x))
    attr(x, which = "na.action") <- NULL
    return(x)
  })
  
  # Reorder longest to shortest
  merged <- merged[order(unlist(lapply(merged, length)), decreasing = TRUE)]
  
  # Get key
  if("key" %in% names(node.redundant.original))
  {
    key <- node.redundant.original$key
    names(key) <- names(node.redundant.original$key)
  }else{
    key <- colnames(node.redundant.original$data)
    names(key) <- key
  }
  
  # Remove missing reundancies
  merged <- lapply(merged, function(x){
    if(length(x) == 0){
      NULL
    }else{x}
  })
  
  nulls <- unlist(lapply(merged, is.null))
  
  if(any(nulls)){
    merged <- merged[!nulls]
  }
  
  # Loop through to make new variables
  for(i in 1:length(merged)){
    
    # Combine into target index
    if(reduce.method == "latent"){
      
      # Get indices from key
      idx <- names(key[match(merged[[i]], key)])
  
      # Create model
      mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,idx]), collapse = " + "))
      
      # Replace arguments
      lavaan.args$model <- mod
      lavaan.args$data <- new.data
      ## Get default estimator
      categories <- apply(new.data[,idx], 2, function(x){
        length(unique(x))
      })
      
      ## Check categories
      if(any(categories < 6)){# Not all continuous
        lavaan.args$estimator <- "WLSMV"
        lavaan.args$missing <- "pairwise"
      }else{# All can be considered continuous
        lavaan.args$estimator <- "MLR"
        lavaan.args$missing <- "fiml"
      }
      
      ## get CFA function from lavaan
      FUN <- lavaan::cfa
      
      ## fit model
      fit <- suppressWarnings(
        suppressMessages(
          do.call(what = "FUN", args = as.list(lavaan.args))
        )
      )
      
      ## identify cases
      cases <- lavaan::inspect(fit, "case.idx")
      
      ## compute latent variable score
      latent <- as.numeric(lavaan::lavPredict(fit))
      
      ## check for missing cases and handle
      if(length(cases) != nrow(new.data)){
        new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
        new.vec[cases] <- latent
      }else{new.vec <- latent}
      
      # Remove variables
      new.data <- new.data[,-match(idx, colnames(new.data))]
      
      # Tack on latent variable
      new.data <- cbind(new.data, new.vec)
      
      # Rename latent variable
      colnames(new.data)[ncol(new.data)] <- names(merged)[i]
      
    }
    
  }
  
  # Transform merged list to matrix
  if(length(merged) != 0){
    
    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))
    
    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))
    
    # Input into merged matrix
    for(i in 1:length(merged)){
      
      diff <- m.rows - length(merged[[i]])
      
      m.mat[,i] <- c(merged[[i]], rep("", diff))
      
    }
    
    colnames(m.mat) <- names(merged)
    
  }
  
  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key))){
    
    # Target names
    target.names <- which(colnames(new.data) %in% names(key))
    
    # new.data names
    new.data.names <- colnames(new.data)[target.names]
    
    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]
  }
  
  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)
    
    if(reduce.method == "latent"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else if(reduce.method == "remove" | reduce.method == "sum"){
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }
  
  # Check for "sum"
  if(reduce.method == "sum"){
    
    # Reinstate new.data
    new.data <- node.redundant.original$data
    
    # Collapse across rows
    for(i in 1:nrow(m.mat)){
      
      # Collapse
      collapse <- row.names(m.mat)[i]
      
      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]
      
      # Collapse and insert into matrix
      new.data[,collapse] <- rowSums(new.data[,c(collapse, redunds)], na.rm = TRUE)
      
      # Remove redundant terms
      new.data <- new.data[,-match(redunds, colnames(new.data))]
    }
    
  }
  
  
  
  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat
  
  return(res)
  
}

#' @noRd
# Item-total correlations
# For UVA
# Updated 15.02.2021
item.total <- function (data.sub, corr)
{
  # Get correlations
  corrs <- switch(corr,
                  "cor_auto" = suppressMessages(qgraph::cor_auto(data.sub)),
                  "pearson" = suppressMessages(cor(data.sub, use = "pairwise.complete.obs")),
                  "spearman" = suppressMessages(cor(data.sub, method = "spearman", use = "pairwise.complete.obs"))
  )
  
  # Check for negatives (reverse if so)
  for(i in 1:nrow(corrs)){
    
    if(any(sign(corrs[i,-i]) == -1)){
      data.sub[,i] <- (max(data.sub[,i]) + min(data.sub[,i])) - data.sub[,i]
      
      if(any(sign(corrs[i,-i]) == 1)){
        target <- colnames(corrs)[(sign(corrs[i,]) == 1)[-i]]
        
        for(j in 1:length(target)){
          data.sub[,target[j]] <- (max(data.sub[,target[j]]) + min(data.sub[,target[j]])) - data.sub[,target[j]]
        }
        
      }
      
      corrs <- switch(corr,
                      "cor_auto" = suppressMessages(qgraph::cor_auto(data.sub)),
                      "pearson" = suppressMessages(cor(data.sub, use = "pairwise.complete.obs")),
                      "spearman" = suppressMessages(cor(data.sub, method = "spearman", use = "pairwise.complete.obs"))
      )
    }
    
  }
  
  # Compute corrected item-total correlations
  cor.corr <- numeric(ncol(data.sub))
  
  # Loop through
  for(i in 1:ncol(data.sub)){
    
    cor.corr[i] <- switch(corr,
                          "cor_auto" = suppressMessages(qgraph::cor_auto(cbind(data.sub[,i], rowSums(data.sub[,-i]))))[1,2],
                          "pearson" = suppressMessages(cor(cbind(data.sub[,i],
                                                                 rowSums(data.sub[,-i])),
                                                           method = "pearson",
                                                           use = "pairwise.complete.obs"))[1,2],
                          "spearman" = suppressMessages(cor(cbind(data.sub[,i],
                                                                  rowSums(data.sub[,-i])),
                                                            method = "spearman",
                                                            use = "pairwise.complete.obs"))[1,2]
    )
    
  }
  
  return(cor.corr)
  
}

#' @noRd
# Menu for redundancy
# For UVA
# Updated 15.02.2021
redundancy.menu <- function (redund, reduce.method, pot, target.item, weights,
                             plot.args, key, node.redundant.obj)
{
  # Configure into list
  pot <- list(pot)
  names(pot) <- target.item
  possible <- unname(unlist(pot))
  
  # Loop through potential redundancies
  count2 <- 1
  
  # Check names
  for(i in 1:length(possible)){
    
    if(possible[i] %in% names(redund)){
      count2 <- count2 + 1
      pot[count2] <- redund[possible[i]]
      names(pot)[count2] <- possible[i]
    }
    
  }
  
  # Check elements
  for(i in 1:length(possible)){
    elements <- redund[sapply(redund, function(x){possible[i] %in% x})]
    pot[(count2 + 1):(count2 + length(elements))] <- elements
    names(pot)[(count2 + 1):(count2 + length(elements))] <- names(elements)
    count2 <- count2 + length(elements)
  }
  
  # Get unique lists
  pot <- redund[unique(names(pot))]
  
  # Possible options
  poss <- unique(c(unname(unlist(pot)), names(pot)[-1]))
  
  # Organize plot of redundancy connections
  mat <- matrix(0, nrow = length(poss) + 1, ncol = length(poss) + 1)
  colnames(mat) <- c(paste("Target"), 1:length(poss))
  row.names(mat) <- colnames(mat)
  
  mat["Target",paste(1:length(unlist(pot[[1]])))] <- weights[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]
  mat[paste(1:length(unlist(pot[[1]]))),"Target"] <- weights[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]
  
  if(length(pot) != 1)
  {
    # Remove first element
    ext <- pot[-1]
    
    # Loop through rest of extended
    for(i in 1:length(ext))
    {
      # Target extended
      target.ext <- ext[[i]]
      
      # Loop through target
      for(j in 1:length(target.ext))
      {
        # Single out each element
        single <- target.ext[[j]]
        
        # Get element in possible redundancies
        elem <- match(names(ext)[i], poss)
        
        # Get elements redundant with element
        red.elem <- match(single, poss)
        
        # Put into matrix
        mat[paste(elem),paste(red.elem)] <- weights[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
        mat[paste(red.elem),paste(elem)] <- weights[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
      }
    }
  }
  
  # Print target and potential options
  cat(paste("Target variable: '", target.item, "'", sep = ""))
  cat("\n\nPotential redundancies:\n\n")
  if(reduce.method == "latent"){
    cat("0. Do not combine with any")
  }else if(reduce.method == "remove" | reduce.method == "sum"){
    cat("0. None")
  }
  
  cat(paste("\n", 1:length(poss), ". ", "'", poss, "'", sep = ""),"\n")
  
  # Plot
  if(node.redundant.obj$model == "tmfg"){
    
    plot.args$title <- "Zero-order Correlations"
    
  }else{
    
    plot.args$title <- switch(node.redundant.obj$method,
                              "wto" = "Regularized Partial Correlations",
                              "pcor" = "Partial Correlations",
                              "cor" = "Zero-order Correlations",
                              "irt" = "Correlated Residuals"
    )
    
  }
  
  if(length(poss) > 1){
    # Global suppress warnings (need a better work around)
    warn <- options("warn")[[1]]
    options(warn = -1)
    plot(redund.plot(plot.matrix = mat, plot.args = plot.args, plot.reduce = TRUE))
    options(warn = warn)
  }else{
    
    if(node.redundant.obj$model == "tmfg"){
      
      plot.args$title <- "Zero-order Correlation"
      
    }else{
      
      plot.args$title <- switch(node.redundant.obj$method,
                                "wto" = "Regularized Partial Correlation",
                                "pcor" = "Partial Correlation",
                                "cor" = "Zero-order Correlation",
                                "irt" = "Correlated Residuals"
      )
      
    }
    
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("There is only one redundant variable with the target variable.\nTheir ",
                                 tolower(plot.args$title), " = ", round(mat[1,2], 3),
                                 sep = ""),
         cex = 1.6, col = "black")
    par(mar = c(5, 4, 4, 2) + 0.1)
    
  }
  
  message("\nPress 'B' to go back")
  
  return(poss)
}

#' @noRd
# Input check for redundancy
# For UVA
# Updated 21.12.2020
input.check <- function (poss, type = c("redund", "remove"))
{
  
  if(type == "redund"){
    
    message("\nEnter numbers of variables redundant with the target variable (separate by commas)")
    input <- readline(prompt = "Selection: ")
    
    # Redo input check
    re.input <- in.check(input, poss = poss)
    
    while(re.input){
      # Print message to try again
      message("Inappropriate input. Try again. 'B' can be used to go back.\n")
      
      # Get input
      message("Enter numbers of variables redundant with the target variable (separate by commas)")
      input <- readline(prompt = "Selection: ")
      
      # Redo input check
      re.input <- in.check(input, poss)
    }
    
  }else if (type == "remove"){
    
    cat(
      paste("\n0. ", "'", poss[1], "'", " (Target)", sep = ""),
      paste("\n", 1:(length(poss) - 1), ". ", "'", poss[-1], "'", sep = ""),
      "\n\n"
    )
    
    input <- readline(prompt = "Select variable to KEEP: ")
    
    # Redo input check
    re.input <- in.check(input, poss = poss)
    
    while(re.input)
    {
      # Print message to try again
      message("Inappropriate input. Try again.\n")
      
      # Get input
      input <- readline(prompt = "Select variable to KEEP: ")
      
      # Redo input check
      re.input <- in.check(input, poss)
    }
    
  }
  
  return(input)
}

#' @noRd
# Input check for redundancy
# For UVA
# Updated 21.12.2020
in.check <- function(input, poss)
{
  if(tolower(input) == "b"){
    ret.val <- FALSE
  }else{
    inp <- suppressWarnings(as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ","))))
    
    ret.val <- FALSE
    
    if(any(is.na(inp)))
    {ret.val <- TRUE}
    
    if(length(inp) == 0)
    {ret.val <- TRUE}
    
    if(length(setdiff(inp, 0:length(poss))) != 0)
    {ret.val <- TRUE}
  }
  
  return(ret.val)
}

#' @noRd 
# Creates a nice looking line break
# For UVA
# Updated 26.02.2021
linebreak <- function(){cat("\n", colortext(paste(rep("-", getOption("width")), collapse = ""), defaults = "message"), "\n\n")}

#' @noRd 
# Changes names that have leading characters that are numbers
# This creates an error in lavaan's formulas
# Numbers are moved to the end of the name
# For UVA
# Updated 24.03.2021
lavaan.formula.names <- function (data){
  
  # Original column names
  original.names <- colnames(data)
  
  # Move leading numeric values to the end of the variable name
  colnames(data) <- unlist(
    lapply(strsplit(colnames(data),
                    split = ""), function(x){
                      
                      ind <- grepl("[[:digit:]]", x)
                      
                      if(isTRUE(rle(ind)$values[1])){
                        rm.ind <- x[1:rle(ind)$lengths[1]]
                        new.name <- paste(
                          paste(x[-c(1:rle(ind)$lengths[1])], collapse = ""),
                          paste(rm.ind, collapse = ""),
                          sep = "_"
                        )
                        return(new.name)
                      }else{return(paste(x, collapse = ""))}
                      
                    })
  )
  
  # Message user so they know
  if(any(colnames(data) != original.names)){
    
    message(paste("Some variable names begin with a number.",
                  "This creates errors in 'lavaan' formulas.",
                  styletext("\nAll numbers have been moved to the end of the variable name.", defaults = "bold"),
                  "\nTo avoid this message, make sure all variable names begin with a character.",
                  sep = "\n"))
    
  }
  
  return(data)
  
}

#%%%%%%%%%%%%%%%%%%%
# itemStability ----
#%%%%%%%%%%%%%%%%%%%

#' @noRd
# Converts memberships to numeric
# For itemStability
# Updated 26.02.2021
numeric.membership <- function(membership){
  
  # Get unique membership
  unique.membership <- unique(membership)
  
  # Check for named memberships
  if(all(is.character(unique.membership))){
    
    # Assign numeric values to memberships
    membership.numbers <- membership
    
    for(i in 1:length(unique.membership)){
      membership.numbers[which(membership.numbers == unique.membership[i])] <- i
    }
    
  }else{membership.numbers <- membership}
  
  # Make sure the memberships are numeric
  membership.numbers <- as.numeric(membership.numbers)
  
  # Name elements
  names(membership.numbers) <- names(membership)
  
  return(membership.numbers)
  
}

#' Computes the mode
#'
#' @param v Numeric vector.
#' Vector of values to find mode in
#'
#' @param fin.vec Alphanumeric vector.
#' Vector of current state of \code{v}
#'
#' @return The mode of a vector
#'
#' @noRd
#'
# Mode
# Updated 15.06.2020
mode <- function(v, fin.vec)
{
  #unique values
  uniqv <- unique(v)
  
  #find mode
  uniq.val <- uniqv[which.max(tabulate(match(v, uniqv)))]
  
  #do not overwrite already identified dimension
  while(uniq.val %in% fin.vec)
  {
    #remove unique value
    uniqv <- uniqv[-which(uniq.val==uniqv)]
    
    if(length(uniqv)==0)
    {
      uniq.val <- NA
      break
    }
    
    #find mode
    uniq.val <- uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  return(uniq.val)
}

#' Converts membership vector into a target membership vector
#'
#' @param target.wc Numeric vector.
#' Target membership vector
#'
#' @param convert.wc Numeric vector, matrix, or data frame.
#' Memberships to be converted to the target membership vector
#'
#' @return Vector or matrix of homogenized memberships
#'
#' @noRd
#'
# Homogenize Membership
# For itemStability
# Updated 14.04.2021
homogenize.membership <- function (target.wc, convert.wc)
{
  # Obtain whether vector or matrix is input for 'convert.wc'
  if(is.vector(convert.wc))
  {
    
    # Error if not same length
    if(length(convert.wc) != length(target.wc))
    {stop("'convert.wc' must be same length as 'target.wc'")}
    
    # Number of comparisons
    n <- 1
    
  }else{
    
    # Make sure 'convert.wc' is matrix
    convert.wc <- as.matrix(convert.wc)
    
    # Number of comparisons
    n <- ncol(convert.wc)
    
  }
  
  # Initialize conversion matrix
  convert.mat <- matrix(NA, nrow = length(target.wc), ncol = n)
  ## Get node names
  if(!is.null(names(target.wc)))
  {row.names(convert.mat) <- names(target.wc)}
  
  # Identify target membership within bootstrapped memberships
  for(i in 1:n)
  {
    # New membership vector
    new.vec <- convert.wc[,i]
    
    # Unique new membership
    new.uniq <- unique(new.vec)
    
    # Check if dimensionality solution was reached
    if(max(new.vec, na.rm = TRUE) == length(new.vec)){
      
      # Submit NAs
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)
      
    }else if(max(target.wc, na.rm = TRUE) > max(new.vec, na.rm = TRUE)){
      # Converge based on maximum number of dimensions
      
      # Initialize rand and length vector
      rand <- vector("numeric", length = max(new.vec, na.rm = TRUE))
      names(rand) <- na.omit(new.uniq)
      len <- rand
      
      for(j in new.uniq)
      {
        # Target nodes
        target <- which(new.vec==j)
        
        # Lengths of target
        len[paste(j)] <- length(target)
        
        # Compute rand index
        if(length(target) >= 3){
          rand[paste(j)] <- igraph::compare(new.vec[target],target.wc[target],method="rand")
        }else{
          # Comparisons give error when variables are less than three
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }
    
      }
      
      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))
      
      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE)]
      
      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)
      
      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord)))
      {
        # Identify target
        new.target <- which(new.vec==j)
        
        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)
        
        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }
      
    }else if(max(target.wc, na.rm = TRUE) < max(new.vec, na.rm = TRUE)){
      
      # Initialize rand and length vector
      rand <- vector("numeric", length = max(new.vec, na.rm = TRUE))
      names(rand) <- na.omit(new.uniq)
      len <- rand
      
      for(j in new.uniq){
        
        # Target nodes
        target <- which(new.vec==j)
        
        # Lengths of target
        len[paste(j)] <- length(target)
        
        # Compute rand index
        if(length(target) >= 3){
          rand[paste(j)] <- igraph::compare(new.vec[target],target.wc[target],method="rand")
        }else{
          # Comparisons give error when variables are less than three
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }
      
      }
      
      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))
      
      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE)]
      
      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)
      
      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord))){
        
        # Identify target
        new.target <- which(new.vec==j)
        
        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)
        
        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }
      
      # Identify number of extra dimensions
      extra.dim <- unique(new.vec[which(is.na(final.vec))])
      
      # Initialize extra dimension length vector
      extra.len <- rep(NA, length = length(extra.dim))
      names(extra.len) <- extra.dim
      
      # Initialize count
      count <- 0
      
      # Order length of extra dimensions
      for(j in extra.dim){
        
        # Increase count
        count <- count + 1
        
        # Length of extra dimensions
        extra.len[count] <- length(which(new.vec==j))
      }
      
      el.ord <- extra.len[order(extra.len, decreasing = TRUE)]
      
      # Reset count
      count <- 0
      
      # Insert extra dimensions into final vector
      for(j in 1:length(el.ord)){
        
        # Increase count
        count <- count + 1
        
        # Target extra dimension
        target.ed <- as.numeric(names(el.ord)[j])
        
        # Insert dimensions into final vector
        final.vec[which(new.vec==target.ed)] <- (max(target.wc) + count)
      }
      
    }else{
      
      # Initialize rand and length vector
      rand <- vector("numeric", length = max(new.vec, na.rm = TRUE))
      names(rand) <- na.omit(new.uniq)
      len <- rand
      
      for(j in new.uniq)
      {
        # Target nodes
        target <- which(new.vec==j)
        
        # Lengths of target
        len[paste(j)] <- length(target)
        
        # Compute rand index
        if(length(target) >= 3){
          rand[paste(j)] <- igraph::compare(new.vec[target],target.wc[target],method="rand")
        }else{
          # Comparisons give error when variables are less than three
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }
        
      }
      
      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))
      
      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE, na.last = FALSE)]
      
      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)
      
      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord)))
      {
        # Identify target
        new.target <- which(new.vec==j)
        
        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)
        
        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }
    }
    
    # Insert final vector into final matrix
    convert.mat[,i] <- final.vec
  }
  
  # Add row names
  row.names(convert.mat) <- names(target.wc)
  
  return(convert.mat)
}

#' @noRd
# Checks for missing dimensions
# For itemStability
# Updated 26.02.2021
missing.dimension.check <- function (proportion, membership, bootstrap)
{
  # Get names of proportion
  proportion.names <- as.numeric(colnames(proportion))
  
  # Check for missing dimensions
  missing.dimensions <- setdiff(membership, proportion.names)
  
  # If there are missing dimensions, then add them (all zeroes)
  if(length(missing.dimensions) != 0){
    
    for(i in 1:length(missing.dimensions)){
      
      # Append missing dimension
      proportion <- cbind(proportion, rep(0, nrow(proportion)))
      colnames(proportion)[ncol(proportion)] <- paste(missing.dimensions[i])
      
    }
    
  }
  
  # Check for NA in bootstrap membership
  if(any(is.na(bootstrap))){
    
    # Get NA proportions
    NA.proportions <- colMeans(apply(bootstrap, 1, is.na))
    
    # Append proportion
    proportion <- cbind(proportion, NA.proportions)
    colnames(proportion)[ncol(proportion)] <- "NA"
    
  }
  
  return(proportion)
  
}

#' @noRd
# Plot configuration for itemStability
# For itemStability
# Updated 17.07.2021
itemStability.plot <- function (res, bootega.obj)
{
  # Obtain empirical membership
  empirical.membership <- res$membership$empirical
  
  # Obtain unique membership
  unique.membership <- res$membership$unique
  
  # Obtain empirical stability
  empirical.stability <- res$item.stability$empirical.dimensions
  
  # Organize plot
  organize.plot <- data.frame(Item = names(empirical.membership),
                              Replication = empirical.stability,
                              Community = factor(empirical.membership, 
                                                 unique.membership[order(
                                                   unique.membership
                                                 )]))
  
  # Item stability plot
  IS.plot <- ggpubr::ggdotchart(organize.plot, x = "Item", y = "Replication",
                                group = "Community", color = "Community",
                                legend.title = "Empirical EGA Communities",
                                add = "segments",
                                rotate = TRUE,
                                dot.size = 6,
                                label = round(organize.plot$Replication, 2),
                                font.label = list(color = "black", size = 8,
                                                  vjust = 0.5),
                                ggtheme = ggpubr::theme_pubr()
  )
  
  # Adjust y-axis and legend title
  IS.plot <- IS.plot + ggplot2::ylim(c(0,1)) + ggplot2::theme(
    legend.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold")
  )
  
  # Manually change alpha
  IS.plot$layers[[2]]$aes_params$alpha <- 0.7
  
  # Adjust item label sizes based on
  sizes <- seq(6,12,.25)
  ## Number of nodes
  nodes <- rev(seq(0, 200, length.out = length(sizes)))
  n.size <- min(which(length(empirical.membership) > nodes))
  ## Number of characters in item name
  chars <- rev(seq(0,100, length.out = length(sizes)))
  ### Maximum characters in item name
  max.chars <- max(unlist(lapply(row.names(organize.plot),nchar)))
  c.size <- min(which(max.chars > chars))
  # Text size
  text.size <- sizes[min(c(n.size,c.size))]
  
  # Change text size
  IS.plot <- IS.plot + 
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=text.size))
  
  # Change color.palette (if necessary)
  if(!ggplot2::is.ggplot(bootega.obj$plot.typical.ega)){
    
    if(!is.null(bootega.obj$color.palette)){
      
      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA(bootega.obj$color.palette,
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership))
      )
      
    }else{
      
      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA("rainbow",
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership))
      )
      
    }
    
  }else{
    if(bootega.obj$color.palette != "Set1"){
      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA(bootega.obj$color.palette,
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership, na.last = TRUE)
                                              )
      )
    }
  }
  
  # Reverse ordering
  IS.plot <- IS.plot + ggplot2::scale_x_discrete(limits = rev(IS.plot$data$Item))
  
  # Insert plot into results
  res$plot <- IS.plot
  
  return(res)
}

#' @noRd
# Average network loadings for itemStability
# For itemStability
# Updated 14.04.2021
itemStability.loadings <- function(res, bootega.obj)
{
  # Get graphs
  graphs <- bootega.obj$bootGraphs
  
  # Get memberships
  memberships <- res$membership$bootstrap
  
  # Maximum number of dimensions
  max.dimensions <- max(memberships, na.rm = TRUE)
  
  # Get number of iterations
  iterations <- ncol(memberships)
  
  # Make list of indices
  index.list <- as.list(1:iterations)
  
  # Loop through (for loop is slightly faster than lapply)
  loadings <- list()
  
  for(i in 1:iterations){
    
    if(all(is.na(memberships[,i]))){
      loadings[[i]] <- NULL
    }else{
      loadings[[i]] <- net.loads(A = graphs[[i]],
                                 wc = memberships[,i])$std
    }
    
  }
  
  # Remove NULL loadings
  if(any(unlist(lapply(loadings, is.null)))){
    loadings <- loadings[-which(unlist(lapply(loadings, is.null)))]
  }
  
  # Initialize final loadings array
  loadings.array <- array(NA,
                          dim = c(nrow(memberships), # number of items
                                  max(memberships, na.rm = TRUE), # number of dimensions
                                  iterations), # number of bootstrap replicates
                          dimnames = list(row.names(memberships),
                                          1:max.dimensions,
                                          NULL))
  
  # Loop through
  for(i in 1:length(loadings)){
    
    loadings.array[row.names(loadings[[i]]), # get available loadings
                   colnames(loadings[[i]]), # get available dimensions
                   i] <- as.matrix(loadings[[i]]) # insert loadings into array
    
  }
  
  # Obtain average loadings
  mean.loadings <- apply(loadings.array, 1:2, mean, na.rm = TRUE)
  
  return(mean.loadings)
}



#%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' Error report
#' 
#' @description Gives necessary information for user reporting error
#' 
#' @param result Character.
#' The error from the result
#' 
#' @param SUB_FUN Character.
#' Sub-routine the error occurred in
#' 
#' @param FUN Character.
#' Main function the error occurred in
#' 
#' @return Error and message to send to GitHub
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @noRd
#' 
#' @importFrom utils packageVersion
#' 
# Error Report
# Updated 26.02.2021
error.report <- function(result, SUB_FUN, FUN)
{
  # Let user know that an error has occurred
  message(paste("\nAn error has occurred in the '", SUB_FUN, "' function of '", FUN, "':\n", sep =""))
  
  # Give them the error to send to you
  cat(paste(result))
  
  # Tell them where to send it
  message("\nPlease open a new issue on GitHub (bug report): https://github.com/hfgolino/EGAnet/issues/new/choose")
  
  # Give them information to fill out the issue
  OS <- as.character(Sys.info()["sysname"])
  OSversion <- paste(as.character(Sys.info()[c("release", "version")]), collapse = " ")
  Rversion <- paste(R.version$major, R.version$minor, sep = ".")
  EGAversion <- paste(unlist(packageVersion("EGAnet")), collapse = ".")
  
  # Let them know to provide this information
  message(paste("\nBe sure to provide the following information:\n"))
  
  # To reproduce
  message(styletext("To Reproduce:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " Function error occurred in: ", SUB_FUN, " function of ", FUN, sep = ""))
  
  # R, SemNetCleaner, and SemNetDictionaries
  message(styletext("\nR and EGAnet versions:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " R version: ", Rversion, sep = ""))
  message(paste(" ", textsymbol("bullet"), " EGAnet version: ", EGAversion, sep = ""))
  
  # Desktop
  message(styletext("\nOperating System:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " OS: ", OS, sep = ""))
  message(paste(" ", textsymbol("bullet"), " Version: ", OSversion, sep = ""))
}

#' Read in Common Data File Extensions (from \code{\link{SemNetCleaner}})
#'
#' @description A single function to read in common data file extensions.
#' Note that this function is specialized for reading in text data in the
#' format necessary for functions in SemNetCleaner
#'
#' File extensions supported:
#' \itemize{
#' \item{.Rdata} \item{.rds} \item{.csv} \item{.xlsx}
#' \item{.xls} \item{.sav} \item{.txt} \item{.mat}
#' }
#'
#' @param file Character.
#' A path to the file to load.
#' Defaults to interactive file selection using \code{\link{file.choose}}
#'
#' @param header Boolean.
#' A logical value indicating whether the file contains the
#' names of the variables as its first line.
#' If missing, the value is determined from the file format:
#' header is set to \code{TRUE} if and only if the first row
#' contains one fewer field than the number of columns
#'
#' @param sep Character.
#' The field separator character.
#' Values on each line of the file are separated by this character.
#' If sep = "" (the default for \code{\link{read.table}}) the separator
#' is a 'white space', that is one or more spaces, tabs, newlines or
#' carriage returns
#'
#' @param ... Additional arguments.
#' Allows for additional arguments to be passed onto
#' the respective read functions. See documentation in the list below:
#'
#' \itemize{
#' \item{.Rdata}
#' {\code{\link{load}}}
#' \item{.rds}
#' {\code{\link{readRDS}}}
#' \item{.csv}
#' {\code{\link[utils]{read.table}}}
#' \item{.xlsx}
#' {\code{\link[readxl]{read_excel}}}
#' \item{.xls}
#' {\code{\link[readxl]{read_excel}}}
#' \item{.sav}
#' {\code{\link[foreign]{read.spss}}}
#' \item{.txt}
#' {\code{\link[utils]{read.table}}}
#' \item{.mat}
#' {\code{\link[R.matlab]{readMat}}}
#' }
#'
#' @return A data frame containing a representation of the data in the file.
#' If file extension is ".Rdata", then data will be read to the global environment
#'
#' @examples
#' # Use this example for your data
#' if(interactive())
#' {read.data()}
#'
#' # Example for CRAN tests
#' ## Create test data
#' test1 <- c(1:5, "6,7", "8,9,10")
#'
#' ## Path to temporary file
#' tf <- tempfile()
#'
#' ## Create test file
#' writeLines(test1, tf)
#'
#' ## Read in data
#' read.data(tf)
#'
#' # See documentation of respective R functions for specific examples
#'
#' @references
#' # R Core Team
#'
#' R Core Team (2019). R: A language and environment for
#' statistical computing. R Foundation for Statistical Computing,
#' Vienna, Austria. URL https://www.R-project.org/.
#'
#' # readxl
#'
#' Hadley Wickham and Jennifer Bryan (2019). readxl: Read Excel
#' Files. R package version 1.3.1.
#' https://CRAN.R-project.org/package=readxl
#'
#' # R.matlab
#'
#' Henrik Bengtsson (2018). R.matlab: Read and Write MAT Files
#' and Call MATLAB from Within R. R package version 3.6.2.
#' https://CRAN.R-project.org/package=R.matlab
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom utils read.table read.csv
#' @importFrom tools file_ext
#'
#' @noRd
# Read data
# Updated 15.04.2020
## FROM SemNeT package (utils-SemNeT.R)
read.data <- function (file = file.choose(), header = TRUE, sep = ",", ...)
{
  # Grab extension
  ext <- tolower(file_ext(file))
  
  # Report error
  if(!ext %in% c("rdata", "rds", "csv", "xlsx",
                 "xls", "sav", "txt", "mat", ""))
  {stop("File extension not supported")}
  
  # Determine data load
  if(ext != "")
  {
    switch(ext,
           rdata = load(file, envir = .GlobalEnv),
           rds = readRDS(file),
           csv = read.csv(file, header = header, sep = sep, as.is = TRUE, ...),
           xlsx = as.data.frame(readxl::read_xlsx(file, col_names = header, ...)),
           xls = as.data.frame(readxl::read_xls(file, col_names = header, ...)),
           sav = foreign::read.spss(file, to.data.frame = TRUE, stringAsFactors = FALSE, ...),
           txt = read.table(file, header = header, sep = sep, ...),
           mat = as.data.frame(R.matlab::readMat(file, ...))
    )
  }else{read.table(file, header = header, sep = sep, ...)}
}

#' System check for OS and RSTUDIO
#'
#' @description Checks for whether text options are available
#'
#' @param ... Additional arguments
#'
#' @return \code{TRUE} if text options are available and \code{FALSE} if not
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# System Check
# Updated 08.09.2020
system.check <- function (...)
{
  OS <- unname(tolower(Sys.info()["sysname"]))
  
  RSTUDIO <- ifelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)
  
  TEXT <- TRUE
  
  if(!RSTUDIO){if(OS != "linux"){TEXT <- FALSE}}
  
  res <- list()
  
  res$OS <- OS
  res$RSTUDIO <- RSTUDIO
  res$TEXT <- TEXT
  
  return(res)
}

#' Colorfies Text
#'
#' Makes text a wide range of colors (8-bit color codes)
#'
#' @param text Character.
#' Text to color
#'
#' @return Colorfied text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Color text
# Updated 08.09.2020
colortext <- function(text, number = NULL, defaults = NULL)
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    # Defaults for number (white text)
    if(is.null(number) || number < 0 || number > 231)
    {number <- 15}
    
    # Check for default color
    if(!is.null(defaults))
    {
      # Adjust highlight color based on background color
      if(defaults == "highlight")
      {
        if(sys.check$RSTUDIO)
        {
          
          if(rstudioapi::getThemeInfo()$dark)
          {number <- 226
          }else{number <- 208}
          
        }else{number <- 208}
      }else{
        
        number <- switch(defaults,
                         message = 204,
                         red = 9,
                         orange = 208,
                         yellow = 11,
                         "light green" = 10,
                         green = 34,
                         cyan = 14,
                         blue = 12,
                         magenta = 13,
                         pink = 211,
        )
        
      }
      
    }
    
    return(paste("\033[38;5;", number, "m", text, "\033[0m", sep = ""))
    
  }else{return(text)}
}

#' Stylizes Text
#'
#' Makes text bold, italics, underlined, and strikethrough
#'
#' @param text Character.
#' Text to stylized
#'
#' @return Sytlized text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Style text
# Updated 08.09.2020
styletext <- function(text, defaults = c("bold", "italics", "highlight",
                                         "underline", "strikethrough"))
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    if(missing(defaults))
    {number <- 0
    }else{
      
      # Get number code
      number <- switch(defaults,
                       bold = 1,
                       italics = 3,
                       underline = 4,
                       highlight = 7,
                       strikethrough = 9
      )
      
    }
    
    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))
  }else{return(text)}
}

#' Text Symbols
#'
#' Makes text symbols (star, checkmark, square root)
#'
#' @param symbol Character.
#'
#' @return Outputs symbol
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Symbols
# Updated 24.04.2020
textsymbol <- function(symbol = c("alpha", "beta", "chi", "delta",
                                  "eta", "gamma", "lambda", "omega",
                                  "phi", "pi", "rho", "sigma", "tau",
                                  "theta", "square root", "infinity",
                                  "check mark", "x", "bullet")
)
{
  # Get number code
  sym <- switch(symbol,
                alpha = "\u03B1",
                beta = "\u03B2",
                chi = "\u03C7",
                delta = "\u03B4",
                eta = "\u03B7",
                gamma = "\u03B3",
                lambda = "\u03BB,",
                omega = "\u03C9",
                phi = "\u03C6",
                pi = "\u03C0",
                rho = "\u03C1",
                sigma = "\u03C3",
                tau = "\u03C4",
                theta = "\u03B8",
                "square root" = "\u221A",
                infinity = "\u221E",
                "check mark" = "\u2713",
                x = "\u2717",
                bullet = "\u2022"
  )
  
  return(sym)
}

