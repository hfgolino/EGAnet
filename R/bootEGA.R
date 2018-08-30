#'  Investigates the stability of EGA's estimation via bootstrap.
#'
#' \code{bootEGA} Estimates the number of dimensions of n bootstraps from the empirical correlation matrix,
#'  and returns a typical network (i.e. the network formed by the median or mean pairwise correlations over the n bootstraps) and its dimensionality.
#'
#' @param data A dataframe with the variables to be used in the analysis
#' @param n An integer value representing the number of bootstraps
#' @param typicalStructure Logical. If true, returns the typical network of partial correlations (estimated via graphical lasso or via TMFG) and estimates its dimensions. The "typical network" is the median of all pairwise correlations over the n bootstraps.
#' @param plot.typicalStructure Logical. If true, returns a plot of the typical network (partial correlations), which is the median of all pairwise correlations over the n bootstraps, and its estimated dimensions.
#' @param model A string indicating the method to use. Current options are:
#' -\code{GGM}:
#' {Gaussian Markov random field estimation using graphical LASSO with extended Bayesian information criterion to select optimal regularization parameter. Using \code{\link[qgraph]{EBICglasso}} from the qgraph package.}
#' \code{TMFG}:
#' {Estimates a Triangulated Maximally Filtered Graph, using the function \code{TMFG} of the NetworkToolbox package}
#' @param type A string indicating the type of bootstrap to use. Current options are:
#' -\code{parametric}:
#' {Generates n new datasets (multivariate normal random distributions) based on the
#' original dataset, via the \code{\link[mvtnorm]{rmvnorm}} function of the mvtnorm package}.
#' -\code{resampling}:
#' {Generates n random subsamples of the original data.}
#' @param ncores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#' @param confirm A vector corresponding to the expected community structure.
#' The number of times this structure appears will be counted.
#' Defaults to NULL
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <apchrist at uncg.edu>
#' @examples
#' \dontrun{
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 500, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "GGM", type = "parametric", ncores = 4)
#' boot.intwl <- bootEGA(data = intelligenceBattery[,8:66], n = 500, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "GGM", type = "parametric", ncores = 4)
#'}
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

# Bootstrap EGA:
bootEGA <- function(data, n, typicalStructure = TRUE, plot.typicalStructure = TRUE, ncores = 4,
                    model = c("GGM", "TMFG"), type = c("parametric", "resampling"),
                    confirm = NULL) {
  if(!require(qgraph)) {
    message("installing the 'qgraph' package")
    install.packages("qgraph")
  }

  if(!require(bootnet)) {
    message("installing the 'bootnet' package")
    install.packages("bootnet")
  }

  if(!require(igraph)) {
    message("installing the 'igraph' package")
    install.packages("igraph")
  }

  if(!require(NetworkToolbox)) {
    message("installing the 'NetworkToolbox' package")
    install.packages("NetworkToolbox")
  }

  if(!require(foreach)) {
    message("installing the 'foreach' package")
    install.packages("foreach")
    library(foreach)
  }

  if(!require(doSNOW)) {
    message("installing the 'doSNOW' package")
    install.packages("doSNOW")
  }
  
  #mode function for item confirm
  mode <- function(v)
  {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  #Parallel processing
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)

  #progress bar
  pb <- txtProgressBar(max=n, style = 3)
  progress <- function(num) setTxtProgressBar(pb, num)
  opts <- list(progress = progress)

  boots <- list()

  #nets
  boots <-foreach::foreach(i=1:n,
                           .packages = c("NetworkToolbox","psych","qgraph"),
                           .options.snow = opts) %dopar%
                           {
                             if(model=="GGM")
                             {
                               if(type=="parametric")  # Use a parametric approach:
                               {
                                 g <- -qgraph::EBICglasso(cov(data), n = nrow(data))
                                 diag(g) <- 1
                                 bootData <- mvtnorm::rmvnorm(nrow(data), sigma = corpcor::pseudoinverse(g))
                                 net <- qgraph::EBICglasso(cov(bootData), n = nrow(data))
                               }else if(type=="resampling") # Random subsample with replace
                               {
                                 mat <- data[sample(1:nrow(data), replace=TRUE),]
                                 net <- qgraph::EBICglasso(cov(mat), n=nrow(data))
                               }

                             }else if(model=="TMFG")
                             {
                               if(type=="parametric"){
                                 g <- -LoGo(data, partial=TRUE)
                                 diag(g) <- 1
                                 bootData <- mvtnorm::rmvnorm(nrow(data), sigma = corpcor::pseudoinverse(g))
                                 net <- TMFG(bootData)$A
                               } else if(type=="resampling"){
                                 mat <- data[sample(1:nrow(data), replace=TRUE),]
                                 net <- TMFG(mat)$A
                               }

                             }
                           }

  parallel::stopCluster(cl)

  bootGraphs <- vector("list", n)
  for (i in 1:n) {
    bootGraphs[[i]] <- boots[[i]]
    colnames(bootGraphs[[i]]) <- colnames(data)
    rownames(bootGraphs[[i]]) <- colnames(data)
  }
  boot.igraph <- vector("list", n)
  for (l in 1:n) {
    boot.igraph[[l]] <- as.igraph(qgraph(abs(bootGraphs[[l]]), DoNotPlot = TRUE))
  }
  boot.wc <- vector("list", n)
  for (m in 1:n) {
    boot.wc[[m]] <- walktrap.community(boot.igraph[[m]])
  }
  boot.ndim <- matrix(NA, nrow = n, ncol = 2)
  #Initiate confirm matrix
  if(!is.null(confirm))
    {
      uniq <- unique(confirm)
      confirm.dim <- matrix(NA, nrow = n, ncol = length(uniq))
      item.confirm <- matrix(NA, nrow = n, ncol = ncol(data))
      dim.nmi <- vector("numeric", length = n)
    
      #check if confirm is character
      if(is.character(confirm))
        {
            uni <- unique(confirm)
            num.comm <- confirm
            
            for(i in 1:length(uni))
            {num.comm[which(num.comm==uniq[i])] <- i}
        } else {num.comm <- confirm} 
    }
  
  for (m in 1:n) {
    boot.ndim[m, 2] <- max(boot.wc[[m]]$membership)
        
    #Check if dimension is confirmed
    if(!is.null(confirm))
      {
        #normalized mutual information of community comparisons
        dim.nmi[m] <- igraph::compare(boot.wc[[m]]$membership, num.comm, method="nmi")

        for(i in 1:length(uniq))
          {
            dim.items <- which(confirm==uniq[i])
            target.dim <- boot.wc[[m]]$membership[dim.items]
            uniq.dim <- unique(target.dim)
            if(length(uniq.dim)==1){confirm.dim[m,i] <- 1}else{confirm.dim[m,i] <- 0}
            
            #Check if item is confirmed within dimension
            if(length(uniq.dim)>1)
              {
                target.mode <- mode(target.dim)
                non.con <- dim.items[which(target.dim!=target.mode)]
                con <- setdiff(dim.items,non.con)
                item.confirm[m,non.con] <- 0
                item.confirm[m,con] <- 1
              } else {item.confirm[m,dim.items] <- 1}
          }
      }
  }
  
  if(!is.null(confirm))
    {
      #Proportion of times dimension is confirmed
      con.dim <- (colSums(confirm.dim))/n
      names(con.dim) <- uniq
    
      #Proportion of times item is confirmed
      con.item <- (colSums(item.confirm)/n)
      names(con.item) <- colnames(data)
    
      #Tables for nmi
      dim.nmi.table <- vector("numeric", length = 2)
      dim.nmi.table[1] <- mean(dim.nmi)
      dim.nmi.table[2] <- sd(dim.nmi)
      names(dim.nmi.table) <- c("mean","sd")
    }
  
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")
  boot.ndim[, 1] <- seq_len(n)
  if (typicalStructure == TRUE) {
    if(model=="GGM")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, median)
    }else if(model=="TMFG")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, mean)}
    typical.igraph <- as.igraph(qgraph(abs(typical.Structure), DoNotPlot = TRUE), attributes = FALSE)
    typical.wc <- walktrap.community(typical.igraph)
    typical.ndim <- max(typical.wc$membership)
    dim.variables <- data.frame(items = colnames(data), dimension = typical.wc$membership)
  }
  if (plot.typicalStructure == TRUE) {
    plot.typical.ega <- qgraph(typical.Structure, layout = "spring",
                               vsize = 6, groups = as.factor(typical.wc$membership))
  }
  Median <- median(boot.ndim[, 2])
  sd.boot <- sd(boot.ndim[, 2])
  se.boot <- (1.253 * sd.boot)/sqrt(nrow(boot.ndim))
  ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
  ci <- se.boot * ciMult
  summary.table <- data.frame(n.Boots = n, median.dim = Median,
                              SD.dim = sd.boot, SE.dim = se.boot, CI.dim = ci, Lower = Median -
                                ci, Upper = Median + ci)

  #compute likelihood
  dim.range <- range(boot.ndim[,2])
  lik <- matrix(0, nrow = diff(dim.range)+1, ncol = 2)
  colnames(lik) <- c("# of Factors", "Likelihood")
  count <- 0

  for(i in seq(from=min(dim.range),to=max(dim.range),by=1))
  {
    count <- count + 1
    lik[count,1] <- i
    lik[count,2] <- length(which(boot.ndim[,2]==i))/n
  }

  result <- list()
  result$n <- n
  result$boot.ndim <- boot.ndim
  result$bootGraphs <- bootGraphs
  result$summary.table <- summary.table
  result$likelihood <- lik
  if(!is.null(confirm))
  {
    result$dim.confirm <- con.dim
    result$item.confirm <- con.item
    result$dim.nmi <- dim.nmi.table
  }
  typicalGraph <- list()
  typicalGraph$graph <- typical.Structure
  typicalGraph$typical.dim.variables <- dim.variables[order(dim.variables[,2]), ]
  result$typicalGraph <- typicalGraph
  class(result) <- "bootEGA"
  return(result)
}
