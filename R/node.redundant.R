#' Detects Redundant Nodes in a Network
#' 
#' @description Identifies redundant nodes in the network based on several
#' measures. Computes the weighted topological overlap between
#' each node and every other node in the network. The weighted topological
#' overlap is implemented using the method from Nowick et al. (2009; see references)
#' and the function \link[wTO]{wTO} from the wTO package.
#'
#' @param data Matrix or data frame.
#' Input can either be data or a correlation matrix
#' 
#' @param n Numeric.
#' If input in \code{data} is a correlation matrix and
#' \code{method = "wTO"}, then sample size is required.
#' Defaults to NULL
#'
#' @param sig Numeric.
#' \emph{p}-value for significance of overlap (defaults to \code{.05}).
#' If more than 200 connections, then \code{\link[fdrtool]{fdrtool}}
#' is used to correct for false positives. In these instances, \code{sig}
#' sets the \emph{q}-value for significance of overlap (defaults to \code{.10})
#'
#' @param method Character.
#' Computes weighted topological overlap (\code{"wTO"} using \code{\link[qgraph]{EBICglasso}}),
#' partial correlations (\code{"pcor"}), and correlations (\code{"cor"})
#' 
#' @param thresh Boolean.
#' Should a threshold be applied?
#' Defaults to \code{FALSE}.
#' If \code{TRUE}, then based on a certain threshold only redundancies
#' above that value will be returned.
#' Uses argument \code{"sig"} to input the desired threshold.
#' Defaults for each method:
#' 
#' \itemize{
#' 
#' \item{\code{"wTO"}}
#' {.20}
#' 
#' \item{\code{"pcor"}}
#' {.20}
#' 
#' \item{\code{"cor"}}
#' {.70}
#' 
#' } 
#'
#' @param type Character.
#' Computes significance using the standard \emph{p}-value (\code{"alpha"}),
#' bonferroni corrected \emph{p}-value (\code{"bonferroni"}),
#' false-discovery rate corrected \emph{p}-value (\code{"FDR"}),
#' or adaptive alpha \emph{p}-value (\code{\link[NetworkToolbox]{adapt.a}}).
#' Defaults to \code{"adapt"}
#' 
#' @param plot Boolean.
#' Should redundancies be plotted in a network plot?
#' Defaults to \code{FALSE}
#'
#' @return Returns a list:
#'
#' \item{redundant}{Vectors nested within the list corresponding
#' to redundant nodes with the name of object in the list}
#'
#' \item{data}{Returns original data}
#'
#' \item{weights}{Returns weights determine by weighted topological overlap
#' or partial correlations}
#'
#' \item{network}{The network compute by \code{\link[qgraph]{EBICglasso}}}
#' 
#' \item{descriptives}{
#' 
#' \itemize{
#' 
#' \item{basic}{A vector containing the mean, standard deviation,
#' median, median absolute deviation (MAD), 3 times the MAD, 6 times the MAD,
#' minimum, maximum, and critical value for the overlap measure
#' (i.e., weighted topological overlap, partial correlation, or threshold)}
#' 
#' \item{centralTendency}{A matrix for all (aboslute) non-zero values and their
#' respective standard deviation from the mean and median absolute deviation
#' from the median}
#' 
#' }
#' }
#' 
#' \item{distribution}{Distribution that was used to determine significance}
#'
#' @examples
#' # obtain SAPA items
#' items <- psychTools::spi[,c(11:20)]
#' 
#' # weighted topological overlap
#' redund <- node.redundant(items, method = "wTO", type = "adapt", plot = TRUE)
#'
#' # partial correlation
#' redund <- node.redundant(items, method = "pcor", type = "adapt", plot = TRUE)
#' 
#' # threshold
#' redund <- node.redundant(items, method = "pcor", thresh = TRUE, sig = .20, plot = TRUE)
#'
#' @references
#' # Simulation using node.redundant \cr
#' Christensen, A. P. (2020).
#' Towards a network psychometrics approach to assessment: Simulations for redundancy, dimensionality, and loadings
#' (Unpublished doctoral dissertation). University of North Carolina at Greensboro, Greensboro, NC, USA.
#' doi: \href{https://doi.org/10.31234/osf.io/84kgd}{10.31234/osf.io/84kgd}
#' 
#' # Implementation of node.redundant \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (in press).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}.
#' doi: \href{https://doi.org/10.1002/per.2265}{10.1002/per.2265}
#' 
#' # wTO measure \cr
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' doi: \href{https://doi.org/10.1073/pnas.0911376106}{10.1073/pnas.0911376106}
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom stats pgamma pnorm qgamma qnorm cov2cor mad
#'
#' @export
#
# Redundant Nodes Function
# Updated 25.11.2020
node.redundant <- function (data, n = NULL, sig, method = c("wTO", "pcor", "cor"),
                            thresh = FALSE, type = c("alpha", "bonferroni", "FDR", "adapt"),
                            plot = FALSE)
{
  #### missing arguments handling ####
  
  if(missing(type)){
    type <- "adapt"
  }else{type <- match.arg(type)}
  
  #### missing arguments handling ####
  
  # change arguments to lower
  type <- tolower(type)
  method <- tolower(method)
  
  # check for correlation matrix
  if(ncol(data) == nrow(data)){
    
    A <- data
    
    if(method == "wto"){# check for number of cases ("wTO" only)
      
      if(is.null(n))
      {stop('Argument \'n\' is NULL. Number of cases must be specified when a correlation matrix is input into the \'data\' argument and method = "wTO"')}
    }
    
  }else{# correlate data
    A <- qgraph::cor_auto(data)
    # number of cases
    n <- nrow(data)
  }

  #compute redundant method
  if(method == "wto"){# compute network
    net <- EBICglasso.qgraph(A, n = n)
    tom <- wTO::wTO(net, sign = "sign")
  }else if(method == "pcor"){# compute precision matrix
    tom <- -cov2cor(solve(A))
  }else{tom <- A}

  #number of nodes
  nodes <- ncol(tom)
  #make diagonal zero
  diag(tom) <- 0
  #absolute values
  tom <- abs(tom)

  #lower triangle of TOM
  lower <- tom[lower.tri(tom)]

  #grab names for lower triangle
  name1 <- colnames(tom)
  name2 <- colnames(tom)

  #initialize name matrix
  name.mat <- tom

  #produce name matrix
  for(i in 1:nodes)
    for(j in 1:nodes)
    {name.mat[i,j] <- paste(name1[j],name2[i],sep="--")}

  #grab lower triangle names
  names(lower) <- name.mat[lower.tri(name.mat)]

  #obtain only positive values
  pos.vals <- na.omit(ifelse(lower==0,NA,lower))
  attr(pos.vals, "na.action") <- NULL
  
  if(thresh){
    
    #get defaults
    if(missing(sig)){
      sig <- switch(method,
                    "wto" = .20,
                    "pcor" = .20,
                    "cor" = .70
                    )
    }
    
    #threshold values
    res <- pos.vals[which(abs(pos.vals) > sig)]
    
  }else{
    
    #determine distribution
    ##distributions
    distr <- c("norm","gamma")
    aic <- vector("numeric", length = length(distr))
    names(aic) <- c("normal","gamma")
    
    for(i in 1:length(distr))
    {aic[i] <- fitdistrplus::fitdist(pos.vals,distr[i],method="mle")$aic}
    
    #obtain distribution parameters
    g.dist <- suppressWarnings(MASS::fitdistr(pos.vals, names(aic)[which.min(aic)]))
    
    #compute significance values
    pval <- switch(names(aic)[which.min(aic)],
                   
                   normal = 1 - unlist(lapply(pos.vals, #positive wTo
                                              pnorm, #probability in normal distribution
                                              mean = g.dist$estimate["mean"], #mean of normal
                                              sd = g.dist$estimate["sd"]) #sd of normal
                   ),
                   
                   gamma = 1 - unlist(lapply(pos.vals, #positive wTo
                                             pgamma, #probability in gamma distribution
                                             shape = g.dist$estimate["shape"], #shape of gamma
                                             rate = g.dist$estimate["rate"]) #rate of gamma
                   ),
    )
    
    #switch for missing arguments
    if(missing(sig))
    {
      sig <- switch(type,
                    fdr = .10,
                    bonferroni = .05,
                    adapt = .05,
                    alpha = .05
      )
    }else{sig <- sig}
    
    #switch to compute pvals
    if(type == "fdr")
    {
      pval <- suppressWarnings(fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE, verbose = FALSE)$qval)
    }else{
      sig <- switch(type,
                    bonferroni = sig / length(pos.vals),
                    adapt = NetworkToolbox::adapt.a("cor", alpha = sig, n = length(pos.vals), efxize = "medium")$adapt.a,
                    alpha = sig
      )
    }
    
    #identify q-values less than significance
    res <- pos.vals[which(pval<=sig)]
    
  }

  #if no redundant, then print no redundant nodes
  if(length(res)==0)
  {
    message(paste("No redundant nodes identified. Increase 'sig' arugment to detect more nodes."))
    res.list <- NA
  }else{

    #create result matrix
    split.res <- unlist(strsplit(names(res), split = "--"))
    res.mat <- t(simplify2array(sapply(names(res), strsplit, split = "--")))

    #initialize result list
    res.list <- list()

    #initialize count
    count <- 0

    while(nrow(res.mat)!=0)
    {
      #increase count
      count <- count + 1

      #get variable counts
      var.counts <- sort(table(split.res), decreasing = TRUE)

      if(!all(var.counts==1))
      {
        #identify targets
        target <- which(res.mat == names(var.counts[1]), arr.ind = TRUE)[,"row"]

        #insert values into list
        res.list[[names(var.counts[1])]] <- setdiff(unlist(strsplit(names(target),split="--")),names(var.counts[1]))

        #remove rows from result matrix
        res.mat <- res.mat[-target,]

        #remove variables from split result
        split.res <- as.vector(res.mat)

        #force matrix
        if(is.vector(res.mat))
        {res.mat <- t(as.matrix(res.mat))}

      }else
      {
        for(i in 1:nrow(res.mat))
        {res.list[[res.mat[i,1]]] <- unname(res.mat[i,2])}

        res.mat <- res.mat[-c(1:nrow(res.mat)),]
      }
    }
  }
  
  # Revert tom to matrix
  tom <- as.matrix(tom)
  
  # Plot
  if(plot)
  {
    # Initialize plot matrix
    plot.mat <- matrix(0, nrow = nrow(tom), ncol = ncol(tom))
    colnames(plot.mat) <- colnames(tom)
    row.names(plot.mat) <- colnames(tom)
    
    for(i in 1:length(res.list))
    {
      plot.mat[names(res.list)[i],res.list[[i]]] <- tom[names(res.list)[i],res.list[[i]]]
      plot.mat[res.list[[i]],names(res.list)[i]] <- tom[res.list[[i]],names(res.list)[i]]
    }
    
    rm.mat <- which(colSums(plot.mat) == 0)
    
    plot.mat <- plot.mat[-rm.mat, -rm.mat]
    
    if(method == "wto")
    {title <- "Weighted Topological Overlaps"
    }else{title <- "Partial Correlations"}
    
    qgraph::qgraph(plot.mat, layout = "spring", title = title,
                   edge.labels = TRUE)
  }
  
  # Initialize descriptives matrix
  desc <- matrix(0, nrow = 1, ncol = 9)
  
  # Row name
  if(method == "wto")
  {row.names(desc) <- "wTO"
  }else{row.names(desc) <- "pcor"}
  
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
  if(thresh)
  {desc[,"Critical Value"] <- sig
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
  
  if(method == "wto")
  {colnames(pos.output)[1] <- "wTO"
  }else{colnames(pos.output)[1] <- "pcor"}
  
  colnames(pos.output)[2:3] <- c("SD from Mean", "MAD from Median")
  
  
  full.res <- list()
  full.res$redundant <- res.list
  full.res$data <- data
  full.res$weights <- tom
  if(exists("net")){full.res$network <- net}
  full.res$descriptives$basic <- round(desc, 3)
  full.res$descriptives$centralTendency <- pos.output
  if(!thresh){full.res$distribution <- names(aic)[which.min(aic)]}
  
  class(full.res) <- "node.redundant"

  return(full.res)
}
#----
