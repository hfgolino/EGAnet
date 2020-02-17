#' Network Loadings
#'
#' @description Computes the between- and within-community
#' \code{\link[NetworkToolbox]{strength}} of each item
#' for each community. This function uses the
#' \code{\link[NetworkToolbox]{comcat}} and
#' \code{\link[NetworkToolbox]{stable}} functions to calculate
#' the between- and within-community strength of each item, respectively.
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' An adjacency matrix of network data
#'
#' @param wc Numeric.
#' A vector of community assignments.
#' Not necessary if an \code{\link[EGAnet]{EGA}} object
#' is input for argument \code{A}
#'
#' @param rm.zero Should zeros be removed from the resulting matrix?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to reduce the noise in the results
#' 
#' @param pos.manifold Boolean.
#' Should a positive manifold be applied (i.e., should
#' all dimensions be positively correlated)?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for a positive manifold
#'
#' @param plot Boolean.
#' Should proportional loadings be plotted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for plot with pie charts
#' visualizing the proportion of loading associated with
#' each dimension
#'
#' @return Returns a list containing:
#'
#' \item{unstd}{A matrix of the unstandardized within- and between-community
#' strength values for each node}
#'
#' \item{std}{A matrix of the standardized within- and between-community
#' strength values for each node}
#'
#' @details Simulation studies have demonstrated that a node's strength
#' centrality is roughly equivalent to factor loadings
#' (Christensen, Golino, & Silvia, 2019; Hallquist, Wright, & Molenaar, in press).
#' Hallquist and colleagues (in press) found that node strength represented a
#' combination of dominant and cross-factor loadings. This function computes
#' each node's strength within each specified dimension, providing a rough
#' equivalent to factor loadings (including cross-loadings).
#'
#' For more details, type \code{vignette("Network_Scores")}
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(wmt)
#'
#' }
#'
#' # Network loadings
#' net.loads(ega.wmt, rm.zero = TRUE)
#'
#' @references
#' Christensen, A. P., Golino, H. F., & Silvia, P. (2019).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/ktejp}{10.31234/osf.io/ktejp}
#'
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. (in press).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}.
#' doi:\href{https://doi.org/10.31234/osf.io/pg4mf}{10.31234/osf.io/pg4mf}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings
# Updated 17.02.2020
net.loads <- function(A, wc, pos.manifold = FALSE, rm.zero = FALSE, plot = FALSE)
{
  # Add signs to loadings function
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
  
  # Detect if input is an 'EGA' object
  if(any(class(A) == "EGA"))
  {
    # Order
    ord <- match(A$dim.variables$items, names(A$wc))
    
    # Grab communities
    wc <- A$wc[ord]
    
    # Replace 'A' with 'EGA' network
    A <- A$network[ord,ord]
  }
  
  # Dimensions
  ## Check for single item dimensions
  dim.uniq <- na.omit(unique(wc))
  
  for(i in 1:length(dim.uniq))
  {
    len <- length(wc[which(wc==dim.uniq[i])])
    
    if(len == 1)
    {wc[which(wc==dim.uniq[i])] <- NA}
  }
  
  dims <- sort(na.omit(unique(wc)))
  attr(dims, "na.action") <- NULL
  
  mat.func <- function(A, wc, metric = "each", absolute, diagonal)
  {
    comc <- NetworkToolbox::comcat(A = A, comm = wc, metric = metric,
                                   absolute = absolute, diagonal = diagonal)
    stab <- NetworkToolbox::stable(A = A, comm = wc, absolute = absolute, diagonal = diagonal)
    
    for(q in 1:nrow(comc))
    {comc[q,which(is.na(comc[q,]))] <- stab[q]}
    
    if(ncol(comc)!=1)
    {
      comm.str <- comc[,order(colnames(comc))]
      comm.str <- round(comm.str,3)
    }else{comm.str <- stab}
    
    return(comm.str)
  }
  
  # Compute aboslute loadings
  comm.str <- mat.func(A = A, wc = wc, absolute = TRUE, diagonal = 0)
  
  if(any(colnames(comm.str)=="NA"))
  {comm.str <- comm.str[,-which(colnames(comm.str) == "NA")]}
  
  # Check for reverse signs
  res.rev <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
  comm.str <- res.rev$comm.str
  A <- res.rev$A
  
  #result list
  res <- list()
  
  #unstandardized loadings
  if(rm.zero)
  {
    comm.str <- as.data.frame(ifelse(comm.str==0,"",round(comm.str,3)))
    unstd <- apply(as.matrix(comm.str),2,as.numeric)
    unstd <- ifelse(is.na(unstd),0,unstd)
    row.names(comm.str) <- colnames(A)
    
    if(is.null(colnames(comm.str)))
    {colnames(comm.str) <- 1:ncol(comm.str)}
    
    comm.str <- comm.str[,order(colnames(comm.str))]
    res$unstd <- comm.str
    
    #standardized loadings
    if(length(dims)!=1)
    {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
    }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
    
    std <- round(std,3)
    std <- as.data.frame(ifelse(std==0,"",std))
    
    row.names(std) <- colnames(A)
    
    res$std <- std
    
  }else{
    unstd <- apply(as.matrix(comm.str),2,as.numeric)
    row.names(unstd) <- colnames(A)
    
    if(is.null(colnames(unstd)))
    {colnames(unstd) <- 1:ncol(unstd)}
    
    unstd <- unstd[,order(colnames(unstd))]
    res$unstd <- round(unstd,3)
    
    #standardized loadings
    if(length(dims)!=1)
    {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
    }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
    
    row.names(std) <- colnames(A)
    
    res$std <- round(std,3)
  }
  
  #Plot?
  if(plot)
  {
    #Set to absolute for multidimensional
    std.res <- abs(res$std)
    
    #Standardize by maximum rspbc
    std.res <- std.res / rowSums(std.res)
    
    #Ensure that pie value is not greater than 1
    std.res <- std.res - .001
    std.res <- ifelse(std.res==-.001,0,std.res)
    
    #Split results to list for each node
    pies <- split(std.res, rep(1:nrow(std.res)))
    
    #Plot
    qgraph::qgraph(A, layout = "spring",
                   groups = as.factor(wc),
                   label.prop = 1.5,
                   pie = pies,
                   vTrans = 200,
                   negDashed = TRUE)
  }
  
  if(rm.zero)
  {message("Argument 'rm.zero = TRUE': Output is provided in factors. Set argument 'rm.zero = FALSE' to provide numeric output for summarizing results")}
  
  return(res)
}
#----
