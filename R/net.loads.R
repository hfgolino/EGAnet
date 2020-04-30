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
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param min.load Numeric.
#' Sets the minimum loading allowd in the standardized
#' network loading matrix. Values equal or greater than
#' the minimum loading are kept in the output. Values
#' less than the minimum loading are removed. This matrix can
#' be viewed using \code{print()} or \code{summary()}
#' Defaults to \code{0}
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
#' \item{minLoad}{The minimum loading to appear in summary of network loadings.
#' Use \code{print()} or \code{summary()} to view}
#' 
#' \item{plot}{A \code{\link[qgraph]{qgraph}} plot of the network loadings.
#' Use \code{plot} to view} 
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
#' net.loads(ega.wmt)
#'
#' @references
#' Christensen, A. P., & Golino, H. (2020).
#' Statistical equivalency of factor and network loadings.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/xakez}{10.31234/osf.io/xakez}
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. (2019).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/ktejp}{10.31234/osf.io/ktejp}
#'
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. M. (in press).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}.
#' doi:\href{https://doi.org/10.31234/osf.io/pg4mf}{10.31234/osf.io/pg4mf}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings
# Updated 27.04.2020
net.loads <- function(A, wc, pos.manifold = FALSE, min.load = 0, plot = FALSE)
{
  ###############################
  #### START DATA MANAGEMENT ####
  ###############################
  
  #------------------------------------------#
  ## DETECT EGA INPUT AND VARIABLE ORDERING ##
  #------------------------------------------#
  
  if(any(class(A) == "EGA"))
  {
    # Order
    ord <- match(A$dim.variables$items, names(A$wc))
    
    # Grab communities
    wc <- A$wc
    
    # Replace 'A' with 'EGA' network
    A <- A$network
  }else{ord <- order(wc)} # Reorder by communities
  
  #----------------------#
  ## REORDER FOR OUTPUT ##
  #----------------------#
  
  # Reorder communities
  wc <- wc[ord]
  
  # Reorder network
  A <- A[ord,ord]
  
  #---------------------------#
  ## DIMENSION QUALITY CHECK ##
  #---------------------------#
  
  # Remove NA dimensions
  dim.uniq <- na.omit(unique(wc))
  
  # Check for single item dimensions
  for(i in 1:length(dim.uniq))
  {
    len <- length(wc[which(wc==dim.uniq[i])])
    
    if(len == 1)
    {wc[which(wc==dim.uniq[i])] <- NA}
  }
  
  # Remove single item dimensions
  dims <- na.omit(unique(wc))
  
  # Remove NA attribute
  attr(dims, "na.action") <- NULL
  
  #############################
  #### END DATA MANAGEMENT ####
  #############################
  
  # Make sure that there are actual dimensions
  if(all(!is.na(wc)) && length(unique(wc)) != 1)
  {
    
    ################################
    #### START COMPUTE LOADINGS ####
    ################################
    
    # Compute aboslute loadings
    comm.str <- mat.func(A = A, wc = wc, absolute = TRUE, diagonal = 0)
    
    # Check for missing dimensions
    if(any(colnames(comm.str)=="NA"))
    {
      # Target dimension
      target <- which(colnames(comm.str) == "NA")
      
      # Remove from matrix
      comm.str <- comm.str[,-target]
    }
    
    # Reorder loading matrix
    comm.str <- comm.str[,paste(dims)]
    
    # Add signs to loadings
    res.rev <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
    comm.str <- res.rev$comm.str
    A <- res.rev$A
    
    ##############################
    #### END COMPUTE LOADINGS ####
    ##############################
    
    #################################
    #### START OUTPUT MANAGEMENT ####
    #################################
    
    # Initialize result list
    res <- list()
    
    # Unstandardized loadings
    unstd <- as.data.frame(round(comm.str,3))
    row.names(unstd) <- colnames(A)
    res$unstd <- unstd
    
    # Standardized loadings
    if(length(dims)!=1)
    {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
    }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
    res$std <- as.data.frame(round(std,3))
    
    #####################
    #### PLOT SET UP ####
    #####################
    
    #Set to absolute for multidimensional
    std.res <- as.matrix(abs(res$std))
    
    #Standardize by maximum rspbc
    std.res <- std.res / rowSums(std.res)
    
    #Ensure that pie value is not greater than 1
    std.res <- std.res - .001
    std.res <- ifelse(std.res==-.001,0,std.res)
    
    #Split results to list for each node
    pies <- split(std.res, rep(1:nrow(std.res)))
    
    # Plot (or not)
    nl.plot <- qgraph::qgraph(A, layout = "spring", groups = as.factor(wc),
                              label.prop = 1.5, pie = pies, vTrans = 200,
                              negDashed = TRUE, DoNotPlot = ifelse(plot,FALSE,TRUE))
    
    # Remove loadings (added as attribute)
    ## S3Methods summary and print
    res$MinLoad <- min.load
    ## S3Methods plot
    res$plot <- nl.plot
    
  }else if(all(is.na(wc)))
  {
    # Create matrix of NAs
    comm.str <- matrix(NA, nrow = ncol(A), ncol = ncol(A))
    
    # Set up dimensions for all 
    dims <- 1:ncol(A)
    
    # Assign column names
    colnames(comm.str) <- dims
    
    # Set up return
    res <- list()
    res$std <- comm.str
    res$unstd <- comm.str
    
  }else{ # One dimension
    
    # Create matrix of NAs
    comm.str <- matrix(NetworkToolbox::strength(A, absolute = TRUE), nrow = ncol(A), ncol = 1)
    
    # Set up dimensions for all 
    dims <- 1
    
    # Assign column names
    colnames(comm.str) <- dims
    
    # Add signs to loadings
    res.rev <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
    comm.str <- res.rev$comm.str
    A <- res.rev$A
    
    # Initialize result list
    res <- list()
    
    # Unstandardized loadings
    unstd <- as.data.frame(round(comm.str,3))
    row.names(unstd) <- colnames(A)
    res$unstd <- unstd
    
    # Standardized loadings
    if(length(dims)!=1)
    {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
    }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
    res$std <- as.data.frame(round(std,3))
    
  }
  
  class(res) <- "NetLoads"
  
  return(res)
}
#----
