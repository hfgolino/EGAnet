#' Network Loadings
#'
#' @description Computes the between- and within-community
#' \code{strength} of each item
#' for each community. This function uses the
#' \code{comcat} and
#' \code{stable} functions to calculate
#' the between- and within-community strength of each item, respectively.
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' A network adjacency matrix
#'
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param min.load Numeric.
#' Sets the minimum loading allowed in the standardized
#' network loading matrix. Values equal or greater than
#' the minimum loading are kept in the output. Values
#' less than the minimum loading are removed. This matrix can
#' be viewed using \code{print()} or \code{summary()}
#' Defaults to \code{0}
#' 
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure. For a list of rotations,
#' see \link{GPArotation}
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
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#'
#' # Network loadings
#' net.loads(ega.wmt)
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   ega.wmt,
#'   stats = "net.loads"
#' )}
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#'
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. M. (2019).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}, 1-25.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings
# Updated 09.11.2022
# Signs updated 18.10.2022
# Rotations added 20.10.2022
net.loads <- function(
    A, wc, rotation = "oblimin",
    min.load = 0,
    ...
)
{
  
  # Deprecated arguments
  if("pos.manifold" %in% names(list(...))){
    message("Argument 'pos.manifold' has been deprecated.")
  }
  
  #------------------------------------------#
  ## DETECT EGA INPUT AND VARIABLE ORDERING ##
  #------------------------------------------#
  
  if(any(class(A) == "EGA")){
    
    # Order
    ord <- match(colnames(A$network), names(A$wc))
    
    # Grab communities
    wc <- A$wc
    
    # Replace 'A' with 'EGA' network
    A <- A$network
    
  }else{ord <- order(wc)} # Reorder by communities
  
  # Make sure membership is named
  names(wc) <- colnames(A)
  
  # Ensure matrix object
  A <- as.matrix(A)
  
  # Ensure data is matrix
  if(nrow(A) != ncol(A)){
    stop("Input for 'A' must be an n x n matrix.")
  }
  
  # Ensure data is symmetric
  row.names(A) <- colnames(A)
  
  if(!isSymmetric(A)){
    stop("'A' is not a symmetric matrix. Network loadings can only be computed with undirected networks.")
  }
  
  # Check if there are actual dimensions
  if(length(wc) == length(unique(wc))){
    
    # Initialize result list
    res <- list()
    
    res$unstd <- matrix(NA, nrow = ncol(A), ncol = ncol(A))
    
    res$std <- matrix(NA, nrow = ncol(A), ncol = ncol(A))
    
  }else{
  
    #### START DATA MANAGEMENT
    
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
    for(i in 1:length(dim.uniq)){
      
      len <- length(wc[which(wc==dim.uniq[i])])
      
      if(len == 1)
      {wc[which(wc==dim.uniq[i])] <- NA}
    }
    
    # Remove single item dimensions
    dims <- na.omit(unique(wc))
    
    # Remove NA attribute
    attr(dims, "na.action") <- NULL
    
    # Make sure that there are actual dimensions
    if(length(dims) != 1)
    {

      #### START COMPUTE LOADINGS
      
      # Compute absolute loadings
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
      # res.rev <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
      # comm.str <- res.rev$comm.str
      comm.str <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
      
      #### START OUTPUT MANAGEMENT
      
      # Initialize result list
      res <- list()
      
      # Unstandardized loadings
      unstd <- as.data.frame(round(comm.str,3))
      row.names(unstd) <- colnames(A)
      res$unstd <- descend.ord(unstd, wc)
      
      # Standardized loadings
      if(length(dims)!=1)
      {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
      }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
      res$std <- as.data.frame(round(descend.ord(std, wc),3))
      
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
      comm.str <- matrix(strength(A, absolute = TRUE), nrow = ncol(A), ncol = 1)
      
      # Assign column names
      colnames(comm.str) <- dims
      
      # Add signs to loadings
      comm.str <- add.signs(comm.str = comm.str, A = A, wc = wc, dims = dims, pos.manifold = pos.manifold)
      # comm.str <- res.rev$comm.str
      
      # Initialize result list
      res <- list()
      
      # Unstandardized loadings
      unstd <- as.data.frame(round(comm.str,3))
      row.names(unstd) <- colnames(A)
      res$unstd <- descend.ord(unstd, wc)
      
      # Standardized loadings
      if(length(dims)!=1)
      {std <- t(t(unstd) / sqrt(colSums(abs(unstd))))
      }else{std <- t(t(unstd) / sqrt(sum(abs(unstd))))}
      res$std <- as.data.frame(round(descend.ord(std, wc),3))
      
    }
    
    class(res) <- "NetLoads"
    
  }
  
  # Check for more than one dimension
  if(ncol(std) > 1){
    
    # Check for {GPArotation}
    check_package("GPArotation")
    
    # Obtain rotation from GPArotation package
    rotation_names <- ls(asNamespace("GPArotation"))
    
    # Check if rotation exists
    rotation <- tolower(rotation)
    rotation_names_lower <- tolower(rotation_names)
    
    # Obtain rotation arguments
    rot_arguments <- list(...)
    
    if(rotation %in% rotation_names_lower){
      
      # Obtain rotation function arguments
      # rotation_function <- get(
      #   rotation_names[which(rotation == rotation_names_lower)],
      #   envir = asNamespace("GPArotation")
      # )
      
      # Obtain arguments
      # rotation_arguments <- obtain.arguments(
      #   FUN = rotation_function, FUN.args = list(...)
      # )
      
      if(rotation != "oblimin"){
        
        # Obtain arguments
        rotation_arguments <- obtain.arguments(
          FUN = psych::faRotations, FUN.args = list(rotate = rotation)
        )
        
        # Check for arguments
        rotation_arguments$loadings <- std
        rotation_arguments$n.rotations <- ifelse(
          "n.rotations" %in% names(rot_arguments),
          rot_arguments$n.rotations,
          10
        )
        rotation_arguments$maxit <- ifelse(
          "maxit" %in% names(rot_arguments),
          rot_arguments$maxit,
          1000
        )
        
        # Set loadings
        # rotation_arguments$L <- as.matrix(res$std)
        # rotation_arguments$Tmat <- diag(ncol(rotation_arguments$L))
        rotation_arguments$loadings <- as.matrix(res$std)
        
        # Obtain rotated loadings
        # res$rotated <- do.call(
        #   what = rotation_function,
        #   args = as.list(rotation_arguments)
        # )
        res$rotated <- do.call(
          what = psych::faRotations,
          args = as.list(rotation_arguments)
        )
        
      }else{
        
        # Obtain arguments
        rotation_arguments <- obtain.arguments(
          FUN = oblimin_rotate, FUN.args = list(...)
        )
        
        # Check for arguments
        rotation_arguments$loadings <- std
        rotation_arguments$n.rotations <- ifelse(
          "n.rotations" %in% names(rot_arguments),
          rot_arguments$n.rotations,
          10
        )
        rotation_arguments$maxit <- ifelse(
          "maxit" %in% names(rot_arguments),
          rot_arguments$maxit,
          1000
        )
        
        # Set loadings
        # rotation_arguments$L <- as.matrix(res$std)
        # rotation_arguments$Tmat <- diag(ncol(rotation_arguments$L))
        rotation_arguments$loadings <- as.matrix(res$std)
        
        # Obtain rotated loadings
        # res$rotated <- do.call(
        #   what = rotation_function,
        #   args = as.list(rotation_arguments)
        # )
        res$rotated <- do.call(
          what = oblimin_rotate,
          args = as.list(rotation_arguments)
        )
        
      }
      
    }
    
  }
  
  
  # Add minimum loading
  res$minLoad <- min.load
  
  return(res)
}
#----
