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
#Communicating----
#Updated 18.03.2020
comcat <- function (A, comm = c("walktrap","louvain"),
                    cent = c("strength","degree"),
                    absolute = TRUE,
                    metric = c("across","each"),
                    diagonal = 0, ...)
{
  ###########################
  #### MISSING ARGUMENTS ####
  ###########################
  
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
  
  #######################
  #### MAIN FUNCTION ####
  #######################
  
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
#----

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
####Stabilizing####
#Updated 18.03.2020
stable <- function (A, comm = c("walktrap","louvain"),
                    cent = c("betweenness","rspbc","closeness",
                             "strength","degree","hybrid"), 
                    absolute = TRUE, diagonal = 0, ...)
{
  ###########################
  #### MISSING ARGUMENTS ####
  ###########################
  
  if(missing(comm))
  {comm <- "walktrap"
  }else{comm <- comm}
  
  if(missing(diagonal))
  {diagonal <- 0
  }else{diagonal <- diagonal}
  
  if(missing(cent))
  {cent <- "strength"
  }else{cent <- match.arg(cent)}
  
  #######################
  #### MAIN FUNCTION ####
  #######################
  
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
#----

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
# Add signs----
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
#----

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
# Add signs----
# Updated 18.03.2020
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
    comc[which(is.na(comc))] <- stab
    
    # Round to 3 decimal places
    comm.str <- round(comc, 3)
    
  }
  
  return(comm.str)
}
#----