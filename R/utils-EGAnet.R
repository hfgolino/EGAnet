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
# Unstandardized Network Loadings----
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

#' Mode
#' @description Computes the mode
#' 
#' @param v Vector
#'
#' @param fin.vec Vector
#' 
#' @return The mode of a vector
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @noRd
#' 
# Mode----
# Updated 18.03.2020
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
#----

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
#----
