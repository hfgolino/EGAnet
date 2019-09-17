#' Network Scores
#'
#' @description This function computes network scores for
#' factor analysis models. Network scores are computed based on 
#' each node's \code{\link[NetworkToolbox]{strength}} within each
#' community (i.e., factor) in the network. These values are used
#' as network "factor loadings" for the weights of each item. Notably,
#' network analysis allows nodes to load onto more than one factor.
#' These loadings are considered in the factor scores. In addition,
#' if the construct is a hierarchy (e.g., personality questionnaire;
#' items in facet scales in a trait domain), then an overall
#' score can be computed (see argument \code{general}). These overall
#' scores are computed using \code{\link[NetworkToolbox]{comm.close}}
#' as weights, which are roughly similar to general factor loadings in a
#' CFA model (see Christensen, Golino, & Silvia, 2019). The score
#' estimates are roughly equivalent to the Maximum Likelihood method in
#' \code{lavaan}'s \code{\link[lavaan]{cfa}} function. An important difference
#' is that the network scores account for cross-loadings in their
#' estimation of scores.
#'
#' @param data Matrix or data frame.
#' Must be a dataset
#' 
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object. 
#' An adjacency matrix of network data
#' 
#' @param wc Numeric. 
#' A vector of community assignments.
#' Not necessary if an \code{\link[EGAnet]{EGA}} object
#' is input for argument \code{A}
#' 
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#' 
#' @return Returns a list containing:
#' 
#' \item{unstd.scores}{The unstandardized network scores for each participant
#' and community (including the overall score)}
#' 
#' \item{std.scores}{The standardized network scores for each participant
#' and community (including the overall score)}
#' 
#' \item{commCor}{Partial correlations between the specified or identified communities}
#' 
#' \item{loads}{Standardized network loadings for each item in each dimension
#' (computed using \code{\link[EGAnet]{net.loads}})}
#' 
#' @details For more details, type \code{vignette("Network_Scores")}
#' 
#' @examples 
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#'  # Estimate EGA
#'  ega.wmt <- EGA(wmt)
#' }
#'  
#'  # Network scores
#'  net.scores(data = wmt, A = ega.wmt)
#' 
#' @references 
#' Christensen, A. P. (2018).
#' NetworkToolbox: Methods and measures for brain, cognitive, and psychometric network analysis in R.
#' \emph{The R Journal}, \emph{10}, 422-439.
#' doi: \href{https://doi.org/10.32614/RJ-2018-065}{10.32614/RJ-2018-065}
#' 
#' Christensen, A. P., Golino, H. F., & Silvia, P. J. (2019).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/ktejp}{10.31234/osf.io/ktejp}
#' 
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#' 
#Network Scores
net.scores <- function (data, A, wc, ...)
{
    ####Missing arguments checks####
    if(missing(data))
    {stop("Argument 'data' is required for analysis")}
    
    # Detect if input is an 'EGA' object
    if(class(A) == "EGA")
    {
        #Compute network loadings
        P <- net.loads(A)$std
        
        # Grab communities
        wc <- A$wc
        
        # Replace 'A' with 'EGA' network
        A <- A$network
        
    }else if(missing(A))
    {stop("Adjacency matrix is required for analysis")
    }else if(missing(wc)) #Default to single  variable
    {wc <- rep(1,ncol(data))
    }else{
        #Compute network loadings
        {P <- net.loads(A=A,wc=wc)$std}
    }
    ####Missing arguments checks####
    
    #Number of factors
    nfacts <- length(unique(wc))
    
    #Initialize factor result matrix
    if(nfacts > 1)
    {fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = (nfacts + 1)))
    }else{fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = nfacts))}
    
    ####NETWORK SCORE FUNCTION####
    net.score.fxn <- function(loads, data)
    {
        #Initialize participant  scores
        net.sco <- matrix(0, nrow = nrow(data), ncol = ncol(loads))
        
        #Compute  factor scores (ML)
        for(i in 1:ncol(loads))
        {
            #Network loadings for each factor
            f.load <- loads[which(loads[,i]!=0),i]
            
            #Grab items associated with factor
            dat <- data[,names(f.load)]
            
            #Grab std dev of items associated with factor 
            f.sds <- apply(dat,2,sd,na.rm = TRUE)
            
            #Obtain relative weights
            rel <- f.load / f.sds
            rel.wei <- rel / sum(rel)
            
            #Compute scores
            net.sco[,i] <- as.vector(rowSums(t(t(dat) * rel.wei)))
        }
        
        colnames(net.sco) <- colnames(loads)
        
        return(net.sco)
    }
    ####NETWORK SCORE FUNCTION####
    
    #Populate factor result matrix
    net.sco <- net.score.fxn(P, data)
    fact.res[,1:nfacts] <- net.sco
    
    if(nfacts > 1)
    {colnames(fact.res)[1:nfacts] <- colnames(P)
    }else{colnames(fact.res) <- "1"}
    
    #Compute partial correlations between factors
    invS <- -cov2cor(solve(cov(net.sco, use = "pairwise.complete.obs")))
    diag(invS) <- 1
    C <- invS
    
    if(nfacts > 1)
    {
        #Compute general network loadings
        Pg <- NetworkToolbox::comm.close(A = A, comm = wc)
        
        #Overall score
        G <- rowSums(t(t(net.sco) * Pg))
        fact.res[,(nfacts + 1)] <- G
        colnames(fact.res)[nfacts + 1] <- "Overall"
        
        #Re-compute partial correlations between factors
        invS <- -cov2cor(solve(cov(net.sco, use = "pairwise.complete.obs")))
        diag(invS) <- 1
        C <- invS
    }
    
    #Results
    res <- list()
    res$unstd.scores <- round(fact.res,3)
    res$std.scores <- round(apply(fact.res,2,scale),3)
    res$commCor <- C
    res$loads <- P
    
    return(res)
}
