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
    res$scores <- round(apply(fact.res,2,scale),3)
    res$commCor <- C
    res$loads <- P
    
    return(res)
}
