#' Item Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, this function
#' computes and plots the number of times an item (variable) is estimated
#' in the same factor/dimension as originaly estimated by \code{\link[EGAnet]{EGA}} (\code{item.replication}).
#' The output also contains each item's replication frequency (i.e., proprotion of
#' bootstraps that an item appeared in each dimension; \code{item.dim.rep}) as well as the average
#' network loading for each item in each dimension (\code{item.loadings}).
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param orig.wc Numeric or character.
#' A vector with community numbers or labels for each item.
#' Typically uses community results (\code{wc}) from \code{\link[EGAnet]{EGA}}
#'
#' @param item.freq A value for lowest frequency allowed in \code{item.dim.rep} output.
#' Removes noise from table to allow for easier interpretation.
#' Defaults to \code{.10}
#'
#' @param plot.item.rep Should the plot be produced for \code{item.replication}?
#' If \code{TRUE}, then a plot for the \code{item.replication} output will be produced.#'
#' Defaults to \code{TRUE}
#'
#' @return Returns a list containing:
#'
#' \item{item.replication}{The proporton of times each item replicated
#' within the defined dimension}
#'
#' \item{item.dim.rep}{The proportion of times each item replicated
#' within each possible dimension. Dimensions greater than the maximum
#' number used in the \code{orig.wc} argument are labeled based on the
#' largest remaining components after the dimensions used to \code{orig.wc}}
#'
#' \item{item.loadings}{Matrix of the average standardied network loading
#' (computed using \code{\link[EGAnet]{net.loads}}) for each item in each dimension}
#'
#' \item{wc}{A matrix containing the community membership values for
#' each bootstrapped sample. The values correspond to the values input
#' for the \code{orig.wc} argument}
#' 
#' \item{plot.itemStability}{A plot of the number of times each item
#' replicates in its original community membership (\code{orig.wc})}
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA network
#' ega.wmt <- EGA(data = wmt, model = "glasso")
#'
#' # Estimate dimension stability
#' boot.wmt <- bootEGA(data = wmt, n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso",
#' type = "parametric", ncores = 4)
#'
#' }
#'
#' # Estimate item stability statistics
#' itemStability(boot.wmt, orig.wc = ega.wmt$wc)
#'
#' @references
#' Danon, L., Diaz-Guilera, A., Duch, J., & Arenas, A. (2005).
#' Comparing community structure identification.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{9}, P09008.
#' <doi:10.1088/1742-5468/2005/09/P09008>
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Item Stability function
itemStability <- function(bootega.obj, orig.wc, item.freq = .10, plot.item.rep = TRUE){
    
    if(class(bootega.obj) != "bootEGA")
    {stop("Input for 'bootega.obj' is not a 'bootEGA' object")}
    
    #mode function
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
    
    #number of bootstrapped networks
    n <- length(bootega.obj$bootGraphs)
    
    #original EGA network
    net <- bootega.obj$EGA$network
    
    #grab membership vectors
    wc.mat <- matrix(NA, nrow = nrow(net), ncol = n)
    
    for(i in 1:n)
    {wc.mat[,i] <- bootega.obj$boot.wc[[i]]$membership}
    
    #grab item names
    row.names(wc.mat) <- row.names(net)
    
    #check if orig.wc is character
    if(is.character(orig.wc))
    {
        uni <- unique(orig.wc)
        num.comm <- orig.wc
        
        for(i in 1:length(uni))
        {num.comm[which(num.comm==uniq[i])] <- i}
    } else {num.comm <- orig.wc}
    
    #unique original cimensions
    uniq <- unique(num.comm)
    uniq <- uniq[order(uniq)]
    
    #initialize final matrix
    final.mat <- matrix(0, nrow = length(num.comm), ncol = n)
    
    ##########################################################
    #### ITEM FREQUENCY AND DIMENSION REPLICATION RESULTS ####
    ##########################################################
    
    #identify confirm membership within bootstrapped memberships
    for(i in 1:n)
    {
        #new membership vector
        new.vec <- wc.mat[,i]
        
        #unique new membership
        new.uniq <- unique(new.vec)
        
        #converge based on maximum number of dimensions
        if(max(num.comm) > max(new.vec))
        {
            #initialize rand and length vector
            rand <- vector("numeric", length = max(new.vec))
            names(rand) <- new.uniq
            len <- rand
            
            for(j in new.uniq)
            {
                #target nodes
                target <- which(new.vec==j)
                
                #lengths of target
                len[paste(j)] <- length(target)
                
                #compute rand index
                rand[paste(j)] <- igraph::compare(new.vec[target],num.comm[target],method="rand")
            }
            
            #order rand by highest rand index and then number of items
            rand.ord <- rand[order(rand, len, decreasing = TRUE)]
            
            #initialize final vector
            final.vec <- vector("numeric", length = length(num.comm))
            names(final.vec) <- names(num.comm)
            
            #insert new values into final vector
            for(j in as.numeric(names(rand.ord)))
            {
                #identify target
                new.target <- which(new.vec==j)
                
                #identify mode
                target.mode <- mode(num.comm[new.target], final.vec)
                
                #insert into final vector
                final.vec[new.target] <- rep(target.mode)
            }
            
        }else if(max(num.comm) < max(new.vec))
        {
            #initialize rand and length vector
            rand <- vector("numeric", length = max(new.vec))
            names(rand) <- new.uniq
            len <- rand
            
            for(j in new.uniq)
            {
                #target nodes
                target <- which(new.vec==j)
                
                #lengths of target
                len[paste(j)] <- length(target)
                
                #compute rand index
                rand[paste(j)] <- igraph::compare(new.vec[target],num.comm[target],method="rand")
            }
            
            #order rand by highest rand index and then number of items
            rand.ord <- rand[order(rand, len, decreasing = TRUE)]
            
            #initialize final vector
            final.vec <- vector("numeric", length = length(num.comm))
            names(final.vec) <- names(num.comm)
            
            #insert new values into final vector
            for(j in as.numeric(names(rand.ord)))
            {
                #identify target
                new.target <- which(new.vec==j)
                
                #identify mode
                target.mode <- mode(num.comm[new.target], final.vec)
                
                #insert into final vector
                final.vec[new.target] <- rep(target.mode)
            }
            
            #identify number of extra dimensions
            extra.dim <- unique(new.vec[which(is.na(final.vec))])
            
            #initialize extra dimension length vector
            extra.len <- vector("numeric", length = length(extra.dim))
            names(extra.len) <- extra.dim
            
            #initialize count
            count <- 0
            
            #order length of extra dimensions
            for(j in extra.dim)
            {
                #increase count
                count <- count + 1
                
                #length of extra dimensions
                extra.len[count] <- length(which(new.vec==j))
            }
            
            el.ord <- extra.len[order(extra.len, decreasing = TRUE)]
            
            #reset count
            count <- 0
            
            #insert extra dimensions into final vector
            for(j in 1:length(el.ord))
            {
                #increase count
                count <- count + 1
                
                #target extra dimension
                target.ed <- as.numeric(names(el.ord)[j])
                
                #insert dimensions into final vector
                final.vec[which(new.vec==target.ed)] <- (max(num.comm) + count)
            }
            
        }else{
            
            #initialize rand and length vector
            rand <- vector("numeric", length = max(new.vec))
            names(rand) <- new.uniq
            len <- rand
            
            for(j in new.uniq)
            {
                #target nodes
                target <- which(new.vec==j)
                
                #lengths of target
                len[paste(j)] <- length(target)
                
                #compute rand index
                rand[paste(j)] <- igraph::compare(new.vec[target],num.comm[target],method="rand")
            }
            
            #order rand by highest rand index and then number of items
            rand.ord <- rand[order(rand, len, decreasing = TRUE)]
            
            #initialize final vector
            final.vec <- vector("numeric", length = length(num.comm))
            names(final.vec) <- names(num.comm)
            
            #insert new values into final vector
            for(j in as.numeric(names(rand.ord)))
            {
                #identify target
                new.target <- which(new.vec==j)
                
                #identify mode
                target.mode <- mode(num.comm[new.target], final.vec)
                
                #insert into final vector
                final.vec[new.target] <- rep(target.mode)
            }
        }
        
        #insert final vector into final matrix
        final.mat[,i] <- final.vec
    }
    
    #get frequency tables
    freq.list <- apply(final.mat,1,table)
    
    #change to matrix
    if(is.list(freq.list))
    {
        #initialize new matrix
        new.mat <- matrix(0, nrow = max(final.mat,na.rm=TRUE), ncol = length(freq.list))
        
        #name rows and columns
        row.names(new.mat) <- paste(1:max(final.mat,na.rm=TRUE))
        colnames(new.mat) <- colnames(net)
        
        #insert values
        for(i in 1:ncol(new.mat))
        {new.mat[names(freq.list[[i]]),i] <- freq.list[[i]]}
        
        freq.list <- new.mat
    }
    
    #initialize item frequency table
    item.tab <- matrix(0,nrow=nrow(net),ncol=max(final.mat,na.rm=TRUE))
    
    #name columns and rows
    colnames(item.tab) <- paste(seq(1,max(final.mat,na.rm=TRUE),1))
    row.names(item.tab) <- colnames(net)
    
    #insert proportion values into item likelihod table
    for(i in 1:ncol(freq.list))
    {
        prop.int <- freq.list[,i]/n
        
        item.tab[i,match(names(prop.int),colnames(item.tab))] <- prop.int
    }
    
    #initialize item confirmation vector
    con.item <- vector("numeric",length=nrow(item.tab))
    
    #grab confirmation value from proportion table
    for(i in 1:nrow(item.tab))
    {con.item[i] <- item.tab[i,orig.wc[i]]}
    
    #name item confirmation vector
    names(con.item) <- colnames(net)
    
    #Item Confirm
    itemCon <- con.item
    
    #Item Frequency
    item.tab[which(item.tab<=item.freq)] <- 0
    
    if(any(apply(item.tab,2,sum)==0))
    {item.tab <- item.tab[,-which(apply(item.tab,2,sum)==0)]}
    
    item.tab[which(item.tab<=item.freq)] <- ""
    
    result <- list()
    
    #Plot
    comm <- orig.wc
    rain <- rev(grDevices::rainbow(max(comm)))
    
    item.repl <- data.frame(Item = names(itemCon),
                            Replication = itemCon,
                            Comm = factor(comm,7:1))
    
    
    ic.plot <- ggpubr::ggdotchart(item.repl, x = "Item", y = "Replication",
                                  group = "Comm", color = "Comm",
                                  palette = rain,
                                  legend.title = "EGA Communities",
                                  add = "segments",
                                  rotate = TRUE,
                                  dot.size = 6,
                                  label = round(item.repl$Replication, 2),
                                  font.label = list(color = "black", size = 8,
                                                    vjust = 0.5),
                                  ggtheme = ggpubr::theme_pubr()
    )
    
    ic.plot <- ic.plot + ggplot2::ylim(c(0,1))
    
    if(plot.item.rep)
    {result$plot.itemStability <- ic.plot}
    
    #match row names to plot output
    itemCon <- itemCon[rev(match(ic.plot$data$Item,names(itemCon)))]
    
    #match row names to itemCon output
    itemLik <- as.data.frame(item.tab[match(names(itemCon),row.names(item.tab)),])
    
    #message for additional item likelihoods
    if(ncol(itemLik)<max(final.mat,na.rm=TRUE))
    {message("Lower the item.freq argument to view item frequencies for additional dimensions")}
    
    ##########################################################
    #### ITEM FREQUENCY AND DIMENSION REPLICATION RESULTS ####
    ##########################################################
    
    ##################################
    #### NETWORK LOADINGS RESULTS ####
    ##################################
    
    item.lik <- itemLik
    
    col <- ncol(item.lik)
    
    item.id.samps <- list()
    
    max.wc <- max(final.mat, na.rm = TRUE)
    
    for(m in 1:n)
    {
        item.id.samps[[m]] <- net.loads(bootega.obj$bootGraphs[[m]], final.mat[,m])$std
        
        dims <- length(unique(final.mat[,m]))
        
        if(dims!=max.wc)
        {
            if(any(colnames(item.id.samps[[m]])=="NA"))
            {item.id.samps[[m]] <- item.id.samps[[m]][,-which(colnames(item.id.samps[[m]])=="NA")]}
            
            diff <- max.wc - ncol(item.id.samps[[m]])
            
            diff.wc <- setdiff(seq(1,max.wc,1),unique(final.mat[,m]))
            
            col.names <- c(colnames(item.id.samps[[m]]),paste(diff.wc))
            
            for(i in 1:diff)
            {item.id.samps[[m]] <- cbind(item.id.samps[[m]],rep(NA,nrow(item.id.samps[[m]])))}
            
            colnames(item.id.samps[[m]]) <- na.omit(col.names)
            
            item.id.samps[[m]] <- item.id.samps[[m]][,order(as.numeric(colnames(item.id.samps[[m]])))]
        }
    }
    
    #Unstandardized
    unstd.item.id <- round(apply(simplify2array(item.id.samps),1:2, mean, na.rm=TRUE),3)
    colnames(unstd.item.id) <- paste(seq(1,max(final.mat, na.rm = TRUE),1))
    
    if(ncol(unstd.item.id)!=col)
    {
        rm.col <- setdiff(colnames(unstd.item.id),names(item.lik))
        
        target.col <- match(rm.col,colnames(unstd.item.id))
        
        unstd.item.id <- unstd.item.id[,-target.col]
    }
    
    item.lik.ord <- item.lik[match(row.names(unstd.item.id),row.names(item.lik)),]
    
    unstd.item.id[which(item.lik.ord=="")] <- ""
    
    unstd.item.id <- as.data.frame(unstd.item.id)
    
    #Unstandardize
    unstd.item.ident <- as.data.frame(cbind(orig.wc,unstd.item.id))
    colnames(unstd.item.ident) <- c("Dimension",colnames(unstd.item.id))
    itemLoads <- unstd.item.ident[match(names(itemCon),row.names(unstd.item.ident)),]
    
    ##################################
    #### NETWORK LOADINGS RESULTS ####
    ##################################
    
    
    result$item.replication <- itemCon
    result$item.dim.rep <- itemLik
    result$item.loadings <- itemLoads
    result$wc <- final.mat
    
    return(result)
}
#----
