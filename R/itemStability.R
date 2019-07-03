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
#' boot.wmt <- bootEGA(data = wmt2, n = 100, typicalStructure = TRUE,
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
#Item Confirm function
itemStability <- function(bootega.obj, orig.wc, item.freq = .10, plot.item.rep = TRUE){
    
    if(class(bootega.obj) != "bootEGA")
    {stop("Input for 'bootega.obj' is not a 'bootEGA' object")}
    
    #mode function
    mode <- function(v, numeric = FALSE)
    {
        #unique values
        uniqv <- unique(v)
        
        if(!numeric)
        {
            #identify letters in values
            uniqv <- uniqv[which(is.na(suppressWarnings(as.numeric(uniqv))))]
        }
        
        #find mode of letters
        uniqv <- uniqv[which.max(tabulate(match(v, uniqv)))]
        
        return(uniqv)
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
    
    #number of original dimensions
    uniq.len <- length(num.comm)
    
    #original dimensions list
    orig.wc.list <- list()
    
    for(i in uniq)
    {orig.wc.list[[i]] <- which(orig.wc==uniq[i])}
    
    names(orig.wc.list) <- uniq
    
    #letter membership vectors
    let.wc.mat <- matrix(NA, nrow = nrow(net), ncol = n)
    
    for(i in 1:n)
    {let.wc.mat[,i] <- letters[wc.mat[,i]]}
    
    ##########################################################
    #### ITEM FREQUENCY AND DIMENSION REPLICATION RESULTS ####
    ##########################################################
    
    #identify confirm membership within bootstrapped memberships
    for(i in 1:n)
    {
        #iteration letter membership vector
        let.vec <- let.wc.mat[,i]
        #number membership vector
        num.vec <- wc.mat[,i]
        #unique int.let membership
        num.uniq <- unique(num.vec)
        #number of unique int.let membership
        num.len <- length(num.uniq)
        
        #initialize vector for closest dimension closest to original
        close.vec <- vector("numeric",length=max(uniq))
        names(close.vec) <- uniq
        
        for(j in uniq)
        {
            #original dimension list names
            orig.wc.names <- names(orig.wc.list[[paste(j)]])
            
            #match to iteration membership vector and construct comparison vector
            comp.vec <- num.vec[match(orig.wc.names,names(num.vec))]
            
            #original dimension vector
            con.vec <- num.comm[match(orig.wc.names,names(num.comm))]
            
            #compare communities via the Rand method
            close.vec[uniq[j]] <- igraph::compare(con.vec,comp.vec,method="rand")
        }
        
        #order close.vec by most items
        ##identify rand with largest number of items
        close.ord <- close.vec
        close.clone <- close.vec
        
        for(y in 1:length(close.vec))
        {
            target.close <- mode(num.comm[which(!is.na(match(num.comm,names(close.clone)[which(close.clone==max(close.clone))])))],TRUE)
            
            target.dim <- which(names(close.clone)==target.close)
            
            dim <- names(close.clone)[target.dim]
            
            close.ord[y] <- dim
            
            close.clone <- close.clone[-target.dim]
        }
        
        close.vec <- close.vec[order(close.ord)]
        
        if(all(close.vec==1)) #perfect match
        {
            #if number of dimensions is less than confirmatory
            if(max(num.comm)!=max(num.vec))
            {
                #create new vector of iteration vector
                new.vec <- num.vec
                
                #identify overlapping dimensions
                target <- intersect(num.comm,new.vec)
                #number of overlapping dimensions
                target.len <- length(target)
                
                #initialize Rand index vector
                rand.vec <- vector("numeric",length=target.len)
                
                #initilize count
                count <- 0
                
                #compute Rand index for each iteration dimension
                for(o in 1:target.len)
                {
                    #increase count
                    count <- count + 1
                    
                    #find iteration dimensions matching overlapping dimensions
                    target.val <- which(num.vec==target[o])
                    
                    #compute Rand index for overlapping dimensions
                    rand.vec[count] <- igraph::compare(new.vec[target.val],num.comm[target.val],method="rand")
                }
                
                #name Rand index vector with iteration dimensions
                names(rand.vec) <- paste(target)
                
                #identify Rand index less than 1
                target.rand <- as.numeric(names(rand.vec)[which(rand.vec!=1)])
                
                #number of Rand index less than 1
                rand.len <- length(target.rand)
                
                #identify confirmatory dimension that most represents the iteration dimension
                for(p in 1:rand.len)
                {
                    #compute mode for confirmatory dimension that matches iteration dimension's elements
                    mode.val <- mode(num.comm[which(new.vec==target.rand[p])], numeric = TRUE)
                    
                    #apply confirmatory dimension to iteration dimension's elements
                    num.vec[which(new.vec == target.rand[p])] <- mode.val
                }
                
                #if the above fails, then different strategy
                if(all(num.vec==new.vec))
                {
                    #create new vector of iteration vector
                    new.vec <- num.vec
                    
                    #identify dimensions that are different between iteration and confirmatory dimensions
                    target <- setdiff(new.vec,num.comm)
                    
                    #number that are different
                    target.len <- length(target)
                    
                    #for the different dimensions
                    #identify confirmatory dimension that most represents the iteration dimension
                    for(o in 1:target.len)
                    {
                        #find iteration dimensions matching different dimensions
                        target.val <- which(new.vec==target[o])
                        
                        #compute mode for confirmatory dimension that matches iteration dimension's elements
                        mode.val <- mode(num.comm[target.val], numeric = TRUE)
                        
                        #apply confirmatory dimension to iteration dimension's elements
                        num.vec[which(new.vec == target[o])] <- mode.val
                    }
                }
                
                let.wc.mat[,i] <- num.vec
                
            }else{let.wc.mat[,i] <- num.comm}
        }else{
            
            while(length(close.vec)!=0)
            {
                #determine dimension closest to original
                max.pos <- as.numeric(names(close.vec)[which.max(close.vec)])
                
                #positons of closest
                close.pos <- which(num.comm==max.pos)
                
                #target letter
                target.let <- mode(let.vec[close.pos])
                #target vector
                target.vec <- which(let.vec==target.let)
                
                #insert appropriate original EGA community number
                let.vec[target.vec] <- rep(max.pos,length(target.vec))
                
                #remove value from closest dimension vector
                close.vec <- close.vec[-which(names(close.vec)==max.pos)]
            }
            
            #check for remaining dimensions
            rem.dim <- intersect(let.vec,let.wc.mat[,i])
            
            if(max(num.uniq)>max(uniq))
            {
                for(u in (max(uniq)+1):max(num.uniq))
                {
                    #identify largest remaining dimension
                    target.rem <- mode(let.vec)
                    #target remaming values
                    target.rem.vec <- which(let.vec==target.rem)
                    
                    #insert next dimension number
                    let.vec[target.rem.vec] <- rep(u,length(target.rem.vec))
                }
            }else if(any(is.na(suppressWarnings(as.numeric(let.vec)))))
            {
                #identify missing numeric value
                uniq.vals <- as.numeric(unique(let.vec[which(!is.na(suppressWarnings(as.numeric(let.vec))))]))
                
                #target value
                target.val <- setdiff(uniq,uniq.vals)
                
                #target remaming values
                target.rem.vec <- which(let.vec==rem.dim)
                
                #insert missing dimension number
                let.vec[target.rem.vec] <- rep(target.val,length(target.rem.vec))
            }
            
            #check for remaining dimensions
            rem.dim <- intersect(let.vec,let.wc.mat[,i])
            
            if(length(rem.dim)!=0)
            {
                target.max <- max(as.numeric(let.vec[suppressWarnings(!is.na(as.numeric(let.vec)))]))
                
                for(r in 1:length(rem.dim))
                {
                    target.rem <- which(let.vec==rem.dim[r])
                    
                    let.vec[target.rem] <- rep(target.max + r, length(target.rem))
                }
            }
            
            #insert values into letter membership matrix
            let.wc.mat[,i] <- let.vec
        }
    }
    
    
    #get letter dimension matrix into numbers
    num.wc.mat <- apply(let.wc.mat,1:2,as.numeric)
    
    #get frequency tables
    freq.list <- apply(num.wc.mat,1,table)
    
    #initialize item frequency table
    item.tab <- matrix(0,nrow=nrow(net),ncol=max(num.wc.mat))
    
    #name columns and rows
    colnames(item.tab) <- paste(seq(1,max(num.wc.mat,1)))
    row.names(item.tab) <- colnames(net)
    
    #insert proportion values into item likelihod table
    for(i in 1:length(freq.list))
    {
        prop.int <- freq.list[[i]]/n
        
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
    
    if(plot.item.rep)
    {result$plot.itemStability <- ic.plot}
    
    #match row names to plot output
    itemCon <- itemCon[rev(match(ic.plot$data$Item,names(itemCon)))]
    
    #match row names to itemCon output
    itemLik <- as.data.frame(item.tab[match(names(itemCon),row.names(item.tab)),])
    
    #message for additional item likelihoods
    if(ncol(itemLik)<max(num.wc.mat))
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
    
    max.wc <- max(num.wc.mat)
    
    for(m in 1:n)
    {
        item.id.samps[[m]] <- net.loads(bootega.obj$bootGraphs[[m]], num.wc.mat[,m])$std
        
        dims <- length(unique(num.wc.mat[,m]))
        
        if(dims!=max.wc)
        {
            diff <- max.wc - ncol(item.id.samps[[m]])
            
            diff.wc <- setdiff(seq(1,max.wc,1),unique(num.wc.mat[,m]))
            
            col.names <- c(colnames(item.id.samps[[m]]),paste(diff.wc))
            
            for(i in 1:diff)
            {item.id.samps[[m]] <- cbind(item.id.samps[[m]],rep(NA,nrow(item.id.samps[[m]])))}
            
            colnames(item.id.samps[[m]]) <- col.names
            
            item.id.samps[[m]] <- item.id.samps[[m]][,order(as.numeric(colnames(item.id.samps[[m]])))]
        }
    }
    
    #Unstandardized
    unstd.item.id <- round(apply(simplify2array(item.id.samps),1:2, mean, na.rm=TRUE),3)
    colnames(unstd.item.id) <- paste(seq(1,max(num.wc.mat),1))
    
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
    result$wc <- num.wc.mat
    
    return(result)
}
#----
