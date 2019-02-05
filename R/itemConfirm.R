#' Replicability of Items from EGA
#'
#' @description Based on the \code{\link{bootEGA}} results, this function 
#' computes and plots the number of times an item (variable) is estimated
#' in the same factor/dimension as originaly estimated by EGA. The output
#' also each item's likelihood (i.e., proprotion of bootstraps that it
#' appeared in) for every possible dimension (defined by the maximum number
#' of dimensions across the bootstrap samples).
#' 
#' @param bootega.obj A \code{\link[EGA]{bootEGA}} object
#' 
#' @param confirm A vector with community numbers or labels for each item
#' 
#' @param plot.ic Should the plot be produced?
#' Defaults to TRUE
#' 
#' @param item.rep A value for lowest likelihood allowed in \code{item.likelihood} output.
#' Removes noise from table to allow for easier interpretation.
#' Defaults to .10
#' 
#' @return Returns a list containing:
#' 
#' \item{item.confirm}{The proporton of times each item replicated
#' within the defined dimension}
#' 
#' \item{item.likelihood}{The proportion of times each item replicated
#' within each possible dimension. Dimensions greater than the maximum
#' number used in the \code{confirm} argument are labeled based on the
#' largest remaining components after the dimensions used to confirm.}
#' 
#' \item{wc}{A matrix containing the community membership values for
#' each bootstrapped sample. The values correspond to the values input
#' for the \code{confrim} argument.}
#' 
#' 
#' @examples
#' \dontrun{
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' 
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso",
#' type = "parametric", ncores = 4)
#' 
#' itemConfirm(boot.wmt, confirm = ega.wmt$wc)
#' }
#'
#' @references
#' Danon, L., Diaz-Guilera, A., Duch, J., & Arenas, A. (2005).
#' Comparing community structure identification.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{9}, P09008.
#' doi:\href{https://doi.org/10.1088/1742-5468/2005/09/P09008}{10.1088/1742-5468/2005/09/P09008}
#'
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' 
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Item Confirm function
itemConfirm <- function(bootega.obj, confirm, item.rep = .10, plot.ic = TRUE){
    
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
    
    #check if confirm is character
    if(is.character(confirm))
    {
        uni <- unique(confirm)
        num.comm <- confirm
        
        for(i in 1:length(uni))
        {num.comm[which(num.comm==uniq[i])] <- i}
    } else {num.comm <- confirm}
    
    #unique confirm cimensions
    uniq <- unique(num.comm)
    
    #number of confirm dimensions
    uniq.len <- length(num.comm)
    
    #confirm dimensions list
    confirm.list <- list()
    
    for(i in uniq)
    {confirm.list[[i]] <- which(confirm==uniq[i])}
    
    #letter membership vectors
    let.wc.mat <- matrix(NA, nrow = nrow(net), ncol = n)
    
    for(i in 1:n)
    {let.wc.mat[,i] <- letters[wc.mat[,i]]}
    
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
            #confirm list names
            conlist.names <- names(confirm.list[[j]])
            
            #match to iteration membership vector and construct comparison vector
            comp.vec <- num.vec[match(conlist.names,names(num.vec))]
            
            #confirm vector
            con.vec <- num.comm[match(conlist.names,names(num.comm))]
            
            #compare communities via the Rand method
            close.vec[uniq[j]] <- igraph::compare(con.vec,comp.vec,method="rand")
        }
        
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
            
            #insert values into letter membership matrix
            let.wc.mat[,i] <- let.vec
        }
    }
    
    
    #get letter dimension matrix into numbers
    num.wc.mat <- apply(let.wc.mat,1:2,as.numeric)
    
    #get frequency tables
    freq.list <- apply(num.wc.mat,1,table)
    
    #initialize item likelihood table
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
    {con.item[i] <- item.tab[i,confirm[i]]}
    
    #name item confirmation vector
    names(con.item) <- colnames(net)
    
    #Item Confirm
    itemCon <- con.item
    
    #Item Likelihood
    item.tab[which(item.tab<=item.rep)] <- 0
    
    if(any(apply(item.tab,2,sum)==0))
    {item.tab <- item.tab[,-which(apply(item.tab,2,sum)==0)]}
    
    item.tab[which(item.tab<=item.rep)] <- ""
    
    result <- list()
    
    #Plot
    comm <- confirm
    rain <- grDevices::rainbow(max(comm))
    
    item.repl <- data.frame(Item = names(itemCon),
                            Replication = itemCon,
                            Comm = factor(comm))
    
    
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
    
    if(plot.ic)
    {result$ic.plot <- ic.plot}
    
    #match row names to plot output
    itemCon <- itemCon[rev(match(ic.plot$data$Item,names(itemCon)))]
    
    #match row names to itemCon output
    itemLik <- as.data.frame(item.tab[match(names(itemCon),row.names(item.tab)),])
    
    #message for additional item likelihoods
    if(ncol(itemLik)<max(num.wc.mat))
    {message("Lower the item.rep argument to view item likelihoods for additional dimensions")}
    
    
    result$item.confirm <- itemCon
    result$item.likelihood <- itemLik
    result$wc <- num.wc.mat
    
    return(result)
}
#----
