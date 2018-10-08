#' Replicability of items from EGA
#' @description Based on the \code{\link{bootEGA}} results, this function 
#' computes and plots the number of times an item (variable) is estimated
#' in the same factor/dimension as originaly estimated by EGA. The output
#' will also contain how often a dimension is replicated.
#' 
#' 
#' @param bootega.obj A \code{\link{bootEGA}} object
#' 
#' @param confirm A vector with community numbers or labels for each item
#' 
#' @param plot.ic Should the plot be produced?
#' Defaults to TRUE
#' 
#' @return Returns a list containing:
#' 
#' \item{dim.confirm}{The proportion of times each dimension replicated
#' (i.e., all items were in the exact same dimension)}
#' 
#' \item{item.confirm}{The proporton of times each item replicated
#' within the defined dimension}
#' 
#' \item{dim.nmi.table}{Computes the Normalized Mutual Informaiton
#' for community replication. Larger numbers mean dimensions and
#' items are replicating closer to the input (\code{confirm})
#' dimensions (see Danon et al. 2005 for more details)}
#' 
#' @examples
#' \dontrun{
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' 
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso",
#' type = "parametric", ncores = 4)
#' 
#' item.confirm(boot.wmt, confirm = ega.wmt$wc)
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
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
item.confirm <- function(bootega.obj, confirm, plot.ic = TRUE){
    require(ggpubr)
    
    #mode function for item confirm
    mode <- function(v)
    {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    n <- length(bootega.obj$bootGraphs) 
    
    #Initiate confirm matrix
    uniq <- unique(confirm)
    confirm.dim <- matrix(NA, nrow = n, ncol = length(uniq))
    dim.nmi <- vector("numeric", length = n)
    item.confirm <- matrix(NA, nrow = n, ncol = length(confirm))
    
    #check if confirm is character
    if(is.character(confirm))
    {
        uni <- unique(confirm)
        num.comm <- confirm
        
        for(i in 1:length(uni))
        {num.comm[which(num.comm==uniq[i])] <- i}
    } else {num.comm <- confirm}
    
    for (m in 1:n)
    {
        #normalized mutual information of community comparisons
        dim.nmi[m] <- igraph::compare(bootega.obj$boot.wc[[m]]$membership, num.comm, method="nmi")
        
        
        for(i in 1:length(uniq))
        {
            dim.items <- which(confirm==uniq[i])
            target.dim <- bootega.obj$boot.wc[[m]]$membership[dim.items]
            
            #count code
            uniq.dim <- unique(target.dim)
            if(length(uniq.dim)==1){confirm.dim[m,i] <- 1}else{confirm.dim[m,i] <- 0}
            
            #Check if item is confirmed within dimension
            if(length(uniq.dim)>1)
            {
                target.mode <- mode(target.dim)
                non.con <- dim.items[which(target.dim!=target.mode)]
                con <- setdiff(dim.items,non.con)
                item.confirm[m,non.con] <- 0
                item.confirm[m,con] <- 1
            } else {item.confirm[m,dim.items] <- 1}
        }
    }
    
    #Proportion of times dimension is confirmed
    con.dim <- (colSums(confirm.dim))/n
    names(con.dim) <- uniq
    
    #Proportion of times item is confirmed
    con.item <- (colSums(item.confirm)/n)
    names(con.item) <- as.character(bootega.obj$EGA$dim.variables$items)
    
    #Tables for nmi
    dim.nmi.table <- vector("numeric", length = 3)
    dim.nmi.table[1] <- mean(dim.nmi)
    dim.nmi.table[2] <- sd(dim.nmi)
    dim.nmi.se <- (1.253 * sd(dim.nmi))/sqrt(n)
    dim.nmi.table[3] <- mean(dim.nmi) - dim.nmi.se
    dim.nmi.table[4] <- mean(dim.nmi) + dim.nmi.se
    names(dim.nmi.table) <- c("mean","sd","lower.ci","upper.ci")
    
    result <- list()
    
    #Plot
    if(plot.ic == TRUE)
    {
        comm <- bootega.obj$EGA$wc
        rain <- grDevices::rainbow(max(comm))
        
        item.rep <- data.frame(Item = bootega.obj$EGA$dim.variables$items,
                               Rep = con.item,
                               Comm = factor(comm))
        
        
        ic.plot <- ggdotchart(item.rep, x = "Item", y = "Rep",
                              group = "Comm", color = "Comm",
                              palette = rain,
                              legend.title = "EGA Communities",
                              sorting = NULL,
                              add = "segments",
                              rotate = TRUE,
                              dot.size = 6,
                              label = round(item.rep$Rep, 2),
                              font.label = list(color = "black", size = 8,
                                                vjust = 0.5),
                              ggtheme = theme_pubr()
        )
        
        result$ic.plot <- ic.plot
    }
    
    result$dim.confirm <- con.dim
    result$dim.nmi.table <- dim.nmi.table
    result$item.confirm <- con.item
    
    return(result)
}
