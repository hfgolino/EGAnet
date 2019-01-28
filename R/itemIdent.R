#' Item Dimension Identification
#' 
#' @description Computes the within- and between-community strength of each item
#' for each community. Based on \code{\link[EGA]{itemConfirm}}, researchers can flag
#' items that are not replicating well to identify which communities they are
#' falling between.
#' 
#' @param bootega.obj A \code{\link[EGA]{bootEGA}} object
#' 
#' @param confirm A vector with community numbers or labels for each item
#' 
#' @param rep.val A replication value between 0 and 1. It's recommended
#' to first run \code{\link[EGA]{itemConfirm}} to determine appropriate cut-off
#' value.
#' Defaults to .80
#' 
#' @return Returns a matrix of the unstandardized within- and between-community
#' strength values for each node
#' 
#' @examples
#' \dontrun{
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' 
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso",
#' type = "parametric", ncores = 4, confirm = ega.wmt$wc)
#' 
#' itemConfirm(boot.wmt, confirm = ega.wmt$wc, plot.ic = TRUE)
#' 
#' itemIdent(boot.wmt, confirm = ega.wmt$wc, rep.val = .80)
#' }
#' 
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Item Identification function
itemIdent <- function (bootega.obj, confirm, rep.val = .80)
{
    n <- length(bootega.obj$bootGraphs)
    
    #Strength of within- and between-communities
    net <- bootega.obj$EGA$network
        
        ident.item <- function (net, comm)
        {
            comc <- NetworkToolbox::comcat(net,comm,metric="each")
            stab <- NetworkToolbox::stable(net,comm)
            
            for(q in 1:nrow(comc))
            {comc[q,which(is.na(comc[q,]))] <- stab[q]}
            
            comm.str <- comc[,order(colnames(comc))]
            comm.str <- round(comm.str,3)
            
            return(comm.str)
        }
        
        item.con <- itemConfirm(bootega.obj,confirm,plot.ic=FALSE)
        
        col <- ncol(item.con$item.likelihood)
        
        item.id.samps <- list()
        
        for(m in 1:n)
        {item.id.samps[[m]] <- ident.item(bootega.obj$bootGraphs[[m]], confirm)}
        
        #Unstandardized
        unstd.item.id <- round(apply(simplify2array(item.id.samps),1:2, mean),3)
        
        #Standardized (proportion)
        #std.item.id <- round(unstd.item.id/rowSums(unstd.item.id),3)
        
        
        #less reliable items
        con.item <- item.con$item.confirm[row.names(unstd.item.id)]
        prob.item <- which(con.item<=rep.val)
        
        prob.vec <- con.item
        
        prob.vec[prob.item] <- "X"
        prob.vec[-prob.item] <- ""
        
        #Unstandardize
        unstd.item.ident <- as.data.frame(cbind(confirm,unstd.item.id,prob.vec))
        colnames(unstd.item.ident) <- c("Dimension",colnames(unstd.item.id),paste("Rep<=",rep.val,sep=""))
        unstd <- unstd.item.ident[match(names(item.con$item.confirm),row.names(unstd.item.ident)),]
        
        #Standaridize
        #std.item.ident <- as.data.frame(cbind(confirm,std.item.id,prob.vec))
        #colnames(std.item.ident) <- c("Dimension",colnames(std.item.id),paste("Rep<=",rep.val,sep=""))
        #std <- std.item.ident[order(std.item.ident$Dimension,decreasing=FALSE),]
        
        #result <- list()
        
        #result$unstd <- unstd
        #result$std <- std
        
        #return(result)
        
        return(unstd)
}
