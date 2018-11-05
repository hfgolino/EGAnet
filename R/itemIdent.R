#' Item Dimension Identification
#' @description Computes the within- and between-community strength of each item
#' for each community. Based on \code{\link{itemConfirm}}, researchers can flag
#' items that are not replicating well to identify which communities they are
#' falling between.
#' 
#' @param bootega.obj A \code{\link{bootEGA}} object
#' 
#' @param confirm A vector with community numbers or labels for each item
#' 
#' @param rep.val A replication value between 0 and 1. It's recommended
#' to first run \code{\link{itemConfirm}} to determine appropriate cut-off
#' value.
#' Defaults to .80
#' 
#' @return Returns a list containing:
#' 
#' \item{unstd}{The unstandardized within- and between-community
#' strength values for each node}
#' 
#' \item{std}{The standardized (proportion of strength)
#' within- and between-community strength values for each node}
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
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
itemIdent <- function (bootega.obj, confirm, rep.val = .80)
{
    n <- length(bootega.obj$bootGraphs)
    
    #Strength of within- and between-communities
    net <- bootega.obj$EGA$network
        
        ident.item <- function (net, comm)
        {
            comcat <- NetworkToolbox::comcat(net,comm,metric="each")
            stable <- NetworkToolbox::stable(net,comm)
            
            for(q in 1:nrow(comcat))
            {comcat[q,which(is.na(comcat[q,]))] <- stable[q]}
            
            comm.str <- comcat[,order(colnames(comcat))]
            comm.str <- round(comm.str,3)
            
            return(comm.str)
        }
        
        item.con <- itemConfirm(bootega.obj,confirm,plot.ic=FALSE)
        
        item.id.samps <- list()
        
        for(m in 1:n)
        {item.id.samps[[m]] <- ident.item(bootega.obj$bootGraphs[[m]], confirm)}
        
        #Unstandardized
        unstd.item.id <- round(apply(simplify2array(item.id.samps),1:2, mean),3)
        
        #Standardized (proportion)
        std.item.id <- round(unstd.item.id/rowSums(unstd.item.id),3)
        
        
        #less reliable items
        prob.item <- which(item.con$item.confirm<=rep.val)
        
        prob.vec <- item.con$item.confirm
        
        prob.vec[prob.item] <- "X"
        prob.vec[-prob.item] <- ""
        
        #Unstandardize
        unstd.item.ident <- as.data.frame(cbind(confirm,unstd.item.id,prob.vec))
        colnames(unstd.item.ident) <- c("Dimension",colnames(unstd.item.id),paste("Rep<=",rep.val,sep=""))
        unstd <- unstd.item.ident[order(unstd.item.ident$Dimension,decreasing=FALSE),]
        
        #Standaridize
        std.item.ident <- as.data.frame(cbind(confirm,std.item.id,prob.vec))
        colnames(std.item.ident) <- c("Dimension",colnames(std.item.id),paste("Rep<=",rep.val,sep=""))
        std <- std.item.ident[order(std.item.ident$Dimension,decreasing=FALSE),]
        
        result <- list()
        
        result$unstd <- unstd
        result$std <- std
        
        return(result)
}
