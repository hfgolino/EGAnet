#' Item Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, this function
#' computes and plots the number of times an item (variable) is estimated
#' in the same factor/dimension as originally estimated by \code{\link[EGAnet]{EGA}} (\code{item.replication}).
#' The output also contains each item's replication frequency (i.e., proportion of
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
#' If \code{TRUE}, then a plot for the \code{item.replication} output will be produced.
#' Defaults to \code{TRUE}
#'
#' @return Returns a list containing:
#'
#' \item{item.replication}{The proportion of times each item replicated
#' within the defined dimension}
#' 
#' \item{mean.dim.rep}{The average replication proportion of items replicating
#' in each dimension. More simply, the average of the \code{item.replication}
#' output for each dimension}
#'
#' \item{item.dim.rep}{The proportion of times each item replicated
#' within each possible dimension. Dimensions greater than the maximum
#' number used in the \code{orig.wc} argument are labeled based on the
#' largest remaining components after the dimensions used to \code{orig.wc}}
#'
#' \item{item.loadings}{Matrix of the average standardized network loading
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
#' \donttest{# Estimate EGA network
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' ega.wmt <- EGA(data = wmt, model = "glasso", plot.type = "qgraph")
#'
#' # Estimate dimension stability
#' boot.wmt <- bootEGA(data = wmt, iter = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso", plot.type = "qgraph",
#' type = "parametric", ncores = 2)
#' }
#' 
#' # Estimate item stability statistics
#' res <- itemStability(boot.wmt, orig.wc = ega.wmt$wc)
#' 
#' # Changing plot features (ggplot2)
#' ## Changing colors (ignore warnings)
#' ### qgraph Defaults
#' res$plot.itemStability + 
#'     ggplot2::scale_color_manual(values = rainbow(max(res$uniq.num)))
#' 
#' ### Pastel
#' res$plot.itemStability + 
#'     ggplot2::scale_color_brewer(palette = "Pastel1")
#'     
#' ## Changing Legend (ignore warnings)
#' res$plot.itemStability + 
#'     ggplot2::scale_color_discrete(labels = "Intelligence")
#'
#' @references
#' Christensen, A. P., & Golino, H. (2019).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/9deay}
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}.
#' \doi{10.1002/per.2265}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Item Stability function
#Updated 15.02.2021
itemStability <- function(bootega.obj, orig.wc, item.freq = .10, plot.item.rep = TRUE){
  
  # Check for 'bootEGA' object
  if(class(bootega.obj) != "bootEGA")
  {stop("Input for 'bootega.obj' is not a 'bootEGA' object")}
  
  # Number of bootstrapped networks
  n <- length(bootega.obj$bootGraphs)
  
  # Original EGA network
  net <- bootega.obj$EGA$network
  
  # Grab membership vectors
  wc.mat <- matrix(NA, nrow = length(orig.wc), ncol = n)
  
  for(i in 1:length(bootega.obj$bootGraphs))
  {wc.mat[,i] <- bootega.obj$boot.wc[[i]]}
  
  # Grab item names
  row.names(wc.mat) <- row.names(net)
  
  # Check if 'orig.wc' is character
  uni <- unique(orig.wc)
  
  if(is.character(orig.wc))
  {
    num.comm <- orig.wc
    
    for(i in 1:length(uni))
    {num.comm[which(num.comm==uni[i])] <- i}
  } else {num.comm <- orig.wc}
  
  # Convert to numeric vector
  num.comm <- as.numeric(num.comm)
  
  # Unique original dimensions
  uniq <- unique(num.comm)
  
  # Convert membership to target membership
  final.mat <- homogenize.membership(num.comm, wc.mat)
  
  ##########################################################
  #### ITEM FREQUENCY AND DIMENSION REPLICATION RESULTS ####
  ##########################################################
  
  #let user know results are being computed has started
  message("Computing results...", appendLF = FALSE)
  
  # Get proportion table
  item.tab <- proportion.table(final.mat)
  row.names(item.tab) <- colnames(net)
  
  if(is.character(uni))
  {colnames(item.tab)[1:length(uni)] <- uni}
  
  #initialize item confirmation vector
  con.item <- vector("numeric",length=nrow(item.tab))
  
  #grab confirmation value from proportion table
  for(i in 1:nrow(item.tab))
  {con.item[i] <- item.tab[i,paste(orig.wc[i])]}
  
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
  #rain <- rev(RColorBrewer::brewer.pal(max(num.comm), "Set1"))
  
  item.repl <- data.frame(Item = names(itemCon),
                          Replication = itemCon,
                          Comm = factor(comm, uni[order(uni)]))
  
  
  ic.plot <- ggpubr::ggdotchart(item.repl, x = "Item", y = "Replication",
                                group = "Comm", color = "Comm",
                                legend.title = "EGA Communities",
                                add = "segments",
                                rotate = TRUE,
                                dot.size = 6,
                                label = round(item.repl$Replication, 2),
                                font.label = list(color = "black", size = 8,
                                                  vjust = 0.5),
                                ggtheme = ggpubr::theme_pubr()
  )
  
  # Adjust y-axis
  ic.plot <- ic.plot + ggplot2::ylim(c(0,1))
  
  # Manually change alpha
  ic.plot$layers[[2]]$aes_params$alpha <- 0.7
  
  # Bold legend title
  ic.plot <- ic.plot + ggplot2::theme(
    legend.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold")
  )
  
  # Adjust item label sizes based on
  sizes <- seq(6,12,.25)
  ## Number of nodes
  nodes <- rev(seq(0, 200, length.out = length(sizes)))
  n.size <- min(which(length(orig.wc) > nodes))
  ## Number of characters in item name
  chars <- rev(seq(0,100, length.out = length(sizes)))
  ### Maximum characters in item name
  max.chars <- max(unlist(lapply(row.names(item.repl),nchar)))
  c.size <- min(which(max.chars > chars))
  # Text size
  text.size <- sizes[min(c(n.size,c.size))]
  
  ic.plot <- ic.plot + ggplot2::theme(axis.text.y = ggplot2::element_text(size=text.size))
  
  # Change color.palette (if necessary)
  if(!ggplot2::is.ggplot(bootega.obj$plot.typical.ega)){
    ic.plot <- suppressMessages(
      ic.plot + ggplot2::scale_color_manual(values = color_palette_EGA("rainbow", orig.wc),
                                            breaks = sort(orig.wc))
    )
  }else{
    if(bootega.obj$color.palette != "Set1"){
      ic.plot <- suppressMessages(
        ic.plot + ggplot2::scale_color_manual(values = color_palette_EGA(bootega.obj$color.palette, orig.wc),
                                              breaks = sort(orig.wc))
      )
    }
  }
  
  # Reverse ordering
  ic.plot <- ic.plot + ggplot2::scale_x_discrete(limits = rev(ic.plot$data$Item))
  
  
  if(plot.item.rep)
  {result$plot.itemStability <- ic.plot}
  
  #match row names to plot output
  itemCon <- itemCon[rev(match(ic.plot$data$Item,names(itemCon)))]
  
  #match row names to itemCon output
  itemLik <- as.data.frame(item.tab[match(names(itemCon),row.names(item.tab)),])
  
  #catch unidimensional structures (ugly band-aid fix)
  if(length(colnames(itemLik)) == 1){
    colnames(itemLik) <- 1
  }
  
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
    
    dims <- ncol(item.id.samps[[m]])
    
    if(dims!=max.wc)
    {
      if(any(colnames(item.id.samps[[m]])=="NA"))
      {item.id.samps[[m]] <- item.id.samps[[m]][,-which(colnames(item.id.samps[[m]])=="NA")]}
      
      diff <- max.wc - ncol(item.id.samps[[m]])
      
      diff.wc <- setdiff(seq(1,max.wc,1),unique(colnames(item.id.samps[[m]])))
      
      col.names <- c(colnames(item.id.samps[[m]]),paste(diff.wc))
      
      for(i in 1:diff)
      {item.id.samps[[m]] <- cbind(item.id.samps[[m]],rep(NA,nrow(item.id.samps[[m]])))}
      
      colnames(item.id.samps[[m]]) <- na.omit(col.names)
      
      item.id.samps[[m]] <- item.id.samps[[m]][,order(as.numeric(colnames(item.id.samps[[m]])))]
    }
  }
  
  #Unstandardized
  arr.func <- function(data)
  {
    # Get maximum dimensions
    r.dims <- max(sapply(data, dim)[1,])
    n.dims <- max(sapply(data, dim)[2,])
    
    arr <- array(NA, dim = c(r.dims, n.dims, length(data)))
    
    for(i in 1:length(data))
    {
      # Target matrix
      target.mat <- as.matrix(data[[i]])
      
      # Reorder based on bootega network
      target.mat <- target.mat[match(colnames(bootega.obj$EGA$network), row.names(target.mat)),]
      
      # Check for NAs
      if(any(is.na(target.mat))){
        NAs <- which(is.na(target.mat))
        
        names(target.mat)[NAs] <- colnames(bootega.obj$EGA$network)[NAs]
        
      }
      
      # Insert into array
      arr[,,i] <- target.mat
    }
    
    return(arr)
  }
  
  
  arr <- arr.func(item.id.samps)
  unstd.item.id <- round(apply(arr,1:2, mean, na.rm=TRUE),3)
  colnames(unstd.item.id) <- paste(seq(1,max(final.mat, na.rm = TRUE),1))
  row.names(unstd.item.id) <- colnames(bootega.obj$EGA$network)
  unstd.item.id[,1:length(uniq)] <- unstd.item.id[,uniq]
  colnames(unstd.item.id)[1:length(uni)] <- uni
  
  #let user know results are computed has ended
  message("done", appendLF = TRUE)
  
  if(ncol(unstd.item.id) > 1){
    
    unstd.item.id <- unstd.item.id[row.names(item.lik),]
    
    unstd.item.id <- unstd.item.id[,colnames(item.lik)]
  }
  
  unstd.item.id[which(item.lik=="")] <- ""
  
  unstd.item.id <- as.data.frame(unstd.item.id)
  
  #Unstandardize
  unstd.item.ident <- as.data.frame(cbind(orig.wc[row.names(unstd.item.id)], unstd.item.id))
  colnames(unstd.item.ident) <- c("Dimension",colnames(unstd.item.id))
  itemLoads <- unstd.item.ident[match(names(itemCon),row.names(unstd.item.ident)),]
  
  # Average replication in each dimension
  org <- orig.wc[match(names(itemCon),names(orig.wc))]
  
  # Initialize dimension replication vector
  dimRep <- numeric(length(uniq))
  
  # Compute average replication in each dimension
  for(i in 1:length(uniq))
  {dimRep[i] <- mean(itemCon[which(org == uni[i])])}
  
  # Name dimension replication vector
  names(dimRep) <- uni
  
  #message for additional item likelihoods
  if(ncol(itemLik)<max(final.mat,na.rm=TRUE))
  {message("\nLower the item.freq argument to view item frequencies for additional dimensions")}
  
  ##################################
  #### NETWORK LOADINGS RESULTS ####
  ##################################
  
  # Remove extra columns
  itemLoads[,-1] <- itemLoads[,colnames(itemLik)]
  
  # Remove NA from item loadings
  blank.itemLoads <- ifelse(as.matrix(itemLoads) == "NaN", "", as.matrix(itemLoads))
  blank.cols <- which(apply(blank.itemLoads, 2, function(x){all(x == "")}))
  if(length(blank.cols) != 0)
  {itemLoads <- itemLoads[,-blank.cols]}
  
  result$item.replication <- itemCon
  result$mean.dim.rep <- dimRep[order(names(dimRep))]
  result$item.dim.rep <- itemLik[,order(colnames(itemLik))]
  result$item.loadings <- itemLoads
  result$wc <- final.mat
  result$uniq.name <- uni
  result$uniq.num <- unique(num.comm)
  
  return(result)
}
#----
