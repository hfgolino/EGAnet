#' Visual Compares \code{\link{EGAnet}} plots
#' 
#' @description Identifies redundant nodes in the network based on several
#' measures. Computes the weighted topological overlap between
#' each node and every other node in the network. The weighted topological
#' overlap is implemented using the method from Nowick et al. (2009; see references)
#' and the function \link[wTO]{wTO} from the wTO package.
#'
#' @param ... \code{\link{EGAnet}} objects
#' 
#' @param labels Character vector.
#' Labels for each \code{\link{EGAnet}} object
#' 
#' @param rows Numeric.
#' Number of rows to spread plots across
#' 
#' @param columns Numeric.
#' Number of columns to spread plots down
#'
#' @return Visual comparison of \code{\link{EGAnet}} objects
#' 
#' @examples
#' # obtain SAPA items
#' items <- psychTools::spi[,c(11:20)]
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#
# Compare EGA plots function
# Updated 10.02.2022
compare.EGA.plots <- function(..., labels, rows, columns)
{
  # Obtain object list
  object.list <- list(...)
  
  # Identify EGA objects
  EGA.idx <- grep("EGA", unlist(lapply(object.list, class)))
  
  # Obtain EGA objects only
  object.list <- object.list[EGA.idx]
  
  # Obtain names
  name <- names(object.list)
  
  # Missing arguments
  if(missing(rows)){
    rows <- 1
  }
  
  if(missing(columns)){
    columns <- length(object.list)
  }
  
  if(missing(labels)){
    labels <- name
  }
  
  # Check for at least two objects
  if(length(object.list) < 2){
    stop("EGA plot comparisons require two or more EGA objects")
  }
  
  # Organize memberships
  for(i in 2:length(object.list)){
    
    # Target membership
    target.wc <- switch(
      class(object.list[[1]]),
      "EGA" = object.list[[1]]$wc,
      "bootEGA" = object.list[[1]]$typicalGraph$wc,
      "dynEGA" = object.list[[1]]$dynEGA$wc
    )
    
    # Covert membership
    convert.wc <- as.matrix(
      switch(
        class(object.list[[i]]),
        "EGA" = object.list[[i]]$wc,
        "bootEGA" = object.list[[i]]$typicalGraph$wc,
        "dynEGA" = object.list[[i]]$dynEGA$wc
      )
    )
    
    # Homogenize membership
    homogenized.wc <- as.vector(homogenize.membership(target.wc, convert.wc))
    
    # Replace membership
    if(class(object.list[[i]]) == "EGA"){
      object.list[[i]]$wc <- homogenized.wc
      names(object.list[[i]]$wc) <- colnames(object.list[[i]]$network)
    }else if(class(object.list[[i]]) == "bootEGA"){
      object.list[[i]]$typicalGraph$wc <- homogenized.wc
    }else if(class(object.list[[i]]$dynEGA$wc)){
      object.list[[i]]$typicalGraph$wc <- homogenized.wc
    }
    
  }
  
  # Initialize plot list
  plots.ega <- list()
  plots.ega <- suppressPackageStartupMessages(
    compare.plot.fix.EGA(object.list)
  )
  
  # Loop through matching 
  for(i in 2:length(object.list)){
    plots.ega[[i]] <- compare.EGA(plots.ega[[1]], plots.ega[[i]])[[2]]
  }
  
  # Set up grid return
  compare.plot <- ggpubr::ggarrange(plotlist = plots.ega, ncol = columns,
                    nrow = rows, labels = labels, label.x = 0.3)
  
  # Name plots
  names(plots.ega) <- labels
  
  # Results
  result <- list()
  result$comparison.plot <- compare.plot
  result$individual.plots <- plots.ega
  
  # Plot
  plot(compare.plot)
  
  return(result)
  
}