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
# Updated 07.06.2021
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
  
  # Initialize plot list
  plots.ega <- list()
  plots.ega[[1]] <- plot(object.list[[1]], produce = FALSE)
  
  # Loop through matching 
  for(i in 2:length(object.list)){
    plots.ega[[i]] <- compare.EGA(object.list[[1]], object.list[[i]])[[2]]
  }
  
  # Set up grid return
  ggpubr::ggarrange(plotlist = plots.ega, ncol = columns, nrow = rows, labels = labels, label.x = 0.3)
  
}