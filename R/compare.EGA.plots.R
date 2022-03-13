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
#' @param input_list List.
#' Bypasses \code{...} argument in favor of using a list
#' as an input
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
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}
#' 
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#'
#'  \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#'
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for
#' more details and examples}
#'
#' }
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
# Updated 03.03.2022
compare.EGA.plots <- function(
  ..., input_list = NULL, labels, rows, columns,
  plot.type = c("GGally", "qgraph"), plot.args = list()
)
{
  # Check for input list
  if(is.null(input_list)){
    object.list <- list(...)
  }else{
    object.list <- input_list
  }
  
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
  
  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{
    plot.type <- match.arg(plot.type)
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
    
    # Ensure names
    names(homogenized.wc) <- names(target.wc)
    
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
    compare.plot.fix.EGA(
      object.list,
      plot.type = plot.type,
      plot.args = plot.args
    )
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
