#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## S3Methods plot() // Updated 07.07.2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' S3Methods for Plotting
#'
#' @name plots
#'
#' @aliases
#' plot.CFA
#'
#' @usage
#' \method{plot}{CFA}(x, layout = "spring", vsize = 6, ...)
#'
#' @description Plots for \code{EGAnet} objects
#'
#' @param x Object from \code{EGAnet} package
#' 
#' @param vsize Numeric.
#' Size of vertices in \code{\link[EGAnet]{CFA}} plots.
#' Defaults to \code{6}
#'
#' @param layout Character.
#' Layout of plot (see \code{\link[semPlot]{semPaths}}).
#' Defaults to "spring"
#'
#' @param ncol Numeric.
#' Number of columns
#'
#' @param nrow Numeric.
#' Number of rows
#'
#' @param title Character.
#' Title of the plot.
#' Defaults to \code{""}
#'
#' @param id Numeric.
#' An integer or character indicating the ID of the individual to plot
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot. See 
#' \code{\link[GGally]{ggnet2}} for full list of arguments:
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
#' \item{\strong{\code{legend.names}}}
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
#' @param produce Boolean.
#' This argument is used internally.
#' Should plot be produced?
#' Defaults to \code{TRUE}
#' 
#'
#' @param ... Arguments passed on to
#'
#' \itemize{
#' 
#' \item{\code{\link[semPlot]{semPaths}}}
#' {Functions: CFA}
#'
#' }
#'
#' @return Plots of \code{EGAnet} object
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom graphics plot
#'
#' @export
# Plot invariance----
# Updated 07.07.2022
#' @export
plot.invariance <- function(
  x, title = "", labels = NULL,
  rows, columns,
  plot.args = list(), ...
)
{
  # Obtain structure
  structure <- x$memberships
  
  # Prepare EGA results for plots
  input_EGA <- lapply(x$groups$EGA, function(x){
    
    # Make class 'EGA'
    class(x) <- "EGA"
    
    # Return list
    return(x)
  
  })
  
  # Set structure
  input_EGA <- lapply(input_EGA, function(x){
    
    # Set structure
    x$wc <- structure
    
    # Return list
    return(x)
    
  })
  
  # Check for labels
  if(is.null(labels)){
    labels <- names(input_EGA)
  }
  
  # Check for rows
  if(missing(rows)){
    rows <- 1
  }
  
  # Check for columns
  if(missing(columns)){
    columns <- length(input_EGA)
  }
  
  # Set up plot arguments
  if(any(x$results$p <= .05)){
    
    # Check for "alpha" in plot.args
    if(!"alpha" %in% names(plot.args)){
      plot.args$alpha <- ifelse(
        x$results$p <= .05, .8, .2
      )
    }
    
  }
  
  # Obtain plots
  plots <- compare_EGA_plots(
    input.list = input_EGA,
    labels = labels, rows = rows,
    columns = columns,
    plot.args = plot.args
  )
  
}


#Plot CFA----
# Updated 02.05.2020
#' @export
plot.CFA <- function(x, layout = "spring", vsize = 6, ...) {
  semPlot::semPaths(x$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                    sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                    "std", cut = 0.5, ...)
}
