#' Visually Compares \code{\link{EGAnet}} plots
#' 
#' @description Organizes EGA plots for comparison. Ensures that nodes are
#' placed in the same layout to maximize comparison. Community memberships
#' are also homogenized across EGA outputs to enhance interpretation
#'
#' @param ... \code{\link{EGAnet}} objects
#' 
#' @param input_list List.
#' Bypasses \code{...} argument in favor of using a list
#' as an input
#' 
#' @param base_plot Numeric.
#' Plot to be used as the base for the configuration of the networks.
#' Uses the number of the order in which the plots are input.
#' Defaults to \code{1} or the first plot
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
# Updated 28.03.2022
compare.EGA.plots <- function(
  ..., input_list = NULL,
  base_plot = 1,
  labels, rows, columns,
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
  
  # Obtain base EGA
  base_EGA <- object.list[[base_plot]]
  
  # Comparison EGAs
  comparison_EGA <- object.list[-base_plot]
  
  # Set up number of communities (for legend later)
  num_wc <- numeric(length(object.list))
  num_wc[1] <- switch(
    class(base_EGA),
    "EGA" = length(na.omit(unique(base_EGA$wc))),
    "bootEGA" = length(na.omit(unique(base_EGA$typicalGraph$wc))),
    "dynEGA" = length(na.omit(unique(base_EGA$dynEGA$wc)))
  )
  
  # Organize memberships
  for(i in 1:length(comparison_EGA)){
    
    # Target membership
    target.wc <- switch(
      class(base_EGA),
      "EGA" = base_EGA$wc,
      "bootEGA" = base_EGA$typicalGraph$wc,
      "dynEGA" = base_EGA$dynEGA$wc
    )
    
    # Covert membership
    convert.wc <- as.matrix(
      switch(
        class(comparison_EGA[[i]]),
        "EGA" = comparison_EGA[[i]]$wc,
        "bootEGA" = comparison_EGA[[i]]$typicalGraph$wc,
        "dynEGA" = comparison_EGA[[i]]$dynEGA$wc
      )
    )
    
    # Homogenize membership
    homogenized.wc <- as.vector(homogenize.membership(target.wc, convert.wc))
    
    # Ensure names
    names(homogenized.wc) <- names(target.wc)
    
    # Replace membership
    if(class(comparison_EGA[[i]]) == "EGA"){
      comparison_EGA[[i]]$wc <- homogenized.wc
      names(comparison_EGA[[i]]$wc) <- colnames(comparison_EGA[[i]]$network)
    }else if(class(comparison_EGA[[i]]) == "bootEGA"){
      comparison_EGA[[i]]$typicalGraph$wc <- homogenized.wc
    }else if(class(comparison_EGA[[i]]$dynEGA$wc)){
      comparison_EGA[[i]]$typicalGraph$wc <- homogenized.wc
    }
    
    # Obtain number of communities
    num_wc[i + 1] <- length(na.omit(unique(homogenized.wc)))
  }
  
  # Reset object list
  for(i in 1:length(object.list)){
    
    if(i == 1){
      object.list[[1]] <- base_EGA
    }else{
      object.list[[i]] <- comparison_EGA[[i-1]]
    }
    
  }
  
  # Initialize plot list
  plots_ega <- list()
  plots_ega <- suppressPackageStartupMessages(
    compare.plot.fix.EGA(
      object.list,
      plot.type = plot.type,
      plot.args = plot.args
    )
  )

  # Loop through matching 
  for(i in 2:length(plots_ega)){
    plots_ega[[i]] <- compare.EGA(plots_ega[[1]], plots_ega[[i]])[[2]]
  }
  
  # Re-organize plot list with reference to base plot
  set_number <- 1:length(plots_ega) # obtain number of plots
  set_diff <- setdiff(set_number, base_plot) # remove base plot
  set_org <- c(base_plot, set_diff) # get set organization
  plots_ega <- plots_ega[set_org] # organize to original inputs
  
  # Re-organize number of communities for legend
  num_wc <- num_wc[set_org]
  
  # Obtain maximum number of communities
  if(!all(num_wc == max(num_wc, na.rm = TRUE))){
    max_wc <- which.max(num_wc)
  }else{
    max_wc <- base_plot
  }
  
  # Set up grid return
  compare.plot <- ggpubr::ggarrange(
    plotlist = plots_ega, ncol = columns,
    nrow = rows, labels = labels, label.x = 0.3,
    legend.grob = ggpubr::get_legend(plots_ega[[max_wc]]),
    common.legend = TRUE, legend = "right"
  )
  
  # Name plots
  names(plots_ega) <- labels
  
  # Results
  result <- list()
  result$comparison.plot <- compare.plot
  result$individual.plots <- plots_ega
  
  # Plot
  plot(compare.plot)
  
  return(result)
  
}
