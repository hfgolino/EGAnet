#' Visually Compares \code{\link{EGAnet}} plots
#' 
#' @description Organizes EGA plots for comparison. Ensures that nodes are
#' placed in the same layout to maximize comparison. Community memberships
#' are also homogenized across EGA outputs to enhance interpretation
#'
#' @param ... \code{\link{EGAnet}} objects
#' 
#' @param input.list List.
#' Bypasses \code{...} argument in favor of using a list
#' as an input
#' 
#' @param base.plot Numeric.
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
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' }
#' 
#' (see \code{\link[GGally]{ggnet2}} for
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
#' # Obtain SAPA items
#' items <- psychTools::spi[,c(11:20)]
#' 
#' # Draw random samples
#' sample1 <- items[sample(1:nrow(items), 1000),]
#' sample2 <- items[sample(1:nrow(items), 1000),]
#' 
#' \dontrun{
#' # Estimate EGAs
#' ega1 <- EGA(sample1)
#' ega2 <- EGA(sample2)
#' 
#' # Compare EGAs via plot
#' compare_EGA_plots(
#'   ega1, ega2,
#'   base.plot = 1, # use "ega1" as base for comparison
#'   labels = c("Sample 1", "Sample 2"),
#'   rows = 1, columns = 2
#' )}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#
# Compare EGA plots ----
# Updated 14.07.2023
compare.EGA.plots <- function(
  ..., input.list = NULL, base.plot = 1,
  labels = NULL, rows = NULL, columns = NULL
)
{
  
  # Start with ellipse objects/arguments
  ellipse <- list(...)
  
  # Determine whether it's necessary to search for
  # plots in the ellipse
  if(is.null(input.list)){
    
    # Get classes for all input
    classes <- lapply(ellipse, class)
    
    # Determine which inputs are `EGA`
    ega_object <- lvapply(classes, function(x){grepl("EGA", x)})
    
    # Check for no `EGA` objects
    if(!any(ega_object)){
      stop(
        "The 'input.list' was `NULL` and there were no `EGA` class objects detected.",
        call. = FALSE
      )
    }
    
    # Extract `EGA` objects from list
    input.list <- ellipse[ega_object]
    ellipse <- ellipse[!ega_object]
    
  }
  
  # With the input list, there could be different types of
  # `EGA` objects... let's figure that out
  ega_classes <- cvapply(input.list, class)
  
  # Separate `dynEGA` from rest of `EGA` functions
  dynega_classes <- grepl("dynEGA", ega_classes)
  
  # Separate input list
  ega_input <- input.list[!dynega_classes]
  dynega_input <- input.list[dynega_classes]
  
  # For `EGA` input, extract `EGA` objects
  if(length(ega_input) != 0){
    ega_list <- lapply(ega_input, get_EGA_object)
  }else{ # Set up empty list
    ega_list <- list()
  }
  
  # For `dynEGA` input, extract `EGA` objects
  if(length(dynega_input) != 0){
    
    # Extract `dynEGA` objects
    dynega_objects <- ulapply(
      dynega_input, get_EGA_object, recursive = FALSE
    )
    
    # Check for `population` objects
    population_objects <- names(dynega_objects) == "population"
    
    # Separate `population` objects from the rest
    dynega_population <- dynega_objects[population_objects]
    dynega_other <- unlist( # Extract `EGA` objects from rest
      dynega_objects[!population_objects], recursive = FALSE
    )
      
    # Put all objects together
    dynega_list <- c(
      dynega_population, dynega_other
    )
    
  }else{ # Set up empty list
    dynega_list <- list()
  }
  
  # Compile all `EGA` objects
  input.list <- c(ega_list, dynega_list)
  
  # Organize input list with base plot
  if(is.character(base.plot)){
    base.plot <- which(names(input.list) == base.plot)
  }
  
  # Re-assemble input list
  input.list <- c(
    input.list[base.plot], input.list[-base.plot]
  )
  
  # Assign "labels"
  if(!is.null(labels)){
    names(input.list) <- labels
  }
  
  # Set up as `dynEGA.Individual` object to
  # leverage the code already implemented there
  class(input.list) <- "dynEGA.Individual"
  
  # For `ellipse`, get legacy arguments
  ellipse <- legacy_EGA_args(ellipse)
  
  # Add labels, rows, and columns
  ellipse[c("nrow", "ncol")] <- list(rows, columns)
  
  # Plot using all "IDs"
  return(
    do.call(
      plot,
      args = c(
        list(
          x = input.list,
          id = seq_along(input.list)
        ),
        ellipse
      )
    )
  )

}

